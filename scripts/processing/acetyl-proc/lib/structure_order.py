"""
Structural ordering of the acetyltransferases plotted in figure 3.2 panel B — the 69
K-locus candidates plus the reference proteins passed as `extra_protein_ids` (currently
PROTEIN01/02/03/05, so 73 in total).

Runs all-vs-all TM-align over the AlphaFold3 models in S4_Data, clusters the resulting
TM-score matrix, and emits a dendrogram leaf order. Figure 3.2 panel B uses that order
for its x-axis, so the panel is arranged by an axis independent of the FoldSeek
annotations it plots — structure orders the proteins, annotation composition is what
the reader then judges against it.

The tmscore table is a cache, not a record of the plotted set: it may span more proteins
than are currently ordered (it still holds PROTEIN04/KL30, dropped from the figure after
it was computed). Only pairs within the requested set feed the clustering, and the
checkpoint recomputes when a *new* protein appears rather than when one is removed.

TM-align itself runs in the `tmtools` environment (which has no pandas/scipy), invoked
as a subprocess via scripts/helpers/tmalign_worker.py. Clustering happens here, in
`jkoszucki`.

Similarity: `tm_max`, i.e. max(tm_norm_1, tm_norm_2) — the TM-score normalised by the
shorter chain. This treats "a domain of A matches all of B" as similar, which suits a
set spanning 122-400+ residues. Distance = 1 - tm_max, average linkage (UPGMA).
Clusters are cut at TM-score 0.5, the conventional same-fold threshold.

Outputs (both under cps_acetylases/tmalign/):
    acetylases_kloci_tmscore.tsv         — long-format pairwise TM-scores (checkpointed)
    acetylases_kloci_structure_order.tsv — protein_id, order_index, structure_cluster
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

TM_SAME_FOLD = 0.5   # conventional TM-score threshold for "same fold"
_MODEL_SUFFIX = "_model_0.cif"


def _resolve_cifs(protein_ids: list[str], s4_data_dir: Path) -> tuple[dict[str, Path], list[str]]:
    """Return ({protein_id: CIF path}, protein_ids with no CIF).

    K-locus candidates are filed under their exact id (`KL4_03_wcuG_model_0.cif`) while
    the S1_Table proteins are filed lowercased (`protein02_mod_ac_k2_model_0.cif`), so
    both spellings are tried. The canonical (requested) id is what the caller gets back —
    filenames must not leak into the output tables.
    """
    found, missing = {}, []
    for pid in protein_ids:
        for candidate in (pid, pid.lower()):
            cif = s4_data_dir / f"{candidate}{_MODEL_SUFFIX}"
            if cif.is_file():
                found[pid] = cif
                break
        else:
            missing.append(pid)
    return found, missing


def _run_tmalign(cif_paths: list[Path], out_tsv: Path) -> None:
    worker = Path(__file__).resolve().parents[3] / "helpers" / "tmalign_worker.py"
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("\n".join(str(p) for p in cif_paths))
        list_path = Path(fh.name)
    try:
        subprocess.run(
            ["conda", "run", "-n", "tmtools", "python", str(worker),
             "--cif-list", str(list_path), "--out", str(out_tsv)],
            check=True,
        )
    finally:
        list_path.unlink(missing_ok=True)


def _similarity_matrix(pairs: pd.DataFrame, protein_ids: list[str]) -> pd.DataFrame:
    sim = pd.DataFrame(0.0, index=protein_ids, columns=protein_ids)
    wanted = set(protein_ids)
    # A tmscore table may span more proteins than are being ordered (e.g. after a
    # protein is dropped); ignore pairs outside the requested set rather than growing
    # the matrix with them.
    for row in pairs.itertuples(index=False):
        if row.id1 not in wanted or row.id2 not in wanted:
            continue
        sim.at[row.id1, row.id2] = row.tm_max
        sim.at[row.id2, row.id1] = row.tm_max
    for pid in protein_ids:
        sim.at[pid, pid] = 1.0
    return sim


def build_structure_order(
    acetylases_kloci_tsv: Path,
    s4_data_dir: Path,
    out_dir: Path,
    run_tmalign: bool = True,
    extra_protein_ids: list[str] | None = None,
) -> pd.DataFrame:
    """
    Compute the structural ordering of the K-locus acetyltransferase candidates.

    Args:
        acetylases_kloci_tsv: cps_acetylases/acetylases_kloci.tsv (source of protein_id)
        s4_data_dir:          input_dir/.../S4_Data (read-only AF3 models)
        out_dir:              cps_acetylases/ — both TSVs land in its tmalign/ subfolder
        run_tmalign:          False skips the TM-align subprocess and reuses the
                              existing tmscore TSV (fails if it is absent)
        extra_protein_ids:    reference proteins to place alongside the K-locus
                              candidates — the experimental and GWAS acetyltransferases
                              of S1_Table. They are ordered by the same clustering, so the
                              figure can show where a characterised enzyme falls among the
                              candidates. Any without an AF3 model in S4_Data are reported
                              and skipped. Callers pass PROTEIN01/02/03/05; PROTEIN04
                              (KL30) is deliberately omitted — see acetyl-proc/main.py.

    Returns:
        DataFrame with columns protein_id, order_index, structure_cluster.
    """
    protein_ids = pd.read_csv(acetylases_kloci_tsv, sep="\t")["protein_id"].tolist()
    if extra_protein_ids:
        protein_ids += [p for p in extra_protein_ids if p not in protein_ids]
    cif_map, missing = _resolve_cifs(protein_ids, s4_data_dir)
    if missing:
        print(f"  [warn] no AF3 model — excluded from the ordering: {', '.join(missing)}")
    print(f"  {len(cif_map)}/{len(protein_ids)} proteins have an AF3 model")

    resolved = list(cif_map)                       # canonical ids, in request order
    # The worker names rows after the CIF stem; map them back to canonical ids so the
    # lowercased S1_Table filenames never reach the output tables.
    stem_to_id = {p.name[: -len(_MODEL_SUFFIX)]: pid for pid, p in cif_map.items()}

    tmalign_dir = out_dir / "tmalign"
    tmalign_dir.mkdir(parents=True, exist_ok=True)
    tmscore_tsv = tmalign_dir / "acetylases_kloci_tmscore.tsv"

    def _covers_all(path: Path) -> bool:
        """True if an existing tmscore table already spans exactly this protein set."""
        done = pd.read_csv(path, sep="\t", usecols=["id1", "id2"])
        seen = {stem_to_id.get(s, s) for s in set(done.id1) | set(done.id2)}
        stale = set(resolved) - seen
        if stale:
            print(f"  {path.name} predates {len(stale)} protein(s) "
                  f"({', '.join(sorted(stale))}) — recomputing TM-align")
        return not stale

    has_tmscore = tmscore_tsv.exists() and tmscore_tsv.stat().st_size > 0
    if run_tmalign and not (has_tmscore and _covers_all(tmscore_tsv)):
        _run_tmalign(list(cif_map.values()), tmscore_tsv)
    elif has_tmscore:
        print(f"  {tmscore_tsv.name} covers all {len(resolved)} proteins — skipping TM-align")
    else:
        raise FileNotFoundError(f"{tmscore_tsv} absent and run_tmalign=False")

    pairs = pd.read_csv(tmscore_tsv, sep="\t")
    pairs["id1"] = pairs["id1"].map(lambda s: stem_to_id.get(s, s))
    pairs["id2"] = pairs["id2"].map(lambda s: stem_to_id.get(s, s))

    sim = _similarity_matrix(pairs, resolved)
    dist = 1.0 - sim
    for pid in resolved:
        dist.at[pid, pid] = 0.0

    z = linkage(squareform(dist.values, checks=False), method="average")
    leaves = dendrogram(z, no_plot=True)["leaves"]
    clusters = fcluster(z, t=1.0 - TM_SAME_FOLD, criterion="distance")

    order = pd.DataFrame({
        "protein_id":        [resolved[i] for i in leaves],
        "structure_cluster": [int(clusters[i]) for i in leaves],
    })
    order.insert(1, "order_index", range(len(order)))

    dest = tmalign_dir / "acetylases_kloci_structure_order.tsv"
    order.to_csv(dest, sep="\t", index=False)

    n_clusters = order["structure_cluster"].nunique()
    sizes = order["structure_cluster"].value_counts()
    n_singletons = int((sizes == 1).sum())
    print(f"  {n_clusters} structural clusters at TM ≥ {TM_SAME_FOLD} "
          f"({n_singletons} singletons, largest {int(sizes.max())})")
    print(f"  → {dest}")
    return order
