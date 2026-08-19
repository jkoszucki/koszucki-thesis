"""
Figure 4.2, Panel C (structure half) — acetyltransferase structural similarity network.

Supporting output, not a manuscript panel. Panel C itself is the 18-representative TM-score
matrix (`figure4_2_panelC_representatives.py` → `figure4_2_panelC_heatmap.py`); this module
links the same 73 proteins by *structural* similarity as a Cytoscape network, which is
useful for inspecting the two clusters but is not what the caption describes. It writes to
its own `tmalign/` subfolder so it keeps the conventional `node.tsv` / `edge.tsv` names
alongside the representatives tables.

It once paired with a sequence-similarity half (`figure4_2_panelC.py`, BLASTP identity),
which was cut from the thesis on 2026-07-27 and now lives at
`other/unused/figure4_2_panelC_sequence_network.py`.

Nothing is recomputed here. All-vs-all TM-align already runs in `processing/acetyl-proc`
(Step 4, see `structure_order.py`), which writes the pairwise scores and the clustering
this module reads. Edges are simply the pairs passing the TM-score cutoff.

Similarity: `tm_max` = max(tm_norm_1, tm_norm_2), the TM-score normalised by the shorter
chain — the same measure the clustering uses, so network and x-axis order agree.

Node set: the 73 proteins in the structural ordering — 69 K-locus candidates plus the four
S1_Table reference proteins. PROTEIN04_GWAS_AC_K30 is absent because it is excluded from
the ordering as a false-positive acetyltransferase.

Edge cutoff: TM-score >= 0.75 by default. This is well above the 0.5 same-fold threshold
used for clustering, so edges mark close structural correspondence rather than shared fold.

Outputs:
    plots_dir/figure4_2-panelC/tmalign/node.tsv
    plots_dir/figure4_2-panelC/tmalign/edge.tsv
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from enzymes_table import load_enzymes_table

TM_SCORE_MIN = 0.75

# Edge styling bands: two intervals above the cutoff, thin/thick in Cytoscape.
_TM_INTERVALS = [(0.90, "high (0.90-1.00)"), (TM_SCORE_MIN, "moderate (0.75-0.90)")]

_REFERENCE_TYPES = (
    (re.compile(r"_MOD_AC_",  re.IGNORECASE), "EXPERIMENTAL"),
    (re.compile(r"_GWAS_AC_", re.IGNORECASE), "GWAS_BEST"),
)


def _tm_interval(tm: float) -> str:
    for threshold, label in _TM_INTERVALS:
        if tm >= threshold:
            return label
    return "below cutoff"


def _node_type(protein_id: str) -> str:
    for pattern, label in _REFERENCE_TYPES:
        if pattern.search(protein_id):
            return label
    return "KLOCI"


def _specificity(protein_id: str) -> str:
    """K-locus the protein belongs to, from either naming scheme."""
    m = re.match(r"^(KL\d+)", protein_id)
    if m:
        return m.group(1)
    m = re.search(r"_K(\d+)$", protein_id)          # PROTEIN02_MOD_AC_K2 -> KL2
    return f"KL{m.group(1)}" if m else protein_id


def _label(protein_id: str, node_type: str) -> str:
    """Short display label: K-locus for candidates, active_AT_/gwas_AT_ for references."""
    if node_type == "EXPERIMENTAL":
        return f"active_AT_{_specificity(protein_id).replace('KL', 'K')}"
    if node_type == "GWAS_BEST":
        return f"gwas_AT_{_specificity(protein_id).replace('KL', 'K')}"
    return _specificity(protein_id)


def _build_nodes(order_df: pd.DataFrame, kloci_tsv: Path, enzymes_xlsx: Path | None) -> pd.DataFrame:
    kloci = pd.read_csv(kloci_tsv, sep="\t").set_index("protein_id")

    # Reference proteins are not in acetylases_kloci.tsv, so their lengths come from
    # S1_Table; without this every reference node has an empty protein_length and cannot
    # be sized alongside the candidates in Cytoscape.
    ref_lengths: dict[str, int] = {}
    if enzymes_xlsx is not None:
        enzymes = load_enzymes_table(enzymes_xlsx, verbose=False)
        for _, row in enzymes.iterrows():
            seq = str(row.get("sequence", ""))
            if seq and seq != "nan":
                ref_lengths[str(row["proteinid"])] = len(seq)

    rows = []
    for _, row in order_df.iterrows():
        pid       = row["protein_id"]
        node_type = _node_type(pid)
        evidence = [
            name for name, col in [
                ("experimental", "source_experimental"),
                ("blastp",       "source_blastp"),
                ("keyword",      "source_keyword"),
                ("foldseek",     "source_foldseek"),
            ] if bool(kloci.at[pid, col]) if pid in kloci.index
        ] if pid in kloci.index else [
            "experimental" if node_type == "EXPERIMENTAL" else "gwas_best"
        ]
        rows.append({
            "name":              pid,
            "proteinID":         pid,
            "type":              node_type,
            "specificity":       _specificity(pid),
            "label":             _label(pid, node_type),
            "evidence":          ",".join(evidence),
            "structure_cluster": int(row["structure_cluster"]),
            "order_index":       int(row["order_index"]),
            "protein_length":    int(kloci.at[pid, "protein_length"]) if pid in kloci.index
                                 else ref_lengths.get(pid, pd.NA),
        })
    return pd.DataFrame(rows)


def _build_edges(tmscore_tsv: Path, keep: set[str], tm_score_min: float) -> pd.DataFrame:
    pairs = pd.read_csv(tmscore_tsv, sep="\t")
    pairs = pairs[pairs["id1"].isin(keep) & pairs["id2"].isin(keep)]
    edges = pairs[pairs["tm_max"] >= tm_score_min].copy()
    edges["tm_interval"] = edges["tm_max"].apply(_tm_interval)
    edges["interaction"] = "structural_similarity"
    return edges.rename(columns={"id1": "source", "id2": "target"})[
        ["source", "target", "interaction", "tm_max", "tm_norm_1", "tm_norm_2",
         "rmsd", "tm_interval"]
    ].sort_values("tm_max", ascending=False)


def plot_figure4_2_panelC_tmalign(
    tmscore_tsv: Path,
    structure_order_tsv: Path,
    acetylases_kloci_tsv: Path,
    plots_dir: Path,
    enzymes_xlsx: Path | None = None,
    tm_score_min: float = TM_SCORE_MIN,
) -> None:
    """
    Build the Cytoscape node/edge tables for the structural similarity network.

    Args:
        tmscore_tsv:          cps_acetylases/tmalign/acetylases_kloci_tmscore.tsv
        structure_order_tsv:  cps_acetylases/tmalign/acetylases_kloci_structure_order.tsv
                              — defines the node set and supplies structure_cluster
        acetylases_kloci_tsv: cps_acetylases/acetylases_kloci.tsv — node metadata
        plots_dir:            scripts/figures/chapter4/plots
        enzymes_xlsx:         S1_Table.xlsx — supplies protein_length for the reference
                              proteins, which are absent from acetylases_kloci.tsv
        tm_score_min:         edge cutoff on tm_max (default 0.75)
    """
    order_df = pd.read_csv(structure_order_tsv, sep="\t").sort_values("order_index")
    nodes = _build_nodes(order_df, acetylases_kloci_tsv, enzymes_xlsx)
    edges = _build_edges(tmscore_tsv, set(nodes["name"]), tm_score_min)

    out_dir = Path(plots_dir) / "figure4_2-panelC" / "tmalign"
    out_dir.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(out_dir / "node.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "edge.tsv", sep="\t", index=False)

    connected = set(edges["source"]) | set(edges["target"])
    print(f"  {len(nodes)} nodes ({nodes['type'].value_counts().to_dict()})")
    print(f"  {len(edges)} edges at TM >= {tm_score_min} "
          f"({edges['tm_interval'].value_counts().to_dict()})")
    print(f"  {len(connected)} connected, {len(nodes) - len(connected)} isolated")
    print(f"  → {out_dir}/node.tsv, edge.tsv")
