"""
Rebuild summary.tsv, UPLOAD batches, and 2_DRAWN_AND_ORGANISED symlinks.

S4_Data is now a single flat folder of {canonical_name}_model_0.cif files —
one real, finalised model per protein, independent of which category it
belongs to. It lives under
input_dir/supplementary-thesis/supplementary-data/S4_Data/ — **read-only**,
same rule as the rest of input_dir. This script only ever reads from there;
it never creates or modifies anything inside S4_Data. New AF3 predictions
must be added to S4_Data directly (as {canonical_name}_model_0.cif) outside
this pipeline before rerunning.

Only the derived, regenerable outputs remain writable, under output_dir/other/alphafold3/:
1_DRAWN/ (renders), 2_DRAWN_AND_ORGANISED/ (renders grouped by set), UPLOAD/ (batch JSONs for
proteins still missing from S4_Data), status_per_protein.tsv, summary.tsv.

Universe: all 12 canonical protein sets.
Status (based on presence in S4_Data):
  DONE     — {canonical_name}_model_0.cif exists in S4_Data and is non-empty
  MISSING  — no (non-empty) file in S4_Data yet
  EXCLUDED — in table but out of scope for AF3

Flow:
  Step A — resolve canonical_name per entry (cl-prefix fallback for the rare
           protein-cluster ID reused across clustering levels)
  Step B — check S4_Data membership (read-only; no writes — S4_Data is frozen)
  Step C — assign DONE / MISSING / EXCLUDED from that membership
  Step D — write status_per_protein.tsv and summary.tsv
  Step E — write UPLOAD/ batch JSONs for MISSING proteins
  Step F — update 2_DRAWN_AND_ORGANISED/ symlinks from 1_DRAWN/
  Step G — render nterminal structures not yet in 1_DRAWN/
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from enzymes_table import load_enzymes_table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import Config

cfg = Config()

BASE      = cfg.output_dir / "other" / "alphafold3"
RENDERED  = BASE / "1_DRAWN"
ORGANISED = BASE / "2_DRAWN_AND_ORGANISED"
MISSING   = BASE / "UPLOAD"

# Finalised thesis deliverable — read-only, never written to by this script.
# Flat: every protein is {canonical_name}_model_0.cif directly under S4_DATA.
S4_DATA = cfg.input_dir / "supplementary-thesis" / "supplementary-data" / "S4_Data"

DEPOL_DIR         = cfg.output_dir / "rbp_depolymerases"
DEACETYL_DIR      = cfg.output_dir / "rbp_deacetylases"
ACETYL_DIR        = cfg.output_dir / "cps_acetylases"
BESTPRED_DIR      = cfg.output_dir / "rbp_best_predictors"
MMSEQS            = cfg.gwas_path / "3_GWAS" / "1_INTERMEDIATE" / "2_MMSEQS"
ENZYMES_XLSX      = cfg.input_dir / "supplementary-thesis" / "supplementary-tables" / "S1_Table.xlsx"
NTERMINAL_RBP_DIR = cfg.output_dir / "other" / "n-terminal" / "prophages_rbps"

NTERMINAL_ANCHOR_COLOR = "#9d5a16"   # N-terminal anchor (residues 1-160)
NTERMINAL_ANCHOR_RESI  = "1-160"

BATCH_SIZE = 10

# Proteins from depolymerases_manualsearch_notproduced.tsv retained for AF3
NOTPRODUCED_RETAIN = frozenset({
    "1406_54", "1409_62", "876_57", "1248_60", "1190_71", "1344_72", "1668_49"
})


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _make_job(name: str, sequence: str, count: int = 3) -> dict:
    sequence = sequence.replace(" ", "")
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": sequence, "count": count}}],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def _seq_from_mmseqs(clustering_level: str, pc: str) -> str | None:
    fasta = MMSEQS / clustering_level / "alignments" / f"{pc}.fasta"
    if not fasta.exists():
        return None
    for record in SeqIO.parse(fasta, "fasta"):
        return str(record.seq).replace("-", "")
    return None


# ---------------------------------------------------------------------------
# Source table loaders
# Each returns list[dict] with keys:
#   protein_id, assigned_set, kl, ecod_type, seq, af3_count, excluded,
#   canonical_hint, cl (optional)
# ---------------------------------------------------------------------------

def _load_depol_set(
    fname: str,
    assigned_set: str,
    kl_col: str,
    exclude_mode: str = "none",
) -> list[dict]:
    """Load a depolymerase table (S1/S4A schema).

    exclude_mode:
      'none'       — all proteins included
      'all'        — all EXCLUDED
      'notproduced' — only NOTPRODUCED_RETAIN included; rest EXCLUDED
    """
    path = DEPOL_DIR / fname
    if not path.exists():
        print(f"  [warn] missing: {path}")
        return []
    df = pd.read_csv(path, sep="\t")
    seq_col = "protein seq" if "protein seq" in df.columns else "protein_seq"
    rows = []
    for _, row in df.iterrows():
        pid = str(row["proteinID"]).strip()
        kl  = str(row.get(kl_col, row.get("K_locus_host", "-"))).strip()
        kl  = kl if kl and kl != "nan" else "-"
        seq = str(row.get(seq_col, "")).strip()
        seq = seq if seq and seq != "nan" else None
        if exclude_mode == "all":
            excluded = True
        elif exclude_mode == "notproduced":
            excluded = pid not in NOTPRODUCED_RETAIN
        else:
            excluded = False
        rows.append({
            "protein_id":     pid,
            "assigned_set":   assigned_set,
            "kl":             kl,
            "ecod_type":      None,
            "seq":            seq,
            "af3_count":      3,
            "excluded":       excluded,
            "canonical_hint": (f"{pid}_{kl}" if kl != "-" else pid).replace("/", "-"),
            "cl":             "",
        })
    return rows


def _load_gwas_pc_set(
    fname: str,
    assigned_set: str,
    src_dir: Path,
    exclude_all: bool = False,
    af3_count: int = 3,
    seq_col: str | None = None,
    canonical_suffix: str = "",
) -> list[dict]:
    """Load a GWAS PC-based table (best_predictors, deacetylases_gwas, acetylases_gwas).

    `canonical_suffix` disambiguates a PC that is modelled twice at different oligomeric
    states. PC0675/KL30 and PC0915/KL111 are each both a best GWAS predictor (modelled as
    a homotrimer, like every other RBP predictor) and a GWAS acetyltransferase (modelled
    as a monomer, like every other acetyltransferase). Canonical names are `{PC}_{locus}`,
    so without a suffix the acetyltransferase entry for KL30 would resolve to the
    predictor's trimer CIF and the two would be silently conflated.
    """
    path = src_dir / fname
    if not path.exists():
        print(f"  [warn] missing: {path}")
        return []
    df = pd.read_csv(path, sep="\t")
    seen: set[tuple] = set()
    rows = []
    for _, row in df.iterrows():
        pc    = str(row["PC"]).strip()
        locus = str(row.get("locus", "-")).strip()
        cl    = str(row.get("clustering_level", "")).strip()
        ecod  = str(row.get("ecod_type", "")).strip() or None
        seq   = None
        if seq_col and seq_col in df.columns:
            s = str(row.get(seq_col, "")).strip()
            seq = s if s and s != "nan" else None
        if seq is None and cl:
            seq = _seq_from_mmseqs(cl, pc)
        key = (pc, locus)
        if key in seen:
            continue
        seen.add(key)
        rows.append({
            "protein_id":     pc,
            "assigned_set":   assigned_set,
            "kl":             locus,
            "ecod_type":      ecod,
            "seq":            seq,
            "af3_count":      af3_count,
            "excluded":       exclude_all,
            "canonical_hint": f"{pc}_{locus}{canonical_suffix}",
            "cl":             cl,
        })
    return rows


def _load_s1_table_set(
    enzymes_xlsx: Path,
    assigned_set: str,
    infix: str,
    af3_count: int = 1,
) -> list[dict]:
    """Load PROTEIN0X rows straight from S1_Table by proteinid infix.

    enzymes-proc splits S1_Table into the literature TSVs by `_MOD_AC_` / `_DAC_`, so rows
    with any other infix reach no output table and would be invisible to this index. The
    GWAS-predicted acetyltransferases (`_GWAS_AC_`) are exactly that case: real proteins
    with real AF3 models that belong in the render tree, but no literature table to load
    them from.
    """
    if not enzymes_xlsx.exists():
        print(f"  [warn] missing: {enzymes_xlsx}")
        return []
    df = load_enzymes_table(enzymes_xlsx, verbose=False)
    rows = []
    for _, row in df.iterrows():
        pid = str(row["proteinid"]).strip()
        if infix.lower() not in pid.lower():
            continue
        seq = str(row.get("sequence", "")).strip()
        kl  = str(row.get("ktype", "-")).strip()
        rows.append({
            "protein_id":     pid,
            "assigned_set":   assigned_set,
            "kl":             kl if kl and kl != "nan" else "-",
            "ecod_type":      None,
            "seq":            seq if seq and seq != "nan" else None,
            "af3_count":      af3_count,
            "excluded":       False,
            "canonical_hint": pid.lower(),
            "cl":             "",
        })
    return rows


def _load_literature_set(
    fname: str,
    assigned_set: str,
    src_dir: Path,
    af3_count: int = 1,
) -> list[dict]:
    """Load a literature enzyme table (PROTEIN0X IDs from enzymes.xlsx)."""
    path = src_dir / fname
    if not path.exists():
        print(f"  [warn] missing: {path}")
        return []
    df = pd.read_csv(path, sep="\t")
    rows = []
    for _, row in df.iterrows():
        pid = str(row["proteinid"]).strip()
        seq = str(row.get("sequence", "")).strip()
        seq = seq if seq and seq != "nan" else None
        kl  = str(row.get("ktype", "-")).strip()
        kl  = kl if kl and kl != "nan" else "-"
        rows.append({
            "protein_id":     pid,
            "assigned_set":   assigned_set,
            "kl":             kl,
            "ecod_type":      None,
            "seq":            seq,
            "af3_count":      af3_count,
            "excluded":       False,
            "canonical_hint": pid.lower(),
            "cl":             "",
        })
    return rows


def _load_kloci_set(fname: str, assigned_set: str, exclude_all: bool = False) -> list[dict]:
    """Load the K-loci acetyltransferases table."""
    path = ACETYL_DIR / fname
    if not path.exists():
        print(f"  [warn] missing: {path}")
        return []
    df = pd.read_csv(path, sep="\t")
    rows = []
    for _, row in df.iterrows():
        pid = str(row["protein_id"]).strip()
        seq = str(row.get("sequence", "")).strip()
        seq = seq if seq and seq != "nan" else None
        rows.append({
            "protein_id":     pid,
            "assigned_set":   assigned_set,
            "kl":             "-",
            "ecod_type":      None,
            "seq":            seq,
            "af3_count":      1,
            "excluded":       exclude_all,
            "canonical_hint": pid,
            "cl":             "",
        })
    return rows


def _load_nterminal_set() -> list[dict]:
    """Load nterminal RBP proteins from other/n-terminal/prophages_rbps/.

    Each *.fasta file is one homotrimer job, except rbps.fasta — a bundled
    multi-sequence file (all RBPs concatenated), not a per-protein job; excluded.
    """
    if not NTERMINAL_RBP_DIR.exists():
        print(f"  [warn] missing: {NTERMINAL_RBP_DIR}")
        return []
    rows = []
    for fasta_path in sorted(NTERMINAL_RBP_DIR.glob("*.fasta")):
        if fasta_path.stem.lower() == "rbps":
            continue  # bundled multi-sequence file, not a single protein job
        pid = fasta_path.stem
        record = next(SeqIO.parse(fasta_path, "fasta"))
        seq = str(record.seq).replace("-", "")
        rows.append({
            "protein_id":     pid,
            "assigned_set":   "nterminal",
            "kl":             "-",
            "ecod_type":      None,
            "seq":            seq,
            "af3_count":      3,
            "excluded":       False,
            "canonical_hint": pid.lower(),
            "cl":             "",
        })
    print(f"  nterminal: {len(rows)} RBP proteins")
    return rows


# ---------------------------------------------------------------------------
# S4_Data lookup — flat, case-insensitive, with cl-prefix fallback for the
# rare protein-cluster ID reused identically across clustering levels
# (e.g. PC0915_KL111 -> PCI80C80_PC0915_KL111).
# ---------------------------------------------------------------------------

def _find_in_s4data(entry: dict, s4data_lower: dict[str, Path]) -> tuple[str, Path | None]:
    """Return (canonical_name, cif_path | None)."""
    hint = entry["canonical_hint"]
    cl = entry.get("cl", "")

    stem = f"{hint}_model_0"
    hit = s4data_lower.get(stem.lower())
    if hit is not None:
        return hint, hit

    if cl:
        prefixed = f"{cl}_{hint}"
        hit = s4data_lower.get(f"{prefixed}_model_0".lower())
        if hit is not None:
            return prefixed, hit

    return hint, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading source tables…")

    all_entries: list[dict] = []

    # Depolymerase sets
    all_entries += _load_depol_set("depolymerases_virulent_active.tsv",           "depolymerases_virulent_active",           "K_locus_specificity")
    all_entries += _load_depol_set("depolymerases_manualsearch_active.tsv",        "depolymerases_manualsearch_active",        "K_locus_host")
    all_entries += _load_depol_set("depolymerases_gwas_active.tsv",                "depolymerases_gwas_active",                "K_locus_host")
    all_entries += _load_depol_set("depolymerases_manualsearch_inactive.tsv",      "depolymerases_manualsearch_inactive",      "K_locus_host",       exclude_mode="all")
    all_entries += _load_depol_set("depolymerases_manualsearch_notproduced.tsv",   "depolymerases_manualsearch_notproduced",   "K_locus_host",       exclude_mode="notproduced")
    all_entries += _load_depol_set("depolymerases_gwas.tsv",                       "depolymerases_gwas",                       "K_locus_host")

    # Best predictors (mix of SSRBH + SGNH; color resolved per ecod_type)
    all_entries += _load_gwas_pc_set("best_predictors_gwas.tsv",  "best_predictors_gwas",  BESTPRED_DIR,  af3_count=3)

    # Deacetylases (SGNH → trimers)
    all_entries += _load_gwas_pc_set("deacetylases_gwas.tsv",  "deacetylases_gwas",  DEACETYL_DIR, af3_count=3, seq_col="representative_sequence")
    all_entries += _load_literature_set("deacetylases_literature.tsv",  "deacetylases_literature",  DEACETYL_DIR, af3_count=3)

    # Acetyltransferases (SSLBH → monomers)
    all_entries += _load_literature_set("acetylases_literature_active.tsv",  "acetylases_literature_active",  ACETYL_DIR, af3_count=1)
    # GWAS-predicted acetyltransferases (KL30, KL111) — monomers, straight from S1_Table:
    # enzymes-proc keeps `_GWAS_AC_` rows out of the literature tables, so without this
    # they have AF3 models but no place in the index or the organised render tree.
    all_entries += _load_s1_table_set(ENZYMES_XLSX, "acetylases_gwas_putative", "_GWAS_AC_", af3_count=1)
    all_entries += _load_kloci_set("acetylases_kloci.tsv",   "acetylases_kloci",  exclude_all=False)
    # acetylases_gwas is excluded — indexed for completeness, no AF3 models requested.
    # The `_monomer` suffix is kept so the set stays correct if it is ever re-enabled:
    # PC0675/KL30 and PC0915/KL111 are also best_predictors_gwas entries, modelled as
    # homotrimers under `{PC}_{locus}`, and would otherwise claim those trimer CIFs.
    all_entries += _load_gwas_pc_set("acetylases_gwas.tsv",  "acetylases_gwas",  ACETYL_DIR,
                                     exclude_all=True, af3_count=1, canonical_suffix="_monomer")

    # N-terminal anchor RBPs (homotrimers; domain-colored on render)
    all_entries += _load_nterminal_set()

    print(f"  Total entries loaded: {len(all_entries)}")

    # -----------------------------------------------------------------------
    # Step A/B — resolve canonical_name for each entry against S4_Data (flat,
    # case-insensitive lookup with cl-prefix fallback)
    # -----------------------------------------------------------------------
    print(f"\nReading {S4_DATA} (read-only, S4_Data is frozen)…")
    s4data_files = list(S4_DATA.glob("*_model_0.cif"))
    s4data_lower: dict[str, Path] = {f.stem.lower(): f for f in s4data_files}
    print(f"  S4_Data: {len(s4data_files)} CIF files present")

    for entry in all_entries:
        if entry["excluded"]:
            entry["canonical_name"] = entry["canonical_hint"]
            entry["source_cif"] = None
        else:
            canonical, cif = _find_in_s4data(entry, s4data_lower)
            entry["canonical_name"] = canonical
            entry["source_cif"] = cif if (cif is not None and cif.stat().st_size > 0) else None

    # -----------------------------------------------------------------------
    # Step C — assign status from S4_Data presence
    # -----------------------------------------------------------------------
    rows: list[dict] = []
    tbd_jobs: list[tuple[str, str, int]] = []

    for entry in all_entries:
        canonical = entry["canonical_name"]
        if entry["excluded"]:
            status = "EXCLUDED"
        elif entry["source_cif"] is not None:
            status = "DONE"
        else:
            status = "MISSING"
        rows.append({
            "proteinID":      entry["protein_id"],
            "assigned_set":   entry["assigned_set"],
            "status":         status,
            "canonical_name": canonical,
            "ecod_type":      entry.get("ecod_type"),
        })
        if status == "MISSING" and entry.get("seq"):
            tbd_jobs.append((canonical, entry["seq"], entry["af3_count"]))

    # -----------------------------------------------------------------------
    # Step D — write status_per_protein.tsv and summary.tsv
    # -----------------------------------------------------------------------
    df = pd.DataFrame(rows, columns=["proteinID", "assigned_set", "status", "canonical_name", "ecod_type"])
    df.to_csv(BASE / "status_per_protein.tsv", sep="\t", index=False)
    done_n     = (df["status"] == "DONE").sum()
    missing_n  = (df["status"] == "MISSING").sum()
    excluded_n = (df["status"] == "EXCLUDED").sum()
    print(f"\nstatus_per_protein.tsv → {len(df)} rows  ({done_n} DONE, {missing_n} MISSING, {excluded_n} EXCLUDED)")

    summary_rows = []
    for s, g in df.groupby("assigned_set"):
        summary_rows.append({
            "assigned_set": s,
            "DONE":     int((g["status"] == "DONE").sum()),
            "MISSING":  int((g["status"] == "MISSING").sum()),
            "EXCLUDED": int((g["status"] == "EXCLUDED").sum()),
        })
    summary_df = pd.DataFrame(summary_rows, columns=["assigned_set", "DONE", "MISSING", "EXCLUDED"])
    summary_df.to_csv(BASE / "summary.tsv", sep="\t", index=False)
    for _, r in summary_df.iterrows():
        print(f"  {r['assigned_set']:<45}  DONE={r['DONE']}  MISSING={r['MISSING']}  EXCLUDED={r['EXCLUDED']}")

    # -----------------------------------------------------------------------
    # Step E — write UPLOAD/ batch JSONs (cleared and rewritten)
    # -----------------------------------------------------------------------
    print("\nWriting UPLOAD/ batches…")
    if MISSING.exists():
        for old in MISSING.glob("batch_*.json"):
            old.unlink()
    MISSING.mkdir(exist_ok=True, parents=True)

    seen_missing: set[str] = set()
    unique_jobs: list[tuple[str, str, int]] = []
    for canonical, seq, count in tbd_jobs:
        if canonical not in seen_missing:
            seen_missing.add(canonical)
            unique_jobs.append((canonical, seq, count))

    batches = [unique_jobs[i:i + BATCH_SIZE] for i in range(0, len(unique_jobs), BATCH_SIZE)]
    for n, batch in enumerate(batches, start=1):
        jobs = [_make_job(name, seq, count) for name, seq, count in batch]
        (MISSING / f"batch_{n:03d}.json").write_text(json.dumps(jobs, indent=2))
    print(f"  {len(unique_jobs)} unique MISSING proteins → {len(batches)} batch files")

    # -----------------------------------------------------------------------
    # Step F — update 2_DRAWN_AND_ORGANISED (recreated from scratch from 1_DRAWN)
    # -----------------------------------------------------------------------
    print("\nUpdating 2_DRAWN_AND_ORGANISED…")
    if ORGANISED.exists():
        shutil.rmtree(ORGANISED)
    ORGANISED.mkdir()

    organised_count = 0
    for _, row in df[df["status"] == "DONE"].iterrows():
        set_dir = ORGANISED / row["assigned_set"]
        set_dir.mkdir(exist_ok=True)
        png = RENDERED / f"{row['canonical_name']}_model_0.png"
        if png.exists():
            link = set_dir / png.name
            if not link.exists():
                link.symlink_to(png)
                organised_count += 1
    print(f"  {organised_count} symlinks created")

    # -----------------------------------------------------------------------
    # Step G — render nterminal structures with N-terminal anchor coloring
    #           residues 1-160: #9d5a16 (anchor); rest: gray
    #           Skips proteins whose PNG already exists in 1_DRAWN/.
    # -----------------------------------------------------------------------
    print("\nRendering nterminal structures …")
    nterminal_done = df[(df["assigned_set"] == "nterminal") & (df["status"] == "DONE")]
    if nterminal_done.empty:
        print("  No nterminal CIFs available yet (add to S4_Data, rerun index).")
    else:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from structure_renderer import StructureRenderer
        renderer = StructureRenderer(
            render_mode="accurate",
            color_mode="domain",
            domain_colors=[(NTERMINAL_ANCHOR_RESI, NTERMINAL_ANCHOR_COLOR)],
            gray_color="#cfcfce",
            orient_nterm=True,
        )
        cif_paths = [S4_DATA / f"{row['canonical_name']}_model_0.cif"
                     for _, row in nterminal_done.iterrows()]
        renderer.render_to_dir(cif_paths, RENDERED)

    print(f"\nDone. All outputs written to {BASE}/")


if __name__ == "__main__":
    main()
