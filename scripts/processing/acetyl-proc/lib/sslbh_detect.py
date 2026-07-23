"""
Steps 1 & 2: SSLBH/acetyltransferase detection at PC80 level.

Step 1 — Detection (three evidence sources from raw_hhsuite.tsv):
    ECOD:   db == 'ECOD', name contains 'Single-stranded left-handed beta-helix', prob >= 0.7
    PHROGs: db == 'PHROGS', name in acetyltransferase keywords, prob >= 0.9
    PFAM:   db == 'PFAM', name contains any of {PF01757, PF24894, PF21711, PF10129, PF07779}, prob >= 0.9
    detected_by column: comma-separated list of sources (e.g. 'ECOD,PHROGs')

Step 2 — Quality filter and K-locus assignment:
    Protein → genome mapping chain:
        proteinID → strip _PROTEIN_\\d+ suffix → prophageID
        prophageID → prophages_metadata.tsv (prophageID, genomeID) → genomeID
        genomeID → bacteria_metadata.tsv (MGG_SC, K_locus)
    Quality filter: n_sequences > 10 AND n_scs > 3
    K-locus assignment:
        most-frequent K_locus accounts for >= 50% of sequences → assign it
        otherwise → MANY-LOCI(K1,K2)  (two most frequent)
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ECOD_PROB_MIN   = 0.7
PHROGS_PROB_MIN = 0.9
PFAM_PROB_MIN   = 0.9

ECOD_NAME_SUBSTR  = "Single-stranded left-handed beta-helix"
PHROGS_NAMES      = {
    "acetyltransferase",
    "o-acetyltransferase",
    "o-antigen acetylase",
    "acyl-coa n-acyltransferase",
    "antimicrobial peptide resistance and lipid a acylation protein pagp",
    "fabg-like 3-oxoacyl-(acyl-carrier-protein) reductase",
}
PFAM_ACCESSIONS = {
    "PF01757",
    "PF24894",
    "PF21711",
    "PF10129",
    "PF07779",
}

MIN_SEQUENCES = 10
MIN_SCS       = 3

_PROPHAGE_ID_PATTERN = re.compile(r"_PROTEIN_\d+$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Step 1 — detect
# ---------------------------------------------------------------------------

def detect_sslbh_pc80(
    raw_hhsuite_tsv: Path,
    pc2proteins_tsv: Path,
) -> pd.DataFrame:
    """
    Detect PC80s with SSLBH / acetyltransferase evidence.

    Returns DataFrame with columns:
        pc80        — PC identifier at PCI80C50 level
        detected_by — comma-separated source list, e.g. 'ECOD', 'PHROGs,PFAM'
        protein_ids — list of protein IDs belonging to this PC80
    """
    hhsuite = pd.read_csv(raw_hhsuite_tsv, sep="\t")
    hhsuite.columns = [c.lower() for c in hhsuite.columns]

    _require_cols(hhsuite, ["query", "db", "name", "prob"], raw_hhsuite_tsv)

    ecod_pcs   = _filter_ecod(hhsuite)
    phrogs_pcs = _filter_phrogs(hhsuite)
    pfam_pcs   = _filter_pfam(hhsuite)

    all_pcs = ecod_pcs | phrogs_pcs | pfam_pcs
    print(
        f"  ECOD: {len(ecod_pcs)}  PHROGs: {len(phrogs_pcs)}  "
        f"PFAM: {len(pfam_pcs)}  union: {len(all_pcs)} PC80s"
    )

    rows = []
    for pc in sorted(all_pcs):
        sources = []
        if pc in ecod_pcs:   sources.append("ECOD")
        if pc in phrogs_pcs: sources.append("PHROGs")
        if pc in pfam_pcs:   sources.append("PFAM")
        rows.append({"pc80": pc, "detected_by": ",".join(sources)})
    detected = pd.DataFrame(rows)

    # Map PC80 → protein IDs
    pc2prot = _load_pc2proteins(pc2proteins_tsv)
    detected["protein_ids"] = detected["pc80"].map(pc2prot).apply(
        lambda v: v if isinstance(v, list) else []
    )

    n_no_proteins = (detected["protein_ids"].apply(len) == 0).sum()
    if n_no_proteins:
        print(f"  [warn] {n_no_proteins} PC80s have no protein IDs in {pc2proteins_tsv.name}")

    return detected


# ---------------------------------------------------------------------------
# Step 2 — K-locus assignment
# ---------------------------------------------------------------------------

def assign_klocus(
    detected: pd.DataFrame,
    prophages_metadata_tsv: Path,
    bacteria_metadata_tsv: Path,
) -> pd.DataFrame:
    """
    Apply quality filter and assign K-locus to each PC80.

    Adds columns: n_sequences, n_scs, k_locus
    Drops PC80s failing the quality filter (n_sequences <= 10 or n_scs <= 3).

    Args:
        detected: output of detect_sslbh_pc80()
        prophages_metadata_tsv: columns prophageID, genomeID
        bacteria_metadata_tsv:  columns genomeID, MGG_SC, K_locus
    """
    prophage_meta = pd.read_csv(prophages_metadata_tsv, sep="\t")
    prophage_meta.columns = [c.lower() for c in prophage_meta.columns]
    _require_cols(prophage_meta, ["prophageid", "genomeid"], prophages_metadata_tsv)
    prophage_to_genome: dict[str, str] = dict(
        zip(prophage_meta["prophageid"], prophage_meta["genomeid"])
    )

    bact_meta = pd.read_csv(bacteria_metadata_tsv, sep="\t")
    bact_meta.columns = [c.lower() for c in bact_meta.columns]
    _require_cols(bact_meta, ["genomeid", "k_locus"], bacteria_metadata_tsv)
    sc_col = next((c for c in bact_meta.columns if "sc" in c.lower() and "mg" in c.lower()), None)
    if sc_col is None:
        raise ValueError(f"{bacteria_metadata_tsv.name}: MGG_SC column not found; columns: {list(bact_meta.columns)}")
    genome_to_sc: dict[str, str]     = dict(zip(bact_meta["genomeid"], bact_meta[sc_col]))
    genome_to_klocus: dict[str, str] = dict(zip(bact_meta["genomeid"], bact_meta["k_locus"]))

    results = []
    for _, row in detected.iterrows():
        protein_ids = row["protein_ids"]
        n_sequences = len(protein_ids)

        # protein → genome via prophage
        genomes = []
        for pid in protein_ids:
            phage_id = _PROPHAGE_ID_PATTERN.sub("", pid)
            genome_id = prophage_to_genome.get(phage_id)
            if genome_id:
                genomes.append(genome_id)

        # quality filter
        sc_set = {genome_to_sc[g] for g in genomes if g in genome_to_sc and genome_to_sc[g]}
        n_scs = len(sc_set)

        if n_sequences <= MIN_SEQUENCES or n_scs <= MIN_SCS:
            continue

        # K-locus assignment
        kloci = [genome_to_klocus[g] for g in genomes
                 if g in genome_to_klocus and _is_valid_klocus(genome_to_klocus[g])]
        k_locus = _assign_dominant_klocus(kloci, n_sequences)

        results.append({
            "pc80":        row["pc80"],
            "detected_by": row["detected_by"],
            "protein_ids": protein_ids,
            "n_sequences": n_sequences,
            "n_scs":       n_scs,
            "k_locus":     k_locus,
        })

    out = pd.DataFrame(results)
    print(
        f"  Quality filter: {len(detected)} → {len(out)} PC80s "
        f"(dropped {len(detected) - len(out)} failing n_seq>{MIN_SEQUENCES} or n_sc>{MIN_SCS})"
    )
    return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _filter_ecod(df: pd.DataFrame) -> set[str]:
    mask = (
        (df["db"] == "ECOD") &
        df["name"].str.contains(ECOD_NAME_SUBSTR, na=False) &
        (df["prob"] >= ECOD_PROB_MIN)
    )
    return set(df.loc[mask, "query"].unique())


def _filter_phrogs(df: pd.DataFrame) -> set[str]:
    # PHROGs functional annotation is in the 'annot' column, not 'name'
    annot_col = "annot" if "annot" in df.columns else "name"
    mask = (
        (df["db"] == "PHROGS") &
        df[annot_col].str.lower().isin(PHROGS_NAMES) &
        (df["prob"] >= PHROGS_PROB_MIN)
    )
    return set(df.loc[mask, "query"].unique())


def _filter_pfam(df: pd.DataFrame) -> set[str]:
    pattern = "|".join(re.escape(acc) for acc in PFAM_ACCESSIONS)
    mask = (
        (df["db"] == "PFAM") &
        df["name"].str.contains(pattern, na=False) &
        (df["prob"] >= PFAM_PROB_MIN)
    )
    return set(df.loc[mask, "query"].unique())


def _load_pc2proteins(pc2proteins_tsv: Path) -> dict[str, list[str]]:
    """
    Load pc2proteins.tsv → dict {pc80: [protein_id, ...]}.
    Expects columns: PC, proteinID  (additional columns like 'repr' are ignored).
    """
    df = pd.read_csv(pc2proteins_tsv, sep="\t", dtype=str)
    df.columns = [c.lower() for c in df.columns]
    # Accept 'pc' or first column as PC identifier; 'proteinid' or second as protein
    pc_col  = "pc"       if "pc"       in df.columns else df.columns[0]
    pid_col = "proteinid" if "proteinid" in df.columns else df.columns[1]

    result: dict[str, list[str]] = {}
    for _, row_data in df.iterrows():
        result.setdefault(str(row_data[pc_col]), []).append(str(row_data[pid_col]))
    return result


def _is_valid_klocus(klocus: str | None) -> bool:
    if not klocus:
        return False
    s = str(klocus).strip().lower()
    return s not in ("", "nan", "none", "unknown", "na")


def _assign_dominant_klocus(kloci: list[str], total_sequences: int) -> str:
    if not kloci:
        return "UNKNOWN"
    counts = Counter(kloci)
    top = counts.most_common(2)
    top_locus, top_count = top[0]
    if top_count / total_sequences >= 0.5:
        return top_locus
    names = [k for k, _ in top]
    # Strip KL prefix for compact display: KL2 → K2, KL13 → K13
    short = [re.sub(r"^KL(\d+)$", r"K\1", n) for n in names]
    return f"MANY-LOCI({','.join(short)})"


def _require_cols(df: pd.DataFrame, cols: list[str], path: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: expected columns {missing}, got {list(df.columns)}")
