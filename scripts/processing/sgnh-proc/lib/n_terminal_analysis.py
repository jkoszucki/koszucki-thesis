"""
Analyse BLAST hits for the shared N-terminal fragment.

Reads sgnh-hydrolases/n-terminal/raw_blast.tsv and produces four tables plus
prophage files and AF3 upload:

  1_high_identity.tsv               — evalue <= 1e-3, qcov >= 0.90, pident >= 0.90;
                                       protein + prophageID + wgrr clustering columns
  2_representatives.tsv             — 2–3 per wgrr50 cluster, preferring distinct wgrr90
  3_one_per_cluster.tsv             — single phage per wgrr50 cluster
  4_genome_metadata.tsv             — table 3 + genomeID + isolation/environment metadata
  prophages/{ID}.gb                 — GenBank files for SELECTED_PROPHAGES
  prophages_rbps/{protein_id}.fasta — one FASTA per RBP defined in SELECTED_PROTEINS

Checkpoints: each step skips if its output already exists.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

EVALUE_MAX  = 1e-3
QCOV_MIN    = 0.90
PIDENT_MIN  = 0.90

CLUSTERING_COLS = ["prophageID", "wgrr50", "wgrr70", "wgrr80", "wgrr90", "wgrr95", "wgrr99"]

# Four prophages selected from distinct wGRR50 clusters to illustrate that the
# shared N-terminal anchor is carried by tail fibre proteins with structurally
# distinct receptor-binding domains across diverse host species and sequence types.
SELECTED_PROPHAGES = [
    "KPN_B2_PHAGE002_M",   # K. pneumoniae ST200,     KL13, wGRR50=PV013
    "KVV_B3_PHAGE169_L",   # K. variicola  ST7019,    KL19, wGRR50=PV004
    "KVV_B1_PHAGE242_M",   # K. variicola  ST207-1LV, KL64, wGRR50=PV027
    "KPN_B18_PHAGE103_M",  # K. pneumoniae ST37,      KL14, wGRR50=PV018
]

# Tail fibre proteins per selected prophage carrying the shared N-terminal anchor.
# Protein numbers are 1-based CDS indices within the prophage GenBank record.
SELECTED_PROTEINS: dict[str, list[int]] = {
    "KPN_B2_PHAGE002_M":  [55, 58],
    "KVV_B3_PHAGE169_L":  [71, 72],
    "KVV_B1_PHAGE242_M":  [14, 11],
    "KPN_B18_PHAGE103_M": [13],
}

GENOME_META_COLS = [
    "genomeID", "collection", "kaspah_infection", "kaspah_specimen",
    "klebpavia_group_summary", "species", "ST", "K_locus", "K_type",
]


def _load_filtered(raw_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(raw_tsv, sep="\t")
    df["qcov"]        = (df["qend"] - df["qstart"] + 1) / df["qlen"]
    df["scov"]        = (df["send"] - df["sstart"] + 1) / df["slen"]
    df["pident_frac"] = df["pident"] / 100
    return df[df["evalue"] <= EVALUE_MAX].copy()


def _high_identity(df: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    hi = df[(df["qcov"] >= QCOV_MIN) & (df["pident_frac"] >= PIDENT_MIN)].copy()
    hi["prophageID"] = hi["sseqid"].str.replace(r"_PROTEIN_\d+$", "", regex=True)
    hi = hi[["sseqid", "prophageID"]].rename(columns={"sseqid": "protein"})
    return hi.merge(meta[CLUSTERING_COLS], on="prophageID", how="left")


def _select_representatives(hi: pd.DataFrame) -> pd.DataFrame:
    selected = []
    for _, grp in hi.groupby("wgrr50"):
        picked_grp90 = set()
        reps = []
        for _, row in grp.iterrows():
            if row["wgrr90"] not in picked_grp90:
                reps.append(row)
                picked_grp90.add(row["wgrr90"])
            if len(reps) == 3:
                break
        if len(reps) < 2:
            for _, row in grp.iterrows():
                if row["protein"] not in [r["protein"] for r in reps]:
                    reps.append(row)
                if len(reps) == 2:
                    break
        selected.extend(reps)
    return pd.DataFrame(selected)


def _build_prophage_genbanks(genbank_dir: Path, out_dir: Path) -> None:
    """Copy GenBank files for SELECTED_PROPHAGES to out_dir. Skips existing files."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for prophage_id in SELECTED_PROPHAGES:
        src = genbank_dir / f"{prophage_id}.gb"
        dst = out_dir / f"{prophage_id}.gb"
        if dst.exists():
            continue
        if not src.exists():
            print(f"  [warn] GenBank not found: {src}")
            continue
        shutil.copy2(src, dst)
        print(f"  {prophage_id}.gb → prophages/")


def _build_rbps_fastas(raw_tsv: Path, genbank_dir: Path, out_dir: Path) -> None:
    """
    Write one FASTA per RBP in SELECTED_PROTEINS to out_dir ({protein_id}.fasta).

    Sequences are taken from raw_blast.tsv (sseq column) where available;
    proteins not in BLAST hits are extracted from GenBank by 1-based CDS index.
    Checkpoint: skips proteins whose .fasta file already exists.
    """
    from Bio import SeqIO

    out_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(raw_tsv, sep="\t")
    seq_map: dict[str, str] = raw.set_index("sseqid")["sseq"].to_dict()

    for prophage_id, protein_nums in SELECTED_PROTEINS.items():
        gb_rec = None
        cdss = None
        for n in protein_nums:
            pid = f"{prophage_id}_PROTEIN_{n}"
            fasta_path = out_dir / f"{pid}.fasta"
            if fasta_path.exists():
                print(f"  {pid}.fasta exists — skipping")
                continue
            if pid in seq_map:
                seq = seq_map[pid]
            else:
                if gb_rec is None:
                    gb_path = genbank_dir / f"{prophage_id}.gb"
                    gb_rec  = SeqIO.read(gb_path, "genbank")
                    cdss    = [f for f in gb_rec.features if f.type == "CDS"]
                seq = cdss[n - 1].qualifiers.get("translation", [""])[0]
            fasta_path.write_text(f">{pid}\n{seq}\n")
            print(f"  {pid}.fasta — {len(seq)} aa")


def analyse_n_terminal(
    prophage_metadata_tsv: Path,
    bacteria_metadata_tsv: Path,
    sgnh_dir: Path,
    genbank_dir: Path | None = None,
) -> None:
    """
    Run the N-terminal BLAST analysis.

    Args:
        prophage_metadata_tsv: path to prophages_metadata.tsv
        bacteria_metadata_tsv: path to bacteria_metadata.tsv
        sgnh_dir:              sgnh-hydrolases/ root (or any dir containing n-terminal/)
        genbank_dir:           directory of prophage .gb files; required for GenBank
                               copying and rbps.fasta extraction
    """
    out_dir = sgnh_dir / "n-terminal"
    checkpoint = out_dir / "4_genome_metadata.tsv"
    if checkpoint.exists():
        print(f"  4_genome_metadata.tsv exists — skipping (delete to recompute)")
        return

    raw_tsv = out_dir / "raw_blast.tsv"
    if not raw_tsv.exists():
        raise FileNotFoundError(f"raw_blast.tsv not found: {raw_tsv}")

    meta = pd.read_csv(prophage_metadata_tsv, sep="\t")

    # Step 1 — filter by evalue, qcov, pident; map clustering
    df = _load_filtered(raw_tsv)
    hi = _high_identity(df, meta)
    hi.to_csv(out_dir / "1_high_identity.tsv", sep="\t", index=False)
    print(f"  1_high_identity.tsv — {len(hi)} hits (evalue <= {EVALUE_MAX}, "
          f"qcov >= {QCOV_MIN}, pident >= {PIDENT_MIN})")

    # Step 2 — 2–3 representatives per wgrr50, preferring distinct wgrr90
    sel = _select_representatives(hi)
    sel.to_csv(out_dir / "2_representatives.tsv", sep="\t", index=False)
    print(f"  2_representatives.tsv — {len(sel)} entries across "
          f"{sel['wgrr50'].nunique()} wgrr50 clusters")

    # Step 3 — one phage per wgrr50
    one = sel.groupby("wgrr50", sort=False).first().reset_index()
    one.to_csv(out_dir / "3_one_per_cluster.tsv", sep="\t", index=False)
    print(f"  3_one_per_cluster.tsv — {len(one)} phages (one per wgrr50)")

    # Step 4 — map genomeID + isolation/environment metadata
    bact = pd.read_csv(bacteria_metadata_tsv, sep="\t")
    # prophages_metadata already loaded as meta; get genomeID per prophageID
    prophage_to_genome = meta[["prophageID", "genomeID"]].drop_duplicates()
    out4 = one.merge(prophage_to_genome, on="prophageID", how="left")
    out4 = out4.merge(bact[GENOME_META_COLS], on="genomeID", how="left")
    out4.to_csv(checkpoint, sep="\t", index=False)
    print(f"  4_genome_metadata.tsv — {len(out4)} rows")
    print(out4[["wgrr50", "prophageID", "genomeID", "collection",
                "kaspah_specimen", "species", "ST", "K_locus"]].to_string(index=False))

    # Step 5 — copy GenBank files and build per-protein FASTAs for selected prophages
    if genbank_dir is not None:
        prophages_dir = out_dir / "prophages"
        rbps_dir      = out_dir / "prophages_rbps"
        print(f"  Copying GenBank files for {len(SELECTED_PROPHAGES)} selected prophages …")
        _build_prophage_genbanks(genbank_dir, prophages_dir)
        print(f"  Building prophages_rbps/ (one .fasta per RBP) …")
        _build_rbps_fastas(raw_tsv, genbank_dir, rbps_dir)
    else:
        print(f"  [skip] GenBank copy + prophages_rbps/ (genbank_dir not provided)")
