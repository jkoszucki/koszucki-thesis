"""
Figure 3.1, Panel B — SGNH hydrolase similarity network (Cytoscape).

All-vs-all BLASTP between:
  - SGNH GWAS best predictor representative sequences (12 loci, precision ≥ 0.50:
    KL2, KL6, KL8, KL11, KL16, KL24, KL30, KL35, KL55, KL64, KL111, KL125)
  - 4 experimentally characterised deacetylases from enzymes.xlsx
    (modification == "deacetylation": PROTEIN06–PROTEIN09)

Node prediction_strength:
  - good   (precision ≥ 0.8): KL2, KL11, KL16, KL24, KL55, KL111
  - likely (precision < 0.8): KL6, KL8, KL30, KL35, KL64, KL125

Outputs (Cytoscape import files):
    plots_dir/figure3_1-panelB/node.tsv
    plots_dir/figure3_1-panelB/edge.tsv

Checkpoint: skips if both files already exist.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

# Coverage interval boundaries
EVALUE_MAX  = 1e-3
PIDENT_MIN  = 30.0   # %
QCOV_MIN    = 0.4    # minimum bidirectional coverage to retain edge

GOOD_PRECISION_MIN = 0.8


# ---------------------------------------------------------------------------
# Node assembly
# ---------------------------------------------------------------------------

def _load_gwas_nodes(gwas_sgnh_best_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(gwas_sgnh_best_tsv, sep="\t")
    df = df[df["representative_sequence"].notna() & (df["representative_sequence"] != "")]
    rows = []
    for _, row in df.iterrows():
        rows.append({
            "proteinID":   row["representative_protein_id"],
            "source":      "PREDICTION",
            "specificity": row["locus"],
            "label":       row["locus"],
            "seq":         row["representative_sequence"],
            "precision":   row["precision"],
            "recall":      row["recall"],
        })
    return pd.DataFrame(rows)


def _load_experimental_nodes(enzymes_xlsx: Path) -> pd.DataFrame:
    df = pd.read_excel(enzymes_xlsx, sheet_name="enzymes")
    df = df[df["modification"] == "deacetylation"].copy()
    rows = []
    for _, row in df.iterrows():
        protein_id = row["proteinid"]
        seq        = str(row["sequence"]).strip()
        if not seq:
            print(f"  [warn] {protein_id}: empty sequence — skipping")
            continue
        parts       = protein_id.split("_")
        specificity = parts[-1]
        category    = "_".join(parts[1:3])
        rows.append({
            "proteinID":   protein_id,
            "source":      f"EXPERIMENTAL_{category}",
            "specificity": specificity,
            "label":       protein_id,
            "seq":         seq,
        })
    return pd.DataFrame(rows)


def _build_nodes(gwas_sgnh_best_tsv: Path, enzymes_xlsx: Path) -> pd.DataFrame:
    gwas = _load_gwas_nodes(gwas_sgnh_best_tsv)
    exp  = _load_experimental_nodes(enzymes_xlsx)
    nodes = pd.concat([gwas, exp], ignore_index=True)
    nodes.index = nodes.index + 1
    nodes["node_id"] = nodes.index.astype(str) + "_" + nodes["label"]
    nodes["clean_id"] = nodes["node_id"].str.replace(r"^\d+_", "", regex=True)
    nodes["full_label"] = nodes.apply(
        lambda r: f"{r['proteinID']}_{r['label']}" if r["source"] == "PREDICTION" else r["proteinID"],
        axis=1,
    )
    nodes["prediction_strength"] = nodes.apply(
        lambda r: (
            "strong" if r["source"] == "PREDICTION"
                        and r.get("precision", 0) >= GOOD_PRECISION_MIN
            else "likely" if r["source"] == "PREDICTION"
            else "rbp"   if r["source"] == "EXPERIMENTAL_RBP_DAC"
            else pd.NA
        ),
        axis=1,
    )
    nodes["type"] = nodes.apply(
        lambda r: (
            "PREDICTION-STRONG"  if str(r["prediction_strength"]) == "strong"
            else "PREDICTION-LIKELY" if str(r["prediction_strength"]) == "likely"
            else "MOD" if r["source"] == "EXPERIMENTAL_MOD_DAC"
            else "RBP"
        ),
        axis=1,
    )
    return nodes[["node_id", "clean_id", "full_label", "proteinID", "type", "source",
                  "specificity", "label", "prediction_strength", "seq"]].copy()


# ---------------------------------------------------------------------------
# Edge generation (all-vs-all BLASTP)
# ---------------------------------------------------------------------------

def _sanitize(s: str) -> str:
    return s.replace(".", "_DOT_")


def _restore(s: str) -> str:
    return s.replace("_DOT_", ".")


def _generate_edges(nodes: pd.DataFrame, tmp_dir: Path) -> pd.DataFrame:
    nodes = nodes.copy()
    nodes["blast_id"] = nodes["node_id"].apply(_sanitize)
    blast_to_node = dict(zip(nodes["blast_id"], nodes["node_id"]))

    fasta_path = tmp_dir / "sequences.fasta"
    blast_out  = tmp_dir / "blastp_raw.tsv"

    with open(fasta_path, "w") as fh:
        for _, row in nodes.iterrows():
            fh.write(f">{row['blast_id']}\n{row['seq']}\n")

    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "prot"],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["blastp",
         "-query", str(fasta_path), "-db", str(fasta_path),
         "-out", str(blast_out),
         "-outfmt", "6 qseqid sseqid evalue pident bitscore nident length mismatch gapopen gaps qlen qstart qend slen sstart send"],
        check=True, capture_output=True,
    )

    cols = ["source", "target", "evalue", "pident", "bitscore",
            "nident", "length", "mismatch", "gapopen", "gaps",
            "qlen", "qstart", "qend", "slen", "sstart", "send"]
    df = pd.read_csv(blast_out, sep="\t", header=None, names=cols)

    df["source"] = df["source"].map(blast_to_node).fillna(df["source"].apply(_restore))
    df["target"] = df["target"].map(blast_to_node).fillna(df["target"].apply(_restore))

    df["qcov"] = np.round((df["qend"] - df["qstart"] + 1) / df["qlen"], 3)
    df["scov"] = np.round((df["send"] - df["sstart"] + 1) / df["slen"], 3)

    # Deduplicate symmetric hits; retain self-hits so singletons appear in Cytoscape
    df["_pair"] = df.apply(lambda r: "-".join(sorted([r["source"], r["target"]])), axis=1)
    df = df.drop_duplicates(subset="_pair")

    df = df[
        (df["qcov"] >= QCOV_MIN) &
        (df["scov"] >= QCOV_MIN) &
        (df["pident"] >= PIDENT_MIN) &
        (df["evalue"] <= EVALUE_MAX)
    ].copy()

    # All retained hits have qcov >= QCOV_MIN (0.4), so fragment-sharing is always "large"
    df["fragment-sharing"] = "large"

    df["seq-identity"] = pd.cut(
        df["pident"],
        bins=[0, 60, 100.001],
        labels=["low", "high"],
        right=False, include_lowest=True,
    )

    df["interaction"] = "pp"
    df = df.drop(columns=["_pair"])

    return df[["source", "target", "evalue", "pident", "bitscore",
               "nident", "length", "mismatch", "gapopen", "gaps",
               "qcov", "scov", "qlen", "slen",
               "fragment-sharing", "seq-identity", "interaction"]]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_figure3_1_panelB(
    gwas_sgnh_best_tsv: Path,
    enzymes_xlsx: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Build SGNH similarity network (GWAS best predictors + experimental deacetylases)
    and write Cytoscape node/edge TSV files.

    Args:
        gwas_sgnh_best_tsv: sgnh-hydrolases/gwas_sgnh_best.tsv
        enzymes_xlsx:       input/enzymes/enzymes.xlsx (rows with modification=="deacetylation")
        plots_dir:          output directory; files go into plots_dir/figure3_1-panelB/
        style:              cfg.style (optional, unused but kept for API consistency)
    """
    out_dir  = plots_dir / "figure3_1-panelB"
    node_out = out_dir / "node.tsv"
    edge_out = out_dir / "edge.tsv"

    if (node_out.exists() and node_out.stat().st_size > 0 and
            edge_out.exists() and edge_out.stat().st_size > 0):
        print("  figure3_1-panelB/node.tsv + edge.tsv exist — skipping (delete to recompute)")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    nodes = _build_nodes(gwas_sgnh_best_tsv, enzymes_xlsx)
    n_gwas = (nodes["source"] == "PREDICTION").sum()
    n_exp  = nodes["source"].str.startswith("EXPERIMENTAL").sum()
    print(f"  {len(nodes)} nodes: {n_gwas} GWAS predictors, {n_exp} experimental")

    with tempfile.TemporaryDirectory() as _tmp:
        edges = _generate_edges(nodes, Path(_tmp))

    print(f"  {len(edges)} edges after filtering (evalue ≤ {EVALUE_MAX})")

    nodes.to_csv(node_out, sep="\t", index=False)
    edges.to_csv(edge_out, sep="\t", index=False)
    print(f"  → figure3_1-panelB/node.tsv ({len(nodes)} rows)")
    print(f"  → figure3_1-panelB/edge.tsv ({len(edges)} rows)")
