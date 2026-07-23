"""
Figure 3.2, Panel B — Acetyltransferase sequence similarity network (Cytoscape).

All-vs-all BLASTP between three node sets:
  - 69 candidate K-locus acetyltransferases (acetylases_kloci.tsv)
  - 3 experimentally confirmed CPS O-acetyltransferases (acetylases_literature_active.tsv):
      K1 WcsT/NeuO-like, K2 WcsU, K57 orf13
  - 2 GWAS best-performing acetyltransferase predictions (FoldSeek-identified):
      KL30 (PC0675, PCI80C80) and KL111 (PC0915, PCI80C80)

Edge thresholds: identity >= 30%, query coverage >= 20%, e-value <= 0.001
Edge weight:     thin  = qcov 20-40%
                 thick = qcov >= 40%

Outputs (Cytoscape import files):
    plots_dir/figure3_2-panelB/node.tsv
    plots_dir/figure3_2-panelB/edge.tsv

Checkpoint: skips if both files already exist.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

EVALUE_MAX = 1e-3
PIDENT_MIN = 30.0   # %
QCOV_MIN   = 0.20   # minimum bidirectional coverage to retain edge


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_fasta_single(path: Path) -> tuple[str, str]:
    """Return (seq_id, sequence) for the first record in a FASTA file."""
    seq_id, parts = "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if not seq_id:
                    seq_id = line[1:].split()[0]
            else:
                parts.append(line)
    return seq_id, "".join(parts)


def _sanitize(s: str) -> str:
    return s.replace(".", "_DOT_").replace("|", "_PIPE_").replace(" ", "_SPC_")


def _restore(s: str) -> str:
    return s.replace("_DOT_", ".").replace("_PIPE_", "|").replace("_SPC_", " ")


# ---------------------------------------------------------------------------
# Node loaders
# ---------------------------------------------------------------------------

def _load_kloci_nodes(acetylases_kloci_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(acetylases_kloci_tsv, sep="\t")
    df["locus"] = df["protein_id"].str.extract(r"^(KL\d+)")
    rows = []
    for _, row in df.iterrows():
        seq = str(row["sequence"]).replace(" ", "").strip()
        if not seq or seq == "nan":
            print(f"  [warn] {row['protein_id']}: empty sequence — skipping")
            continue
        evidence = [
            name for name, col in [
                ("experimental", "source_experimental"),
                ("blastp",       "source_blastp"),
                ("keyword",      "source_keyword"),
                ("foldseek",     "source_foldseek"),
            ] if row.get(col)
        ]
        rows.append({
            "proteinID":   row["protein_id"],
            "type":        "KLOCI",
            "specificity": row["locus"],
            "label":       row["locus"],
            "evidence":    ",".join(evidence),
            "seq":         seq,
        })
    return pd.DataFrame(rows)


def _load_experimental_nodes(acetylases_literature_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(acetylases_literature_tsv, sep="\t")
    rows = []
    for _, row in df.iterrows():
        ktype = str(row["ktype"]).strip()
        rest  = ktype[1:].strip()
        locus = f"KL{rest}" if ktype.startswith("K") and rest.isdigit() else ktype
        gene  = str(row.get("protein_name", "")).strip()
        seq   = str(row["sequence"]).replace(" ", "").strip()
        rows.append({
            "proteinID":   row["proteinid"],
            "type":        "EXPERIMENTAL",
            "specificity": locus,
            "label":       gene if gene and gene != "nan" else row["proteinid"],
            "evidence":    "experimental",
            "seq":         seq,
        })
    return pd.DataFrame(rows)


def _load_gwas_best_nodes(gwas_best_at_fastas: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for locus, fasta_path in gwas_best_at_fastas.items():
        seq_id, seq = _read_fasta_single(fasta_path)
        rows.append({
            "proteinID":   seq_id,
            "type":        "GWAS_BEST",
            "specificity": locus,
            "label":       locus,
            "evidence":    "gwas_best",
            "seq":         seq,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Node assembly
# ---------------------------------------------------------------------------

def _build_nodes(
    acetylases_kloci_tsv: Path,
    acetylases_literature_tsv: Path,
    gwas_best_at_fastas: dict[str, Path],
) -> pd.DataFrame:
    kloci = _load_kloci_nodes(acetylases_kloci_tsv)
    exp   = _load_experimental_nodes(acetylases_literature_tsv)
    gwas  = _load_gwas_best_nodes(gwas_best_at_fastas)

    nodes = pd.concat([kloci, exp, gwas], ignore_index=True)
    nodes.index = nodes.index + 1
    nodes["node_id"]   = nodes.index.astype(str) + "_" + nodes["label"]
    nodes["clean_id"]  = nodes["node_id"].str.replace(r"^\d+_", "", regex=True)
    nodes["full_label"] = nodes.apply(
        lambda r: f"{r['proteinID']}_{r['label']}" if r["type"] == "KLOCI" else r["proteinID"],
        axis=1,
    )
    return nodes[["node_id", "clean_id", "full_label", "proteinID",
                  "type", "specificity", "label", "evidence", "seq"]].copy()


# ---------------------------------------------------------------------------
# Edge generation (all-vs-all BLASTP)
# ---------------------------------------------------------------------------

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
         "-query", str(fasta_path),
         "-db",    str(fasta_path),
         "-out",   str(blast_out),
         "-outfmt", "6 qseqid sseqid evalue pident bitscore nident length "
                    "mismatch gapopen gaps qlen qstart qend slen sstart send"],
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
        (df["qcov"]   >= QCOV_MIN)   &
        (df["scov"]   >= QCOV_MIN)   &
        (df["pident"] >= PIDENT_MIN) &
        (df["evalue"] <= EVALUE_MAX)
    ].copy()

    df["coverage-sharing"] = pd.cut(
        df["qcov"],
        bins=[QCOV_MIN, 0.40, 1.001],
        labels=["thin", "thick"],
        right=False, include_lowest=True,
    )
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
               "coverage-sharing", "seq-identity", "interaction"]]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_figure3_2_panelB(
    acetylases_kloci_tsv: Path,
    acetylases_literature_tsv: Path,
    gwas_best_at_fastas: dict[str, Path],
    plots_dir: Path,
    style=None,
) -> None:
    """
    Build acetyltransferase similarity network and write Cytoscape node/edge TSVs.

    Args:
        acetylases_kloci_tsv:      cps_acetylases/acetylases_kloci.tsv (69 K-locus candidates)
        acetylases_literature_tsv: cps_acetylases/acetylases_literature_active.tsv (3 experimental anchors)
        gwas_best_at_fastas:       {locus: Path} for the 2 GWAS best AT predictions
                                   (KL30/PC0675 and KL111/PC0915, FoldSeek-identified)
        plots_dir:                 output directory; files go into plots_dir/figure3_2-panelB/
        style:                     cfg.style (optional, unused)
    """
    out_dir  = plots_dir / "figure3_2-panelB"
    node_out = out_dir / "node.tsv"
    edge_out = out_dir / "edge.tsv"

    if (node_out.exists() and node_out.stat().st_size > 0 and
            edge_out.exists() and edge_out.stat().st_size > 0):
        print("  figure3_2-panelB/node.tsv + edge.tsv exist — skipping (delete to recompute)")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    nodes   = _build_nodes(acetylases_kloci_tsv, acetylases_literature_tsv, gwas_best_at_fastas)
    n_kloci = (nodes["type"] == "KLOCI").sum()
    n_exp   = (nodes["type"] == "EXPERIMENTAL").sum()
    n_gwas  = (nodes["type"] == "GWAS_BEST").sum()
    print(f"  {len(nodes)} nodes: {n_kloci} K-locus candidates, {n_exp} experimental, {n_gwas} GWAS best")

    with tempfile.TemporaryDirectory() as _tmp:
        edges = _generate_edges(nodes, Path(_tmp))

    n_self     = (edges["source"] == edges["target"]).sum()
    n_non_self = len(edges) - n_self
    print(f"  {len(edges)} edges after filtering ({n_non_self} non-self, {n_self} self-hits)")

    nodes.to_csv(node_out, sep="\t", index=False)
    edges.to_csv(edge_out, sep="\t", index=False)
    print(f"  → figure3_2-panelB/node.tsv ({len(nodes)} rows)")
    print(f"  → figure3_2-panelB/edge.tsv ({len(edges)} rows)")
