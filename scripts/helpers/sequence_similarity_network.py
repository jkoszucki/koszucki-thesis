"""
Shared sequence similarity network builder — importable by any chapter script.

Runs all-vs-all BLASTP on a set of protein sequences, filters hits by
coverage / identity / e-value thresholds, assigns coverage interval labels,
and writes Cytoscape-ready edge.tsv and node.tsv files.

Node metadata (category, locus, label columns, …) is passed in by the caller,
keeping this class generic across all use cases.

Usage (in any figures/chapterX/main.py or lib/)
------------------------------------------------
    sys.path.insert(0, str(Path(__file__).resolve().parents[N] / "helpers"))
    from sequence_similarity_network import SequenceSimilarityNetwork

    net = SequenceSimilarityNetwork(
        min_pident          = 50.0,          # % identity threshold
        min_cov             = 0.30,          # query AND subject coverage (0–1)
        max_evalue          = 1e-3,          # e-value threshold
        coverage_intervals  = [30, 50, 80],  # bin boundaries in %; auto-labels low/medium/high
        include_singletons  = True,          # include unconnected nodes as self-loops
    )

    nodes, edges = net.build(
        sequences  = sequences_fasta,        # Path to FASTA  OR  dict {id: seq}
        output_dir = plots_dir,              # writes edge.tsv + node.tsv here
        tmp_dir    = plots_dir / "blast_tmp",
        node_attrs = my_metadata_df,         # optional DataFrame with extra node columns
                                             # must have a column matching the sequence IDs
                                             # (default join key: "name")
        node_id_col = "name",               # column in node_attrs to join on
    )

Coverage intervals
------------------
coverage_intervals defines the bin boundary values in percent, e.g.:
    [30, 50, 80]  →  3 bins: low (30–50%), medium (50–80%), high (80–100%)
    [50, 80]      →  2 bins: low (50–80%), high (80–100%)
    [80, 90]      →  2 bins: low (80–90%), high (90–100%)

The last bin always extends to 100%. Labels: 2 bins → low/high,
3 bins → low/medium/high, 4+ bins → bin1/bin2/…/binN.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_INTERVALS = [30, 50, 80]   # bin boundaries in percent
_BIN_LABEL_PRESETS = {
    1: ["all"],
    2: ["low", "high"],
    3: ["low", "medium", "high"],
}


# ---------------------------------------------------------------------------
# SequenceSimilarityNetwork
# ---------------------------------------------------------------------------

class SequenceSimilarityNetwork:
    """
    All-vs-all BLASTP similarity network with Cytoscape export.

    Parameters
    ----------
    min_pident : float
        Minimum % identity to retain an edge. Default: 50.0.
    min_cov : float
        Minimum query AND subject coverage to retain an edge (0–1 fraction).
        Default: 0.30.
    max_evalue : float
        Maximum e-value to retain an edge. Default: 1e-3.
    coverage_intervals : list[int]
        Bin boundary values in percent (ascending). The final bin always
        extends to 100%. Labels are auto-assigned: 2 bins → low/high,
        3 bins → low/medium/high, 4+ → bin1/bin2/…
        Default: [30, 50, 80].
    include_singletons : bool
        If True, nodes with no edges are added as self-loop rows so that
        Cytoscape displays them as isolated nodes. Default: True.
    """

    def __init__(
        self,
        min_pident:         float      = 50.0,
        min_cov:            float      = 0.30,
        max_evalue:         float      = 1e-3,
        coverage_intervals: list[int]  = None,
        include_singletons: bool       = True,
    ) -> None:
        self.min_pident         = min_pident
        self.min_cov            = min_cov
        self.max_evalue         = max_evalue
        self.coverage_intervals = list(coverage_intervals) if coverage_intervals else _DEFAULT_INTERVALS
        self.include_singletons = include_singletons

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(
        self,
        sequences:   dict[str, str] | Path,
        output_dir:  Path,
        tmp_dir:     Path,
        node_attrs:  pd.DataFrame | None = None,
        node_id_col: str                 = "name",
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Build the network and write edge.tsv + node.tsv to output_dir.

        Args:
            sequences:   FASTA path or dict {seq_id: sequence}.
            output_dir:  destination for edge.tsv and node.tsv.
            tmp_dir:     scratch directory for BLAST intermediates.
            node_attrs:  optional DataFrame with extra per-node columns.
                         Merged onto the node table by node_id_col.
            node_id_col: column in node_attrs that matches the sequence IDs.

        Returns:
            (nodes_df, edges_df) — full DataFrames (seq column included).
        """
        output_dir = Path(output_dir)
        tmp_dir    = Path(tmp_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        seqs  = self._load_sequences(sequences)
        nodes = self._make_nodes(seqs, node_attrs, node_id_col)
        raw   = self._run_blast(nodes, tmp_dir)
        edges = self._build_edges(raw, nodes)

        # Write outputs — drop seq from node table
        nodes_out = nodes.drop(columns=["seq", "blast_id"], errors="ignore")
        nodes_out.to_csv(output_dir / "node.tsv", sep="\t", index=False)
        edges.to_csv(output_dir / "edge.tsv", sep="\t", index=False)

        n_edges = len(edges[edges["source"] != edges["target"]])
        print(f"  {len(nodes)} nodes, {n_edges} edges → {output_dir / 'edge.tsv'}")
        return nodes, edges

    # ------------------------------------------------------------------
    # Sequence loading
    # ------------------------------------------------------------------

    @staticmethod
    def _load_sequences(sequences: dict[str, str] | Path) -> dict[str, str]:
        if isinstance(sequences, (str, Path)):
            from Bio import SeqIO
            return {rec.id: str(rec.seq) for rec in SeqIO.parse(sequences, "fasta")}
        return dict(sequences)

    # ------------------------------------------------------------------
    # Node assembly
    # ------------------------------------------------------------------

    @staticmethod
    def _make_nodes(
        seqs: dict[str, str],
        node_attrs: pd.DataFrame | None,
        node_id_col: str,
    ) -> pd.DataFrame:
        nodes = pd.DataFrame([
            {"name": sid, "seq": seq, "seq_len": len(seq)}
            for sid, seq in seqs.items()
        ])

        if node_attrs is not None and not node_attrs.empty:
            nodes = nodes.merge(
                node_attrs, left_on="name", right_on=node_id_col, how="left"
            )

        # Sanitize IDs for BLAST (dots truncate query IDs)
        nodes["blast_id"] = nodes["name"].str.replace(".", "_DOT_", regex=False)
        return nodes

    # ------------------------------------------------------------------
    # BLAST
    # ------------------------------------------------------------------

    @staticmethod
    def _run_blast(nodes: pd.DataFrame, tmp_dir: Path) -> pd.DataFrame:
        fasta_path = tmp_dir / "sequences.fasta"
        blast_out  = tmp_dir / "blastp_raw.tsv"

        with open(fasta_path, "w") as f:
            for _, row in nodes.iterrows():
                f.write(f">{row['blast_id']}\n{row['seq']}\n")

        subprocess.run(
            ["makeblastdb", "-in", str(fasta_path), "-dbtype", "prot"],
            check=True, capture_output=True,
        )
        subprocess.run(
            ["blastp",
             "-query", str(fasta_path), "-db", str(fasta_path),
             "-out",   str(blast_out),
             "-outfmt", "6 qseqid sseqid evalue pident bitscore "
                        "qlen qstart qend slen sstart send"],
            check=True, capture_output=True,
        )

        cols = ["source", "target", "evalue", "pident", "bitscore",
                "qlen", "qstart", "qend", "slen", "sstart", "send"]
        df = pd.read_csv(blast_out, sep="\t", header=None, names=cols)

        blast_to_name = dict(zip(nodes["blast_id"], nodes["name"]))
        for col in ("source", "target"):
            df[col] = (df[col]
                       .map(blast_to_name)
                       .fillna(df[col].str.replace("_DOT_", ".", regex=False)))
        return df

    # ------------------------------------------------------------------
    # Edge filtering and annotation
    # ------------------------------------------------------------------

    def _build_edges(self, df: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()

        # Coverage fractions
        df["qcov"] = np.round((df["qend"] - df["qstart"] + 1) / df["qlen"], 2)
        df["scov"] = np.round((df["send"] - df["sstart"] + 1) / df["slen"], 2)

        # Deduplicate symmetric hits
        df["_pair"] = df.apply(
            lambda r: "-".join(sorted([r["source"], r["target"]])), axis=1
        )
        df = df.drop_duplicates(subset="_pair").drop(columns=["_pair"])

        # Filter
        keep = (
            (df["qcov"]   >= self.min_cov)
            & (df["scov"]  >= self.min_cov)
            & (df["pident"] >= self.min_pident)
            & (df["evalue"] <= self.max_evalue)
            & (df["source"] != df["target"])
        )
        df = df.loc[keep].copy()

        # Coverage interval labels
        df["coverage_interval"] = self._label_coverage(df)

        # Singletons
        if self.include_singletons:
            connected  = set(df["source"]) | set(df["target"])
            singletons = [n for n in nodes["name"] if n not in connected]
            if singletons:
                sing = pd.DataFrame({"source": singletons, "target": singletons})
                df   = pd.concat([df, sing], ignore_index=True)

        df["interaction"] = "A"
        return df

    def _label_coverage(self, df: pd.DataFrame) -> pd.Series:
        """Assign coverage interval labels using self.coverage_intervals."""
        boundaries = sorted(self.coverage_intervals)
        if boundaries[-1] < 100:
            boundaries = boundaries + [100]

        n_bins = len(boundaries) - 1
        prefix_list = _BIN_LABEL_PRESETS.get(
            n_bins,
            [f"bin{i + 1}" for i in range(n_bins)]
        )
        labels = [
            f"{prefix_list[i]} ({boundaries[i]}–{boundaries[i + 1]}%)"
            for i in range(n_bins)
        ]

        cov = df[["qcov", "scov"]].min(axis=1) * 100
        return pd.cut(
            cov,
            bins=boundaries,
            labels=labels,
            right=True,
            include_lowest=True,
        )
