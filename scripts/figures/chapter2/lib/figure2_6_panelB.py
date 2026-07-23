"""
Chapter 2 — Figure 2.6 Panel B: Cytoscape network of distinct-specificity depolymerase pairs.

Reads BLASTP hits from figure6_panelA analysis outputs, filters to distinct-specificity
pairs, and writes node/edge tables for Cytoscape.

Node colour encodes source: magenta = prophage, purple = virulent.

Reads:
    active_enzymes_tsv   — analysis_dir/active_enzymes.tsv
    hits_tsv             — analysis_dir/figure6_panelA_table.tsv

Writes:
    plots_dir/figure2_6-panelB/node.tsv
    plots_dir/figure2_6-panelB/edge.tsv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_PROPHAGE_COLOR = "#e000e0"   # magenta
_VIRULENT_COLOR = "#6a0dad"   # purple


def _normalize_label(spec: str) -> str:
    """Normalize each /‑separated part to KL prefix for display."""
    parts = []
    for p in str(spec).split("/"):
        p = p.strip()
        if p.upper().startswith("KL") or p.upper().startswith("KN"):
            parts.append(p)
        elif p.upper().startswith("K"):
            parts.append("KL" + p[1:])
        else:
            parts.append(p)
    return "/".join(parts)


def plot_figure2_6_panelB(
    active_enzymes_tsv: Path,
    hits_tsv: Path,
    plots_dir: Path,
    style=None,
) -> None:
    """
    Write Cytoscape node/edge tables for distinct-specificity depolymerase pairs.

    Args:
        active_enzymes_tsv: analysis_dir/active_enzymes.tsv (proteinID, source, K_locus_specificity)
        hits_tsv:           analysis_dir/figure6_panelA_table.tsv (clean hit table)
        plots_dir:          output directory for plots
        style:              cfg.style (optional, unused — colours are hardcoded)
    """
    enzymes = pd.read_csv(active_enzymes_tsv, sep="\t")
    source_map = dict(zip(enzymes["proteinID"], enzymes["source"]))
    spec_map   = dict(zip(enzymes["proteinID"], enzymes["K_locus_specificity"]))

    hits = pd.read_csv(hits_tsv, sep="\t")

    # Unique distinct-specificity pairs (use "all" category to avoid duplicates)
    pairs = hits[hits["category"] == "distinct-specificity (all)"].copy()
    pairs = pairs.drop_duplicates(subset=["protein1", "protein2"])
    print(f"  [figure2_6-panelB] Distinct-specificity pairs: {len(pairs)}")

    # Collect all node IDs from these pairs
    node_ids = pd.unique(pairs[["protein1", "protein2"]].values.ravel())

    nodes = pd.DataFrame({
        "node_id":     node_ids,
        "label":       [_normalize_label(spec_map.get(n, n)) for n in node_ids],
        "source":      [source_map.get(n, "unknown") for n in node_ids],
        "color":       [
            _PROPHAGE_COLOR if source_map.get(n) == "prophage" else _VIRULENT_COLOR
            for n in node_ids
        ],
        "specificity": [spec_map.get(n, "") for n in node_ids],
    })

    edges = pairs.rename(columns={
        "protein1":     "source_node",
        "protein2":     "target_node",
        "seq_similarity": "seq_similarity",
        "C-terminus":   "C_terminus",
    })[["source_node", "target_node", "seq_similarity", "C_terminus"]].copy()

    def _source_pair(q, s):
        sq, ss = source_map.get(q, "unknown"), source_map.get(s, "unknown")
        if sq == "prophage" and ss == "prophage":
            return "prophage-prophage"
        if sq == "virulent" and ss == "virulent":
            return "virulent-virulent"
        if {sq, ss} == {"prophage", "virulent"}:
            return "prophage-virulent"
        return "other"

    edges["source_pair"] = [_source_pair(q, s) for q, s in zip(edges["source_node"], edges["target_node"])]

    out_dir = Path(plots_dir) / "figure2_6-panelB"
    out_dir.mkdir(parents=True, exist_ok=True)

    nodes.to_csv(out_dir / "node.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "edge.tsv", sep="\t", index=False)

    print(f"  [figure2_6-panelB] {len(nodes)} nodes, {len(edges)} edges")
    print(f"  [figure2_6-panelB] → {out_dir}/node.tsv")
    print(f"  [figure2_6-panelB] → {out_dir}/edge.tsv")
