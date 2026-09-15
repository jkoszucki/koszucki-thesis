"""
S11 Figure, Panel B (support) — representative-level structural similarity network.

The all-vs-all network (`figureS11_panelB_tmalign.py`) is dense: at TM >= 0.75 its 73 nodes
carry 345 edges, and near-identical proteins dominate the picture. This module collapses
that redundancy — cluster at TM >= 0.75, keep one representative per cluster, and draw the
relationships between representatives — so each node stands for a distinct structure rather
than a family of near-copies.

Pipeline:
  1. Cluster the 73 proteins by TM-score, average linkage (UPGMA), cut at `cluster_tm`
     (0.75). The same clustering machinery as `processing/acetyl-proc/structure_order.py`,
     at a stricter cut: 0.5 groups by fold, 0.75 groups near-identical structures.
  2. One representative per cluster. The default is the medoid — the member with the
     highest mean TM-score to the rest of its cluster, i.e. the most typical structure.
     With `prefer_reference=True` a cluster containing an S1_Table reference protein
     (active_AT_*/gwas_AT_*) is represented by it instead, so the characterised enzymes
     stay visible in the figure rather than being absorbed into a K-locus candidate.
  3. Edges between representatives, each classified by `edge_class_tm` (0.6) into
     `high (>= 0.6)` / `low (< 0.6)`.

Note on density: the current 19 representatives give 171 possible pairs, of which 59 pass
0.6. Drawing them all is a complete graph and unreadable; `min_edge_tm` sets a floor on which
edges are written at all (default `edge_class_tm`, i.e. only the `high` class). Set it to
0.0 to emit every pair and filter inside Cytoscape instead.

Outputs:
    plots_dir/supplement/figureS11-panelB-representatives/node.tsv
    plots_dir/supplement/figureS11-panelB-representatives/edge.tsv
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from figureS11_panelB_tmalign import _label, _node_type, _specificity

CLUSTER_TM    = 0.75   # cut for grouping near-identical structures
EDGE_CLASS_TM = 0.60   # high/low split on edges between representatives


def _similarity_matrix(pairs: pd.DataFrame, ids: list[str]) -> np.ndarray:
    idx = {k: i for i, k in enumerate(ids)}
    n = len(ids)
    sim = np.ones((n, n))
    for row in pairs.itertuples(index=False):
        if row.id1 in idx and row.id2 in idx:
            sim[idx[row.id1], idx[row.id2]] = sim[idx[row.id2], idx[row.id1]] = row.tm_max
    return sim


def _pick_representative(members: list[int], sim: np.ndarray, ids: list[str],
                         prefer_reference: bool) -> int:
    if prefer_reference:
        refs = [i for i in members if _node_type(ids[i]) != "KLOCI"]
        if refs:
            # More than one reference in a cluster: take the most typical of them.
            return max(refs, key=lambda i: np.mean([sim[i, j] for j in members if j != i])) \
                if len(refs) > 1 else refs[0]
    if len(members) == 1:
        return members[0]
    return max(members, key=lambda i: np.mean([sim[i, j] for j in members if j != i]))


def plot_figureS11_panelB_representatives(
    tmscore_tsv: Path,
    structure_order_tsv: Path,
    plots_dir: Path,
    cluster_tm: float = CLUSTER_TM,
    edge_class_tm: float = EDGE_CLASS_TM,
    min_edge_tm: float | None = None,
    prefer_reference: bool = True,
) -> None:
    """
    Cluster at `cluster_tm`, pick representatives, and write their similarity network.

    Args:
        tmscore_tsv:         cps_acetylases/tmalign/acetylases_kloci_tmscore.tsv
        structure_order_tsv: defines the protein set (73 rows)
        plots_dir:           scripts/figures/chapter4/plots
        cluster_tm:          TM-score cut defining a cluster (default 0.75)
        edge_class_tm:       high/low classification of representative edges (default 0.6)
        min_edge_tm:         floor on which edges are written; None → `edge_class_tm`.
                             0.0 emits every pair (complete graph) for filtering elsewhere.
        prefer_reference:    represent a cluster by its S1_Table protein when it has one
    """
    min_edge_tm = edge_class_tm if min_edge_tm is None else min_edge_tm

    ids = pd.read_csv(structure_order_tsv, sep="\t").sort_values("order_index")["protein_id"].tolist()
    sim = _similarity_matrix(pd.read_csv(tmscore_tsv, sep="\t"), ids)

    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    z = linkage(squareform(dist, checks=False), method="average")
    labels = fcluster(z, t=1.0 - cluster_tm, criterion="distance")

    rows = []
    for cluster_id in sorted(set(labels)):
        members = [i for i in range(len(ids)) if labels[i] == cluster_id]
        rep = _pick_representative(members, sim, ids, prefer_reference)
        member_ids = [ids[i] for i in members]
        rows.append({
            "name":              ids[rep],
            "proteinID":         ids[rep],
            "type":              _node_type(ids[rep]),
            "specificity":       _specificity(ids[rep]),
            "label":             _label(ids[rep], _node_type(ids[rep])),
            "cluster_id":        int(cluster_id),
            "cluster_size":      len(members),
            "contains_reference": any(_node_type(m) != "KLOCI" for m in member_ids),
            "members":           ",".join(member_ids),
        })
    nodes = pd.DataFrame(rows).sort_values("cluster_size", ascending=False)

    pos = {k: i for i, k in enumerate(ids)}
    edges = []
    for a, b in itertools.combinations(nodes["name"], 2):
        tm = sim[pos[a], pos[b]]
        if tm < min_edge_tm:
            continue
        edges.append({
            "source": a, "target": b,
            "interaction": "structural_similarity",
            "tm_max": round(float(tm), 4),
            "tm_class": f"high (>= {edge_class_tm})" if tm >= edge_class_tm
                        else f"low (< {edge_class_tm})",
        })
    edges = pd.DataFrame(edges).sort_values("tm_max", ascending=False)

    out_dir = Path(plots_dir) / "supplement" / "figureS11-panelB-representatives"
    out_dir.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(out_dir / "node.tsv", sep="\t", index=False)
    edges.to_csv(out_dir / "edge.tsv", sep="\t", index=False)

    n_possible = len(nodes) * (len(nodes) - 1) // 2
    connected = set(edges["source"]) | set(edges["target"]) if len(edges) else set()
    print(f"  {len(nodes)} representatives from {len(ids)} proteins "
          f"(clusters at TM >= {cluster_tm}; sizes {sorted(nodes['cluster_size'], reverse=True)})")
    print(f"  {len(edges)}/{n_possible} possible edges written (floor TM >= {min_edge_tm}); "
          f"{edges['tm_class'].value_counts().to_dict() if len(edges) else {}}")
    print(f"  density {len(edges) / n_possible:.0%}; "
          f"{len(nodes) - len(connected)} representative(s) isolated")
    print(f"  → {out_dir}/node.tsv, edge.tsv")
