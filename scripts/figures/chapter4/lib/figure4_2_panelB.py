from pathlib import Path

import pandas as pd


def plot_figure4_2_panelB(
    similarity_csv: Path,
    ktypes_csv: Path,
    modifications_csv: Path,
    output_dir: Path,
    style=None,
) -> None:
    pyr_color = style.pyruvylation_color if style else "#3283bf"
    ace_color = style.acetylation_color if style else "#d45515"
    both_color = style.both_color if style else "#7A3FD4"
    gray_color = style.gray_color if style else "#bfbfbf"

    sim = pd.read_csv(similarity_csv)
    ktypes = pd.read_csv(ktypes_csv)
    mods = pd.read_csv(modifications_csv)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- EDGES ----
    core_thr = 0.6
    branch_thr = 0.6

    c = sim["path_jaccard_core"].fillna(0)
    b = sim["path_jaccard_branch"].fillna(0)

    mask_total  = (c > core_thr)   & (b > branch_thr)
    mask_core   = (c > core_thr)   & (b > 0) & ~mask_total
    mask_branch = (c > 0)          & (b > branch_thr) & ~mask_total

    edges = sim[mask_total | mask_core | mask_branch].copy()
    edges = edges.rename(columns={
        "structure_id_1": "source",
        "structure_id_2": "target",
    })

    c2 = edges["path_jaccard_core"].fillna(0)
    b2 = edges["path_jaccard_branch"].fillna(0)
    edges["interaction"] = "sim_core"
    edges.loc[(c2 > 0) & (b2 > branch_thr), "interaction"] = "sim_branch"
    edges.loc[(c2 > core_thr) & (b2 > branch_thr), "interaction"] = "sim_total"

    # Add self-loop rows for singletons so Cytoscape creates a node for each structure
    connected = set(edges["source"]) | set(edges["target"])
    all_ids = ktypes["structure_id"].tolist()
    singletons = [sid for sid in all_ids if sid not in connected]
    if singletons:
        singleton_rows = pd.DataFrame({
            "source": singletons,
            "target": singletons,
            "path_jaccard_core": float("nan"),
            "path_jaccard_branch": float("nan"),
            "path_jaccard_total": float("nan"),
            "interaction": "singleton",
        })
        edges = pd.concat([edges, singleton_rows], ignore_index=True)

    edge_cols = [
        "source", "target",
        "path_jaccard_core", "path_jaccard_branch", "path_jaccard_total",
        "interaction",
    ]
    edges[[c for c in edge_cols if c in edges.columns]].to_csv(
        output_dir / "edge.tsv", sep="\t", index=False
    )

    # ---- NODES ----
    has_ace = set(mods.loc[mods["modification"] == "acetylation", "structure_id"].dropna())
    has_pyr = set(mods.loc[mods["modification"] == "pyruvylation", "structure_id"].dropna())

    def mod_category(sid):
        a, p = sid in has_ace, sid in has_pyr
        if a and p:
            return "both"
        if a:
            return "acetylation"
        if p:
            return "pyruvylation"
        return "none"

    color_map = {
        "both": both_color,
        "acetylation": ace_color,
        "pyruvylation": pyr_color,
        "none": gray_color,
    }

    meta_cols = ["structure_id", "K_type", "Source"]
    for col in ["mw_struct_kda", "Year"]:
        if col in ktypes.columns:
            meta_cols.append(col)

    nodes = ktypes[meta_cols].copy()
    nodes["modification"] = nodes["structure_id"].apply(mod_category)
    nodes["color"] = nodes["modification"].map(color_map)
    nodes["label"] = nodes["structure_id"]
    nodes.to_csv(output_dir / "node.tsv", sep="\t", index=False)

    counts = edges["interaction"].value_counts().to_dict()
    print(
        f"Cytoscape export: {len(edges)} edges "
        f"(sim_total={counts.get('sim_total', 0)}, "
        f"sim_core={counts.get('sim_core', 0)}, "
        f"sim_branch={counts.get('sim_branch', 0)}, "
        f"singleton={counts.get('singleton', 0)}), "
        f"{len(nodes)} nodes → {output_dir}"
    )
