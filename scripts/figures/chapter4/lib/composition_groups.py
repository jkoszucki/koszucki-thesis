import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_set(val):
    if pd.isna(val):
        return frozenset()
    s = str(val).strip()
    tokens = re.split(r"[;,|]\s*|\s+", s)
    toks = [t.strip() for t in tokens if t and t.strip()]
    return frozenset(toks)


def compute_jaccard(s1: frozenset, s2: frozenset) -> float:
    if not s1 and not s2:
        return 1.0
    union = s1 | s2
    if len(union) == 0:
        return 0.0
    inter = s1 & s2
    return len(inter) / len(union)


def plot_composition_groups(input_csv: Path, outdir: Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    if "monos_unique_all" not in df.columns:
        raise KeyError("'monos_unique_all' column missing in input")

    df["composition_set"] = df["monos_unique_all"].apply(parse_set)

    comp_counts = df["composition_set"].value_counts().sort_values(ascending=False)

    total = comp_counts.sum()
    comp_df = pd.DataFrame({
        "composition": [";".join(sorted(x)) if x else "" for x in comp_counts.index],
        "count": comp_counts.values,
        "frequency": comp_counts.values / total,
    })
    comp_df.to_csv(outdir / "composition_freq.csv", index=False)

    unique_sets = list(comp_counts.index)
    labels = [";".join(sorted(s)) if s else "" for s in unique_sets]

    records = []
    for i, s1 in enumerate(unique_sets):
        for j, s2 in enumerate(unique_sets):
            rec = {
                "composition_1": labels[i],
                "composition_2": labels[j],
                "jaccard": compute_jaccard(s1, s2),
            }
            records.append(rec)

    sim_long_df = pd.DataFrame(records)
    sim_long_df = sim_long_df.query("jaccard>=0.8")
    sim_long_df.to_csv(outdir / "composition_sim.csv", index=False)

    plt.figure(figsize=(5, 4))
    sim_long_df["jaccard"].hist(bins=30, edgecolor="black")
    plt.xlabel("Jaccard similarity", fontsize=12, fontweight="bold")
    plt.ylabel("Count", fontsize=12, fontweight="bold")
    plt.title("Distribution of Jaccard similarities", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(outdir / "composition_jaccard_hist.png", dpi=300)
    plt.close()

    print("Wrote:")
    print("  ", outdir / "composition_freq.csv")
    print("  ", outdir / "composition_sim.csv")
