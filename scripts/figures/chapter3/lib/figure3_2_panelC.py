"""Styled PNG/PDF table for Figure 3.2C — CPS structural pair comparison."""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# Manually curated CPS pair comparison table (14 rows).
#
# Last validated against: Chapter 3 Results, thesis-manuscript Doc, 2026-08-25.
# (The Results section is unchanged from 2026-08-23; only the Discussion was revised.)
# Every cell below traces to a sentence in that section; the fourteen pairs are the
# ones the text says are "compared in detail" (14 of the 35 network edges, 40%):
#   - eight near-identical pairs   -> Figure 3.3
#   - three core-shared pairs      -> Figure 3.4
#   - three branch-shared pairs    -> Figure 3.5
# Do not edit a row without editing the corresponding manuscript sentence.
_CPS_PAIRS_DATA = [
    # --- Near-identical (both core and branch path Jaccard > 0.6), Figure 3.3 ---
    {"cps_pair": "K2/K2",     "serotype": "same",       "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "acetylation",                     "modification_location": None,                  "structures_comparison": "Fig 3.3A"},
    {"cps_pair": "K22/K37",   "serotype": "distinct",   "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "acetylation",                     "modification_location": None,                  "structures_comparison": "Fig 3.3D"},
    {"cps_pair": "K30/K33",   "serotype": "distinct",   "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "acetyl occupancy",                "modification_location": "33% vs 100%",         "structures_comparison": "Fig 3.3C"},
    {"cps_pair": "K30/K69",   "serotype": "distinct",   "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "pyruvylation",                    "modification_location": "ketal bond position", "structures_comparison": "Fig 3.3C"},
    {"cps_pair": "K33/K69",   "serotype": "distinct",   "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "acetyl occupancy, pyruvylation",  "modification_location": "33% vs 100%; ketal bond", "structures_comparison": "Fig 3.3C"},
    {"cps_pair": "K82/SK1",   "serotype": "unassigned", "similarity_location": "near-identical", "dissimilarity_location": None,     "structural_change": None,                          "modification_change": "acetylation, glutamylation",      "modification_location": None,                  "structures_comparison": "Fig 3.3B"},
    {"cps_pair": "K8/K82",    "serotype": "distinct",   "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition (Gal)",      "modification_change": "glutamylation",                   "modification_location": None,                  "structures_comparison": "Fig 3.3B"},
    {"cps_pair": "SK1/K8",    "serotype": "unassigned", "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition (Gal)",      "modification_change": "acetylation",                     "modification_location": None,                  "structures_comparison": "Fig 3.3B"},
    # --- Core-shared (core > 0.6 only), Figure 3.4 ---
    {"cps_pair": "K2/K13",    "serotype": "distinct",   "similarity_location": "core",           "dissimilarity_location": "branch", "structural_change": "residue addition (Gal)",      "modification_change": "acetylation, pyruvylation",       "modification_location": "20% occupancy",       "structures_comparison": "Fig 3.4A"},
    {"cps_pair": "K74/K80",   "serotype": "distinct",   "similarity_location": "core",           "dissimilarity_location": "branch", "structural_change": "residue swap (Gal to Rha)",   "modification_change": "pyruvylation",                    "modification_location": "distinct bond",       "structures_comparison": "Fig 3.4B"},
    {"cps_pair": "K16/K58",   "serotype": "distinct",   "similarity_location": "core",           "dissimilarity_location": "branch", "structural_change": "branch relocation, anomeric", "modification_change": "acetylation, pyruvylation",       "modification_location": "core",                "structures_comparison": "Fig 3.4C"},
    # --- Branch-shared (branch > 0.6 only), Figure 3.5 ---
    {"cps_pair": "K8/K8",     "serotype": "same",       "similarity_location": "branch",         "dissimilarity_location": "core",   "structural_change": "residue displacement",        "modification_change": "pyruvylation",                    "modification_location": "core",                "structures_comparison": "Fig 3.5A"},
    {"cps_pair": "K21a/K21b", "serotype": "same",       "similarity_location": "branch",         "dissimilarity_location": "core",   "structural_change": "residue swap (Man to Rha)",   "modification_change": "acetylation",                     "modification_location": None,                  "structures_comparison": "Fig 3.5B"},
    {"cps_pair": "K27a/K27b", "serotype": "same",       "similarity_location": "branch",         "dissimilarity_location": "core",   "structural_change": "anomeric configuration",      "modification_change": None,                              "modification_location": None,                  "structures_comparison": "Fig 3.5C"},
]

# The five pairs whose *monosaccharide composition* differs (a residue is added or
# swapped). A narrower set than `structural_change`, which also covers relocation,
# displacement and anomeric changes that leave composition untouched.
_COMPOSITION_DIFF_PAIRS = {"K8/K82", "SK1/K8", "K2/K13", "K74/K80", "K21a/K21b"}

# Counts asserted by the Chapter 3 Summary section. Three are quoted verbatim from
# the manuscript; the acetylation count is NOT — see the note below.
_EXPECTED_COUNTS = {
    "differ_in_any_modification": 13,   # manuscript: "thirteen differ in at least one modification (92.9%)"
    "differ_in_modification_only": 6,   # manuscript: "six differ in modifications only (42.9%)"
    "differ_in_composition": 5,         # manuscript: "differs in only five pairs and never alone (35.7%)"
    "differ_in_acetylation": 9,         # manuscript says TEN (71.4%) -- see note
}

# NOTE ON THE ACETYLATION COUNT
# The manuscript Summary states acetylation "was most often identified to be distinct
# between similar pairs of structures (71.4%; 10 out of 14 pairs analysed in detail)".
# Only nine of the fourteen pairs can differ in acetylation, given the manuscript's own
# description of the K8/K82/SK1 component: K8 "carries ... no modifications at all",
# and K82 and SK1 "differ in the presence of acetylation and of N-L-glutamate" -- so
# acetylation sits on exactly one of K82/SK1, and therefore separates only one of the
# two pairs K8/K82 and SK1/K8 from K8, not both. Counting both would also require SK1/K8
# and K8/K82 each to differ in acetylation, which contradicts K8 having no modifications.
# The nine are: K2/K2, K22/K37, K30/K33, K33/K69, K82/SK1, SK1/K8, K2/K13, K16/K58,
# K21a/K21b. The five that do not: K30/K69 (pyruvate bond), K8/K82 (glutamate),
# K74/K80 (pyruvylation), K8/K8 (pyruvylation), K27a/K27b (anomeric only).
# Flagged to the author; if the manuscript is revised to 9/14 (64.3%) this table already
# agrees. If instead a pair is reclassified, update `_EXPECTED_COUNTS` with the reason.


def _verify_manuscript_statistics(rows):
    """Cross-check the table against the counts quoted in the Chapter 3 Summary.

    Raises if a future edit to `_CPS_PAIRS_DATA` silently breaks agreement with the
    manuscript -- the failure mode this table has already been through once.
    """
    def _has(row, key):
        val = row.get(key)
        return val is not None and str(val).strip() not in ("", "-")

    actual = {
        "differ_in_any_modification": sum(_has(r, "modification_change") for r in rows),
        "differ_in_modification_only": sum(
            _has(r, "modification_change") and not _has(r, "structural_change") for r in rows
        ),
        "differ_in_composition": sum(r["cps_pair"] in _COMPOSITION_DIFF_PAIRS for r in rows),
        "differ_in_acetylation": sum(
            "acetyl" in str(r.get("modification_change") or "").lower() for r in rows
        ),
    }
    mismatches = {k: (v, actual[k]) for k, v in _EXPECTED_COUNTS.items() if actual[k] != v}
    if mismatches:
        detail = "; ".join(
            f"{k}: expected {exp}, table gives {got}" for k, (exp, got) in mismatches.items()
        )
        raise ValueError(
            f"Figure 3.2C table no longer matches the Chapter 3 Summary statistics ({detail}). "
            "Reconcile the table with the manuscript before regenerating the figure."
        )
    if len(rows) != 14:
        raise ValueError(f"Figure 3.2C must have 14 rows, got {len(rows)}.")


# Column display names for the 6-column format
COL_LABELS = {
    "cps_pair": "CPS pair",
    "similarity_location": "Relationship",
    "serotype": "Serotype",
    "structural_difference": "Structural\ndifference",
    "modification_difference": "Modification\ndifference",
    "structures_comparison": "Figure",
}

# Colour scheme by similarity_location
GROUP_COLORS = {
    "near-identical": "#dce9f5",
    "branch": "#d5e8d4",
    "core": "#fff2cc",
}

# Serotype colours
SEROTYPE_COLORS = {
    "same": "#e8f4e8",
    "distinct": "#fce8e8",
    "unassigned": "#eeeeee",   # SK1 carries no serological assignment
}

# Relationship display labels
RELATIONSHIP_LABELS = {
    "near-identical": "Near-identical",
    "branch": "Branch-shared",
    "core": "Core-shared",
}


def _merge_with_location(change, location):
    """Merge a change description with its location in parentheses."""
    if pd.isna(change) or str(change).strip() in ("", "—"):
        return "—"
    change = str(change).strip()
    if pd.isna(location) or str(location).strip() in ("", "—"):
        return change
    return f"{change}\n({location})"


def _wrap(text, width=18):
    """Wrap long strings at spaces."""
    if not isinstance(text, str) or len(text) <= width:
        return text
    words = text.split()
    lines, line = [], []
    for w in words:
        if sum(len(x) for x in line) + len(line) + len(w) > width:
            lines.append(" ".join(line))
            line = [w]
        else:
            line.append(w)
    if line:
        lines.append(" ".join(line))
    return "\n".join(lines)


def plot_figure3_2_panelC(output_path: Path, style=None) -> None:
    _verify_manuscript_statistics(_CPS_PAIRS_DATA)
    raw = pd.DataFrame(_CPS_PAIRS_DATA)

    # Build merged 6-column dataframe
    df = pd.DataFrame()
    df["cps_pair"] = raw["cps_pair"]
    df["similarity_location"] = raw["similarity_location"]
    df["serotype"] = raw["serotype"]
    df["structural_difference"] = [
        _merge_with_location(row["structural_change"], row["dissimilarity_location"])
        for _, row in raw.iterrows()
    ]
    df["modification_difference"] = [
        _merge_with_location(row["modification_change"], row["modification_location"])
        for _, row in raw.iterrows()
    ]
    df["structures_comparison"] = raw["structures_comparison"]

    cols = list(COL_LABELS.keys())
    headers = [COL_LABELS[c] for c in cols]

    # Column widths (relative): pair, relationship, serotype, structural, modification, figure
    col_widths = [0.95, 1.2, 0.9, 2.2, 3.2, 0.75]

    # Per-column char wrap widths — wide columns wrap late to avoid unnecessary line breaks
    col_wrap_widths = [12, 16, 12, 28, 42, 10]

    # Build display data: replace relationship codes with labels, wrap text per column
    display_data = []
    for _, row in df.iterrows():
        cell_row = []
        for c, ww in zip(cols, col_wrap_widths):
            val = row[c]
            if c == "similarity_location":
                val = RELATIONSHIP_LABELS.get(str(val), str(val))
            cell_row.append(_wrap(str(val), ww) if isinstance(val, str) else str(val))
        display_data.append(cell_row)

    n_rows = len(display_data)

    fig_width = (sum(col_widths) * 1.25 + 0.4) * 1.5
    row_height = 0.55
    header_height = 0.60
    fig_height = header_height + n_rows * row_height + 0.3

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis("off")

    # Normalise col widths to 0–1 range for positioning
    total_w = sum(col_widths)
    cum_x = [0.0]
    for w in col_widths[:-1]:
        cum_x.append(cum_x[-1] + w / total_w)
    norm_widths = [w / total_w for w in col_widths]

    margin_l = 0.01
    margin_r = 0.01
    usable = 1.0 - margin_l - margin_r
    x_starts = [margin_l + cx * usable for cx in cum_x]
    x_widths = [nw * usable for nw in norm_widths]

    total_h = header_height + n_rows * row_height
    y_top = 1.0
    y_header_bottom = y_top - header_height / total_h

    def y_row_top(i):
        return y_top - (header_height + i * row_height) / total_h

    def y_row_bottom(i):
        return y_top - (header_height + (i + 1) * row_height) / total_h

    # Draw header
    for lbl, x0, xw in zip(headers, x_starts, x_widths):
        rect = mpatches.FancyBboxPatch(
            (x0, y_header_bottom), xw, (y_top - y_header_bottom),
            boxstyle="square,pad=0", linewidth=0.5,
            edgecolor="#555555", facecolor="#3a3a3a",
            transform=ax.transAxes, clip_on=False,
        )
        ax.add_patch(rect)
        ax.text(
            x0 + xw / 2, (y_top + y_header_bottom) / 2,
            lbl, ha="center", va="center",
            fontsize=9.5, fontweight="bold", color="white",
            transform=ax.transAxes, clip_on=False,
            linespacing=1.3,
        )

    # Draw data rows
    for i, row in enumerate(display_data):
        rel = df["similarity_location"].iloc[i]
        row_bg = GROUP_COLORS.get(rel, "#ffffff")
        row_top = y_row_top(i)
        row_bot = y_row_bottom(i)
        row_h = row_top - row_bot

        for j, (val, x0, xw) in enumerate(zip(row, x_starts, x_widths)):
            if cols[j] == "serotype":
                bg = SEROTYPE_COLORS.get(str(df["serotype"].iloc[i]), row_bg)
            else:
                bg = row_bg

            rect = mpatches.FancyBboxPatch(
                (x0, row_bot), xw, row_h,
                boxstyle="square,pad=0", linewidth=0.4,
                edgecolor="#aaaaaa", facecolor=bg,
                transform=ax.transAxes, clip_on=False,
            )
            ax.add_patch(rect)
            ax.text(
                x0 + xw / 2, (row_top + row_bot) / 2,
                val, ha="center", va="center",
                fontsize=9, color="#222222", fontweight="bold",
                transform=ax.transAxes, clip_on=False,
                linespacing=1.3,
            )

    # Outer border
    outer = mpatches.FancyBboxPatch(
        (margin_l, y_row_bottom(n_rows - 1)),
        usable, y_top - y_row_bottom(n_rows - 1),
        boxstyle="square,pad=0", linewidth=1.0,
        edgecolor="#333333", facecolor="none",
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(outer)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=GROUP_COLORS["near-identical"], edgecolor="#555", label="Near-identical"),
        mpatches.Patch(facecolor=GROUP_COLORS["branch"], edgecolor="#555", label="Branch-shared"),
        mpatches.Patch(facecolor=GROUP_COLORS["core"], edgecolor="#555", label="Core-shared"),
        mpatches.Patch(facecolor=SEROTYPE_COLORS["same"], edgecolor="#555", label="Same serotype"),
        mpatches.Patch(facecolor=SEROTYPE_COLORS["distinct"], edgecolor="#555", label="Distinct serotype"),
        mpatches.Patch(facecolor=SEROTYPE_COLORS["unassigned"], edgecolor="#555", label="No serotype assigned"),
    ]
    ax.legend(
        handles=legend_elements, loc="lower center",
        bbox_to_anchor=(0.5, -0.07), ncol=6,
        fontsize=8.5, frameon=False,
    )

    dpi = style.dpi if style else 200
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {output_path}")
