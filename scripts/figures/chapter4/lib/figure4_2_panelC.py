"""Styled PNG/PDF table for Figure 4.2C — CPS structural pair comparison."""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# Manually curated CPS pair comparison table (14 rows).
# Formerly stored in the chapter3_figure2 sheet of cps.xlsx (archive).
_CPS_PAIRS_DATA = [
    {"cps_pair": "K2/K2",    "serotype": "same",     "similarity_location": "near-identical", "dissimilarity_location": None,    "structural_change": None,              "modification_change": "acetylation",  "modification_location": "core",   "structures_comparison": "Fig3A"},
    {"cps_pair": "SK1/K8",   "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition", "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3B"},
    {"cps_pair": "K82/K8",   "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition", "modification_change": "deacetylation","modification_location": "core",   "structures_comparison": "Fig3B"},
    {"cps_pair": "K40/K2",   "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition", "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3C"},
    {"cps_pair": "K31/K2",   "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "branch", "structural_change": "residue addition", "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3C"},
    {"cps_pair": "K47/K2",   "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "core",   "structural_change": "residue addition", "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3D"},
    {"cps_pair": "K54/K43",  "serotype": "distinct",  "similarity_location": "near-identical", "dissimilarity_location": "core",   "structural_change": "residue addition", "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3E"},
    {"cps_pair": "K48/K15",  "serotype": "distinct",  "similarity_location": "branch",         "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3F"},
    {"cps_pair": "K48/K1",   "serotype": "distinct",  "similarity_location": "branch",         "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3F"},
    {"cps_pair": "K57/K43",  "serotype": "distinct",  "similarity_location": "branch",         "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3G"},
    {"cps_pair": "K3/K57",   "serotype": "distinct",  "similarity_location": "branch",         "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3G"},
    {"cps_pair": "K32/K43",  "serotype": "distinct",  "similarity_location": "core",           "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3H"},
    {"cps_pair": "K72/K2",   "serotype": "distinct",  "similarity_location": "core",           "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3I"},
    {"cps_pair": "K64/K63",  "serotype": "distinct",  "similarity_location": "core",           "dissimilarity_location": None,    "structural_change": None,              "modification_change": None,           "modification_location": None,     "structures_comparison": "Fig3J"},
]


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


def plot_figure4_2_panelC(output_path: Path, style=None) -> None:
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
    col_widths = [0.85, 1.2, 0.85, 2.0, 3.2, 0.65]

    # Per-column char wrap widths — wide columns wrap late to avoid unnecessary line breaks
    col_wrap_widths = [12, 16, 12, 28, 42, 8]

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
    ]
    ax.legend(
        handles=legend_elements, loc="lower center",
        bbox_to_anchor=(0.5, -0.07), ncol=5,
        fontsize=8.5, frameon=False,
    )

    dpi = style.dpi if style else 200
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {output_path}")
