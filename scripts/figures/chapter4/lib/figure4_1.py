"""Plotting API exposing streamlined helpers for the K-type datasets.

This module consolidates the three visualisation utilities that remain from the
legacy plotting scripts:

* :meth:`KTypePlotAPI.plot_cumulative_structures` plots the cumulative number of
  K-types discovered per year.
* :meth:`KTypePlotAPI.plot_monosaccharide_prevalence` summarises how frequently
  monosaccharides and their acetylation (OAc) / pyruvation (OPy) modifications
  occur across the catalogue.
* :meth:`KTypePlotAPI.plot_similarity_histogram`,
  :meth:`KTypePlotAPI.plot_similarity_scatter` and
  :meth:`KTypePlotAPI.export_similarity_network` operate on the pairwise
  similarity table used to create Cytoscape networks.

The module is intentionally self-contained so that notebooks or command line
scripts can import :class:`KTypePlotAPI` directly from the repository root.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


@dataclass
class PrevalenceResult:
    """Return value for :meth:`KTypePlotAPI.plot_monosaccharide_prevalence`.

    ``data`` stores a tidy representation of the prevalence table.  Each row
    corresponds to a category (``Py``/``Ac``/monosaccharide) and one modification
    state.  ``figure`` contains the rendered stacked bar chart.  The
    ``presence`` array mirrors the order of categories in ``data`` and exposes
    the raw prevalence ratios for convenience.
    """

    data: pd.DataFrame
    presence: np.ndarray
    figure: plt.Figure


class KTypePlotAPI:
    """High level plotting helpers for the K-type tables."""

    def __init__(
        self,
        *,
        ktypes_csv: Optional[Path | str] = None,
        similarity_csv: Optional[Path | str] = None,
        modifications_csv: Optional[Path | str] = None,
        style=None,
    ) -> None:
        self._ktypes_csv = Path(ktypes_csv) if ktypes_csv else None
        self._similarity_csv = Path(similarity_csv) if similarity_csv else None
        self._modifications_csv = Path(modifications_csv) if modifications_csv else None
        self._ktypes_df: Optional[pd.DataFrame] = None
        self._similarity_df: Optional[pd.DataFrame] = None
        self._modifications_df: Optional[pd.DataFrame] = None
        self._style = style


    # ------------------------------------------------------------------
    # Public helpers operating on the K-type catalogue
    # ------------------------------------------------------------------
    def plot_cumulative_structures(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        year_col: str = "Year",
        source_col: str = "Source",
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Plot cumulative K-type structures grouped by ``source_col``.

        The dataset is split according to the values contained in ``source_col``.
        Separate cumulative curves are drawn for each source so that the
        contribution from KPAM and the literature can be tracked over time.  The
        returned dataframe contains one row per source/year combination with the
        raw counts and cumulative totals.
        """

        df = self._load_ktypes_dataframe(input_csv)
        if year_col not in df.columns:
            raise ValueError(f"Column '{year_col}' not found in the dataset")


        if source_col not in df.columns:
            raise ValueError(f"Column '{source_col}' not found in the dataset")

        working = df[[year_col, source_col]].copy()
        working[year_col] = pd.to_numeric(working[year_col], errors="coerce").astype("Int64")
        working = working.dropna(subset=[year_col])
        working[year_col] = working[year_col].astype(int)

        grouped = (
            working.groupby([source_col, year_col])
            .size()
            .rename("count")
            .reset_index()
            .sort_values([source_col, year_col])
            .reset_index(drop=True)
        )

        grouped["cumulative"] = grouped.groupby(source_col)["count"].cumsum()

        st = self._style
        dpi = st.dpi if st else 200
        label_fs = st.axis_label_fontsize if st else 12
        label_fw = st.axis_label_fontweight if st else "bold"
        tick_fs = st.tick_fontsize if st else 8
        tick_fw = st.tick_fontweight if st else "bold"
        kpam_color = st.kpam_color if st else "#1f77b4"
        lit_color = st.lit_color if st else "#d62728"

        fig, ax = plt.subplots(figsize=(4, 4))
        colours = {"KPAM": kpam_color, "LIT": lit_color}
        labels = {"KPAM": "k-pam database", "LIT": "literature review"}

        for source, data in grouped.groupby(source_col):
            color = colours.get(source, None)
            label = labels.get(source, source)
            ax.plot(
                data[year_col],
                data["cumulative"],
                marker="o",
                linestyle="--",
                label=label,
                color=color,
            )

        ax.set_xlabel("Year", fontsize=label_fs, fontweight=label_fw)
        ax.set_ylabel("Cumulative # of CPS structures", fontsize=label_fs, fontweight=label_fw)
        ax.grid(True, linestyle=":")
        ax.legend()

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(tick_fs)
            label.set_fontweight(tick_fw)

        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
        if show:
            plt.show()
        else:
            plt.close(fig)

        return grouped



    def plot_modification_and_monosaccharide_distribution(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        modifications_csv: Optional[Path | str] = None,
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> PrevalenceResult:
        """Render the two panel summary for modifications and monosaccharides.

        Parameters
        ----------
        input_csv:
            Optional override for the processed K-type catalogue.  When omitted,
            the dataframe supplied at construction time is reused.
        modifications_csv:
            Optional override for the modifications table.  If neither the
            parameter nor the constructor attribute are provided the plot falls
            back to an empty modifications table (i.e. only grey / white
            segments are rendered).
        """

        dpi = self._style.dpi if self._style else 200

        data = self._prepare_prevalence_data(
            input_csv=input_csv,
            modifications_csv=modifications_csv,
        )

        figure = self._render_modification_distribution_chart(data)

        if output_path:
            figure.savefig(output_path, bbox_inches="tight", dpi=dpi)
        if show:
            plt.show()
        else:
            plt.close(figure)

        tidy = self._build_prevalence_table(data)
        return PrevalenceResult(data=tidy, presence=data.presence, figure=figure)


    def plot_monosaccharide_combinations(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        backbone_col: str = "backbone_linkages_coarse",
        branch_col: str = "branch",
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Plot the unique monosaccharide combinations and their frequencies."""

        df = self._load_ktypes_dataframe(input_csv)

        if backbone_col not in df.columns:
            raise ValueError(
                f"Column '{backbone_col}' not found in the K-type dataset"
            )
        if branch_col not in df.columns:
            raise ValueError(f"Column '{branch_col}' not found in the K-type dataset")

        combinations: List[Tuple[str, ...]] = []

        for _, row in df[[backbone_col, branch_col]].iterrows():
            backbone_monos = self._extract_backbone_monos(row[backbone_col])
            branch_monos = self._extract_branch_monos(row[branch_col])
            combo = tuple(sorted(set(backbone_monos + branch_monos)))
            if combo:
                combinations.append(combo)

        if not combinations:
            raise ValueError("No monosaccharide combinations could be extracted")

        counts = Counter(combinations)
        ordered = sorted(counts.items(), key=lambda item: item[1], reverse=True)

        unique_monos: List[str] = sorted({mono for combo, _ in ordered for mono in combo})
        mono_index = {mono: idx for idx, mono in enumerate(unique_monos)}

        matrix = np.zeros((len(ordered), len(unique_monos)), dtype=int)
        for row_idx, (combo, _) in enumerate(ordered):
            for mono in combo:
                matrix[row_idx, mono_index[mono]] = 1

        labels = [" + ".join(combo) if combo else "" for combo, _ in ordered]
        frequencies = [freq for _, freq in ordered]

        height = max(4.5, 0.35 * len(labels))
        fig, (ax_matrix, ax_count) = plt.subplots(
            1,
            2,
            figsize=(12, height),
            sharey=True,
            gridspec_kw={"width_ratios": [max(2.5, len(unique_monos) * 0.6), 1.0]},
        )

        ax_matrix.imshow(matrix, aspect="auto", cmap="Blues")
        ax_matrix.set_xticks(range(len(unique_monos)))
        ax_matrix.set_xticklabels(unique_monos, rotation=45, ha="right")
        ax_matrix.set_yticks(range(len(labels)))
        ax_matrix.set_yticklabels(labels)
        ax_matrix.set_xlabel("Monosaccharides")
        ax_matrix.set_ylabel("Monosaccharide combinations")
        ax_matrix.set_title("Combination composition")

        ax_count.barh(range(len(labels)), frequencies, color="#4c72b0")
        ax_count.set_xlabel("Frequency")
        ax_count.set_title("Combination frequency")
        ax_count.set_ylim(-0.5, len(labels) - 0.5)
        ax_count.grid(axis="x", linestyle=":", alpha=0.5)
        ax_count.set_yticks(range(len(labels)))
        ax_count.set_yticklabels([])

        for spine in ["top", "right"]:
            ax_matrix.spines[spine].set_visible(False)
            ax_count.spines[spine].set_visible(False)

        fig.tight_layout()

        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=300)
        if show:
            plt.show()
        else:
            plt.close(fig)

        summary = pd.DataFrame(
            {
                "combination": labels,
                "monosaccharides": [list(combo) for combo, _ in ordered],
                "count": frequencies,
            }
        )
        return summary


    # ------------------------------------------------------------------
    # Public helpers operating on the similarity table
    # ------------------------------------------------------------------

    def plot_composition_combinations(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        backbone_col: str = "backbone_linkages_coarse",
        branch_col: str = "branch",
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Count unique monosaccharide combinations for the K-type catalogue.

        Monosaccharides are extracted with straightforward ``split`` / ``strip``
        operations.  ``-OUT`` suffixes are removed from ``backbone_col`` entries
        while anchors (the final residues) are dropped from ``branch_col``
        values.  Edge descriptors inside brackets are ignored which yields the
        bare monosaccharide names.  The method returns the frequency table used
        for plotting for convenient inspection.
        """

        df = self._load_ktypes_dataframe(input_csv)
        for column in (backbone_col, branch_col):
            if column not in df.columns:
                raise ValueError(f"Column '{column}' not found in the dataset")

        combinations: List[str] = []
        for backbone_value, branch_value in zip(df[backbone_col], df[branch_col]):
            backbone_monos = self._extract_backbone_monos(backbone_value)
            branch_monos = self._extract_branch_monos(branch_value)
            unique = sorted({*backbone_monos, *branch_monos})
            if not unique:
                continue
            combinations.append(
                ", ".join(unique)
            )

        if not combinations:
            raise ValueError("No monosaccharide combinations could be parsed")

        counts = (
            pd.Series(combinations)
            .value_counts()
            .rename_axis("composition")
            .reset_index(name="frequency")
        )

        fig_height = max(4, 0.45 * len(counts))
        fig, ax = plt.subplots(figsize=(10, fig_height))
        ax.barh(counts["composition"], counts["frequency"], color="#4c72b0")
        ax.set_xlabel("Count")
        ax.set_ylabel("K-type composition")
        ax.set_title("Unique monosaccharide combinations")
        ax.invert_yaxis()
        ax.grid(axis="x", linestyle=":", alpha=0.6)

        fig.tight_layout()
        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=300)
        if show:
            plt.show()
        else:
            plt.close(fig)

        return counts


    def plot_similarity_histogram(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        column: str = "jaccard_similarity",
        bins: int = 20,
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> pd.Series:
        """Plot a histogram for ``column`` from the similarity table."""

        df = self._load_similarity_dataframe(input_csv)
        if column not in df.columns:
            raise ValueError(f"Column '{column}' is missing from the similarity table")

        st = self._style
        dpi = st.dpi if st else 200
        label_fs = st.axis_label_fontsize if st else 12
        label_fw = st.axis_label_fontweight if st else "bold"
        tick_fs = st.tick_fontsize if st else 8
        tick_fw = st.tick_fontweight if st else "bold"

        series = pd.to_numeric(df[column], errors="coerce").dropna()
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.hist(series, bins=bins, edgecolor="black")
        ax.set_xlabel(column, fontsize=label_fs, fontweight=label_fw)
        ax.set_ylabel("count", fontsize=label_fs, fontweight=label_fw)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(tick_fs)
            label.set_fontweight(tick_fw)

        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=dpi)
        if show:
            plt.show()
        else:
            plt.close(fig)

        return series

    def plot_similarity_scatter(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        x_column: str = "weighted_core",
        y_column: str = "weighted_branch",
        output_path: Optional[Path | str] = None,
        show: bool = True,
    ) -> pd.DataFrame:
        """Scatter plot of two similarity metrics against each other."""

        df = self._load_similarity_dataframe(input_csv)
        missing = [c for c in (x_column, y_column) if c not in df.columns]
        if missing:
            raise ValueError(
                "Columns {} are missing from the similarity table".format(
                    ", ".join(missing)
                )
            )

        x = pd.to_numeric(df[x_column], errors="coerce")
        y = pd.to_numeric(df[y_column], errors="coerce")

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(x, y, alpha=0.6, edgecolor="k")
        ax.set_xlabel(x_column)
        ax.set_ylabel(y_column)
        ax.set_title(f"{x_column} vs {y_column}")
        ax.grid(True, linestyle=":")
        ax.axhline(0, color="gray", linewidth=0.7)
        ax.axvline(0, color="gray", linewidth=0.7)

        if output_path:
            fig.savefig(output_path, bbox_inches="tight", dpi=200)
        if show:
            plt.show()
        else:
            plt.close(fig)

        return pd.DataFrame({x_column: x, y_column: y})

    def export_similarity_network(
        self,
        *,
        input_csv: Optional[Path | str] = None,
        output_path: Optional[Path | str] = None,
        min_weighted_core: Optional[float] = None,
        max_weighted_core: Optional[float] = None,
        min_weighted_branch: Optional[float] = None,
        max_weighted_branch: Optional[float] = None,
    ) -> pd.DataFrame:
        """Filter similarity edges based on weighted score thresholds."""

        df = self._load_similarity_dataframe(input_csv)
        for col in ("weighted_core", "weighted_branch"):
            if col not in df.columns:
                raise ValueError(f"Column '{col}' is required in the similarity table")

        if min_weighted_core is not None:
            df = df[df["weighted_core"] >= min_weighted_core]
        if max_weighted_core is not None:
            df = df[df["weighted_core"] <= max_weighted_core]
        if min_weighted_branch is not None:
            df = df[df["weighted_branch"] >= min_weighted_branch]
        if max_weighted_branch is not None:
            df = df[df["weighted_branch"] <= max_weighted_branch]

        df = df.reset_index(drop=True)

        near_identical = df['weighted_total'] >= 0.9
        very_similar = (df['weighted_total'] >= 0.8) & (df['weighted_total'] < 0.9)
        somewhat_similar = df['weighted_total'] < 0.8

        df['similarity_cat'] = 'unassigned'
        df.loc[near_identical, 'similarity_cat'] = 'near-identical (total ≥ 0.9)'
        df.loc[very_similar, 'similarity_cat'] = 'very-similar with branch variation (total 0.8–0.9)'
        df.loc[somewhat_similar, 'similarity_cat'] = 'somewhat-similar with branch variation (total < 0.8)'

        

        if output_path:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, index=False)
        return df

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    @dataclass
    class _PrevalenceData:
        column_names: List[str]
        monos: List[str]
        presence: np.ndarray
        type_splits: List[Dict[object, float]]
        total_ktypes: int
        mono_counts: List[int]

    def _prepare_prevalence_data(
        self,
        *,
        input_csv: Optional[Path | str],
        modifications_csv: Optional[Path | str],
    ) -> "KTypePlotAPI._PrevalenceData":
        df = self._load_ktypes_dataframe(input_csv)
        if df.empty:
            raise ValueError("K-type table is empty")

        ktype_col = self._resolve_k_type_column(df)
        if ktype_col is None:
            raise ValueError("Could not identify the K-type column in the dataset")

        monos_by_type: Dict[str, set[str]] = {}
        for _, row in df.iterrows():
            ktype = str(row[ktype_col]).strip()
            if not ktype:
                continue
            monos = monos_by_type.setdefault(ktype, set())
            for column in ("core_monos_all", "branch_monos_all", "core_monos_unique", "branch_monos_unique"):
                if column in df.columns:
                    monos.update(self._split_processed_monos(row[column]))

        unique_ktypes = sorted(monos_by_type.keys())
        total_ktypes = len(unique_ktypes)
        if total_ktypes == 0:
            raise ValueError("No K-types with monosaccharide information were found")

        mono_counter: Counter[str] = Counter()
        for monos in monos_by_type.values():
            for mono in monos:
                mono_counter[mono] += 1
                
        monos = [mono for mono, _ in sorted(mono_counter.items(), key=lambda item: (-item[1], item[0]))]

        # Track which ktypes have pyruvylation / acetylation globally and per mono
        has_pyro: set[str] = set()
        has_ac: set[str] = set()
        pyro_ktypes_by_mono: Dict[str, set[str]] = {}
        ac_ktypes_by_mono: Dict[str, set[str]] = {}

        mods_df = self._load_modifications_dataframe(modifications_csv)
        if not mods_df.empty:
            mod_ktype_col = self._resolve_k_type_column(mods_df)
            if mod_ktype_col is None:
                raise ValueError(
                    "Could not identify the K-type column in the modifications table"
                )

            for _, row in mods_df.iterrows():
                ktype = str(row[mod_ktype_col]).strip()
                if not ktype or ktype not in monos_by_type:
                    continue
                mod_name = str(row.get("modification", "")).strip().lower()
                if mod_name not in {"pyruvylation", "acetylation"}:
                    continue
                mono = self._normalise_monosaccharide_name(row.get("monosaccharide"))

                if mod_name == "pyruvylation":
                    has_pyro.add(ktype)
                    if mono:
                        pyro_ktypes_by_mono.setdefault(mono, set()).add(ktype)
                else:
                    has_ac.add(ktype)
                    if mono:
                        ac_ktypes_by_mono.setdefault(mono, set()).add(ktype)

        # Build stacked ratios for the global modifications bar
        column_names = ["Pyruvylation & Acetylation"] + monos
        presence = np.zeros(len(column_names), dtype=float)
        type_splits: List[Dict[str, float]] = []

        pyro_only = has_pyro - has_ac
        ac_only = has_ac - has_pyro
        both = has_pyro & has_ac
        n_unmodified = total_ktypes - len(has_pyro | has_ac)

        combined_split: Dict[str, float] = {}
        if pyro_only:
            combined_split["pyruvylation"] = len(pyro_only) / total_ktypes
        if ac_only:
            combined_split["acetylation"] = len(ac_only) / total_ktypes
        if both:
            combined_split["_both"] = len(both) / total_ktypes
        if n_unmodified > 0:
            combined_split["_unmodified"] = n_unmodified / total_ktypes

        presence[0] = 1.0 - n_unmodified / total_ktypes
        type_splits.append(combined_split)

        # Per-monosaccharide: exclusive breakdown by which ktypes carry each modification
        for idx, mono in enumerate(monos, start=1):
            presence_ratio = mono_counter[mono] / total_ktypes
            mono_ktypes = {k for k, ms in monos_by_type.items() if mono in ms}
            py_here = pyro_ktypes_by_mono.get(mono, set()) & mono_ktypes
            ac_here = ac_ktypes_by_mono.get(mono, set()) & mono_ktypes
            both_here = py_here & ac_here
            py_only_here = py_here - both_here
            ac_only_here = ac_here - both_here
            unmod_here = mono_ktypes - py_here - ac_here

            split: Dict[str, float] = {}
            if py_only_here:
                split["pyruvylation"] = len(py_only_here) / total_ktypes
            if ac_only_here:
                split["acetylation"] = len(ac_only_here) / total_ktypes
            if both_here:
                split["_both"] = len(both_here) / total_ktypes
            if unmod_here:
                split["_unmodified"] = len(unmod_here) / total_ktypes

            presence[idx] = presence_ratio
            type_splits.append(split)

        mono_counts = [mono_counter[mono] for mono in monos]

        return KTypePlotAPI._PrevalenceData(
            column_names=column_names,
            monos=monos,
            presence=presence,
            type_splits=type_splits,
            total_ktypes=total_ktypes,
            mono_counts=mono_counts,
        )

    def _build_prevalence_table(
        self, data: "KTypePlotAPI._PrevalenceData"
    ) -> pd.DataFrame:
        records = []
        for name, presence, split in zip(
            data.column_names, data.presence, data.type_splits
        ):
            for modification, share in split.items():
                if modification == "_unmodified":
                    label = "unmodified"
                elif modification == "_both":
                    label = "both modifications"
                else:
                    label = str(modification)
                records.append(
                    {
                        "category": name,
                        "modification": label,
                        "share": share,
                        "presence": presence,
                    }
                )
        return pd.DataFrame(records)


    def _render_modification_distribution_chart(
        self, data: "KTypePlotAPI._PrevalenceData"
    ) -> plt.Figure:
        st = self._style
        pyro_color = st.pyruvylation_color if st else "#4f81bd"
        ac_color = st.acetylation_color if st else "#f2c94c"
        both_color = st.both_color if st else "#DC143C"
        gray_color = st.gray_color if st else "#bfbfbf"
        count_fs = st.tick_fontsize if st else 8

        width_mod = 0.6
        width_mono = 0.8
        fig_width = max(10.0, len(data.monos) * 0.6 + 4.0)
        fig, (ax_mod, ax_mono) = plt.subplots(
            1, 2, sharey=True, figsize=(fig_width, 5.5), width_ratios=[1, 4]
        )

        def _render_split_bar(ax, x, split, width):
            bottom = 0.0
            color_map = {
                "pyruvylation": pyro_color,
                "acetylation": ac_color,
                "_both": both_color,
                "_unmodified": gray_color,
            }
            for key in ("pyruvylation", "acetylation", "_both", "_unmodified"):
                share = split.get(key, 0.0)
                if share <= 0:
                    continue
                ax.bar(x, share, bottom=bottom, width=width,
                       color=color_map[key], edgecolor="none")
                bottom += share
            if bottom < 1.0:
                ax.bar(x, 1.0 - bottom, bottom=bottom, width=width,
                       color="white", edgecolor="none")
            ax.bar(x, 1.0, bottom=0.0, width=width,
                   color="none", edgecolor="#333333", linewidth=0.8)

        ax_mod.set_title("Prevalence across\nK-types")
        ax_mod.set_xticks([0], ["pyruvylation / acetylation"])
        ax_mod.set_ylim(0.0, 1.18)
        ax_mod.set_ylabel("Ratio of CPS structures")

        combined_split = data.type_splits[0]
        _render_split_bar(ax_mod, 0, combined_split, width_mod)
        ax_mod.text(0, 1.03, f"n={data.total_ktypes}", ha="center", va="bottom",
                    fontsize=count_fs, fontweight="bold")

        ax_mod.set_xlim(-0.75, 0.75)
        ax_mod.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.7)

        ax_mono.set_title(
            "Monosaccharide prevalence\nand modification distribution"
        )
        ax_mono.set_xticks(range(len(data.monos)))
        ax_mono.set_xticklabels(data.monos, rotation=45, ha="right")
        ax_mono.set_ylim(0.0, 1.18)
        ax_mono.set_ylabel("")

        for idx, (mono, count) in enumerate(zip(data.monos, data.mono_counts)):
            split = data.type_splits[idx + 1]
            _render_split_bar(ax_mono, idx, split, width_mono)
            ax_mono.text(idx, 1.03, f"n={count}", ha="center", va="bottom",
                         fontsize=count_fs, fontweight="bold")

        ax_mono.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.7)
        ax_mono.set_xlim(-0.6, len(data.monos) - 0.4)

        legend_handles: List[Patch] = []
        has_pyro = combined_split.get("pyruvylation", 0.0) > 0 or combined_split.get("_both", 0.0) > 0
        has_ac = combined_split.get("acetylation", 0.0) > 0 or combined_split.get("_both", 0.0) > 0
        both_present = combined_split.get("_both", 0.0) > 0
        if has_pyro:
            legend_handles.append(Patch(facecolor=pyro_color, edgecolor="none", label="Pyruvylation"))
        if has_ac:
            legend_handles.append(Patch(facecolor=ac_color, edgecolor="none", label="Acetylation"))
        if both_present:
            legend_handles.append(Patch(facecolor=both_color, edgecolor="none", label="Both"))
        legend_handles.append(Patch(facecolor=gray_color, edgecolor="none", label="No modification"))
        legend_handles.append(
            Patch(facecolor="white", edgecolor="#333333", label="Residue absent", linewidth=0.8)
        )
        ax_mono.legend(
            legend_handles,
            [h.get_label() for h in legend_handles],
            loc=(1.06, 0.1),
            frameon=False,
        )

        fig.tight_layout()
        return fig

    # ----------------------------- parsing helpers -----------------------------
    @staticmethod
    def _resolve_k_type_column(df: pd.DataFrame) -> Optional[str]:
        # Prefer structure_id as the primary unique identifier
        for column in df.columns:
            if str(column).strip().lower() == "structure_id":
                return column
        for column in df.columns:
            normalized = str(column).replace("\ufeff", "").strip().lower()
            normalized = normalized.replace("-", "_").replace(" ", "_")
            if normalized in {"k_type", "ktype"}:
                return column
        return None

    @staticmethod
    def _split_processed_monos(value: object) -> List[str]:
        if pd.isna(value):
            return []
        text = str(value).strip()
        if not text or text == "0":
            return []
        text = text.replace("~~", ";").replace(",", ";")
        tokens: List[str] = []
        for chunk in text.split(";"):
            token = chunk.strip()
            if not token:
                continue
            token = token.split("@", 1)[0]
            token = token.split("(", 1)[0]
            token = token.strip()
            if not token or token.lower() == "sug":
                continue
            tokens.append(token)
        return tokens

    @staticmethod
    def _normalise_bond_type(value: object) -> str:
        if pd.isna(value):
            return "unknown"
        text = str(value).strip()
        if not text or text == "0":
            return "unknown"
        lowered = text.lower().replace(" ", "")
        lowered = lowered.replace("−", "-").replace("–", "-").replace("—", "-")
        if lowered == "unknown":
            return "unknown"
        cleaned = lowered.replace("0-", "o-")
        cleaned = cleaned.upper()
        return cleaned

    @staticmethod
    def _normalise_monosaccharide_name(value: object) -> Optional[str]:
        if pd.isna(value):
            return None
        text = str(value).strip()
        if not text or text == "0":
            return None
        name = text.split("@", 1)[0].strip()
        if not name or name.lower() == "unknown":
            return None
        return name

    @staticmethod
    def _extract_backbone_monos(value: object) -> List[str]:
        if pd.isna(value):
            return []
        text = str(value).strip()
        if not text:
            return []
        if "-OUT" in text:
            text = text.split("-OUT")[0].strip()
        return KTypePlotAPI._extract_monos_from_chain(text)

    @staticmethod
    def _extract_branch_monos(value: object) -> List[str]:
        if pd.isna(value):
            return []
        text = str(value).strip()
        if not text:
            return []
        monos = KTypePlotAPI._extract_monos_from_chain(text)
        if monos:
            return monos[:-1]
        return monos

    @staticmethod
    def _extract_monos_from_chain(text: str) -> List[str]:
        monos: List[str] = []
        for segment in text.split(";"):
            remaining = segment.strip()
            while "(" in remaining:
                prefix, remaining = remaining.split("(", 1)
                prefix = prefix.strip()
                prefix = prefix.strip(";,:")
                prefix = prefix.strip()
                if prefix:
                    monos.append(prefix)
                if ")" in remaining:
                    remaining = remaining.split(")", 1)[1].strip()
                else:
                    remaining = ""
            if remaining:
                tail = remaining.strip()
                tail = tail.strip(";,:")
                tail = tail.strip()
                if tail:
                    monos.append(tail)
        return monos


    def _load_ktypes_dataframe(self, input_csv: Optional[Path | str]) -> pd.DataFrame:
        if input_csv is not None:
            return pd.read_csv(input_csv)
        if self._ktypes_df is None:
            if self._ktypes_csv is None:
                raise ValueError("No K-type CSV provided")
            self._ktypes_df = pd.read_csv(self._ktypes_csv)
        return self._ktypes_df.copy()

    def _load_similarity_dataframe(
        self, input_csv: Optional[Path | str]
    ) -> pd.DataFrame:
        if input_csv is not None:
            return pd.read_csv(input_csv)
        if self._similarity_df is None:
            if self._similarity_csv is None:
                raise ValueError("No similarity CSV provided")
            self._similarity_df = pd.read_csv(self._similarity_csv)
        return self._similarity_df.copy()

    def _load_modifications_dataframe(
        self, input_csv: Optional[Path | str]
    ) -> pd.DataFrame:
        if input_csv is not None:
            return pd.read_csv(input_csv)
        if self._modifications_df is None:
            if self._modifications_csv is None:
                return pd.DataFrame()
            self._modifications_df = pd.read_csv(self._modifications_csv)
        return self._modifications_df.copy()


__all__ = ["KTypePlotAPI", "PrevalenceResult"]
