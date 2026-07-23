"""API for accessing K-type modification metadata."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ktypes_process import KTypeTablesAPI


class KTypeModificationsAPI(KTypeTablesAPI):
    """Expose helper methods focused on the modifications worksheet."""

    MODIFICATIONS_FILENAME = "ktypes_modifications.csv"

    def build_modifications_table(self) -> pd.DataFrame:
        """Return the modifications table with structure_id at the front.

        If the modifications sheet already contains a structure_id column it is used directly —
        each row maps 1-to-1 to a structure.  The fallback K-type join is only applied when
        structure_id is absent from the modifications sheet.
        """
        mods = self.get_modifications()
        raw = self.get_sheet(self.KTYPE_SHEET)

        if "structure_id" in mods.columns:
            # structure_id already present — just move it to the front.
            cols = ["structure_id"] + [c for c in mods.columns if c != "structure_id"]
            return mods[cols]

        if "structure_id" not in raw.columns:
            return mods

        k_col = next(
            (c for c in raw.columns if str(c).strip().lower().replace("-", "_") in {"k_type", "ktype"}),
            None,
        )
        mod_k_col = next(
            (c for c in mods.columns if str(c).strip().lower().replace("-", "_") in {"k_type", "ktype"}),
            None,
        )
        if not k_col or not mod_k_col:
            return mods

        id_lookup = raw[[k_col, "structure_id"]].copy()
        id_lookup[k_col] = id_lookup[k_col].astype(str).str.strip()
        mods[mod_k_col] = mods[mod_k_col].astype(str).str.strip()

        mods = mods.merge(id_lookup, left_on=mod_k_col, right_on=k_col, how="left")
        if k_col != mod_k_col and k_col in mods.columns:
            mods = mods.drop(columns=[k_col])

        cols = ["structure_id"] + [c for c in mods.columns if c != "structure_id"]
        return mods[cols]

    def export_modifications_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        table = df if df is not None else self.build_modifications_table()
        return self.save_dataframe(table, filename or self.MODIFICATIONS_FILENAME)
