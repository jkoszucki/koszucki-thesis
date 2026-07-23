"""API for computing K-type similarity metrics."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ktypes_process import KTypeTablesAPI


class KTypeSimilarityAPI(KTypeTablesAPI):
    """Convenience wrapper dedicated to similarity analyses."""

    def build_similarity_table(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        return self.compute_similarity(df)

    def export_similarity_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        return super().export_similarity_table(df=df, filename=filename)
