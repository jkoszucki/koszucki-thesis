"""Common base utilities for K-type analysis APIs."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

_DEFAULT_INPUT = Path(__file__).resolve().parent / "../capsule_structures/analysis/ktypes.xlsx"
_DEFAULT_OUTPUT = Path(__file__).resolve().parent / "../capsule_structures/analysis"


class BaseKTypeAPI:
    """Base class providing workbook access for K-type analyses."""

    def __init__(
        self,
        input_xlsx: Optional[Path | str] = None,
        output_dir: Optional[Path | str] = None,
    ) -> None:
        self.input_xlsx = Path(input_xlsx or _DEFAULT_INPUT).expanduser().resolve()
        self.output_dir = Path(output_dir or _DEFAULT_OUTPUT).expanduser().resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._sheet_cache: Dict[str, pd.DataFrame] = {}
        self._load_workbook()

    def _load_workbook(self) -> None:
        if not self.input_xlsx.exists():
            raise FileNotFoundError(f"Input workbook not found: {self.input_xlsx}")
        # Pre-load sheet names for quick validation. DataFrames are cached lazily.
        self._sheet_names = set(pd.ExcelFile(self.input_xlsx).sheet_names)

    def get_sheet(self, sheet_name: str) -> pd.DataFrame:
        if sheet_name not in self._sheet_cache:
            if sheet_name not in getattr(self, "_sheet_names", set()):
                raise KeyError(f"Sheet '{sheet_name}' not found in {self.input_xlsx}")
            self._sheet_cache[sheet_name] = pd.read_excel(
                self.input_xlsx, sheet_name=sheet_name
            )
        return self._sheet_cache[sheet_name].copy()

    def save_dataframe(self, df: pd.DataFrame, filename: str) -> Path:
        path = self.output_dir / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)
        return path
