"""
Single reader for the experimental-enzyme table (S1_Table.xlsx).

Three modules consume this spreadsheet — enzymes-proc (literature tables), acetyl-proc
(source_experimental flag) and figures/chapter4 (figure 4.1B network nodes) — and it is
hand-edited between pipeline runs, so both of its quirks are handled once, here:

* **Sheet name drift.** The sheet has been called both `enzymes` and `S1_Table`; saving
  from Excel has renamed it at least once, which broke all three consumers
  simultaneously. Candidate names are tried in turn, then a single-sheet workbook is
  accepted whatever its sheet is called.
* **Whitespace inside sequences.** PROTEIN02_MOD_AC_K2 carried a space at position 130,
  inflating its reported length by one and corrupting anything consuming the sequence
  verbatim (AF3 submissions, BLAST, FoldSeek queries). The spreadsheet is read-only
  input, so the cleanup lives here.

Callers get a ready-to-use DataFrame:

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
    from enzymes_table import load_enzymes_table

    df = load_enzymes_table(enzymes_xlsx)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

_SHEET_CANDIDATES = ("enzymes", "S1_Table")


def _resolve_sheet(xlsx: Path) -> str:
    sheets = pd.ExcelFile(xlsx).sheet_names
    for candidate in _SHEET_CANDIDATES:
        if candidate in sheets:
            return candidate
    if len(sheets) == 1:
        return sheets[0]
    raise KeyError(
        f"{xlsx.name}: no sheet named {' / '.join(_SHEET_CANDIDATES)} and the workbook "
        f"has several sheets ({', '.join(sheets)}) — cannot guess which one to read."
    )


def load_enzymes_table(xlsx: Path, verbose: bool = True) -> pd.DataFrame:
    """Read S1_Table.xlsx and return its enzyme rows with sequences cleaned."""
    sheet = _resolve_sheet(Path(xlsx))
    df = pd.read_excel(xlsx, sheet_name=sheet, engine="openpyxl")

    if "sequence" in df.columns:
        cleaned = df["sequence"].astype(str).str.replace(r"\s+", "", regex=True)
        # astype(str) turns blank cells into the string "nan"; put them back.
        cleaned = cleaned.where(df["sequence"].notna(), df["sequence"])
        changed = df.loc[cleaned.ne(df["sequence"]) & df["sequence"].notna(), "proteinid"]
        if verbose and len(changed):
            print(f"  [fix] whitespace stripped from sequence of: {', '.join(changed)}")
        df["sequence"] = cleaned

    return df
