"""
Build acetylases_kloci.tsv from S3_Table.xlsx.

Source: the `acetyltrasferases` sheet (69 rows) is a pre-merged summary of the
three detection methods (keyword / blastp / foldseek), each already flagged
as a boolean column, plus the candidate's sequence. Product annotation is
looked up from the `keyword` and `foldseek` sheets (not present on the merged
sheet itself).

  acetyltrasferases sheet (69 rows, pre-merged)  → source_keyword/source_blastp/source_foldseek flags + sequence
  keyword / foldseek sheets                      → product_kaptive lookup by protein id
  S1_Table.xlsx enzymes sheet (experimental K-types) → source_experimental=True

Output columns:
  protein_id, product_kaptive, protein_length, sequence,
  source_experimental, source_blastp, source_keyword, source_foldseek
"""

from __future__ import annotations
from pathlib import Path

import pandas as pd
import re


def _ktype_to_klocus(ktype: str) -> str:
    """'K2' → 'KL2', 'K57' → 'KL57'."""
    m = re.match(r"^K(\d+)$", str(ktype).strip())
    return f"KL{m.group(1)}" if m else ""


def build_acetylases_kloci(
    s3_table_xlsx: Path,
    enzymes_xlsx: Path,
) -> pd.DataFrame:
    """
    Return the 69-row acetylases_kloci DataFrame.

    Parameters
    ----------
    s3_table_xlsx
        Path to S3_Table.xlsx (sheets: acetyltrasferases, keyword, blastp, foldseek).
    enzymes_xlsx
        Path to S1_Table.xlsx; the enzymes sheet lists experimental K-types.
    """
    # ------------------------------------------------------------------
    # 1. Pre-merged candidate table (69 rows, source flags already computed)
    # ------------------------------------------------------------------
    at_df = pd.read_excel(s3_table_xlsx, sheet_name="acetyltrasferases")

    # ------------------------------------------------------------------
    # 2. Product annotation lookup — keyword sheet, falling back to foldseek
    # ------------------------------------------------------------------
    kw_df = pd.read_excel(s3_table_xlsx, sheet_name="keyword")
    kw_product = kw_df.set_index("locus_tag")["product"]

    fs_df = pd.read_excel(s3_table_xlsx, sheet_name="foldseek")
    fs_product = (
        fs_df[["query", "product"]]
        .dropna(subset=["product"])
        .drop_duplicates(subset=["query"])
        .set_index("query")["product"]
    )

    def _lookup_product(pid: str) -> str:
        if pid in kw_product.index:
            return str(kw_product.at[pid]).strip()
        if pid in fs_product.index:
            return str(fs_product.at[pid]).strip()
        return ""

    rows = []
    for _, row in at_df.iterrows():
        pid = str(row["proteinid"]).strip()
        seq = str(row["sequence"]).strip()
        rows.append({
            "protein_id": pid,
            "product_kaptive": _lookup_product(pid),
            "protein_length": len(seq),
            "sequence": seq,
            "source_blastp": bool(row["blastp"]),
            "source_keyword": bool(row["keyword"]),
            "source_foldseek": bool(row["foldseek"]),
        })
    combined = pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # 3. source_experimental flag — from S1_Table.xlsx enzymes sheet
    # ------------------------------------------------------------------
    exp_df = pd.read_excel(enzymes_xlsx, sheet_name="enzymes")
    # Only acetylation enzymes (not deacetylases) set source_experimental
    acetyl_exp_df = exp_df[exp_df["modification"].str.lower() == "acetylation"]
    # Derive K-locus prefix (e.g. "K2" → "KL2") from the ktype column
    exp_kloci = set()
    for ktype in acetyl_exp_df["ktype"].dropna():
        kl = _ktype_to_klocus(str(ktype))
        if kl:
            exp_kloci.add(kl)

    def _is_experimental(pid: str) -> bool:
        m = re.match(r"^(KL\d+)_", pid)
        return m.group(1) in exp_kloci if m else False

    combined["source_experimental"] = combined["protein_id"].apply(_is_experimental)

    # ------------------------------------------------------------------
    # 4. Final column order
    # ------------------------------------------------------------------
    return combined[[
        "protein_id", "product_kaptive", "protein_length", "sequence",
        "source_experimental", "source_blastp", "source_keyword", "source_foldseek",
    ]]
