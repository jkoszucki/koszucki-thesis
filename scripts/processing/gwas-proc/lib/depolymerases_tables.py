"""
Build rbp_depolymerases/ output tables from PLOS Biology supplement tables.

Outputs
-------
depolymerases_virulent_active.tsv       n=58 — S1_Table.xlsx (all virulent-phage depolymerases)
depolymerases_gwas.tsv                  n=26 — S3_Table.xlsx (GWAS-predicted depolymerases)
depolymerases_manualsearch_active.tsv   n=14 — S4A manual section: expression level != 'none' AND active==True
depolymerases_manualsearch_inactive.tsv n=2  — S4A manual section: expression level != 'none' AND active==False
depolymerases_manualsearch_notproduced.tsv
                                        n=34 — S4A manual section: expression level == 'none'
depolymerases_gwas_active.tsv           n=6  — S4A GWAS section (after 'GWAS proteins:' marker): active==True

S4A split:
  manual section = rows BEFORE 'GWAS proteins:' marker; NaN spacer excluded.
    expression level == 'none'             → not produced
    expression level != 'none', active=T   → produced, active
    expression level != 'none', active=F   → produced, inactive
  GWAS section = rows AFTER 'GWAS proteins:' marker; active==True → gwas_active
"""

from pathlib import Path
import pandas as pd


def build_depolymerases_tables(plos_tables_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    _write_virulent(plos_tables_dir, output_dir)
    _write_gwas(plos_tables_dir, output_dir)
    _write_manualsearch(plos_tables_dir, output_dir)


def _write_virulent(plos_dir: Path, out_dir: Path) -> None:
    df = pd.read_excel(plos_dir / "S1_Table.xlsx")
    path = out_dir / "depolymerases_virulent_active.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"  depolymerases_virulent_active.tsv — {len(df)} rows → {path}")


def _write_gwas(plos_dir: Path, out_dir: Path) -> None:
    df = pd.read_excel(plos_dir / "S3_Table.xlsx")
    path = out_dir / "depolymerases_gwas.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"  depolymerases_gwas.tsv — {len(df)} rows → {path}")


def _write_manualsearch(plos_dir: Path, out_dir: Path) -> None:
    xlsx = plos_dir / "S4_Table.xlsx"

    # Locate "GWAS proteins:" marker in column 0 (row 0 = header in raw read)
    df_raw = pd.read_excel(xlsx, sheet_name="S4A_Produced_depolymerases", header=None)
    gwas_marker_idx = df_raw[
        df_raw.iloc[:, 0].astype(str).str.strip() == "GWAS proteins:"
    ].index[0]

    df = pd.read_excel(xlsx, sheet_name="S4A_Produced_depolymerases")
    # Manual section: rows before the GWAS marker, drop NaN spacer rows
    # gwas_marker_idx (0-indexed, row 0 = header) → position gwas_marker_idx-1 in header-df
    df_manual = df.iloc[: gwas_marker_idx - 1].dropna(subset=["proteinID"]).copy()

    produced = df_manual[df_manual["expression level"] != "none"]
    not_produced = df_manual[df_manual["expression level"] == "none"]

    active = produced[produced["active"] == True].copy()
    inactive = produced[produced["active"] == False].copy()

    # GWAS section: rows from the marker onward, active==True
    df_gwas_section = df.iloc[gwas_marker_idx:].dropna(subset=["proteinID"]).copy()
    gwas_active = df_gwas_section[df_gwas_section["active"] == True].copy()

    for fname, subset in [
        ("depolymerases_manualsearch_active.tsv", active),
        ("depolymerases_manualsearch_inactive.tsv", inactive),
        ("depolymerases_manualsearch_notproduced.tsv", not_produced),
        ("depolymerases_gwas_active.tsv", gwas_active),
    ]:
        path = out_dir / fname
        subset.to_csv(path, sep="\t", index=False)
        print(f"  {fname} — {len(subset)} rows → {path}")
