"""
Build literature tables from S1_Table.xlsx.

Outputs
-------
rbp_deacetylases/deacetylases_literature.tsv  n=4  — proteins with _DAC_ in proteinID
cps_acetylases/acetylases_literature_active.tsv n=3  — proteins with _MOD_AC_ in proteinID
"""

from pathlib import Path

from enzymes_table import load_enzymes_table


def build_literature_tables(
    enzymes_xlsx: Path,
    rbp_deacetylases_dir: Path,
    cps_acetylases_dir: Path,
) -> None:
    df = load_enzymes_table(enzymes_xlsx)

    rbp_deacetylases_dir.mkdir(parents=True, exist_ok=True)
    cps_acetylases_dir.mkdir(parents=True, exist_ok=True)

    deacetylases = df[df["proteinid"].str.contains("_DAC_", na=False)].copy()
    path_deac = rbp_deacetylases_dir / "deacetylases_literature.tsv"
    deacetylases.to_csv(path_deac, sep="\t", index=False)
    print(f"  deacetylases_literature.tsv — {len(deacetylases)} rows → {path_deac}")

    acetylases = df[df["proteinid"].str.contains("_MOD_AC_", na=False)].copy()
    path_ac = cps_acetylases_dir / "acetylases_literature_active.tsv"
    acetylases.to_csv(path_ac, sep="\t", index=False)
    print(f"  acetylases_literature_active.tsv — {len(acetylases)} rows → {path_ac}")
