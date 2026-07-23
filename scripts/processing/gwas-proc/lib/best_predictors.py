"""
Build rbp_best_predictors/best_predictors_gwas.tsv from gwas_hits.tsv.

Selection logic (from manuscript):
  1. Start from all GWAS hits in gwas_hits.tsv (pvalue_corr already ≤ 0.05, across all 6 clustering levels)
  2. Filter: precision ≥ 0.8, F1_score ≥ 0.5, MCC ≥ 0.5
  3. Per K-locus: select the row with the highest F1_score × MCC product
  Result: n=16 best-performing predictors across 16 K-loci.

Note: gwas_hits.tsv uses precision ≥ 0.60 as its baseline threshold; here
the stricter 0.80 threshold is applied to identify the best per-locus predictors.
"""

from pathlib import Path
import pandas as pd


def build_best_predictors(gwas_hits_tsv: Path, output_dir: Path) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(gwas_hits_tsv, sep="\t")

    filtered = df[
        (df["precision"] >= 0.8) &
        (df["F1_score"] >= 0.5) &
        (df["MCC"] >= 0.5)
    ].copy()

    filtered["_f1_mcc"] = filtered["F1_score"] * filtered["MCC"]
    best = filtered.loc[filtered.groupby("locus")["_f1_mcc"].idxmax()].copy()
    best = best.drop(columns="_f1_mcc").reset_index(drop=True)

    out_path = output_dir / "best_predictors_gwas.tsv"
    best.to_csv(out_path, sep="\t", index=False)
    print(f"  best_predictors_gwas.tsv — {len(best)} rows → {out_path}")

    return best
