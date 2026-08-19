"""
Classify the GWAS acetyltransferase predictions in `acetylases_gwas.tsv`.

Step 2 of acetyl-proc selects GWAS hits whose PC80 carries SSLBH evidence (ECOD
single-stranded left-handed β-helix / PFAM / PHROGs), which is how nearly every
acetyltransferase prediction in the table was found. One prediction matters for a
different reason: it is the **best-performing GWAS predictor for its K-locus**
(highest F1 × MCC), and the manual finding behind figure 3.2 is that it turns out to be an
acetyltransferase despite carrying no ECOD annotation at all (`no-ecod-reported-topology`).

That distinction is what the `prediction_class` column records:

    sslbh-carrying   the PC80 carries SSLBH/PFAM/PHROGs acetyltransferase evidence
    best-performing  hand-curated: the top GWAS predictor for KL111

**KL111 / PC0915** is already in the SSLBH subset (its PC80 hits PFAM + PHROGs), so it is
only relabelled. The mechanism for pulling a missing best performer out of the flagged
pyseer table is kept below — it is what a future addition would need — but no current
entry uses it.

> **KL30 / PC0675 was listed here until 2026-07-27** as the second best-performing
> acetyltransferase (PROTEIN04_GWAS_AC_K30). Both structural methods refuted the call:
> FoldSeek returns 234 hits of which exactly one passes PROB. ≥ 0.7 — the protein itself —
> and its best TM-align score against the other 73 acetyltransferases is 0.347. Removing
> it from this list also removes its row from `acetylases_gwas.tsv`, since it carried no
> SSLBH evidence and was present only by virtue of this hand-curated list.

This PC is the same one hardcoded as `_FOLDSEEK_ACETYLTRANSFERASES` in
`figures/chapter4/lib/figure4_2_panelA.py` and as `gwas_best_at_fastas` in
`figures/chapter4/main.py`, and corresponds to PROTEIN05_GWAS_AC_K111 in S1_Table.xlsx.
All three places encode the same manual result; if it is ever revised, revise it in all
three.
"""

from __future__ import annotations

import pandas as pd

SSLBH_CARRYING = "sslbh-carrying"
BEST_PERFORMING = "best-performing"

# (locus, clustering_level, PC) — the best-performing GWAS predictor per K-locus.
BEST_PERFORMING_PREDICTORS = [
    ("KL111", "PCI80C80", "PC0915"),   # PROTEIN05_GWAS_AC_K111
]

_KEY = ["locus", "clustering_level", "PC"]


def add_prediction_class(sslbh_hits: pd.DataFrame, all_flagged: pd.DataFrame) -> pd.DataFrame:
    """
    Label every GWAS acetyltransferase prediction and add any missing best performer.

    Args:
        sslbh_hits:  the `acetyl_hhsearch == True` subset (the current table contents)
        all_flagged: the full flagged pyseer table (all lasso hits, pvalue_corr ≤ 0.05),
                     used to recover best-performing rows that carry no SSLBH evidence

    Returns:
        sslbh_hits plus any missing best performer, with a `prediction_class` column.
    """
    out = sslbh_hits.copy()
    out["prediction_class"] = SSLBH_CARRYING

    present = set(map(tuple, out[_KEY].itertuples(index=False, name=None)))

    missing_rows = []
    for key in BEST_PERFORMING_PREDICTORS:
        if key in present:
            mask = (out["locus"] == key[0]) & (out["clustering_level"] == key[1]) & (out["PC"] == key[2])
            out.loc[mask, "prediction_class"] = BEST_PERFORMING
            print(f"  [best-performing] {key[0]} {key[1]} {key[2]} — already SSLBH-flagged, relabelled")
            continue

        locus, level, pc = key
        row = all_flagged[
            (all_flagged["locus"] == locus)
            & (all_flagged["clustering_level"] == level)
            & (all_flagged["PC"] == pc)
        ]
        if row.empty:
            raise ValueError(
                f"best-performing predictor {locus} {level} {pc} not found in the flagged "
                f"pyseer table — the GWAS input changed and this hardcoded list is stale."
            )
        row = row.copy()
        row["prediction_class"] = BEST_PERFORMING
        missing_rows.append(row)
        print(f"  [best-performing] {locus} {level} {pc} — no SSLBH evidence, row added from pyseer")

    if missing_rows:
        out = pd.concat([out, *missing_rows], ignore_index=True)

    return out.sort_values(["locus", "clustering_level", "PC"]).reset_index(drop=True)
