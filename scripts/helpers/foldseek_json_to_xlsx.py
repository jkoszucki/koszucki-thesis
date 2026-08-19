"""
Convert a FoldSeek Search JSON download into the 21-column .xlsx layout used by S3_Data.

The 69 K-locus acetyltransferase result files in
`input_dir/supplementary-thesis/supplementary-data/S3_Data/` were exported from the
FoldSeek web server's results table. This script reproduces that exact table from the
raw JSON the server also offers, so results retrieved as JSON (e.g. for the experimental
acetyltransferases) can be filed in the same format and read by
`figures/chapter4/lib/figure4_2_panelB.py` without special-casing.

    conda run -n jkoszucki python scripts/helpers/foldseek_json_to_xlsx.py \
        --json ~/Downloads/Foldseek_2026_07_26_16_33_18.json \
        --out  <somewhere>/PROTEIN02_MOD_AC_K2.xlsx

**The search must have been run with alignment mode = TM-align.** In that mode the
server reports a TM-score where the default 3Di+AA mode reports an e-value, and carries
it in the JSON's `eval` field (0-1, and the sort key of the hit list). Panel B filters
on `Tmscore`, so a 3Di+AA JSON has nothing to filter on and is rejected outright rather
than silently converted into a file whose Tmscore column holds e-values.

A structure uploaded as a multimer is searched chain by chain, giving one job per chain
in a single download; `--job` selects which to keep (default `job_A`), matching the
single-chain searches behind the 69 reference files.

The protein's identity is not in the JSON -- the server names every query `job_A` -- so
it lives in the output filename, exactly as it does for the existing 69 files. `--query`
overwrites the QUERY column if you want it in the table as well.

Note the TM-score is rounded to 3 decimals in the JSON against the table export's 4;
irrelevant except for a hit sitting exactly on a cutoff.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Column order and dtype of the S3_Data tables, in order.
_COLUMNS: list[tuple[str, str, type]] = [
    # (xlsx column,  JSON field,     dtype)
    ("QUERY",        "query",        str),
    ("TARGET",       "",             str),   # target + " " + description
    ("SEQ. ID.",     "seqId",        float),
    ("ALENLEN",      "alnLength",    int),
    ("MISMATCH",     "missmatches",  int),
    ("GAPOPEN",      "gapsopened",   int),
    ("QSTART",       "qStartPos",    int),
    ("QEND",         "qEndPos",      int),
    ("TSTART",       "dbStartPos",   int),
    ("TEND",         "dbEndPos",     int),
    ("PROB.",        "prob",         float),
    ("Tmscore",      "eval",         float),  # TM-align mode: TM-score, not an e-value
    ("SCORE",        "score",        int),
    ("QLEN",         "qLen",         int),
    ("TLEN",         "dbLen",        int),
    ("QSEQ",         "qAln",         str),
    ("TSEQ",         "dbAln",        str),
    ("COORDS",       "tCa",          str),
    ("TSEQ_FULL",    "tSeq",         str),
    ("taxid",        "taxId",        int),
    ("organism",     "taxName",      str),
]


def _select_job(payload, job: str):
    """Pick one search job from a FoldSeek JSON payload.

    A structure uploaded as a multimer is searched chain by chain, so the download holds
    one job per chain (`job_A`, `job_B`, …) — same query sequence, different hit lists.
    The 69 K-locus reference searches are single-chain, so the convention here is to keep
    chain A and discard the rest rather than invent a way to merge them.
    """
    if not isinstance(payload, list):
        return payload
    if len(payload) == 1:
        return payload[0]

    headers = [b["queries"][0]["header"] for b in payload]
    for block in payload:
        if block["queries"][0]["header"] == job:
            print(f"  multimer search: {len(payload)} chains ({', '.join(headers)}) — keeping {job}")
            return block
    raise ValueError(
        f"{len(payload)} search jobs ({', '.join(headers)}) and none is named '{job}'"
    )


def _alignments(payload, job: str = "job_A") -> tuple[list[dict], str, str]:
    """Return (hits, alignment mode, database) from a FoldSeek Search JSON payload."""
    payload = _select_job(payload, job)

    results = payload["results"]
    if len(results) != 1:
        raise ValueError(f"expected a single target database, found {len(results)}")
    block = results[0]

    # `alignments` is a dict keyed by stringified rank, each value a list of hits.
    hits = [hit for group in block["alignments"].values() for hit in group]
    return hits, payload.get("mode", ""), block.get("db", "")


def _target(hit: dict) -> str:
    description = (hit.get("description") or "").strip()
    return f"{hit['target']} {description}".strip()


def convert(json_path: Path, out_path: Path, query_name: str | None = None,
            job: str = "job_A") -> pd.DataFrame:
    payload = json.loads(json_path.read_text())
    hits, mode, db = _alignments(payload, job)

    if mode != "tmalign":
        raise ValueError(
            f"alignment mode is '{mode}', not 'tmalign' — this search has no TM-score "
            f"(its 'eval' field is a real e-value). Re-run on the FoldSeek server with "
            f"alignment mode set to TM-align."
        )

    print(f"  {json_path.name}: {len(hits)} hits, mode={mode}, db={db}")

    frame = {}
    for column, field, dtype in _COLUMNS:
        values = [_target(h) for h in hits] if column == "TARGET" else [h[field] for h in hits]
        frame[column] = pd.Series(values, dtype=object).astype(dtype)

    df = pd.DataFrame(frame)
    if query_name:
        df["QUERY"] = query_name

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(out_path, index=False)
    print(f"  Tmscore {df['Tmscore'].min():.3f}-{df['Tmscore'].max():.3f}, "
          f"PROB. {df['PROB.'].min():.2f}-{df['PROB.'].max():.2f}")
    print(f"  → {out_path}")
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, type=Path, help="FoldSeek Search JSON download")
    ap.add_argument("--out", required=True, type=Path,
                    help="output .xlsx; its stem identifies the protein")
    ap.add_argument("--query", default=None,
                    help="overwrite the QUERY column (the server writes 'job_A')")
    ap.add_argument("--job", default="job_A",
                    help="which chain to keep from a multimer search (default: job_A)")
    args = ap.parse_args()

    try:
        convert(args.json, args.out, args.query, args.job)
    except (ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
