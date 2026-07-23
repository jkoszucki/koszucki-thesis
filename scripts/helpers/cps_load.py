"""
Quick script to load cps.xlsx and produce a flat table with a human-readable
unique identifier column for manual refinement.

Output: cps_ids.csv  (written next to this script)

identifier format:
  {structure_id}_{K-type}              — when strain_lit is empty
  {structure_id}_{K-type}_(strain_lit) — when strain_lit is present
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts" / "helpers"))
from config import Config

cfg = Config()

df = pd.read_excel(cfg.input_dir / "cps.xlsx", sheet_name="ktypes")

# Build identifier
def make_id(row):
    parts = [str(row["structure_id"]), str(row["K-type"])]
    if pd.notna(row["strain_lit"]) and str(row["strain_lit"]).strip():
        parts.append(str(row["strain_lit"]))
    return "_".join(parts)

df.insert(0, "identifier", df.apply(make_id, axis=1))

out_path = Path(__file__).resolve().parent / "cps_ids.csv"
df.to_csv(out_path, index=False)
print(f"Written {len(df)} rows → {out_path}")
print()
print(df[["identifier", "structure_id", "K-type", "strain_lit"]].to_string(index=False))
