from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))

from config import Config
from ktypes_draw import KTypeStructureDrawer
from ktypes_modifications import KTypeModificationsAPI
from ktypes_process import KTypeTablesAPI
from ktypes_similarity import KTypeSimilarityAPI

SAMPLE_DIR = Path(__file__).resolve().parent / "sample"


def main() -> None:
    cfg = Config()
    input_xlsx = cfg.input_dir / "supplementary-thesis" / "supplementary-tables" / "S2_Table.xlsx"
    chapter4_dir = cfg.output_dir / "cps_structures"

    ### BUILD TABLES
    tables_api = KTypeTablesAPI(input_xlsx=input_xlsx, output_dir=chapter4_dir)
    processed_df = tables_api.build_processed_table()
    processed_csv = tables_api.export_processed_table(processed_df)

    modifications_api = KTypeModificationsAPI(input_xlsx=input_xlsx, output_dir=chapter4_dir)
    modifications_df = modifications_api.build_modifications_table()
    modifications_csv = modifications_api.export_modifications_table(modifications_df)

    similarity_api = KTypeSimilarityAPI(input_xlsx=input_xlsx, output_dir=chapter4_dir)
    similarity_df = similarity_api.build_similarity_table(processed_df)
    similarity_api.export_similarity_table(df=similarity_df)

    ### DRAW STRUCTURES (checkpoint: skip already-generated)
    drawer = KTypeStructureDrawer(
        output_dir=cfg.output_dir / "other",
        processed_csv=processed_csv,
        modifications_csv=modifications_csv,
        plots_subdir="cps_drawn",
    )
    ktypes_df = pd.read_csv(processed_csv)
    skipped = []
    for i, row in ktypes_df.iterrows():
        kname = str(row.get("structure_id", "") or row.get("K_type", "")).strip()
        stem = (kname if kname else f"row_{i}").replace("/", "_").replace(" ", "_")
        if (drawer.png_dir / f"{stem}.png").exists() and (drawer.svg_dir / f"{stem}.svg").exists():
            skipped.append(kname)
        else:
            drawer.draw_single(row, stem_name=stem)
    if skipped:
        print(f"Skipped {len(skipped)} already-generated structures: {', '.join(skipped)}")

    ### SAMPLES
    if not processed_df.empty:
        sample_size = min(10, len(processed_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        processed_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes.csv", index=False
        )

    if not modifications_df.empty:
        sample_size = min(10, len(modifications_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        modifications_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes_modifications.csv", index=False
        )

    if not similarity_df.empty:
        sample_size = min(10, len(similarity_df))
        SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
        similarity_df.sample(n=sample_size, random_state=0).to_csv(
            SAMPLE_DIR / "ktypes_sim.csv", index=False
        )


if __name__ == "__main__":
    main()
