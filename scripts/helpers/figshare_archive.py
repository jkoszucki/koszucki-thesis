"""Create a Zenodo/Figshare-ready compressed archive of cfg.input_dir + cfg.output_dir.

Archive layout: everything is nested under an `input/` or `output/` top-level
prefix so both roots can share one deposited record.

Excludes:
    both:       `.DS_Store` (anywhere)
    input_dir:  `archive/`, `data-gwas/1_BACTERIA/`, `supplementary-thesis/thesis-manuscript/`
                (large and/or not needed for the deposited record)
    output_dir: `other/prophage_blastp_db/` (regenerable BLASTP DB, not a result),
                `other/alphafold3/1_DRAWN/` (regenerable pLDDT renders, superseded by
                `2_DRAWN_AND_ORGANISED/`),
                `other/zenodo/` (the archive's own output location — avoids nesting
                a previous archive inside a new one)

Symlinks (e.g. other/alphafold3/2_DRAWN_AND_ORGANISED/*, which links into
1_DRAWN/*) are always dereferenced — the real file content is embedded, never
the symlink itself, so the archive is self-contained regardless of format.

Read-only with respect to input_dir and output_dir — only ever reads from them,
writes the archive to output_dir/other/zenodo/.

Usage:
    conda run -n jkoszucki python scripts/helpers/figshare_archive.py [--format zip|tar.gz] [--out <path>]
"""

import argparse
import sys
import tarfile
import zipfile
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import Config

INPUT_EXCLUDED_DIRS = {
    Path("archive"),
    Path("data-gwas/1_BACTERIA"),
    Path("supplementary-thesis/thesis-manuscript"),
}

OUTPUT_EXCLUDED_DIRS = {
    Path("other/prophage_blastp_db"),
    Path("other/alphafold3/1_DRAWN"),
    Path("other/zenodo"),
}

EXCLUDED_FILENAMES = {".DS_Store"}


def _is_excluded(relative_path: Path, excluded_dirs: set[Path]) -> bool:
    if relative_path.name in EXCLUDED_FILENAMES:
        return True
    return any(
        relative_path == excluded or excluded in relative_path.parents
        for excluded in excluded_dirs
    )


def _iter_included_files(root_dir: Path, excluded_dirs: set[Path], arc_prefix: str):
    for path in sorted(root_dir.rglob("*")):
        if path.is_dir():
            continue
        if not path.exists():
            print(f"  [warn] broken symlink, skipping: {path}")
            continue
        relative_path = path.relative_to(root_dir)
        if _is_excluded(relative_path, excluded_dirs):
            continue
        yield path, Path(arc_prefix) / relative_path


def build_figshare_archive(input_dir: Path, output_dir: Path, out_path: Path, archive_format: str) -> Path:
    files = list(_iter_included_files(input_dir, INPUT_EXCLUDED_DIRS, "input"))
    files += list(_iter_included_files(output_dir, OUTPUT_EXCLUDED_DIRS, "output"))

    if archive_format == "zip":
        with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for path, arcname in files:
                zf.write(path, arcname=arcname)
    elif archive_format == "tar.gz":
        # dereference=True: embed the real file content for symlinks (e.g.
        # other/alphafold3/2_DRAWN_AND_ORGANISED/* → 1_DRAWN/*) instead of
        # storing the symlink itself, which would point at an absolute local
        # path that doesn't exist once the archive is extracted elsewhere.
        with tarfile.open(out_path, "w:gz", dereference=True) as tf:
            for path, arcname in files:
                tf.add(path, arcname=arcname)
    else:
        raise ValueError(f"Unsupported format: {archive_format}")

    return out_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=["zip", "tar.gz"], default="tar.gz")
    parser.add_argument("--out", type=Path, default=None, help="Output archive path")
    args = parser.parse_args()

    cfg = Config()
    input_dir = cfg.input_dir
    output_dir = cfg.output_dir

    if args.out is not None:
        out_path = args.out
    else:
        extension = "zip" if args.format == "zip" else "tar.gz"
        out_path = cfg.output_dir / "other" / "zenodo" / f"submission_archive_{date.today():%Y-%m-%d}.{extension}"

    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(
        f"Archiving {input_dir} + {output_dir} -> {out_path}\n"
        f"  excluding input:  {sorted(str(p) for p in INPUT_EXCLUDED_DIRS)}\n"
        f"  excluding output: {sorted(str(p) for p in OUTPUT_EXCLUDED_DIRS)}"
    )
    build_figshare_archive(input_dir, output_dir, out_path, args.format)

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"Done: {out_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
