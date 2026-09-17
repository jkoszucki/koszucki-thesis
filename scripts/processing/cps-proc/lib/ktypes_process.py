"""API for preparing K-type analytical tables."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from ktypes_base import BaseKTypeAPI

ALLOWED = {"Glc", "GlcA", "Man", "Rha", "Fuc", "Gal", "GalA", "Galf", "Sug"}
ANCHOR_POS_RE = re.compile(r"@P(\d+)")
BOND_RE = re.compile(r"^[ab]1-[1-6]$")          # normalised bond token, e.g. a1-3, b1-4
BRANCH_ANCHOR_RE = re.compile(r"\(P(\d+)\)$")   # trailing (P#) on a branch path


class StructureFormatError(ValueError):
    """Raised when a repeating-unit string in the input workbook is malformed."""


def _reverse_bond(bond: str) -> str:
    """Reverse the direction of a glycosidic bond token: a1-3 -> a3-1.

    A bond token reads donor-anomeric-position -> acceptor-position. When a
    path is read backwards the acceptor becomes the first residue, so the
    positions swap. Tokens that do not match the standard form are returned
    with a trailing apostrophe so they can never collide with a forward token.
    """
    m = re.fullmatch(r"([ab])(\d)-(\d)", bond)
    return f"{m.group(1)}{m.group(3)}-{m.group(2)}" if m else bond + "'"


def validate_structure(structure_id: str, backbone: str, branch) -> List[str]:
    """Return a list of format problems for one repeating unit (empty if clean).

    Checks, on the normalised strings:
      * backbone ends with -OUT and every segment is Residue(bond)
      * every residue is in ALLOWED, every bond matches BOND_RE
      * the core has as many bonds as residues (it is a ring)
      * every branch path ends with a (P#) anchor whose residue is the one
        actually at core position #, and every branch bond matches BOND_RE
    """
    problems: List[str] = []
    bb = str(backbone)
    if not bb.strip().endswith("-OUT"):
        problems.append(f"{structure_id}: backbone does not end with '-OUT': {bb!r}")
    core_monos, core_bonds = [], []
    for part in (x for x in re.sub(r"-?OUT$", "", _norm_text(bb)).split(")") if x):
        mono, sep, bond = part.partition("(")
        if not sep or not mono or not bond:
            problems.append(f"{structure_id}: backbone segment is not Residue(bond): {part!r}")
            continue
        if mono not in ALLOWED:
            problems.append(f"{structure_id}: unknown core residue {mono!r}")
        if not BOND_RE.match(bond):
            problems.append(f"{structure_id}: malformed core bond {bond!r} (expected e.g. a1-3)")
        core_monos.append(mono)
        core_bonds.append(bond)
    if core_monos and len(core_bonds) != len(core_monos):
        problems.append(f"{structure_id}: core has {len(core_monos)} residues but {len(core_bonds)} bonds")

    if _is_missing_branch_text(branch):
        return problems
    for path in str(branch).split(";"):
        ps = _norm_text(path)
        if not ps:
            problems.append(f"{structure_id}: empty branch path (stray ';')")
            continue
        m = BRANCH_ANCHOR_RE.search(ps)
        if not m:
            problems.append(f"{structure_id}: branch path lacks a (P#) anchor: {path.strip()!r}")
            continue
        segs = []
        for part in (x for x in ps.split(")") if x):
            mono, sep, right = part.partition("(")
            if not sep or not mono or not right:
                problems.append(f"{structure_id}: branch segment is not Residue(bond): {part!r}")
                continue
            segs.append((mono, right))
        if not segs:
            continue
        pos = int(m.group(1))
        anchor = segs[-1][0]
        if pos < 1 or pos > len(core_monos):
            problems.append(f"{structure_id}: anchor P{pos} out of range, core has {len(core_monos)} residues: {path.strip()!r}")
        elif core_monos[pos - 1] != anchor:
            problems.append(f"{structure_id}: anchor says {anchor}(P{pos}) but core position {pos} is {core_monos[pos - 1]}: {path.strip()!r}")
        for mono, bond in segs[:-1]:
            if mono not in ALLOWED:
                problems.append(f"{structure_id}: unknown branch residue {mono!r}")
            if not BOND_RE.match(bond):
                problems.append(f"{structure_id}: malformed branch bond {bond!r} in {path.strip()!r}")
    return problems


def _pick_col(df: pd.DataFrame, candidates_lower: Sequence[str]) -> Optional[str]:
    for column in df.columns:
        lowered = str(column).replace("\ufeff", "").strip().lower()
        if any(token in lowered for token in candidates_lower):
            return column
    return None


def _norm_text(text: str) -> str:
    t = str(text) if text is not None else ""
    t = t.replace("α", "a").replace("β", "b").replace("ALPHA", "a").replace("BETA", "b")
    t = t.replace("→", "-").replace("–", "-").replace("—", "-").replace("->", "-")
    t = re.sub(r"\s+", "", t)
    return t


def _canon_mono(mono: str) -> str:
    return mono if mono in ALLOWED else "Other"


def _parse_core_segments(backbone: str) -> tuple[list[str], list[str]]:
    s = _norm_text(backbone)
    s = re.sub(r"-?OUT$", "", s)
    monos, bonds = [], []
    for part in (p for p in s.split(")") if p):
        mono, sep, bond = part.partition("(")
        mono, bond = mono.strip(), bond.strip()
        if not sep or not mono or not bond:
            continue
        monos.append(mono)
        bonds.append(bond)
    return monos, bonds


def _assemble_core_fields(monos: Sequence[str], bonds: Sequence[str]) -> Dict[str, str]:
    if not monos:
        return {
            "core_positions": "",
            "core_bonds_all": "",
            "core_bonds_unique": "",
            "core_edges": "",
            "core_pairs": "",
            "core_monos_all": "",
            "core_monos_unique": "",
        }

    canon_monos = [_canon_mono(m) for m in monos]
    core_positions = [f"P{i + 1}={mono}" for i, mono in enumerate(canon_monos)]

    edges, pairs = [], set()
    for idx in range(len(canon_monos) - 1):
        left = canon_monos[idx]
        right = canon_monos[idx + 1]
        bond = bonds[idx] if idx < len(bonds) else ""
        edges.append(f"{left}({bond})->{right}")
        pairs.add("~~".join(sorted([left, right])))

    bonds_all = list(bonds)
    bonds_unique = sorted(set(bonds_all))

    seen = set()
    unique_monos: List[str] = []
    for mono in canon_monos:
        if mono not in seen:
            seen.add(mono)
            unique_monos.append(mono)

    return {
        "core_positions": "; ".join(core_positions),
        "core_bonds_all": "; ".join(bonds_all),
        "core_bonds_unique": "; ".join(bonds_unique),
        "core_edges": "; ".join(edges),
        "core_pairs": "; ".join(sorted(pairs)),
        "core_monos_all": "; ".join(canon_monos),
        "core_monos_unique": "; ".join(unique_monos),
    }


def _reverse_path(path: tuple) -> tuple:
    """Read a (mono, bond, mono, ...) path backwards, reversing each bond token."""
    rev = []
    for i, elem in enumerate(reversed(path)):
        # In the reversed tuple, odd indices hold bond tokens.
        rev.append(_reverse_bond(elem) if i % 2 == 1 else elem)
    return tuple(rev)


def _path_fingerprint(monos: list[str], bonds: list[str], circular: bool = False) -> frozenset[tuple]:
    """All contiguous subpaths of a residue/bond sequence, in both directions.

    A path is the tuple (m1, b1, m2, b2, ..., mk). Each path is added together
    with its reverse; the reverse carries reversed bond tokens (a1-3 -> a3-1),
    because a glycosidic bond is directional and a path read backwards must
    not collide with a genuinely opposite linkage read forwards.

    If circular=True and len(bonds) == len(monos), the sequence is a ring of n
    residues joined by n bonds (the n-th closing the ring). Subpaths crossing
    the ring closure are captured by doubling the sequence to length 2n and
    taking every window of length 1..n from each of the n start positions.
    """
    paths: set[tuple] = set()
    n = len(monos)
    if n == 0:
        return frozenset()

    if circular and len(bonds) == n:
        monos_ext, bonds_ext = list(monos) * 2, list(bonds) * 2
        windows = ((start, start + length) for start in range(n) for length in range(1, n + 1))
    else:
        monos_ext, bonds_ext = list(monos), list(bonds)
        windows = ((start, end) for start in range(n) for end in range(start + 1, n + 1))

    for start, end in windows:
        seg_m = monos_ext[start:end]
        seg_b = bonds_ext[start:end - 1]
        elems: list = []
        for i, m in enumerate(seg_m):
            elems.append(m)
            if i < len(seg_b):
                elems.append(seg_b[i])
        fwd = tuple(elems)
        paths.add(fwd)
        paths.add(_reverse_path(fwd))
    return frozenset(paths)


def _branch_path_fingerprint(branch_text) -> frozenset[tuple]:
    """Path fingerprint from a branch string (';'-separated paths)."""
    if _is_missing_branch_text(branch_text):
        return frozenset()
    all_paths: set[tuple] = set()
    for path_str in (p for p in str(branch_text).split(";") if p.strip()):
        s = _norm_text(path_str)
        monos, bonds = [], []
        for part in (p for p in s.split(")") if p):
            mono, sep, right = part.partition("(")
            mono, right = mono.strip(), right.strip()
            if not sep or not mono or not right:
                continue
            if re.match(r"P\d+", right, re.IGNORECASE):
                monos.append(_canon_mono(mono))
                break
            monos.append(_canon_mono(mono))
            bonds.append(right)
        if len(monos) >= 1:
            all_paths |= set(_path_fingerprint(monos, bonds))
    return frozenset(all_paths)


def _rotate_core_by_anchor(
    monos: Sequence[str], bonds: Sequence[str], anchor_positions: Sequence[int]
) -> tuple[list[str], list[str], Dict[int, int]]:
    count = len(monos)
    if count == 0:
        return list(monos), list(bonds), {}

    valid_positions = [pos for pos in anchor_positions if 1 <= pos <= count]
    if not valid_positions:
        mapping = {i + 1: i + 1 for i in range(count)}
        return list(monos), list(bonds), mapping

    anchor = min(valid_positions)
    anchor_idx = anchor - 1
    order = list(range(count))
    rotated_indices = order[anchor_idx:] + order[:anchor_idx]

    rotated_monos = [monos[i] for i in rotated_indices]
    rotated_bonds = [bonds[i] for i in rotated_indices] if bonds else []
    mapping = {old_idx + 1: new_idx + 1 for new_idx, old_idx in enumerate(rotated_indices)}
    return rotated_monos, rotated_bonds, mapping


def _extract_anchor_positions(anchor_sites_text: str) -> list[int]:
    if not isinstance(anchor_sites_text, str) or not anchor_sites_text.strip():
        return []
    return [int(match.group(1)) for match in ANCHOR_POS_RE.finditer(anchor_sites_text)]


def _remap_anchor_positions(text: str, mapping: Dict[int, int]) -> str:
    if not isinstance(text, str) or not text:
        return text

    def repl(match: re.Match[str]) -> str:
        old_pos = int(match.group(1))
        new_pos = mapping.get(old_pos, old_pos)
        return f"@P{new_pos}"

    return ANCHOR_POS_RE.sub(repl, text)


def _remap_branch_anchor_positions(branch_dict: Dict[str, str], mapping: Dict[int, int]) -> Dict[str, str]:
    if not mapping:
        return dict(branch_dict)
    updated = dict(branch_dict)
    updated["branch_anchor_sites"] = _remap_anchor_positions(
        updated.get("branch_anchor_sites", ""), mapping
    )
    updated["branch_anchor_edges"] = _remap_anchor_positions(
        updated.get("branch_anchor_edges", ""), mapping
    )
    return updated


def _is_missing_branch_text(branch_text) -> bool:
    if branch_text is None:
        return True
    if isinstance(branch_text, float):
        try:
            if pd.isna(branch_text):
                return True
        except Exception:
            pass
    text = str(branch_text).strip()
    if not text:
        return True
    lowered = text.lower()
    if lowered in {"nan", "na", "n/a", "none", "null"}:
        return True
    return False


def _parse_branch_simple(branch_text: str) -> Dict[str, str | int]:
    if _is_missing_branch_text(branch_text):
        return {
            "n_branch": 0,
            "branch_monos_all": "",
            "branch_monos_unique": "",
            "branch_edges": "",
            "branch_pairs": "",
            "branch_bonds_all": "",
            "branch_bonds_unique": "",
            "branch_anchor_sites": "",
            "branch_anchor_edges": "",
        }

    n_paths = 0
    monos_all: List[str] = []
    monos_set = set()
    edges_set = set()
    pairs_set = set()
    bonds_all: List[str] = []
    bonds_set = set()
    anchor_sites = set()
    anchor_edges = set()

    for path in (p for p in str(branch_text).split(";") if p.strip()):
        s = _norm_text(path)
        segments = []
        for part in (p for p in s.split(")") if p):
            mono, sep, right = part.partition("(")
            mono, right = mono.strip(), right.strip()
            if not sep or not mono or not right:
                continue
            segments.append((mono, right))
        if len(segments) < 2:
            continue

        n_paths += 1
        anchor_mono, anchor_token = segments[-1]
        anchor_sites.add(f"{anchor_mono}@{anchor_token}")
        pre_segments = segments[:-1]

        for mono, bond in pre_segments:
            canon = _canon_mono(mono)
            monos_all.append(canon)
            monos_set.add(canon)
            bonds_all.append(bond)
            bonds_set.add(bond)

        for idx in range(len(pre_segments) - 1):
            left, bond = pre_segments[idx]
            right, _ = pre_segments[idx + 1]
            left_c, right_c = _canon_mono(left), _canon_mono(right)
            edges_set.add(f"{left_c}({bond})->{right_c}")
            pairs_set.add("~~".join(sorted([left_c, right_c])))

        left, bond = pre_segments[-1]
        left_c, right_c = _canon_mono(left), _canon_mono(anchor_mono)
        edges_set.add(f"{left_c}({bond})->{right_c}")
        pairs_set.add("~~".join(sorted([left_c, right_c])))
        anchor_edges.add(f"{left_c}({bond})->{right_c}@{anchor_token}")

    seen = set()
    unique_monos: List[str] = []
    for mono in monos_all:
        if mono not in seen:
            seen.add(mono)
            unique_monos.append(mono)

    return {
        "n_branch": n_paths,
        "branch_monos_all": "; ".join(monos_all),
        "branch_monos_unique": "; ".join(unique_monos),
        "branch_edges": "; ".join(sorted(edges_set)),
        "branch_pairs": "; ".join(sorted(pairs_set)),
        "branch_bonds_all": "; ".join(bonds_all),
        "branch_bonds_unique": "; ".join(sorted(bonds_set)),
        "branch_anchor_sites": "; ".join(sorted(anchor_sites)),
        "branch_anchor_edges": "; ".join(sorted(anchor_edges)),
    }


class KTypeTablesAPI(BaseKTypeAPI):
    """High-level API for preparing K-type analysis tables."""

    KTYPE_SHEET = "capsule_structures"
    MODIFICATIONS_SHEET = "capsule_modifications"
    PROCESSED_FILENAME = "ktypes.csv"
    SIMILARITY_FILENAME = "ktypes_sim.csv"

    def get_raw_ktypes(self) -> pd.DataFrame:
        """Return the raw K-type worksheet."""
        return self.get_sheet(self.KTYPE_SHEET)

    def get_modifications(self) -> pd.DataFrame:
        """Return the modifications worksheet."""
        return self.get_sheet(self.MODIFICATIONS_SHEET)

    def process_dataframe(self, df_raw: pd.DataFrame) -> pd.DataFrame:
        k_col = _pick_col(df_raw, ["k-type", "k_type", "ktype", "serotype", "k "])
        backbone_col = _pick_col(df_raw, ["backbone_linkages_coarse", "backbone", "linkage"])
        branches_col = _pick_col(df_raw, ["branch", "branches"])
        oac_col = _pick_col(df_raw, ["oac"])
        opy_col = _pick_col(df_raw, ["opy", "o py", "o_py"])

        df = pd.DataFrame()
        df["K_type"] = df_raw[k_col].astype(str).str.replace(r"\s+", "", regex=True)
        df["_backbone"] = df_raw[backbone_col].astype(str)
        df["_branch"] = df_raw[branches_col].astype(str) if branches_col else ""
        df["_OAc"] = df_raw[oac_col].astype(str) if oac_col else ""
        df["_Opy"] = df_raw[opy_col].astype(str) if opy_col else ""

        problems: List[str] = []
        for idx, record in df.iterrows():
            sid = str(df_raw.loc[idx, "structure_id"]) if "structure_id" in df_raw.columns else record["K_type"]
            problems += validate_structure(sid, record["_backbone"], df_raw.loc[idx, branches_col] if branches_col else None)
        if problems:
            raise StructureFormatError(
                f"{len(problems)} problem(s) in the input structures — fix the workbook, "
                "nothing was written:\n  " + "\n  ".join(problems)
            )

        rows = []
        for idx, record in df.iterrows():
            branch_parsed = _parse_branch_simple(record["_branch"])
            anchor_positions = _extract_anchor_positions(
                branch_parsed.get("branch_anchor_sites", "")
            )

            core_parsed, position_map = self._parse_core_simple(
                record["_backbone"], anchor_positions=anchor_positions
            )
            branch_adjusted = _remap_branch_anchor_positions(branch_parsed, position_map)

            row_dict: Dict[str, object] = {
                "K_type": record["K_type"],
            }
            if backbone_col:
                row_dict[backbone_col] = df_raw.loc[idx, backbone_col]
            if branches_col:
                row_dict[branches_col] = df_raw.loc[idx, branches_col]
            if oac_col:
                row_dict[oac_col] = df_raw.loc[idx, oac_col]
            if opy_col:
                row_dict[opy_col] = df_raw.loc[idx, opy_col]

            for extra_col in ["structure_id", "Source", "mw_struct_kda", "Year", "DOI"]:
                if extra_col in df_raw.columns:
                    row_dict[extra_col] = df_raw.loc[idx, extra_col]

            row_dict.update(core_parsed)
            row_dict.update(branch_adjusted)

            combined_unique: List[str] = []
            seen_combined: set[str] = set()
            for key in ("core_monos_unique", "branch_monos_unique"):
                tokens = str(row_dict.get(key, "")).split(";")
                for token in tokens:
                    mono = token.strip()
                    if not mono or mono in seen_combined:
                        continue
                    seen_combined.add(mono)
                    combined_unique.append(mono)

            row_dict["monos_unique_all"] = "; ".join(combined_unique)
            rows.append(row_dict)

        out_df = pd.DataFrame(rows)

        orig_cols = [c for c in ["K_type", backbone_col, branches_col, oac_col, opy_col] if c]
        for extra_col in ["Source", "mw_struct_kda", "Year", "DOI"]:
            if extra_col in out_df.columns and extra_col not in orig_cols:
                orig_cols.insert(1, extra_col)
        if "structure_id" in out_df.columns and "structure_id" not in orig_cols:
            orig_cols.insert(0, "structure_id")

        core_cols = sorted(c for c in out_df.columns if c.startswith("core_"))
        branch_cols = sorted(
            c for c in out_df.columns if c.startswith("branch_") or c == "n_branch"
        )
        ordered_cols = orig_cols + core_cols + branch_cols
        leftover = [c for c in out_df.columns if c not in ordered_cols]

        # Molecular-weight category: optional, only when the workbook carries the column.
        if "mw_struct_kda" in out_df.columns:
            mw = pd.to_numeric(out_df["mw_struct_kda"], errors="coerce")
            out_df["mw_struct_kda_cat"] = "unassigned"
            out_df.loc[mw < 0.7, "mw_struct_kda_cat"] = "low"
            out_df.loc[(mw >= 0.7) & (mw < 1), "mw_struct_kda_cat"] = "medium"
            out_df.loc[mw >= 1, "mw_struct_kda_cat"] = "high"
            leftover.append("mw_struct_kda_cat")

        return out_df[ordered_cols + leftover]

    def _parse_core_simple(
        self, backbone: str, anchor_positions: Optional[Sequence[int]] = None
    ) -> tuple[Dict[str, str], Dict[int, int]]:
        monos, bonds = _parse_core_segments(backbone)
        mapping = {i + 1: i + 1 for i in range(len(monos))}
        core_data = _assemble_core_fields(monos, bonds)
        return core_data, mapping

    def build_processed_table(self, df_raw: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        raw = df_raw if df_raw is not None else self.get_raw_ktypes()
        return self.process_dataframe(raw)

    def export_processed_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        processed = df if df is not None else self.build_processed_table()
        return self.save_dataframe(processed, filename or self.PROCESSED_FILENAME)

    def compute_similarity(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """All-vs-all path-based Jaccard similarity between repeating units.

        For every unordered pair of structures the score is the Jaccard
        coefficient of their path fingerprints (see ``_path_fingerprint``),
        computed separately for the cyclic core, the branches, and the union
        of the two.  The branch score is NaN when neither structure carries a
        branch.  Structures are keyed by ``structure_id`` so that serotypes
        with two resolved structures are scored as separate entries.  Only
        the identifiers and the three scores are written; the fingerprint
        sets themselves are intermediates and are not exported.
        """
        processed = df if df is not None else self.build_processed_table()
        required = ["K_type", "core_monos_all", "core_bonds_all"]
        missing = [col for col in required if col not in processed.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        def split(value: object) -> list[str]:
            if pd.isna(value) or str(value).strip() == "":
                return []
            return [token for token in str(value).split("; ") if token]

        def jaccard(a: frozenset, b: frozenset) -> float:
            if not a and not b:
                return np.nan
            if not a or not b:
                return 0.0
            return len(a & b) / len(a | b)

        id_col = "structure_id" if "structure_id" in processed.columns else "K_type"
        branch_col = next((c for c in processed.columns if c.lower() in {"branch", "branches"}), None)

        ids: list[str] = []
        id_to_ktype: Dict[str, str] = {}
        core_fps: Dict[str, frozenset] = {}
        branch_fps: Dict[str, frozenset] = {}
        for _, row in processed.iterrows():
            sid = str(row[id_col])
            ids.append(sid)
            id_to_ktype[sid] = str(row["K_type"])
            core_fps[sid] = _path_fingerprint(
                split(row["core_monos_all"]), split(row["core_bonds_all"]), circular=True
            )
            branch_fps[sid] = _branch_path_fingerprint(row[branch_col] if branch_col else "")
        total_fps = {sid: core_fps[sid] | branch_fps[sid] for sid in ids}

        rows = []
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a, b = ids[i], ids[j]
                pj_core = jaccard(core_fps[a], core_fps[b])
                pj_branch = jaccard(branch_fps[a], branch_fps[b])
                pj_total = jaccard(total_fps[a], total_fps[b])
                rows.append(
                    {
                        "structure_id_1": a,
                        "structure_id_2": b,
                        "K_type_1": id_to_ktype[a],
                        "K_type_2": id_to_ktype[b],
                        "path_jaccard_core": round(float(pj_core), 3),
                        "path_jaccard_branch": None if pd.isna(pj_branch) else round(float(pj_branch), 3),
                        "path_jaccard_total": round(float(pj_total), 3),
                    }
                )
        return pd.DataFrame(rows)

    def export_similarity_table(
        self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None
    ) -> Path:
        similarity = df if df is not None else self.compute_similarity()
        return self.save_dataframe(similarity, filename or self.SIMILARITY_FILENAME)
