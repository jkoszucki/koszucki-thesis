"""Tests for the path-fingerprint similarity used in chapter 3.

Run:  conda run -n jkoszucki python scripts/processing/cps-proc/tests/test_fingerprint.py
(also collectable by pytest if installed).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))

from ktypes_process import (  # noqa: E402
    _branch_path_fingerprint,
    _path_fingerprint,
    _reverse_bond,
    _reverse_path,
    validate_structure,
)


def jaccard(a, b):
    return len(a & b) / len(a | b)


def test_reverse_bond():
    assert _reverse_bond("a1-3") == "a3-1"
    assert _reverse_bond("b1-4") == "b4-1"
    assert _reverse_bond("b1,2") == "b1,2'"  # malformed never collides with a forward token


def test_reverse_path_reverses_bonds():
    assert _reverse_path(("Glc", "b1-4", "Man")) == ("Man", "b4-1", "Glc")
    assert _reverse_path(("Glc",)) == ("Glc",)


def test_ring_rotation_invariance():
    m, b = ["Man", "Glc", "Glc"], ["b1-4", "a1-3", "b1-4"]
    fp = _path_fingerprint(m, b, circular=True)
    rot = _path_fingerprint(m[1:] + m[:1], b[1:] + b[:1], circular=True)
    assert fp == rot


def test_direction_aware_no_spurious_match():
    # A(b1-4)B forwards must not match B(b1-4)A forwards.
    x = _path_fingerprint(["GlcA", "Fuc"], ["b1-4", "a1-3"], circular=True)
    y = _path_fingerprint(["Fuc", "GlcA"], ["b1-4", "a1-3"], circular=True)
    assert ("Fuc", "b1-4", "GlcA") not in x
    assert ("Fuc", "b1-4", "GlcA") in y
    # single-residue paths are shared, two-residue paths are not
    assert {p for p in x & y if len(p) > 1} == set()


def test_ring_closure_captured():
    fp = _path_fingerprint(["A", "B", "C"], ["x1-1", "y1-1", "z1-1"], circular=True)
    assert ("C", "z1-1", "A") in fp  # crosses the closing bond


def test_branch_anchor_included():
    fp = _branch_path_fingerprint("GlcA(α1-3)Man(P1)")
    assert fp == frozenset({("GlcA",), ("Man",), ("GlcA", "a1-3", "Man"), ("Man", "a3-1", "GlcA")})


def test_k2_k13_worked_example():
    # K2 branch: GlcA(a1-3)Man ; K13 branch: Gal(b1-4)GlcA(a1-3)Man  -> 4 shared of 9 = 0.444
    k2 = _branch_path_fingerprint("GlcA(α1-3)Man(P1)")
    k13 = _branch_path_fingerprint("Gal(β1-4)GlcA(α1-3)Man(P1)")
    assert round(jaccard(k2, k13), 3) == 0.444


def test_identical_scores_one():
    core = _path_fingerprint(["Gal", "Glc"], ["b1-3", "b1-4"], circular=True)
    assert jaccard(core, core) == 1.0


def test_validate_clean():
    assert validate_structure("X", "Man(β1-4)Glc(α1-3)Glc(β1-4)-OUT", "GlcA(α1-3)Man(P1)") == []
    assert validate_structure("X", "Man(β1-4)Glc(α1-3)-OUT", None) == []


def test_validate_catches_known_errors():
    assert any("malformed core bond" in p for p in validate_structure("K3", "GalA(α1-3)Gal(β1,2)-OUT", None))
    assert any("lacks a (P#) anchor" in p for p in validate_structure("K34", "Rha(α1-2)GalA(α1-2)-OUT", "Rha(α1-4)GalA(α1-4)"))
    assert any("anchor says" in p for p in validate_structure("K12", "Gal(α1-2)Galf(β1-6)-OUT", "GlcA(β1-3)Gal(P2)"))
    assert any("malformed branch bond" in p for p in validate_structure("K18", "Rha(α1-3)Glc(α1-3)-OUT", "GlcA(β1-2)Rha(α2-3)Glc(P2)"))
    assert any("unknown core residue" in p for p in validate_structure("X", "Xyl(α1-3)Glc(α1-3)-OUT", None))
    assert any("does not end with '-OUT'" in p for p in validate_structure("X", "Glc(α1-3)Glc(α1-3)", None))


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print("ok ", t.__name__)
    print(f"{len(tests)} tests passed")
