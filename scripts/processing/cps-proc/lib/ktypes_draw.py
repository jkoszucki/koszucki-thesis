"""API for rendering K-type structures as static plots."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib.patches import Circle, RegularPolygon


FIG_SIZE = (12, 8)
ALLOWED_BRANCH_MONOS = {"Glc", "GlcA", "Man", "Rha", "Fru", "Gal", "GalA", "Galf", "Sug"}
_ANCHOR_TOKEN_RE = re.compile(r"^P(\d+)", re.IGNORECASE)

BRANCH_EDGE_LENGTH = 1.0
SIBILINGS_SPREAD = 0.6

def _rotate_vector(vec: tuple[float, float], angle_radians: float) -> tuple[float, float]:
    """Rotate a 2D vector by the given angle (in radians)."""
    x, y = vec
    cos_a = math.cos(angle_radians)
    sin_a = math.sin(angle_radians)
    return (x * cos_a - y * sin_a, x * sin_a + y * cos_a)


def parse_core_positions(val: str):
    if not isinstance(val, str) or not val.strip():
        return []
    items = [p.strip() for p in val.split(";") if p.strip()]
    out = []
    for item in items:
        match = re.match(r"P(\d+)\s*[:=]?\s*([A-Za-z0-9]+)", item)
        if match:
            out.append((int(match.group(1)), match.group(2)))
        else:
            tokens = item.replace("=", " ").replace(":", " ").split()
            if len(tokens) >= 2 and tokens[0].startswith("P"):
                try:
                    idx = int(re.sub(r"\D", "", tokens[0]))
                    out.append((idx, tokens[-1]))
                except Exception:
                    pass
    out.sort(key=lambda x: x[0])
    return out


def parse_branch_anchor_sites(val: str):
    if not isinstance(val, str) or not val.strip():
        return []
    parts = [p.strip() for p in val.split(";") if p.strip()]
    out = []
    for part in parts:
        match = re.match(r"P(\d+)\s*[:=]?\s*([A-Za-z0-9@]+)", part)
        if match:
            out.append((int(match.group(1)), match.group(2)))
    return out


def parse_edges_core(val: str):
    if not isinstance(val, str) or not val.strip():
        return []
    parts = [p.strip() for p in val.split(";") if p.strip()]
    out = []
    for part in parts:
        match = re.match(r"([A-Za-z0-9]+)\(([^)]+)\)->([A-Za-z0-9]+)$", part)
        if match:
            out.append((match.group(1), match.group(2), match.group(3)))
    return out


def parse_bond_sequence(val: str):
    if not isinstance(val, str) or not val.strip():
        return []
    parts = re.split(r"[;,]", val)
    return [p.strip() for p in parts if p.strip()]


def _canon_branch_mono(mono: str) -> str:
    mono = (mono or "").strip()
    if not mono:
        return ""
    return mono if mono in ALLOWED_BRANCH_MONOS else "Other"

def _parse_mod_location(token: str) -> tuple[str, Optional[int]]:
    """Parse the location suffix of a modification token.

    The suffix indicates whether the target residue is in the core (``C``) or in a
    branch (``B``) and may include a positional index (e.g. ``B2``).  Bare numbers
    default to core positions, while an empty string returns an empty prefix.
    """

    if not isinstance(token, str):
        return "", None

    cleaned = token.strip()
    if not cleaned:
        return "", None

    match = re.match(r"([A-Za-z]*)(\d*)", cleaned)
    if not match:
        return "", None

    prefix = match.group(1).upper()
    number = match.group(2)

    if not prefix:
        prefix = "C" if number else ""

    idx = int(number) if number else None
    return prefix, idx


def _alpha_for_frequency(freq: Optional[float]) -> float:
    """Return an alpha value in [0, 1] based on the mean frequency."""

    if freq is None or (isinstance(freq, float) and math.isnan(freq)):
        return 0.6
    if freq < 30:
        return 0.4
    if freq < 60:
        return 0.65
    if freq < 90:
        return 0.85
    return 1.0


def _color_for_modification(name: str) -> str:
    """Return the colour for a modification name."""

    key = (name or "").strip().lower()
    if not key:
        return "#6c757d"
    if key.startswith("pyruv"):
        return "#1f77b4"
    if key.startswith("acetyl"):
        return "#d62728"
    return "#6c757d"


def _label_for_modification(name: str) -> str:
    mapping = {
        "pyruvylation": "Pyr",
        "pyruv": "Pyr",
        "acetylation": "Ac",
        "acetyl": "Ac",
        "formylation": "For",
        "formyl": "For",
        "lactylation": "Lac",
        "lactyl": "Lac",
        "glutamate": "Glu",
        "glutemate": "Glu",
        "glutamylation": "Glu",
        "glutam": "Glu",
    }

    key = (name or "").strip().lower()
    if not key:
        return "?"
    for prefix, label in mapping.items():
        if key.startswith(prefix):
            return label
    if "pyruv" in key:
        return "Pyr"
    if "acetyl" in key:
        return "Ac"
    if "formyl" in key:
        return "For"
    if "lactyl" in key:
        return "Lac"
    if "glut" in key:
        return "Glu"
    return "?"


def _normalise_k_type(name: Optional[str]) -> str:
    text = str(name or "")
    text = re.sub(r"\s+", "", text)
    return text.strip()


@dataclass
class ModificationTarget:
    kind: str
    key: str
    coords: tuple[float, float]
    parent_coords: Optional[tuple[float, float]] = None
    mono: Optional[str] = None
    order: Optional[int] = None


def _parse_mod_location(token: str) -> tuple[str, Optional[int]]:
    """Parse the location suffix of a modification token.

    The suffix indicates whether the target residue is in the core (``C``) or in a
    branch (``B``) and may include a positional index (e.g. ``B2``).  Bare numbers
    default to core positions, while an empty string returns an empty prefix.
    """

    if not isinstance(token, str):
        return "", None

    cleaned = token.strip()
    if not cleaned:
        return "", None

    match = re.match(r"([A-Za-z]*)(\d*)", cleaned)
    if not match:
        return "", None

    prefix = match.group(1).upper()
    number = match.group(2)

    if not prefix:
        prefix = "C" if number else ""

    idx = int(number) if number else None
    return prefix, idx


def _alpha_for_frequency(freq: Optional[float]) -> float:
    """Return an alpha value in [0, 1] based on the mean frequency."""

    if freq is None or (isinstance(freq, float) and math.isnan(freq)):
        return 0.6
    if freq < 30:
        return 0.4
    if freq < 60:
        return 0.65
    if freq < 90:
        return 0.85
    return 1.0


def _color_for_modification(name: str) -> str:
    """Return the colour for a modification name."""

    key = (name or "").strip().lower()
    if not key:
        return "#6c757d"
    if key.startswith("pyruv"):
        return "#1f77b4"
    if key.startswith("acetyl"):
        return "#d62728"
    return "#6c757d"


def _label_for_modification(name: str) -> str:
    mapping = {
        "pyruvylation": "Pyr",
        "pyruv": "Pyr",
        "acetylation": "Ac",
        "acetyl": "Ac",
        "formylation": "For",
        "formyl": "For",
        "lactylation": "Lac",
        "lactyl": "Lac",
        "glutamate": "Glu",
        "glutemate": "Glu",
        "glutamylation": "Glu",
        "glutam": "Glu",
    }

    key = (name or "").strip().lower()
    if not key:
        return "?"
    for prefix, label in mapping.items():
        if key.startswith(prefix):
            return label
    if "pyruv" in key:
        return "Pyr"
    if "acetyl" in key:
        return "Ac"
    if "formyl" in key:
        return "For"
    if "lactyl" in key:
        return "Lac"
    if "glut" in key:
        return "Glu"
    return "?"


def _normalise_k_type(name: Optional[str]) -> str:
    text = str(name or "")
    text = re.sub(r"\s+", "", text)
    return text.strip()


@dataclass
class ModificationTarget:
    kind: str
    key: str
    coords: tuple[float, float]
    parent_coords: Optional[tuple[float, float]] = None
    mono: Optional[str] = None
    order: Optional[int] = None
def parse_branch_paths_text(branch_text: str):
    if not isinstance(branch_text, str):
        return []

    paths = []
    seen = set()
    for raw_path in (p.strip() for p in branch_text.split(";") if p.strip()):
        s = raw_path
        s = s.replace("α", "a").replace("β", "b").replace("ALPHA", "a").replace("BETA", "b")
        s = s.replace("→", "-").replace("–", "-").replace("—", "-").replace("->", "-")
        s = re.sub(r"\s+", "", s)

        segments = []
        for part in (p for p in s.split(")") if p):
            mono, sep, right = part.partition("(")
            mono, right = mono.strip(), right.strip()
            if not sep or not mono or not right:
                continue
            segments.append((_canon_branch_mono(mono), right))

        if len(segments) < 2:
            continue

        anchor_mono, anchor_tok = segments[-1]
        match = _ANCHOR_TOKEN_RE.match(anchor_tok)
        if not match:
            continue
        anchor_pos = int(match.group(1))

        branch_segments = []
        for mono, bond in segments[:-1]:
            branch_segments.append({"mono": mono, "bond": bond})

        # The parsed segments are collected in the order they appear in the
        # textual representation (from the outermost monosaccharide towards the
        # anchor).  Downstream rendering expects the sequence to start at the
        # anchor and move outwards so that shared prefixes between different
        # branch strings collapse into a single path.  Reverse the collected
        # segments so that the first element describes the residue that attaches
        # directly to the anchor.
        branch_segments.reverse()

        if not branch_segments:
            continue

        key = (
            anchor_pos,
            anchor_mono,
            tuple((seg["mono"], seg["bond"]) for seg in branch_segments),
        )
        if key in seen:
            continue
        seen.add(key)

        paths.append(
            {
                "anchor_pos": anchor_pos,
                "anchor_mono": anchor_mono,
                "segments": branch_segments,
            }
        )

    return paths


def format_bond_label(label: Optional[str]):
    if label is None:
        return None
    text = str(label).strip()
    if not text:
        return None
    greek = {"a": "α", "b": "β"}
    first = text[0]
    repl = greek.get(first.lower())
    if repl:
        text = repl + text[1:]
    return text


def layout_core_positions(core_list, anchor_positions, core_edges, fallback_bonds=None):
    if not core_list:
        return [], {}, [], None

    ordered = sorted(core_list, key=lambda x: x[0])
    anchors = set(anchor_positions or [])
    seq_labels = list(fallback_bonds or [])
    seq_idx = 0

    def first_anchor_index(seq):
        for idx, (pos, _) in enumerate(seq):
            if pos in anchors:
                return idx
        return None

    forward = ordered
    backward = list(reversed(ordered))
    f_idx = first_anchor_index(forward)
    b_idx = first_anchor_index(backward)

    if f_idx is None and b_idx is None:
        display = forward
        rotate_idx = 0
    else:
        if f_idx is None:
            display, rotate_idx = backward, b_idx
        elif b_idx is None:
            display, rotate_idx = forward, f_idx
        else:
            if b_idx < f_idx:
                display, rotate_idx = backward, b_idx
            else:
                display, rotate_idx = forward, f_idx

    if rotate_idx and rotate_idx > 0:
        display = display[rotate_idx:] + display[:rotate_idx]

    core_nodes = []
    pos_to_x = {}
    for idx, (pos, mono) in enumerate(display):
        x = float(idx)
        entry = {"pos": pos, "mono": mono, "x": x}
        core_nodes.append(entry)
        pos_to_x[pos] = x

    core_bonds = []
    edge_lookup = {}
    for u, bond, v in core_edges:
        key = (u.strip(), v.strip())
        cleaned = bond.strip() if isinstance(bond, str) else bond
        edge_lookup.setdefault(key, []).append(cleaned)

    for i in range(len(display) - 1):
        from_pos, from_mono = display[i]
        to_pos, to_mono = display[i + 1]
        pair_key = (from_mono, to_mono)
        reverse_key = (to_mono, from_mono)

        labels = edge_lookup.get(pair_key)
        if labels:
            bond = labels.pop(0)
        else:
            rev_labels = edge_lookup.get(reverse_key)
            if rev_labels:
                bond = rev_labels.pop(0)
            else:
                bond = None

        seq_label = seq_labels[seq_idx] if seq_idx < len(seq_labels) else None
        if seq_idx < len(seq_labels):
            seq_idx += 1

        if bond is None:
            bond = seq_label

        core_bonds.append(
            {
                "from_pos": from_pos,
                "to_pos": to_pos,
                "bond": format_bond_label(bond),
            }
        )

    terminal_bond = format_bond_label(seq_labels[seq_idx]) if seq_idx < len(seq_labels) else None
    return core_nodes, pos_to_x, core_bonds, terminal_bond


class LabelPlacer:
    def __init__(self, pad=0.05):
        self.boxes = []
        self.pad = pad

    def _overlaps(self, box):
        x0, y0, x1, y1 = box
        for a0, b0, a1, b1 in self.boxes:
            if not (x1 < a0 or a1 < x0 or y1 < b0 or b1 < y0):
                return True
        return False

    def _measure(self, ax, x, y, text, **kwargs):
        fig = ax.figure
        label = ax.text(x, y, text, **kwargs)
        fig.canvas.draw()
        renderer = getattr(fig.canvas, "get_renderer", None)
        if renderer is not None:
            renderer = renderer()
        bb_disp = label.get_window_extent(renderer=renderer)
        label.remove()
        inv = ax.transData.inverted()
        bb_data = inv.transform(bb_disp.get_points())
        (x0, y0), (x1, y1) = bb_data
        pad = self.pad
        return (x0 - pad, y0 - pad, x1 + pad, y1 + pad)

    def place(self, ax, x, y, text, **kwargs):
        box = self._measure(ax, x, y, text, **kwargs)
        tries = 0
        max_tries = 100
        dx = 0.0
        dy = 0.0
        step = self.pad
        while self._overlaps(box) and tries < max_tries:
            dx += step
            dy += step
            x_try = x + dx
            y_try = y + dy
            box = self._measure(ax, x_try, y_try, text, **kwargs)
            tries += 1
        self.boxes.append(box)
        ax.text((box[0] + box[2]) / 2, (box[1] + box[3]) / 2, text, **kwargs)



# symbols formatting
R_OUT = 0.20
R_IN = 0.05
R_TRI = 0.25
LW_OUT = 3
LW_IN = 2
BRANCH_PERP_OFFSETS = [0.0, 0.14, -0.14, 0.26, -0.26, 0.38, -0.38]
CORE_ANGLES_DEG = [90, 135, 45, 180, 0, 225, -45, 270]
MOD_EDGE_LENGTH = 0.55
MOD_NODE_RADIUS = R_OUT / 3.0
MOD_EDGE_LINEWIDTH = LW_OUT / 2.0
MOD_NODE_EDGEWIDTH = LW_OUT / 2.0
UNASSIGNED_OFFSET_X = 2
UNASSIGNED_OFFSET_Y = -1

def _modification_angle_deg(index: int, *, anchored: bool) -> float:
    """Return the angular offset (in degrees) for a modification marker."""

    if not anchored:
        return 0.0
    return 60.0 if index % 2 == 0 else -60.0


def _angle_offset_deg(index: int) -> float:
    if index <= 0:
        return 0.0
    step = 20
    group = (index + 1) // 2
    sign = 1 if index % 2 else -1
    return float(step * group * sign)
def _draw_ring(ax, x, y, outer_color, inner_color="white"):
    circ = Circle((x, y), R_OUT, facecolor=outer_color, edgecolor="black", linewidth=LW_OUT)
    ax.add_patch(circ)
    inner = Circle((x, y), R_IN, facecolor=inner_color, edgecolor=inner_color, linewidth=LW_IN)
    ax.add_patch(inner)


def draw_node(ax, x, y, mono, label_position="above"):
    if mono == "Glc":
        circ = Circle((x, y), R_OUT, facecolor="#4ea3ff", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(circ)
    elif mono == "GlcA":
        _draw_ring(ax, x, y, outer_color="#4ea3ff")
    elif mono == "Man":
        circ = Circle((x, y), R_OUT, facecolor="#3cb371", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(circ)
    elif mono == "Rha":
        tri = RegularPolygon((x, y), numVertices=3, radius=R_TRI, facecolor="#3cb371", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(tri)
    elif mono == "Fru":
        tri = RegularPolygon((x, y), numVertices=3, radius=R_TRI, facecolor="#ff6b6b", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(tri)
    elif mono == "Gal":
        circ = Circle((x, y), R_OUT, facecolor="#ffbf00", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(circ)
    elif mono == "GalA":
        _draw_ring(ax, x, y, outer_color="#ffbf00")
    elif mono == "Galf":
        circ = Circle((x, y), R_OUT, facecolor="#ffbf00", edgecolor="black", linewidth=LW_OUT)
        ax.add_patch(circ)
        dot = Circle((x, y), R_IN, facecolor="white", edgecolor="white", linewidth=LW_IN)
        ax.add_patch(dot)
    else:
        circ = Circle((x, y), R_OUT, facecolor="lightgray", edgecolor="gray", linewidth=LW_OUT)
        ax.add_patch(circ)

    label_position = (label_position or "above").lower()
    if label_position == "below":
        ax.text(x, y - 0.25, mono, ha="center", va="top", fontsize=14, fontweight="bold")
    elif label_position == "left":
        ax.text(x - 0.25, y, mono, ha="right", va="center", fontsize=14, fontweight="bold")
    else:
        ax.text(x, y + 0.25, mono, ha="center", va="bottom", fontsize=14, fontweight="bold")


def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=3, mutation_scale=15))


def place_bond_label(ax, labeler, x1, y1, x2, y2, text, offset=0.14):
    if not text:
        return
    dx, dy = (x2 - x1), (y2 - y1)
    length = math.hypot(dx, dy) or 1.0
    ux, uy = dx / length, dy / length
    nx, ny = -uy, ux

    vertical = abs(dx) < 1e-6 and abs(dy) > 0
    if vertical:
        nx, ny = 1.0, 0.0
    else:
        if abs(ny) < 1e-9:
            nx, ny = 0.0, 1.0
        elif ny < 0:
            nx, ny = -nx, -ny
    norm = math.hypot(nx, ny) or 1.0
    nx, ny = nx / norm, ny / norm

    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    tx, ty = mx + nx * offset, my + ny * offset
    angle = math.degrees(math.atan2(dy, dx))
    if vertical:
        angle = 0.0
    else:
        if angle > 90:
            angle -= 180
        elif angle < -90:
            angle += 180

    text_kwargs = dict(
        fontsize=14,
        color="black",
        fontweight="bold",
        rotation=angle,
        rotation_mode="anchor",
        ha="center",
        va="center",
    )
    if vertical:
        text_kwargs["ha"] = "left"

    labeler.place(ax, tx, ty, text, **text_kwargs)


class BranchNode:
    __slots__ = ("name", "children", "depth", "layout_x")

    def __init__(self, name):
        self.name = name
        self.children = []
        self.depth = 0
        self.layout_x = 0.0

    def add_child(self, bond, child):
        for existing_bond, existing_child in self.children:
            if existing_bond == bond and existing_child.name == child.name:
                return existing_child
        self.children.append((bond, child))
        return child


def _sort_branch_children(node):
    node.children.sort(key=lambda item: ((item[0] or ""), item[1].name))
    for _, child in node.children:
        _sort_branch_children(child)


def compute_branch_depth(node):
    if not node.children:
        return 1
    return 1 + max(compute_branch_depth(child) for _, child in node.children)


def assign_branch_layout(node, depth=0, next_leaf=None):
    if next_leaf is None:
        next_leaf = [0]
    node.depth = depth
    if not node.children:
        node.layout_x = float(next_leaf[0])
        next_leaf[0] += 1
    else:
        child_positions = []
        for bond, child in node.children:
            assign_branch_layout(child, depth + 1, next_leaf)
            child_positions.append(child.layout_x)
        node.layout_x = sum(child_positions) / len(child_positions) if child_positions else float(next_leaf[0])
    return node


def _cluster_branch_entries(entries):
    if len(entries) <= 1:
        return [entries]

    parents = list(range(len(entries)))

    def find(idx):
        while parents[idx] != idx:
            parents[idx] = parents[parents[idx]]
            idx = parents[idx]
        return idx

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parents[rb] = ra

    def shared_prefix_len(first, second):
        length = 0
        for seg_a, seg_b in zip(first["segments"], second["segments"]):
            if seg_a["mono"] == seg_b["mono"] and seg_a["bond"] == seg_b["bond"]:
                length += 1
            else:
                break
        return length

    for i in range(len(entries)):
        for j in range(i + 1, len(entries)):
            if shared_prefix_len(entries[i], entries[j]) > 0:
                union(i, j)

    clusters: dict[int, list[dict[str, object]]] = {}
    for idx, entry in enumerate(entries):
        root_idx = find(idx)
        clusters.setdefault(root_idx, []).append(entry)
    return list(clusters.values())


def build_branch_trees(branch_paths):
    branches = []
    grouped: dict[tuple[int, str], list[dict[str, object]]] = {}
    for entry in branch_paths:
        key = (entry["anchor_pos"], entry["anchor_mono"])
        grouped.setdefault(key, []).append(entry)

    for (anchor_pos, anchor_mono), entries in grouped.items():
        for cluster in _cluster_branch_entries(entries):
            root = BranchNode(anchor_mono)
            for entry in cluster:
                node = root
                for segment in entry["segments"]:
                    child = BranchNode(segment["mono"])
                    node = node.add_child(segment["bond"], child)
            _sort_branch_children(root)
            assign_branch_layout(root)
            max_depth = compute_branch_depth(root)
            branches.append(
                {
                    "anchor_pos": anchor_pos,
                    "anchor_mono": anchor_mono,
                    "root": root,
                    "depth": max_depth,
                    "path_count": len(cluster),
                }
            )
    return branches


def render_branch(ax, branch, labeler, record: Optional[Callable] = None):
    root = branch["root"]
    anchor_x = branch.get("anchor_x", 0.0)
    fan_index = branch.get("fan_index", 0)
    side = branch.get("side", 1)
    if not side:
        side = 1

    spread = 0.8
    offset = fan_index * spread
    branch_x_offset = side * (SIBILINGS_SPREAD + offset)

    positions: dict[BranchNode, tuple[float, float]] = {}

    start_pos = (anchor_x + branch_x_offset, 0.0)
    positions[root] = start_pos

    def _angles_for_children(count: int):
        if count <= 0:
            return []
        if count == 1:
            return [math.pi / 2]
        min_angle = math.pi / 4
        max_angle = 3 * math.pi / 4
        if count == 2:
            return [min_angle, max_angle]
        step = (max_angle - min_angle) / (count - 1)
        return [min_angle + step * idx for idx in range(count)]

    def assign_positions(node):
        parent_pos = positions[node]
        children = node.children
        if not children:
            return
        angles = _angles_for_children(len(children))
        for idx, (bond, child) in enumerate(children):
            angle = angles[idx]
            dx = math.cos(angle) * BRANCH_EDGE_LENGTH * side
            dy = math.sin(angle) * BRANCH_EDGE_LENGTH
            positions[child] = (parent_pos[0] + dx, parent_pos[1] + dy)
            assign_positions(child)

    assign_positions(root)

    xs = [anchor_x, start_pos[0]]
    ys = [0.0, start_pos[1]]

    def draw_subtree(node, parent_coords, incoming_bond):
        x, y = positions[node]
        xs.append(x)
        ys.append(y)
        draw_arrow(ax, x, y, parent_coords[0], parent_coords[1])
        place_bond_label(
            ax,
            labeler,
            x,
            y,
            parent_coords[0],
            parent_coords[1],
            format_bond_label(incoming_bond),
            offset=0.14,
        )
        draw_node(ax, x, y, node.name, label_position="left")
        if record is not None:
            record(
                node,
                (x, y),
                parent_coords,
                branch,
                False,
            )
        for bond, child in node.children:
            draw_subtree(child, (x, y), bond)

    for bond, child in root.children:
        draw_subtree(child, (anchor_x, 0.0), bond)

    ax.plot([anchor_x], [0], marker="o", markersize=2.8, color="black")

    xs.append(anchor_x)
    ys.append(0.0)
    
    return (min(xs), max(xs), max(ys))


class KTypeStructureDrawer:
    """Generate structural diagrams for K-types."""

    def __init__(
        self,
        output_dir: Path,
        processed_csv: Path,
        modifications_csv: Optional[Path] = None,
        plots_subdir: str = "cps_drawn",
    ) -> None:
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.processed_csv = Path(processed_csv)
        self.modifications_csv = Path(modifications_csv) if modifications_csv else None
        self._modifications_df: Optional[pd.DataFrame] = None
        self.plots_subdir = plots_subdir
        self.png_dir = (self.output_dir / plots_subdir / "png").resolve()
        self.svg_dir = (self.output_dir / plots_subdir / "svg").resolve()
        self.png_dir.mkdir(parents=True, exist_ok=True)
        self.svg_dir.mkdir(parents=True, exist_ok=True)

    def _ensure_processed(self, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        if df is not None:
            return df
        return pd.read_csv(self.processed_csv)

    def _ensure_modifications_table(self) -> pd.DataFrame:
        if self._modifications_df is None:
            if self.modifications_csv is not None and self.modifications_csv.exists():
                df = pd.read_csv(self.modifications_csv)
            else:
                df = pd.DataFrame()
            df = df.copy()
            if not df.empty:
                df.columns = [str(c).strip() for c in df.columns]
                if "mean_frequency_percent" in df.columns:
                    df["mean_frequency_percent"] = pd.to_numeric(
                        df["mean_frequency_percent"], errors="coerce"
                    )
                if "K-type" in df.columns and "_K_norm" not in df.columns:
                    df["_K_norm"] = df["K-type"].apply(_normalise_k_type)
                elif "_K_norm" not in df.columns:
                    df["_K_norm"] = [_normalise_k_type(None)] * len(df)
                for col in ["monosaccharide", "modification", "bond_type"]:
                    if col in df.columns:
                        df[col] = df[col].fillna("").astype(str).str.strip()
                    else:
                        df[col] = [""] * len(df)
            else:
                df = pd.DataFrame(
                    columns=[
                        "K-type",
                        "monosaccharide",
                        "mean_frequency_percent",
                        "modification",
                        "bond_type",
                        "_K_norm",
                    ]
                )
            self._modifications_df = df
        return self._modifications_df

    def _get_modifications_for_k_type(self, k_type: Optional[str]) -> list[dict[str, object]]:
        df = self._ensure_modifications_table()
        if df.empty:
            return []
        key = _normalise_k_type(k_type)
        filtered = df[df["_K_norm"] == key]
        if filtered.empty:
            return []
        return filtered.to_dict("records")

    def _resolve_modification_target(
        self,
        token: str,
        core_node_info: dict[str, dict[str, object]],
        core_by_mono: dict[str, list[str]],
        branch_records_by_name: dict[str, list[dict[str, object]]],
    ) -> Optional[ModificationTarget]:
        if not isinstance(token, str) or not token.strip():
            return None

        mono_part, sep, location_part = token.partition("@")
        mono = mono_part.strip()
        location = location_part.strip() if sep else ""
        if not mono:
            return None

        if not location:
            candidates: list[ModificationTarget] = []
            for core_key in core_by_mono.get(mono, []):
                info = core_node_info.get(core_key)
                if info:
                    candidates.append(
                        ModificationTarget(
                            kind="core",
                            key=f"core:{core_key}",
                            coords=info["coords"],
                            mono=info.get("mono"),
                            order=info.get("index"),
                        )
                    )
            for entry in branch_records_by_name.get(mono, []):
                candidates.append(
                    ModificationTarget(
                        kind="branch",
                        key=f"branch:{mono}:{entry.get('order', 0)}",
                        coords=entry.get("coords", (0.0, 0.0)),
                        parent_coords=entry.get("parent"),
                        mono=mono,
                        order=entry.get("order"),
                    )
                )
            if len(candidates) == 1:
                return candidates[0]
            return None

        loc_type, loc_index = _parse_mod_location(location)
        if loc_type == "C" and loc_index is not None:
            core_key = f"C{loc_index}"
            info = core_node_info.get(core_key)
            if info:
                return ModificationTarget(
                    kind="core",
                    key=f"core:{core_key}",
                    coords=info["coords"],
                    mono=info.get("mono"),
                    order=info.get("index"),
                )
            # fallback to the first core entry of this monosaccharide if unique
            matches = [k for k in core_by_mono.get(mono, []) if k == core_key]
            if len(matches) == 1:
                info = core_node_info.get(matches[0])
                if info:
                    return ModificationTarget(
                        kind="core",
                        key=f"core:{matches[0]}",
                        coords=info["coords"],
                        mono=info.get("mono"),
                        order=info.get("index"),
                    )
            return None

        if loc_type == "B":
            entries = branch_records_by_name.get(mono, [])
            if not entries:
                return None
            chosen = None
            if loc_index is not None:
                # indices in the spreadsheet are 1-based
                idx = loc_index - 1
                for entry in entries:
                    order = entry.get("order")
                    if order is not None and order == idx:
                        chosen = entry
                        break
                if chosen is None and 0 <= idx < len(entries):
                    chosen = entries[idx]
            else:
                chosen = entries[0]
            if chosen:
                return ModificationTarget(
                    kind="branch",
                    key=f"branch:{mono}:{chosen.get('order', 0)}",
                    coords=chosen.get("coords", (0.0, 0.0)),
                    parent_coords=chosen.get("parent"),
                    mono=mono,
                    order=chosen.get("order"),
                )
            return None

        if loc_type == "C" and loc_index is None:
            entries = core_by_mono.get(mono, [])
            if len(entries) == 1:
                core_key = entries[0]
                info = core_node_info.get(core_key)
                if info:
                    return ModificationTarget(
                        kind="core",
                        key=f"core:{core_key}",
                        coords=info["coords"],
                        mono=info.get("mono"),
                        order=info.get("index"),
                    )
            return None

        return None

    def _draw_modifications(
        self,
        ax,
        labeler: LabelPlacer,
        core_node_info: dict[str, dict[str, object]],
        core_by_mono: dict[str, list[str]],
        branch_records_by_name: dict[str, list[dict[str, object]]],
        k_type: Optional[str],
        base_limits: tuple[float, float, float, float],
    ) -> Optional[tuple[float, float, float, float]]:
        mods = self._get_modifications_for_k_type(k_type)
        if not mods:
            return None

        extents: list[tuple[float, float, float, float]] = []
        counts: dict[str, int] = {}
        pending_unassigned: list[dict[str, object]] = []

        for record in mods:
            location = record.get("monosaccharide", "")
            if not location:
                pending_unassigned.append(record)
                continue
            target = self._resolve_modification_target(
                location, core_node_info, core_by_mono, branch_records_by_name
            )
            if not target:
                continue
            key = target.key
            idx = counts.get(key, 0)
            counts[key] = idx + 1
            extent = self._draw_modification_marker(ax, labeler, target, record, idx)
            if extent:
                extents.append(extent)

        if pending_unassigned:
            extent = self._draw_unassigned_modifications(
                ax, labeler, pending_unassigned, base_limits
            )
            if extent:
                extents.append(extent)

        if not extents:
            return None

        min_x = min(e[0] for e in extents)
        max_x = max(e[1] for e in extents)
        min_y = min(e[2] for e in extents)
        max_y = max(e[3] for e in extents)
        return (min_x, max_x, min_y, max_y)

    def _draw_modification_marker(
        self,
        ax,
        labeler: LabelPlacer,
        target: ModificationTarget,
        record: dict[str, object],
        index: int,
    ) -> Optional[tuple[float, float, float, float]]:
        freq_val = record.get("mean_frequency_percent")
        freq = None if pd.isna(freq_val) else float(freq_val)
        color = _color_for_modification(record.get("modification"))
        if freq is None:
            color = "#6c757d"
        alpha = _alpha_for_frequency(freq)
        mod_label = _label_for_modification(record.get("modification"))

        base = target.coords
        anchored = False
        base_direction = (1.0, 0.0)
        if target.kind == "core":
            anchored = True
        elif target.kind == "branch":
            parent = target.parent_coords
            if parent and (parent[0] != base[0] or parent[1] != base[1]):
                anchored = True
        else:
            return None

        angle_deg = _modification_angle_deg(index, anchored=anchored)
        if anchored:
            direction = _rotate_vector(base_direction, math.radians(angle_deg))
        else:
            direction = base_direction

        norm = math.hypot(*direction) or 1.0
        direction = (direction[0] / norm, direction[1] / norm)

        start_x = base[0] + direction[0] * R_OUT
        start_y = base[1] + direction[1] * R_OUT
        end_x = base[0] + direction[0] * (R_OUT + MOD_EDGE_LENGTH)
        end_y = base[1] + direction[1] * (R_OUT + MOD_EDGE_LENGTH)

        ax.plot([start_x, end_x], [start_y, end_y], color="black", linewidth=MOD_EDGE_LINEWIDTH)

        face_color = mcolors.to_rgba(color, alpha)
        edge_color = mcolors.to_rgba(color, min(alpha + 0.2, 1.0))
        marker = Circle(
            (end_x, end_y),
            MOD_NODE_RADIUS,
            facecolor=face_color,
            edgecolor='black',
            linewidth=MOD_NODE_EDGEWIDTH,
        )
        ax.add_patch(marker)

        label_x = end_x + direction[0] * 0.12
        label_y = end_y + direction[1] * 0.12
        ax.text(
            label_x - 0.02,
            label_y + 0.08,
            mod_label,
            fontsize=12,
            fontweight="bold",
            color="black",
            ha="center",
            va="center",
        )

        min_x = min(start_x, end_x) - MOD_NODE_RADIUS - 0.1
        max_x = max(start_x, end_x) + MOD_NODE_RADIUS + 0.1
        min_y = min(start_y, end_y) - MOD_NODE_RADIUS - 0.1
        max_y = max(start_y, end_y) + MOD_NODE_RADIUS + 0.2
        return (min_x, max_x, min_y, max_y)



    def _draw_unassigned_modifications(
        self,
        ax,
        labeler: LabelPlacer,
        records: list[dict[str, object]],
        base_limits: tuple[float, float, float, float],
    ) -> Optional[tuple[float, float, float, float]]:
        if not records:
            return None

        base_x = base_limits[1] - UNASSIGNED_OFFSET_X
        base_y = base_limits[3] + UNASSIGNED_OFFSET_Y

        extents: list[tuple[float, float, float, float]] = []
        vertical_spacing = MOD_EDGE_LENGTH + 0.25

        for idx, record in enumerate(records):
            y_center = base_y + idx * vertical_spacing
            end_x = base_x
            end_y = y_center
            start_y = end_y - MOD_EDGE_LENGTH
            start_x = end_x

            freq_val = record.get("mean_frequency_percent")
            freq = None if pd.isna(freq_val) else float(freq_val)
            color = _color_for_modification(record.get("modification"))
            if freq is None:
                color = "#6c757d"
            alpha = _alpha_for_frequency(freq)
            mod_label = _label_for_modification(record.get("modification"))

            ax.plot([start_x, end_x], [start_y, end_y], color="black", linewidth=MOD_EDGE_LINEWIDTH)
            face_color = mcolors.to_rgba(color, alpha)
            edge_color = mcolors.to_rgba(color, min(alpha + 0.2, 1.0))
            marker = Circle(
                (end_x, end_y),
                MOD_NODE_RADIUS,
                facecolor=face_color,
                edgecolor='black',
                linewidth=MOD_NODE_EDGEWIDTH,
            )
            ax.add_patch(marker)

            label_x = end_x
            label_y = end_y + 0.14
            ax.text(
                label_x,
                label_y,
                mod_label,
                fontsize=12,
                fontweight="bold",
                color="black",
                ha="center",
                va="center",
            )

            min_x = min(start_x, end_x) - MOD_NODE_RADIUS - 0.1
            max_x = max(start_x, end_x) + MOD_NODE_RADIUS + 0.3
            min_y = min(start_y, end_y) - MOD_NODE_RADIUS - 0.1
            max_y = max(start_y, end_y) + MOD_NODE_RADIUS + 0.1
            extents.append((min_x, max_x, min_y, max_y))

        if not extents:
            return None

        min_x = min(e[0] for e in extents)
        max_x = max(e[1] for e in extents)
        min_y = min(e[2] for e in extents)
        max_y = max(e[3] for e in extents)
        return (min_x, max_x, min_y, max_y)
    def draw_all(self, df: Optional[pd.DataFrame] = None, *, overwrite: bool = True, test: bool = False) -> list[Path]:
        processed = self._ensure_processed(df)
        written: list[Path] = []
        if test: print('WARNING! testing.')
        for _, row in processed.iterrows():
            kname = str(row.get("structure_id", "") or row.get("K_type", "")).strip()
            stem = (kname if kname else f"row_{_}").replace("/", "_").replace(" ", "_")

            if not test: pass
            elif test and kname in ['K3', 'K7', 'K24', 'K35', 'K43','K59', 'K71', 'K8_STRAIN_A4']: pass
            else: continue
            
            png_path = self.png_dir / f"{stem}.png"
            svg_path = self.svg_dir / f"{stem}.svg"
            if not overwrite and png_path.exists() and svg_path.exists():
                continue
            success, outputs = self.draw_single(row, stem_name=stem)
            if success:
                written.extend(outputs)
        return written

    def draw_single(self, row: pd.Series, stem_name: Optional[str] = None) -> tuple[bool, list[Path]]:
        stem = (stem_name or str(row.get("structure_id", row.get("K_type", "K?")))).replace("/", "_").replace(" ", "_")
        core_list = parse_core_positions(row.get("core_positions", ""))
        if not core_list:
            return False, []

        anchor_sites = parse_branch_anchor_sites(row.get("branch_anchor_sites", ""))
        branch_text = row.get("branch") or row.get("branches") or ""
        branch_paths = parse_branch_paths_text(branch_text)

        core_edges_list = parse_edges_core(row.get("core_edges", ""))
        core_bond_sequence = parse_bond_sequence(row.get("core_bonds_all", ""))

        anchor_positions = sorted({pos for pos, _ in anchor_sites})
        core_nodes, pos_map, core_bonds, core_terminal_bond = layout_core_positions(
            core_list,
            anchor_positions,
            core_edges_list,
            fallback_bonds=core_bond_sequence,
        )

        branches = build_branch_trees(branch_paths)
        core_node_info = {}
        core_by_mono: dict[str, list[str]] = {}
        for idx, node in enumerate(core_nodes):
            pos_key = f"C{node['pos']}"
            coords = (node["x"], 0.0)
            core_node_info[pos_key] = {"mono": node["mono"], "coords": coords, "index": idx}
            core_by_mono.setdefault(node["mono"], []).append(pos_key)

        branch_records_by_name: dict[str, list[dict[str, object]]] = {}

        def record_branch_node(node, coords, parent_coords, branch_data, is_anchor):
            name = getattr(node, "name", "")
            branch_records_by_name.setdefault(name, []).append(
                {
                    "coords": coords,
                    "parent": parent_coords,
                    "is_anchor": is_anchor,
                    "anchor_pos": branch_data.get("anchor_pos"),
                    "anchor_x": branch_data.get("anchor_x"),
                }
            )

        for branch in branches:
            branch["anchor_x"] = pos_map.get(branch["anchor_pos"], 0.0)
            if branch["root"].children:
                first_child = branch["root"].children[0]
                if first_child:
                    branch["anchor_bond"] = first_child[0]

        fig, ax = plt.subplots(figsize=FIG_SIZE)
        ax.set_title(f"{row.get('structure_id', row.get('K_type', 'K?'))}", fontsize=12, pad=12)
        labeler = LabelPlacer(pad=0.06)

        for node in core_nodes:
            draw_node(ax, node["x"], 0, node["mono"], label_position="below")
        terminal_arrow_end = None

        for bond_info in core_bonds:
            x1 = pos_map.get(bond_info["from_pos"], 0.0)
            x2 = pos_map.get(bond_info["to_pos"], 0.0)
            draw_arrow(ax, x1, 0, x2, 0)
            bond = bond_info.get("bond")
            if bond:
                place_bond_label(ax, labeler, x1, 0, x2, 0, bond, offset=0.20)

        if core_terminal_bond:
            last_x = core_nodes[-1]["x"] if core_nodes else 0.0
            terminal_arrow_end = last_x + 1.0
            draw_arrow(ax, last_x, 0, terminal_arrow_end, 0)
            place_bond_label(ax, labeler, last_x, 0, terminal_arrow_end, 0, core_terminal_bond, offset=0.20)

        by_anchor = {}
        for branch in branches:
            by_anchor.setdefault(branch["anchor_x"], []).append(branch)

        branch_extents = []
        for anchor_x, brs in by_anchor.items():
            brs = sorted(brs, key=lambda data: data.get("depth", 0), reverse=True)
            left, right = [], []
            for idx, branch in enumerate(brs):
                (left if idx % 2 == 0 else right).append(branch)
            for idx, branch in enumerate(left):
                branch["side"] = -1
                branch["fan_index"] = idx
                branch_extents.append(
                    render_branch(ax, branch, labeler, record=record_branch_node)
                )
            for idx, branch in enumerate(right):
                branch["side"] = +1
                branch["fan_index"] = idx
                branch_extents.append(
                    render_branch(ax, branch, labeler, record=record_branch_node)
                )

        for entries in branch_records_by_name.values():
            entries.sort(
                key=lambda item: (
                    item.get("anchor_pos", 0),
                    item.get("coords", (0.0, 0.0))[0],
                    item.get("coords", (0.0, 0.0))[1],
                )
            )
            for idx, entry in enumerate(entries):
                entry["order"] = idx

        core_xs = [node["x"] for node in core_nodes] or [0.0]
        if terminal_arrow_end is not None:
            core_xs.append(terminal_arrow_end)
        xmin = min(core_xs) - 1.5
        xmax = max(core_xs) + 1.5
        ymin = -1.5
        ymax = 2.5
        for mn, mx, my in branch_extents:
            xmin = min(xmin, mn - 0.8)
            xmax = max(xmax, mx + 0.8)
            ymax = max(ymax, my + 1.0)

        mod_extents = self._draw_modifications(
            ax,
            labeler,
            core_node_info,
            core_by_mono,
            branch_records_by_name,
            row.get("K_type"),
            (xmin, xmax, ymin, ymax),
        )
        if mod_extents:
            mod_xmin, mod_xmax, mod_ymin, mod_ymax = mod_extents
            xmin = min(xmin, mod_xmin)
            xmax = max(xmax, mod_xmax)
            ymin = min(ymin, mod_ymin)
            ymax = max(ymax, mod_ymax)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[:].set_visible(False)

        fig.tight_layout()
        png_path = self.png_dir / f"{stem}.png"
        svg_path = self.svg_dir / f"{stem}.svg"
        fig.savefig(png_path, dpi=300)
        fig.savefig(svg_path)
        plt.close(fig)
        return True, [png_path, svg_path]
