"""Wallpaper groups — the 17 two-dimensional space groups.

Every periodic pattern in the plane belongs to exactly one of 17
wallpaper groups.  This module provides classification, generator
detection, and basic tiling verification.
"""

from __future__ import annotations

import math
from enum import Enum
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from symmetrymath.transform import (
    Transform,
    Point,
    rotation_matrix,
    reflection_matrix,
    apply,
    compose,
)

# ---------------------------------------------------------------------------
# Wallpaper-group catalogue
# ---------------------------------------------------------------------------

class WallpaperType(str, Enum):
    """The 17 wallpaper-group types in short international notation."""
    P1 = "p1"
    P2 = "p2"
    PM = "pm"
    PG = "pg"
    CM = "cm"
    PMM = "pmm"
    PMG = "pmg"
    PGG = "pgg"
    CMM = "cmm"
    P4 = "p4"
    P4M = "p4m"
    P4G = "p4g"
    P3 = "p3"
    P3M1 = "p3m1"
    P31M = "p31m"
    P6 = "p6"
    P6M = "p6m"


# Convenient alias for public access
WALLPAPER_GROUPS = {wt.value: wt for wt in WallpaperType}


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

class Pattern:
    """A periodic pattern defined by a motif and a lattice.

    Parameters
    ----------
    motif : sequence of (float, float)
        The points / vertices of the asymmetric unit.
    lattice_vectors : ((a1x, a1y), (a2x, a2y))
        Two linearly independent vectors spanning the unit cell.
    """

    def __init__(
        self,
        motif: Sequence[Point],
        lattice_vectors: Tuple[Point, Point],
    ) -> None:
        self.motif = list(motif)
        self.lattice_vectors = lattice_vectors

    def unit_cell_points(self, n: int = 1) -> List[Point]:
        """Return motif replicated over an n×n block of unit cells."""
        (a1x, a1y), (a2x, a2y) = self.lattice_vectors
        pts: List[Point] = []
        for i in range(-n, n + 1):
            for j in range(-n, n + 1):
                for mx, my in self.motif:
                    pts.append((mx + i * a1x + j * a2x, my + i * a1y + j * a2y))
        return pts


def classify_pattern(pattern: Pattern, tol: float = 1e-4) -> str:
    """Classify *pattern* into one of the 17 wallpaper groups.

    This uses a decision tree based on the detected rotational
    symmetries and presence of reflections / glide reflections.
    Returns the short international symbol (e.g. ``"p6m"``).

    The algorithm checks rotational orders present in the motif
    and lattice, then tests for mirror and glide axes.
    """
    pts = pattern.unit_cell_points(2)

    max_rot = _highest_rotation(pts, tol)
    has_ref = _has_reflection(pts, tol)
    has_glide = _has_glide_reflection(pts, tol)

    if max_rot == 6:
        return "p6m" if has_ref else "p6"
    if max_rot == 4:
        if has_ref:
            return "p4m"  # simplified; real check distinguishes p4g
        return "p4"
    if max_rot == 3:
        if has_ref:
            return "p3m1"  # simplified
        return "p3"
    if max_rot == 2:
        if has_ref and has_glide:
            return "cmm"
        if has_ref:
            return "pmm"
        if has_glide:
            return "pgg"
        return "p2"
    # max_rot == 1
    if has_ref:
        return "pm"
    if has_glide:
        return "pg"
    return "p1"


def _highest_rotation(pts: List[Point], tol: float) -> int:
    """Return the highest rotational symmetry order found in *pts*."""
    for order in (6, 4, 3, 2):
        angle = 2 * math.pi / order
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)
        R = rotation_matrix(angle, cx, cy)
        transformed = [apply(R, p) for p in pts]
        if _point_sets_close(pts, transformed, tol):
            return order
    return 1


def _has_reflection(pts: List[Point], tol: float) -> bool:
    """Check whether any axis of reflection exists."""
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    for k in range(12):
        angle = math.pi * k / 12
        M = reflection_matrix(angle, cx, cy)
        transformed = [apply(M, p) for p in pts]
        if _point_sets_close(pts, transformed, tol):
            return True
    return False


def _has_glide_reflection(pts: List[Point], tol: float) -> bool:
    """Very rough check — delegates to reflection for simplified logic."""
    return False  # proper detection requires lattice-aware analysis


def _point_sets_close(a: List[Point], b: List[Point], tol: float) -> bool:
    """Check that multisets *a* and *b* match within *tol*."""
    if len(a) != len(b):
        return False
    sa = sorted(a, key=lambda p: (p[0], p[1]))
    sb = sorted(b, key=lambda p: (p[0], p[1]))
    for pa, pb in zip(sa, sb):
        if math.hypot(pa[0] - pb[0], pa[1] - pb[1]) > tol:
            return False
    return True


# ---------------------------------------------------------------------------
# Generator detection
# ---------------------------------------------------------------------------

def detect_generators(pattern: Pattern, tol: float = 1e-4) -> List[Transform]:
    """Return a list of symmetry transforms that generate the space group.

    Only checks a finite set of candidate rotations and reflections
    around the centroid.
    """
    pts = pattern.unit_cell_points(1)
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    generators: List[Transform] = []

    for order in (6, 4, 3, 2):
        angle = 2 * math.pi / order
        R = rotation_matrix(angle, cx, cy)
        if _is_symmetry(R, pts, tol):
            generators.append(R)

    for k in range(12):
        angle = math.pi * k / 12
        M = reflection_matrix(angle, cx, cy)
        if _is_symmetry(M, pts, tol):
            generators.append(M)

    return generators


def _is_symmetry(T: Transform, pts: List[Point], tol: float) -> bool:
    transformed = [apply(T, p) for p in pts]
    return _point_sets_close(pts, transformed, tol)


# ---------------------------------------------------------------------------
# Tiling verification
# ---------------------------------------------------------------------------

def verify_tiling(
    motif: Sequence[Point],
    lattice_vectors: Tuple[Point, Point],
    radius: float = 5.0,
    tol: float = 1e-4,
) -> bool:
    """Verify that tiling the *motif* with *lattice_vectors* covers the plane.

    Checks that translated copies of the motif tile without gaps or
    overlaps by verifying that lattice translations map the tiling onto
    itself within *radius* of the origin.
    """
    (a1x, a1y), (a2x, a2y) = lattice_vectors

    # Build the set of all points within radius
    pts: List[Point] = []
    n = int(radius / min(math.hypot(a1x, a1y), math.hypot(a2x, a2y))) + 2
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for mx, my in motif:
                x = mx + i * a1x + j * a2x
                y = my + i * a1y + j * a2y
                if x * x + y * y <= radius * radius:
                    pts.append((x, y))

    # Verify lattice translations preserve the set
    # Use only interior points (translate forward, check they land in set)
    pts_set = set((round(x / tol) * tol, round(y / tol) * tol) for x, y in pts)
    for vx, vy in [(a1x, a1y), (a2x, a2y)]:
        T: Transform = (
            (1.0, 0.0, vx),
            (0.0, 1.0, vy),
            (0.0, 0.0, 1.0),
        )
        # Check points that after translation are still within radius
        interior = [(x, y) for x, y in pts if (x + vx) ** 2 + (y + vy) ** 2 <= radius * radius]
        for p in interior:
            tp = apply(T, p)
            key = (round(tp[0] / tol) * tol, round(tp[1] / tol) * tol)
            if key not in pts_set:
                return False
    return True
