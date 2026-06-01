"""2-D geometric transformations and symmetry checks.

All transforms are represented as 3×3 affine matrices suitable for
homogeneous-coordinate multiplication.  Only rigid motions and glide
reflections are covered (no arbitrary scaling or shearing).
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence, Tuple

# Type alias — a 3×3 matrix stored as a tuple of 3 tuples of 3 floats.
Transform = Tuple[
    Tuple[float, float, float],
    Tuple[float, float, float],
    Tuple[float, float, float],
]

Point = Tuple[float, float]


def _mat_mul(A: Transform, B: Transform) -> Transform:
    """Multiply two 3×3 affine matrices."""
    result: list[list[float]] = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                result[i][j] += A[i][k] * B[k][j]
    return (
        (result[0][0], result[0][1], result[0][2]),
        (result[1][0], result[1][1], result[1][2]),
        (result[2][0], result[2][1], result[2][2]),
    )


def _identity() -> Transform:
    return (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


def rotation_matrix(angle: float, cx: float = 0.0, cy: float = 0.0) -> Transform:
    """Return a 3×3 matrix for rotation by *angle* radians about (*cx*, *cy*).

    >>> R = rotation_matrix(math.pi / 2)
    >>> round(R[0][0], 6)
    0.0
    """
    c, s = math.cos(angle), math.sin(angle)
    R: Transform = (
        (c, -s, cx * (1 - c) + cy * s),
        (s, c, cy * (1 - c) - cx * s),
        (0.0, 0.0, 1.0),
    )
    return R


def reflection_matrix(angle: float, px: float = 0.0, py: float = 0.0) -> Transform:
    """Return a reflection matrix across a line through (*px*, *py*) at *angle*.

    The line makes *angle* radians with the x-axis.
    """
    c2 = math.cos(2 * angle)
    s2 = math.sin(2 * angle)
    return (
        (c2, s2, px * (1 - c2) - py * s2),
        (s2, -c2, -px * s2 + py * (1 + c2)),
        (0.0, 0.0, 1.0),
    )


def glide_reflection_matrix(
    angle: float,
    glide: float,
    px: float = 0.0,
    py: float = 0.0,
) -> Transform:
    """Return a glide-reflection: reflect then translate *glide* along the axis.

    Parameters
    ----------
    angle : float
        Orientation of the mirror line (radians from x-axis).
    glide : float
        Distance to translate along the mirror line after reflecting.
    px, py : float
        A point the mirror line passes through.
    """
    refl = reflection_matrix(angle, px, py)
    tx = glide * math.cos(angle)
    ty = glide * math.sin(angle)
    T: Transform = (
        (1.0, 0.0, tx),
        (0.0, 1.0, ty),
        (0.0, 0.0, 1.0),
    )
    return _mat_mul(T, refl)


def compose(t1: Transform, t2: Transform) -> Transform:
    """Compose two transforms: first apply *t2*, then *t1*.

    Returns the matrix product *t1 · t2*.
    """
    return _mat_mul(t1, t2)


def apply(transform: Transform, point: Point) -> Point:
    """Apply a 3×3 affine *transform* to a 2-D *point*."""
    x, y = point
    nx = transform[0][0] * x + transform[0][1] * y + transform[0][2]
    ny = transform[1][0] * x + transform[1][1] * y + transform[1][2]
    return (nx, ny)


def is_symmetry_of(
    transform: Transform,
    shape: Iterable[Point],
    tol: float = 1e-6,
) -> bool:
    """Return ``True`` if *transform* maps *shape* onto itself.

    *shape* is an iterable of 2-D points.  Two shapes are considered
    equal when, after sorting, every corresponding point is within
    *tol* distance.
    """
    pts = list(shape)
    if not pts:
        return True
    transformed = [apply(transform, p) for p in pts]
    return _point_sets_equal(pts, transformed, tol)


def _point_sets_equal(
    a: list[Point], b: list[Point], tol: float
) -> bool:
    """Check that sorted point sets *a* and *b* coincide within *tol*."""
    if len(a) != len(b):
        return False
    sa = sorted(a, key=lambda p: (round(p[0] / tol) * tol, round(p[1] / tol) * tol))
    sb = sorted(b, key=lambda p: (round(p[0] / tol) * tol, round(p[1] / tol) * tol))
    for pa, pb in zip(sa, sb):
        if math.hypot(pa[0] - pb[0], pa[1] - pb[1]) > tol:
            return False
    return True
