"""Crystallographic point groups — the 32 three-dimensional point groups.

Provides look-up tables, Laue-class mapping, and systematic-absence
rules for the 32 crystallographic point groups organised by crystal
system.
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, FrozenSet, List, Optional, Set, Tuple


class CrystalSystem(str, Enum):
    """The seven crystal systems."""
    TRICLINIC = "triclinic"
    MONOCLINIC = "monoclinic"
    ORTHORHOMBIC = "orthorhombic"
    TETRAGONAL = "tetragonal"
    TRIGONAL = "trigonal"
    HEXAGONAL = "hexagonal"
    CUBIC = "cubic"


class PointGroup:
    """A crystallographic point group.

    Attributes
    ----------
    symbol : str
        Short Hermann–Mauguin symbol (e.g. ``"4/mmm"``).
    system : CrystalSystem
        Parent crystal system.
    order : int
        Group order.
    laue : str
        Name of the corresponding Laue class.
    rotations : tuple of int
        Rotational symmetry orders present (besides 1).
    mirrors : bool
        Whether the group contains mirror planes.
    """

    __slots__ = ("symbol", "system", "order", "laue", "rotations", "mirrors")

    def __init__(
        self,
        symbol: str,
        system: CrystalSystem,
        order: int,
        laue: str,
        rotations: Tuple[int, ...] = (),
        mirrors: bool = False,
    ) -> None:
        self.symbol = symbol
        self.system = system
        self.order = order
        self.laue = laue
        self.rotations = rotations
        self.mirrors = mirrors

    def __repr__(self) -> str:
        return f"PointGroup({self.symbol!r}, system={self.system.value})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, PointGroup):
            return self.symbol == other.symbol
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.symbol)


# ---------------------------------------------------------------------------
# The 32 crystallographic point groups
# ---------------------------------------------------------------------------

_POINT_GROUPS: List[PointGroup] = [
    # Triclinic
    PointGroup("1", CrystalSystem.TRICLINIC, 1, "1"),
    PointGroup("-1", CrystalSystem.TRICLINIC, 2, "-1"),
    # Monoclinic
    PointGroup("2", CrystalSystem.MONOCLINIC, 2, "2/m", rotations=(2,)),
    PointGroup("m", CrystalSystem.MONOCLINIC, 2, "2/m", mirrors=True),
    PointGroup("2/m", CrystalSystem.MONOCLINIC, 4, "2/m", rotations=(2,), mirrors=True),
    # Orthorhombic
    PointGroup("222", CrystalSystem.ORTHORHOMBIC, 4, "mmm", rotations=(2,)),
    PointGroup("mm2", CrystalSystem.ORTHORHOMBIC, 4, "mmm", rotations=(2,), mirrors=True),
    PointGroup("mmm", CrystalSystem.ORTHORHOMBIC, 8, "mmm", rotations=(2,), mirrors=True),
    # Tetragonal
    PointGroup("4", CrystalSystem.TETRAGONAL, 4, "4/m", rotations=(4,)),
    PointGroup("-4", CrystalSystem.TETRAGONAL, 4, "4/m", rotations=()),
    PointGroup("4/m", CrystalSystem.TETRAGONAL, 8, "4/m", rotations=(4,), mirrors=True),
    PointGroup("422", CrystalSystem.TETRAGONAL, 8, "4/mmm", rotations=(4, 2)),
    PointGroup("4mm", CrystalSystem.TETRAGONAL, 8, "4/mmm", rotations=(4,), mirrors=True),
    PointGroup("-42m", CrystalSystem.TETRAGONAL, 8, "4/mmm", rotations=(2,), mirrors=True),
    PointGroup("4/mmm", CrystalSystem.TETRAGONAL, 16, "4/mmm", rotations=(4, 2), mirrors=True),
    # Trigonal
    PointGroup("3", CrystalSystem.TRIGONAL, 3, "3", rotations=(3,)),
    PointGroup("-3", CrystalSystem.TRIGONAL, 6, "-3", rotations=(3,), mirrors=True),
    PointGroup("32", CrystalSystem.TRIGONAL, 6, "-3m", rotations=(3, 2)),
    PointGroup("3m", CrystalSystem.TRIGONAL, 6, "-3m", rotations=(3,), mirrors=True),
    PointGroup("-3m", CrystalSystem.TRIGONAL, 12, "-3m", rotations=(3, 2), mirrors=True),
    # Hexagonal
    PointGroup("6", CrystalSystem.HEXAGONAL, 6, "6/m", rotations=(6,)),
    PointGroup("-6", CrystalSystem.HEXAGONAL, 6, "6/m", rotations=()),
    PointGroup("6/m", CrystalSystem.HEXAGONAL, 12, "6/m", rotations=(6,), mirrors=True),
    PointGroup("622", CrystalSystem.HEXAGONAL, 12, "6/mmm", rotations=(6, 2)),
    PointGroup("6mm", CrystalSystem.HEXAGONAL, 12, "6/mmm", rotations=(6,), mirrors=True),
    PointGroup("-62m", CrystalSystem.HEXAGONAL, 12, "6/mmm", rotations=(3, 2), mirrors=True),
    PointGroup("6/mmm", CrystalSystem.HEXAGONAL, 24, "6/mmm", rotations=(6, 2), mirrors=True),
    # Cubic
    PointGroup("23", CrystalSystem.CUBIC, 12, "m-3", rotations=(3, 2)),
    PointGroup("m-3", CrystalSystem.CUBIC, 24, "m-3", rotations=(3, 2), mirrors=True),
    PointGroup("432", CrystalSystem.CUBIC, 24, "m-3m", rotations=(4, 3, 2)),
    PointGroup("-43m", CrystalSystem.CUBIC, 24, "m-3m", rotations=(3, 2), mirrors=True),
    PointGroup("m-3m", CrystalSystem.CUBIC, 48, "m-3m", rotations=(4, 3, 2), mirrors=True),
]

# Index for fast look-up
_PG_BY_SYMBOL: Dict[str, PointGroup] = {pg.symbol: pg for pg in _POINT_GROUPS}

CRYSTAL_POINT_GROUPS: Dict[str, PointGroup] = dict(_PG_BY_SYMBOL)


def point_group_from_hermann_mauguin(symbol: str) -> PointGroup:
    """Look up a point group by its Hermann–Mauguin *symbol*."""
    pg = _PG_BY_SYMBOL.get(symbol)
    if pg is None:
        raise ValueError(f"Unknown point-group symbol: {symbol!r}")
    return pg


# ---------------------------------------------------------------------------
# Laue classes
# ---------------------------------------------------------------------------

_LAUE_CLASSES: Dict[str, List[str]] = {
    "-1": ["1", "-1"],
    "2/m": ["2", "m", "2/m"],
    "mmm": ["222", "mm2", "mmm"],
    "4/m": ["4", "-4", "4/m"],
    "4/mmm": ["422", "4mm", "-42m", "4/mmm"],
    "-3": ["3", "-3"],
    "-3m": ["32", "3m", "-3m"],
    "6/m": ["6", "-6", "6/m"],
    "6/mmm": ["622", "6mm", "-62m", "6/mmm"],
    "m-3": ["23", "m-3"],
    "m-3m": ["432", "-43m", "m-3m"],
}


def laue_class(pg: PointGroup) -> str:
    """Return the Laue class (Friedel class) of *pg*.

    The Laue class is the point group augmented by a centre of
    symmetry; it determines the diffraction symmetry.
    """
    return pg.laue


# ---------------------------------------------------------------------------
# Systematic absences
# ---------------------------------------------------------------------------

def systematic_absences(pg: PointGroup) -> Dict[str, str]:
    """Return a summary of systematic-absence rules for *pg*.

    Returns a dict mapping condition descriptions to brief explanations.
    This is a simplified pedagogical summary, not a full space-group
    treatment.
    """
    rules: Dict[str, str] = {}
    sym = pg.symbol

    # Centre-of-symmetry groups: Friedel's law
    if sym.startswith("-") or sym in ("2/m", "mmm", "4/m", "4/mmm",
                                       "6/m", "6/mmm", "m-3", "m-3m"):
        rules["hkl: h+k+l odd"] = "I-centring or body-diagonal glide"

    # Screw axes
    if "4" in pg.rotations:
        rules["00l: l ≠ 4n"] = "4₁ screw axis"
    if "6" in pg.rotations:
        rules["00l: l ≠ 6n"] = "6₁ screw axis"
    if "3" in pg.rotations:
        rules["00l: l ≠ 3n"] = "3₁ screw axis"
    if "2" in pg.rotations:
        rules["0k0: k ≠ 2n"] = "2₁ screw axis"

    # Glide planes
    if pg.mirrors:
        rules["h0l: h ≠ 2n"] = "a-glide perpendicular to b"
        rules["0kl: k ≠ 2n"] = "b-glide perpendicular to a"

    return rules
