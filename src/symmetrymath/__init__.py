"""symmetry-math — Mathematics of symmetry groups.

Provides group-theoretic primitives, 2D transforms, wallpaper-group
classification, orbit-stabiliser calculations, and crystallographic
point-group utilities.  Pure Python, no external dependencies.
"""

from symmetrymath.groups import (
    Group,
    SymmetricGroup,
    CyclicGroup,
    DihedralGroup,
)
from symmetrymath.transform import (
    rotation_matrix,
    reflection_matrix,
    glide_reflection_matrix,
    compose,
    is_symmetry_of,
    Transform,
)
from symmetrymath.wallpaper import (
    WALLPAPER_GROUPS,
    classify_pattern,
    detect_generators,
    verify_tiling,
)
from symmetrymath.orbit import orbit, stabilizer, burnside_lemma
from symmetrymath.crystal import (
    CRYSTAL_POINT_GROUPS,
    laue_class,
    systematic_absences,
    point_group_from_hermann_mauguin,
)

__all__ = [
    "Group",
    "SymmetricGroup",
    "CyclicGroup",
    "DihedralGroup",
    "rotation_matrix",
    "reflection_matrix",
    "glide_reflection_matrix",
    "compose",
    "is_symmetry_of",
    "Transform",
    "WALLPAPER_GROUPS",
    "classify_pattern",
    "detect_generators",
    "verify_tiling",
    "orbit",
    "stabilizer",
    "burnside_lemma",
    "CRYSTAL_POINT_GROUPS",
    "laue_class",
    "systematic_absences",
    "point_group_from_hermann_mauguin",
]
__version__ = "0.1.0"
