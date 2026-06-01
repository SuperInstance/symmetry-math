# symmetry-math

> Mathematics of symmetry groups for Python — group theory, crystallography, wallpaper groups, and 2D transforms.

## What This Does

`symmetry-math` provides group theory and symmetry operations with practical applications in crystallography, pattern classification, and computational geometry. It builds cyclic, symmetric, and dihedral groups; computes orbits and stabilizers; applies 2D affine transforms (rotation, reflection, glide reflection); classifies periodic patterns into the 17 wallpaper groups; and provides crystallographic point group data. Use it for mathematical art generation, crystallography tools, or teaching abstract algebra.

## The Cultural Root

Symmetry is universal — from Islamic geometric patterns to the hexagonal symmetry of snowflakes to the bilateral symmetry of the human body. The mathematical formalization of symmetry as group theory (by Galois, then Klein's Erlangen program) revealed that **symmetry is not a property of objects but a structure: a group of transformations that leave the object invariant**. The 17 wallpaper groups were proven to be the only possible symmetries of 2D periodic patterns — a theorem with deep connections to crystallography and material science.

## Install

```bash
pip install symmetry-math
```

## Quick Start

```python
from symmetrymath import cyclic_group, symmetric_group, dihedral_group
from symmetrymath import orbit, stabilizer, burnside_lemma
from symmetrymath.transform import rotation_matrix, reflection_matrix, compose, apply
from symmetrymath.wallpaper import classify_pattern, Pattern
from symmetrymath.crystal import point_group_from_hermann_mauguin, laue_class

# Group theory
C4 = cyclic_group(4)
S3 = symmetric_group(3)
D6 = dihedral_group(6)

# Orbits and stabilizers
orb = orbit(element=1, group=D6)
stab = stabilizer(element=1, group=D6)

# Burnside's lemma — count distinct colorings
# For a necklace of 4 beads colored with 3 colors under C4:
distinct = burnside_lemma(group=C4, n_colors=3)

# 2D transforms
import math
R = rotation_matrix(math.pi / 4)        # 45° rotation
M = reflection_matrix(angle=math.pi/4)  # Reflect across diagonal
combined = compose(R, M)
point = apply(combined, (1.0, 0.0))

# Wallpaper group classification
motif = [(0, 0), (1, 0), (0.5, 0.5)]
pattern = Pattern(motif=motif, lattice_vectors=((2, 0), (0, 2)))
result = classify_pattern(pattern)
print(f"Wallpaper group: {result.type}")  # e.g., "p4m"

# Crystallography
pg = point_group_from_hermann_mauguin("4mm")
laue = laue_class(pg)
absences = pg.systematic_absences()
```

## API Reference

### Group Constructors

#### `cyclic_group(n) → Group`
Cₙ: rotations of an n-gon. Order n.

#### `symmetric_group(n) → Group`
Sₙ: all permutations of n elements. Order n!.

#### `dihedral_group(n) → Group`
Dₙ: symmetries of a regular n-gon (rotations + reflections). Order 2n.

### Group Operations

#### `orbit(element, group) → set`
All elements reachable from `element` under the group action.

#### `stabilizer(element, group) → Group`
Subgroup of elements that fix `element`.

#### `burnside_lemma(group, n_colors) → int`
Count distinct colorings: (1/|G|) Σ |Fix(g)|.

### `transform` module

#### `rotation_matrix(angle, cx=0, cy=0) → ndarray`
3×3 affine rotation matrix.

#### `reflection_matrix(angle, px=0, py=0) → ndarray`
Reflection across a line at `angle` through (px, py).

#### `glide_reflection_matrix(angle, glide, px=0, py=0) → ndarray`
Reflect then translate along the axis.

#### `compose(t1, t2) → ndarray`
Matrix product t1 · t2 (apply t2 first, then t1).

#### `apply(transform, point) → tuple`
Apply 3×3 affine transform to a 2D point.

#### `is_symmetry_of(transform, shape, tol=1e-6) → bool`
True if transform maps the shape onto itself.

### `wallpaper` module

#### `Pattern(motif, lattice_vectors)`
A periodic pattern defined by a motif and lattice basis vectors.

#### `classify_pattern(pattern) → WallpaperResult`
Classify into one of the 17 wallpaper groups using a decision tree on rotational symmetry, reflection, and glide reflection.

### `crystal` module

#### `CrystalSystem` (enum)
The seven crystal systems: `TRICLINIC`, `MONOCLINIC`, `ORTHORHOMBIC`, `TETRAGONAL`, `TRIGONAL`, `HEXAGONAL`, `CUBIC`.

#### `PointGroup`
Crystallographic point group with Hermann–Mauguin symbol, crystal system, and rotation orders.

#### `point_group_from_hermann_mauguin(symbol) → PointGroup`
Look up a point group by its HM symbol (e.g., "4mm", "m3m").

#### `laue_class(pg) → PointGroup`
Friedel class — point group augmented by inversion symmetry.

#### `systematic_absences(pg) → dict`
X-ray diffraction extinction rules for the point group.

## How It Works

**Groups** are represented as Cayley tables with identity, inverse, and compose operations. Cyclic groups use modular arithmetic; symmetric groups enumerate all permutations; dihedral groups compose rotations and reflections.

**Burnside's Lemma** counts distinct colorings: |X/G| = (1/|G|) Σ_{g∈G} |Fix(g)|, where Fix(g) is the set of colorings unchanged by g.

**Wallpaper Classification** uses a decision tree: detect highest rotational symmetry order (1, 2, 3, 4, or 6), then check for reflection axes and glide reflections. The combination uniquely determines one of the 17 groups.

**2D Transforms** use homogeneous coordinates (3×3 affine matrices) for unified rotation, reflection, translation, and composition.

## The Math

**Burnside's Lemma:** |Orbits| = (1/|G|) Σ_{g∈G} |Fix(g)|.

**Euler's Formula:** For periodic tilings, V − E + F = χ (Euler characteristic).

**Crystallographic Restriction Theorem:** Only rotations of order 1, 2, 3, 4, and 6 are compatible with 2D periodic lattices — this is why there are exactly 17 wallpaper groups and 230 space groups.

**References:**
- Armstrong, "Groups and Symmetry" (Springer)
- Conway et al., "The Symmetries of Things" (A.K. Peters)
- IUCr International Tables for Crystallography

## License

MIT
