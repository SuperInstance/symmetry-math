# symmetry-math

Mathematics of symmetry groups — the bridge between art and physics.

Symmetry appears everywhere: Islamic tile patterns, molecular structure, particle physics.

## Installation

```bash
pip install symmetry-math
```

## Modules

- **groups** — `Group` base class, `SymmetricGroup(n)`, `CyclicGroup(n)`, `DihedralGroup(n)`, `is_abelian()`, `order()`, `subgroups()`
- **transform** — 2D rotation, reflection, glide reflection matrices, `compose()`, `is_symmetry_of()`
- **wallpaper** — 17 wallpaper groups, `classify_pattern()`, generator detection, tiling verification
- **orbit** — Orbit-stabilizer theorem, `orbit()`, `stabilizer()`, `burnside_lemma()`
- **crystal** — 32 crystallographic point groups, Laue classes, systematic absences

No external dependencies. Type hints, docstrings throughout.
