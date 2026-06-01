"""Orbit-stabiliser theorem, Burnside's lemma, and related counting tools."""

from __future__ import annotations

from typing import Callable, Hashable, Set, TYPE_CHECKING

if TYPE_CHECKING:
    from symmetrymath.groups import Group


def orbit(group: "Group", element: Hashable) -> Set[Hashable]:
    """Return the orbit of *element* under the *group* action.

    The action is the group operation itself (left multiplication).
    """
    return {group.operation(g, element) for g in group.elements}


def stabilizer(group: "Group", element: Hashable) -> Set[Hashable]:
    """Return the stabiliser subgroup of *element*.

    Elements *g* ∈ *G* satisfying g·x = x.
    """
    return {g for g in group.elements if group.operation(g, element) == element}


def burnside_lemma(
    group: "Group",
    colorings: Callable[[Hashable], int],
) -> int:
    """Count distinct colorings using Burnside's lemma.

    Parameters
    ----------
    group : Group
        The symmetry group acting on the positions.
    colorings : callable
        ``colorings(g)`` must return the number of colorings fixed by
        group element *g*.  Equivalently, the number of colorings of the
        elements left unchanged by *g*.

    Returns
    -------
    int
        The number of orbits (distinct colourings) under the group
        action: ``|X/G| = (1/|G|) Σ_{g∈G} |X^g|``.
    """
    total = sum(colorings(g) for g in group.elements)
    return total // group.order()
