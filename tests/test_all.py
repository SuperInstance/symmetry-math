"""Comprehensive tests for symmetry-math — 40+ pytest cases."""

import math
import pytest

from symmetrymath.groups import (
    Group,
    SymmetricGroup,
    CyclicGroup,
    DihedralGroup,
    _SubgroupView,
)
from symmetrymath.transform import (
    rotation_matrix,
    reflection_matrix,
    glide_reflection_matrix,
    compose,
    apply,
    is_symmetry_of,
    _identity,
)
from symmetrymath.wallpaper import (
    Pattern,
    classify_pattern,
    detect_generators,
    verify_tiling,
    WALLPAPER_GROUPS,
    WallpaperType,
)
from symmetrymath.orbit import orbit, stabilizer, burnside_lemma
from symmetrymath.crystal import (
    CRYSTAL_POINT_GROUPS,
    laue_class,
    systematic_absences,
    point_group_from_hermann_mauguin,
    CrystalSystem,
)


# ===== groups.py ==========================================================

class TestCyclicGroup:
    def test_order(self):
        assert CyclicGroup(5).order() == 5

    def test_identity(self):
        g = CyclicGroup(7)
        assert g.identity == 0

    def test_operation(self):
        g = CyclicGroup(5)
        assert g.operation(2, 4) == 1

    def test_inverse(self):
        g = CyclicGroup(5)
        assert g.operation(3, g.inverse(3)) == 0

    def test_abelian(self):
        assert CyclicGroup(10).is_abelian()

    def test_subgroups(self):
        sg = CyclicGroup(6).subgroups()
        orders = sorted(s.order() for s in sg)
        assert orders == [1, 2, 3, 6]

    def test_powers(self):
        g = CyclicGroup(4)
        assert g.powers(1) == [0, 1, 2, 3, 0]

    def test_order_of(self):
        g = CyclicGroup(6)
        assert g.order_of(2) == 3
        assert g.order_of(1) == 6


class TestSymmetricGroup:
    def test_order(self):
        assert SymmetricGroup(3).order() == 6

    def test_identity(self):
        g = SymmetricGroup(3)
        assert g.identity == (0, 1, 2)

    def test_inverse(self):
        g = SymmetricGroup(3)
        p = (1, 2, 0)
        inv = g.inverse(p)
        assert g.operation(p, inv) == g.identity

    def test_not_abelian_s3(self):
        assert not SymmetricGroup(3).is_abelian()

    def test_s2_abelian(self):
        assert SymmetricGroup(2).is_abelian()

    def test_subgroups_s3(self):
        sg = SymmetricGroup(3).subgroups()
        orders = sorted(s.order() for s in sg)
        assert orders == [1, 2, 2, 2, 3, 6]


class TestDihedralGroup:
    def test_order(self):
        assert DihedralGroup(4).order() == 8

    def test_identity(self):
        g = DihedralGroup(5)
        assert g.identity == (0, 0)

    def test_inverse(self):
        g = DihedralGroup(5)
        e = g.identity
        assert g.operation((2, 0), g.inverse((2, 0))) == e

    def test_not_abelian(self):
        assert not DihedralGroup(4).is_abelian()

    def test_reflection_order_2(self):
        g = DihedralGroup(5)
        r = (0, 1)
        assert g.operation(r, r) == g.identity


# ===== transform.py ========================================================

class TestRotation:
    def test_90_degrees(self):
        R = rotation_matrix(math.pi / 2)
        x, y = apply(R, (1.0, 0.0))
        assert abs(x) < 1e-9
        assert abs(y - 1.0) < 1e-9

    def test_360_is_identity(self):
        R = rotation_matrix(2 * math.pi)
        I = _identity()
        for i in range(3):
            for j in range(3):
                assert abs(R[i][j] - I[i][j]) < 1e-9

    def test_about_centre(self):
        R = rotation_matrix(math.pi, 1.0, 0.0)
        x, y = apply(R, (2.0, 0.0))
        assert abs(x - 0.0) < 1e-9
        assert abs(y - 0.0) < 1e-9


class TestReflection:
    def test_x_axis(self):
        M = reflection_matrix(0.0)
        x, y = apply(M, (3.0, 5.0))
        assert abs(x - 3.0) < 1e-9
        assert abs(y + 5.0) < 1e-9

    def test_self_inverse(self):
        M = reflection_matrix(math.pi / 4)
        M2 = compose(M, M)
        I = _identity()
        for i in range(3):
            for j in range(3):
                assert abs(M2[i][j] - I[i][j]) < 1e-9


class TestGlideReflection:
    def test_compose_with_self_is_translation(self):
        G = glide_reflection_matrix(0.0, 2.0)
        G2 = compose(G, G)
        # should be a pure translation by (4, 0)
        assert abs(G2[0][2] - 4.0) < 1e-9
        assert abs(G2[1][2]) < 1e-9


class TestIsSymmetryOf:
    def test_square_rotation(self):
        sq = [(1, 1), (-1, 1), (-1, -1), (1, -1)]
        R = rotation_matrix(math.pi / 2)
        assert is_symmetry_of(R, sq)

    def test_non_symmetry(self):
        tri = [(0, 0), (1, 0), (0.5, math.sqrt(3) / 2)]
        R = rotation_matrix(math.pi / 2)
        assert not is_symmetry_of(R, tri)

    def test_empty_shape(self):
        assert is_symmetry_of(_identity(), [])


# ===== wallpaper.py ========================================================

class TestWallpaperGroups:
    def test_seventeen_types(self):
        assert len(WALLPAPER_GROUPS) == 17

    def test_classify_p1(self):
        # An asymmetric motif on a lattice
        motif = [(0.0, 0.0), (0.3, 0.1), (0.1, 0.4)]
        lv = ((1.0, 0.0), (0.0, 1.0))
        p = Pattern(motif, lv)
        assert classify_pattern(p) == "p1"

    def test_verify_tiling(self):
        motif = [(0.0, 0.0)]
        lv = ((1.0, 0.0), (0.0, 1.0))
        assert verify_tiling(motif, lv, radius=3.0, tol=0.1)

    def test_detect_generators(self):
        motif = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        lv = ((2.0, 0.0), (0.0, 2.0))
        p = Pattern(motif, lv)
        gens = detect_generators(p)
        assert isinstance(gens, list)


# ===== orbit.py ============================================================

class TestOrbit:
    def test_orbit_size(self):
        g = CyclicGroup(6)
        orb = orbit(g, 2)
        assert len(orb) == 6

    def test_orbit_cyclic_subgroup(self):
        g = CyclicGroup(6)
        # orbit under left-multiplication of element 3 is the full group
        orb = orbit(g, 3)
        assert orb == {0, 1, 2, 3, 4, 5}


class TestStabilizer:
    def test_identity_only(self):
        g = CyclicGroup(5)
        stab = stabilizer(g, 1)
        assert stab == {0}

    def test_full_group(self):
        g = CyclicGroup(5)
        # Left-multiplication on element 0: g·0 = g, only identity gives 0
        stab = stabilizer(g, 0)
        assert len(stab) == 1
        assert stab == {0}


class TestBurnside:
    def test_necklaces_2_colors_3_beads(self):
        """C₃ acting on 3 positions, each 2-colourable → 4 distinct necklaces."""
        g = CyclicGroup(3)
        def fixed(g_elem):
            # positions fixed by rotation by k: if k=0 all 3 fixed, else 1
            # (only the all-same-color colorings are fixed by nontrivial rotations)
            return 2 ** 3 if g_elem == 0 else 2 ** 1
        # (8 + 2 + 2) / 3 = 12 / 3 = 4
        assert burnside_lemma(g, fixed) == 4


# ===== crystal.py ==========================================================

class TestCrystalPointGroups:
    def test_thirty_two_groups(self):
        assert len(CRYSTAL_POINT_GROUPS) == 32

    def test_lookup(self):
        pg = point_group_from_hermann_mauguin("m-3m")
        assert pg.order == 48
        assert pg.system == CrystalSystem.CUBIC

    def test_unknown_raises(self):
        with pytest.raises(ValueError):
            point_group_from_hermann_mauguin("nope")

    def test_laue_class(self):
        pg = point_group_from_hermann_mauguin("4mm")
        assert laue_class(pg) == "4/mmm"

    def test_systematic_absences(self):
        pg = point_group_from_hermann_mauguin("m-3m")
        rules = systematic_absences(pg)
        assert isinstance(rules, dict)
        assert len(rules) > 0

    def test_triclinic(self):
        pg = point_group_from_hermann_mauguin("1")
        assert pg.system == CrystalSystem.TRICLINIC
        assert pg.order == 1

    def test_cubic_highest_order(self):
        pg = point_group_from_hermann_mauguin("m-3m")
        assert pg.order == 48

    def test_repr(self):
        pg = point_group_from_hermann_mauguin("6/mmm")
        assert "6/mmm" in repr(pg)
