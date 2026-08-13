"""Solid-body mass and inertia from geometry × an editable density (PHYSICS.md 4.8)."""

import pytest

from mochi.kinematics import RotaryGeometry
from mochi.mass_properties import (
    BUSH_DENSITY_KG_M3,
    JOURNAL_RADIUS_M,
    ROTOR_DENSITY_KG_M3,
    bush_mass_properties,
    bush_piece_mass_properties,
    rotor_mass_properties,
)

GRID_STEP_M = 8.0e-5  # coarser than the default for a fast test


@pytest.fixture(scope="module")
def geometry():
    return RotaryGeometry.default()


@pytest.fixture(scope="module")
def rotor(geometry):
    return rotor_mass_properties(geometry, grid_step_m=GRID_STEP_M)


@pytest.fixture(scope="module")
def bush(geometry):
    return bush_mass_properties(geometry, grid_step_m=GRID_STEP_M)


def test_rotor_mass_is_physical_cast_iron(rotor) -> None:
    """A 68 mm-OD × 21 mm cast-iron rotor (14.2 mm bore, mouth cut) is a few 100 g."""

    assert 0.2 < rotor.mass_kg < 0.8
    assert rotor.mass_kg == pytest.approx(rotor.area_m2 * rotor.height_m * rotor.density_kg_m3)
    assert rotor.density_kg_m3 == ROTOR_DENSITY_KG_M3


def test_rotor_inertia_matches_the_annulus_estimate(rotor, geometry) -> None:
    """Polar inertia about O_r is near the disc-minus-bore annulus (mouth is a small cut)."""

    annulus = 0.5 * rotor.mass_kg * (geometry.rotor_radius_m**2 + JOURNAL_RADIUS_M**2)
    assert rotor.inertia_kg_m2 == pytest.approx(annulus, rel=0.05)
    assert rotor.inertia_kg_m2 < annulus  # the mouth removes off-axis material


def test_bush_mass_is_small_and_positive(bush) -> None:
    """The two swing-bush pieces are light (tens of grams) and run axially short."""

    assert 0.0 < bush.mass_kg < 0.1
    assert bush.height_m < 0.021  # cylinder height minus 2× axial clearance
    assert bush.density_kg_m3 == BUSH_DENSITY_KG_M3


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0])
def test_density_scales_mass_and_inertia_linearly(geometry, factor) -> None:
    """Editing the density rescales mass and inertia linearly (the whole point)."""

    base = rotor_mass_properties(geometry, grid_step_m=GRID_STEP_M)
    scaled = rotor_mass_properties(
        geometry, density_kg_m3=factor * ROTOR_DENSITY_KG_M3, grid_step_m=GRID_STEP_M
    )
    assert scaled.mass_kg == pytest.approx(factor * base.mass_kg, rel=1.0e-12)
    assert scaled.inertia_kg_m2 == pytest.approx(factor * base.inertia_kg_m2, rel=1.0e-12)
    assert scaled.area_m2 == pytest.approx(base.area_m2, rel=1.0e-12)  # geometry unchanged


def test_rotor_rejects_bad_inputs(geometry) -> None:
    with pytest.raises(ValueError, match="positive"):
        rotor_mass_properties(geometry, density_kg_m3=-1.0, grid_step_m=GRID_STEP_M)
    with pytest.raises(ValueError, match="bore"):
        rotor_mass_properties(geometry, bore_radius_m=1.0, grid_step_m=GRID_STEP_M)


def test_bush_rejects_bad_density(geometry) -> None:
    with pytest.raises(ValueError, match="positive"):
        bush_mass_properties(geometry, density_kg_m3=0.0, grid_step_m=GRID_STEP_M)


def test_two_pieces_sum_to_the_whole_bush() -> None:
    """The IN and OUT piece masses add up to the two-piece bush mass (no gap/overlap)."""

    geometry = RotaryGeometry.default()
    both = bush_mass_properties(geometry, grid_step_m=GRID_STEP_M)
    in_piece = bush_piece_mass_properties(geometry, 1.0, grid_step_m=GRID_STEP_M)
    out_piece = bush_piece_mass_properties(geometry, -1.0, grid_step_m=GRID_STEP_M)
    assert in_piece.mass_kg + out_piece.mass_kg == pytest.approx(both.mass_kg, rel=1e-9)
    assert in_piece.mass_kg == pytest.approx(out_piece.mass_kg, rel=0.02)  # symmetric


def test_piece_inertia_and_cg_offset() -> None:
    geometry = RotaryGeometry.default()
    piece = bush_piece_mass_properties(geometry, 1.0, grid_step_m=GRID_STEP_M)
    assert piece.inertia_about_cg_kg_m2 > 0.0
    # the CM sits off the groove centre, toward the piece (away from the vane)
    assert piece.cg_offset_m[0] > 0.0


def test_piece_invalid_side_raises() -> None:
    with pytest.raises(ValueError):
        bush_piece_mass_properties(RotaryGeometry.default(), 0.0)
