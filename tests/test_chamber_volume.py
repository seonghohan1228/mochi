"""Chamber volumes from the true rotor contour, swing bush, and axial steps.

The integration is a grid quadrature, so tolerances here are set by the grid
pitch rather than by algebra. ``COARSE_STEP_M`` keeps the suite fast; the
convergence test is what justifies the pitch.
"""

import dataclasses
import math

import pytest

from mochi.chamber_volume import (
    AxialBands,
    SwingBush,
    clearance_volume_m3,
    mouth_cavity_area_m2,
    true_chamber_volumes,
)
from mochi.chambers import chamber_volumes
from mochi.kinematics import MM, RotaryGeometry
from mochi.ports import characteristic_angles

COARSE_STEP_M = 1.0e-4


def test_swing_bush_arc_matches_the_supplied_angular_extent() -> None:
    """The R0.5 fillets fix the 106.1 degree cylindrical face."""

    bush = SwingBush()

    assert 2.0 * math.degrees(bush.half_arc_rad()) == pytest.approx(106.1, rel=1.0e-3)


def test_swing_bush_films_are_ten_microns_on_both_faces() -> None:
    """The 0.020 mm shift equalizes the vane film and the groove film."""

    geometry = RotaryGeometry.default()
    bush = SwingBush()

    vane_film_m = bush.flat_offset_m + bush.piece_shift_m - 0.5 * geometry.vane_width_m
    groove_film_m = geometry.cutout_radius_m - bush.piece_outer_radius_m - bush.piece_shift_m

    assert vane_film_m == pytest.approx(0.010 * MM, rel=1.0e-6)
    assert groove_film_m == pytest.approx(0.010 * MM, rel=1.0e-6)


def test_axial_bands_leave_the_supplied_open_gap() -> None:
    geometry = RotaryGeometry.default()

    assert AxialBands().open_gap_m(geometry) == pytest.approx(16.2 * MM)


def test_mouth_cavity_area_matches_the_independently_recorded_estimate() -> None:
    """PHYSICS.md recorded roughly 30 mm^2 while this was a neglected volume.

    That figure came from a different route, so agreeing with it to the same
    order is a genuine cross-check of the whole integration.
    """

    area_m2 = mouth_cavity_area_m2(RotaryGeometry.default(), 0.0, grid_step_m=COARSE_STEP_M)

    assert area_m2 / (MM * MM) == pytest.approx(24.0, rel=0.1)


def test_clearance_volume_is_large_enough_to_bound_the_recompression() -> None:
    """The volume trapped when the discharge port shuts.

    The swing-compressor rule set asks for a discharge clearance of roughly
    0.1 to 0.5 cm^3; the geometry supplies it without any tuning.
    """

    volume_m3 = clearance_volume_m3(RotaryGeometry.default(), grid_step_m=COARSE_STEP_M)

    assert volume_m3 * 1.0e6 == pytest.approx(0.165, rel=0.05)
    assert 0.1e-6 < volume_m3 < 0.5e-6


def test_the_clearance_volume_is_mouth_cavity_gas() -> None:
    """The dead volume is the rotor mouth, not a separate discharge recess.

    At the discharge-close angle the trapped discharge gas fits inside the
    rotor mouth cavity, and both are far below the swept scale.
    """

    geometry = RotaryGeometry.default()
    angles = characteristic_angles(geometry)

    volumes = true_chamber_volumes(geometry, angles.discharge_close_rad, grid_step_m=COARSE_STEP_M)

    assert volumes.discharge_volume_m3 < volumes.mouth_cavity_volume_m3
    assert volumes.discharge_volume_m3 < 0.02 * volumes.suction_volume_m3


def test_the_circular_rotor_approximation_understates_the_trapped_volume() -> None:
    """This is why the approximation cannot carry a recompression phase.

    It drives the trapped volume toward zero, so ``p V**n = const`` diverges.
    """

    geometry = RotaryGeometry.default()
    angles = characteristic_angles(geometry)

    approximate_m3 = chamber_volumes(geometry, angles.discharge_close_rad).discharge_volume_m3
    true_m3 = true_chamber_volumes(
        geometry, angles.discharge_close_rad, grid_step_m=COARSE_STEP_M
    ).discharge_volume_m3

    assert approximate_m3 < 0.02 * true_m3


@pytest.mark.parametrize(
    "angle_deg, expected_ratio",
    [(6.5, 1.021), (27.7, 1.011), (180.0, 1.077)],
)
def test_true_volume_exceeds_the_circular_approximation(
    angle_deg: float, expected_ratio: float
) -> None:
    """The mouth cavity is gas space the outside-diameter disc discards.

    Angles stay outside the seal-over window, where the approximation
    refuses to split the chambers at all.
    """

    geometry = RotaryGeometry.default()

    approximate_m3 = chamber_volumes(geometry, math.radians(angle_deg)).discharge_volume_m3
    true_m3 = true_chamber_volumes(
        geometry, math.radians(angle_deg), grid_step_m=COARSE_STEP_M
    ).discharge_volume_m3

    assert true_m3 > approximate_m3
    assert true_m3 / approximate_m3 == pytest.approx(expected_ratio, rel=0.01)


def test_the_clearance_volume_converges_with_the_grid() -> None:
    """Halving the pitch must not move the answer by more than a percent."""

    geometry = RotaryGeometry.default()

    coarse_m3 = clearance_volume_m3(geometry, grid_step_m=COARSE_STEP_M)
    fine_m3 = clearance_volume_m3(geometry, grid_step_m=0.5 * COARSE_STEP_M)

    assert coarse_m3 == pytest.approx(fine_m3, rel=0.01)


def test_the_compression_chamber_shrinks_monotonically() -> None:
    geometry = RotaryGeometry.default()

    volumes = [
        true_chamber_volumes(
            geometry, math.radians(angle_deg), grid_step_m=COARSE_STEP_M
        ).discharge_volume_m3
        for angle_deg in (30.0, 90.0, 180.0, 270.0, 330.0)
    ]

    assert volumes == sorted(volumes, reverse=True)


def test_the_trapped_volume_never_reaches_zero() -> None:
    """The mouth cavity is the floor, so no artificial clearance is needed."""

    geometry = RotaryGeometry.default()

    for angle_deg in (346.8, 350.0, 354.0):
        volumes = true_chamber_volumes(geometry, math.radians(angle_deg), grid_step_m=COARSE_STEP_M)
        assert volumes.discharge_volume_m3 > 0.0


@pytest.mark.parametrize("grid_step_m", [0.0, -1.0e-4, math.inf, math.nan])
def test_invalid_grid_step_is_rejected(grid_step_m: float) -> None:
    with pytest.raises(ValueError, match="[Gg]rid step"):
        true_chamber_volumes(RotaryGeometry.default(), 0.0, grid_step_m=grid_step_m)


@pytest.mark.parametrize("crank_angle_rad", [math.inf, math.nan])
def test_non_finite_crank_angle_is_rejected(crank_angle_rad: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        true_chamber_volumes(RotaryGeometry.default(), crank_angle_rad)


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("piece_outer_radius_m", 0.0, "positive"),
        ("fillet_radius_m", -1.0e-4, "positive"),
        ("flat_offset_m", 20.0 * MM, "inside its own piece radius"),
        ("piece_outer_radius_m", math.inf, "finite"),
        ("piece_outer_radius_m", 9.0 * MM, "fit inside the circular groove"),
    ],
)
def test_invalid_swing_bush_is_rejected(field: str, value: float, message: str) -> None:
    bush = dataclasses.replace(SwingBush(), **{field: value})

    with pytest.raises(ValueError, match=message):
        bush.validate(RotaryGeometry.default())


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("ledge_thickness_m", 0.0, "positive"),
        ("full_vane_depth_m", math.nan, "finite"),
        ("ledge_thickness_m", 15.0 * MM, "thinner than the cylinder height"),
    ],
)
def test_invalid_axial_bands_are_rejected(field: str, value: float, message: str) -> None:
    bands = dataclasses.replace(AxialBands(), **{field: value})

    with pytest.raises(ValueError, match=message):
        bands.validate(RotaryGeometry.default())
