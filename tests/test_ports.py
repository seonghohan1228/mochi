"""Port angular windows and the characteristic angles they define."""

import dataclasses
import math

import pytest

from mochi.kinematics import RotaryGeometry
from mochi.ports import (
    FULL_TURN_RAD,
    PortWindow,
    characteristic_angles,
    contains,
    discharge_port_radius_m,
    discharge_window,
    port_open_area_m2,
    suction_window,
)


def test_supplied_port_angles_give_the_expected_characteristic_angles() -> None:
    """The supplied phi, beta, gamma, delta reproduce the cycle events."""

    angles = characteristic_angles(RotaryGeometry.default())

    assert math.degrees(angles.suction_open_rad) == pytest.approx(10.4)
    assert math.degrees(angles.compression_start_rad) == pytest.approx(27.7)
    assert math.degrees(angles.discharge_open_rad) == pytest.approx(339.6)
    assert math.degrees(angles.discharge_close_rad) == pytest.approx(346.8)


def test_suction_window_spans_from_phi_to_beta() -> None:
    """The suction port opens at phi and its far edge is beta."""

    window = suction_window(RotaryGeometry.default())

    assert math.degrees(window.start_rad) == pytest.approx(10.4)
    assert math.degrees(window.end_rad) == pytest.approx(27.7)
    assert math.degrees(window.span_rad) == pytest.approx(17.3)


def test_the_two_ports_sit_on_opposite_sides_of_top_dead_center() -> None:
    """Suction opens just after top dead center; discharge just before it."""

    suction = suction_window(RotaryGeometry.default())
    discharge = discharge_window(RotaryGeometry.default())

    assert math.sin(suction.centre_rad) > 0.0
    assert math.sin(discharge.centre_rad) < 0.0


def test_discharge_window_span_matches_the_supplied_port_angle() -> None:
    window = discharge_window(RotaryGeometry.default())

    assert math.degrees(window.span_rad) == pytest.approx(7.2)
    assert math.degrees(window.centre_rad) == pytest.approx(343.2)


def test_discharge_port_radius_is_half_the_arc_it_subtends() -> None:
    """A 7.2 degree span on the 38.5 mm bore is a 4.84 mm wide port."""

    geometry = RotaryGeometry.default()

    radius_m = discharge_port_radius_m(geometry)

    assert 2.0 * radius_m == pytest.approx(4.838e-3, rel=1.0e-3)


@pytest.mark.parametrize(
    "angle_deg, expected",
    [
        (0.0, True),
        (5.0, True),
        (350.0, True),
        (355.0, True),
        (10.0, False),
        (180.0, False),
        (340.0, False),
    ],
)
def test_contains_handles_a_window_wrapping_through_zero(angle_deg: float, expected: bool) -> None:
    window = PortWindow(start_rad=math.radians(350.0), end_rad=math.radians(5.0))

    assert contains(window, math.radians(angle_deg)) is expected


def test_contains_accepts_both_window_bounds() -> None:
    window = PortWindow(start_rad=math.radians(339.6), end_rad=math.radians(346.8))

    assert contains(window, math.radians(339.6))
    assert contains(window, math.radians(346.8))
    assert not contains(window, math.radians(339.5))
    assert not contains(window, math.radians(346.9))


def test_contains_is_periodic() -> None:
    window = discharge_window(RotaryGeometry.default())
    angle = math.radians(343.2)

    assert contains(window, angle) is contains(window, angle + FULL_TURN_RAD)


def test_port_stays_fully_open_through_compression() -> None:
    """The contact only reaches the port near the end of the revolution.

    The compression chamber runs clockwise from the contact to the vane, so
    the port belongs to it for the whole compression stroke.
    """

    geometry = RotaryGeometry.default()
    full_area_m2 = math.pi * discharge_port_radius_m(geometry) ** 2

    for angle_deg in (0.0, 27.7, 90.0, 180.0, 300.0, 339.6):
        assert port_open_area_m2(geometry, math.radians(angle_deg)) == pytest.approx(full_area_m2)


def test_port_sweeps_shut_and_closes_exactly_at_the_recompression_angle() -> None:
    """The contact wipes the port closed from full area to zero.

    The port does not switch shut at one angle; it closes over its own
    width, which is what a later orifice flow model needs.
    """

    geometry = RotaryGeometry.default()
    angles = characteristic_angles(geometry)
    full_area_m2 = math.pi * discharge_port_radius_m(geometry) ** 2

    areas = [
        port_open_area_m2(geometry, math.radians(angle_deg))
        for angle_deg in (340.0, 341.5, 343.2, 345.0, 346.0)
    ]

    assert areas == sorted(areas, reverse=True)
    assert areas[2] == pytest.approx(0.5 * full_area_m2, rel=1.0e-6)
    assert port_open_area_m2(geometry, angles.discharge_close_rad) == pytest.approx(0.0)
    assert port_open_area_m2(geometry, math.radians(350.0)) == pytest.approx(0.0)


def test_the_supplied_ports_are_not_symmetric() -> None:
    """The old plus/minus 30 degree markers implied a symmetry that is gone."""

    angles = characteristic_angles(RotaryGeometry.default())
    suction_centre_deg = 0.5 * math.degrees(angles.suction_open_rad + angles.compression_start_rad)
    discharge_centre_deg = math.degrees(discharge_window(RotaryGeometry.default()).centre_rad)

    assert suction_centre_deg != pytest.approx(360.0 - discharge_centre_deg, abs=1.0)


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("suction_seal_angle_deg", 0.0, "positive"),
        ("discharge_port_span_deg", -1.0, "positive"),
        ("recompression_angle_deg", math.inf, "finite"),
        ("compression_start_angle_deg", math.nan, "finite"),
        ("compression_start_angle_deg", 5.0, "before compression"),
        ("recompression_angle_deg", 340.0, "after compression starts"),
    ],
)
def test_invalid_port_timing_is_rejected(field: str, value: float, message: str) -> None:
    geometry = dataclasses.replace(RotaryGeometry.default(), **{field: value})

    with pytest.raises(ValueError, match=message):
        characteristic_angles(geometry)


def test_port_open_area_rejects_a_non_finite_crank_angle() -> None:
    with pytest.raises(ValueError, match="finite"):
        port_open_area_m2(RotaryGeometry.default(), math.nan)
