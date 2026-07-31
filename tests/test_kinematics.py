import dataclasses
import math

import pytest

from mochi.kinematics import (
    MM,
    RotaryGeometry,
    distance,
    port_position,
    prescribed_state,
    vane_fillet_geometry,
)


def test_supplied_geometry_uses_radial_difference_as_eccentricity() -> None:
    geometry = RotaryGeometry.default()

    assert geometry.radial_clearance_m == pytest.approx(4.5 * MM)
    assert geometry.eccentricity_m == pytest.approx(geometry.radial_clearance_m)
    assert geometry.speed_rpm == pytest.approx(1800.0)
    assert geometry.cutout_offset_m == pytest.approx(25.0 * MM)
    assert geometry.cutout_radius_m == pytest.approx(8.0 * MM)
    assert geometry.vane_width_m == pytest.approx(8.0 * MM)
    assert geometry.vane_tip_distance_at_top_m == pytest.approx(9.0 * MM)
    assert geometry.cylinder_height_m == pytest.approx(21.0 * MM)


@pytest.mark.parametrize("angle", [0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi])
def test_rotor_remains_tangent_and_cutout_stays_on_vane_axis(angle: float) -> None:
    geometry = RotaryGeometry.default()
    state = prescribed_state(geometry, angle)

    assert distance((0.0, 0.0), state.rotor_center_m) == pytest.approx(geometry.eccentricity_m)
    assert distance((0.0, 0.0), state.rotor_center_m) + geometry.rotor_radius_m == pytest.approx(
        geometry.cylinder_radius_m
    )
    assert state.cutout_center_m[0] == pytest.approx(0.0)
    assert distance(state.rotor_center_m, state.cutout_center_m) == pytest.approx(
        geometry.cutout_offset_m
    )


def test_cutout_center_moves_up_and_down_by_twice_the_eccentricity() -> None:
    geometry = RotaryGeometry.default()
    top = prescribed_state(geometry, 0.0)
    bottom = prescribed_state(geometry, math.pi)

    travel = top.cutout_center_m[1] - bottom.cutout_center_m[1]
    assert travel == pytest.approx(2.0 * geometry.eccentricity_m)


def test_positive_angle_rotates_clockwise_from_top() -> None:
    geometry = RotaryGeometry.default()
    top = prescribed_state(geometry, 0.0)
    after_top = prescribed_state(geometry, 0.25)

    assert top.rotor_center_m == pytest.approx((0.0, geometry.eccentricity_m))
    assert after_top.rotor_center_m[0] > 0.0
    assert after_top.rotor_center_m[1] < top.rotor_center_m[1]


def test_vane_tip_is_fixed_and_calibrated_at_top_dead_center() -> None:
    geometry = RotaryGeometry.default()
    fixed_tip_y = geometry.eccentricity_m + geometry.vane_tip_distance_at_top_m
    fixed_length = geometry.cylinder_radius_m - fixed_tip_y

    top = prescribed_state(geometry, 0.0)
    assert distance(top.rotor_center_m, top.vane_tip_m) == pytest.approx(
        geometry.vane_tip_distance_at_top_m
    )
    assert fixed_length == pytest.approx(25.0 * MM)

    for angle in (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi):
        state = prescribed_state(geometry, angle)
        assert state.vane_tip_m[0] == pytest.approx(0.0)
        assert state.vane_tip_m[1] == pytest.approx(fixed_tip_y)
        assert state.vane_length_m == pytest.approx(fixed_length)


def test_port_position_mirrors_about_the_positive_y_axis() -> None:
    """Opposite angles land on mirrored bore points.

    This pins the angle convention of the helper, not the port layout: the
    real ports are angular windows and are not symmetric (see
    ``tests/test_ports.py``).
    """

    radius = RotaryGeometry.default().cylinder_radius_m
    inlet = port_position(radius, 30.0)
    outlet = port_position(radius, -30.0)

    assert inlet[0] == pytest.approx(-outlet[0])
    assert inlet[1] == pytest.approx(outlet[1])
    assert inlet[0] > 0.0
    assert inlet[1] > 0.0


def test_vane_root_blend_is_tangent_to_the_flank_and_the_bore() -> None:
    """The R2.1 blend touches the flank at its edge and the bore internally."""

    geometry = RotaryGeometry.default()
    (centre_x, centre_y), flank_tangent, bore_tangent = vane_fillet_geometry(geometry)
    fillet = geometry.vane_cylinder_fillet_m
    half_width = 0.5 * geometry.vane_width_m

    # Flank tangent lies on the vane flank, one blend radius from the centre.
    assert flank_tangent[0] == pytest.approx(half_width)
    assert distance((centre_x, centre_y), flank_tangent) == pytest.approx(fillet)
    # Bore tangent lies on the bore, one blend radius inside from the centre.
    assert math.hypot(*bore_tangent) == pytest.approx(geometry.cylinder_radius_m)
    assert distance((centre_x, centre_y), bore_tangent) == pytest.approx(fillet)
    # Centre sits one blend radius outside the flank and inside the bore.
    assert centre_x == pytest.approx(half_width + fillet)
    assert math.hypot(centre_x, centre_y) == pytest.approx(geometry.cylinder_radius_m - fillet)


@pytest.mark.parametrize(
    "value, message",
    [
        (0.0, "positive"),
        (-1.0 * MM, "positive"),
        (math.inf, "finite"),
        (math.nan, "finite"),
        (20.0 * MM, "too large to meet the bore"),
    ],
)
def test_invalid_vane_root_blend_is_rejected(value: float, message: str) -> None:
    geometry = dataclasses.replace(RotaryGeometry.default(), vane_cylinder_fillet_m=value)

    with pytest.raises(ValueError, match=message):
        geometry.validate()


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (0.0, "positive"),
        (-1.0 * MM, "positive"),
        (math.inf, "finite"),
        (math.nan, "finite"),
        (5.0 * MM, "too large to leave a tip flat"),
    ],
)
def test_invalid_vane_tip_round_is_rejected(value: float, message: str) -> None:
    geometry = dataclasses.replace(RotaryGeometry.default(), vane_tip_fillet_m=value)

    with pytest.raises(ValueError, match=message):
        geometry.validate()


def test_invalid_intersecting_rotor_is_rejected() -> None:
    geometry = RotaryGeometry(
        cylinder_id_m=77.0 * MM,
        rotor_od_m=68.0 * MM,
        eccentricity_m=5.0 * MM,
        cutout_radius_m=8.0 * MM,
        cutout_offset_m=25.0 * MM,
        vane_width_m=8.0 * MM,
        vane_tip_distance_at_top_m=9.0 * MM,
        frequency_hz=30.0,
    )

    with pytest.raises(ValueError, match="intersect"):
        geometry.validate()


def test_nonpositive_cylinder_height_is_rejected() -> None:
    geometry = RotaryGeometry(
        cylinder_id_m=77.0 * MM,
        rotor_od_m=68.0 * MM,
        eccentricity_m=4.5 * MM,
        cutout_radius_m=8.0 * MM,
        cutout_offset_m=25.0 * MM,
        vane_width_m=8.0 * MM,
        vane_tip_distance_at_top_m=9.0 * MM,
        frequency_hz=30.0,
        cylinder_height_m=0.0,
    )

    with pytest.raises(ValueError, match="height"):
        geometry.validate()
