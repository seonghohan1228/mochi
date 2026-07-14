import math

import pytest

from mochi.kinematics import MM, RotaryGeometry, distance, port_position, prescribed_state


def test_supplied_geometry_uses_radial_difference_as_eccentricity() -> None:
    geometry = RotaryGeometry.default()

    assert geometry.radial_clearance_m == pytest.approx(4.5 * MM)
    assert geometry.eccentricity_m == pytest.approx(geometry.radial_clearance_m)
    assert geometry.speed_rpm == pytest.approx(1800.0)
    assert geometry.cutout_offset_m == pytest.approx(25.0 * MM)
    assert geometry.cutout_radius_m == pytest.approx(8.0 * MM)
    assert geometry.vane_width_m == pytest.approx(8.0 * MM)
    assert geometry.vane_tip_distance_at_top_m == pytest.approx(9.0 * MM)


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


def test_vane_tip_calibration_and_length_change() -> None:
    geometry = RotaryGeometry.default()
    top = prescribed_state(geometry, 0.0)
    side = prescribed_state(geometry, 0.5 * math.pi)
    bottom = prescribed_state(geometry, math.pi)

    assert distance(top.rotor_center_m, top.vane_tip_m) == pytest.approx(
        geometry.vane_tip_distance_at_top_m
    )
    assert top.vane_tip_m[0] == pytest.approx(0.0)
    assert side.vane_tip_m[0] == pytest.approx(0.0)
    assert side.vane_length_m != pytest.approx(top.vane_length_m)
    assert bottom.vane_length_m - top.vane_length_m == pytest.approx(2.0 * geometry.eccentricity_m)


def test_ports_are_symmetric_about_positive_y_axis() -> None:
    radius = RotaryGeometry.default().cylinder_radius_m
    inlet = port_position(radius, 30.0)
    outlet = port_position(radius, -30.0)

    assert inlet[0] == pytest.approx(-outlet[0])
    assert inlet[1] == pytest.approx(outlet[1])
    assert inlet[0] > 0.0
    assert inlet[1] > 0.0


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
