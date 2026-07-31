import math
from dataclasses import replace

import pytest

from mochi.kinematics import MM, RotaryGeometry
from mochi.rotor_profile import MouthGeometry, rotor_contour

FULL_TURN = 2.0 * math.pi


def _radius(geometry: RotaryGeometry, crank_angle_rad: float, point: tuple[float, float]) -> float:
    eccentricity = geometry.eccentricity_m
    rotor = (eccentricity * math.sin(crank_angle_rad), eccentricity * math.cos(crank_angle_rad))
    return math.hypot(point[0] - rotor[0], point[1] - rotor[1])


def test_supplied_mouth_dimensions() -> None:
    mouth = MouthGeometry.default()

    assert math.degrees(mouth.flat_start_rad) == pytest.approx(34.8)
    assert math.degrees(mouth.flat_end_rad) == pytest.approx(13.4)
    assert mouth.flat_end_radius_m == pytest.approx(33.657 * MM)
    assert math.degrees(mouth.outlet_tangent_rad) == pytest.approx(13.0)
    assert math.degrees(mouth.inlet_sweep_rad) == pytest.approx(102.2)
    assert math.degrees(mouth.outlet_sweep_rad) == pytest.approx(92.6)
    assert mouth.lip_radius_m == pytest.approx(1.5 * MM)
    assert mouth.blend_radius_m == pytest.approx(1.0 * MM)


def test_flat_endpoints_sit_at_the_supplied_angles_and_radii() -> None:
    geometry = RotaryGeometry.default()
    mouth = MouthGeometry.default()
    contour = rotor_contour(geometry, 0.0)

    start, end = contour.inlet_flat
    assert _radius(geometry, 0.0, start) == pytest.approx(geometry.rotor_radius_m)
    assert _radius(geometry, 0.0, end) == pytest.approx(mouth.flat_end_radius_m)
    rotor_y = geometry.eccentricity_m
    assert math.atan2(start[0], start[1] - rotor_y) == pytest.approx(mouth.flat_start_rad)
    assert math.atan2(end[0], end[1] - rotor_y) == pytest.approx(mouth.flat_end_rad)
    assert start[0] > 0.0 and end[0] > 0.0  # the flat is on the inlet side only


def test_contour_edges_join_without_gaps() -> None:
    geometry = RotaryGeometry.default()

    for angle in (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi):
        contour = rotor_contour(geometry, angle)
        assert contour.inlet_flat[0] == pytest.approx(contour.od_arc[0])
        assert contour.inlet_flat[1] == pytest.approx(contour.mouth_path[0])
        assert contour.mouth_path[-1] == pytest.approx(contour.od_arc[-1])


def test_contour_never_exceeds_the_rotor_outside_diameter() -> None:
    geometry = RotaryGeometry.default()

    for index in range(72):
        angle = FULL_TURN * index / 72
        contour = rotor_contour(geometry, angle, od_points=121, groove_points=41)
        widest = max(_radius(geometry, angle, p) for p in contour.material)
        assert widest <= geometry.rotor_radius_m + 1.0e-12


def test_contour_stays_inside_the_bore_over_a_revolution() -> None:
    geometry = RotaryGeometry.default()

    for index in range(72):
        angle = FULL_TURN * index / 72
        contour = rotor_contour(geometry, angle, od_points=121, groove_points=41)
        reach = max(math.hypot(*p) for p in contour.material)
        assert reach <= geometry.cylinder_radius_m + 1.0e-12


def test_mouth_reaches_the_groove_on_both_sides() -> None:
    geometry = RotaryGeometry.default()
    contour = rotor_contour(geometry, 0.0)
    groove = (0.0, geometry.eccentricity_m + geometry.cutout_offset_m)

    on_groove = [
        p
        for p in contour.mouth_path
        if math.hypot(p[0] - groove[0], p[1] - groove[1])
        == pytest.approx(geometry.cutout_radius_m, abs=1.0e-9)
    ]
    assert any(p[0] > 0.0 for p in on_groove)
    assert any(p[0] < 0.0 for p in on_groove)


def test_the_two_lips_differ_because_the_mouth_is_asymmetric() -> None:
    geometry = RotaryGeometry.default()
    contour = rotor_contour(geometry, 0.0)
    rotor_y = geometry.eccentricity_m

    inlet_start = contour.mouth_path[0]
    outlet_start = contour.mouth_path[-1]
    inlet_angle = abs(math.atan2(inlet_start[0], inlet_start[1] - rotor_y))
    outlet_angle = abs(math.atan2(outlet_start[0], outlet_start[1] - rotor_y))
    assert math.degrees(inlet_angle) == pytest.approx(13.4)
    assert math.degrees(outlet_angle) == pytest.approx(13.0)
    assert _radius(geometry, 0.0, outlet_start) == pytest.approx(geometry.rotor_radius_m)


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("flat_end_rad", math.radians(40.0), "larger to a smaller"),
        ("flat_end_radius_m", 40.0 * MM, "inside the rotor"),
        ("lip_radius_m", 0.0, "radii must be positive"),
        ("inlet_sweep_rad", -1.0, "sweeps must be positive"),
        ("outlet_tangent_rad", 2.0, "between 0 and 90"),
        ("flat_start_rad", math.nan, "must be finite"),
    ],
)
def test_invalid_mouth_dimensions_are_rejected(field: str, value: float, message: str) -> None:
    geometry = RotaryGeometry.default()
    mouth = replace(MouthGeometry.default(), **{field: value})

    with pytest.raises(ValueError, match=message):
        rotor_contour(geometry, 0.0, mouth)
