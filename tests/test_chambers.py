import math
from collections.abc import Callable

import pytest

from mochi.chambers import (
    DEFAULT_POLYTROPIC_EXPONENT,
    DISCHARGE_PORT_PRESSURE_PA,
    DISCHARGE_VALVE_RISE_FRACTION,
    FULL_TURN_RAD,
    RECESS_PRESSURE_PA,
    SUCTION_PORT_PRESSURE_PA,
    ChamberSeparationError,
    CycleTrace,
    _valve_opening_angle_rad,
    build_cycle_trace,
    chamber_areas,
    chamber_pressures,
    chamber_volumes,
    in_seal_over_window,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.kinematics import MM, RotaryGeometry, prescribed_state

FULL_TURN = 2.0 * math.pi
VALVE_OPENING_PA = DISCHARGE_PORT_PRESSURE_PA * (1.0 + DISCHARGE_VALVE_RISE_FRACTION)


def _separated_angles(geometry: RotaryGeometry, count: int = 24) -> list[float]:
    """Sample crank angles strictly inside the separated range."""

    margin = seal_over_half_angle_rad(geometry) + 1.0e-6
    return [margin + (FULL_TURN - 2.0 * margin) * index / (count - 1) for index in range(count)]


def _simpson(function: Callable[[float], float], lower: float, upper: float) -> float:
    samples = 8001
    step = (upper - lower) / (samples - 1)
    total = function(lower) + function(upper)
    for index in range(1, samples - 1):
        total += function(lower + index * step) * (4.0 if index % 2 else 2.0)
    return total * step / 3.0


def _numeric_suction_area(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Independent quadrature of the suction cross-section for validation."""

    eccentricity = geometry.eccentricity_m
    rotor_radius = geometry.rotor_radius_m
    cylinder_radius = geometry.cylinder_radius_m

    def rotor_boundary_distance(ray_angle_rad: float) -> float:
        relative = ray_angle_rad - crank_angle_rad
        chord = eccentricity * math.sin(relative)
        return eccentricity * math.cos(relative) + math.sqrt(rotor_radius**2 - chord**2)

    swept = _simpson(
        lambda ray: 0.5 * (cylinder_radius**2 - rotor_boundary_distance(ray) ** 2),
        0.0,
        crank_angle_rad,
    )
    state = prescribed_state(geometry, crank_angle_rad)
    rotor_x, rotor_y = state.rotor_center_m

    def strip_height(x_m: float) -> float:
        return (
            math.sqrt(cylinder_radius**2 - x_m**2)
            - rotor_y
            - math.sqrt(rotor_radius**2 - (x_m - rotor_x) ** 2)
        )

    strip = _simpson(strip_height, 0.0, 0.5 * geometry.vane_width_m)
    return swept - strip


def test_chamber_areas_conserve_the_crescent() -> None:
    geometry = RotaryGeometry.default()
    crescent = math.pi * (geometry.cylinder_radius_m**2 - geometry.rotor_radius_m**2)

    for angle in _separated_angles(geometry):
        areas = chamber_areas(geometry, angle)
        combined = areas.suction_area_m2 + areas.discharge_area_m2 + areas.vane_area_in_crescent_m2
        assert combined == pytest.approx(crescent, rel=1.0e-12)
        assert areas.crescent_area_m2 == pytest.approx(crescent, rel=1.0e-12)
        assert areas.suction_area_m2 >= 0.0
        assert areas.discharge_area_m2 >= 0.0


@pytest.mark.parametrize("angle_deg", [10.0, 30.0, 90.0, 180.0, 270.0, 350.0])
def test_suction_area_matches_independent_quadrature(angle_deg: float) -> None:
    geometry = RotaryGeometry.default()
    angle = math.radians(angle_deg)

    areas = chamber_areas(geometry, angle)
    numeric = _numeric_suction_area(geometry, angle)
    assert areas.suction_area_m2 == pytest.approx(numeric, rel=1.0e-9, abs=1.0e-15)


@pytest.mark.parametrize("angle_deg", [45.0, 120.0, 200.0, 300.0])
def test_mirrored_angles_swap_suction_and_discharge(angle_deg: float) -> None:
    geometry = RotaryGeometry.default()
    angle = math.radians(angle_deg)

    direct = chamber_areas(geometry, angle)
    mirrored = chamber_areas(geometry, FULL_TURN - angle)
    assert direct.suction_area_m2 == pytest.approx(mirrored.discharge_area_m2, rel=1.0e-12)
    assert direct.discharge_area_m2 == pytest.approx(mirrored.suction_area_m2, rel=1.0e-12)


def test_suction_grows_from_zero_over_one_revolution() -> None:
    geometry = RotaryGeometry.default()
    half_window = seal_over_half_angle_rad(geometry)

    at_seal = chamber_areas(geometry, half_window)
    assert at_seal.suction_area_m2 == pytest.approx(0.0, abs=1.0e-12)

    before_seal = chamber_areas(geometry, FULL_TURN - half_window)
    assert before_seal.discharge_area_m2 == pytest.approx(0.0, abs=1.0e-12)

    samples = [half_window + 0.05, 0.5 * math.pi, math.pi, 1.5 * math.pi, FULL_TURN - half_window]
    values = [chamber_areas(geometry, angle).suction_area_m2 for angle in samples]
    assert values == sorted(values)
    assert values[0] > 0.0


def test_contact_under_the_vane_is_reported_not_forced() -> None:
    geometry = RotaryGeometry.default()
    half_window = seal_over_half_angle_rad(geometry)

    for angle in (0.0, 0.5 * half_window, FULL_TURN - 0.5 * half_window, -0.25 * half_window):
        with pytest.raises(ChamberSeparationError, match="merge"):
            chamber_areas(geometry, angle)


def test_rotor_cylinder_gap_breaks_separation() -> None:
    geometry = RotaryGeometry(
        cylinder_id_m=77.0 * MM,
        rotor_od_m=68.0 * MM,
        eccentricity_m=4.0 * MM,
        cutout_radius_m=8.0 * MM,
        cutout_offset_m=25.0 * MM,
        vane_width_m=8.0 * MM,
        vane_tip_distance_at_top_m=9.0 * MM,
        frequency_hz=30.0,
    )
    geometry.validate()

    with pytest.raises(ChamberSeparationError, match="gap"):
        chamber_areas(geometry, math.pi)


def test_volumes_default_to_the_supplied_21_mm_height() -> None:
    geometry = RotaryGeometry.default()
    angle = 0.5 * math.pi

    areas = chamber_areas(geometry, angle)
    volumes = chamber_volumes(geometry, angle)
    assert geometry.cylinder_height_m == pytest.approx(21.0 * MM)
    assert volumes.cylinder_height_m == pytest.approx(21.0 * MM)
    assert volumes.suction_volume_m3 == pytest.approx(areas.suction_area_m2 * 21.0 * MM)
    assert volumes.discharge_volume_m3 == pytest.approx(areas.discharge_area_m2 * 21.0 * MM)


def test_volumes_accept_an_explicit_height_override() -> None:
    geometry = RotaryGeometry.default()
    angle = 0.5 * math.pi
    height = 30.0 * MM

    areas = chamber_areas(geometry, angle)
    volumes = chamber_volumes(geometry, angle, height)
    assert volumes.suction_volume_m3 == pytest.approx(areas.suction_area_m2 * height)
    assert volumes.discharge_volume_m3 == pytest.approx(areas.discharge_area_m2 * height)
    assert volumes.cylinder_height_m == pytest.approx(height)


@pytest.mark.parametrize("height_m", [0.0, -1.0e-2, math.inf, math.nan])
def test_nonphysical_height_is_rejected(height_m: float) -> None:
    geometry = RotaryGeometry.default()

    with pytest.raises(ValueError, match="height"):
        chamber_volumes(geometry, 0.5 * math.pi, height_m)


def test_supplied_case_pressures_are_absolute_and_ordered() -> None:
    assert SUCTION_PORT_PRESSURE_PA == pytest.approx(0.82e6)
    assert DISCHARGE_PORT_PRESSURE_PA == pytest.approx(3.24e6)
    assert DEFAULT_POLYTROPIC_EXPONENT == pytest.approx(1.07)
    assert DISCHARGE_VALVE_RISE_FRACTION == pytest.approx(0.05)
    assert RECESS_PRESSURE_PA == pytest.approx(4.0e6)
    assert DISCHARGE_PORT_PRESSURE_PA > SUCTION_PORT_PRESSURE_PA


def test_suction_chamber_always_holds_the_suction_port_pressure() -> None:
    geometry = RotaryGeometry.default()

    for angle in (seal_over_half_angle_rad(geometry), 0.5 * math.pi, math.pi, 1.5 * math.pi):
        pressures = chamber_pressures(geometry, angle)
        assert pressures.suction_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)


def test_compression_starts_continuously_at_the_suction_pressure() -> None:
    geometry = RotaryGeometry.default()

    at_seal = chamber_pressures(geometry, seal_over_half_angle_rad(geometry))
    assert at_seal.discharge_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)
    assert at_seal.seal_over_mixing is False
    assert at_seal.discharge_valve_open is False


def test_compression_follows_the_polytropic_area_ratio() -> None:
    geometry = RotaryGeometry.default()
    angle = 0.5 * math.pi

    sealed = chamber_areas(geometry, seal_over_half_angle_rad(geometry)).discharge_area_m2
    current = chamber_areas(geometry, angle).discharge_area_m2
    expected = SUCTION_PORT_PRESSURE_PA * (sealed / current) ** DEFAULT_POLYTROPIC_EXPONENT
    assert expected < VALVE_OPENING_PA

    pressures = chamber_pressures(geometry, angle)
    assert pressures.discharge_pressure_pa == pytest.approx(expected, rel=1.0e-12)
    assert pressures.discharge_valve_open is False


def test_compression_pressure_rises_monotonically_while_the_valve_is_closed() -> None:
    geometry = RotaryGeometry.default()
    half_window = seal_over_half_angle_rad(geometry)

    angles = [half_window + 0.05, 0.5 * math.pi, math.pi]
    results = [chamber_pressures(geometry, angle) for angle in angles]
    values = [result.discharge_pressure_pa for result in results]
    assert values == sorted(values)
    assert all(result.discharge_valve_open is False for result in results)
    assert all(value < VALVE_OPENING_PA for value in values)


def test_delivery_pressure_declines_linearly_to_the_port_pressure() -> None:
    geometry = RotaryGeometry.default()
    half_window = seal_over_half_angle_rad(geometry)
    window_entry = 2.0 * math.pi - half_window

    early = chamber_pressures(geometry, 4.2)
    middle = chamber_pressures(geometry, 1.5 * math.pi)
    late = chamber_pressures(geometry, window_entry - 1.0e-4)
    for result in (early, middle, late):
        assert result.discharge_valve_open is True
        assert result.seal_over_mixing is False
        assert result.valve_opening_pressure_pa == pytest.approx(VALVE_OPENING_PA)
    assert VALVE_OPENING_PA >= early.discharge_pressure_pa
    assert early.discharge_pressure_pa > middle.discharge_pressure_pa
    assert middle.discharge_pressure_pa > late.discharge_pressure_pa
    assert late.discharge_pressure_pa == pytest.approx(DISCHARGE_PORT_PRESSURE_PA, rel=1.0e-3)


def test_delivery_pressure_is_continuous_at_the_valve_opening() -> None:
    geometry = RotaryGeometry.default()

    flip_angle = None
    angle = 3.7
    while angle < 4.1:
        if chamber_pressures(geometry, angle).discharge_valve_open:
            flip_angle = angle
            break
        angle += 1.0e-3
    assert flip_angle is not None
    before = chamber_pressures(geometry, flip_angle - 1.0e-3)
    after = chamber_pressures(geometry, flip_angle)
    assert before.discharge_pressure_pa == pytest.approx(VALVE_OPENING_PA, rel=2.0e-3)
    assert after.discharge_pressure_pa == pytest.approx(VALVE_OPENING_PA, rel=2.0e-3)


def test_zero_valve_rise_caps_at_the_port_pressure() -> None:
    geometry = RotaryGeometry.default()

    late = chamber_pressures(geometry, 1.5 * math.pi, valve_rise_fraction=0.0)
    assert late.discharge_pressure_pa == pytest.approx(DISCHARGE_PORT_PRESSURE_PA)
    assert late.discharge_valve_open is True


@pytest.mark.parametrize("angle_factor", [0.0, 0.5])
def test_leaving_seal_over_window_holds_the_suction_pressure(angle_factor: float) -> None:
    geometry = RotaryGeometry.default()
    angle = angle_factor * seal_over_half_angle_rad(geometry)

    assert in_seal_over_window(geometry, angle)
    pressures = chamber_pressures(geometry, angle)
    assert pressures.suction_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)
    assert pressures.discharge_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)
    assert pressures.seal_over_mixing is True
    assert pressures.discharge_valve_open is False


def test_entering_seal_over_window_ramps_linearly_to_the_suction_pressure() -> None:
    geometry = RotaryGeometry.default()
    half_window = seal_over_half_angle_rad(geometry)

    halfway = chamber_pressures(geometry, 2.0 * math.pi - 0.5 * half_window)
    expected = DISCHARGE_PORT_PRESSURE_PA + 0.5 * (
        SUCTION_PORT_PRESSURE_PA - DISCHARGE_PORT_PRESSURE_PA
    )
    assert halfway.seal_over_mixing is True
    assert halfway.discharge_valve_open is False
    assert halfway.suction_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)
    assert halfway.discharge_pressure_pa == pytest.approx(expected, rel=1.0e-12)

    near_tdc = chamber_pressures(geometry, 2.0 * math.pi - 1.0e-5)
    assert near_tdc.discharge_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA, rel=1.0e-2)


def test_pressures_reject_a_rotor_cylinder_gap() -> None:
    geometry = RotaryGeometry(
        cylinder_id_m=77.0 * MM,
        rotor_od_m=68.0 * MM,
        eccentricity_m=4.0 * MM,
        cutout_radius_m=8.0 * MM,
        cutout_offset_m=25.0 * MM,
        vane_width_m=8.0 * MM,
        vane_tip_distance_at_top_m=9.0 * MM,
        frequency_hz=30.0,
    )

    with pytest.raises(ChamberSeparationError, match="gap"):
        chamber_pressures(geometry, math.pi)


@pytest.mark.parametrize("pressure_pa", [0.0, -1.0e5, math.inf, math.nan])
def test_nonphysical_port_pressures_are_rejected(pressure_pa: float) -> None:
    geometry = RotaryGeometry.default()

    with pytest.raises(ValueError, match="[Pp]ressure"):
        chamber_pressures(geometry, 0.5 * math.pi, pressure_pa, 9.0e5)
    with pytest.raises(ValueError, match="[Pp]ressure"):
        chamber_pressures(geometry, 0.5 * math.pi, 1.0e5, pressure_pa)


def test_pressure_ordering_and_exponent_are_validated() -> None:
    geometry = RotaryGeometry.default()

    with pytest.raises(ValueError, match="exceed"):
        chamber_pressures(geometry, 0.5 * math.pi, 3.24e6, 0.82e6)
    for exponent in (0.0, -1.3, math.inf, math.nan):
        with pytest.raises(ValueError, match="[Ee]xponent"):
            chamber_pressures(geometry, 0.5 * math.pi, polytropic_exponent=exponent)
    for fraction in (-0.05, math.inf, math.nan):
        with pytest.raises(ValueError, match="fraction"):
            chamber_pressures(geometry, 0.5 * math.pi, valve_rise_fraction=fraction)


# --- Port-timed cycle on true-geometry volumes (PHYSICS.md section 3.4) ---

TRACE_SAMPLES = 72
TRACE_STEP_M = 2.0e-4


@pytest.fixture(scope="module")
def cycle_trace() -> CycleTrace:
    """One shared trace: each sample runs an area integration."""

    return build_cycle_trace(
        RotaryGeometry.default(), samples=TRACE_SAMPLES, grid_step_m=TRACE_STEP_M
    )


@pytest.mark.parametrize(
    "angle_deg, expected_phase",
    [
        (2.0, "suction-connected"),
        (20.0, "suction-connected"),
        (60.0, "compression"),
        (180.0, "compression"),
        (300.0, "discharge"),
        (345.0, "discharge"),
        (348.0, "recompression"),
        (352.0, "recompression"),
    ],
)
def test_port_timing_selects_the_expected_phase(
    cycle_trace: CycleTrace, angle_deg: float, expected_phase: str
) -> None:
    result = port_timed_pressures(
        RotaryGeometry.default(), math.radians(angle_deg), trace=cycle_trace
    )

    assert result.phase == expected_phase


def test_the_discharge_valve_opens_from_a_pressure_condition(
    cycle_trace: CycleTrace,
) -> None:
    """The opening angle must move with the operating condition.

    Raising the discharge pressure has to delay the opening; a hardcoded
    angle would not.
    """

    geometry = RotaryGeometry.default()

    def opening_angle_deg(discharge_pa: float) -> float:
        for step in range(1, 3600):
            angle_deg = 0.1 * step
            result = port_timed_pressures(
                geometry,
                math.radians(angle_deg),
                trace=cycle_trace,
                discharge_port_pressure_pa=discharge_pa,
            )
            if result.discharge_valve_open:
                return angle_deg
        raise AssertionError("the valve never opened")

    assert opening_angle_deg(3.6e6) > opening_angle_deg(3.24e6)


def test_suction_holds_the_port_pressure_once_the_port_opens(
    cycle_trace: CycleTrace,
) -> None:
    """The machine has no suction valve, so no suction loss is modelled."""

    geometry = RotaryGeometry.default()

    for angle_deg in (10.4, 90.0, 180.0, 270.0, 350.0):
        result = port_timed_pressures(geometry, math.radians(angle_deg), trace=cycle_trace)
        assert result.suction_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)


def test_re_expansion_falls_to_the_suction_pressure_before_the_port_opens(
    cycle_trace: CycleTrace,
) -> None:
    """Residual gas expands from the trapped pressure down to suction."""

    geometry = RotaryGeometry.default()

    early = port_timed_pressures(geometry, math.radians(1.0), trace=cycle_trace)
    late = port_timed_pressures(geometry, math.radians(9.0), trace=cycle_trace)
    opened = port_timed_pressures(geometry, math.radians(10.4), trace=cycle_trace)

    assert early.suction_pressure_pa > late.suction_pressure_pa
    assert late.suction_pressure_pa > SUCTION_PORT_PRESSURE_PA
    assert opened.suction_pressure_pa == pytest.approx(SUCTION_PORT_PRESSURE_PA)


def test_recompression_is_finite_and_overshoots_the_discharge_pressure(
    cycle_trace: CycleTrace,
) -> None:
    """The over-compression spike the circular approximation could not give.

    With no leakage this peak is an upper bound, not a prediction.
    """

    geometry = RotaryGeometry.default()

    start = port_timed_pressures(geometry, math.radians(346.8), trace=cycle_trace)
    peak = port_timed_pressures(geometry, math.radians(354.0), trace=cycle_trace)

    assert start.compression_pressure_pa == pytest.approx(DISCHARGE_PORT_PRESSURE_PA, rel=1e-3)
    assert math.isfinite(peak.compression_pressure_pa)
    assert peak.compression_pressure_pa > DISCHARGE_PORT_PRESSURE_PA


def test_the_residual_pressure_closes_the_cycle(cycle_trace: CycleTrace) -> None:
    """The recompression end pressure is what the next re-expansion starts from.

    It is derived from the geometry and the port pressures, never hardcoded.
    """

    geometry = RotaryGeometry.default()

    result = port_timed_pressures(geometry, 0.0, trace=cycle_trace)
    window_entry = FULL_TURN_RAD - seal_over_half_angle_rad(geometry)
    at_entry = port_timed_pressures(geometry, window_entry, trace=cycle_trace)

    assert result.residual_pressure_pa == pytest.approx(at_entry.compression_pressure_pa)
    assert result.residual_pressure_pa > DISCHARGE_PORT_PRESSURE_PA


def test_moving_the_compression_start_barely_moves_the_valve_angle() -> None:
    """The compression chamber volume is flat near top dead center.

    So the reading of beta, which is not settled, cannot change the
    compression stroke appreciably.
    """

    geometry = RotaryGeometry.default()
    sealed_early = chamber_areas(geometry, seal_over_half_angle_rad(geometry)).discharge_area_m2
    sealed_late = chamber_areas(geometry, math.radians(27.7)).discharge_area_m2

    def valve_angle_deg(sealed_area_m2: float) -> float:
        return math.degrees(
            _valve_opening_angle_rad(
                geometry,
                SUCTION_PORT_PRESSURE_PA,
                DISCHARGE_PORT_PRESSURE_PA * (1.0 + DISCHARGE_VALVE_RISE_FRACTION),
                DEFAULT_POLYTROPIC_EXPONENT,
                sealed_area_m2,
                seal_over_half_angle_rad(geometry),
                FULL_TURN_RAD - seal_over_half_angle_rad(geometry),
            )
        )

    assert abs(valve_angle_deg(sealed_early) - valve_angle_deg(sealed_late)) < 0.5


@pytest.mark.parametrize("samples", [0, 7])
def test_a_cycle_trace_needs_enough_samples(samples: int) -> None:
    with pytest.raises(ValueError, match="samples"):
        build_cycle_trace(RotaryGeometry.default(), samples=samples)
