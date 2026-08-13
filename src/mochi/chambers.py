"""Chamber cross-section areas for the prescribed rotary mechanism.

The rotor is approximated by its full outside-diameter disc; the rotor slot
and circular cutout are ignored and the drawn geometry is unchanged. The
crescent between the cylinder bore and this disc is divided by the union of
the vane and the rotor disc into a suction chamber on the inlet (+x) side and
a discharge chamber on the outlet (-x) side. Volumes multiply a cross-section
area by the axial cylinder height, and the accepted constant chamber pressure
rule assigns each chamber its port pressure with seal-over mixing. See
PHYSICS.md sections 3.1 and 3.2 for the equations and validity conditions.
"""

from dataclasses import dataclass
from functools import lru_cache
from math import asin, degrees, hypot, inf, isfinite, pi, sin, sqrt

from mochi.chamber_volume import (
    DEFAULT_GRID_STEP_M,
    AxialBands,
    SwingBush,
    true_chamber_volumes,
)
from mochi.kinematics import PrescribedState, RotaryGeometry, prescribed_state
from mochi.ports import characteristic_angles

FULL_TURN_RAD = 2.0 * pi

# Supplied absolute R410A port pressures (2026-07-16); they correspond to
# 0.85 degC evaporating and 52.5 degC condensing saturation temperatures.
SUCTION_PORT_PRESSURE_PA = 0.82e6
DISCHARGE_PORT_PRESSURE_PA = 3.24e6

# Effective R410A polytropic exponent for p V**n = const along the real-gas
# isentrope from saturated vapour at 0.82 MPa to 3.24 MPa (CoolProp HEOS:
# 1.064 saturated to 1.075 at 10 K superheat).
DEFAULT_POLYTROPIC_EXPONENT = 1.07

# Discharge reed valves open only above the line pressure; the accepted
# approximation raises the delivery pressure by this fraction of the
# discharge-port pressure (typical small hermetic compressor loss 2-6%).
DISCHARGE_VALVE_RISE_FRACTION = 0.05

# Fixed absolute pressure of the sealed recess spaces above and below the
# rotor's thin 14.7 mm channel, held by the full-height rotor, the swing
# bush, and the vane (PHYSICS.md section 3.3). Not used in any computation
# yet; revisit with the swing-bush pressure boundary conditions.
RECESS_PRESSURE_PA = 4.0e6

# Up to one nanometre of rotor-to-cylinder gap is treated as the ideal
# tangency seal; a larger gap connects the two chambers around the rotor.
SEALING_GAP_TOLERANCE_M = 1.0e-9

# Negative round-off no larger than this fraction of the crescent area is
# clamped to zero; a more negative area indicates a formula error.
RELATIVE_AREA_TOLERANCE = 1.0e-9


class ChamberSeparationError(ValueError):
    """The vane and rotor disc do not split the crescent into two chambers."""


@dataclass(frozen=True, slots=True)
class ChamberAreas:
    """Suction/discharge cross-section split at one crank angle."""

    crank_angle_rad: float
    suction_area_m2: float
    discharge_area_m2: float
    vane_area_in_crescent_m2: float
    crescent_area_m2: float


@dataclass(frozen=True, slots=True)
class ChamberVolumes:
    """Chamber volumes for the axial cylinder height."""

    crank_angle_rad: float
    cylinder_height_m: float
    suction_volume_m3: float
    discharge_volume_m3: float


@dataclass(frozen=True, slots=True)
class ChamberPressures:
    """Prescribed chamber pressures at one crank angle (absolute pascals)."""

    crank_angle_rad: float
    suction_pressure_pa: float
    discharge_pressure_pa: float
    valve_opening_pressure_pa: float
    seal_over_mixing: bool
    discharge_valve_open: bool


def seal_over_half_angle_rad(geometry: RotaryGeometry) -> float:
    """Half-width of the crank window where the contact hides under the vane.

    Inside ``(-half, +half)`` around crank angle zero the rotor-cylinder
    contact point lies within the vane width, the vane covers the tangency
    region, and the two chambers form one connected region. This is the
    physical seal-over instant where the finished suction chamber becomes the
    next compression chamber.
    """

    return asin(0.5 * geometry.vane_width_m / geometry.cylinder_radius_m)


def _circle_area_integral(x_m: float, radius_m: float) -> float:
    """Antiderivative of sqrt(radius**2 - x**2) used for circle-bounded strips."""

    # Clipping guards round-off at the circle edge only; callers keep the
    # argument inside the circle through the separation checks.
    clipped_ratio = min(max(x_m / radius_m, -1.0), 1.0)
    return 0.5 * (x_m * sqrt(max(radius_m**2 - x_m**2, 0.0)) + radius_m**2 * asin(clipped_ratio))


def _swept_crescent_area(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Crescent area swept clockwise from the vane centreline to the contact.

    Closed form of the integral of (R_c**2 - rho(phi)**2) / 2 for rays from
    the cylinder centre, where rho(phi) is the distance to the rotor-disc
    boundary and the rotor sits at the same crank angle, so the sweep ends at
    the rotor-cylinder tangency point.
    """

    eccentricity = geometry.eccentricity_m
    rotor_radius = geometry.rotor_radius_m
    cylinder_radius = geometry.cylinder_radius_m
    chord = eccentricity * sin(crank_angle_rad)
    return (
        0.5 * (cylinder_radius**2 - rotor_radius**2) * crank_angle_rad
        - 0.25 * eccentricity**2 * sin(2.0 * crank_angle_rad)
        - 0.5 * chord * sqrt(rotor_radius**2 - chord**2)
        - 0.5 * rotor_radius**2 * asin(chord / rotor_radius)
    )


def _vane_strip_area(
    geometry: RotaryGeometry,
    state: PrescribedState,
    inlet_side: bool,
) -> float:
    """Crescent area occupied by the vane on one side of its centreline.

    Integrates cylinder-arc minus upper-rotor-arc height over half the vane
    width. Valid because the separation checks keep the vane bottom inside
    the rotor disc, so the vane covers this strip completely.
    """

    half_width = 0.5 * geometry.vane_width_m
    rotor_x, rotor_y = state.rotor_center_m
    lower = 0.0 if inlet_side else -half_width
    upper = half_width if inlet_side else 0.0
    cylinder_part = _circle_area_integral(upper, geometry.cylinder_radius_m)
    cylinder_part -= _circle_area_integral(lower, geometry.cylinder_radius_m)
    rotor_part = _circle_area_integral(upper - rotor_x, geometry.rotor_radius_m)
    rotor_part -= _circle_area_integral(lower - rotor_x, geometry.rotor_radius_m)
    return cylinder_part - rotor_part - half_width * rotor_y


def in_seal_over_window(geometry: RotaryGeometry, crank_angle_rad: float) -> bool:
    """Return True when the rotor-cylinder contact hides under the vane."""

    normalized = crank_angle_rad % FULL_TURN_RAD
    half_window = seal_over_half_angle_rad(geometry)
    return normalized < half_window or normalized > FULL_TURN_RAD - half_window


def _require_sealing_geometry(geometry: RotaryGeometry, state: PrescribedState) -> None:
    """Raise ``ChamberSeparationError`` when the geometry can never seal."""

    gap_m = geometry.radial_clearance_m - geometry.eccentricity_m
    if gap_m > SEALING_GAP_TOLERANCE_M:
        raise ChamberSeparationError(
            "chambers stay connected: the rotor-cylinder gap is "
            f"{gap_m * 1.0e3:.6g} mm; the circular-rotor approximation needs "
            "tangency, e = (cylinder ID - rotor OD) / 2"
        )
    if geometry.eccentricity_m >= geometry.rotor_radius_m:
        raise ChamberSeparationError(
            "eccentricity reaches the rotor radius; the cylinder centre must "
            "stay inside the approximated rotor disc"
        )
    rotor_x, rotor_y = state.rotor_center_m
    vane_tip_y = state.vane_tip_m[1]
    half_width = 0.5 * geometry.vane_width_m
    for corner_x in (half_width, -half_width):
        if hypot(corner_x - rotor_x, vane_tip_y - rotor_y) >= geometry.rotor_radius_m:
            raise ChamberSeparationError(
                "chambers stay connected: the vane bottom corner at "
                f"x = {corner_x * 1.0e3:+.3g} mm does not reach inside the "
                "approximated circular rotor"
            )


def _require_separated_chambers(geometry: RotaryGeometry, state: PrescribedState) -> None:
    """Raise ``ChamberSeparationError`` unless exactly two chambers exist."""

    _require_sealing_geometry(geometry, state)
    if in_seal_over_window(geometry, state.crank_angle_rad):
        half_window = seal_over_half_angle_rad(geometry)
        raise ChamberSeparationError(
            "chambers merge into one: the rotor-cylinder contact lies under "
            f"the vane within +/-{degrees(half_window):.2f} deg of crank angle "
            "zero (suction seal-over)"
        )


def _clamped_area(area_m2: float, crescent_area_m2: float) -> float:
    if area_m2 < -RELATIVE_AREA_TOLERANCE * crescent_area_m2:
        raise ValueError(f"internal chamber-area inconsistency: computed {area_m2} m^2")
    return max(area_m2, 0.0)


def chamber_areas(geometry: RotaryGeometry, crank_angle_rad: float) -> ChamberAreas:
    """Split the crescent cross-section at one crank angle.

    The suction chamber runs clockwise from the vane's inlet-side face to the
    rotor-cylinder contact; the discharge chamber continues clockwise back to
    the vane's outlet-side face. Raises ``ChamberSeparationError`` when the
    union of the vane and the approximated circular rotor does not divide the
    crescent into exactly these two regions.
    """

    state = prescribed_state(geometry, crank_angle_rad)
    _require_separated_chambers(geometry, state)

    crescent = pi * (geometry.cylinder_radius_m**2 - geometry.rotor_radius_m**2)
    normalized = crank_angle_rad % FULL_TURN_RAD
    swept = _swept_crescent_area(geometry, normalized)
    inlet_strip = _vane_strip_area(geometry, state, inlet_side=True)
    outlet_strip = _vane_strip_area(geometry, state, inlet_side=False)
    return ChamberAreas(
        crank_angle_rad=crank_angle_rad,
        suction_area_m2=_clamped_area(swept - inlet_strip, crescent),
        discharge_area_m2=_clamped_area(crescent - swept - outlet_strip, crescent),
        vane_area_in_crescent_m2=inlet_strip + outlet_strip,
        crescent_area_m2=crescent,
    )


def chamber_volumes(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    cylinder_height_m: float | None = None,
) -> ChamberVolumes:
    """Multiply the chamber cross-sections by the axial cylinder height.

    The supplied 21 mm height stored on the geometry is used unless a
    different height is passed explicitly.
    """

    height_m = geometry.cylinder_height_m if cylinder_height_m is None else cylinder_height_m
    if not isfinite(height_m) or height_m <= 0.0:
        raise ValueError("Cylinder height must be a positive, finite length in metres.")
    areas = chamber_areas(geometry, crank_angle_rad)
    return ChamberVolumes(
        crank_angle_rad=crank_angle_rad,
        cylinder_height_m=height_m,
        suction_volume_m3=areas.suction_area_m2 * height_m,
        discharge_volume_m3=areas.discharge_area_m2 * height_m,
    )


def _compressed_pressure_pa(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    suction_port_pressure_pa: float,
    polytropic_exponent: float,
    sealed_area_m2: float,
) -> float:
    current_area_m2 = chamber_areas(geometry, crank_angle_rad).discharge_area_m2
    if current_area_m2 <= 0.0:
        return inf
    return float(
        suction_port_pressure_pa * (sealed_area_m2 / current_area_m2) ** polytropic_exponent
    )


def _valve_opening_angle_rad(
    geometry: RotaryGeometry,
    suction_port_pressure_pa: float,
    opening_pressure_pa: float,
    polytropic_exponent: float,
    sealed_area_m2: float,
    low_rad: float,
    high_rad: float,
) -> float:
    """Bisect the crank angle where the polytropic pressure reaches the valve."""

    for _ in range(60):
        mid = 0.5 * (low_rad + high_rad)
        mid_pa = _compressed_pressure_pa(
            geometry, mid, suction_port_pressure_pa, polytropic_exponent, sealed_area_m2
        )
        if mid_pa < opening_pressure_pa:
            low_rad = mid
        else:
            high_rad = mid
    return 0.5 * (low_rad + high_rad)


def chamber_pressures(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    suction_port_pressure_pa: float = SUCTION_PORT_PRESSURE_PA,
    discharge_port_pressure_pa: float = DISCHARGE_PORT_PRESSURE_PA,
    polytropic_exponent: float = DEFAULT_POLYTROPIC_EXPONENT,
    valve_rise_fraction: float = DISCHARGE_VALVE_RISE_FRACTION,
) -> ChamberPressures:
    """Assign the accepted chamber pressure boundary conditions.

    All pressures are absolute. The suction chamber always takes the
    suction-port pressure. The compression (discharge) chamber starts at the
    suction-port pressure when it seals over and follows the polytropic
    relation p V**n = const (the volume ratio equals the cross-section area
    ratio at constant height) until it reaches the valve opening pressure
    (the discharge-port pressure raised by the valve rise fraction). During
    delivery the pressure then falls linearly in crank angle from the
    opening pressure to the discharge-port pressure at the seal-over window
    entry, approximating the decaying valve/flow overpressure. Across the
    entering half of the seal-over window the merged region follows a linear
    mixing ramp from the discharge-port pressure down to the suction-port
    pressure, reaching it exactly at 360 degrees; the leaving half of the
    window stays at the suction-port pressure. The discharge-chamber trace
    is therefore continuous over the whole revolution (PHYSICS.md section
    3.2). These ramps are prescribed boundary conditions, not valve or flow
    models; leakage, heat transfer, and clearance re-expansion are
    neglected, and pressure forces are not yet implemented. Raises
    ``ChamberSeparationError`` when the geometry cannot seal at any angle
    (for example a rotor-cylinder gap).
    """

    for pressure_pa in (suction_port_pressure_pa, discharge_port_pressure_pa):
        if not isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise ValueError(
                "Port pressures must be positive, finite absolute pressures in pascals."
            )
    if discharge_port_pressure_pa <= suction_port_pressure_pa:
        raise ValueError("Discharge-port pressure must exceed the suction-port pressure.")
    if not isfinite(polytropic_exponent) or polytropic_exponent <= 0.0:
        raise ValueError("Polytropic exponent must be a positive, finite number.")
    if not isfinite(valve_rise_fraction) or valve_rise_fraction < 0.0:
        raise ValueError("Valve rise fraction must be a non-negative, finite number.")
    opening_pressure_pa = discharge_port_pressure_pa * (1.0 + valve_rise_fraction)
    state = prescribed_state(geometry, crank_angle_rad)
    _require_sealing_geometry(geometry, state)

    half_window = seal_over_half_angle_rad(geometry)
    window_entry = FULL_TURN_RAD - half_window
    normalized = crank_angle_rad % FULL_TURN_RAD
    if normalized < half_window:
        # leaving half of the seal-over window: mixing is complete
        return ChamberPressures(
            crank_angle_rad=crank_angle_rad,
            suction_pressure_pa=suction_port_pressure_pa,
            discharge_pressure_pa=suction_port_pressure_pa,
            valve_opening_pressure_pa=opening_pressure_pa,
            seal_over_mixing=True,
            discharge_valve_open=False,
        )
    if normalized > window_entry:
        # entering half: linear mixing ramp from the port pressure to suction
        fraction = (normalized - window_entry) / half_window
        mixed_pa = (
            discharge_port_pressure_pa
            + (suction_port_pressure_pa - discharge_port_pressure_pa) * fraction
        )
        return ChamberPressures(
            crank_angle_rad=crank_angle_rad,
            suction_pressure_pa=suction_port_pressure_pa,
            discharge_pressure_pa=mixed_pa,
            valve_opening_pressure_pa=opening_pressure_pa,
            seal_over_mixing=True,
            discharge_valve_open=False,
        )
    sealed_area_m2 = chamber_areas(geometry, half_window).discharge_area_m2
    compressed_pa = _compressed_pressure_pa(
        geometry, normalized, suction_port_pressure_pa, polytropic_exponent, sealed_area_m2
    )
    if compressed_pa < opening_pressure_pa:
        return ChamberPressures(
            crank_angle_rad=crank_angle_rad,
            suction_pressure_pa=suction_port_pressure_pa,
            discharge_pressure_pa=compressed_pa,
            valve_opening_pressure_pa=opening_pressure_pa,
            seal_over_mixing=False,
            discharge_valve_open=False,
        )
    valve_angle = _valve_opening_angle_rad(
        geometry,
        suction_port_pressure_pa,
        opening_pressure_pa,
        polytropic_exponent,
        sealed_area_m2,
        half_window,
        normalized,
    )
    fraction = (normalized - valve_angle) / (window_entry - valve_angle)
    delivery_pa = (
        opening_pressure_pa + (discharge_port_pressure_pa - opening_pressure_pa) * fraction
    )
    return ChamberPressures(
        crank_angle_rad=crank_angle_rad,
        suction_pressure_pa=suction_port_pressure_pa,
        discharge_pressure_pa=delivery_pa,
        valve_opening_pressure_pa=opening_pressure_pa,
        seal_over_mixing=False,
        discharge_valve_open=True,
    )


# --- Port-timed cycle on true-geometry volumes (PHYSICS.md section 3.4) ---
#
# The rule above is the accepted constant-pressure boundary condition of
# section 3.2 and stays the regression baseline. What follows adds the phase
# structure the supplied port angles make possible: a re-expansion phase
# before the suction port opens and a recompression phase after the discharge
# port closes. Both need a clearance volume that the circular-rotor
# approximation cannot provide, so these run on the true-geometry volumes.

CYCLE_TRACE_SAMPLES = 121


@dataclass(frozen=True, slots=True)
class CycleTrace:
    """Sampled chamber volumes over one revolution.

    The volume integration is far too slow to run per query, so a trace is
    built once per geometry and interpolated. Angles are equally spaced and
    ascending over ``[0, 2*pi)``.
    """

    angles_rad: tuple[float, ...]
    suction_volume_m3: tuple[float, ...]
    compression_volume_m3: tuple[float, ...]

    def suction_at(self, crank_angle_rad: float) -> float:
        """Interpolate the suction-chamber volume."""

        return _interpolate(self.angles_rad, self.suction_volume_m3, crank_angle_rad)

    def compression_at(self, crank_angle_rad: float) -> float:
        """Interpolate the compression-chamber volume."""

        return _interpolate(self.angles_rad, self.compression_volume_m3, crank_angle_rad)


@dataclass(frozen=True, slots=True)
class PortTimedPressures:
    """Chamber pressures from the port-timed phase rule.

    ``phase`` names the branch that produced ``compression_pressure_pa`` so
    callers and tests can assert the phase boundaries directly.
    """

    crank_angle_rad: float
    suction_pressure_pa: float
    compression_pressure_pa: float
    phase: str
    valve_opening_pressure_pa: float
    residual_pressure_pa: float
    discharge_valve_open: bool
    seal_over_mixing: bool


def _interpolate(
    nodes_rad: tuple[float, ...],
    values: tuple[float, ...],
    angle_rad: float,
) -> float:
    """Linear interpolation on the periodic sampled trace."""

    normalized = angle_rad % FULL_TURN_RAD
    spacing = FULL_TURN_RAD / len(nodes_rad)
    lower = int(normalized / spacing)
    fraction = (normalized - lower * spacing) / spacing
    upper = (lower + 1) % len(nodes_rad)
    return values[lower] + (values[upper] - values[lower]) * fraction


@lru_cache(maxsize=8)
def build_cycle_trace(
    geometry: RotaryGeometry,
    samples: int = CYCLE_TRACE_SAMPLES,
    grid_step_m: float = DEFAULT_GRID_STEP_M,
    bush: SwingBush | None = None,
    bands: AxialBands | None = None,
) -> CycleTrace:
    """Sample the true-geometry chamber volumes over one revolution.

    Cached because each sample runs an area integration; a full trace costs
    seconds, so callers must not rebuild it per crank angle.
    """

    if samples < 8:
        raise ValueError("A cycle trace needs at least eight samples.")
    angles = tuple(FULL_TURN_RAD * index / samples for index in range(samples))
    suction: list[float] = []
    compression: list[float] = []
    for angle in angles:
        volumes = true_chamber_volumes(
            geometry, angle, bush=bush, bands=bands, grid_step_m=grid_step_m
        )
        suction.append(volumes.suction_volume_m3)
        compression.append(volumes.discharge_volume_m3)
    return CycleTrace(
        angles_rad=angles,
        suction_volume_m3=tuple(suction),
        compression_volume_m3=tuple(compression),
    )


def _polytropic_pa(
    reference_pressure_pa: float,
    reference_volume_m3: float,
    volume_m3: float,
    polytropic_exponent: float,
) -> float:
    """Return ``p V**n = const`` with a guard against a vanishing volume."""

    if volume_m3 <= 0.0:
        return inf
    return float(reference_pressure_pa * (reference_volume_m3 / volume_m3) ** polytropic_exponent)


def port_timed_pressures(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    trace: CycleTrace | None = None,
    suction_port_pressure_pa: float = SUCTION_PORT_PRESSURE_PA,
    discharge_port_pressure_pa: float = DISCHARGE_PORT_PRESSURE_PA,
    polytropic_exponent: float = DEFAULT_POLYTROPIC_EXPONENT,
    valve_rise_fraction: float = DISCHARGE_VALVE_RISE_FRACTION,
) -> PortTimedPressures:
    """Assign chamber pressures from the supplied port timing.

    Phases follow the swing-compressor rule set, all on the instantaneous
    view over ``[0, 2*pi)`` with the characteristic angles of
    :func:`mochi.ports.characteristic_angles`:

    * ``[0, phi)`` re-expansion: the new suction chamber is sealed from both
      ports and its residual gas expands from the previous cycle's trapped
      pressure, clamped once it reaches the suction port.
    * ``[phi, 2*pi)`` suction at the suction-port pressure. The machine has
      no suction valve, so no suction pressure loss is modelled.
    * ``[beta, theta_vo)`` compression: polytropic from the volume at
      ``beta``, the angle where the chamber loses the suction port.
    * ``[theta_vo, 2*pi - delta)`` delivery at the discharge-port pressure.
      ``theta_vo`` is found from the pressure condition, never hardcoded, so
      it moves when the operating condition moves.
    * ``[2*pi - delta, 2*pi - alpha)`` recompression: the discharge port is
      shut and the trapped gas is compressed further. This is the phase the
      circular-rotor approximation cannot represent, because it drives the
      trapped volume to zero.

    The recompression phase ends at the seal-over window entry, where the
    two chambers merge and no separate compression chamber exists. Its end
    pressure is the residual that opens the next cycle's re-expansion, so
    the loop closes without hardcoding a residual pressure.

    Leakage is neglected, which makes the recompression peak an upper bound.
    """

    for pressure_pa in (suction_port_pressure_pa, discharge_port_pressure_pa):
        if not isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise ValueError(
                "Port pressures must be positive, finite absolute pressures in pascals."
            )
    if discharge_port_pressure_pa <= suction_port_pressure_pa:
        raise ValueError("Discharge-port pressure must exceed the suction-port pressure.")
    if not isfinite(polytropic_exponent) or polytropic_exponent <= 0.0:
        raise ValueError("Polytropic exponent must be a positive, finite number.")
    if not isfinite(valve_rise_fraction) or valve_rise_fraction < 0.0:
        raise ValueError("Valve rise fraction must be a non-negative, finite number.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    angles = characteristic_angles(geometry)
    opening_pressure_pa = discharge_port_pressure_pa * (1.0 + valve_rise_fraction)
    half_window = seal_over_half_angle_rad(geometry)
    window_entry = FULL_TURN_RAD - half_window
    normalized = crank_angle_rad % FULL_TURN_RAD

    # Recompression closes the loop: the trapped gas at the discharge close
    # angle is compressed to the window entry, and that is the residual the
    # next re-expansion starts from.
    trapped_volume_m3 = trace.compression_at(angles.discharge_close_rad)
    residual_pa = _polytropic_pa(
        discharge_port_pressure_pa,
        trapped_volume_m3,
        trace.compression_at(window_entry),
        polytropic_exponent,
    )

    sealed_volume_m3 = trace.compression_at(angles.compression_start_rad)
    valve_angle_rad = _port_timed_valve_angle_rad(
        trace,
        angles.compression_start_rad,
        angles.discharge_close_rad,
        sealed_volume_m3,
        suction_port_pressure_pa,
        opening_pressure_pa,
        polytropic_exponent,
    )

    suction_pa = suction_port_pressure_pa
    if normalized < angles.suction_open_rad:
        expanded_pa = _polytropic_pa(
            residual_pa,
            trace.compression_at(window_entry),
            trace.suction_at(normalized),
            polytropic_exponent,
        )
        suction_pa = max(min(expanded_pa, residual_pa), suction_port_pressure_pa)

    if normalized >= window_entry:
        phase = "recompression"
        compression_pa = _polytropic_pa(
            discharge_port_pressure_pa,
            trapped_volume_m3,
            trace.compression_at(window_entry),
            polytropic_exponent,
        )
        merged = True
        valve_open = False
    elif normalized >= angles.discharge_close_rad:
        phase = "recompression"
        compression_pa = _polytropic_pa(
            discharge_port_pressure_pa,
            trapped_volume_m3,
            trace.compression_at(normalized),
            polytropic_exponent,
        )
        merged = False
        valve_open = False
    elif normalized >= valve_angle_rad:
        phase = "discharge"
        compression_pa = discharge_port_pressure_pa
        merged = False
        valve_open = True
    elif normalized >= angles.compression_start_rad:
        phase = "compression"
        compression_pa = _polytropic_pa(
            suction_port_pressure_pa,
            sealed_volume_m3,
            trace.compression_at(normalized),
            polytropic_exponent,
        )
        merged = False
        valve_open = False
    else:
        # Before beta the chamber still sees the suction port.
        phase = "suction-connected"
        compression_pa = suction_port_pressure_pa
        merged = normalized < half_window
        valve_open = False

    return PortTimedPressures(
        crank_angle_rad=crank_angle_rad,
        suction_pressure_pa=suction_pa,
        compression_pressure_pa=compression_pa,
        phase=phase,
        valve_opening_pressure_pa=opening_pressure_pa,
        residual_pressure_pa=residual_pa,
        discharge_valve_open=valve_open,
        seal_over_mixing=merged,
    )


def _port_timed_valve_angle_rad(
    trace: CycleTrace,
    low_rad: float,
    high_rad: float,
    sealed_volume_m3: float,
    suction_port_pressure_pa: float,
    opening_pressure_pa: float,
    polytropic_exponent: float,
) -> float:
    """Bisect the angle where compression reaches the valve opening pressure.

    Kept a pressure condition rather than a constant so it tracks the
    operating condition.
    """

    for _ in range(60):
        mid = 0.5 * (low_rad + high_rad)
        mid_pa = _polytropic_pa(
            suction_port_pressure_pa,
            sealed_volume_m3,
            trace.compression_at(mid),
            polytropic_exponent,
        )
        if mid_pa < opening_pressure_pa:
            low_rad = mid
        else:
            high_rad = mid
    return 0.5 * (low_rad + high_rad)
