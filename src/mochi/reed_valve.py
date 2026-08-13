"""Discharge reed valve on the chamber mass balance (quasi-static check valve).

The port-timed rule (Section 3.4) and the leakage model (Section 3.7) hold the
compression chamber at the discharge-port pressure once the delivery phase
starts. A real discharge is gated by a **reed valve**: it lifts only when the
chamber pressure exceeds the discharge line, and it flows through the geometric
port opening, which the rotor-cylinder contact closes over the last few degrees.
This module adds that valve to the leakage mass balance, letting the chamber
pressure **float** during delivery so the finite, shrinking port area produces a
physical **discharge overpressure**. See PHYSICS.md Section 3.8.

**Level: quasi-static check valve.** The reed opens at ``p_dis*(1+rise)`` (the
existing 5% valve-rise as the spring preload) and flows out through the geometric
port area ``ports.port_open_area_m2(theta)`` (full ~18.4 mm^2 through most of
delivery, tapering to zero over the [339.6, 346.8] deg discharge window). No reed
inertia is modelled. Because the chamber pressure stays at or above ``p_dis``
throughout the open window, there is **no back-flow** in this quasi-static limit
(back-flow needs the reed to stay open past pressure reversal — a valve-dynamics,
Section 9, effect). The deliverable here is the **discharge overpressure loss**;
after the geometric port closes (``theta > discharge_close``) the port area is
zero, so the valve cannot vent the recompression pocket — that spike is limited
only by leakage (Section 3.7), as before.

Mass-scaled polytropic pressure and the compressible orifice are reused from
Section 3.7 (``mochi.leakage``); the tangency-gap leakage to suction runs
alongside the valve.
"""

from dataclasses import dataclass
from math import isfinite, pi

from mochi.chambers import (
    DEFAULT_POLYTROPIC_EXPONENT,
    DISCHARGE_PORT_PRESSURE_PA,
    DISCHARGE_VALVE_RISE_FRACTION,
    SUCTION_PORT_PRESSURE_PA,
    CycleTrace,
    build_cycle_trace,
    seal_over_half_angle_rad,
)
from mochi.indicated_work import indicated_work_j
from mochi.kinematics import RotaryGeometry
from mochi.leakage import (
    DISCHARGE_COEFFICIENT,
    ISENTROPIC_EXPONENT,
    LEAKAGE_GAP_M,
    SUCTION_DENSITY_KG_M3,
    orifice_mass_flow,
)
from mochi.ports import characteristic_angles, port_open_area_m2

FULL_TURN_RAD = 2.0 * pi


@dataclass(frozen=True, slots=True)
class ValvedCycle:
    """Reed-valve compression-chamber pressure history and losses (SI units)."""

    angles_rad: tuple[float, ...]
    compression_pressure_pa: tuple[float, ...]
    mass_kg: tuple[float, ...]
    delivery_peak_pa: float
    recompression_peak_pa: float
    delivered_mass_kg: float
    leaked_mass_kg: float
    volumetric_efficiency: float
    overpressure_work_j: float
    overpressure_power_w: float
    baseline_indicated_work_j: float
    valve_indicated_work_j: float
    valve_indicated_power_w: float
    frequency_hz: float
    samples: int

    def compression_at(self, crank_angle_rad: float) -> float:
        """Interpolate the valved compression pressure over ``[beta, window_entry]``."""

        angles = self.angles_rad
        if crank_angle_rad <= angles[0]:
            return self.compression_pressure_pa[0]
        if crank_angle_rad >= angles[-1]:
            return self.compression_pressure_pa[-1]
        position = (crank_angle_rad - angles[0]) / (angles[-1] - angles[0]) * (len(angles) - 1)
        lower = int(position)
        fraction = position - lower
        pressures = self.compression_pressure_pa
        return pressures[lower] + (pressures[lower + 1] - pressures[lower]) * fraction


def valved_cycle(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
    gap_m: float = LEAKAGE_GAP_M,
    discharge_coefficient: float = DISCHARGE_COEFFICIENT,
    valve_area_scale: float = 1.0,
    suction_density_kg_m3: float = SUCTION_DENSITY_KG_M3,
    isentropic_exponent: float = ISENTROPIC_EXPONENT,
) -> ValvedCycle:
    """Integrate the compression-chamber mass balance with a discharge reed valve.

    From the seal angle ``beta`` the chamber is marched in ``samples`` crank-angle
    steps with the pressure floating from a mass-scaled polytropic reference at
    the seal state. Once the pressure reaches the reed-opening pressure the valve
    vents to the discharge line through the geometric port area
    (:func:`mochi.ports.port_open_area_m2`); the tangency-gap leak to suction runs
    throughout. As the port area shrinks to zero over the discharge window the
    valve can no longer keep pace and the pressure rises above ``p_dis`` — the
    discharge overpressure. Returns the pressure/mass history, the delivery and
    recompression peaks, the delivered mass, the volumetric efficiency, and the
    valve overpressure loss (extra indicator-diagram work above ``p_dis``).
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("A valved cycle needs at least eight crank-angle samples.")
    if not isfinite(gap_m) or gap_m <= 0.0:
        raise ValueError("Leakage gap must be a positive, finite length in metres.")
    if not isfinite(suction_density_kg_m3) or suction_density_kg_m3 <= 0.0:
        raise ValueError("Suction density must be a positive, finite value in kg/m^3.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    angles = characteristic_angles(geometry)
    beta = angles.compression_start_rad
    discharge_close = angles.discharge_close_rad
    window_entry = FULL_TURN_RAD - seal_over_half_angle_rad(geometry)

    suction_pa = SUCTION_PORT_PRESSURE_PA
    discharge_pa = DISCHARGE_PORT_PRESSURE_PA
    opening_pa = discharge_pa * (1.0 + DISCHARGE_VALVE_RISE_FRACTION)
    exponent = DEFAULT_POLYTROPIC_EXPONENT
    omega = geometry.angular_speed_rad_s
    leak_area = gap_m * geometry.cylinder_height_m

    seal_volume = trace.compression_at(beta)
    seal_mass = suction_density_kg_m3 * seal_volume
    clearance_volume = trace.compression_at(discharge_close)

    def density(pressure_pa: float) -> float:
        return suction_density_kg_m3 * (pressure_pa / suction_pa) ** (1.0 / exponent)

    def flow(up_pa: float, down_pa: float, area_m2: float) -> float:
        if area_m2 <= 0.0 or up_pa <= down_pa:
            return 0.0
        return orifice_mass_flow(
            up_pa,
            down_pa,
            area_m2,
            discharge_coefficient=discharge_coefficient,
            isentropic_exponent=isentropic_exponent,
            upstream_density_kg_m3=density(up_pa),
        )

    span = window_entry - beta
    step = span / samples
    angles_out: list[float] = []
    pressures: list[float] = []
    masses: list[float] = []

    mass = seal_mass
    leaked = 0.0
    delivered = 0.0
    overpressure_work = 0.0
    valve_open = False
    delivery_peak = discharge_pa
    recompression_peak = discharge_pa
    previous_volume = seal_volume

    for index in range(samples + 1):
        theta = beta + step * index
        volume = trace.compression_at(theta)
        pressure = suction_pa * (mass / seal_mass) * (seal_volume / volume) ** exponent

        if not valve_open and pressure >= opening_pa:
            valve_open = True

        valve_out = 0.0
        if valve_open and theta < discharge_close:
            valve_area = valve_area_scale * port_open_area_m2(geometry, theta)
            valve_out = flow(pressure, discharge_pa, valve_area)
            if pressure > delivery_peak:
                delivery_peak = pressure
            overpressure_work += max(pressure - discharge_pa, 0.0) * max(
                previous_volume - volume, 0.0
            )
        elif theta >= discharge_close and pressure > recompression_peak:
            recompression_peak = pressure

        leak_out = flow(pressure, suction_pa, leak_area)
        mass -= (valve_out + leak_out) * step / omega
        delivered += valve_out * step / omega
        leaked += leak_out * step / omega

        angles_out.append(theta)
        pressures.append(pressure)
        masses.append(mass)
        previous_volume = volume

    displacement_volume = seal_volume - clearance_volume
    ideal_mass = suction_density_kg_m3 * displacement_volume
    volumetric_efficiency = max(delivered, 0.0) / ideal_mass

    # Valve-aware indicated work: the no-valve baseline (Section 3.5, discharge
    # clamped at p_dis) plus the discharge overpressure area the reed adds.
    frequency = geometry.frequency_hz
    baseline_work = indicated_work_j(geometry, trace=trace).net_work_j
    valve_work = baseline_work + overpressure_work

    return ValvedCycle(
        angles_rad=tuple(angles_out),
        compression_pressure_pa=tuple(pressures),
        mass_kg=tuple(masses),
        delivery_peak_pa=delivery_peak,
        recompression_peak_pa=recompression_peak,
        delivered_mass_kg=delivered,
        leaked_mass_kg=leaked,
        volumetric_efficiency=volumetric_efficiency,
        overpressure_work_j=overpressure_work,
        overpressure_power_w=overpressure_work * frequency,
        baseline_indicated_work_j=baseline_work,
        valve_indicated_work_j=valve_work,
        valve_indicated_power_w=valve_work * frequency,
        frequency_hz=frequency,
        samples=samples,
    )
