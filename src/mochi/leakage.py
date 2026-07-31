"""Gap leakage in a chamber mass balance, and volumetric efficiency.

The no-leakage port-timed rule (Section 3.4) over-compresses the trapped
clearance gas during recompression to a ~10 MPa spike, a strict upper bound
because nothing bleeds the over-pressurised gas out. This module adds a **single
equivalent leakage orifice** from the compression/recompression chamber to the
suction side and integrates the chamber **mass balance** over one revolution. The
leak caps the recompression spike and yields the **volumetric efficiency**. See
PHYSICS.md Section 3.7.

It is a *parallel* model: the accepted `mochi.chambers.port_timed_pressures`
baseline is untouched, as with the two side-by-side pressure rules of Sections
3.2 and 3.4.

**Mass-scaled polytropic pressure (no temperature).** The rest of `mochi` is pure
polytropic `p Vⁿ = const` with no absolute temperature; the mass balance keeps
that form. For a chamber sealed at ``(p_ref, V_ref, m_ref)`` the pressure with a
changed mass is ``p = p_ref · (m/m_ref) · (V_ref/V)ⁿ`` (the ideal-gas identity
``p = m R T / V`` with ``T`` on the fixed-mass polytropic ``T V^{n-1}`` path; the
``m = m_ref`` case is the existing `_polytropic_pa`). Density along the process is
``ρ = ρ_suc (p/p_suc)^{1/n}``.

**Introduced R410A properties** (hard-coded from CoolProp exactly as ``n = 1.07``
was, not a runtime dependency): the suction saturated-vapour density and the
isentropic (nozzle) exponent. **Introduced leakage geometry:** the rotor-cylinder
tangency is an ideal zero-gap seal in the model, so an effective leakage gap is
assumed (editable constant), all paths lumped into one equivalent orifice.
"""

from dataclasses import dataclass
from math import isfinite, pi, sqrt

from mochi.chambers import (
    DEFAULT_POLYTROPIC_EXPONENT,
    DISCHARGE_PORT_PRESSURE_PA,
    DISCHARGE_VALVE_RISE_FRACTION,
    SUCTION_PORT_PRESSURE_PA,
    CycleTrace,
    build_cycle_trace,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.kinematics import RotaryGeometry
from mochi.ports import characteristic_angles

FULL_TURN_RAD = 2.0 * pi

# R410A saturated-vapour properties at the 0.82 MPa suction port (0.85 degC
# evaporating), from CoolProp HEOS. Hard-coded like DEFAULT_POLYTROPIC_EXPONENT,
# not a runtime dependency.
SUCTION_DENSITY_KG_M3 = 31.4
# Isentropic (nozzle) exponent k = a**2 rho / p = 1.10 at the suction state;
# distinct from the compression pv-index n = 1.07.
ISENTROPIC_EXPONENT = 1.10

# Assumed effective leakage gap (the rotor-cylinder tangency is an ideal
# zero-gap line seal in the model, and no end-face path is defined, so a single
# equivalent orifice is introduced; override as needed). The orifice area is this
# gap times the cylinder height (the tangency-line length).
LEAKAGE_GAP_M = 5.0e-6
DISCHARGE_COEFFICIENT = 0.6


@dataclass(frozen=True, slots=True)
class LeakyCycle:
    """Leakage-coupled compression-chamber pressure history and efficiency (SI)."""

    angles_rad: tuple[float, ...]
    compression_pressure_pa: tuple[float, ...]
    mass_kg: tuple[float, ...]
    no_leak_peak_pa: float
    capped_peak_pa: float
    delivered_mass_kg: float
    leaked_mass_kg: float
    volumetric_efficiency: float
    frequency_hz: float
    samples: int

    def compression_at(self, crank_angle_rad: float) -> float:
        """Interpolate the leaky compression pressure over ``[beta, window_entry]``."""

        angles = self.angles_rad
        normalized = crank_angle_rad
        if normalized <= angles[0]:
            return self.compression_pressure_pa[0]
        if normalized >= angles[-1]:
            return self.compression_pressure_pa[-1]
        span = angles[-1] - angles[0]
        position = (normalized - angles[0]) / span * (len(angles) - 1)
        lower = int(position)
        fraction = position - lower
        pressures = self.compression_pressure_pa
        return pressures[lower] + (pressures[lower + 1] - pressures[lower]) * fraction


def orifice_mass_flow(
    upstream_pressure_pa: float,
    downstream_pressure_pa: float,
    area_m2: float,
    *,
    discharge_coefficient: float = DISCHARGE_COEFFICIENT,
    isentropic_exponent: float = ISENTROPIC_EXPONENT,
    upstream_density_kg_m3: float,
) -> float:
    """Compressible isentropic-nozzle mass flow through a gap (kg/s), choked-capable.

    Subsonic: ``m = Cd A sqrt(2k/(k-1) · p_up ρ_up · (r^{2/k} - r^{(k+1)/k}))`` with
    ``r = p_down/p_up``; choked (``r <= (2/(k+1))^{k/(k-1)}``):
    ``m = Cd A sqrt(k p_up ρ_up) (2/(k+1))^{(k+1)/(2(k-1))}``. Returns 0 for a
    non-positive pressure difference (no back-flow modelled here).
    """

    for name, value in (
        ("upstream pressure", upstream_pressure_pa),
        ("upstream density", upstream_density_kg_m3),
        ("area", area_m2),
    ):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Orifice {name} must be a positive, finite value.")
    if not isfinite(isentropic_exponent) or isentropic_exponent <= 1.0:
        raise ValueError("Isentropic exponent must be a finite value greater than 1.")
    if upstream_pressure_pa <= downstream_pressure_pa:
        return 0.0

    k = isentropic_exponent
    ratio = max(downstream_pressure_pa, 0.0) / upstream_pressure_pa
    critical = (2.0 / (k + 1.0)) ** (k / (k - 1.0))
    coefficient = discharge_coefficient * area_m2
    product = upstream_pressure_pa * upstream_density_kg_m3
    if ratio <= critical:
        flux = sqrt(k * product) * (2.0 / (k + 1.0)) ** ((k + 1.0) / (2.0 * (k - 1.0)))
    else:
        expansion = ratio ** (2.0 / k) - ratio ** ((k + 1.0) / k)
        flux = sqrt(2.0 * k / (k - 1.0) * product * max(expansion, 0.0))
    return coefficient * flux


def leaky_cycle(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
    gap_m: float = LEAKAGE_GAP_M,
    discharge_coefficient: float = DISCHARGE_COEFFICIENT,
    suction_density_kg_m3: float = SUCTION_DENSITY_KG_M3,
    isentropic_exponent: float = ISENTROPIC_EXPONENT,
) -> LeakyCycle:
    """Integrate the compression chamber's leakage-coupled mass balance over a cycle.

    From the seal angle ``beta`` to the seal-over entry the chamber is marched in
    ``samples`` crank-angle steps: compression (closed, leaking to suction) until
    the pressure reaches the valve-opening pressure, delivery at the discharge
    pressure, then recompression (closed, leaking) of the trapped clearance gas.
    The leak bleeds mass out so the recompression pressure is capped below the
    no-leakage ~10 MPa bound. Returns the pressure/mass history, the capped and
    no-leak peaks, the delivered mass, and the volumetric efficiency.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("A leaky cycle needs at least eight crank-angle samples.")
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
    area = gap_m * geometry.cylinder_height_m
    discharge_density = suction_density_kg_m3 * (discharge_pa / suction_pa) ** (1.0 / exponent)

    seal_volume = trace.compression_at(beta)
    seal_mass = suction_density_kg_m3 * seal_volume
    clearance_volume = trace.compression_at(discharge_close)
    clearance_mass = discharge_density * clearance_volume

    def leak(pressure_pa: float) -> float:
        density = suction_density_kg_m3 * (pressure_pa / suction_pa) ** (1.0 / exponent)
        return orifice_mass_flow(
            pressure_pa,
            suction_pa,
            area,
            discharge_coefficient=discharge_coefficient,
            isentropic_exponent=isentropic_exponent,
            upstream_density_kg_m3=density,
        )

    span = window_entry - beta
    step = span / samples
    angles_out: list[float] = []
    pressures: list[float] = []
    masses: list[float] = []

    mass = seal_mass
    leaked = 0.0
    delivered = 0.0
    valve_open = False
    recompressing = False
    reference_mass = seal_mass
    reference_volume = seal_volume
    reference_pressure = suction_pa
    previous_volume = seal_volume

    for index in range(samples + 1):
        theta = beta + step * index
        volume = trace.compression_at(theta)

        if theta >= discharge_close:
            # Recompression: closed, trapped clearance gas, leaking out.
            if not recompressing:
                recompressing = True
                mass = clearance_mass
                reference_mass = clearance_mass
                reference_volume = clearance_volume
                reference_pressure = discharge_pa
            pressure = (
                reference_pressure
                * (mass / reference_mass)
                * (reference_volume / volume) ** exponent
            )
            outflow = leak(pressure)
            mass -= outflow * step / omega
            leaked += outflow * step / omega
        elif not valve_open:
            # Compression: closed, leaking, until the valve-opening pressure.
            pressure = (
                reference_pressure
                * (mass / reference_mass)
                * (reference_volume / volume) ** exponent
            )
            if pressure >= opening_pa:
                valve_open = True
                pressure = discharge_pa
                mass = discharge_density * volume
            else:
                outflow = leak(pressure)
                mass -= outflow * step / omega
                leaked += outflow * step / omega
        if valve_open and not recompressing:
            # Delivery at the discharge pressure; volume shrink pushes mass out.
            pressure = discharge_pa
            delivered += discharge_density * max(previous_volume - volume, 0.0)
            outflow = leak(pressure)
            mass = discharge_density * volume
            leaked += outflow * step / omega
            delivered -= outflow * step / omega

        angles_out.append(theta)
        pressures.append(pressure)
        masses.append(mass)
        previous_volume = volume

    recompression_pressures = [
        pressures[i] for i, theta in enumerate(angles_out) if theta >= discharge_close
    ]
    capped_peak = max(recompression_pressures) if recompression_pressures else discharge_pa
    no_leak_peak = port_timed_pressures(geometry, window_entry, trace=trace).residual_pressure_pa

    displacement_volume = seal_volume - clearance_volume
    ideal_mass = suction_density_kg_m3 * displacement_volume
    volumetric_efficiency = max(delivered, 0.0) / ideal_mass

    return LeakyCycle(
        angles_rad=tuple(angles_out),
        compression_pressure_pa=tuple(pressures),
        mass_kg=tuple(masses),
        no_leak_peak_pa=no_leak_peak,
        capped_peak_pa=capped_peak,
        delivered_mass_kg=delivered,
        leaked_mass_kg=leaked,
        volumetric_efficiency=volumetric_efficiency,
        frequency_hz=geometry.frequency_hz,
        samples=samples,
    )
