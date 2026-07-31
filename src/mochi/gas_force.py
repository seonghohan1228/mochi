"""Gas pressure force on the rotor and the net gas torque about the crank axis.

The chamber gas (Section 3.2/3.4) presses on the rotor outside diameter. This
module turns that pressure boundary condition into the ``F_gas`` contribution of
PHYSICS.md Section 4.4: the net gas force on the rotor, the reaction on the fixed
vane, and the net gas torque about the cylinder/crank axis ``O`` -- the indicated
torque. See PHYSICS.md Section 4.5 for the derivation and validity.

Two facts from the circular-rotor approximation make the result a closed form:

* The rotor OD is a circle about the rotor centre ``O_r = e (sin theta,
  cos theta)``. A chamber at uniform pressure ``p`` acts radially on its arc, so
  the net force on an arc is fixed by its endpoints, and the whole-rotor
  resultant passes through ``O_r``.
* Because the resultant passes through ``O_r``, the torque about ``O`` is simply
  ``T_gas = O_r x F_rotor`` (z-component).

The vane is part of the ground body (its tip is fixed, Section 4.2), so gas load
on the vane reacts to the cylinder and carries no shaft torque; the vane force is
reported only to close the force balance.

Sign discipline (Section 4.2): moments are written about ``+z`` (counter-
clockwise positive). With ``T_gas = O_r x F_rotor`` the crank-angle integral
returns the indicated work, ``oint T_gas d(theta) = -oint p dV = W`` (Section
3.5), on the circular-rotor volume basis (``chamber_areas``). This is the
independent cross-check binding the force model to the P-V work.
"""

from dataclasses import dataclass
from math import asin, cos, isfinite, pi, sin, sqrt

from mochi.chambers import (
    CycleTrace,
    build_cycle_trace,
    chamber_areas,
    in_seal_over_window,
    port_timed_pressures,
    seal_over_half_angle_rad,
)
from mochi.kinematics import RotaryGeometry

FULL_TURN_RAD = 2.0 * pi


@dataclass(frozen=True, slots=True)
class GasLoad:
    """Gas pressure loads at one crank angle (global frame, SI units)."""

    crank_angle_rad: float
    rotor_force_n: tuple[float, float]
    rotor_torque_nm: float
    vane_force_n: tuple[float, float]
    suction_pressure_pa: float
    compression_pressure_pa: float


@dataclass(frozen=True, slots=True)
class GasTorqueWork:
    """Gas-torque cross-check against the indicated work (SI units)."""

    torque_work_j: float
    reference_pv_work_j: float
    power_w: float
    frequency_hz: float
    samples: int


def _foot_normal_angles(geometry: RotaryGeometry, crank_angle_rad: float) -> tuple[float, float]:
    """Rotor-OD normal angles of the two vane feet, measured from ``+y`` clockwise.

    The rotor-OD point with outward normal ``(sin phi, cos phi)`` sits at
    ``x = e sin(theta) + R_r sin(phi)``. Setting ``x = +/- w/2`` gives the inlet
    and outlet feet where the vane flanks meet the rotor; the near-top
    (principal) root is the physical contact.
    """

    rotor_x = geometry.eccentricity_m * sin(crank_angle_rad)
    half_width = 0.5 * geometry.vane_width_m
    rotor_radius = geometry.rotor_radius_m
    inlet_foot = asin((half_width - rotor_x) / rotor_radius)
    outlet_foot = asin((-half_width - rotor_x) / rotor_radius)
    return inlet_foot, outlet_foot


def _arc_force(
    geometry: RotaryGeometry,
    pressure_pa: float,
    phi_start: float,
    phi_end: float,
) -> tuple[float, float]:
    """Closed-form net gas force on a uniform-pressure rotor-OD arc.

    Integrates ``dF = -p n R_r H d(phi)`` with the outward normal
    ``n = (sin phi, cos phi)`` from ``phi_start`` to ``phi_end``. The result
    depends only on the endpoints:
    ``F = -p R_r H (cos phi_s - cos phi_e, sin phi_e - sin phi_s)``.
    """

    scale = pressure_pa * geometry.rotor_radius_m * geometry.cylinder_height_m
    force_x = -scale * (cos(phi_start) - cos(phi_end))
    force_y = -scale * (sin(phi_end) - sin(phi_start))
    return force_x, force_y


def gas_load(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    trace: CycleTrace | None = None,
) -> GasLoad:
    """Return the gas pressure force on the rotor and the net gas torque about ``O``.

    The suction chamber wets the rotor OD from the inlet vane foot clockwise to
    the rotor-cylinder contact (bore angle ``theta``); the compression chamber
    wets it from the contact clockwise to the outlet vane foot. The vane
    footprint over the top is excluded. Pressures are the port-timed rule of
    Section 3.4. The torque is ``O_r x F_rotor`` about ``+z`` (Section 4.2).

    Raises ``ValueError`` inside the seal-over window, where the two chambers
    merge and the contact hides under the vane, so no separated arcs exist.
    """

    geometry.validate()
    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    if in_seal_over_window(geometry, crank_angle_rad):
        raise ValueError(
            "Gas load is undefined in the seal-over window: the rotor-cylinder "
            "contact hides under the vane and the chambers merge."
        )

    trace = build_cycle_trace(geometry) if trace is None else trace
    pressures = port_timed_pressures(geometry, crank_angle_rad, trace=trace)
    suction_pa = pressures.suction_pressure_pa
    compression_pa = pressures.compression_pressure_pa

    normalized = crank_angle_rad % FULL_TURN_RAD
    inlet_foot, outlet_foot = _foot_normal_angles(geometry, normalized)
    contact = normalized  # the OD normal points along the contact at bore angle theta

    # Suction arc: inlet foot -> contact. Compression arc: contact -> outlet
    # foot, taken the long way round (through the bottom), so its end angle is
    # lifted by a full turn to keep phi increasing.
    suction_force = _arc_force(geometry, suction_pa, inlet_foot, contact)
    compression_force = _arc_force(geometry, compression_pa, contact, outlet_foot + FULL_TURN_RAD)
    rotor_force = (
        suction_force[0] + compression_force[0],
        suction_force[1] + compression_force[1],
    )

    rotor_x = geometry.eccentricity_m * sin(normalized)
    rotor_y = geometry.eccentricity_m * cos(normalized)
    rotor_torque = rotor_x * rotor_force[1] - rotor_y * rotor_force[0]

    vane_force = _vane_force(
        geometry, normalized, suction_pa, compression_pa, inlet_foot, outlet_foot
    )

    return GasLoad(
        crank_angle_rad=crank_angle_rad,
        rotor_force_n=rotor_force,
        rotor_torque_nm=rotor_torque,
        vane_force_n=vane_force,
        suction_pressure_pa=suction_pa,
        compression_pressure_pa=compression_pa,
    )


def _vane_force(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    suction_pa: float,
    compression_pa: float,
    inlet_foot: float,
    outlet_foot: float,
) -> tuple[float, float]:
    """Net gas force on the two exposed vane flanks (first-cut, balance only).

    Each vertical flank at ``x = +/- w/2`` is wetted from the rotor foot up to
    the bore. The suction chamber pushes the inlet flank toward ``-x`` and the
    compression chamber pushes the outlet flank toward ``+x``; only the
    horizontal component survives. The vane tip and root fillets are neglected.
    """

    half_width = 0.5 * geometry.vane_width_m
    rotor_y = geometry.eccentricity_m * cos(crank_angle_rad)
    rotor_radius = geometry.rotor_radius_m
    height = geometry.cylinder_height_m
    bore_y = sqrt(geometry.cylinder_radius_m**2 - half_width**2)
    inlet_flank = bore_y - (rotor_y + rotor_radius * cos(inlet_foot))
    outlet_flank = bore_y - (rotor_y + rotor_radius * cos(outlet_foot))
    force_x = (compression_pa * outlet_flank - suction_pa * inlet_flank) * height
    return (force_x, 0.0)


def gas_torque_work_j(
    geometry: RotaryGeometry,
    *,
    samples: int = 720,
    trace: CycleTrace | None = None,
) -> GasTorqueWork:
    """Cross-check the gas torque against the indicated work over one revolution.

    ``torque_work_j`` integrates the gas torque, ``oint T_gas d(theta)``, and
    ``reference_pv_work_j`` integrates ``-oint p dV`` on the same circular-rotor
    volumes (``chamber_areas`` times the height) with the same port-timed
    pressures. The identity ``oint T_gas d(theta) = -oint p dV = W`` (PHYSICS.md
    Section 4.5) makes the two equal; ``power_w`` multiplies the torque work by
    the shaft frequency for comparison with the Section 3.5 headline.

    Both integrals run over the separated range ``[alpha, 2*pi - alpha]`` (the
    seal-over half-angle ``alpha`` on each end), where two distinct chambers
    exist. The excluded merge window carries negligible work, exactly as in
    :func:`mochi.indicated_work.indicated_work_j`.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("The gas-torque cross-check needs at least eight samples.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    half_angle = seal_over_half_angle_rad(geometry)
    start = half_angle
    span = FULL_TURN_RAD - 2.0 * half_angle
    angles = tuple(start + span * index / samples for index in range(samples + 1))

    torques = tuple(gas_load(geometry, angle, trace=trace).rotor_torque_nm for angle in angles)
    height = geometry.cylinder_height_m
    suction_v = tuple(chamber_areas(geometry, angle).suction_area_m2 * height for angle in angles)
    compression_v = tuple(
        chamber_areas(geometry, angle).discharge_area_m2 * height for angle in angles
    )
    pressures = tuple(port_timed_pressures(geometry, angle, trace=trace) for angle in angles)
    suction_p = tuple(p.suction_pressure_pa for p in pressures)
    compression_p = tuple(p.compression_pressure_pa for p in pressures)

    torque_work = 0.0
    reference_work = 0.0
    for index in range(samples):
        torque_work += (
            0.5 * (torques[index] + torques[index + 1]) * (angles[index + 1] - angles[index])
        )
        reference_work -= (
            0.5
            * (suction_p[index] + suction_p[index + 1])
            * (suction_v[index + 1] - suction_v[index])
        )
        reference_work -= (
            0.5
            * (compression_p[index] + compression_p[index + 1])
            * (compression_v[index + 1] - compression_v[index])
        )

    return GasTorqueWork(
        torque_work_j=torque_work,
        reference_pv_work_j=reference_work,
        power_w=torque_work * geometry.frequency_hz,
        frequency_hz=geometry.frequency_hz,
        samples=samples,
    )
