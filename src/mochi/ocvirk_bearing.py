"""Short-bearing (Ocvirk) crank-pin journal force -- eccentricity-coupled film.

Section 4.7 models the crank-pin journal friction with **Petroff's concentric
law**, which is load-*independent* and, at the peak reaction (Sommerfeld
``S ~ 0.07``, a heavily loaded bearing running well off centre), *underestimates*
the friction. This rung lifts the concentric assumption to the **short-bearing
(Ocvirk) hydrodynamic film**: a closed-form oil-film force that depends on the
journal eccentricity ratio ``epsilon = e / c_j`` (0 = centred, 1 = metal contact),
its rate ``epsilon_dot``, and the film **entrainment speed** -- with no lubrication
PDE solved (that is the later L2 Reynolds rung).

This module provides only a **quasi-static force law**: ``eccentricity_cycle`` finds
``epsilon(theta)`` from a *steady force balance* at each crank angle
(``epsilon_dot = 0``) -- it writes **no equation of motion and integrates nothing in
time**. The ``epsilon_dot`` squeeze term is exposed for the future D2 rotordynamics
EOM that would integrate ``F(epsilon, epsilon_dot, attitude)``, but that dynamic rung
is not built here. See PHYSICS.md Section 4.9.

**Derivation (short bearing / Ocvirk, pi-film).** With ``beta`` measured
circumferentially from the point of maximum film thickness, the film is
``h = c(1 + epsilon cos beta)``. The short-bearing approximation keeps only the
**axial** pressure-flow term of the Reynolds equation (valid for ``L/D <~ 0.5``;
here ``L_j/2 r_j = 0.74`` -- an order-of-magnitude, not asymptotic, application),
so the pressure integrates in closed form,

    p(beta, y) = (6 mu / c^2) (y^2 - L^2/4) (epsilon_dot cos beta - epsilon Omega
                 sin beta) / (1 + epsilon cos beta)^3,

with ``Omega = (omega_j + omega_b)/2 - psi_dot`` the entrainment (mean surface)
speed less the attitude whirl. Integrating the pressure over the journal with the
**pi-film (Gumbel/half-Sommerfeld)** cavitation boundary gives the closed-form
force, resolved into a component along the line of centres (``F_e``, the
load-supporting direction toward maximum film) and perpendicular to it (``F_t``):

    F_e = (mu r L^3 / c^2) [ epsilon_dot * P - epsilon Omega * Q ],
    F_t = (mu r L^3 / c^2) [ epsilon_dot * Q - epsilon Omega * S ],

    P = pi (1 + 2 eps^2) / (2 (1 - eps^2)^{5/2}),
    Q = -2 eps / (1 - eps^2)^2,
    S = pi / (2 (1 - eps^2)^{3/2}).

For **pure rotation** (``epsilon_dot = 0``) this reduces to the textbook Ocvirk
load capacity ``|W| = (mu Omega r L^3 / c^2) * eps sqrt(16 eps^2 + pi^2(1-eps^2)) /
(2(1-eps^2)^2)`` (with ``Omega = omega/2`` for a fixed sleeve this is the familiar
``mu omega r L^3 / (4 c^2) * ...``) and the attitude angle
``tan(psi) = pi sqrt(1-eps^2) / (4 eps)``.

**Crank-pin kinematics.** The journal (crank pin) turns at the shaft speed
``omega``; the sleeve (rotor bore) turns at the rotor spin ``omega dphi/dtheta``.
The **entrainment** speed is their mean ``Omega = omega - |omega_rel|/2`` (using
``omega_rel = omega(1 - dphi/dtheta)`` of Section 4.7, so ``Omega = omega/2`` when
the rotor does not spin), while the **shear/friction** is driven by their
*difference* ``omega_rel`` -- the same speed the Petroff model uses.

**Eccentric friction.** The Couette shear over the eccentric film integrates to
``T_f = 2 pi mu |omega_rel| r^3 L / (c sqrt(1 - eps^2))`` = Petroff's torque
divided by ``sqrt(1 - eps^2)`` -- so as ``eps -> 0`` it **recovers Petroff exactly**
(the reduction check of the validation ladder), and at the peak load it is larger,
capturing what the concentric model misses. The pressure-flow shear term is
second order and omitted at this (L1) fidelity.

Geometry and lubricant are shared with Section 4.7 (``r_j = 14.2`` mm CAD,
``L_j = 21`` mm, assumed ``c_j = 15`` micron, POE VG68 ``mu = 0.010`` Pa*s).
"""

from dataclasses import dataclass
from math import atan2, hypot, isfinite, pi, sqrt

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S
from mochi.chambers import CycleTrace, build_cycle_trace, seal_over_half_angle_rad
from mochi.journal_bearing import (
    JOURNAL_CLEARANCE_M,
    JOURNAL_LENGTH_M,
    JOURNAL_RADIUS_M,
    journal_relative_speed_rad_s,
)
from mochi.kinematics import RotaryGeometry
from mochi.true_gas_force import FULL_TURN_RAD, true_gas_load

# Largest eccentricity ratio the equilibrium solver will report; ``eps -> 1`` is
# metal-to-metal contact (zero film), so a load the film cannot support at any
# eps < 1 is clamped here and flagged by a near-zero minimum film thickness.
_MAX_ECCENTRICITY_RATIO = 1.0 - 1.0e-9


@dataclass(frozen=True, slots=True)
class OcvirkForce:
    """Short-bearing oil-film force at one journal state (SI units).

    ``radial_force_n`` (``F_e``) is along the line of centres toward *maximum* film
    thickness (the load-supporting direction); ``tangential_force_n`` (``F_t``) is
    perpendicular. ``attitude_angle_rad`` is the angle from the load line to the
    line of centres.
    """

    eccentricity_ratio: float
    eccentricity_rate_per_s: float
    entrainment_speed_rad_s: float
    radial_force_n: float
    tangential_force_n: float
    magnitude_n: float
    attitude_angle_rad: float
    min_film_thickness_m: float


def _check_dimensions(radius_m: float, length_m: float, clearance_m: float) -> None:
    for name, value in (("radius", radius_m), ("length", length_m), ("clearance", clearance_m)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Journal {name} must be a positive, finite length in metres.")


def _film_integrals(eccentricity_ratio: float) -> tuple[float, float, float]:
    """The three pi-film Sommerfeld integrals ``P, Q, S`` at this eccentricity."""

    one_minus = 1.0 - eccentricity_ratio * eccentricity_ratio
    p_int = pi * (1.0 + 2.0 * eccentricity_ratio**2) / (2.0 * one_minus**2.5)
    q_int = -2.0 * eccentricity_ratio / one_minus**2
    s_int = pi / (2.0 * one_minus**1.5)
    return p_int, q_int, s_int


def crank_pin_entrainment_speed_rad_s(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Film entrainment (mean surface) speed at the crank-pin journal.

    ``Omega = (omega_journal + omega_sleeve)/2 = omega - |omega_rel|/2`` with
    ``omega_rel`` the Section 4.7 shear speed. The crank pin turns at ``omega`` and
    the rotor bore at ``omega dphi/dtheta``; ``Omega`` reduces to ``omega/2`` when
    the rotor does not spin, matching the classic fixed-sleeve half-speed.
    """

    shaft_speed = geometry.angular_speed_rad_s
    relative_speed = journal_relative_speed_rad_s(geometry, crank_angle_rad)
    return shaft_speed - 0.5 * abs(relative_speed)


def short_bearing_force(
    eccentricity_ratio: float,
    entrainment_speed_rad_s: float,
    *,
    eccentricity_rate_per_s: float = 0.0,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
) -> OcvirkForce:
    """Closed-form Ocvirk short-bearing film force at one journal state.

    ``eccentricity_ratio`` is ``e/c`` in ``[0, 1)``; ``entrainment_speed_rad_s`` is
    ``Omega`` (see :func:`crank_pin_entrainment_speed_rad_s`);
    ``eccentricity_rate_per_s`` is ``d eps/dt`` (the squeeze term, 0 for a steady
    eccentricity). Returns the radial/tangential components, magnitude, attitude
    angle, and minimum film thickness.
    """

    if not (isfinite(eccentricity_ratio) and 0.0 <= eccentricity_ratio < 1.0):
        raise ValueError("Eccentricity ratio must be finite and in [0, 1).")
    if not isfinite(entrainment_speed_rad_s):
        raise ValueError("Entrainment speed must be finite.")
    if not isfinite(eccentricity_rate_per_s):
        raise ValueError("Eccentricity rate must be finite.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    _check_dimensions(radius_m, length_m, clearance_m)

    prefactor = viscosity_pa_s * radius_m * length_m**3 / clearance_m**2
    p_int, q_int, s_int = _film_integrals(eccentricity_ratio)
    radial = prefactor * (
        eccentricity_rate_per_s * p_int - eccentricity_ratio * entrainment_speed_rad_s * q_int
    )
    tangential = prefactor * (
        eccentricity_rate_per_s * q_int - eccentricity_ratio * entrainment_speed_rad_s * s_int
    )
    return OcvirkForce(
        eccentricity_ratio=eccentricity_ratio,
        eccentricity_rate_per_s=eccentricity_rate_per_s,
        entrainment_speed_rad_s=entrainment_speed_rad_s,
        radial_force_n=radial,
        tangential_force_n=tangential,
        magnitude_n=hypot(radial, tangential),
        attitude_angle_rad=atan2(tangential, radial),
        min_film_thickness_m=clearance_m * (1.0 - eccentricity_ratio),
    )


def static_load_capacity_n(
    eccentricity_ratio: float,
    entrainment_speed_rad_s: float,
    *,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
) -> float:
    """Steady (``eps_dot = 0``) load the short-bearing film supports at this ``eps``.

    Monotonic in ``eccentricity_ratio`` (0 at centre, diverging as ``eps -> 1``),
    so it inverts uniquely for the running eccentricity of a given load.
    """

    return short_bearing_force(
        eccentricity_ratio,
        entrainment_speed_rad_s,
        viscosity_pa_s=viscosity_pa_s,
        radius_m=radius_m,
        length_m=length_m,
        clearance_m=clearance_m,
    ).magnitude_n


def equilibrium_eccentricity_ratio(
    load_n: float,
    entrainment_speed_rad_s: float,
    *,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
    tolerance: float = 1.0e-10,
    iterations: int = 100,
) -> float:
    """Eccentricity ratio at which the film supports ``load_n`` under pure rotation.

    Inverts :func:`static_load_capacity_n` by bisection (it is monotonic in ``eps``).
    A load larger than the film can carry for any ``eps < 1`` is clamped to
    :data:`_MAX_ECCENTRICITY_RATIO` (metal contact); zero/negative load gives 0.
    """

    if not (isfinite(load_n) and load_n >= 0.0):
        raise ValueError("Journal load must be a non-negative, finite force in newtons.")
    if not isfinite(entrainment_speed_rad_s):
        raise ValueError("Entrainment speed must be finite.")
    _check_dimensions(radius_m, length_m, clearance_m)
    if load_n == 0.0 or entrainment_speed_rad_s == 0.0:
        return 0.0

    def capacity(eps: float) -> float:
        return static_load_capacity_n(
            eps,
            entrainment_speed_rad_s,
            viscosity_pa_s=viscosity_pa_s,
            radius_m=radius_m,
            length_m=length_m,
            clearance_m=clearance_m,
        )

    if capacity(_MAX_ECCENTRICITY_RATIO) < load_n:
        return _MAX_ECCENTRICITY_RATIO

    low, high = 0.0, _MAX_ECCENTRICITY_RATIO
    for _ in range(iterations):
        if high - low < tolerance:
            break
        mid = 0.5 * (low + high)
        if capacity(mid) < load_n:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def eccentric_friction_torque_nm(
    eccentricity_ratio: float,
    relative_speed_rad_s: float,
    *,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
) -> float:
    """Couette friction torque over the eccentric film, ``T_petroff / sqrt(1-eps^2)``.

    Reduces to the Section 4.7 Petroff torque as ``eps -> 0`` and grows with
    eccentricity. The pressure-flow shear contribution is second order (omitted).
    """

    if not (isfinite(eccentricity_ratio) and 0.0 <= eccentricity_ratio < 1.0):
        raise ValueError("Eccentricity ratio must be finite and in [0, 1).")
    if not isfinite(relative_speed_rad_s):
        raise ValueError("Relative speed must be finite.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    _check_dimensions(radius_m, length_m, clearance_m)

    petroff = (
        2.0 * pi * viscosity_pa_s * abs(relative_speed_rad_s) * radius_m**3 * length_m / clearance_m
    )
    return petroff / sqrt(1.0 - eccentricity_ratio * eccentricity_ratio)


@dataclass(frozen=True, slots=True)
class EccentricityCycle:
    """Journal eccentricity, film thickness, and friction over one revolution (SI)."""

    crank_angle_rad: tuple[float, ...]
    load_n: tuple[float, ...]
    entrainment_speed_rad_s: tuple[float, ...]
    relative_speed_rad_s: tuple[float, ...]
    eccentricity_ratio: tuple[float, ...]
    min_film_thickness_m: tuple[float, ...]
    friction_torque_nm: tuple[float, ...]
    friction_power_w: tuple[float, ...]
    peak_eccentricity_ratio: float
    minimum_film_thickness_m: float
    mean_friction_power_w: float
    petroff_mean_friction_power_w: float
    frequency_hz: float
    samples: int


def eccentricity_cycle(
    geometry: RotaryGeometry,
    *,
    samples: int = 360,
    trace: CycleTrace | None = None,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
) -> EccentricityCycle:
    """Trace the running eccentricity, minimum film, and eccentric friction.

    At each crank angle in the separated range ``[alpha, 2*pi - alpha]`` (outside
    seal-over, as in :func:`mochi.true_gas_force.peak_rotor_force_n`), the film is
    balanced against the true gas reaction ``|R_j(theta)|`` (Section 4.6) at the
    crank-pin entrainment speed to find the eccentricity ratio ``eps(theta)``, the
    minimum film ``c(1-eps)``, and the eccentric Couette friction. The mean
    friction power is reported alongside the Section 4.7 Petroff mean for the
    concentric-vs-eccentric comparison.
    """

    geometry.validate()
    if samples < 8:
        raise ValueError("The eccentricity sweep needs at least eight samples.")
    _check_dimensions(radius_m, length_m, clearance_m)

    trace = build_cycle_trace(geometry) if trace is None else trace
    half_angle = seal_over_half_angle_rad(geometry)
    span = FULL_TURN_RAD - 2.0 * half_angle

    angles: list[float] = []
    loads: list[float] = []
    entrainments: list[float] = []
    relatives: list[float] = []
    eccentricities: list[float] = []
    films: list[float] = []
    torques: list[float] = []
    powers: list[float] = []
    petroff_powers: list[float] = []
    for index in range(samples + 1):
        angle = half_angle + span * index / samples
        load = true_gas_load(geometry, angle, trace=trace).rotor_force_mag_n
        entrainment = crank_pin_entrainment_speed_rad_s(geometry, angle)
        relative = journal_relative_speed_rad_s(geometry, angle)
        eccentricity = equilibrium_eccentricity_ratio(
            load,
            entrainment,
            viscosity_pa_s=viscosity_pa_s,
            radius_m=radius_m,
            length_m=length_m,
            clearance_m=clearance_m,
        )
        torque = eccentric_friction_torque_nm(
            eccentricity,
            relative,
            viscosity_pa_s=viscosity_pa_s,
            radius_m=radius_m,
            length_m=length_m,
            clearance_m=clearance_m,
        )
        petroff_torque = eccentric_friction_torque_nm(
            0.0,
            relative,
            viscosity_pa_s=viscosity_pa_s,
            radius_m=radius_m,
            length_m=length_m,
            clearance_m=clearance_m,
        )
        angles.append(angle)
        loads.append(load)
        entrainments.append(entrainment)
        relatives.append(relative)
        eccentricities.append(eccentricity)
        films.append(clearance_m * (1.0 - eccentricity))
        torques.append(torque)
        powers.append(torque * abs(relative))
        petroff_powers.append(petroff_torque * abs(relative))

    mean_power = _trapezoid_mean(angles, powers)
    petroff_mean_power = _trapezoid_mean(angles, petroff_powers)
    return EccentricityCycle(
        crank_angle_rad=tuple(angles),
        load_n=tuple(loads),
        entrainment_speed_rad_s=tuple(entrainments),
        relative_speed_rad_s=tuple(relatives),
        eccentricity_ratio=tuple(eccentricities),
        min_film_thickness_m=tuple(films),
        friction_torque_nm=tuple(torques),
        friction_power_w=tuple(powers),
        peak_eccentricity_ratio=max(eccentricities),
        minimum_film_thickness_m=min(films),
        mean_friction_power_w=mean_power,
        petroff_mean_friction_power_w=petroff_mean_power,
        frequency_hz=geometry.frequency_hz,
        samples=samples,
    )


def _trapezoid_mean(angles: list[float], values: list[float]) -> float:
    """Trapezoidal average of ``values`` over the sampled ``angles`` span."""

    total = 0.0
    for index in range(len(angles) - 1):
        step = angles[index + 1] - angles[index]
        total += 0.5 * (values[index] + values[index + 1]) * step
    return total / (angles[-1] - angles[0])
