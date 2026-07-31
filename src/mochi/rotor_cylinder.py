"""Rotor-cylinder sealing-line contact: sliding kinematics and boundary friction.

The rotor outer diameter seals against the cylinder bore along a line near the vane.
Unlike a conventional **rolling-piston** compressor -- where the roller spins on the
shaft and *rolls* on the cylinder with little slip -- this is a **swing** compressor:
the rotor is tied to the fixed vane through the swing bush, so it **swings**
(+/-10.37 deg orientation, *net-zero rotation per revolution*) instead of spinning.
A body that does not rotate cannot roll, so the sealing line **slides**. The contact
tangential sliding speed at the point ``P = O_b + R_r n_hat`` (``n_hat = O_b/|O_b|``)
is, with the cylinder stationary,

    v_slide(theta) = O_b_dot . t_hat + phi_orient_dot * R_r
                   = -e omega + omega (dphi_orient/dtheta) R_r,

the first term the orbital sliding (``-e omega = -0.85`` m/s, constant) and the second
the swing modulation (``+/-1.15`` m/s). Over a cycle ``|v_slide|`` averages ~0.95 m/s
(peak ~2.0) -- i.e. the seal genuinely rubs, so this loss is potentially comparable to
the crank-pin journal loss, not a small rolling correction.

**Friction law and the contact force.** In **boundary** lubrication the loss is
``P = mu * N_c * |v_slide|`` integrated over the cycle; ``v_slide`` (kinematics) and
``mu`` (assumed) are direct, and the contact normal force ``N_c`` -- statically
indeterminate on its own (the crank-pin bearing, this contact, and the bush/vane
over-constrain the rotor) -- is resolved **self-consistently** by
:func:`integrate_sealing_contact`: a compliant Hertz line contact
(:func:`hertz_line_contact_force_n`) is added to the Section 4.13 rotor EOM, so the
load splits between the journal film and the contact by their true stiffnesses. That
also removes the ~6 um free-orbit bore-penetration artifact (it collapses to the
physical ~1.4 um Hertz deflection). The centrifugal ``m omega^2 e ~ 44 N`` (always
outward) keeps the seal largely engaged, giving ``N_c`` mean ~140 N and a
rotor-cylinder loss of ~5-7 W -- first-order, comparable to the crank-pin journal.
:func:`contact_normal_force_n` is a cheaper quasi-static (rigid-contact) estimate;
:func:`rotor_cylinder_friction_power_w` / :func:`mixed_ehl_friction_power_w` apply the
boundary / mixed closures to any ``N_c(theta)``.

**Friction coefficient.** Steel-on-steel, refrigerant-oil boundary lubrication:
``mu ~ 0.1`` (dry-ish boundary 0.1-0.15). Mixed-lubrication analyses of rolling-piston
roller-cylinder / vane-roller contacts report **0.04-0.08** for the asperity-borne
component (Greenwood-Tripp asperity load + Patir-Cheng average-flow film). Editable
via :data:`ROTOR_CYLINDER_FRICTION_COEFF`. See PHYSICS.md Section 4.15 and the
references there (Yanagisawa & Shimizu 1985; the 2024 mixed-lubrication review).
"""

from collections.abc import Callable
from dataclasses import dataclass
from math import cos, hypot, isfinite, pi, sin, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from mochi.chambers import CycleTrace, build_cycle_trace, seal_over_half_angle_rad
from mochi.journal_bearing import JOURNAL_CLEARANCE_M
from mochi.kinematics import RotaryGeometry, prescribed_state
from mochi.rotor_dynamics import ROTOR_MASS_KG, _film_force_on_rotor, _prescribe_gas_and_speed
from mochi.true_gas_force import FULL_TURN_RAD, true_gas_load

# Steel-on-steel boundary-lubrication friction coefficient (editable). Compressor
# mixed-lubrication studies of the roller-cylinder / vane-roller contact report
# 0.04-0.08 for the asperity-borne part; 0.1 is a conservative pure-boundary value.
ROTOR_CYLINDER_FRICTION_COEFF = 0.10

# Mixed-EHL defaults for the (nearly conformal) rotor-cylinder line contact.
# Film parameter Lambda = h_min / sigma uses the COMPOSITE RMS roughness (Rq), not Ra:
# from the Ra 0.3 um design spec on both surfaces, Rq = 1.25 Ra = 0.375 um each and the
# composite is sqrt(2) * 0.375 = 0.53 um. (Editable per call via ``roughness_m``.)
COMPOSITE_ROUGHNESS_M = 0.53e-6  # composite RMS, from Ra 0.3 um (design standard)
PRESSURE_VISCOSITY_COEFF_PA = 15.0e-9  # POE oil, alpha ~ 15 GPa^-1
REDUCED_MODULUS_PA = 210.0e9 / (1.0 - 0.3**2)  # steel-steel, E' = E/(1-nu^2)

# Palmgren steel line-contact elastic approach: delta[mm] = coeff * Q[N]^0.9 / L[mm]^0.8.
# This is the compliant contact law that closes the indeterminate N_c self-consistently.
PALMGREN_LINE_CONTACT_COEFF = 3.84e-5

# Central-difference half-step for the orientation swing rate (analytic pose).
_ORIENTATION_STEP_RAD = 1.0e-5


def contact_sliding_speed_m_s(geometry: RotaryGeometry, crank_angle_rad: float) -> float:
    """Signed tangential sliding speed at the rotor-cylinder sealing line (m/s).

    ``v_slide = -e omega + omega (dphi_orient/dtheta) R_r`` (module docstring). Negative
    is the dominant orbital-sliding sense; the swing term modulates it. The cylinder is
    stationary, so this is the rotor surface speed at the contact.
    """

    if not isfinite(crank_angle_rad):
        raise ValueError("Crank angle must be finite.")
    geometry.validate()
    omega = geometry.angular_speed_rad_s
    throw = geometry.eccentricity_m
    rotor_radius = geometry.rotor_radius_m
    step = _ORIENTATION_STEP_RAD

    def orientation(angle: float) -> float:
        return prescribed_state(geometry, angle).rotor_orientation_rad

    d_orientation = (orientation(crank_angle_rad + step) - orientation(crank_angle_rad - step)) / (
        2.0 * step
    )
    return -throw * omega + omega * d_orientation * rotor_radius


def cycle_mean_abs_sliding_speed_m_s(geometry: RotaryGeometry, *, samples: int = 360) -> float:
    """Cycle-mean ``|v_slide|`` at the sealing line (the boundary-friction speed)."""

    if samples < 8:
        raise ValueError("Need at least 8 samples.")
    angles = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    speeds = np.array([contact_sliding_speed_m_s(geometry, float(a)) for a in angles])
    return float(np.mean(np.abs(speeds)))


def rotor_cylinder_friction_power_w(
    geometry: RotaryGeometry,
    normal_force_n: float | Callable[[float], float],
    *,
    friction_coefficient: float = ROTOR_CYLINDER_FRICTION_COEFF,
    samples: int = 360,
) -> float:
    """Cycle-mean boundary friction power ``<mu * N_c * |v_slide|>`` (W).

    ``normal_force_n`` is the contact normal force ``N_c`` -- the input the model does
    **not** yet resolve (statically indeterminate; needs the contact/EHL rung). Pass a
    constant (N) or a callable ``theta -> N`` once an ``N_c(theta)`` estimate exists.
    ``friction_coefficient`` defaults to the editable boundary value.
    """

    if not (isfinite(friction_coefficient) and friction_coefficient >= 0.0):
        raise ValueError("Friction coefficient must be finite and non-negative.")
    if samples < 8:
        raise ValueError("Need at least 8 samples.")
    angles = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    total = 0.0
    for angle in angles:
        theta = float(angle)
        load = normal_force_n(theta) if callable(normal_force_n) else normal_force_n
        if not (isfinite(load) and load >= 0.0):
            raise ValueError("Contact normal force must be finite and non-negative.")
        total += friction_coefficient * load * abs(contact_sliding_speed_m_s(geometry, theta))
    return total / samples


def contact_normal_force_n(
    geometry: RotaryGeometry,
    crank_angle_rad: float,
    *,
    trace: CycleTrace | None = None,
    rotor_mass_kg: float = ROTOR_MASS_KG,
) -> float:
    """Quasi-static rotor-cylinder contact normal force ``N_c(theta)`` (N).

    ``N_c`` is statically indeterminate in general; this is the **rigid-contact /
    journal-tangential** estimate: with the bore pinned so the OD just touches the bore
    (no penetration), the crank-pin journal reacts the *tangential* gas load and the
    contact reacts the net *outward radial* load -- the gas radial component plus the
    orbital centrifugal ``m omega^2 e`` (44 N, always outward). Clamped at zero where
    the net radial pulls the rotor **off** the bore (the seal lifts, ``N_c = 0``). A
    time-accurate value would come from the coupled contact/EHL rung; a stiff dynamic
    penalty contact gives a similar mean and a higher transient peak.

    Returns 0 in the seal-over window (gas load undefined).
    """

    geometry.validate()
    half_angle = seal_over_half_angle_rad(geometry)
    wrapped = crank_angle_rad % FULL_TURN_RAD
    if wrapped < half_angle or wrapped > FULL_TURN_RAD - half_angle:
        return 0.0
    trace = build_cycle_trace(geometry) if trace is None else trace
    force_x, force_y = true_gas_load(geometry, wrapped, trace=trace).rotor_force_n
    omega = geometry.angular_speed_rad_s
    radial = (
        force_x * sin(wrapped)
        + force_y * cos(wrapped)
        + rotor_mass_kg * omega * omega * geometry.eccentricity_m
    )
    return max(radial, 0.0)


def mixed_ehl_friction_power_w(
    geometry: RotaryGeometry,
    normal_force_n: float | Callable[[float], float],
    *,
    roughness_m: float = COMPOSITE_ROUGHNESS_M,
    pressure_viscosity_pa: float = PRESSURE_VISCOSITY_COEFF_PA,
    reduced_modulus_pa: float = REDUCED_MODULUS_PA,
    viscosity_pa_s: float = 0.010,
    boundary_coefficient: float = ROTOR_CYLINDER_FRICTION_COEFF,
    samples: int = 360,
) -> float:
    """Cycle-mean mixed-lubrication friction power at the sealing line (W).

    The rotor-cylinder contact is nearly conformal (equivalent radius
    ``R_eq = R_r R_cyl / (R_cyl - R_r) ~ 291`` mm), so a partial EHL film forms. Per
    crank angle the Dowson-Higginson line-contact minimum film ``h_min`` gives the
    specific film thickness ``Lambda = h_min / sigma`` (entrainment ``ubar = |v|/2``,
    the bore stationary); the load and friction split between the film (viscous shear
    ``mu_0 v / h_min`` over the Hertz zone) and the asperities (Greenwood-Tripp-style
    load fraction ``f_a = clip((3 - Lambda)/2, 0, 1)``, carrying ``mu_boundary``). It
    reduces to pure boundary as ``Lambda -> 0`` and to viscous-only as ``Lambda >= 3``.
    ``normal_force_n`` is ``N_c`` (scalar or ``theta -> N``); see
    :func:`contact_normal_force_n`.
    """

    for name, value in (
        ("roughness", roughness_m),
        ("pressure-viscosity", pressure_viscosity_pa),
        ("reduced modulus", reduced_modulus_pa),
        ("viscosity", viscosity_pa_s),
    ):
        if not (isfinite(value) and value > 0.0):
            raise ValueError(f"The {name} parameter must be positive and finite.")
    if samples < 8:
        raise ValueError("Need at least 8 samples.")

    angles = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    total = 0.0
    for angle in angles:
        theta = float(angle)
        load = normal_force_n(theta) if callable(normal_force_n) else normal_force_n
        speed = abs(contact_sliding_speed_m_s(geometry, theta))
        total += (
            _mixed_ehl_friction_force_n(
                geometry,
                load,
                speed,
                roughness_m=roughness_m,
                pressure_viscosity_pa=pressure_viscosity_pa,
                reduced_modulus_pa=reduced_modulus_pa,
                viscosity_pa_s=viscosity_pa_s,
                boundary_coefficient=boundary_coefficient,
            )
            * speed
        )
    return total / samples


def ehl_film_thickness_m(
    geometry: RotaryGeometry,
    load_n: float,
    sliding_speed_m_s: float,
    *,
    pressure_viscosity_pa: float = PRESSURE_VISCOSITY_COEFF_PA,
    reduced_modulus_pa: float = REDUCED_MODULUS_PA,
    viscosity_pa_s: float = 0.010,
) -> float:
    """Dowson-Higginson minimum EHL film thickness for the rotor-cylinder line contact.

    ``h_min/R = 2.65 U*^0.7 G*^0.54 W*^-0.13`` at equivalent radius ``R_eq = R_r R_cyl /
    (R_cyl - R_r) ~ 291`` mm (nearly conformal), entrainment ``ubar = |v|/2`` (the bore
    stationary). This is the *hard*-EHL film that rides in the elastically deflected
    contact; the film parameter ``Lambda = h_min / sigma`` then sets the lubrication
    regime. Returns ``+inf`` for no load / no sliding (no metal proximity).
    """

    if load_n <= 0.0 or sliding_speed_m_s <= 1.0e-12:
        return float("inf")
    radius_eq = geometry.rotor_radius_m * geometry.cylinder_radius_m / geometry.radial_clearance_m
    load_per_length = load_n / geometry.cylinder_height_m
    u_star = viscosity_pa_s * (0.5 * sliding_speed_m_s) / (reduced_modulus_pa * radius_eq)
    g_star = pressure_viscosity_pa * reduced_modulus_pa
    w_star = load_per_length / (reduced_modulus_pa * radius_eq)
    return radius_eq * 2.65 * u_star**0.7 * g_star**0.54 * w_star**-0.13


def _mixed_ehl_friction_force_n(
    geometry: RotaryGeometry,
    load_n: float,
    sliding_speed_m_s: float,
    *,
    roughness_m: float,
    pressure_viscosity_pa: float,
    reduced_modulus_pa: float,
    viscosity_pa_s: float,
    boundary_coefficient: float,
) -> float:
    """Mixed-lubrication friction *force* at one operating point (shared by the cycle
    average and the self-consistent orbit). Dowson-Higginson ``h_min`` -> film parameter
    -> Greenwood-Tripp-style asperity fraction + viscous shear over the Hertz zone.
    """

    if load_n <= 0.0 or sliding_speed_m_s <= 1.0e-12:
        return 0.0
    radius_eq = geometry.rotor_radius_m * geometry.cylinder_radius_m / geometry.radial_clearance_m
    length = geometry.cylinder_height_m
    load_per_length = load_n / length
    h_min = ehl_film_thickness_m(
        geometry,
        load_n,
        sliding_speed_m_s,
        pressure_viscosity_pa=pressure_viscosity_pa,
        reduced_modulus_pa=reduced_modulus_pa,
        viscosity_pa_s=viscosity_pa_s,
    )
    asperity_fraction = float(np.clip((3.0 - h_min / roughness_m) / 2.0, 0.0, 1.0))
    half_width = sqrt(4.0 * load_per_length * radius_eq / (pi * reduced_modulus_pa))
    viscous_force = viscosity_pa_s * sliding_speed_m_s / h_min * (2.0 * half_width * length)
    return boundary_coefficient * asperity_fraction * load_n + viscous_force


# ---------------------------------------------------------------------------
# Self-consistent sealing-contact rung (compliant contact solved in the rotor EOM)
# ---------------------------------------------------------------------------


def hertz_line_contact_force_n(
    approach_m: float,
    length_m: float,
    *,
    coeff: float = PALMGREN_LINE_CONTACT_COEFF,
) -> float:
    """Steel line-contact load from the elastic approach (Palmgren, inverted).

    ``delta[mm] = coeff * Q[N]^0.9 / L[mm]^0.8`` gives ``Q = (delta L^0.8 / coeff)^{1/0.9}``.
    Zero for a non-positive approach (the surfaces are apart). This is the *physical*
    contact stiffness that replaces the arbitrary penalty and closes ``N_c``.
    """

    if approach_m <= 0.0:
        return 0.0
    return (approach_m * 1.0e3 * (length_m * 1.0e3) ** 0.8 / coeff) ** (1.0 / 0.9)


@dataclass(frozen=True, slots=True)
class SealContactOrbit:
    """Rotor-cylinder sealing contact over the steady revolution (SI units)."""

    crank_angle_rad: tuple[float, ...]
    normal_force_n: tuple[float, ...]
    penetration_m: tuple[float, ...]  # signed: > 0 contact (Hertz), < 0 seal lifts (gap)
    sliding_speed_m_s: tuple[float, ...]
    boundary_friction_power_w: float
    mixed_ehl_friction_power_w: float
    mean_normal_force_n: float
    peak_normal_force_n: float
    max_penetration_m: float
    max_liftoff_m: float
    contact_fraction: float
    # EHL line-contact film in the (elastically deflected) contact -- sets the regime.
    mean_ehl_film_thickness_m: float
    mean_film_parameter: float  # Lambda = h_min / roughness; < 3 mixed, < 1 boundary
    friction_coefficient: float
    roughness_m: float


def integrate_sealing_contact(
    geometry: RotaryGeometry,
    *,
    trace: CycleTrace | None = None,
    revolutions: int = 3,
    samples: int = 360,
    grid_samples: int = 180,
    friction_coefficient: float = ROTOR_CYLINDER_FRICTION_COEFF,
    roughness_m: float = COMPOSITE_ROUGHNESS_M,
    rotor_mass_kg: float = ROTOR_MASS_KG,
    journal_clearance_m: float = JOURNAL_CLEARANCE_M,
    contact_smoothing_m: float = 3.0e-8,
    rtol: float = 1.0e-5,
) -> SealContactOrbit:
    """Solve the rotor-cylinder contact **self-consistently** in the rotor lateral EOM.

    The compliant Hertz line contact (:func:`hertz_line_contact_force_n`) is added to the
    Section 4.13 rotor EOM as a one-sided radial reaction plus sliding friction, and the
    whole system is integrated to its steady orbit. The load then distributes between the
    crank-pin journal film and the contact by their true stiffnesses -- no arbitrary
    penalty -- so ``N_c(theta)`` is determined and the ~6 um free-orbit penetration
    collapses to the physical Hertz deflection (sub-micron). Because the orbital
    **centrifugal** load (``m omega^2 e ~ 44 N``, always outward) presses the rotor onto
    the bore, the seal stays largely engaged; the gas modulates ``N_c`` and can briefly
    lift it (a leakage gap). The contact onset is smoothed over ``contact_smoothing_m``
    for the stiff solver; results are insensitive to it. Returns the boundary and
    mixed-EHL friction over the final revolution.
    """

    geometry.validate()
    if revolutions < 1 or samples < 8 or grid_samples < 16:
        raise ValueError("Need revolutions >= 1, samples >= 8, grid_samples >= 16.")
    if not (rotor_mass_kg > 0.0 and journal_clearance_m > 0.0 and roughness_m > 0.0):
        raise ValueError("Mass, clearance, and roughness must be positive.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    omega = geometry.angular_speed_rad_s
    throw = geometry.eccentricity_m
    length = geometry.cylinder_height_m
    grid, gas_x, gas_y, speed = _prescribe_gas_and_speed(geometry, grid_samples, trace)
    smoothing = contact_smoothing_m

    def contact_load(bore_radius: float) -> float:
        overlap = bore_radius - throw  # signed penetration
        softened = 0.5 * (overlap + hypot(overlap, smoothing))  # smooth one-sided max
        return hertz_line_contact_force_n(softened, length)

    def rhs(t: float, state: np.ndarray) -> list[float]:
        e_x, e_y, v_x, v_y = state
        theta = omega * t
        wrapped = theta % FULL_TURN_RAD
        gas_fx = float(np.interp(wrapped, grid, gas_x))
        gas_fy = float(np.interp(wrapped, grid, gas_y))
        entrainment = float(np.interp(wrapped, grid, speed))
        pin_x, pin_y = throw * sin(theta), throw * cos(theta)
        bore_x, bore_y = pin_x + e_x, pin_y + e_y
        bore_r = hypot(bore_x, bore_y)
        u_x, u_y = bore_x / bore_r, bore_y / bore_r
        load = contact_load(bore_r)
        slide = contact_sliding_speed_m_s(geometry, theta)
        tangent = -np.sign(slide)
        contact_x = -load * u_x - friction_coefficient * load * tangent * (-u_y)
        contact_y = -load * u_y - friction_coefficient * load * tangent * u_x
        film_x, film_y = _film_force_on_rotor(e_x, e_y, v_x, v_y, entrainment, journal_clearance_m)
        accel_x = (gas_fx + film_x + contact_x + rotor_mass_kg * omega**2 * pin_x) / rotor_mass_kg
        accel_y = (gas_fy + film_y + contact_y + rotor_mass_kg * omega**2 * pin_y) / rotor_mass_kg
        return [v_x, v_y, accel_x, accel_y]

    start_theta = grid[1]
    t_start = start_theta / omega
    period = 1.0 / geometry.frequency_hz
    t_end = t_start + revolutions * period
    t_eval = np.linspace(t_end - period, t_end, samples)
    t_eval[0] = max(t_eval[0], t_start)
    atol = np.array([1.0e-11, 1.0e-11, 1.0e-6, 1.0e-6])
    solution = solve_ivp(
        rhs,
        (t_start, t_end),
        [0.0, 0.0, 0.0, 0.0],
        method="BDF",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Sealing-contact integration failed: {solution.message}")

    angles = (omega * t_eval) % FULL_TURN_RAD
    bore_x = throw * np.sin(angles) + solution.y[0]
    bore_y = throw * np.cos(angles) + solution.y[1]
    penetration = np.hypot(bore_x, bore_y) - throw
    normal = np.array([contact_load(float(r)) for r in np.hypot(bore_x, bore_y)])
    slide = np.array([abs(contact_sliding_speed_m_s(geometry, float(a))) for a in angles])

    boundary = float(np.mean(friction_coefficient * normal * slide))
    mixed = float(
        np.mean(
            [
                _mixed_ehl_friction_force_n(
                    geometry,
                    float(n),
                    float(v),
                    roughness_m=roughness_m,
                    pressure_viscosity_pa=PRESSURE_VISCOSITY_COEFF_PA,
                    reduced_modulus_pa=REDUCED_MODULUS_PA,
                    viscosity_pa_s=0.010,
                    boundary_coefficient=friction_coefficient,
                )
                * v
                for n, v in zip(normal, slide, strict=True)
            ]
        )
    )
    films = np.array(
        [
            ehl_film_thickness_m(geometry, float(n), float(v))
            for n, v in zip(normal, slide, strict=True)
        ]
    )
    loaded = np.isfinite(films) & (normal > 1.0)
    mean_film = float(films[loaded].mean()) if loaded.any() else float("inf")
    return SealContactOrbit(
        crank_angle_rad=tuple(float(a) for a in angles),
        normal_force_n=tuple(float(n) for n in normal),
        penetration_m=tuple(float(p) for p in penetration),
        sliding_speed_m_s=tuple(float(v) for v in slide),
        boundary_friction_power_w=boundary,
        mixed_ehl_friction_power_w=mixed,
        mean_normal_force_n=float(normal.mean()),
        peak_normal_force_n=float(normal.max()),
        max_penetration_m=float(penetration.max()),
        max_liftoff_m=float(-penetration.min()),
        contact_fraction=float(np.mean(normal > 1.0)),
        mean_ehl_film_thickness_m=mean_film,
        mean_film_parameter=mean_film / roughness_m,
        friction_coefficient=friction_coefficient,
        roughness_m=roughness_m,
    )
