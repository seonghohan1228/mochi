"""Rotor lateral equation of motion and its time integration (D2 dynamics).

The project's **first real mechanical time integration**: the rotor bore no longer
sits exactly on the prescribed crank pin but moves within the crank-pin journal
clearance under gas load, oil-film reaction, and its own inertia. At constant drive
speed ``omega`` the state is the **bearing eccentricity** vector (Section 4.10)

    e = O_b - O_j = (e_x, e_y),   |e| <= c_j,

the rotor-bore centre relative to the prescribed crank pin ``O_j(theta) =
e_throw (sin theta, cos theta)``, ``theta = omega t``. With the rotor mass ``m_r``
the absolute rotor centre is ``O_b = O_j + e`` and Newton's second law gives

    m_r ( e'' - omega^2 O_j ) = F_gas(theta) + F_film(e, e'),

i.e. ``m_r e'' = F_gas + F_film + m_r omega^2 O_j`` -- the last term is the orbital
**centrifugal load** ``m_r e_throw omega^2`` the quasi-static model dropped (Section
4.8). ``F_gas(theta)`` is the mouth-aware gas force (Section 4.5); ``F_film`` is the
Ocvirk short-bearing reaction (Section 4.9) evaluated at the instantaneous
eccentricity and its rate -- so the squeeze term ``eps_dot`` is now **driven live**.

**Stiffness and the squeeze-film lag.** The squeeze film is strongly damping and
stiff, so the system is integrated with an implicit method (``scipy.integrate.
solve_ivp``, BDF). The centrifugal (inertia) load ``m_r e_throw omega^2 ~ 44 N`` is
only ~2 % of the ~2.5 kN gas reaction, so inertia is a small correction. The
**dominant dynamic effect is the squeeze-film lag**: the film's finite squeeze
response has a time constant that is a fixed fraction of the cycle (``tau/T ~ 0.26``,
independent of speed, mass, and clearance -- a geometric function of eccentricity),
so the rotor cannot follow the fast quasi-static swing. It **attenuates and
phase-lags** it -- peak eccentricity ~0.50 (vs the quasi-static ~0.71) and a larger
minimum film (~7.4 vs ~4.4 um). **Reduction check:** under a constant (frozen) load
the orbit relaxes exactly to the Section 4.9 quasi-static equilibrium; the
operating-speed deviation is the genuine squeeze-film effect, so the quasi-static
model **overestimates** the peak eccentricity.

**Friction.** Feeding the dynamic ``eps(theta)`` into the eccentric-friction law
(Section 4.9) gives the realistic cycle-mean journal loss ``~10.2 W`` -- between the
concentric Petroff (``~9.1 W``) and the quasi-static eccentric estimate (``~10.7 W``),
reported on ``RotorOrbit`` for comparison. See PHYSICS.md Section 4.13.
"""

from dataclasses import dataclass
from math import atan2, cos, hypot, sin

import numpy as np
from scipy.integrate import solve_ivp

from mochi.chambers import build_cycle_trace, seal_over_half_angle_rad
from mochi.journal_bearing import JOURNAL_CLEARANCE_M, journal_relative_speed_rad_s
from mochi.ocvirk_bearing import (
    crank_pin_entrainment_speed_rad_s,
    eccentric_friction_torque_nm,
    equilibrium_eccentricity_ratio,
    short_bearing_force,
)
from mochi.true_gas_force import FULL_TURN_RAD, true_gas_load

# Confirmed real rotor mass (2026-07), overriding the simplified-geometry density
# estimate of Section 4.8; editable per call.
ROTOR_MASS_KG = 0.275

# Largest eccentricity ratio the film force is evaluated at; nearer contact the
# closed form diverges, so the state is clamped here (the run should stay well below).
_MAX_ECCENTRICITY_RATIO = 0.995


@dataclass(frozen=True, slots=True)
class RotorOrbit:
    """Rotor lateral orbit over the final steady revolution (SI units)."""

    crank_angle_rad: tuple[float, ...]
    eccentricity_x_m: tuple[float, ...]
    eccentricity_y_m: tuple[float, ...]
    eccentricity_ratio: tuple[float, ...]
    min_film_thickness_m: tuple[float, ...]
    quasi_static_eccentricity_ratio: tuple[float, ...]
    peak_eccentricity_ratio: float
    minimum_film_thickness_m: float
    quasi_static_peak_eccentricity_ratio: float
    centrifugal_load_n: float
    dynamic_journal_friction_w: float
    quasi_static_journal_friction_w: float
    petroff_journal_friction_w: float
    mass_kg: float
    clearance_m: float
    revolutions: int
    frequency_hz: float


def _film_force_on_rotor(
    e_x: float,
    e_y: float,
    v_x: float,
    v_y: float,
    entrainment_speed_rad_s: float,
    clearance_m: float,
) -> tuple[float, float]:
    """Ocvirk short-bearing oil-film reaction on the rotor, in the global frame.

    Maps the Cartesian state to the ``(eps, Omega, eps_dot)`` the closed form uses:
    the attitude (whirl) rate ``psi_dot`` reduces the entrainment, and ``eps_dot`` is
    the eccentricity-magnitude rate. The force on the journal is
    ``F_e e_hat + F_t (z x e_hat)``; the reaction on the rotor is its negative.
    """

    magnitude = hypot(e_x, e_y)
    if magnitude <= 1.0e-12:
        return 0.0, 0.0
    eccentricity_ratio = min(magnitude / clearance_m, _MAX_ECCENTRICITY_RATIO)
    radial_hat = (e_x / magnitude, e_y / magnitude)
    tangent_hat = (-radial_hat[1], radial_hat[0])  # z x e_hat
    eccentricity_rate = (e_x * v_x + e_y * v_y) / (magnitude * clearance_m)
    attitude_rate = (e_x * v_y - e_y * v_x) / (magnitude * magnitude)

    force = short_bearing_force(
        eccentricity_ratio,
        entrainment_speed_rad_s - attitude_rate,
        eccentricity_rate_per_s=eccentricity_rate,
        clearance_m=clearance_m,
    )
    f_e, f_t = force.radial_force_n, force.tangential_force_n
    on_rotor_x = -(f_e * radial_hat[0] + f_t * tangent_hat[0])
    on_rotor_y = -(f_e * radial_hat[1] + f_t * tangent_hat[1])
    return on_rotor_x, on_rotor_y


def _prescribe_gas_and_speed(
    geometry, samples, trace
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample F_gas(theta) and the entrainment speed on a periodic grid.

    The gas force is undefined in the small seal-over window (chambers merge), so the
    grid runs over the separated range and the seal-over seam is bridged linearly by
    a wrap point at theta = 0 / 2*pi -- adequate because the window is a few degrees
    and the rotor is overdamped.
    """

    half_angle = seal_over_half_angle_rad(geometry)
    inner = np.linspace(half_angle, FULL_TURN_RAD - half_angle, samples)
    gas_x = np.empty(samples)
    gas_y = np.empty(samples)
    for index, angle in enumerate(inner):
        force = true_gas_load(geometry, float(angle), trace=trace).rotor_force_n
        gas_x[index] = force[0]
        gas_y[index] = force[1]
    # Bridge the seal-over seam with a wrap point at 0 and 2*pi (mean of the ends).
    bridge_x = 0.5 * (gas_x[0] + gas_x[-1])
    bridge_y = 0.5 * (gas_y[0] + gas_y[-1])
    grid = np.concatenate(([0.0], inner, [FULL_TURN_RAD]))
    gx = np.concatenate(([bridge_x], gas_x, [bridge_x]))
    gy = np.concatenate(([bridge_y], gas_y, [bridge_y]))

    speed = np.array([crank_pin_entrainment_speed_rad_s(geometry, float(a)) for a in grid])
    return grid, gx, gy, speed


def integrate_rotor_orbit(
    geometry,
    *,
    mass_kg: float = ROTOR_MASS_KG,
    clearance_m: float = JOURNAL_CLEARANCE_M,
    revolutions: int = 4,
    samples: int = 360,
    grid_samples: int = 360,
    trace=None,
    rtol: float = 1.0e-6,
    atol: float = 1.0e-12,
) -> RotorOrbit:
    """Integrate the rotor lateral EOM to its steady periodic orbit.

    Runs ``revolutions`` crank turns from the quasi-static equilibrium at the start
    angle (the overdamped transient decays within ~1 turn) and returns the final
    revolution sampled at ``samples`` points, alongside the quasi-static eccentricity
    for the reduction check. Uses an implicit stiff integrator (BDF).
    """

    geometry.validate()
    if not (mass_kg > 0.0 and clearance_m > 0.0):
        raise ValueError("Mass and clearance must be positive.")
    if revolutions < 1 or samples < 8 or grid_samples < 16:
        raise ValueError("Need revolutions >= 1, samples >= 8, grid_samples >= 16.")

    trace = build_cycle_trace(geometry) if trace is None else trace
    omega = geometry.angular_speed_rad_s
    throw = geometry.eccentricity_m
    grid, gas_x, gas_y, speed = _prescribe_gas_and_speed(geometry, grid_samples, trace)

    def gas_at(theta: float) -> tuple[float, float]:
        wrapped = theta % FULL_TURN_RAD
        return float(np.interp(wrapped, grid, gas_x)), float(np.interp(wrapped, grid, gas_y))

    def speed_at(theta: float) -> float:
        return float(np.interp(theta % FULL_TURN_RAD, grid, speed))

    def rhs(t: float, state: np.ndarray) -> list[float]:
        e_x, e_y, v_x, v_y = state
        theta = omega * t
        gx, gy = gas_at(theta)
        entrainment = speed_at(theta)
        pin_x, pin_y = throw * sin(theta), throw * cos(theta)
        film_x, film_y = _film_force_on_rotor(e_x, e_y, v_x, v_y, entrainment, clearance_m)
        accel_x = (gx + film_x + mass_kg * omega * omega * pin_x) / mass_kg
        accel_y = (gy + film_y + mass_kg * omega * omega * pin_y) / mass_kg
        return [v_x, v_y, accel_x, accel_y]

    # Quasi-static equilibrium at the start angle for the initial condition and the
    # reduction-check reference (film force -F_gas => e at attitude F_gas_angle - psi).
    def quasi_static_state(theta: float) -> tuple[float, float, float]:
        gx, gy = gas_at(theta)
        load = hypot(gx, gy)
        entrainment = speed_at(theta)
        eps = equilibrium_eccentricity_ratio(load, entrainment, clearance_m=clearance_m)
        if eps <= 0.0:
            return 0.0, 0.0, 0.0
        psi = short_bearing_force(eps, entrainment, clearance_m=clearance_m).attitude_angle_rad
        direction = atan2(gy, gx) - psi
        return eps, eps * clearance_m * cos(direction), eps * clearance_m * sin(direction)

    start_theta = grid[1]  # first separated angle
    _, e0_x, e0_y = quasi_static_state(start_theta)
    t_start = start_theta / omega
    period = 1.0 / geometry.frequency_hz
    t_end = t_start + revolutions * period
    t_eval = np.linspace(t_end - period, t_end, samples)

    solution = solve_ivp(
        rhs,
        (t_start, t_end),
        [e0_x, e0_y, 0.0, 0.0],
        method="BDF",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Rotor orbit integration failed: {solution.message}")

    ex = solution.y[0]
    ey = solution.y[1]
    angles = (omega * t_eval) % FULL_TURN_RAD
    ecc_ratio = np.hypot(ex, ey) / clearance_m
    films = clearance_m * (1.0 - np.minimum(ecc_ratio, 1.0))
    quasi = np.array([quasi_static_state(float(a))[0] for a in angles])

    # Journal friction from the eccentric Couette law (Section 4.9). The t_eval grid
    # is uniform over one period, so a simple sample mean is the cycle mean. Compare
    # the dynamic eccentricity to the quasi-static estimate and the concentric Petroff.
    omega_rel = np.array([journal_relative_speed_rad_s(geometry, float(a)) for a in angles])

    def friction_mean(eps_series: np.ndarray) -> float:
        total = 0.0
        for eps, speed in zip(eps_series, omega_rel, strict=True):
            capped = min(float(eps), _MAX_ECCENTRICITY_RATIO)
            total += eccentric_friction_torque_nm(capped, float(speed)) * abs(float(speed))
        return total / len(omega_rel)

    return RotorOrbit(
        crank_angle_rad=tuple(float(a) for a in angles),
        eccentricity_x_m=tuple(float(v) for v in ex),
        eccentricity_y_m=tuple(float(v) for v in ey),
        eccentricity_ratio=tuple(float(v) for v in ecc_ratio),
        min_film_thickness_m=tuple(float(v) for v in films),
        quasi_static_eccentricity_ratio=tuple(float(v) for v in quasi),
        peak_eccentricity_ratio=float(ecc_ratio.max()),
        minimum_film_thickness_m=float(films.min()),
        quasi_static_peak_eccentricity_ratio=float(quasi.max()),
        centrifugal_load_n=mass_kg * omega * omega * throw,
        dynamic_journal_friction_w=friction_mean(ecc_ratio),
        quasi_static_journal_friction_w=friction_mean(quasi),
        petroff_journal_friction_w=friction_mean(np.zeros_like(ecc_ratio)),
        mass_kg=mass_kg,
        clearance_m=clearance_m,
        revolutions=revolutions,
        frequency_hz=geometry.frequency_hz,
    )
