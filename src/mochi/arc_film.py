"""Partial-arc oil-film force (curved swing-bush piece film), axial-uniform 1-D.

Each swing-bush piece is a **partial** journal: its convex OD arc runs in the
concave rotor groove over a limited span (the piece is not a full ring, Section
3.3/4.11). This gives the 2-D film reaction on the piece as a function of the piece
eccentricity ``(e_x, e_y) = O_p - O_g`` (piece centre relative to groove centre) and
its rate.

The pressure is solved with the **axial-uniform (long-bearing) 1-D Reynolds**
reduction: at the bush aspect ratio the film escapes **circumferentially**, not
axially (``L/D > 1``), so ``dp/dz ~ 0`` and the pressure is solved *along the arc*
and multiplied by the axial height ``H``. This is the opposite reduction to the
axial short-bearing (:mod:`mochi.reynolds_1d`, Ocvirk), which is the correct one for
the ``L/D < 1`` crank-pin journal. Film thickness

    h(beta) = c - (e_x cos beta + e_y sin beta),

with ``beta`` the angle from the groove centre; along the arc-length ``x = r beta``
the Reynolds equation ``d/dx (h^3/(12 mu) dp/dx) = ubar dh/dx + dh/dt`` becomes

    d/dbeta ( h^3 dp/dbeta ) = 12 mu r^2 S(beta),
    S(beta) = Omega (e_x sin beta - e_y cos beta) - (edot_x cos beta + edot_y sin beta),

``Omega = ubar / r`` the entrainment (mean surface) angular speed of the rotor groove
and piece. The arc is open to oil at both circumferential edges, so ``p = 0`` there;
the line problem is solved by :func:`mochi.line_reynolds.solve_line_pressure` with
Gumbel (half-Sommerfeld) cavitation. The pressure acts radially inward on the piece
over the wetted area ``dA = r H dbeta``, so

    F_x = - r H integral p cos beta dbeta,
    F_y = - r H integral p sin beta dbeta.

The radial pressure passes through the arc's centre of curvature ``O_p``, so it
contributes **no moment** about ``O_p`` (only the shear does, second order and
omitted here); the piece attitude is driven instead by the flat-film moment
(:mod:`mochi.slider_film`). Keeping a **partial** arc preserves the tangential
asymmetry a full ring averages out. See PHYSICS.md Section 4.11.
"""

from dataclasses import dataclass
from math import isfinite

import numpy as np

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S
from mochi.line_reynolds import poiseuille_bias_pressure, solve_line_pressure


@dataclass(frozen=True, slots=True)
class ArcFilmForce:
    """Partial-arc film reaction on a bush piece (global frame, SI units)."""

    force_x_n: float
    force_y_n: float
    min_film_thickness_m: float
    max_pressure_pa: float
    throughflow_m3_s: float = 0.0  # pressure-driven (Poiseuille) leakage across the arc
    shear_moment_nm: float = 0.0  # Couette shear moment on the piece about O_p (rotor-swing drag)


def arc_film_force(
    ecc_x_m: float,
    ecc_y_m: float,
    ecc_dot_x_m_s: float,
    ecc_dot_y_m_s: float,
    entrainment_speed_rad_s: float,
    *,
    arc_center_rad: float,
    arc_half_span_rad: float,
    radius_m: float,
    length_m: float,
    clearance_m: float,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    n_beta: int = 361,
    pressure_start_pa: float = 0.0,
    pressure_end_pa: float = 0.0,
    cavitation_pressure_pa: float = 0.0,
    shear_speed_rad_s: float = 0.0,
    check_physical: bool = False,
    reference_pressure_pa: float = 0.0,
) -> ArcFilmForce:
    """Axial-uniform 1-D film reaction on a partial-arc (bush) piece.

    ``ecc = (e_x, e_y) = O_p - O_g`` is the piece-centre offset from the groove
    centre; the arc spans ``arc_center +/- arc_half_span``; ``length_m`` is the axial
    height ``H``. The pressure is solved along the arc (axial-uniform) and integrated
    with the height, so the returned 2-D force on the piece is ``O(H)``; its reaction
    ``-F`` acts on the rotor groove.

    ``pressure_start_pa`` / ``pressure_end_pa`` are the gas pressures the arc ends open
    to (the local chamber / crank-bore pressures, gauge; the arc runs from
    ``arc_center - arc_half_span`` to ``arc_center + arc_half_span``). They add the
    Poiseuille gas-bias field (:func:`mochi.line_reynolds.poiseuille_bias_pressure`) to
    the hydrodynamic pressure by superposition; both default to zero (the pure
    hydrodynamic film). ``throughflow_m3_s`` reports the resulting pressure-driven
    leakage ``H (p_start - p_end)/(12 mu integral h^-3 dx)``.

    ``cavitation_pressure_pa`` is the pressure below which the film cavitates, in the same
    gauge as the end pressures; the clamp is applied to the **total** (hydrodynamic + gas)
    field. The default ``0.0`` is classic Gumbel (half-Sommerfeld) cavitation.
    """

    for value in (ecc_x_m, ecc_y_m, ecc_dot_x_m_s, ecc_dot_y_m_s, entrainment_speed_rad_s):
        if not isfinite(value):
            raise ValueError("Eccentricity, its rate, and the speed must be finite.")
    if not (0.0 < arc_half_span_rad <= np.pi):
        raise ValueError("Arc half-span must be in (0, pi].")
    for name, value in (("radius", radius_m), ("length", length_m), ("clearance", clearance_m)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Piece {name} must be a positive, finite length in metres.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    if n_beta < 16:
        raise ValueError("The arc grid needs at least 16 points.")

    beta = np.linspace(
        arc_center_rad - arc_half_span_rad, arc_center_rad + arc_half_span_rad, n_beta
    )
    cos_b, sin_b = np.cos(beta), np.sin(beta)
    film = clearance_m - (ecc_x_m * cos_b + ecc_y_m * sin_b)
    if np.any(film <= 0.0):
        raise ValueError("Piece eccentricity exceeds the clearance (metal contact).")

    # S(beta) = ubar dh/dx + dh/dt with x = r*beta and Omega = ubar/r; the scaled line
    # source is g = 12 mu r^2 S (see module docstring).
    source = entrainment_speed_rad_s * (ecc_x_m * sin_b - ecc_y_m * cos_b) - (
        ecc_dot_x_m_s * cos_b + ecc_dot_y_m_s * sin_b
    )
    line_source = 12.0 * viscosity_pa_s * radius_m**2 * source

    d_beta = 2.0 * arc_half_span_rad / (n_beta - 1)
    gas_present = pressure_start_pa != 0.0 or pressure_end_pa != 0.0
    # Reynolds is linear in p, so the hydrodynamic field is solved UNCLAMPED and the gas bias
    # superposed; cavitation is then applied to the **total**, which is where the physical
    # constraint lives (see mochi.slider_film for the same treatment).
    pressure = solve_line_pressure(film, line_source, d_beta, cavitation=False)

    throughflow = 0.0
    if gas_present:
        gas = poiseuille_bias_pressure(film, d_beta, pressure_start_pa, pressure_end_pa)
        pressure = pressure + gas  # superposition (Reynolds linear in p)
        # Pressure-driven throughflow (Poiseuille slot leakage) along the arc, arc
        # length element dx = r dbeta, so the h^-3 integral carries the r factor.
        inv_h3_integral = float(np.trapezoid(film**-3.0, dx=radius_m * d_beta))
        throughflow = (
            length_m * (pressure_start_pa - pressure_end_pa)
            / (12.0 * viscosity_pa_s * inv_h3_integral)
        )
    np.maximum(pressure, cavitation_pressure_pa, out=pressure)
    if check_physical:
        from mochi.physical_checks import check_film_positive, check_pressure_field

        check_film_positive(film, label="arc film")
        check_pressure_field(
            pressure, reference_pressure_pa=reference_pressure_pa, label="arc film"
        )

    force_x = -radius_m * length_m * float(np.trapezoid(pressure * cos_b, dx=d_beta))
    force_y = -radius_m * length_m * float(np.trapezoid(pressure * sin_b, dx=d_beta))
    # Couette shear moment on the piece about its centre O_p: the groove surface shears
    # the film at the relative rotation ``shear_speed_rad_s``, dragging the piece. The
    # tangential traction mu*(shear_speed*r)/h acts at lever r over dA = r H dbeta, so
    # M = mu * shear_speed * r^3 * H * integral(dbeta / h). This is the (viscous) moment
    # that couples the piece rotation to the rotor swing; zero when shear_speed = 0.
    shear_moment = (
        viscosity_pa_s * shear_speed_rad_s * radius_m**3 * length_m
        * float(np.trapezoid(1.0 / film, dx=d_beta))
    )
    return ArcFilmForce(
        force_x_n=force_x,
        force_y_n=force_y,
        min_film_thickness_m=float(np.min(film)),
        max_pressure_pa=float(np.max(pressure)),
        throughflow_m3_s=throughflow,
        shear_moment_nm=shear_moment,
    )
