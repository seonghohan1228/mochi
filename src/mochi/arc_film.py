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
from mochi.line_reynolds import solve_line_pressure


@dataclass(frozen=True, slots=True)
class ArcFilmForce:
    """Partial-arc film reaction on a bush piece (global frame, SI units)."""

    force_x_n: float
    force_y_n: float
    min_film_thickness_m: float
    max_pressure_pa: float


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
) -> ArcFilmForce:
    """Axial-uniform 1-D film reaction on a partial-arc (bush) piece.

    ``ecc = (e_x, e_y) = O_p - O_g`` is the piece-centre offset from the groove
    centre; the arc spans ``arc_center +/- arc_half_span``; ``length_m`` is the axial
    height ``H``. The pressure is solved along the arc (axial-uniform) and integrated
    with the height, so the returned 2-D force on the piece is ``O(H)``; its reaction
    ``-F`` acts on the rotor groove.
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
    pressure = solve_line_pressure(film, line_source, d_beta)  # Pa, Gumbel-clamped

    force_x = -radius_m * length_m * float(np.trapezoid(pressure * cos_b, dx=d_beta))
    force_y = -radius_m * length_m * float(np.trapezoid(pressure * sin_b, dx=d_beta))
    return ArcFilmForce(
        force_x_n=force_x,
        force_y_n=force_y,
        min_film_thickness_m=float(np.min(film)),
        max_pressure_pa=float(np.max(pressure)),
    )
