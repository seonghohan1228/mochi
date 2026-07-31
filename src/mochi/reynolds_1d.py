"""1-D finite-width (short-bearing) Reynolds solver for the crank-pin journal film.

Step (i) of the Section 4.10 1D->2D ladder: numerically integrate the **axial**
(short-bearing) Reynolds equation on a finite-difference grid with **Gumbel
(half-Sommerfeld) cavitation**, and validate the result against the Section 4.9
Ocvirk closed form. Dropping the circumferential pressure-flow term leaves, at each
circumferential angle ``beta``, a 1-D axial problem

    d/dz ( h^3/(12 mu) dp/dz ) = S(beta),
    S(beta) = c ( eps_dot cos beta - eps Omega sin beta ),   p(+/- L/2) = 0,

with the film thickness ``h(beta) = c (1 + eps cos beta)`` measured, as in Section
4.9, from the point of maximum film thickness. Because ``h`` does not vary along
``z`` for an aligned journal the axial solution is a parabola; solving it by a
tridiagonal finite difference (exact for the quadratic) and integrating the
pressure over the wetted surface **must reproduce the Ocvirk force**. The value of
this rung is that it exercises the numerical machinery -- tridiagonal axial solve,
cavitation clamp, and pressure-force integration -- that the later **2-D**
(circumferential + axial + tilt) Reynolds solver reuses; here that machinery is
checked against a closed form. The force is resolved, like Section 4.9, into a
component along the line of centres (``F_e``, toward maximum film) and perpendicular
(``F_t``):

    F_e = - r oint ( int p dz ) cos beta d beta,
    F_t = - r oint ( int p dz ) sin beta d beta.

Geometry and lubricant default to the crank-pin journal (Section 4.7). Cavitation
is Gumbel here to match the Section 4.9 pi-film; the physical Swift-Stieber
boundary and the asymmetric/4 MPa reference (Section 4.10) belong to the 2-D rung.
See PHYSICS.md Section 4.12.
"""

from dataclasses import dataclass
from math import isfinite

import numpy as np

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S
from mochi.journal_bearing import (
    JOURNAL_CLEARANCE_M,
    JOURNAL_LENGTH_M,
    JOURNAL_RADIUS_M,
)


@dataclass(frozen=True, slots=True)
class ShortBearingSolution:
    """Numerical short-bearing (1-D axial Reynolds) journal force (SI units).

    ``radial_force_n`` (``F_e``) is along the line of centres toward maximum film;
    ``tangential_force_n`` (``F_t``) is perpendicular -- the same convention as
    :func:`mochi.ocvirk_bearing.short_bearing_force`.
    """

    eccentricity_ratio: float
    entrainment_speed_rad_s: float
    eccentricity_rate_per_s: float
    radial_force_n: float
    tangential_force_n: float
    magnitude_n: float
    n_circumferential: int
    n_axial: int


def _simpson(values: np.ndarray, step: float) -> np.ndarray:
    """Composite Simpson integral of ``values`` along axis 0 (odd sample count)."""

    n = values.shape[0]
    weights = np.ones(n)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    weights *= step / 3.0
    return np.tensordot(weights, values, axes=([0], [0]))


def solve_short_bearing_1d(
    eccentricity_ratio: float,
    entrainment_speed_rad_s: float,
    *,
    eccentricity_rate_per_s: float = 0.0,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    radius_m: float = JOURNAL_RADIUS_M,
    length_m: float = JOURNAL_LENGTH_M,
    clearance_m: float = JOURNAL_CLEARANCE_M,
    n_circumferential: int = 720,
    n_axial: int = 41,
) -> ShortBearingSolution:
    """Solve the short-bearing (axial 1-D) Reynolds film numerically.

    Returns the journal force in the ``(F_e, F_t)`` convention of
    :func:`mochi.ocvirk_bearing.short_bearing_force`, which it reproduces (Gumbel
    cavitation). ``n_axial`` must be odd (composite Simpson).
    """

    if not (isfinite(eccentricity_ratio) and 0.0 <= eccentricity_ratio < 1.0):
        raise ValueError("Eccentricity ratio must be finite and in [0, 1).")
    if not isfinite(entrainment_speed_rad_s):
        raise ValueError("Entrainment speed must be finite.")
    if not isfinite(eccentricity_rate_per_s):
        raise ValueError("Eccentricity rate must be finite.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    for name, value in (("radius", radius_m), ("length", length_m), ("clearance", clearance_m)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Journal {name} must be a positive, finite length in metres.")
    if n_circumferential < 16:
        raise ValueError("The circumferential grid needs at least 16 points.")
    if n_axial < 5 or n_axial % 2 == 0:
        raise ValueError("The axial grid needs an odd number of points (>= 5) for Simpson.")

    beta = np.linspace(0.0, 2.0 * np.pi, n_circumferential, endpoint=False)
    cos_b, sin_b = np.cos(beta), np.sin(beta)
    film = clearance_m * (1.0 + eccentricity_ratio * cos_b)
    source = clearance_m * (
        eccentricity_rate_per_s * cos_b - eccentricity_ratio * entrainment_speed_rad_s * sin_b
    )
    # Axial ODE p'' = 12 mu S / h^3 (constant in z at each beta), p = 0 at z = +/- L/2.
    curvature = 12.0 * viscosity_pa_s * source / film**3  # (n_circumferential,)

    n_interior = n_axial - 2
    step = length_m / (n_axial - 1)
    tridiag = (
        np.diag(np.full(n_interior, -2.0))
        + np.diag(np.ones(n_interior - 1), 1)
        + np.diag(np.ones(n_interior - 1), -1)
    )
    rhs = (step * step * curvature)[None, :] * np.ones((n_interior, 1))
    interior = np.linalg.solve(tridiag, rhs)  # (n_interior, n_circumferential)

    pressure = np.zeros((n_axial, n_circumferential))
    pressure[1:-1, :] = interior
    np.maximum(pressure, 0.0, out=pressure)  # Gumbel (half-Sommerfeld) cavitation

    axial_integral = _simpson(pressure, step)  # (n_circumferential,)
    d_beta = 2.0 * np.pi / n_circumferential
    radial = -radius_m * d_beta * float(np.sum(axial_integral * cos_b))
    tangential = -radius_m * d_beta * float(np.sum(axial_integral * sin_b))

    return ShortBearingSolution(
        eccentricity_ratio=eccentricity_ratio,
        entrainment_speed_rad_s=entrainment_speed_rad_s,
        eccentricity_rate_per_s=eccentricity_rate_per_s,
        radial_force_n=radial,
        tangential_force_n=tangential,
        magnitude_n=float(np.hypot(radial, tangential)),
        n_circumferential=n_circumferential,
        n_axial=n_axial,
    )
