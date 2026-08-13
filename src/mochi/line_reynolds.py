"""1-D Reynolds *line* solver: axial-uniform pressure solved along the film line.

Shared machinery for the swing-bush films (:mod:`mochi.arc_film`,
:mod:`mochi.slider_film`). Under the **axial-uniform (long-bearing)** assumption the
pressure does not vary across the axial width -- it is solved as a 1-D boundary-value
problem *along the film line* (circumferential ``beta`` for the curved arc, sliding
``s`` for the flat pad) and then multiplied by the axial height ``H`` for the force.

This is the reduction requested for the bush films: at their aspect ratio (film-line
length short, axial width ``H`` long, ``L/D > 1``) the Poiseuille pressure escapes
**along the film line**, not axially, so ``dp/dz ~ 0`` and the 1-D problem is the
line direction. It is the *opposite* reduction to the axial short-bearing
(:mod:`mochi.reynolds_1d`, Ocvirk) that is correct for the ``L/D < 1`` crank-pin
journal, where the pressure escapes axially instead.

The line problem is the variable-coefficient two-point BVP

    d/dxi ( h^3 dp/dxi ) = g(xi),    p = 0 at both ends,

with ``h(xi)`` the film thickness and ``g(xi)`` the (already-scaled) Reynolds source.
It is discretised in **conservative flux form** on a uniform grid -- fluxes
``h^3 dp/dxi`` evaluated at the cell faces with the midpoint ``h`` -- giving a
tridiagonal system solved by :func:`scipy.linalg.solve_banded`. **Gumbel
(half-Sommerfeld) cavitation** then clamps ``p >= 0`` pointwise. See PHYSICS.md
Section 4.11.
"""

import numpy as np
from scipy.linalg import solve_banded


def solve_line_pressure(
    film: np.ndarray, source: np.ndarray, step: float, *, cavitation: bool = True
) -> np.ndarray:
    """Solve ``d/dxi (film^3 dp/dxi) = source`` on a uniform grid, ``p = 0`` at the ends.

    ``film`` (``h > 0``) and ``source`` (``g``) are equal-length samples on a uniform
    grid of spacing ``step``. Returns the pressure (same length, zero at the two ends).
    The conservative flux discretisation with midpoint ``h^3`` faces is exact for a
    quadratic pressure and second-order in ``step`` for the variable-coefficient problem.

    ``cavitation`` (default ``True``) applies **Gumbel (half-Sommerfeld) cavitation**,
    clamping ``p >= 0`` -- correct when the film ambient is near the cavitation pressure.
    Set it ``False`` for the **full-Sommerfeld** field (no clamp): the physical choice
    when a high gas pressure floods the film and suppresses cavitation (the superposed
    total stays well above the cavitation pressure), so the diverging half carries a real
    sub-ambient pressure rather than cavitating.
    """

    film = np.asarray(film, dtype=float)
    source = np.asarray(source, dtype=float)
    if film.shape != source.shape or film.ndim != 1:
        raise ValueError("film and source must be 1-D arrays of the same length.")
    n = film.size
    if n < 3:
        raise ValueError("The line grid needs at least three points.")
    if not (np.all(np.isfinite(film)) and np.all(film > 0.0)):
        raise ValueError("Film thickness must be finite and positive everywhere.")
    if not np.all(np.isfinite(source)):
        raise ValueError("Source must be finite everywhere.")
    if not (np.isfinite(step) and step > 0.0):
        raise ValueError("Grid step must be a positive, finite length.")

    # h^3 at the interior cell faces (xi_{i+1/2}), length n - 1.
    faces = (0.5 * (film[:-1] + film[1:])) ** 3
    lower = faces[:-1]  # h^3_{i-1/2} for interior nodes i = 1 .. n-2
    upper = faces[1:]  # h^3_{i+1/2}
    diag = -(lower + upper)

    # Banded storage for solve_banded((1, 1), ...): row 0 super-, 1 diag, 2 sub-diagonal.
    ab = np.zeros((3, n - 2))
    ab[0, 1:] = upper[:-1]
    ab[1, :] = diag
    ab[2, :-1] = lower[1:]
    rhs = step * step * source[1:-1]

    pressure = np.zeros(n)
    pressure[1:-1] = solve_banded((1, 1), ab, rhs)
    if cavitation:
        np.maximum(pressure, 0.0, out=pressure)  # Gumbel (half-Sommerfeld) cavitation
    return pressure


def poiseuille_bias_pressure(
    film: np.ndarray, step: float, pressure_start: float, pressure_end: float
) -> np.ndarray:
    """Homogeneous (source-free) Reynolds field driven by unequal end pressures.

    The film ends open to gas at ``pressure_start`` / ``pressure_end`` (the local chamber
    or crank/bore pressures, in gauge). With no wedge/squeeze source the Reynolds equation
    ``d/dxi(h^3 dp/dxi) = 0`` has ``h^3 dp/dxi = Q`` constant, so

        p(xi) = p_start + (p_end - p_start) * C(xi) / C(L),   C(xi) = integral_0^xi h^-3 dxi',

    a pure pressure-driven (Poiseuille) throughflow field. Because the Reynolds equation
    is **linear in pressure**, the full film pressure is the superposition
    ``p = p_hydro (this module's :func:`solve_line_pressure`, p=0 ends) + p_gas`` -- valid
    while the total stays above the cavitation pressure, which the MPa-scale chamber
    pressures ensure. Returns the gas-bias pressure on the same grid; the end values are
    exactly ``pressure_start`` / ``pressure_end``. See PHYSICS.md Section 4.11.
    """

    film = np.asarray(film, dtype=float)
    if film.ndim != 1 or film.size < 3:
        raise ValueError("Film must be a 1-D array with at least three points.")
    if not (np.all(np.isfinite(film)) and np.all(film > 0.0)):
        raise ValueError("Film thickness must be finite and positive everywhere.")
    if not (np.isfinite(step) and step > 0.0):
        raise ValueError("Grid step must be a positive, finite length.")
    if not (np.isfinite(pressure_start) and np.isfinite(pressure_end)):
        raise ValueError("End pressures must be finite.")

    inv_h3 = film**-3.0
    # Cumulative trapezoidal integral C(xi) of h^-3, zero at the first node.
    cumulative = np.concatenate(([0.0], np.cumsum(0.5 * (inv_h3[:-1] + inv_h3[1:]) * step)))
    total = cumulative[-1]
    if total <= 0.0:  # degenerate; fall back to a linear ramp
        fraction = np.linspace(0.0, 1.0, film.size)
    else:
        fraction = cumulative / total
    return pressure_start + (pressure_end - pressure_start) * fraction
