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


def solve_line_pressure(film: np.ndarray, source: np.ndarray, step: float) -> np.ndarray:
    """Solve ``d/dxi (film^3 dp/dxi) = source`` on a uniform grid, ``p = 0`` at the ends.

    ``film`` (``h > 0``) and ``source`` (``g``) are equal-length samples on a uniform
    grid of spacing ``step``. Returns the Gumbel-clamped pressure (``>= 0``), the same
    length, and zero at the two ends. The conservative flux discretisation with
    midpoint ``h^3`` faces is exact for a quadratic pressure and second-order in
    ``step`` for the variable-coefficient problem.
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
    np.maximum(pressure, 0.0, out=pressure)  # Gumbel (half-Sommerfeld) cavitation
    return pressure
