"""1-D Reynolds line solver: manufactured-solution and convergence checks."""

import numpy as np
import pytest

from mochi.line_reynolds import poiseuille_bias_pressure, solve_line_pressure


def _grid(n, length):
    xi = np.linspace(0.0, length, n)
    return xi, length / (n - 1)


def test_constant_coefficient_is_exact_quadratic() -> None:
    """Constant film + constant source: the solution is a parabola, recovered exactly."""

    n, length, h0, g = 401, 4.0e-3, 10.0e-6, -2.0e-3
    xi, step = _grid(n, length)
    film = np.full(n, h0)
    source = np.full(n, g)
    # h0^3 p'' = g, p(0) = p(L) = 0  ->  p = g/(2 h0^3) * xi (xi - L)  (> 0 for g < 0).
    exact = g / (2.0 * h0**3) * xi * (xi - length)

    pressure = solve_line_pressure(film, source, step)
    assert np.all(pressure >= 0.0)
    assert np.allclose(pressure, exact, rtol=1e-9, atol=1e-9 * exact.max())


def _manufactured(n):
    """Variable-h manufactured case with p = sin(pi xi / L) >= 0 (no cavitation)."""

    length, h0 = 4.0e-3, 10.0e-6
    xi, step = _grid(n, length)
    k = np.pi / length
    film = h0 * (1.0 + 0.3 * np.cos(k * xi))
    film_prime = -0.3 * h0 * k * np.sin(k * xi)
    p_exact = np.sin(k * xi)
    dp = k * np.cos(k * xi)
    d2p = -k * k * np.sin(k * xi)
    source = 3.0 * film**2 * film_prime * dp + film**3 * d2p  # g = d/dxi(h^3 dp/dxi)
    p_num = solve_line_pressure(film, source, step)
    return float(np.max(np.abs(p_num - p_exact)))


def test_variable_coefficient_matches_and_converges() -> None:
    """The variable-h solve matches the manufactured p and is ~second-order accurate."""

    err_coarse = _manufactured(251)
    err_fine = _manufactured(501)
    assert err_fine < 1.0e-3
    # halving the step should cut the error by ~4 (second order); require > 3.5.
    assert err_coarse / err_fine > 3.5


def test_all_negative_source_cavitates_to_zero() -> None:
    """A source that drives p < 0 everywhere is Gumbel-clamped to zero (no load)."""

    n, length, h0 = 201, 3.0e-3, 8.0e-6
    _, step = _grid(n, length)
    film = np.full(n, h0)
    pressure = solve_line_pressure(film, np.full(n, +5.0e-3), step)  # g > 0 -> p < 0
    assert np.all(pressure == 0.0)


def test_rejects_bad_inputs() -> None:
    good = np.full(5, 1.0e-5)
    with pytest.raises(ValueError):
        solve_line_pressure(good, good[:4], 1e-4)  # length mismatch
    with pytest.raises(ValueError):
        solve_line_pressure(np.array([1e-5, -1e-5, 1e-5]), np.zeros(3), 1e-4)  # film <= 0
    with pytest.raises(ValueError):
        solve_line_pressure(good, good, 0.0)  # non-positive step


# ---------------------------------------------------------------------------
# Gas-pressure boundary conditions (Poiseuille bias, PHYSICS.md 4.11)
# ---------------------------------------------------------------------------


def test_poiseuille_bias_uniform_film_is_linear() -> None:
    """A uniform gap gives a linear pressure ramp between the two end pressures."""

    n, length = 201, 0.01
    _, step = _grid(n, length)
    film = np.full(n, 1.0e-5)
    p = poiseuille_bias_pressure(film, step, 2.0e6, 5.0e5)
    assert p[0] == pytest.approx(2.0e6)
    assert p[-1] == pytest.approx(5.0e5)
    assert np.allclose(p, np.linspace(2.0e6, 5.0e5, n), rtol=1e-9)


def test_poiseuille_bias_varying_film_is_monotone_with_exact_ends() -> None:
    """A converging gap keeps the ends exact and the field monotone high->low."""

    n, length = 401, 0.01
    _, step = _grid(n, length)
    film = 1.0e-5 * (1.0 + 0.5 * np.sin(np.linspace(0.0, np.pi, n)))
    p = poiseuille_bias_pressure(film, step, 1.0e6, 0.0)
    assert p[0] == pytest.approx(1.0e6)
    assert p[-1] == pytest.approx(0.0, abs=1e-6)
    assert np.all(np.diff(p) <= 1e-6)  # non-increasing from high to low end


def test_poiseuille_bias_equal_ends_is_uniform() -> None:
    """Equal end pressures give a uniform field (no throughflow gradient)."""

    n, length = 101, 0.01
    _, step = _grid(n, length)
    film = 1.0e-5 * (1.0 + 0.3 * np.cos(np.linspace(0.0, np.pi, n)))
    p = poiseuille_bias_pressure(film, step, 3.0e6, 3.0e6)
    assert np.allclose(p, 3.0e6, rtol=1e-12)
