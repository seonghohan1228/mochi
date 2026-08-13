"""Flat slider oil-film force (swing-bush piece flat against the fixed vane)."""

import math

import pytest

from mochi.bush_film import AXIAL_CLEARANCE_M, LUBRICANT_VISCOSITY_PA_S
from mochi.kinematics import RotaryGeometry
from mochi.slider_film import flat_slider_film

H = RotaryGeometry.default().cylinder_height_m - 2 * AXIAL_CLEARANCE_M
LF = 11.94e-3
CF = 10.0e-6


def _slider(approach, tilt, approach_rate, tilt_rate, slide):
    return flat_slider_film(
        approach,
        tilt,
        approach_rate,
        tilt_rate,
        slide,
        length_m=LF,
        height_m=H,
        clearance_m=CF,
    )


def test_parallel_non_approaching_pad_carries_no_load() -> None:
    """A flat parallel pad that only slides generates no hydrodynamic pressure."""

    f = _slider(0.0, 0.0, 0.0, 0.0, 0.86)
    assert f.normal_force_n == pytest.approx(0.0, abs=1e-9)
    assert f.moment_nm == pytest.approx(0.0, abs=1e-12)


def test_matches_the_incline_slider_closed_form() -> None:
    """A pure converging wedge matches Reynolds' analytic fixed-incline slider load.

    For a pad of length ``B`` with inlet/outlet films ``h1``/``h0`` (ratio ``a``) and
    surface speed ``U``, the load per unit width is
    ``w = 6 mu U B^2 / ((a-1)^2 h0^2) * (ln a - 2(a-1)/(a+1))``; the numerical film,
    times the height ``H``, must reproduce it.
    """

    mu, u, gamma = LUBRICANT_VISCOSITY_PA_S, 0.86, -5.0e-4
    f = flat_slider_film(0.0, gamma, 0.0, 0.0, u, length_m=LF, height_m=H, clearance_m=CF, n_s=4001)
    h1 = CF + 0.5 * LF * abs(gamma)  # inlet (thick, s = -L/2)
    h0 = CF - 0.5 * LF * abs(gamma)  # outlet (thin, s = +L/2)
    a = h1 / h0
    w = (
        (6.0 * mu * u * LF**2)
        / ((a - 1.0) ** 2 * h0**2)
        * (math.log(a) - 2.0 * (a - 1.0) / (a + 1.0))
    )
    assert f.normal_force_n == pytest.approx(w * H, rel=2.0e-3)


def test_squeeze_gives_restoring_normal_force() -> None:
    """Approaching the vane (approach_rate > 0) is resisted by a positive normal force."""

    f = _slider(0.0, 0.0, 0.02, 0.0, 0.0)
    assert f.normal_force_n > 0.0


@pytest.mark.parametrize("rate", [0.005, 0.02, 0.05])
def test_matches_the_parallel_squeeze_closed_form(rate) -> None:
    """A parallel-plate squeeze pad matches the analytic load ``W = mu ddot L^3 H / h^3``.

    With no tilt and no sliding the gap is uniform (``h = c_f``) and the Reynolds source
    is the constant squeeze ``-ddot``; solving ``h^3 p'' = -12 mu ddot`` with ``p = 0``
    at both ends gives ``p(s) = 6 mu ddot (L^2/4 - s^2)/h^3`` and, integrated over the
    pad (times the height ``H``), the closed-form normal force below.
    """

    mu = LUBRICANT_VISCOSITY_PA_S
    f = flat_slider_film(
        0.0, 0.0, rate, 0.0, 0.0, length_m=LF, height_m=H, clearance_m=CF, n_s=6001
    )
    w_closed = mu * rate * LF**3 * H / CF**3
    assert f.normal_force_n == pytest.approx(w_closed, rel=1e-4)
    assert f.moment_nm == pytest.approx(0.0, abs=1e-9)  # symmetric squeeze -> no moment


def test_converging_wedge_carries_load_and_moment() -> None:
    """A converging wedge with sliding builds load and a restoring moment."""

    f = _slider(0.0, -5.0e-4, 0.0, 0.0, 0.86)
    assert f.normal_force_n > 0.0
    assert f.moment_nm != pytest.approx(0.0, abs=1e-6)


def test_diverging_wedge_cavitates() -> None:
    """A diverging wedge generates no positive pressure (Gumbel)."""

    f = _slider(0.0, +5.0e-4, 0.0, 0.0, 0.86)
    assert f.normal_force_n == pytest.approx(0.0, abs=1e-9)


def test_load_grows_with_wedge_steepness() -> None:
    steep = _slider(0.0, -5.0e-4, 0.0, 0.0, 0.86)
    shallow = _slider(0.0, -2.0e-4, 0.0, 0.0, 0.86)
    assert steep.normal_force_n > shallow.normal_force_n


def test_contact_raises() -> None:
    with pytest.raises(ValueError):
        _slider(CF, 0.0, 0.0, 0.0, 0.0)  # approach == clearance -> zero film


def test_invalid_inputs_raise() -> None:
    with pytest.raises(ValueError):
        flat_slider_film(0.0, 0.0, 0.0, 0.0, 0.0, length_m=0.0, height_m=H, clearance_m=CF)
    with pytest.raises(ValueError):
        flat_slider_film(0.0, 0.0, 0.0, 0.0, 0.0, length_m=LF, height_m=H, clearance_m=CF, n_s=4)


# ---------------------------------------------------------------------------
# Gas-pressure boundary conditions (chamber/crank pressures at the pad ends)
# ---------------------------------------------------------------------------


def test_gas_bias_equal_ends_add_pressure_times_area() -> None:
    """Both ends at chamber pressure P add a uniform bias P over the pad (no leakage)."""

    p_gas = 2.0e6
    base = flat_slider_film(0.0, 0.0, 0.0, 0.0, 0.0, length_m=LF, height_m=H, clearance_m=CF)
    biased = flat_slider_film(
        0.0, 0.0, 0.0, 0.0, 0.0, length_m=LF, height_m=H, clearance_m=CF,
        pressure_start_pa=p_gas, pressure_end_pa=p_gas,
    )
    assert biased.normal_force_n - base.normal_force_n == pytest.approx(p_gas * LF * H, rel=1e-3)
    assert biased.throughflow_m3_s == 0.0


def test_gas_bias_ramp_shifts_load_and_leaks() -> None:
    """A P->0 pressure drop across the pad adds the mean (P/2) load and a Poiseuille leak."""

    p_gas, mu = 1.0e6, LUBRICANT_VISCOSITY_PA_S
    f = flat_slider_film(
        0.0, 0.0, 0.0, 0.0, 0.0, length_m=LF, height_m=H, clearance_m=CF, n_s=4001,
        pressure_start_pa=p_gas, pressure_end_pa=0.0,
    )
    assert f.normal_force_n == pytest.approx(0.5 * p_gas * LF * H, rel=2e-3)
    assert f.throughflow_m3_s == pytest.approx(H * p_gas * CF**3 / (12.0 * mu * LF), rel=1e-3)


def test_gas_bias_defaults_are_backward_compatible() -> None:
    """Zero end pressures (the default) reproduce the pure hydrodynamic film exactly."""

    a = flat_slider_film(0.0, 1.0e-4, 0.0, 0.0, 0.5, length_m=LF, height_m=H, clearance_m=CF)
    b = flat_slider_film(
        0.0, 1.0e-4, 0.0, 0.0, 0.5, length_m=LF, height_m=H, clearance_m=CF,
        pressure_start_pa=0.0, pressure_end_pa=0.0,
    )
    assert a.normal_force_n == pytest.approx(b.normal_force_n)
    assert a.throughflow_m3_s == 0.0
