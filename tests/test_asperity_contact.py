"""Tests for the Greenwood-Tripp asperity contact (:mod:`mochi.asperity_contact`)."""

import numpy as np
import pytest

from mochi.asperity_contact import AsperityParams, greenwood_tripp_pressure


def test_roughness_parameter_in_literature_range():
    # eta*beta*sigma ~ 0.02-0.06 for ground steel.
    assert 0.02 <= AsperityParams().roughness_parameter <= 0.06


def test_zero_beyond_four_sigma():
    p = AsperityParams()
    # No asperity contact once the surfaces are > 4 sigma apart.
    assert float(greenwood_tripp_pressure(np.array([4.0 * p.roughness_m]), p)[0]) == 0.0
    assert float(greenwood_tripp_pressure(np.array([6.0 * p.roughness_m]), p)[0]) == 0.0


def test_monotonic_decreasing_with_separation():
    p = AsperityParams()
    h = np.linspace(-1.0e-6, 3.0e-6, 40)
    pres = greenwood_tripp_pressure(h, p)
    # Pressure must not increase as the gap opens.
    assert np.all(np.diff(pres) <= 1e-6)
    # Steep rise near contact: p(sigma) >> p(2 sigma).
    p_1s = float(greenwood_tripp_pressure(np.array([p.roughness_m]), p)[0])
    p_2s = float(greenwood_tripp_pressure(np.array([2.0 * p.roughness_m]), p)[0])
    assert p_1s > 5.0 * p_2s > 0.0


def test_deep_contact_capped():
    p = AsperityParams()
    # A deep incursion is capped at the plastic-flow scale (no overflow).
    deep = float(greenwood_tripp_pressure(np.array([-5.0e-6]), p)[0])
    assert deep == pytest.approx(5.0e9)


def test_scales_with_reduced_modulus():
    soft = AsperityParams(reduced_modulus_pa=50e9)
    stiff = AsperityParams(reduced_modulus_pa=150e9)
    h = np.array([0.8e-6])
    assert greenwood_tripp_pressure(h, stiff)[0] == pytest.approx(
        3.0 * greenwood_tripp_pressure(h, soft)[0], rel=1e-6
    )


def test_validate_rejects_nonpositive():
    with pytest.raises(ValueError):
        AsperityParams(summit_radius_m=0.0).validate()
