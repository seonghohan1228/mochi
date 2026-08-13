"""Tests for the internal physical-plausibility checks (:mod:`mochi.physical_checks`).

The first test is the point of the module: it reconstructs the 2026-08 cavitation bug and
asserts the check catches it. The rest guard against the checks firing on valid states.
"""

import numpy as np
import pytest

from mochi.physical_checks import (
    PhysicalityError,
    check_eccentricity_ratio,
    check_energy_sign,
    check_film_positive,
    check_pressure_field,
)
from mochi.slider_film import flat_slider_film

BORE_PA = 4.0e6
CHAMBER_PA = 0.82e6
PAD = dict(length_m=11.94e-3, height_m=20.98e-3, clearance_m=10.0e-6, n_s=401)


def test_catches_the_2026_08_cavitation_bug():
    """The state that ran at -11.6 MPa absolute must be rejected."""
    with pytest.raises(PhysicalityError, match="absolute pressure"):
        flat_slider_film(
            10.0e-6 - 0.906e-6, 1.115e-5, 0.0, 0.0, 0.86,
            pressure_start_pa=0.0, pressure_end_pa=CHAMBER_PA - BORE_PA,
            cavitation_pressure_pa=-1.0e12,  # the old "no clamp under gas" behaviour
            check_physical=True, reference_pressure_pa=BORE_PA, **PAD,
        )


def test_current_cavitation_floor_passes():
    """With the floor at the lowest connected pressure the same state is physical."""
    f = flat_slider_film(
        10.0e-6 - 0.906e-6, 1.115e-5, 0.0, 0.0, 0.86,
        pressure_start_pa=0.0, pressure_end_pa=CHAMBER_PA - BORE_PA,
        cavitation_pressure_pa=CHAMBER_PA - BORE_PA,
        check_physical=True, reference_pressure_pa=BORE_PA, **PAD,
    )
    assert np.isfinite(f.normal_force_n)


def test_pure_hydrodynamic_path_does_not_false_positive():
    f = flat_slider_film(2.0e-6, -1.0e-5, 0.0, 0.0, 0.86, check_physical=True, **PAD)
    assert f.normal_force_n > 0.0


def test_pressure_field_respects_the_reference():
    # -1 MPa gauge is fine when immersed at 4 MPa, and impossible referenced to absolute.
    field = np.array([-1.0e6, 0.0, 1.0e6])
    check_pressure_field(field, reference_pressure_pa=BORE_PA)
    with pytest.raises(PhysicalityError):
        check_pressure_field(field, reference_pressure_pa=0.0)


def test_film_positive_and_eccentricity_and_energy():
    check_film_positive(np.array([1e-6, 2e-6]))
    with pytest.raises(PhysicalityError, match="interpenetrate"):
        check_film_positive(np.array([1e-6, -1e-9]))

    check_eccentricity_ratio(0.98)
    with pytest.raises(PhysicalityError, match="eccentricity ratio"):
        check_eccentricity_ratio(1.33)  # the value the coupled orbit actually reaches

    check_energy_sign(6.15)
    with pytest.raises(PhysicalityError, match="cannot be negative"):
        check_energy_sign(-0.1)
