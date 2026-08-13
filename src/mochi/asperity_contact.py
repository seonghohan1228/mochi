"""Greenwood-Tripp statistical asperity contact for mixed-lubrication bush films.

When a bush film (the curved rotor-groove arc or the flat vane pad) is squeezed into the
surface roughness (film parameter Lambda = h_min/sigma < ~3), the asperities carry part of
the load. This module gives that **asperity contact pressure** as a function of the local
mean separation ``h``, grounded in measurable surface parameters rather than an arbitrary
pressure scale.

Greenwood-Tripp (elastic, Gaussian summit heights):

    p_asp(h) = K * E' * sqrt(sigma / beta) * F_{5/2}(h / sigma)
    K        = (16*sqrt(2)/15) * pi * (eta * beta * sigma)^2

with ``sigma`` the composite (summit) roughness, ``beta`` the summit radius, ``eta`` the
areal summit density, ``E'`` the reduced elastic modulus, and ``F_{5/2}`` the standardized
Gaussian summit integral. ``F_{5/2}`` uses the widely-used engineering fit

    F_{5/2}(H) ~ 4.4086e-5 * (4 - H)^6.804   (H < 4),   0 otherwise,

which extends naturally to ``H < 0`` (deep contact). The result feeds the hydrodynamic
Reynolds film by **load sharing** (the total film + asperity force supports the load),
replacing the hard contact clamps. See PHYSICS.md Section 4.11 / docs/mixed_lubrication.

The dimensionless roughness group ``eta*beta*sigma`` is ~0.02-0.06 for ground steel; the
individual ``beta``/``eta`` below are literature-typical placeholders to be refined by a
surface measurement (the friction is sensitive to them -- run a sensitivity sweep).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np

# Literature-typical ground-steel summit parameters (refine by measurement).
DEFAULT_SUMMIT_RADIUS_M = 10.0e-6  # beta
DEFAULT_SUMMIT_DENSITY_PER_M2 = 1.0e10  # eta  (eta*beta*sigma ~ 0.05 at sigma=0.53 um)
DEFAULT_REDUCED_MODULUS_PA = 115.0e9  # E' for steel-steel: 1/E' = 2(1-nu^2)/E
# Plastic-flow cap on the asperity pressure (~steel hardness scale) so a transient deep
# incursion during the stiff ODE's Newton iterations cannot overflow.
_MAX_ASPERITY_PRESSURE_PA = 5.0e9


@dataclass(frozen=True, slots=True)
class AsperityParams:
    """Surface/material parameters for the Greenwood-Tripp asperity contact (SI units)."""

    roughness_m: float = 0.53e-6  # composite summit roughness sigma
    summit_radius_m: float = DEFAULT_SUMMIT_RADIUS_M  # beta
    summit_density_per_m2: float = DEFAULT_SUMMIT_DENSITY_PER_M2  # eta
    reduced_modulus_pa: float = DEFAULT_REDUCED_MODULUS_PA  # E'

    def validate(self) -> None:
        for name in (
            "roughness_m",
            "summit_radius_m",
            "summit_density_per_m2",
            "reduced_modulus_pa",
        ):
            value = getattr(self, name)
            if not (np.isfinite(value) and value > 0.0):
                raise ValueError(f"AsperityParams.{name} must be a positive, finite number.")

    @property
    def roughness_parameter(self) -> float:
        """The dimensionless roughness group ``eta*beta*sigma`` (~0.02-0.06 for ground steel)."""
        return self.summit_density_per_m2 * self.summit_radius_m * self.roughness_m


#: Module-level default so the parameters are constructed once, not per call.
DEFAULT_ASPERITY = AsperityParams()


def greenwood_tripp_pressure(h_m, params: AsperityParams | None = None):
    """Asperity contact pressure ``p_asp(h)`` (Pa), vectorised over the separation ``h_m``.

    ``h_m`` is the local mean film/separation (can be <= 0 for deep contact). Returns 0 where
    the surfaces are more than ~4 sigma apart (no asperity contact), and is capped at a
    plastic-flow pressure for numerical safety in deep incursions.
    """

    params = DEFAULT_ASPERITY if params is None else params
    sigma = params.roughness_m
    lam = np.asarray(h_m, dtype=float) / sigma  # H = h/sigma
    base = np.maximum(4.0 - lam, 0.0)
    f_five_halves = 4.4086e-5 * base**6.804  # standardized Gaussian summit integral (fit)
    k = (16.0 * sqrt(2.0) / 15.0) * pi * params.roughness_parameter**2
    pressure = k * params.reduced_modulus_pa * sqrt(sigma / params.summit_radius_m) * f_five_halves
    return np.minimum(pressure, _MAX_ASPERITY_PRESSURE_PA)
