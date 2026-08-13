"""Infinitely-long (Sommerfeld) journal-bearing analytical film force.

The **long-bearing** (axial-uniform) limit ``L/D -> inf``: the Poiseuille pressure
escapes only *circumferentially* (``dp/dz ~ 0``), so the pressure is a 1-D function of
the circumferential angle alone. This is the exact reduction the swing-bush **curved**
film uses (:mod:`mochi.arc_film`, Section 4.11), and the analytical counterpart -- at
the opposite aspect ratio -- to the short-bearing Ocvirk model
(:mod:`mochi.ocvirk_bearing`, Section 4.9). It is the **closed-form benchmark** for
``arc_film``.

Steady, pure entrainment (no squeeze). With the film measured from the maximum-film
point, ``h(theta) = c(1 + eps cos theta)`` (thickest at ``theta = 0``, thinnest at
``theta = pi``), the exact **Sommerfeld pressure** is

    p(theta) = (6 mu U R / c^2) * eps sin(theta) (2 + eps cos(theta))
                                / [(2 + eps^2)(1 + eps cos(theta))^2],

``U = 2 R Omega`` the mean surface speed for entrainment ``Omega`` (so ``ubar = U/2 =
R Omega``, matching ``arc_film``'s ``entrainment_speed_rad_s``). Two textbook closures:

* **full-Sommerfeld** -- keep ``p`` over the whole ``0..2pi``; by antisymmetry the load
  is purely tangential (attitude ``90 deg``, radial load zero). Physical only for a
  fully flooded, cavitation-free film.
* **half-Sommerfeld (Gumbel)** -- clamp ``p >= 0``; only the converging half
  ``0..pi`` carries load. The attitude angle is the pure-number invariant
  ``tan(phi) = pi sqrt(1 - eps^2) / (2 eps)`` (independent of ``mu, U, R, L, c``).

The load is the pressure resultant over the wetted arc (radius ``R``, axial length
``L``), reported in the natural frame: ``radial_force_n`` is the film force on the
journal along the maximum-film direction (opposing the eccentricity, restoring, the
same sense as ``arc_film``'s force on the piece); ``tangential_force_n`` is the
cross-load. Because ``p`` is written from the exact closed form and only *integrated*
numerically, this is an analytical reference independent of ``arc_film``'s
finite-difference solve. See PHYSICS.md Section 4.11.
"""

from dataclasses import dataclass
from math import atan2, hypot

import numpy as np

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S


@dataclass(frozen=True, slots=True)
class LongBearingLoad:
    """Infinitely-long journal-bearing film load in the natural frame (SI units)."""

    radial_force_n: float  # along the maximum-film direction (opposes eccentricity)
    tangential_force_n: float
    magnitude_n: float
    attitude_angle_rad: float  # from the line of centres to the load, in (0, pi/2]
    max_pressure_pa: float
    min_film_thickness_m: float


def sommerfeld_pressure(
    eccentricity_ratio: float,
    theta_rad,
    *,
    radius_m: float,
    clearance_m: float,
    entrainment_speed_rad_s: float,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
):
    """Exact full-Sommerfeld long-bearing pressure at ``theta`` (from the max-film point).

    Accepts a scalar or array ``theta_rad`` and returns the same shape. This is the
    fully flooded (non-cavitated) pressure; the half-Sommerfeld field is
    ``max(p, 0)``.
    """

    eps = eccentricity_ratio
    theta = np.asarray(theta_rad, dtype=float)
    surface_speed = 2.0 * radius_m * entrainment_speed_rad_s
    prefactor = 6.0 * viscosity_pa_s * surface_speed * radius_m / clearance_m**2
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    return (
        prefactor * eps * sin_t * (2.0 + eps * cos_t) / ((2.0 + eps**2) * (1.0 + eps * cos_t) ** 2)
    )


def long_bearing_load(
    eccentricity_ratio: float,
    entrainment_speed_rad_s: float,
    *,
    radius_m: float,
    length_m: float,
    clearance_m: float,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    condition: str = "half",
    n_theta: int = 4001,
) -> LongBearingLoad:
    """Closed-form long-bearing film load by integrating the exact Sommerfeld pressure.

    ``condition="half"`` is the Gumbel (half-Sommerfeld) converging half ``0..pi``;
    ``"full"`` is the fully flooded ``0..2pi`` (radial load vanishes, attitude
    ``90 deg``). The pressure acts radially inward over ``dA = R L d theta``, so the
    force on the journal is ``F = -R L integral p (cos theta, sin theta) d theta`` --
    the same convention as :func:`mochi.arc_film.arc_film_force`.
    """

    if not (0.0 <= eccentricity_ratio < 1.0):
        raise ValueError("Eccentricity ratio must be in [0, 1).")
    for name, value in (("radius", radius_m), ("length", length_m), ("clearance", clearance_m)):
        if not (value > 0.0):
            raise ValueError(f"Bearing {name} must be a positive length in metres.")
    if viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive value in Pa*s.")
    if condition not in ("half", "full"):
        raise ValueError("condition must be 'half' (Gumbel) or 'full' (Sommerfeld).")
    if n_theta < 16:
        raise ValueError("Need at least 16 circumferential samples.")

    upper = np.pi if condition == "half" else 2.0 * np.pi
    theta = np.linspace(0.0, upper, n_theta)
    pressure = sommerfeld_pressure(
        eccentricity_ratio,
        theta,
        radius_m=radius_m,
        clearance_m=clearance_m,
        entrainment_speed_rad_s=entrainment_speed_rad_s,
        viscosity_pa_s=viscosity_pa_s,
    )
    if condition == "half":
        pressure = np.maximum(pressure, 0.0)  # Gumbel cavitation (converging half only)

    scale = radius_m * length_m
    radial = -scale * float(np.trapezoid(pressure * np.cos(theta), theta))
    tangential = -scale * float(np.trapezoid(pressure * np.sin(theta), theta))
    return LongBearingLoad(
        radial_force_n=radial,
        tangential_force_n=tangential,
        magnitude_n=hypot(radial, tangential),
        attitude_angle_rad=atan2(abs(tangential), radial),
        max_pressure_pa=float(np.max(pressure)),
        min_film_thickness_m=clearance_m * (1.0 - eccentricity_ratio),
    )


def half_sommerfeld_attitude_rad(eccentricity_ratio: float) -> float:
    """Half-Sommerfeld long-bearing attitude angle ``atan(pi sqrt(1-eps^2)/(2 eps))``.

    The pure-number invariant that any correct long-bearing (Gumbel) solver must
    reproduce, independent of ``mu, U, R, L, c``.
    """

    eps = eccentricity_ratio
    if not (0.0 < eps < 1.0):
        raise ValueError("Eccentricity ratio must be in (0, 1).")
    from math import atan, sqrt

    return atan(np.pi * sqrt(1.0 - eps * eps) / (2.0 * eps))
