"""Flat slider oil-film force (swing-bush piece flat against the fixed vane), 1-D.

Each bush piece's flat face runs on the fixed vane flank -- a flat pad on a flat
wall (Section 3.6/4.11). The pressure is solved with the **axial-uniform
(long-bearing) 1-D Reynolds** reduction: at the pad aspect ratio the sliding length
``L_f`` is short and the axial height ``H`` long (``L/D > 1``), so the film escapes
**along the sliding direction**, not axially -- ``dp/dz ~ 0`` and the pressure is
solved along ``s`` and multiplied by ``H``. This is the correct reduction here (the
opposite of the axial short-bearing used for the ``L/D < 1`` journal).

At each position ``s`` along the flat (``s = 0`` at the pad centre, span ``L_f``,
axial height ``H``) the film is

    h(s) = (c_f - delta) + s * gamma,

with ``delta`` the pad's normal approach toward the vane and ``gamma`` its in-plane
tilt (the wedge = the piece attitude ``phi_k``). Along ``s`` the Reynolds equation
``d/ds (h^3/(12 mu) dp/ds) = ubar dh/ds + dh/dt`` becomes

    d/ds ( h^3 dp/ds ) = 12 mu S(s),
    S(s) = (U/2) gamma - delta_dot + s * gamma_dot,

``U`` the sliding speed along the vane (bush translation, Section 3.6) and
``ubar = U/2`` the mean surface speed. The pad is open to oil at both ends, so
``p = 0`` there; the line problem is solved by
:func:`mochi.line_reynolds.solve_line_pressure` with Gumbel cavitation. Integrating
the pressure over the pad (times ``H``) gives the **normal force** (pushing the piece
off the vane) and the **moment** about the pad centre (the attitude-restoring moment
for ``phi_k``):

    F_n = H integral p ds,   M = H integral p * s ds.

A parallel, non-approaching pad carries no hydrodynamic load (only Couette shear);
load comes from the converging wedge (``gamma`` with sliding) or squeeze
(``delta_dot``). See PHYSICS.md Section 4.11.
"""

from dataclasses import dataclass
from math import isfinite

import numpy as np

from mochi.bush_film import LUBRICANT_VISCOSITY_PA_S
from mochi.line_reynolds import poiseuille_bias_pressure, solve_line_pressure


@dataclass(frozen=True, slots=True)
class SliderFilmForce:
    """Flat slider film reaction on a bush piece (SI units).

    ``normal_force_n`` pushes the piece away from the vane (restoring against the
    approach ``delta``); ``moment_nm`` is about the pad centre (restoring ``phi_k``).
    """

    normal_force_n: float
    moment_nm: float
    min_film_thickness_m: float
    throughflow_m3_s: float = 0.0  # pressure-driven (Poiseuille) leakage along the pad


def flat_slider_film(
    approach_m: float,
    tilt_rad: float,
    approach_rate_m_s: float,
    tilt_rate_rad_s: float,
    slide_speed_m_s: float,
    *,
    length_m: float,
    height_m: float,
    clearance_m: float,
    viscosity_pa_s: float = LUBRICANT_VISCOSITY_PA_S,
    n_s: int = 201,
    pressure_start_pa: float = 0.0,
    pressure_end_pa: float = 0.0,
    cavitation_pressure_pa: float = 0.0,
) -> SliderFilmForce:
    """Axial-uniform 1-D flat-slider film force and moment on a bush piece.

    ``approach_m`` (delta) is the pad's motion toward the vane (gap = ``c_f - delta``),
    ``tilt_rad`` (gamma) its attitude; their rates are ``approach_rate`` /
    ``tilt_rate``, and ``slide_speed`` is the sliding along the vane. ``length_m`` is
    the sliding span ``L_f`` and ``height_m`` the axial height ``H``. Returns the
    normal force (off the vane) and the moment about the pad centre.

    ``pressure_start_pa`` / ``pressure_end_pa`` are the gas pressures the pad ends
    (``s = -L_f/2`` and ``s = +L_f/2``) open to -- the local chamber and crank/bore
    pressures (gauge). They superpose the Poiseuille gas-bias field onto the
    hydrodynamic pressure (Reynolds is linear in ``p``; both default to zero). The
    pressure-driven leakage is reported in ``throughflow_m3_s``.

    ``cavitation_pressure_pa`` is the pressure below which the film cavitates, in the **same
    gauge** as ``pressure_start_pa`` / ``pressure_end_pa``; the clamp is applied to the total
    (hydrodynamic + gas) field. The default ``0.0`` reproduces classic Gumbel
    (half-Sommerfeld) cavitation for the pure-hydrodynamic pad. When the ends are referenced
    to a pressurised ambient (e.g. the rotor-bore gas), pass the ambient's negative
    (``-p_ambient``) so the clamp sits at absolute zero rather than at the reference.
    """

    for value in (approach_m, tilt_rad, approach_rate_m_s, tilt_rate_rad_s, slide_speed_m_s):
        if not isfinite(value):
            raise ValueError("Pad state and its rates must be finite.")
    for name, value in (("length", length_m), ("height", height_m), ("clearance", clearance_m)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"Pad {name} must be a positive, finite length in metres.")
    if not isfinite(viscosity_pa_s) or viscosity_pa_s <= 0.0:
        raise ValueError("Viscosity must be a positive, finite value in Pa*s.")
    if n_s < 16:
        raise ValueError("The slider grid needs at least 16 points.")

    s = np.linspace(-0.5 * length_m, 0.5 * length_m, n_s)
    film = (clearance_m - approach_m) + s * tilt_rad
    if np.any(film <= 0.0):
        raise ValueError("Pad approach/tilt closes the film (metal contact).")

    source = 0.5 * slide_speed_m_s * tilt_rad - approach_rate_m_s + s * tilt_rate_rad_s
    line_source = 12.0 * viscosity_pa_s * source

    d_s = length_m / (n_s - 1)
    gas_present = pressure_start_pa != 0.0 or pressure_end_pa != 0.0
    # Reynolds is linear in p, so the hydrodynamic field is solved UNCLAMPED and the gas
    # bias superposed; cavitation is then applied to the **total**, which is where the
    # physical constraint lives (a sub-cavitation total cannot exist, whatever its parts
    # are). Clamping only the hydrodynamic part -- or skipping the clamp because gas is
    # present -- lets a reciprocating pad keep a huge unphysical suction lobe on the
    # diverging half-stroke and slams it into contact.
    pressure = solve_line_pressure(film, line_source, d_s, cavitation=False)

    throughflow = 0.0
    if gas_present:
        gas = poiseuille_bias_pressure(film, d_s, pressure_start_pa, pressure_end_pa)
        pressure = pressure + gas  # superposition (Reynolds linear in p)
        inv_h3_integral = float(np.trapezoid(film**-3.0, dx=d_s))
        throughflow = (
            height_m * (pressure_start_pa - pressure_end_pa)
            / (12.0 * viscosity_pa_s * inv_h3_integral)
        )
    np.maximum(pressure, cavitation_pressure_pa, out=pressure)

    normal_force = height_m * float(np.trapezoid(pressure, dx=d_s))
    moment = height_m * float(np.trapezoid(pressure * s, dx=d_s))
    return SliderFilmForce(
        normal_force_n=normal_force,
        moment_nm=moment,
        min_film_thickness_m=float(np.min(film)),
        throughflow_m3_s=throughflow,
    )
