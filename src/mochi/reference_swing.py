"""Kinematics and blade gas load of the *reference* swing compressor (Tanaka rig).

This is the machine of the single-case film-thickness validation (see
``single_case_validation_plan.md``): a **Daikin-type swing** compressor, where the blade is
integral with the roller and the bush pivots in a seat in the cylinder. That is the mirror of
the machine this project models (fixed vane, bush riding in the rotor groove), so this module
deliberately lives apart from :mod:`mochi.kinematics` -- it exists to drive the *bush film
module* with the reference machine's motion and load, not to model our compressor.

Geometry from Tanaka et al. (2002) Purdue C10-4 Table 3 (rig) and the bush drawing of
Tanaka et al. (2008) Trans. JSRAE 25(4) Fig.2; operating conditions from the 2008 standard
case. Every value that is *not* in either paper is flagged in :class:`SwingReference` and is
meant to be swept, not trusted.

**Kinematics.** With the cylinder centre at the origin, the bush pivot fixed at distance
``a`` and the roller centre orbiting at the eccentricity ``e``,

    rho(theta)   = sqrt(a^2 + e^2 - 2 a e cos theta)      (roller centre -> bush centre)
    psi(theta)   = atan2(-e sin theta, a - e cos theta)   (blade swing angle)
    ell(theta)   = rho - R_b - R_o                        (blade length exposed to the gas)

``ell`` runs from 0 at top dead centre to ``2e`` at bottom dead centre -- the blade stroke --
which is the check :func:`SwingReference.validate` makes. The blade slides through the bush
slot at ``U = d rho/dt`` and the bush swings with it at ``psi_dot``, so both the entrainment
and the squeeze at the bush-blade film follow from these.

**Gas load.** The chamber volumes come from the same true-geometry sweep this project uses
elsewhere: at azimuth ``phi`` from the cylinder centre the roller surface sits at
``r(phi) = e cos(phi - phi_r) + sqrt(R_o^2 - e^2 sin^2(phi - phi_r))``, so the chamber area is
``(1/2) int (R_c^2 - r^2) dphi``. The compression chamber is then polytropic from suction
until it reaches the discharge pressure. The blade carries the chamber pressure difference
over its exposed length, giving the force and the moment about the bush pivot that the bush
films must react.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, pi, sin, sqrt

import numpy as np

MM = 1.0e-3


@dataclass(frozen=True)
class SwingReference:
    """Reference-machine geometry and operating point (SI units).

    Sourced values carry their citation; the rest are estimates to sweep.
    """

    # --- geometry, Tanaka 2002 Table 3 (rig) ---
    cylinder_radius_m: float = 22.0 * MM  # R_c
    roller_radius_m: float = 16.8 * MM  # R_o
    eccentricity_m: float = 5.2 * MM  # e = R_c - R_o
    bush_radius_m: float = 6.0 * MM  # R_b  (also J08 Fig.2 "R6")
    blade_thickness_m: float = 4.0 * MM  # blade width a
    axial_length_m: float = 25.0 * MM  # L_p (also J08 Fig.2 "25")

    # --- bush pad, J08 Fig.2 (dimensioned drawing) ---
    pad_length_m: float = 10.0 * MM  # bush flat height, along the sliding direction
    sensor_spacing_m: float = 4.5 * MM  # L in J08 Eq.(1); sensors at +/- L/2

    # --- operating point, J08 standard case ---
    suction_pressure_pa: float = 0.4e6
    discharge_pressure_pa: float = 1.6e6
    frequency_hz: float = 28.0
    polytropic_index: float = 1.18  # R410A, typical; sweep

    # --- ESTIMATED, not in either paper -- sweep these ---
    bush_pivot_distance_m: float = 28.0 * MM
    """Distance from the cylinder centre to the bush pivot, ``a``. **Estimated.** Neither
    paper gives it. Taken as ``R_c + R_b`` so the bush's inner surface is flush with the
    bore, which is the only placement that neither intrudes into the swept volume nor leaves
    a step. Sets the swing amplitude (``asin(e/a)`` ~ 10.7 deg), so sweep it."""

    @property
    def angular_speed_rad_s(self) -> float:
        return 2.0 * pi * self.frequency_hz

    @property
    def swing_amplitude_rad(self) -> float:
        """Peak blade/bush swing angle, ``asin(e/a)`` for ``a > e``."""
        return abs(np.arcsin(self.eccentricity_m / self.bush_pivot_distance_m))

    def validate(self) -> None:
        """Check the geometric identities the reference machine must satisfy."""
        if not abs((self.cylinder_radius_m - self.roller_radius_m) - self.eccentricity_m) < 1e-9:
            raise ValueError("Reference geometry: e must equal R_c - R_o.")
        if self.bush_pivot_distance_m <= self.eccentricity_m:
            raise ValueError("Bush pivot distance must exceed the eccentricity.")
        # The blade stroke must come out as 2e (TDC fully retracted, BDC fully extended).
        stroke = self.blade_exposure_m(pi) - self.blade_exposure_m(0.0)
        if abs(stroke - 2.0 * self.eccentricity_m) > 1e-9:
            raise ValueError(f"Blade stroke {stroke:.6e} m should be 2e.")
        if not (self.sensor_spacing_m < self.pad_length_m):
            raise ValueError("The sensors must sit inside the pad.")

    # ------------------------------------------------------------------
    # Kinematics
    # ------------------------------------------------------------------

    def roller_to_bush_m(self, theta: float) -> float:
        """``rho(theta)``: roller centre to bush pivot."""
        a, e = self.bush_pivot_distance_m, self.eccentricity_m
        return sqrt(a * a + e * e - 2.0 * a * e * cos(theta))

    def swing_angle_rad(self, theta: float) -> float:
        """Blade (and bush) swing angle about the bush pivot; zero at TDC/BDC."""
        a, e = self.bush_pivot_distance_m, self.eccentricity_m
        return atan2(-e * sin(theta), a - e * cos(theta))

    def blade_exposure_m(self, theta: float) -> float:
        """Blade length exposed to the gas: 0 at TDC, ``2e`` at BDC."""
        return self.roller_to_bush_m(theta) - self.bush_radius_m - self.roller_radius_m

    def slide_speed_m_s(self, theta: float) -> float:
        """``U = d rho/dt`` -- the blade sliding through the bush slot."""
        a, e = self.bush_pivot_distance_m, self.eccentricity_m
        return self.angular_speed_rad_s * a * e * sin(theta) / self.roller_to_bush_m(theta)

    def swing_rate_rad_s(self, theta: float) -> float:
        """``psi_dot`` -- the bush's own angular velocity in its seat."""
        a, e = self.bush_pivot_distance_m, self.eccentricity_m
        rho = self.roller_to_bush_m(theta)
        return self.angular_speed_rad_s * e * (e - a * cos(theta)) / (rho * rho)

    # ------------------------------------------------------------------
    # Chamber volumes and gas load
    # ------------------------------------------------------------------

    def chamber_area_m2(self, theta: float, n: int = 2001) -> float:
        """Crescent area swept between the blade and the contact point, at shaft angle theta.

        Integrates ``(1/2)(R_c^2 - r(phi)^2)`` over the azimuth from the blade to the roller
        contact, with ``r(phi)`` the roller surface seen from the cylinder centre.
        """
        rc, ro, e = self.cylinder_radius_m, self.roller_radius_m, self.eccentricity_m
        phi = np.linspace(0.0, theta, n) if theta > 0.0 else np.zeros(n)
        off = phi - theta  # azimuth relative to the roller-centre direction
        r = e * np.cos(off) + np.sqrt(np.maximum(ro * ro - e * e * np.sin(off) ** 2, 0.0))
        return 0.5 * float(np.trapezoid(rc * rc - r * r, x=phi)) if theta > 0.0 else 0.0

    def compression_pressure_pa(self, theta: float) -> float:
        """Polytropic compression from suction, capped at the discharge pressure.

        The compression chamber is the one closed off behind the blade; its volume shrinks
        from the full crescent at TDC to zero at TDC of the next turn.
        """
        full = self.chamber_area_m2(2.0 * pi)
        volume = self.chamber_area_m2(2.0 * pi - theta) if theta > 0.0 else full
        if volume <= 0.0:
            return self.discharge_pressure_pa
        ratio = full / volume
        p = self.suction_pressure_pa * ratio**self.polytropic_index
        return min(p, self.discharge_pressure_pa)

    def blade_gas_load(self, theta: float) -> tuple[float, float]:
        """``(force_n, moment_nm)`` the gas puts on the blade, about the bush pivot.

        The pressure difference across the blade acts over its exposed length ``ell`` and the
        axial length; its centroid sits ``R_b + ell/2`` from the bush pivot, which is the
        moment arm the bush films have to react.
        """
        ell = max(self.blade_exposure_m(theta), 0.0)
        d_p = self.compression_pressure_pa(theta) - self.suction_pressure_pa
        force = d_p * ell * self.axial_length_m
        arm = self.bush_radius_m + 0.5 * ell
        return force, force * arm
