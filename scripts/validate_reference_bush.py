"""Add the roller/blade journal freedom to the reference-machine harness.

The previous harness held the blade rigidly on its prescribed path, so the only thing that
could move the bush-blade gap was the bush itself -- and the modelled gap swung 5.8 um against
a measured 9.1 um. In the machine this project models, the equivalent freedom is already
carried: the rotor rides its crank pin on an oil film and its bore centre wanders by ~9 um
over a revolution (``e_j`` in rotor_bush_dynamics). Here the blade is integral with the
roller, so the same journal excursion translates the blade directly against the bush, and its
measured span (8.98 um normal to the pad in our machine) is the size of the amplitude the
harness is missing.

So the roller gets two freedoms of its own, driven exactly as the rotor's are: gas load on the
blade, the crank-pin journal film, the bush reaction, and the centrifugal term.

    m_r e''  =  F_gas(blade)  +  F_journal(e, e')  -  F_flat  +  m_r w^2 r_pin

The flat film gap then follows both bodies: h = (bush offset) - (blade offset along the pad
normal). Unknowns rise to four with the journal clearance, still against 144 measured points
and with the waveform phase to reproduce.
"""

import sys

sys.path.insert(0, "src")

import numpy as np
from scipy.integrate import solve_ivp

from mochi.arc_film import arc_film_force
from mochi.ocvirk_bearing import short_bearing_force
from mochi.reference_swing import SwingReference
from mochi.slider_film import flat_slider_film

REF = SwingReference()
GAUGE = REF.suction_pressure_pa - REF.discharge_pressure_pa
W = REF.angular_speed_rad_s
ETA = 2.8e-3

_R, _H = REF.bush_radius_m, 4.0e-3
_seg = _R**2 * np.arccos((_R - _H) / _R) - (_R - _H) * np.sqrt(2 * _R * _H - _H**2)
M_B = 7850.0 * _seg * REF.axial_length_m
I_B = 0.5 * M_B * _R**2
# Roller: annulus, outer R_o, bore taken as the crank eccentric radius, axial 25 mm.
R_BORE = 0.45 * REF.roller_radius_m
M_R = 7850.0 * np.pi * (REF.roller_radius_m**2 - R_BORE**2) * REF.axial_length_m
J_LEN = REF.axial_length_m

meas = np.load("refs/fig8_grid.npy")
TH_M, HU_M, HL_M = meas


def journal_force(ex, ey, vx, vy, c_j, omega_ent):
    """Ocvirk reaction on the roller, same mapping rotor_dynamics uses."""
    mag = np.hypot(ex, ey)
    if mag <= 1e-12:
        return 0.0, 0.0
    eps = min(mag / c_j, 0.995)
    rx, ry = ex / mag, ey / mag
    tx, ty = -ry, rx
    eps_dot = (ex * vx + ey * vy) / (mag * c_j)
    psi_dot = (ex * vy - ey * vx) / (mag * mag)
    f = short_bearing_force(eps, omega_ent - psi_dot, eccentricity_rate_per_s=eps_dot,
                            viscosity_pa_s=ETA, radius_m=R_BORE, length_m=J_LEN,
                            clearance_m=c_j)
    return -(f.radial_force_n * rx + f.tangential_force_n * tx), \
           -(f.radial_force_n * ry + f.tangential_force_n * ty)


def run(c_flat, c_seat, pad, c_j, revs=5, samples=144):
    def rhs(t, y):
        xb, gam, ex, ey, vxb, vgam, vex, vey = y
        th = W * t
        psi = REF.swing_angle_rad(th)
        # Pad normal in the world; the blade swings with psi, so the journal excursion
        # projects onto it through the swing angle.
        nx, ny = np.cos(psi), np.sin(psi)
        blade = ex * nx + ey * ny
        vblade = vex * nx + vey * ny
        h0 = float(np.clip(0.5 * c_flat + xb - blade, 0.03 * c_flat, 0.97 * c_flat))
        gmax = 0.9 * min(h0, c_flat - h0) / (0.5 * pad)
        g = float(np.clip(gam, -gmax, gmax))
        flat = flat_slider_film(
            c_flat - h0, g, -(vxb - vblade), vgam, -REF.slide_speed_m_s(th),
            length_m=pad, height_m=REF.axial_length_m, clearance_m=c_flat,
            viscosity_pa_s=ETA, n_s=161,
            pressure_start_pa=GAUGE, pressure_end_pa=0.0, cavitation_pressure_pa=GAUGE,
        )
        seat = arc_film_force(
            float(np.clip(-xb, -0.97 * c_seat, 0.97 * c_seat)), 0.0, -vxb, 0.0,
            abs(REF.swing_rate_rad_s(th)),
            arc_center_rad=0.0, arc_half_span_rad=0.5 * np.pi,
            radius_m=_R, length_m=REF.axial_length_m, clearance_m=c_seat,
            viscosity_pa_s=ETA, n_beta=101,
        )
        fgas, _ = REF.blade_gas_load(th)
        jx, jy = journal_force(ex, ey, vex, vey, c_j, 0.5 * W)
        pin = REF.eccentricity_m
        # Gas pushes the blade onto the suction bush (-n); the film pushes back (+n).
        fx = fgas * (-nx) + flat.normal_force_n * nx + jx + M_R * W * W * pin * np.sin(th)
        fy = fgas * (-ny) + flat.normal_force_n * ny + jy + M_R * W * W * pin * np.cos(th)
        return [vxb, vgam, vex, vey,
                (flat.normal_force_n - seat.force_x_n) / M_B,
                flat.moment_nm / I_B,
                fx / M_R, fy / M_R]

    period = 1.0 / REF.frequency_hz
    t_end = revs * period
    t_eval = np.linspace(t_end - period, t_end, samples)
    try:
        sol = solve_ivp(rhs, (0.0, t_end), np.zeros(8), method="BDF", t_eval=t_eval,
                        rtol=1e-6, atol=[1e-11, 1e-9, 1e-11, 1e-11, 1e-6, 1e-4, 1e-6, 1e-6])
    except Exception:
        return None
    if not sol.success or not np.all(np.isfinite(sol.y)):
        return None
    xb, gam, ex, ey = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    psi = np.array([REF.swing_angle_rad(float(W * t)) for t in t_eval])
    blade = ex * np.cos(psi) + ey * np.sin(psi)
    h0 = np.clip(0.5 * c_flat + xb - blade, 0.03 * c_flat, 0.97 * c_flat)
    gmax = 0.9 * np.minimum(h0, c_flat - h0) / (0.5 * pad)
    g = np.clip(gam, -gmax, gmax)
    s = 0.5 * REF.sensor_spacing_m
    ang = np.degrees(W * t_eval) % 360.0
    o = np.argsort(ang)
    hu = np.interp(TH_M, ang[o], ((h0 + s * g) * 1e6)[o], period=360.0)
    hl = np.interp(TH_M, ang[o], ((h0 - s * g) * 1e6)[o], period=360.0)
    return hu, hl, np.ptp(np.hypot(ex, ey)) * 1e6


print(f"roller m={M_R*1e3:.1f} g, bore R={R_BORE*1e3:.1f} mm")
print(f"measured: h_u ptp {np.ptp(HU_M):.2f} um, h_l ptp {np.ptp(HL_M):.2f} um")
print()
print(f"{'C_flat':>7} {'C_seat':>7} {'L_c':>5} {'c_j':>5} | {'RMS':>6} {'h_u ptp':>8} "
      f"{'|e_j| ptp':>10} {'corr_u':>7}")
print("-" * 70)
rows = []
for c_flat in (20e-6, 30e-6):
    for c_seat in (20e-6, 40e-6):
        for pad in (6.0e-3, 8.0e-3):
            for c_j in (10e-6, 20e-6):
                out = run(c_flat, c_seat, pad, c_j)
                if out is None:
                    continue
                hu, hl, ej = out
                ru = np.sqrt(np.mean((hu - HU_M) ** 2))
                rl = np.sqrt(np.mean((hl - HL_M) ** 2))
                rms = float(np.sqrt(0.5 * (ru**2 + rl**2)))
                cu = float(np.corrcoef(hu, HU_M)[0, 1]) if np.ptp(hu) > 1e-6 else 0.0
                rows.append((rms, c_flat, c_seat, pad, c_j, hu, hl, cu))
                print(f"{c_flat*1e6:7.0f} {c_seat*1e6:7.0f} {pad*1e3:5.1f} {c_j*1e6:5.0f} "
                      f"| {rms:6.2f} {np.ptp(hu):8.2f} {ej:10.2f} {cu:7.3f}", flush=True)
rows.sort(key=lambda r: r[0])
if rows:
    rms, cf, cs, pd_, cj, hu, hl, cu = rows[0]
    print()
    print(f"best: C_flat={cf*1e6:.0f} C_seat={cs*1e6:.0f} L_c={pd_*1e3:.1f}mm c_j={cj*1e6:.0f}um"
          f" -> RMS {rms:.2f} um, corr {cu:.3f}")
    np.save("refs/fit_roller.npy", np.vstack([TH_M, hu, hl]))
