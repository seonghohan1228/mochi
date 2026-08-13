# Stage 5 spec — 9-DOF rotor + two-piece swing-bush multibody dynamics

**Status:** spec locked, no code yet. Build this in a fresh, full-context session by
reading this file top to bottom. New module: `src/mochi/rotor_bush_dynamics.py`
(+ `tests/test_rotor_bush_dynamics.py`, PHYSICS §4.14). This extends the rotor
lateral EOM (`rotor_dynamics.py`, §4.13) by making the **two swing-bush pieces
independent bodies** and coupling them to the rotor through the curved (arc) and flat
(slider) oil films rebuilt in the axial-uniform 1-D form (`arc_film.py`,
`slider_film.py`, `line_reynolds.py`, §4.11).

Shaft axis == cylinder axis throughout (the crank pin orbits; the shaft centreline
does not). All films are Gumbel-cavitated, POE VG68, µ = 0.010 Pa·s.

---

## 1. Confirmed data (all SI; from CAD / user 2026-07)

| Quantity | Symbol | Value | Source |
|---|---|---|---|
| Crank throw (pin orbit radius) | e_throw | 4.5 mm | `geometry.eccentricity_m` |
| Drive | ω, f | 188.5 rad/s, 30 Hz | `geometry.angular_speed_rad_s` |
| Rotor mass | m_r | **0.275 kg** | `ROTOR_MASS_KG` (confirmed) |
| Rotor polar inertia about **CG** | I_r | **2.11e-4 kg·m²** | user (confirmed) |
| Rotor CG rel. O_b | — | (0.032, 0.067) mm ≈ (0.032, −4.433)+O_b | user; O_b=(0,−4.5) |
| Rotor gas moment about O_b | M_gas | `true_gas_load(...).rotor_moment_about_center_nm` | §4.5 |
| Crank-pin journal radius / length | r_j, L_j | 14.2 / 21.0 mm | `journal_bearing` |
| Journal radial clearance | c_j | 15 µm | `JOURNAL_CLEARANCE_M` |
| Bush piece mass (each) | m_p | **6.341 g** (confirmed; grid gives 6.33 g) | user |
| Bush piece inertia about **CG** | I_p | ≈ 6.98e-8 kg·m² (rescale to m_p=6.341 g) | `bush_piece_mass_properties` |
| Bush piece CG offset from groove centre | — | (5.66, 0.0) mm (toward its own side) | `bush_piece_mass_properties.cg_offset_m` |
| Bush piece OD radius | r_b | 7.970 mm | `SwingBush.piece_outer_radius_m` |
| Rotor groove radius | r_cut | 8.000 mm | `geometry.cutout_radius_m` |
| Piece shift (sets 10 µm curved film) | — | 0.020 mm | `SwingBush.piece_shift_m` |
| Curved film clearance | c_c | 10 µm | `film_thicknesses_m` |
| Flat film clearance | c_f | 10 µm | `film_thicknesses_m` |
| Curved arc half-span | β_half | 0.9260 rad (53.05°) → full 106.1° | `SwingBush.half_arc_rad()` |
| Curved arc length | 2 r_b β_half | 14.76 mm | derived |
| Flat contact length (full / near-BDC) | L_f | 11.94 mm / ~11.47 mm | `flat_contact_length_m(θ)` |
| Effective axial height | H | 20.983 mm (cyl H − 2×8.5 µm) | `cylinder_height_m − 2·AXIAL_CLEARANCE_M` |
| Groove-centre lever from O_b | ℓ_g | 25.0 mm along rotor orientation | `geometry.cutout_offset_m` |

**Reference kinematics** at crank angle θ (rigid-film / ε=0 limit), from
`prescribed_state(geometry, θ)`:
- crank pin / reference rotor centre `O_j = e_throw·(sinθ, cosθ)` = `rotor_center_m`
- reference rotor orientation `φ_orient,ref = rotor_orientation_rad` (≈ 90° ± 10.37°);
  vane-referenced swing `φ_r,ref = φ_orient,ref − 90°`
- reference groove centre `O_g,ref = cutout_center_m = (0, cutout_y)` — **always on the
  vane axis x = 0**, because the vane guides the bush. Identity to verify in code:
  `O_g = O_b + ℓ_g·(cos φ_orient, sin φ_orient)` reproduces `cutout_center_m` when
  O_b = O_j and φ_orient = φ_orient,ref.

---

## 2. Governing principle — perturb about the prescribed kinematics ✅ DECIDED (b)

The prescribed kinematics (§ above) is the **rigid-film limit**: bush glued to the
vane axis and to the groove, giving the exact ±10.37° rotor swing. With compliant
~10 µm films the true motion is that reference **plus small (µm / mrad) deviations**.
All 9 DOF are therefore modelled as **deviations about the prescribed reference**, the
same way `rotor_dynamics.py` already treats the rotor centre (state ê_j = O_b − O_j is
a ≤15 µm perturbation of the pin orbit). This keeps the system well-conditioned (films
carry only the perturbation load, not the whole swing) and gives a clean reduction
check: **zero all deviations → recover the prescribed kinematics and the §4.11 10 µm
films exactly.**

> **DECIDED (user, 2026-07):** rotor attitude φ_r = **(b) reference + small δφ_r
> deviation**. So `φ_orient = φ_orient,ref(θ) + δφ_r`, and the rotor attitude EOM
> (§6) integrates **δφ_r** (state #3), not the absolute angle. δφ_r = 0 ⇒ prescribed
> swing exactly. (Rejected: (c) free absolute — ill-conditioned; (a) prescribed 8-DOF
> — no attitude response.)

---

## 3. State vector (18 = 9 DOF + 9 rates)

| # | DOF | Meaning | Reference (deviation = 0) |
|---|---|---|---|
| 1–2 | ê_j = (e_jx, e_jy) | rotor bore centre − crank pin, O_b − O_j | 0 (bore on pin) |
| 3 | δφ_r | rotor attitude deviation, φ_orient = φ_orient,ref + δφ_r | 0 |
| 4–5 | (x_IN, y_IN) | IN piece centre O_p,IN (global) | O_p,IN,ref (§4) |
| 6 | φ_IN | IN piece attitude (rel. vane) | φ_IN,ref |
| 7–8 | (x_OUT, y_OUT) | OUT piece centre O_p,OUT (global) | O_p,OUT,ref |
| 9 | φ_OUT | OUT piece attitude | φ_OUT,ref |
| 10–18 | time-derivatives of 1–9 | | 0 (steady orbit ⇒ periodic, not zero) |

Integrate with `scipy.solve_ivp` **BDF** (stiff squeeze films), exactly as §4.13:
start from the reference at θ_start, run ~4 revolutions, keep the last revolution.

---

## 4. Piece reference position (wire carefully — reduction check catches errors)

Each piece sits in the rotor groove and against the vane. In the reference:
- **flat film 10 µm:** piece flat is 10 µm off its side of the vane;
- **curved film 10 µm:** the 0.020 mm `piece_shift` toward the vane side makes the
  minimum curved gap 10 µm (concentric gap would be r_cut − r_b = 30 µm).

So the reference curved eccentricity is **e_c,ref = piece_shift along the piece→vane
normal**, magnitude 20 µm, giving min curved film 10 µm. The reference piece centre is
`O_p,k,ref = O_g,ref − side_k · piece_shift · n̂_flat`, where n̂_flat is the flat-face
outward normal (≈ ±x̂ near TDC, rotates with φ_orient) and side_k = +1 (IN) / −1 (OUT).
**Reduction check:** at deviation 0, `arc_film_force` min film = 10 µm and
`flat_slider_film` gap = 10 µm, matching `film_thicknesses_m`.

---

## 5. Coupling map (states → film inputs), per piece k

Let n̂ = (cos φ_orient, sin φ_orient) be the O_b→O_g direction, τ̂ = (−sin, cos).

**Groove centre and its velocity** (rotor-driven):
```
O_g   = O_j + ê_j + ℓ_g · n̂
Ȯ_g  = Ȯ_j + ê̇_j + ℓ_g · φ̇_orient · τ̂        (Ȯ_j = e_throw·ω·(cosθ, −sinθ))
φ̇_orient = φ̇_orient,ref(θ) + δφ̇_r
```
`φ_orient,ref(θ)` and `φ̇_orient,ref` come from `prescribed_state` (central-difference
the orientation, as `journal_relative_speed_rad_s` does).

**Curved film (arc_film) for piece k** — bush OD in rotor groove:
```
e_c,k     = O_p,k − O_g                        (2-vector; feed as ecc_x, ecc_y)
ė_c,k    = Ȯ_p,k − Ȯ_g
entrainment Ω_c = ½ (rotor-groove + piece) surface angular speed about O_g
              → mean of rotor swing rate and piece spin: use
              crank-pin-style entrainment; leading term ½ φ̇_orient (rotor groove)
              since the piece barely spins. Confirm sign vs curved_slide_velocity.
arc_center_rad = angle of n̂ toward the piece's contact side (per side_k)
arc_half_span_rad = SwingBush().half_arc_rad()  (0.926)
radius_m = r_b, length_m = H, clearance_m = c_c (10 µm)
→ ArcFilmForce(force_x_n, force_y_n)  = F_c,k ON THE PIECE
```

**Flat film (slider) for piece k** — bush flat on fixed vane (ground):
```
approach_m δ_k     = signed piece-flat displacement toward the vane
                     = −side_k · ( component of (O_p,k − O_p,k,ref) along n̂_flat )
approach_rate      = time-derivative of δ_k  (= −side_k · Ȯ_p,k · n̂_flat)
tilt_rad  γ_k      = φ_k − φ_k,ref             (piece attitude deviation)
tilt_rate          = φ̇_k
slide_speed U_k    = piece translation along the vane (τ̂-component of Ȯ_p,k);
                     reference value ≈ flat_slide_velocity(θ)
length_m = flat_contact_length_m(geometry, θ)  (varies near BDC), height_m = H,
clearance_m = c_f (10 µm)
→ SliderFilmForce(normal_force_n F_n,k, moment_nm M_f,k)
   normal acts along +side_k·n̂_flat (off the vane); moment about piece centre
```

Gas pressure end-conditions on the films (chamber vs 4 MPa recess) are the §4.11 BCs;
for the dynamic wedge/squeeze load they enter only as the Gumbel reference (p = 0
gauge at open ends). Keep the gas-pressure normal load separate (as `bush_film.py`
already computes it) if a mean bias is wanted — decide in build.

---

## 6. The nine equations of motion

**Rotor (3):**
```
m_r ë_j   = F_gas(θ) + F_journal(ê_j, ê̇_j) + Σ_k (−F_c,k) + m_r ω² O_j
I_r δφ̈_r = M_gas(θ) + M_journal_fric + Σ_k [ (O_g − O_b) × (−F_c,k) ]_z
```
- `F_gas`, `M_gas` = `true_gas_load(...).rotor_force_n`, `.rotor_moment_about_center_nm`
- `F_journal` = existing `_film_force_on_rotor` (Ocvirk journal reaction on rotor, §4.13)
- `−F_c,k` = reaction of each curved film on the rotor groove, applied at O_g;
  lever (O_g − O_b) = ℓ_g·n̂ gives the swing moment (the 25 mm arm)
- `M_journal_fric` = eccentric journal friction torque (`eccentric_friction_torque_nm`),
  opposing the rotor spin; small.
- **Inertia note:** m_r ω² O_j is the centrifugal (orbital) term (§4.8, ~44 N). If the
  CG offset from O_b matters for I_r φ̈, account for the (0.032, 0.067) mm offset;
  it is tiny — first cut may use I_r about O_b ≈ I_r,CG.

**Each bush piece k (3 each → 6):**
```
m_p (ẍ_k, ÿ_k) = F_c,k (curved, arc_film)  +  F_f,k (flat normal, slider)  [+ gas bias?]
I_p φ̈_k       = M_f,k (flat slider moment about piece centre)
                 + (curved film gives no moment about O_p — radial through O_p)
```
- F_f,k normal is along +side_k·n̂_flat; its line of action + M_f,k set the attitude.
- Weight/Coriolis of the ~6 g piece are negligible vs kN-less film forces but the
  piece loads are small — keep the arithmetic exact, don't drop terms silently.

---

## 7. Outputs (`RotorBushOrbit` dataclass)

Over the final revolution: ê_j(θ), δφ_r(θ), each piece (x,y,φ)(θ); min curved & flat
films per piece; **bush-film friction power** (dynamic) = shear over each film
integrated over the cycle — the Stage 6 deliverable, to compare against the quasi-static
`friction_power_cycle_w` (`bush_film.py`) and the journal loss. Report peak
eccentricities and minimum films as safety indicators.

---

## 8. Validation / reduction checklist (write these tests)

1. **Rigid-film reduction:** freeze films very stiff (or zero deviations as IC with zero
   load perturbation) → recovers prescribed kinematics; O_g stays on x≈0, φ_orient
   matches `prescribed_state`.
2. **Film reduction:** at zero deviation, `arc_film`/`slider` min films = 10 µm
   (matches `film_thicknesses_m`).
3. **φ_b = 0 recovery:** setting piece attitude + curved-relative spin to the §3.6
   fixed-attitude values reproduces `bush_film.curved_slide_velocity` /
   `flat_slide_velocity` (already the §4.11 guarantee).
4. **Journal sub-check:** with the bush pieces frozen at reference, the rotor ê_j orbit
   must match `integrate_rotor_orbit` (§4.13) to tight tolerance — proves the extension
   didn’t perturb the validated rotor rung.
5. **Energy sanity:** total dynamic bush friction is positive and O(quasi-static
   `friction_power_cycle_w`); no negative dissipation.
6. **Periodicity:** last-revolution start ≈ end (steady orbit reached).
7. **Grid/tolerance:** orbit stable under n_beta, n_s, and rtol refinement.

---

## 9. Open decisions to confirm BEFORE coding

1. **φ_r formulation** (§2): ✅ **DECIDED — (b) reference + small δφ_r deviation.**
2. **Curved-film entrainment Ω_c** sign/magnitude: leading ½ φ̇_orient (rotor groove
   swing) vs including piece spin; cross-check against `curved_slide_velocity`.
3. **Gas-pressure normal bias** on the films: include the mean gas load
   (`bush_film._film_face` normal) as a static bias, or model only the hydrodynamic
   wedge/squeeze? First cut: hydrodynamic only, note the omission.
4. **Piece inertia reference** (about CG with the 5.66 mm offset vs about piece centre):
   use CG inertia + parallel-axis to the piece centre for the φ̈ equation.
5. **I_r reference** (CG vs O_b): tiny offset; first cut O_b.

---

## 10. Exact API surface (verified signatures)

```
arc_film_force(ecc_x_m, ecc_y_m, ecc_dot_x_m_s, ecc_dot_y_m_s, entrainment_speed_rad_s,
    *, arc_center_rad, arc_half_span_rad, radius_m, length_m, clearance_m,
    viscosity_pa_s=0.010, n_beta=361) -> ArcFilmForce(force_x_n, force_y_n,
    min_film_thickness_m, max_pressure_pa)

flat_slider_film(approach_m, tilt_rad, approach_rate_m_s, tilt_rate_rad_s,
    slide_speed_m_s, *, length_m, height_m, clearance_m, viscosity_pa_s=0.010,
    n_s=201) -> SliderFilmForce(normal_force_n, moment_nm, min_film_thickness_m)

solve_line_pressure(film, source, step) -> pressure   # shared 1-D BVP core

true_gas_load(geometry, theta, *, trace) -> .rotor_force_n (fx,fy),
    .rotor_moment_about_center_nm, .rotor_force_mag_n
prescribed_state(geometry, theta) -> .rotor_center_m, .cutout_center_m,
    .rotor_orientation_rad
crank_pin_entrainment_speed_rad_s(geometry, theta);  journal_relative_speed_rad_s(...)
_film_force_on_rotor(e_x,e_y,v_x,v_y, entrainment, clearance)  # from rotor_dynamics
eccentric_friction_torque_nm(eps, omega_rel)                   # journal friction
bush_piece_mass_properties(geometry, side, *, density_kg_m3, bush, grid_step_m)
    -> mass_kg, inertia_about_cg_kg_m2, cg_offset_m, height_m
SwingBush(): .piece_outer_radius_m, .half_arc_rad(), .piece_shift_m
film_thicknesses_m(geometry, bush) -> (flat, curved)
flat_contact_length_m(geometry, theta, bush);  full_flat_contact_length_m(bush)
```

Current gate: **335 passed, ruff clean.** Suite before Stage 5: line_reynolds (5),
arc_film (8), slider_film (8) among them.
