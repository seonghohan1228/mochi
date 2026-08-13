# Changes

Substantive changes are recorded here with their reason. Keep entries concise
and append new work under `Unreleased`; Git remains the detailed history.

## Unreleased

### Fixed

- **Cavitation was skipped whenever a film had gas boundary conditions**
  (`mochi.slider_film`, `mochi.arc_film`). The clamp was switched off on the assumption that
  "high gas pressure floods the film and suppresses cavitation"; evaluating the real state
  showed the total **absolute** pressure reaching **-11.6 MPa** on the diverging half-stroke
  of the reciprocating bush pad — worth -2310 N of spurious suction. Cavitation is a
  constraint on the *total* field, so both solvers now solve the hydrodynamic problem
  unclamped, superpose the Poiseuille gas bias (Reynolds is linear in `p`), and clamp the sum
  at the new `cavitation_pressure_pa` (same gauge as the end pressures, default `0.0` =
  classic Gumbel, so the validated pure-hydrodynamic path is bit-for-bit unchanged). The
  floor now follows one rule across all three films — *a film cannot fall below the lowest
  pressure it is connected to* — which degenerates to the journal's existing `p_cav = ambient`
  (PHYSICS.md 4.10). Bush loss in the 3-option regime falls 91 -> 6 W and, more importantly,
  the flat film finally **thickens with speed** instead of thinning. See
  `docs/bush_film_revision_2026-08.md`.

### Added

- **Curved-film gas boundary and sealing land** (`mochi.rotor_bush_dynamics`, opt-in with
  `gas_film_boundary`). The bore->chamber drop was imposed on the flat film only, and the arc
  was integrated over the whole piece span. Per PHYSICS.md 3.3 the curved film lives on the
  **86.5 deg sealing land**, whose ends open to the rotor mouth (the piece's chamber) and to
  the 4.0 MPa recess channel; both are now modelled (mirrored between the pieces). The
  resulting +367 N/piece away from the vane largely cancels the flat film's -400 N squeeze.
- **Greenwood-Tripp asperity contact** (`mochi.asperity_contact`) replacing the arbitrary
  `P0 = 50 MPa` exponential scale, driven by measurable surface parameters (sigma, beta, eta,
  E'), and applied to **both** bush films by load sharing instead of the hard eccentricity
  clamp. Adds `curved_asperity` / `asperity_params` toggles for isolation studies.
- **Bush friction diagnostics**: the flat-film loss is split into its boundary
  (asperity-borne, linear in `FLAT_BOUNDARY_COEFF`) and viscous parts, with the cycle-mean
  asperity load, on `RotorBushOrbit` — the split that identified the boundary term as the
  dominant channel and ruled out a `mu_b` calibration.

- **Raw-data export for research post-processing (`mochi.tecplot`,
  `generate_results.py --data`).** The `results/` PNGs are illustration only; the numbers
  behind the coupled-orbit and validation figures are now exported as Tecplot ASCII `.dat`
  under `results/data/` (1-D point zones and 2-D ordered surface zones). Datasets: the
  coupled orbit's per-crank-angle kinematics/films/seal and its 28 per-part force/moment
  channels (`RotorBushOrbit.sample_channels`, newly retained per angle), the three films'
  thickness fields over the whole cycle (crank-angle × position), and the 1-D-Reynolds
  validation curves. Run with `--data` (figures + data) or `--data-only`; the coupled orbit
  is integrated once and shared. See `docs/data_export.md`.
- **Bush clearance / rotor-mouth animation (`mochi.bush_gui`, `mochi bush-gui`).**
  A separate viewer from `mochi.gui`: animates the coupled 9-DOF orbit over one
  revolution with a two-panel figure — a zoomed-out, true-scale main panel showing the
  real rotor contour at the mouth (`rotor_profile.rotor_contour`, whose groove is part of
  the rotor shape), the fixed vane and the two swing-bush pieces moving; and a x60
  exaggerated inset (the `assembly/bush_clearance` view) with the curved/flat films and the
  `O_g`/`O_p(IN/OUT)` centre offsets. Tk window (matplotlib embedded) with a crank-angle
  slider, play/pause and *Save GIF*; `--gif` renders `bush_clearance_cycle.gif` headless.
  The orbit is computed once (optional `.npz` cache). `mochi.gui` is unchanged.
- **Stage 5 — 9-DOF rotor + two-piece swing-bush multibody dynamics**
  (`mochi.rotor_bush_dynamics`, PHYSICS.md 4.14). The two bush pieces become
  independent bodies coupled to the rotor by the curved (`arc_film`) and flat
  (`slider_film`) oil films; an 18-state stiff BDF integration of deviations about
  the prescribed kinematics. Finding: the rotor gas moment (~22 N*m) is reacted
  against the fixed vane through the bush, driving the curved film near contact at
  peak (~0.4 um). The **dynamic bush-film friction ~0.6 W is the value of record**
  (the 0.2 W quasi-static `bush_film` estimate is kept only for comparison).
- **Rotor-cylinder seal contact coupled into the 9-DOF EOM** (`seal_contact=True`,
  now the default; PHYSICS.md 4.14). The compliant Hertz line contact + boundary
  friction act on the rotor lateral/attitude EOM, so the crank-pin journal, the swing
  bush, and the cylinder share the load by their true stiffnesses. This roughly
  **doubles the sealing load vs the free-rotor model** — N_c mean ~296 N, and the
  rotor-cylinder loss rises to **~11 W (the value of record)**, the largest single
  mechanical loss. Attribution: the bush reaction is *inward*; tying the rotor to the
  fixed vane loads the journal (+323->+709 N outward), which presses the rotor ~1 um
  deeper into the bore. Seal engagement scales with omega^2 (centrifugal), so the
  coupling must stay on for any operating-speed sweep. Validated: revolution-converged
  N_c, closed radial force balance, undistorted attitude. New result figures:
  `bush_film/reynolds_{curved_vs_long_bearing, flat_vs_incline_slider, journal_vs_ocvirk}.png`
  (1-D Reynolds vs the analytical model for each film, error <= 1e-4),
  `bush_film/film_clearance_{journal, bush_curved, bush_flat}.png` (each film's clearance at
  the most-loaded crank angle), and
  `bearing_load/friction_dynamic_vs_quasistatic.png` (per-film quasi-static vs dynamic:
  bush x3.3, journal x1.0, seal x3.7, total x1.6). Note (open item): the bush films use
  p=0 Dirichlet ends, so the chamber discharge/suction/crank gas pressures are not yet
  imposed as boundary conditions (gas Poiseuille bias + throughflow leakage omitted).
- **Gas-pressure film boundary conditions** (`mochi.line_reynolds.poiseuille_bias_pressure`;
  `arc_film_force`/`flat_slider_film` `pressure_start_pa`/`pressure_end_pa`; PHYSICS.md 4.11).
  The 1-D film solvers now impose the chamber/crank gas pressures at the film ends,
  superposing the Poiseuille gas-bias field (Reynolds is linear in p; full-Sommerfeld under
  gas flooding) and reporting the throughflow leakage. Available in the coupled 9-DOF EOM
  as an opt-in finding (`integrate_rotor_bush_orbit` `gas_film_boundary=True`, default
  False; `bore_pressure_pa=4 MPa`):
  the piece is immersed in the bore gas, so referenced to it (divergence theorem) only the
  flat vane-sealing film carries a bore->chamber drop (IN=suction, OUT=compression). **The
  effect is dominant** — the bore gas presses each piece onto the vane, driving the flat
  film from ~1.8 um to the 0.20 um contact floor and raising the bush friction from 0.67 to
  ~5.5 W (x8); the rotor orbit and rotor-cylinder seal are unchanged. The flat film reaches
  metal contact, so this value is contact-clamp-limited: the vane-bush flat interface now
  needs an EHL/contact treatment (as the rotor-cylinder seal got). `gas_film_boundary=False`
  recovers the pure hydrodynamic films.
- **Long-bearing (Sommerfeld) analytical model** (`mochi.long_bearing`, PHYSICS.md
  4.11) — the infinitely-long journal counterpart to the short-bearing Ocvirk model,
  and the closed-form benchmark for the curved bush film: `arc_film` reproduces it to
  ~1e-7 over eps 0.2-0.9, with the half-Sommerfeld attitude invariant
  tan(phi)=pi*sqrt(1-eps^2)/(2 eps). Also added the parallel-plate squeeze closed-form
  check for `slider_film`.
- **Rotor-cylinder sealing contact** (`mochi.rotor_cylinder`, PHYSICS.md 4.15). The
  swing rotor has net-zero rotation per revolution, so the seal *slides* (mean |v|
  ~0.95 m/s), not rolls. `integrate_sealing_contact` resolves the statically
  indeterminate contact force self-consistently via a compliant Hertz line contact in
  the rotor EOM (removing the ~6 um free-orbit bore-penetration artifact), and the
  EHL line-contact film (`ehl_film_thickness_m`) sets the lubrication regime. At the
  Ra 0.3 um design finish (composite RMS sigma = 0.53 um) the film parameter is
  ~0.45 (boundary), giving a **rotor-cylinder loss ~6.7 W** in this **standalone
  4-state free-rotor** model — now a cross-check; the coupled 9-DOF value of record is
  ~11 W (see the seal-coupling entry above). References: Yanagisawa & Shimizu 1985; the
  2024 mixed-lubrication review (Lubricants 12(8):273); Daikin swing compressor.
- Added the first physics layers on top of the prescribed pressure/volume
  models. `mochi.indicated_work` computes the indicated work and power from the
  closed P-V loop of the port-timed cycle (`W = -oint p dV`, `P = W f`),
  integrating to the seal-over entry so the degenerate merge-window volume does
  not corrupt the loop; the supplied geometry gives about 24.6 J per revolution
  and 738 W. `mochi.bush_film` models the two 10 um swing-bush oil films (curved
  bush/groove and flat bush/vane) as incompressible Couette + Poiseuille flow:
  the flat film is driven by the bush translation (up to 0.86 m/s) and the
  curved film by the rotor's +/-10.37 deg orientation swing (up to 0.27 m/s, not
  the shaft speed), the pressure is linear between the gas boundary pressures
  (chamber pressure and the 4.0 MPa recess pressure, the first use of
  `RECESS_PRESSURE_PA`), the shear traction is mu*U/h with POE VG68
  mu = 0.010 Pa*s, and the cycle-mean film friction is about 0.2 W. Added
  `indicated_pv_diagram.png` and `bush_film_pressure.png` figures, PHYSICS.md
  sections 3.5 and 3.6, and tests. Pure Couette carries no load-bearing
  pressure, so the film pressure comes from the gas boundary conditions; a
  hydrodynamic (Reynolds/wedge) model is left for a later rung.
- Added `scripts/generate_results.py`, a reproducible renderer for the
  git-ignored `results/` gallery, plus a `viz` optional dependency
  (matplotlib, pillow) it uses. It regenerates `chamber_pressures_vs_crank.png`
  from the port-timed true-geometry-volume rule (`port_timed_pressures`,
  PHYSICS.md section 3.4) instead of the earlier circular-rotor section 3.2
  trace, so the chart now shows the re-expansion spike before phi, the
  delivery plateau held at the discharge-port pressure, and the recompression
  rise to the about 10 MPa leakage-free peak. It also regenerates
  `rotor_motion.gif` in the latest confirmed geometry — the asymmetric rotor
  mouth, the R2.1 mm vane-root fillets, and both swing-bush pieces, which the
  previous GIF predated (its title read "swing bush excluded"). Recorded the
  convention in AGENTS.md that `results/` holds only current renders, since the
  directory is git-ignored and never shared.
- Made `generate_results.py` a renderer registry (`FIGURES`) so legacy figures
  are ported one at a time, and ported `port_geometry.png` (final rotor, swing
  bush, and the dimensioned suction/discharge port windows at top dead center)
  and `tdc_bdc_definition.png` (top/bottom dead-centre chamber definition on the
  section 3.2 circular-rotor model, with the equal 476 mm² / 10.0 cm³ chambers
  computed from the current geometry rather than hardcoded). `--prune` deletes
  only images outside the manifest; older hand-made figures stay listed in
  `PENDING_LEGACY_FIGURES` and are preserved until their renderers are ported or
  they are confirmed obsolete, so no current figure is lost while the port is in
  progress. Re-created `rotor_mouth_lip_detail.png` from the current asymmetric
  contour, so it now shows the R1.5/R1.0 lips (102.2° arc, 52.6° groove blend)
  the earlier R1.4/R0.9 figure predated, and dropped `rotor_mouth_48deg.png`
  (the superseded symmetric 48°/48° mouth, whose premise the asymmetric mouth
  made obsolete) from the manifest so `--prune` reclaims it.
- Ported the remaining hand-made figures into `FIGURES`, so every kept image is
  now regenerated from the model: `geometry_master`, `dimensioned_top_view`,
  `dimensioned_side_section`, `stepped_vane_structure`, `vane_side_view`,
  `vane_model_comparison`, `bush_placement_clearances`, `chamber_case1_seal_over`,
  and `chamber_case2_gap`. The axial-section and assembly plates are drawn as
  schematics from the `AxialBands`/`SwingBush` constants (the 2-D model has no
  axial geometry); the top-view panels reuse the live `rotor_contour`. Dropped
  angle/area values (476 mm² / 10.0 cm³, the 0.5 mm gap, etc.) are computed from
  the current geometry. `PENDING_LEGACY_FIGURES` is now empty: only
  `rotor_mouth_48deg.png` remains outside the manifest, awaiting `--prune`.

### Fixed

- Drew the swing-bush R0.5 corner fillets in the GUI. `_draw_swing_bush` had
  built each piece as the outer cylindrical arc joined straight to the full
  chord intersection (`sqrt(r^2 - flat^2)`), which skipped the fillets and ran
  the flat about 0.93 mm past the true fillet-flat tangent, so the piece read
  as a plain circular segment. The outline now follows the same profile as
  `SwingBush.occupies` — outer arc, top fillet, flat between the tangent
  points, bottom fillet — matching the volume model to within grid rounding.

### Added

- Blended each vane flank into the cylinder bore with the supplied R2.1 mm
  radius instead of a sharp corner, a stress-relief fillet at the vane root.
  Added `vane_cylinder_fillet_m` to `RotaryGeometry` (with `vane_fillet_geometry`
  returning the tangent points) and drew the rounded roots in the GUI and the
  `results/port_geometry.png` figure. The blend is full axial height, so on
  each side it turns a small pocket between the flank, the bore, and the
  blend arc from gas into cylinder material; the true-geometry volume model
  now excludes it. Its only model effect is that volume: it removes about
  0.028 cm3, lowering the clearance volume from 0.193 to 0.165 cm3 (still in
  the 0.1-0.5 cm3 band) and the leakage-free recompression peak from about
  11.8 to about 10 MPa. It does not touch the seal-over geometry or the
  rotor. Documented in PHYSICS.md sections 3.3 and 3.4.
- Confirmed (A): `phi` is the suction-port opening angle, measured with
  `beta`, `gamma`, `delta` from the vane centerline, matching the reference
  rule set (R2.1) and this implementation. Added `ports.suction_window`
  alongside `discharge_window` so both ports are real `PortWindow` objects
  with start and end angles — suction `[phi, beta] = [10.4, 27.7]` degrees,
  discharge `[339.6, 346.8]` degrees — and pointed the GUI at them. Recorded
  in PHYSICS.md section 3.4 that `phi` not matching the vane band (11.93
  degrees) is not a discrepancy: a machined port location and the seal-over
  width are unrelated quantities, so the two earlier open items on this are
  resolved in section 9. Added `results/port_geometry.png` (git-ignored)
  showing the final rotor contour, both swing-bush pieces, and the two port
  windows dimensioned from top dead center.
- Turned the ports into angular windows and measured the chamber volumes on
  the true rotor boundary, which together give the cycle its missing
  clearance volume. The ports were single display angles (`+30` and `-30`
  degrees from `+y`) read only by the GUI's tick drawing; the supplied
  timing `phi` 10.4, `beta` 27.7, `gamma` 7.2, `delta` 13.2 degrees now
  lives on `RotaryGeometry` and the new `ports` module derives the four
  characteristic angles from it, plus a discharge-port open area that the
  rotor-cylinder contact wipes shut over the port's own width rather than
  switching off at one angle. `gamma` checks out (a 4.84 mm arc on the
  38.5 mm bore radius); `phi` does not follow from our 8 mm vane, which
  subtends 11.93 degrees, so it is carried as an independent supplied angle
  and kept separate from the 5.96 degree seal-over half angle.
- Added `chamber_volume`, which integrates the gas space on the real
  boundary — the rotor contour of `rotor_profile`, both swing-bush pieces,
  and the stepped vane's axial bands — instead of the outside-diameter disc.
  This was needed because the circular-rotor approximation cannot carry a
  recompression phase: at `2*pi - delta` it gives a 0.00138 cm3 trapped
  volume and drives it to zero, so `p V**n = const` diverges. The true
  geometry gives 0.193 cm3 (0.165 cm3 once the vane-root blend below is
  included), which is rotor mouth cavity gas, so the dead volume of this
  machine is the rotor mouth itself and no artificial
  volume floor is required. Cross-checked against the 30 mm2 mouth free area
  that PHYSICS.md had recorded as a known neglected volume, measured here as
  26.5 mm2 by an independent route. Grid convergence is 0.3 % between a
  0.050 mm and a 0.025 mm pitch. The swing bush is assumed to hold a fixed
  attitude and the 0.010 mm films are assumed filled; both are registered as
  open items.
- Added the port-timed pressure phases (`port_timed_pressures`) on those
  volumes: re-expansion before the suction port opens, suction, polytropic
  compression from `beta`, delivery, and recompression after the discharge
  port shuts. The valve opening angle stays a pressure condition rather than
  a constant so it tracks the operating condition. The recompression end
  pressure is the residual the next re-expansion starts from, so the cycle
  closes without hardcoding one; it reaches 11.8 MPa, which is an upper
  bound because leakage is neglected. The accepted constant-pressure rule of
  PHYSICS.md section 3.2 is untouched and remains the regression baseline;
  the two run side by side, and the port-timed compression pressures sit
  about 7 % lower at mid-stroke because the mouth cavity adds volume the
  disc discarded. Moving the compression start from 5.96 to 27.7 degrees
  shifts the valve angle by 0.35 degrees only, so the unsettled reading of
  the supplied angles does not change the compression stroke.
- Drew the ports as arcs and the two swing-bush pieces in the GUI, added the
  four port angles as inputs, and added a port-phase and open-area status
  line. The true-geometry volumes are not shown per frame because each crank
  angle costs an area integration.
- Stopped the GUI drawing from rescaling on every frame. The status text
  changed length as the crank turned, and because the window had no explicit
  size Tk resized it to the label's requested width, which changed the canvas
  size and therefore the drawing scale. The window now opens at a set size
  (still user-resizable), the status label uses a fixed-width font, and every
  status branch is padded to the same character count, so the canvas stays
  constant through a full revolution.
- Added `mochi.rotor_profile`, which turns the confirmed mouth dimensions
  into the rotor's real boundary polylines, and switched the GUI to draw it.
  The GUI previously showed a placeholder narrow slot with lip fillets on a
  plain circular rotor; it now shows the inlet OD flat and both asymmetric
  lips, and strokes only real material edges so the mouth stays open to the
  chamber instead of being closed by a drawing artefact.
- Made the rotor mouth asymmetric and corrected the lip radii to R1.5 and
  R1.0 (the earlier 1.4 and 0.9 ignored the chamfer). On the inlet side a
  straight now replaces the OD from 34.8 degrees to 13.4 degrees at radius
  33.657 mm, cutting the OD by up to 0.773 mm over a 12.566 mm chord, and
  the lip runs R1.5 x 102.2 degrees, 0.399 mm straight, R1.0 x 52.62 degrees
  into the groove at 47.74 degrees. The outlet side keeps the original
  OD-tangent design at 13 degrees, which with the new radii gives R1.5 x
  92.6 degrees, 0.744 mm straight, R1.0 x 52.43 degrees into the groove at
  47.97 degrees. Both blend sweeps match the 52 degrees on the supplied CAD
  view, and a full-revolution sweep confirms neither side touches the bore.
- Completed the swing-bush piece profile: the cylindrical face spans 106.1°
  and meets the flat through R0.5 tangent fillets rather than a sharp
  corner. The fillet radius reproduces the supplied arc exactly
  (arccos((3.990 + 0.5)/(7.970 - 0.5)) x 2 = 106.10 degrees), so the earlier
  clearance scheme is unchanged; the flat contact length becomes 11.940 mm
  and the arc 14.76 mm. The capture argument was restated on the physical
  basis (piece thickness 3.980 mm versus the 2.573 mm mouth opening) and the
  86.5 degree sealing land was noted to narrow to 84.07 degrees near the
  extreme rotor tilt.
- Recorded the supplied swing-bush geometry and clearances: each half-moon
  piece has a 7.970 mm outer radius (Ø15.94) with its flat face 3.990 mm
  from its own center, so a groove-centered piece would overlap the vane by
  0.010 mm; offsetting each piece 0.020 mm toward its own side opens both
  the flat-to-vane and curved-to-groove films to a uniform 0.010 mm (the
  offset is exactly the midpoint of the 0.020 mm free play, confirming the
  numbers are mutually consistent). The axial clearance is 0.0085 mm per
  side, so the bush height is 20.983 mm.
- Placed the swing bush in the documented geometry (two half-moon pieces
  centered on the circular groove, slot aligned with the fixed vane) and
  derived its clearance conditions: the pieces stay captured with an 8.1°
  margin at maximum tilt, the bush must span the full 21 mm height so its
  OD seals the 4.0 MPa recess bands where the channel meets the groove, the
  mouth lip adds no extra clearance condition, and the sliding/oscillation
  demands are 9 mm at up to 0.86 m/s and ±10.37° at up to 34 rad/s (30 Hz).
  The radial, slot-side, and axial clearance values remain to be supplied.
- Made the discharge-chamber pressure trace continuous over the whole
  revolution, as agreed: after the valve opens near 221 degrees the delivery
  pressure now falls linearly in crank angle from the 3.40 MPa opening
  pressure to the 3.24 MPa port pressure at the seal-over window entry, and
  the merged region then follows a linear mixing ramp down to the 0.82 MPa
  suction pressure exactly at 360 degrees (replacing the constant-delivery
  plateau and the sudden drop at the window). The GUI seal-over message now
  reports the ramping mixing pressure.

- Recorded the confirmed rotor mouth opening: the rotor OD is cut 13° to
  each side of the vane (measured at the rotor center) and each lip runs
  tangent-continuously from the OD through an R1.4 arc (92.6° sweep), a
  0.97 mm straight, and an R0.9 arc (52.3° sweep) into the circular groove.
  The exact G1 chain reproduces the CAD values (R0.9 sweep 51.9° computed
  vs 52.3° displayed) and puts the groove blend 48.5° from the groove
  center, matching the earlier 48° nominal and clearing the vane at maximum
  relative tilt (44.6° < 48.5°). The lip touches the rotor OD only at its
  tangent point, so no bore interference exists at any crank angle; earlier
  semicircular-lip variants penetrated up to 0.32 mm and were discarded.
  The GUI still draws the old narrow slot; its update is a planned item,
  and the chamber model notes the mouth gas space (about 3 % of the
  crescent) as a known neglected volume.
- Recorded the confirmed axial structure: the vane is a fixed stepped part of
  the cylinder (full 21 mm thickness to 15.4 mm from the bore, then 2.4 mm
  top/bottom ledges to the 25 mm tip) and the rotor carries a 14.7 mm-thick
  recess channel of width 8*sqrt(2) mm from the circular groove to its
  center, so the previously flagged planar vane-rotor overlap is not a real
  interference. The vane tip is now fixed in the kinematics (constant 25 mm
  length, replacing the invisible reference-line rule), TDC/BDC are defined
  (contact at/opposite the vane), the sealed recess spaces carry a fixed
  4.0 MPa absolute pressure constant for later swing-bush work, and
  PHYSICS.md section 4 records the accepted bodies, frames, sign discipline,
  and the single prescribed crank degree of freedom (7 body DOFs minus 6
  joint constraints).

- Confirmed R410A as the working fluid: the effective polytropic exponent
  became n = 1.07 (fitted to the CoolProp HEOS isentrope between the
  supplied ports, replacing the 1.3 placeholder), and the discharge valve is
  approximated to open 5 % above the discharge-port pressure (about
  3.40 MPa) with delivery held at that opening pressure.
- Fixed the supplied absolute port pressures (0.82 MPa suction, 3.24 MPa
  discharge) and replaced the constant discharge-chamber pressure with a
  capped polytropic rule chosen with the collaborators: the compression
  chamber starts at the suction pressure at seal-over, follows p V^n = const
  (default n = 1.3, refrigerant to be confirmed), and holds the
  discharge-port pressure once reached, removing the pressure jump at the
  seal-over window boundary. The GUI status bar now shows both chamber
  pressures and the valve-open state.
- Added the supplied 21 mm axial cylinder height to the geometry so chamber
  areas become volumes, and adopted the first constant chamber pressure rule:
  each chamber takes its port pressure, and inside the roughly +/-6 degree
  seal-over window both chambers take the suction-port pressure as a mixing
  process. Documented the recommendation to keep the eccentricity locked to
  e = (cylinder ID - rotor OD) / 2 because a smaller value opens a
  rotor-cylinder gap that leaves the chambers connected.
- Added crank-angle suction/discharge chamber cross-section areas using the
  circular-rotor approximation (`mochi.chambers`), as groundwork for later
  suction and discharge pressure conditions. Closed-form areas are validated
  against independent quadrature; when the chambers are not geometrically
  separated (rotor-cylinder gap, unsealed vane, or the roughly +/-6 degree
  seal-over window at top dead center) the evaluation raises a specific error
  instead of forcing a split. The GUI status bar reports the areas per frame
  without changing the drawn geometry; volumes require an explicitly supplied
  axial cylinder height.

- Added an optional stop-at-crank-angle control to the test GUI so the
  animation pauses exactly at a chosen angle for geometry inspection. The
  control affects the displayed animation only, never physical time or speed.

- Corrected the project scope to a rotary compressor
  throughout code and documentation.
- Added prescribed rotor-cutout kinematics for the supplied 77 mm cylinder ID,
  68 mm rotor OD, 4.5 mm eccentricity, 25 mm cutout-center distance, 8 mm
  cutout radius, 8 mm vane width, and a 9 mm top-position vane-tip
  calibration.
- Added a Tkinter test GUI with editable inputs, slowed animation, symmetric
  inlet/outlet markers, and physical versus displayed speed diagnostics.
- Refined the GUI to rotate clockwise from the top position, vary vane length,
  and draw the rotor plus cylinder-vane geometry as clean parametric contours
  without overpainted masking shapes.
- Rounded both rotor-opening lips with parametric tangent fillets, removed the
  lower square notch, and kept the changing-length cylinder-vane contour
  transparent so it can overlap the rotor.
- Added cylinder- and rotor-center markers plus a dashed rotor-center locus
  tied directly to the eccentricity.
- Added the initial Python package, command-line, test, and continuous-
  integration scaffold so development begins from a runnable cross-platform
  baseline.
- Added collaboration guidance for protected `main`, independent experiment
  branches, selective integration, review, and cleanup.
- Added living physics, architecture, and development-plan contracts so
  unverified compressor assumptions are not embedded silently in code.
