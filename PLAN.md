# Development plan

This is a living plan, not a description of implemented features. Retain
completed items so the development path remains reviewable.

## Phase 0 - repository foundation

- [x] Create an installable Python package shell.
- [x] Define formatting, linting, typing, and test commands.
- [x] Add Linux and Windows continuous integration.
- [x] Document the protected-main and selected-integration workflow.
- [x] Establish physics and architecture contracts.

## Phase 1 - agree the physical problem

- [ ] Confirm the cylinder, top vane, cutout, and rotor against compressor CAD.
- [x] Define the prototype global frame and crank-angle convention.
- [x] Document the supplied cylinder, rotor, circular-cutout, and vane geometry.
- [x] Define 9 mm as the rotor-center-to-vane-tip distance at the top position.
- [x] Confirm the vane-tip relation and clearances against CAD: fixed stepped
  vane integral with the cylinder, tip 25 mm from the bore; ledges run in the
  rotor recess channel (PHYSICS.md section 3.3).
- [x] Choose the bodies, degrees of freedom, and constraints: one prescribed
  crank degree of freedom; cylinder+vane ground, crankshaft, rotor, and swing
  bush with six joint constraints (PHYSICS.md section 4).
- [ ] Choose the required outputs.
- [x] Place the swing bush (two half-moon pieces on the groove center,
  slot aligned with the vane) and derive its clearance conditions: capture
  margin 8.1° at maximum tilt, full-height requirement for the 4 MPa seal,
  9 mm / 0.86 m/s sliding and ±10.37° / 34 rad/s oscillation demands.
- [x] Obtain the swing-bush clearance values: piece Ø15.94 with its flat
  face 3.99 mm from the piece center, offset 0.020 mm per side for uniform
  0.010 mm films, and 0.0085 mm axial clearance per side.
- [x] Complete the swing-bush profile: the cylindrical face keeps 106.1° and
  joins the flat through R0.5 tangent fillets, giving an 11.94 mm flat
  contact length and a 14.76 mm arc.
- [x] Add the inlet-side rotor OD flat (34.8° on the OD to 13.4° at radius
  33.657), correct the lip radii to R1.5 / R1.0, and re-solve both mouth
  lips: inlet R1.5 x 102.2° + 0.399 straight + R1.0 x 52.62° (groove
  47.74°), outlet unchanged in form with R1.5 x 92.6° + 0.744 straight +
  R1.0 x 52.43° (groove 47.97°); no bore interference on either side.
- [ ] Decide the swing-bush pressure boundary conditions (recess spaces
  fixed at 4.0 MPa absolute meanwhile).
- [x] Update the GUI rotor contour to the confirmed asymmetric mouth: the
  new `rotor_profile` module builds the inlet OD flat and both
  tangent-continuous lips, and the GUI strokes only the real material edges
  so the mouth stays open to the chamber.
- [ ] Build a parameter table with symbols, SI units, sources, and ranges.
- [ ] Inventory all forces and identify which are prescribed or solved.
- [ ] Select at least one independently verifiable reference case.

Exit criterion: both collaborators can calculate the sign and expected order of
magnitude of each force for a simple shared case without consulting code.

## Phase 2 - reference dynamics solver

- [ ] Implement validated configuration and planar state data structures.
- [ ] Implement modular force assembly with an explicit zero-force case.
- [ ] Add a SciPy reference integrator with visible diagnostics.
- [ ] Pass zero-force, constant-force, and oscillator analytic tests.
- [ ] Add result objects and a minimal versioned case/result format.
- [ ] Demonstrate tolerance/step refinement.

Exit criterion: the reference solver passes analytic cases and produces a fully
traceable result without compressor-specific hidden assumptions.

## Phase 3 - compressor model

- [x] Implement and validate the prescribed rotary-compressor test kinematics.
- [x] Add an interactive test GUI with physical/display speed separation.
- [x] Use clockwise, top-zero motion and continuous parametric rotor and
  cylinder-vane contours.
- [x] Fillet both rotor lips, end the opening at the circular cutout, and keep
  the changing-length vane overlay transparent.
- [x] Show the cylinder center, moving rotor center, and dashed eccentric
  rotor-center locus.
- [x] Compute crank-angle suction/discharge chamber cross-section areas with
  the circular-rotor approximation as groundwork for pressure conditions;
  report instead of forcing a split when the chambers are not separated.
- [x] Obtain the axial cylinder height (21 mm supplied) so chamber areas
  become volumes.
- [x] Adopt the chamber pressure rule: the suction chamber holds the
  suction-port pressure; the compression chamber rises polytropically to the
  valve opening pressure, then declines linearly to the discharge-port
  pressure at the seal-over window entry, and a linear mixing ramp brings
  the merged region down to the suction pressure exactly at 360 degrees, so
  the discharge trace is continuous over the whole revolution.
- [x] Fix the port pressures: 0.82 MPa suction and 3.24 MPa discharge,
  absolute.
- [x] Confirm the refrigerant and polytropic exponent: R410A, effective
  n = 1.07 from the CoolProp HEOS isentrope between the supplied ports.
- [x] Approximate the discharge-valve opening as a 5 % pressure rise over
  the discharge port (valve opens at about 3.40 MPa and delivery holds that
  pressure).
- [ ] Confirm the valve rise fraction against actual valve data.
- [x] Turn the ports from display markers into angular windows and derive
  the four characteristic angles from the supplied phi 10.4, beta 27.7,
  gamma 7.2, delta 13.2 degrees, with a swept discharge-port open area.
- [x] Measure the chamber volumes on the true rotor contour with the swing
  bush and the stepped vane, giving a 0.193 cm3 clearance volume that is
  99 % rotor mouth cavity; the circular-rotor approximation understated it
  by 140x and could not carry a recompression phase at all.
- [x] Add the re-expansion and recompression phases the clearance volume
  makes possible, closing the cycle on a derived residual pressure.
- [ ] Confirm the port-angle reference convention and the source of phi
  against the drawing (PHYSICS.md section 9).
- [ ] Quantify the fixed-attitude swing-bush assumption in the volume model.
- [ ] Add leakage so the recompression peak stops being an upper bound.
- [ ] Add accepted load/force models one at a time.
- [ ] Add contact and friction only after their mathematical treatment is
  selected.
- [ ] Compare competing implementations on identical cases.
- [ ] Assemble selected changes on an integration branch.

Exit criterion: the combined model meets the agreed compressor reference-case
acceptance limits.

## Phase 4 - robustness and release

- [ ] Define parameter sweeps and edge cases.
- [ ] Profile representative workloads before optimizing.
- [ ] Add conservation and failure-mode regression tests.
- [ ] Freeze and document the first stable input schema and public API.
- [ ] Tag a reproducible release.
