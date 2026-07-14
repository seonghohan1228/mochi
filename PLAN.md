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

- [ ] Draw and document the bushing and neighboring compressor components.
- [ ] Define the global/body coordinate frames and crank-angle convention.
- [ ] Choose degrees of freedom, constraints, and required outputs.
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

- [ ] Implement and validate prescribed compressor kinematics.
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
