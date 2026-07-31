# AGENTS.md

This file is the repository-level operating guide for coding agents. It is
intentionally named `AGENTS.md` (uppercase) so agent tools discover it.

## Project overview

`mochi` is a Python simulation program for a rotary compressor. The current GUI
is a supporting prescribed-motion visualization, not the main solver. The wider
simulation is in the model-definition stage. Do not infer missing force laws,
material properties, or boundary conditions from the visualization or project
name.

The hierarchy of truth is:

1. Tests and reviewed reference cases define verified behavior.
2. `PHYSICS.md` defines the accepted mathematical and numerical model.
3. `ARCHITECTURE.md` defines module boundaries and data flow.
4. `PLAN.md` records intended work, not implemented behavior.

If these disagree, report and resolve the mismatch rather than silently
choosing one.

## Setup and checks

```bash
python -m venv .venv
python -m pip install -e ".[dev]"
ruff format --check .
ruff check .
mypy src/mochi
pytest
```

On PowerShell, activate the environment with
`.venv\Scripts\Activate.ps1`. On POSIX shells, use
`source .venv/bin/activate`.

Run the smallest relevant test during development and the full check set before
declaring a task complete. Do not weaken a tolerance or delete a regression
case merely to make a change pass.

## Git discipline

Use one coherent task per branch. The default branch prefixes are:

- `feature/<owner>/<scope>` for accepted-direction product work.
- `experiment/<owner>/<scope>` for competing or uncertain approaches.
- `fix/<owner>/<scope>` for defects.
- `integration/<scope>` for assembling selected changes from other branches.
- `docs/<owner>/<scope>` for documentation-only work.

Rules:

- Keep `main` clean, protected, and reproducible. Do not develop directly on it.
- Start a task from current `main`; do not base new work on an unrelated branch.
- Preserve unrelated user changes in the worktree.
- Keep commits small enough that a selected change can be cherry-picked without
  importing rejected work.
- Never rewrite shared `main` history or force-push a shared branch.
- Do not push, merge, delete remote branches, or open a pull request unless the
  user asked for that external action.
- Before integration, update the branch from `main`, run relevant checks, and
  update affected documentation.
- When choosing between parallel implementations, assemble only the accepted
  commits on an `integration/<scope>` branch. Record source branch names and
  commit hashes in the pull request.
- Delete short-lived branches after their useful changes are integrated and the
  collaborators agree that the comparison history is no longer needed.

See `CONTRIBUTING.md` for the complete human workflow.

## Physics and numerical standards

- Use SI units internally. Put units in names only at input/output boundaries
  when they prevent ambiguity; document every conversion.
- Define coordinate frames, positive directions, angle origins, and force-on-
  body conventions before implementing an equation.
- Use `float64` arrays by default. Any different precision requires a measured
  reason and a regression test.
- Keep state variables explicit. Avoid module-level mutable state and hidden
  unit conversions.
- Separate model equations from time integration, configuration, and file I/O.
- Express forces as independently testable contributions and state which body
  each force acts on.
- Guard invalid geometry and nonphysical inputs early with useful errors.
- For every new model term, add its equation, assumptions, parameter units,
  source or derivation, and validity range to `PHYSICS.md`.
- Validate numerical changes with an analytic solution, manufactured solution,
  trusted reference data, or a documented refinement study.
- Make convergence tolerances and failure visible. Do not silently return an
  unconverged state.

## Python conventions

- Support Python 3.11 and newer on Windows and Linux.
- Use type hints for public and internal interfaces.
- Prefer small pure functions and immutable configuration/state records.
- Use NumPy vector operations where they improve clarity; do not obscure a
  simple equation for premature optimization.
- Keep public APIs narrow. Re-export only deliberate user-facing names from
  `mochi.__init__`.
- Do not add a dependency when the standard library or an existing dependency
  provides a clear solution.
- Use snake_case for modules, functions, and variables; PascalCase for classes;
  uppercase snake case for true constants.
- Comments should explain a physical or numerical reason, not restate syntax.

## Intended architecture

The package should evolve toward these boundaries without creating empty
modules in advance:

```text
configuration -> model/state -> physics/forces -> numerics/integrators
                                           \-> simulation orchestration -> I/O
CLI -> configuration and orchestration only
```

Core physics must not import the CLI or output layer. Integrators operate on a
derivative or residual interface and must not know compressor-specific details.
See `ARCHITECTURE.md` before adding a new top-level module.

## Documentation responsibilities

- Update `README.md` for setup, dependencies, commands, layout, or public usage.
- Update `PHYSICS.md` for equations, assumptions, coordinates, algorithms, or
  tolerances.
- Update `ARCHITECTURE.md` when module ownership or data flow changes.
- Update `PLAN.md` when scope or milestones change; retain completed history.
- Append a meaningful entry to `CHANGES.md` for substantive changes.

## Generated and local files

Do not commit virtual environments, caches, coverage output, local editor
settings, generated `results/`, or large simulation artifacts. Reference data
must be small, intentional, traceable, and stored under `tests/data/`.

### `results/` figures: keep only the latest render

`results/` is git-ignored and never shared, so it must never accumulate a
figure that pictures a superseded geometry. The rule is that every image in
`results/` is a current render, and anything that is not is removed rather than
kept "just in case."

- `scripts/generate_results.py` is the single source of truth for the gallery.
  Each figure is one renderer in its `FIGURES` registry, rendered from the
  current `mochi` geometry and physics, so after any geometry or model change
  the figures are regenerated, not edited by hand. It needs the plotting extra:
  `python -m pip install -e ".[viz]"`.
- The script's manifest (`GENERATED_FIGURES`) is the figures it produces plus a
  short, explicit `PENDING_LEGACY_FIGURES` list of older hand-made figures that
  are kept until their renderers are ported. `python scripts/generate_results.py
  --prune` regenerates the produced set and deletes every image *not* in the
  manifest; without `--prune` those are only listed. Because `results/` is not
  shared, deleting an unmanifested render loses nothing reproducible.
- Port legacy figures into `FIGURES` one at a time: add the renderer, remove the
  filename from `PENDING_LEGACY_FIGURES`. When a legacy figure is obsolete
  rather than portable, drop it from `PENDING_LEGACY_FIGURES` so `--prune`
  reclaims it. A figure with no renderer and not on the pending list is, by
  definition, neither reproducible nor current.

## Definition of done

A change is complete when the intended behavior is implemented, relevant tests
pass, the full local check set passes, physics and user documentation agree
with the code, failure modes are explicit, and the diff contains no unrelated
or generated files.
