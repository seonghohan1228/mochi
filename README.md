# mochi

`mochi` is an early-stage Python simulation program for a rotary compressor.
The currently implemented feature is a small prescribed-motion visualization
for checking the cylinder, top-mounted vane, and rotor geometry. It is intended to be
used by the broader simulation; it is not the main numerical solver.

The broader goal remains a simulation whose equations, assumptions, units, and
numerical error can be reviewed independently of its Python implementation.

## Quick start

Requires Python 3.11 or newer.

PowerShell:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
python -m mochi --version
```

Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
python -m mochi --version
```

Before opening a pull request, run all local checks:

```bash
ruff format --check .
ruff check .
mypy src/mochi
pytest
```

## Prescribed-motion test GUI

Launch the interactive Tkinter prototype:

```bash
python -m mochi gui
```

The initial inputs are:

| Quantity | Value |
|---|---:|
| Cylinder ID | 77.0 mm |
| Rotor OD | 68.0 mm |
| Eccentricity | 4.5 mm |
| Physical speed | 30 Hz = 1800 rpm |
| Display slow factor | 0.01 |
| Circular cutout radius | 8 mm |
| Cutout center distance | 25 mm from rotor center |
| Vane width | 8 mm |
| Vane-tip calibration | 9 mm from the rotor center when the rotor is at the top |

The ports are drawn as angular windows on the bore, not as single markers:
the suction port opens at phi = 10.4 degrees and closes at beta = 27.7, and
the discharge port spans gamma = 7.2 degrees ending delta = 13.2 degrees
before top dead center. All four are editable in the input panel. The status
bar reports which port phase the current crank angle is in and how much of
the discharge port is still open to the compression chamber, which the
rotor-cylinder contact wipes shut over the port's own width. See
`PHYSICS.md` section 3.4. Angle zero is the top rotor
position and positive angle advances clockwise. The vane is integral with the
cylinder and drawn with its fixed 25 mm length; its unfilled outline may
overlap the rotor silhouette below the circular cutout, which the real parts
avoid through the axial stepped structure documented in `PHYSICS.md` section
3.3. The rotor is drawn from `mochi.rotor_profile`, which builds the
confirmed asymmetric mouth: the inlet side of the outside diameter is cut by
a straight and both lips run tangent-continuously into the circular groove. The GUI shows prescribed geometry only and does not yet
calculate pressure, contact force, or torque.
The display marks the cylinder center `C`, the moving rotor center `R`, and
the rotor-center orbit as a dashed circle of radius equal to the eccentricity.

The status bar also reports the suction (IN side) and discharge (OUT side)
chamber cross-section areas and volumes at the current crank angle, using
the supplied 21 mm axial cylinder height. These come from the circular-rotor
approximation documented in `PHYSICS.md` section 3.1: the rotor is treated
as its full outside-diameter disc and the crescent between rotor and
cylinder is split by the vane and the rotor-cylinder tangency. The status
bar also shows the accepted absolute R410A chamber pressures
(`mochi.chambers.chamber_pressures`, `PHYSICS.md` section 3.2): the suction
chamber stays at the 0.82 MPa suction-port pressure while the compression
chamber rises polytropically (effective R410A exponent n = 1.07) from
0.82 MPa at seal-over until the discharge valve opens at about 3.40 MPa
(3.24 MPa port pressure plus a 5 % valve rise), then declines linearly to
the 3.24 MPa port pressure during delivery. Near top dead center (about
+/-6 degrees for the supplied geometry) the two chambers merge into one;
the status bar reports this seal-over mixing window, during which the
merged region ramps linearly from 3.24 MPa down to the suction pressure,
reaching it exactly at 360 degrees. Keep the
`Lock e` checkbox enabled: with any smaller eccentricity the rotor no longer
seals on the cylinder and no chamber split is defined.

`mochi.chamber_volume` measures the same chambers on the real boundary
instead — the rotor contour, both swing-bush pieces, the stepped vane's
axial bands, and the R2.1 vane-root blends — which the circular-rotor
approximation cannot do. It supplies the 0.165 cm3 clearance volume that the
recompression phase of `mochi.chambers.port_timed_pressures` starts from,
which is rotor mouth cavity gas. Each crank angle costs an area integration,
so these volumes are not drawn per frame.

The `Stop at crank angle (deg)` input with its `Stop at the crank angle above`
checkbox pauses the animation exactly at a chosen crank angle. The animation
stops the first time the clockwise motion reaches that angle; pressing `Start`
again completes one more full revolution back to the same angle. This is a
display control only and does not change the physical speed.

## Bush clearance / rotor-mouth animation

A separate viewer animates the **coupled 9-DOF rotor+bush orbit** over one revolution:

```bash
python -m mochi bush-gui                 # interactive Tk window (computes the orbit, ~minutes)
python -m mochi bush-gui --gif out.gif   # render the animation to a GIF and exit
python -m mochi bush-gui --cache orbit.npz   # cache/reuse the orbit so relaunches are instant
```

It shows two panels per crank angle: a zoomed-out, true-scale **main** view of the real
rotor contour at the mouth (`mochi.rotor_profile.rotor_contour`, whose groove is part of the
rotor shape), the fixed vane, and the two swing-bush pieces moving; and a **x60 exaggerated
inset** — the `assembly/bush_clearance` view — with the curved (rotor-groove) and flat (vane)
films and the centre offsets `O_g`, `O_p(IN/OUT)`. The Tk window has a crank-angle slider,
play/pause and a *Save GIF* button. This is `mochi.bush_gui`, independent of the
prescribed-motion GUI above. Needs the plotting extra (`pip install -e ".[viz]"`).

## Collaboration model

`main` is the protected, reproducible project history. It should contain only
code that both collaborators have chosen and that passes the agreed checks.
Neither student develops directly on `main`.

When two implementations are developed separately, do not merge both branches
wholesale. Create a short-lived `integration/<scope>` branch from current
`main`, cherry-pick or manually port only the selected commits, verify the
combined solver, and review that integration branch before it reaches `main`.
The unselected implementations remain available in their original branches or
Git history.

See [CONTRIBUTING.md](CONTRIBUTING.md) for commands, branch names, conflict
handling, and the exact selection workflow.

## Documentation map

- [PHYSICS.md](PHYSICS.md): governing equations, coordinates, units,
  assumptions, numerical method, and validation contract.
- [ARCHITECTURE.md](ARCHITECTURE.md): module boundaries and data flow.
- [CONTRIBUTING.md](CONTRIBUTING.md): the two-person Git and review workflow.
- [PLAN.md](PLAN.md): staged development plan and unresolved decisions.
- [CHANGES.md](CHANGES.md): user-visible changes and why they were made.
- [AGENTS.md](AGENTS.md): repository instructions for coding agents.

If code and documentation disagree, treat that as a defect. Physics changes
must update `PHYSICS.md`; workflow or public-interface changes must update the
corresponding guide in the same pull request.

## Project layout

```text
src/mochi/             installable Python package
  cli.py               command-line entry point
  kinematics.py        prescribed rotary mechanism geometry and motion
  rotor_profile.py     rotor contour: OD flat and asymmetric mouth lips
  ports.py             port angular windows and characteristic angles
  chambers.py          crank-angle chamber cross-sections (circular rotor)
                       and the port-timed pressure phases
  chamber_volume.py    chamber volumes from the true rotor contour, swing
                       bush, and stepped vane
  gui.py               Tkinter test GUI and animation
  __main__.py          `python -m mochi` entry point
scripts/
  generate_results.py  regenerate the git-ignored results/ figures
tests/                 fast unit and regression tests
.github/workflows/     checks run for branches and pull requests
AGENTS.md               instructions for automated coding agents
PHYSICS.md              model and numerical-method contract
ARCHITECTURE.md         intended package boundaries
CONTRIBUTING.md         collaboration and Git workflow
PLAN.md                 development stages and open decisions
CHANGES.md              append-only change history
```

Generated simulations belong under `results/` and are ignored by Git. They are
rebuilt from the current model by `python scripts/generate_results.py` (which
needs the `viz` extra: `python -m pip install -e ".[viz]"`); `--prune` deletes
any figure that is no longer produced, so `results/` only ever holds the latest
render. See AGENTS.md, "Generated and local files". Small, reviewed reference
data may later be committed under `tests/data/` with its origin, units, and
generation procedure documented.

## Current scope

The visualization prescribes the rotor-center orbit and keeps the internal
circular-cutout center on the stationary top-vane axis. It does not calculate
pressure, force, contact, or torque. Those models belong to the broader
simulation and remain outside the accepted implementation until their equations
and validation cases are documented in `PHYSICS.md`.

## License

This project is released under the terms in [LICENSE](LICENSE).
