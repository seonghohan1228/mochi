# Changes

Substantive changes are recorded here with their reason. Keep entries concise
and append new work under `Unreleased`; Git remains the detailed history.

## Unreleased

### Added

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
