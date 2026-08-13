# Raw data export (Tecplot ASCII `.dat`)

The images under `results/` are **illustration only**. The numbers behind them — film
thickness, per-part forces and moments, validation curves — are exported as Tecplot ASCII
`.dat` under `results/data/`, so research post-processing happens in Tecplot on the raw data
rather than off a picture. Emitted by `scripts/generate_results.py` from the `DATASETS`
manifest (parallel to the `FIGURES` manifest).

## Running

```
python scripts/generate_results.py            # figures only (default, unchanged)
python scripts/generate_results.py --data      # figures + raw .dat data
python scripts/generate_results.py --data-only # only the .dat data (skips figures + prune)
```

The coupled 9-DOF orbit (a few minutes) is integrated **once** and shared by every
orbit dataset and figure in a run.

## Format

Tecplot ASCII, `F=POINT`. One file is one dataset with one shared `VARIABLES` list. Two zone
shapes:

- **1-D point zone** — one row per sample (a quantity vs crank angle / eccentricity).
- **Transient strand** (many 1-D line zones) — one line zone per crank angle, all sharing
  `STRANDID=1`, each tagged with `SOLUTIONTIME = crank angle in degrees` (so Tecplot's time
  readout shows the crank angle directly, 0–360). Each frame is `h_um` vs position, so
  Tecplot's Solution-Time animation steps through crank angle (x = position, y = film
  thickness). Used by the film files; the zone title also carries the angle (`theta_XX.Xdeg`).

Angles are in degrees, lengths in µm (films) / mm (part positions), forces in N, moments in
N·m; the crank angle (= solution time) is sorted ascending so the strand plays in order.
Physical time, if needed, is `deg/360 / freq` (`= theta_rad / ω`).

## Manifest (`results/data/`)

| File | Zone | Variables |
|------|------|-----------|
| `orbit/state_timeseries.dat` | 1-D, I=samples | `theta_deg`, rotor `e_jx/e_jy_um`, `dphi_r_mrad`, each piece `x_mm/y_mm/phi_mrad`, four film `*_um`, `seal_normal_n`, `seal_penetration_um` |
| `orbit/force_timeseries.dat` | 1-D, I=samples | `theta_deg` + 28 per-part channels: rotor forces (gas/journal/bush/seal/centrifugal `x,y`), rotor moments (gas/journal-friction/bush-curved/seal), each piece's curved force `x,y`, flat normal, flat/shear/asperity moments |
| `orbit/film_journal.dat` | transient, 1 line/crank angle | `phi_deg`, `h_um` |
| `orbit/film_curved_{in,out}.dat` | transient, 1 line/crank angle | `beta_rel_deg`, `h_um` |
| `orbit/film_flat_{in,out}.dat` | transient, 1 line/crank angle | `s_mm`, `h_um` (actual pad position; its range varies with θ) |
| `validation/reynolds_curved.dat` | 1-D | `eps`, `F_long_bearing_N`, `F_arc_film_N` |
| `validation/reynolds_flat.dat` | 1-D | `tilt_mrad`, `W_incline_N`, `W_slider_N` |
| `validation/reynolds_journal.dat` | 1-D | `eps`, `F_ocvirk_N`, `F_reynolds_1d_N` |

The per-angle force/moment channels come from `RotorBushOrbit.sample_channels` (populated by
`mochi.rotor_bush_dynamics.integrate_rotor_bush_orbit`). The writer lives in
`mochi.tecplot`. Add a dataset by writing a `data_*` function and registering it in
`DATASETS`.
