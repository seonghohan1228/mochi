# Reference-machine bush-film validation (2026-08) — why it does not match

Record of validating this project's bush oil-film module against a **measured** film-thickness
waveform, and of the disagreement that came out of it. The disagreement is the point of the
document: it is quantified, its candidate causes are separated into ones that were ruled out
and ones that remain, and it is stated plainly that **no model parameter was changed to close
it**.

**Sources** (both read in full, not second-hand)

| | |
|---|---|
| **T02** | Tanaka, Nakahara, Kyogoku, Toyama (2002), *Purdue ICEC* Paper 1534 — "Lubrication Characteristics Between Bush And Blade Of Swing Compressor". Analysis + rig geometry. |
| **J08** | Tanaka, Zuo, Hikam, Toyama (2008), *Trans. JSRAE* 25(4) 375–382 — eddy-current **measurement** of the same sliding pair. Same group, same manufacturer. |

**Machine**: a Daikin-type swing compressor — blade integral with the roller, bush pivoting in
a cylinder seat. That is the *mirror* of the machine this project models (fixed vane, bush
riding in the rotor groove), so the validation targets the **bush film module**, not the whole
coupled model.

---

## 1. What was established (independent of the disagreement)

* **Geometry confirmed at source.** T02 Table 3 (rig): `R_c` 22.0, `R_o` 16.8, `R_b` 6.0,
  blade `a` 4.0, `L_p` 25.0 mm. J08 Fig.2 adds the bush drawing: **R6, 25 mm axial, 10 mm
  flat height, 4 mm segment, R0.1 fillets, sensor spacing 4.5 mm**, and the machine is a
  **2.8 kW** air-conditioning unit with a high-speed-tool-steel bush on a grey-cast-iron blade.
  Note T02 **Table 1 ≠ Table 3**: Table 1 (`R_c` 30, `R_o` 26.8, `L_p` 10 mm) is the machine
  T02 *simulated*, not the rig.
* **Kinematics validated against T02 Fig.8.** Running `mochi.reference_swing` on T02's
  *simulation* machine reproduces both published sliding speeds — bush–blade **1.17 vs
  ~1.15 m/s**, bush–cylinder **0.21 vs ~0.20 m/s**. That also pins the one geometric quantity
  neither paper gives, the bush pivot distance: only `a = R_c + R_b` matches the second speed.
* **Our contact model is the field standard.** T02 eq.(17) is
  `k_c = 4.4086e-5 · E' · (4 − h/σ)^6.804` — the same Greenwood–Tripp fit, coefficient and
  exponent included, that `mochi.asperity_contact` implements.
* **Digitisation is trustworthy.** The Fig.8 digitisation reproduces the four values quoted in
  J08's text to ±0.05 µm on three of them, and — the independent check — the tilt recomputed
  from it via J08 eq.(1) lands at **0.017–0.089°** against the **0.02–0.085°** of Fig.9, peak
  at the same angle. That simultaneously confirms `L = 4.5 mm` is the *sensor spacing*, not a
  pad length.
* **Effective pad length ≤ 8 mm.** Purely geometric: with the measured tilt, a 10 mm pad (the
  drawing value) drives `h(±L/2)` negative at 5 of 72 crank angles. 6 mm is the value that
  also balances the load at BDC.

## 2. The disagreement

| source | bush–blade film |
|---|---|
| **T02's own simulation** (Fig.6) | **0.28 – 0.74 µm** |
| **this project's model** (our machine, flat film) | **0.89 – 1.18 µm** |
| **J08 measurement** (Fig.8) | **2.1 – 14.5 µm** |

Two independently written models land in the same sub-micron decade; the measurement is
**10–20× thicker**.

**The tilt, by contrast, agrees.** An earlier draft of this document claimed it did not; that
was read off the *pure-hydrodynamic* configuration, where the pad has no gas ramp to tilt it
and sits at ~10⁻³ deg. With the gas boundary and the cavitation fix in place — the current
model — the IN piece reaches **9.6e-4 rad = 0.055°**, inside the measured **0.02–0.085°**
band. (The OUT piece runs 11× flatter at 8.9e-5 rad; whether that asymmetry is physical is
open.) So the disagreement is confined to the **absolute film thickness**, not to the tilt.

Worth noting: J08 came six years after T02, from the same group, and **does not compare its
measurement to T02's simulation**.

## 3. Causes ruled out

* **Grid convergence — excluded.** At the reference conditions the flat film's force is
  converged to **0.002 %** by `n_s = 401` (thin-film case 0.007 %), and the arc film to
  **0.16 %** by `n_beta = 121`; these are the grids in production use. Structurally it cannot
  be the cause either: film force goes as `F ∝ h^-n` with `n ≈ 2–3`, so `δh/h = −(1/n)·δF/F`
  and a 10× error in `h` would need a **10²–10³×** error in force. A 0.2 % force error moves
  `h` by under 0.1 %. Two independently implemented codes agreeing to a decade also rules out
  a shared discretisation artefact.
* **Seat friction moment — VERDICT RETRACTED.** An earlier pass swept `μ_seat` over
  0 / 0.05 / 0.10 / 0.15 and reported no effect (RMS unchanged, `h_u > h_l` failing 16 of 16).
  That verdict is **withdrawn**: the harness it ran in was solving the pad's moment balance by
  root-finding per crank angle and, whenever the root fell outside the bracket, snapping the
  tilt to `±gmax` picked by whichever endpoint had the smaller residual. That is a
  discontinuous selection, and it made the tilt flip sign at **27 of 71** angle steps and
  saturate at the limit in **29 %** of them, with `h0` jumping up to 11.5 µm between adjacent
  angles. The oscillation dominated the RMS, so the sweep could not have resolved the friction
  term either way. Re-tested with a time-integrated harness (§3.1).

  The saturation was specific to that harness. `mochi.rotor_bush_dynamics` time-integrates the
  piece attitude as an ODE state, so nothing selects between branches; its `min`/`max` clamps
  on approach and tilt are monotone and, measured over a revolution, **never bind**
  (approach peaks at 6.96 / 8.47 µm against a 9.80 µm limit; tilt-clamp activations 0 of 180
  on both pieces). The one clamp that does bind is the arc eccentricity at 0.99 — already
  documented in PHYSICS.md 4.11, and shown there to be what the load path demands rather than
  an artefact.
* **Wrong pressure reference — found and fixed.** An early pass referenced the film to the
  suction chamber rather than to the bush's immersion (back) pressure. Referencing to the
  immersion, as `mochi.rotor_bush_dynamics` already does for our own piece, made the implied
  seat reaction **insensitive to the assumed pad length** (51–60 N across 4.5–8 mm, against
  84–172 N before), and gives a seat contact pressure of **0.1–0.2 MPa** — comfortably
  physical. This corrected a real error in the harness, not in the model.

## 3.1 The three-parameter fit — run properly, and failed

With the oscillating harness replaced by a time-integrated one (two bush freedoms carried as
ODE states, BDF to the periodic orbit — the way `rotor_bush_dynamics` solves ours), the
artefact is gone: **0 sign flips in 72 angles**, against 27 before, and `h_u > h_l` now holds
in every case. So the model could finally be judged on physics.

The three quantities neither paper publishes — flat reference clearance, seat clearance,
effective pad length — were swept over 3x3x3, everything else fixed at published values
(viscosity 2.8 mPa*s and the friction coefficients from T02 Table 2; geometry from T02
Table 3 and J08 Fig.2). Three parameters against **144 measured points** (72 angles x two
edges), with the waveform *shape* having to come out, is a fit with real discriminating power.

**It did not pass.** Best case `C_flat = 20 um, C_seat = 40 um, L_c = 8.0 mm`:

| | measured | model |
|---|---|---|
| `h_u` range | 5.4 – 14.5 um | 11.0 – 16.8 um |
| **`h_u` peak-to-peak** | **9.13 um** | **5.78 um** |
| `h_l` peak-to-peak | 6.43 um | 6.52 um |
| RMS (upper / lower) | — | 6.75 / 2.18 um |
| correlation (upper) | — | 0.639 |

Two things stand out and they are more informative than the RMS:

* **Amplitude, not offset, is the dominant failure.** In **24 of the 27** cases the upper
  edge moves by under 1 um over the whole revolution against a measured 9.13 um — the modelled
  bush sits nearly still. Only opening the seat clearance revives it (0.03 -> 5.78 um going
  from 10 to 40 um), which is the direction J08's own explanation points ("the bush has more
  room, so it tilts more"), but even 40 um falls short.
* **The residual is one-signed.** The model is thicker at 11 of 12 sampled angles, worst
  around 60–90 deg at −9.4 um. That is a systematic error, not scatter.

`C_flat` wants to sit at the bottom of its sweep (RMS degrades monotonically 20 -> 30 -> 40 um)
but cannot go lower without failing to contain the measured 14.5 um peak — the fit is boxed in
structurally, not short of tuning.

So the honest reading is that this is **missing physics, not unfound parameters**: with the
unknowns cut to three and the whole waveform used as the target, it still fails, and it fails
by not moving the bush far enough over the cycle. That is plausibly the same root as the
sub-micron film both models predict.

## 4. Causes still open

1. **Runout residual in the measurement.** J08 reports runout of **±0.04 V ≈ ±4.8 µm** on the
   upper sensor — the same size as the thin part of the signal it is measuring. The authors
   state the geometric correction was impossible because of clearances between parts, and fell
   back on **waveform pattern matching**. The quoted ±0.3/0.2 µm error covers **temperature
   drift only**; residual runout is not in it.
2. **What the sensor measures.** The eddy-current sensors and thermocouple sit **~100 µm below
   the bush surface, in epoxy** (J08 §2.1). Differential thermal expansion between the bush
   steel and the epoxy enters the displacement directly — the authors acknowledge it inside
   their temperature-drift correction — and the rig runs at 61 °C while the calibration
   (Fig.4) is at room temperature.
3. **Physics both models omit.** Thermal wedge, surface texture, and bulk elastic deformation
   are absent from T02's model and from ours alike. Elastic deflection is the largest of these
   by scale: at the pressures involved it is micron-order, comparable to the films themselves.
4. **Whatever drives the bush's stroke.** Section 3.1's failure is specifically an amplitude
   failure — the modelled bush barely moves over the revolution while the measured one sweeps
   ~9 um. Either the forcing that moves it (the gas-ramp variation and squeeze) is too weak or
   the seat film holds it too stiffly. This is the sharpest single lead left, and it may share
   a root with the sub-micron film both models predict.

## 5. What was NOT done

**No parameter in this repository was changed to close the gap.** The seat-friction term was
added to the *validation harness* only, measured to be ineffective, and not carried into
`mochi.rotor_bush_dynamics`. Chasing the measurement through the remaining free quantities —
flat clearance, seat clearance, effective pad length, `μ_seat`, blade tilt, runout residual:
six unknowns against two measured waveforms — would be curve fitting, not validation.

The honest status is: **our bush film module agrees with the published model of the same
sliding pair and disagrees with the published measurement of it, by the same factor that
model does.** That is a defensible position to hold, and the open items above are what would
have to be settled to move off it.

## 6. Data

`results/data/validation/reference_bush_film.dat` (Tecplot ASCII, one point zone over the
revolution) carries, per shaft angle: measured `h_upper`/`h_lower`, the best-fit model's
`h_upper`/`h_lower` (the section 3.1 winner, time-integrated), T02's simulated band, the tilt
from both, and the two residual columns. Plotting `residual_upper_um` alone shows the
one-signed bias directly.

`refs/` holds the two source papers, the page renders used for digitising, and the digitised
Fig.8 (`Fig8.csv` as exported by WebPlotDigitizer, `fig8_grid.npy` resampled onto a 5 deg
grid). The digitisation reproduces the four values J08 quotes in its text to +/-0.05 um on
three of them, and its recomputed tilt matches Fig.9 in both range and peak angle.
