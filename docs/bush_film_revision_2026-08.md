# Bush-film revision (2026-08) — cavitation, groove gas boundary, mixed lubrication

Record of a review of the swing-bush oil films, triggered by comparing the coupled 9-DOF
model against Pan et al. (2022) and ending in one genuine bug fix, two physics additions, and
several negative results worth not repeating.

**Status of the reference figure — Pan 2022 validated nothing.** It is worth saying plainly,
because the framing above ("against Pan et al. 2022") overstates its role. Pan is a *different
machine* (valveless swing compressor, swing rod on the rotor, bush in the cylinder, 10.34 cc)
at different pressures and with a different refrigerant model, and it reports a loss split
over eight channels of which this model covers three. Nothing in it can judge a term in our
model.

What it did was **prompt a question** — "is 91 W reasonable?" — and every finding below then
rested on internal evidence:

| finding | evidence | needed Pan? |
|---|---|---|
| cavitation bug | total **absolute** pressure at −11.6 MPa | no |
| curved-film gas BC missing | our own PHYSICS.md 3.3 geometry | no |
| Greenwood-Tripp vs the hard clamp | clamp sensitivity, section 5 | no |
| thin curved film is a load-path result | the attitude moment balance | no |
| journal runs full-film (Λ≈11) | film thickness vs roughness | no |

Even the prompt was not unique to Pan: "can the absolute pressure be −11.6 MPa" is a check
this repo could have run against itself at any time, and section 7 now does. The real
lesson of the episode is that the model lacked an internal physical-plausibility check, not
that it lacked an external reference.

The reference that *did* validate something is Tanaka (2002/2008) —
`docs/reference_bush_validation.md` — which supplies measured film thickness, passes a
quantitative kinematics check, and pins a quantitative disagreement.

---

## 1. What triggered the review

At the 3-option regime (`curved_shear_moment` + `flat_mixed_lubrication` +
`gas_film_boundary`) the bush loss came out at **91 W**, while the seal (11 W) and the
crank-pin journal (10 W) both sat within ~25 % of any reasonable anchor. One channel being
~7x out while its neighbours agreed is a model smell, not a machine difference.

## 2. The bug: cavitation was not applied when a gas boundary was present

`slider_film` / `arc_film` used to switch the Gumbel clamp **off** whenever end pressures were
imposed:

```python
pressure = solve_line_pressure(film, line_source, d_s, cavitation=not gas_present)
```

with the justification "high gas pressure floods the film and suppresses cavitation". That was
asserted, never checked. Evaluating the real operating state (film 0.84 um, tilt 1.1e-5 rad,
slide 0.86 m/s, ends 4.0 MPa -> 0.82 MPa) shows the **total absolute pressure reaching
-11.6 MPa** on the diverging half-stroke of the reciprocating pad — physically impossible, and
worth -2310 N of spurious suction that slammed the pad into contact.

**Fix.** Cavitation is a constraint on the *total* field, not on the hydrodynamic part alone.
Both film modules now solve the hydrodynamic problem unclamped, superpose the Poiseuille gas
bias (Reynolds is linear in `p`), and clamp the sum:

```python
pressure = solve_line_pressure(film, line_source, d_s, cavitation=False)
if gas_present:
    pressure = pressure + poiseuille_bias_pressure(...)
np.maximum(pressure, cavitation_pressure_pa, out=pressure)
```

`cavitation_pressure_pa` is expressed in the same gauge as the end pressures and defaults to
`0.0`, which reproduces classic Gumbel exactly — **the validated pure-hydrodynamic path is
unchanged** (its regression tests pass untouched).

## 3. Cavitation floor: consistency with the journal film

With a pressurised ambient the floor is not "zero gauge" and not "absolute zero" either. The
journal film in this repo already cavitates at its ambient (PHYSICS.md 4.10: symmetric 4.0 MPa
ends, `p_cav = 4.0 MPa`), which is the right physics for refrigerant-saturated oil — dissolved
refrigerant evolves as soon as the pressure falls below the local saturation, long before
absolute zero.

The bush films have **asymmetric** ends (4.0 MPa recess <-> chamber), so their floor cannot be
the ambient (that would clip the real throughflow ramp). The rule adopted is

> a film cannot fall below the lowest pressure it is connected to

i.e. `cavitation_pressure_pa = min(0, chamber - bore)`. For the journal (both ends at ambient)
this rule degenerates to the existing `p_cav = ambient`, so the three films are now consistent.

**This floor is still the main open uncertainty** — the true value lies between the chamber
pressure and absolute zero and depends on the oil/refrigerant solubility, which is not
modelled. The bush friction is very sensitive to it (see the table in section 6).

## 4. Curved film: missing gas boundary and the wrong land

`gas_film_boundary` imposed the bore->chamber drop on the flat (vane) film only; the curved
(rotor-groove) arc got nothing, and the film was integrated over the whole piece arc
(+/-53.05 deg). PHYSICS.md 3.3 says otherwise:

* the groove wall runs from the mouth blend to the channel opening, giving an **86.5 deg
  sealing land** — the rest of the arc faces gas, not oil;
* one end of that land opens to the rotor **mouth** (the piece's own chamber), the other to
  the **recess channel** at `RECESS_PRESSURE_PA = 4.0 MPa` (the bore reference).

So the curved film carries a genuine pressure drop, mirrored between the two pieces. Both are
now modelled: the film is integrated over the land (offset centre + 43.25 deg half-span,
reproducing 41.5 + 45.0 = 86.5 deg exactly) with the corresponding end pressures.

The resulting curved gas force is **+367 N per piece away from the vane**, largely cancelling
the flat film's -400 N squeeze — the counteracting force that was missing. This closes the
PHYSICS.md open item "the bush-face recess BCs remain to be confirmed" for the film ends
(the recess pressure value itself was already documented).

## 5. Mixed lubrication: Greenwood-Tripp replaces an arbitrary pressure scale

The flat film's asperity contact used `p_asp = P0 exp(-h/sigma)` with `P0 = 50 MPa` chosen to
put the film in the mixed regime — a model-sensitive knob the friction scaled linearly with.
It is replaced by the standard Greenwood-Tripp elastic contact
(`mochi.asperity_contact`), driven by measurable surface parameters:

```
p_asp(h) = (16 sqrt(2)/15) pi (eta beta sigma)^2 E' sqrt(sigma/beta) F_{5/2}(h/sigma)
```

with literature-typical ground-steel values (`beta = 10 um`, `eta = 1e10 /m^2`,
`E' = 115 GPa`, giving `eta beta sigma ~ 0.05`, inside the usual 0.02-0.06 band). The same
contact is now applied to the **curved** film as well, sharing load with the hydrodynamic
film instead of relying on the hard eccentricity clamp.

Isolation run: switching the curved asperity off makes things *worse* (bush 48.7 -> 59.1 W at
1800 rpm), because the piece then runs into the hard clamp at `ecc/gap = 0.99`, mispositions,
and loads the flat film harder. The curved asperity stays on.

**Why the contact model was needed at all** (measured on our machine, no external reference).
There is no measured film for this compressor, so "accuracy" cannot be scored. What can be
measured is how much of the answer a *numerical* constant decides. `_MAX_ARC_ECC_RATIO` is a
solver guard with no physical content, and the curved film sits against it. Running the same
3-option configuration with the contact model off, and moving that guard over its range:

| | bush W | min curved film | min flat film | ecc ratio |
|---|---|---|---|---|
| **Greenwood-Tripp contact** | 6.15 | **0.583 µm** | **0.890 µm** | **1.33** |
| clamp only, 0.99 | 4.90 | 0.300 µm | 0.200 µm | **82.8** |
| clamp only, 0.95 | 4.63 | 1.500 µm | 0.200 µm | **100.7** |
| clamp only, 0.999 | 8.30 | 0.030 µm | 0.200 µm | **79.5** |

Moving an arbitrary constant from 0.95 to 0.999 swings the bush loss by **79 %** and the
minimum curved film by **4900 %** (0.030 → 1.500 µm) — and the minimum film is what decides
wear margin. Worse, with no contact reaction the eccentricity ratio runs to **80–100**: the
piece drives a hundred clearances into the groove, and the clamp only hides it because it is
applied to what the film solver is handed rather than to the state. The flat film meanwhile
pins at exactly 0.200 µm in all three — its `_MIN_FLAT_FILM_FRACTION` floor, not a result.

With Greenwood-Tripp the contact carries a real reaction, the eccentricity settles at 1.33,
the flat film leaves its floor, and none of it depends on the guard. GT has its own spread
(±16 % over the literature `beta` range) but those parameters are measurable with a
profilometer, which the clamp value never was. The friction *totals* are similar either way
(4.6–8.3 vs 6.15 W); the difference is in the kinematics, which is what wear and attitude
predictions rest on.

`beta` sensitivity (5/10/20/40 um) moves the bush loss over 66-91 W non-monotonically — the
remaining parameter uncertainty, to be pinned by a surface measurement.

## 6. Effect of the changes (3-option regime, our machine, our pressures)

| stage | bush @1800 | bush @5400 | F_asp @1800 | min flat film |
|---|---|---|---|---|
| as found (cavitation bug) | 91.3 W | 462.9 W | 1654 N | 0.839 um |
| + cavitation on the total | 68.0 W | 277.0 W | 1226 N | 0.867 um |
| + curved gas boundary & land | 48.7 W | 247.9 W | 732 N | 0.934 um |
| + consistent cavitation floor | 6.2 W | 33.5 W | 24.8 N | 0.890 um |

The last row also changes the *character* of the solution: the boundary/viscous split goes
from 95:5 to 49:51, and the flat film **thickens with speed** (0.89 -> 1.18 um) instead of
thinning — the first physically healthy mixed-lubrication behaviour in this channel.

For scale only (not a target): Pan's swing-rod-bush loss scaled by displacement is ~13.6 W at
1800 rpm, ~19.7 W if also scaled by the pressure difference. Our result brackets that band
from above (48.7 W) before the floor change and from below (6.2 W) after it, which is the
honest statement of where the model stands.

## 6.1 The check that should have caught it — `mochi.physical_checks`

The cavitation bug ran undetected for as long as the gas boundary existed, and surfaced only
because an external comparison prompted someone to look. Yet the evidence was internal the
whole time: a field at −11.6 MPa absolute is impossible on its own terms. A solver returning
such a field should say so, so now it can.

`mochi.physical_checks` holds assertions over quantities whose valid range comes from physics
rather than from tuning: absolute pressure against the lubricant's vapour pressure (referenced
to the immersion, so a −1 MPa *gauge* field is fine at 4 MPa immersion and impossible at
absolute), film thickness against interpenetration, eccentricity ratio against the clearance,
and dissipated power against sign. `flat_slider_film` and `arc_film_force` take
`check_physical=True` plus `reference_pressure_pa`; they default off because a coupled orbit
evaluates these films tens of thousands of times, and are meant to be switched on when a
configuration is new or a number is surprising.

`tests/test_physical_checks.py` reconstructs the 2026-08 state and asserts it is rejected
(`absolute pressure reaches -11.555 MPa`), that the current cavitation floor passes it, and
that the validated pure-hydrodynamic path does not false-positive.

## 7. Negative results (do not repeat)

* **Bore-pressure sweep** (4.0 -> 3.0 MPa): bush loss only 91 -> 81 W and the speed exponent
  got *worse*. The gas load was never the dominant term; the spurious suction was.
* **mu_b calibration**: the boundary coefficient needed to hit a target differs by 9x between
  1800 and 5400 rpm. A material constant cannot do that — the discrepancy was structural.
* **Piezoviscous EHL**: with `alpha ~ 6-13 /GPa` and film pressures of 10-30 MPa, `alpha p`
  is 0.06-0.4, i.e. a 6-50 % viscosity rise. Far too small for the magnitudes in question.
* **Elastic (Winkler) compliance on the curved film**: implemented and measured, then
  reverted. It did **not** move the kinematics (`ecc/gap` 1.328 -> 1.330, `dphi_r` 1.006e-3 ->
  1.027e-3 rad) because the curved eccentricity is set by *geometry* — the piece is pinned to
  the vane by the flat film while the groove rides the rotor — not by film compliance. It
  thinned the curved film (0.583 -> 0.348 um), raised the 5400 rpm loss (33 -> 53 W), and cost
  4x runtime.

## 8. What the review established about the design (not model errors)

* At the 0 deg reference `piece_shift = 20 um` gives **both** films exactly 10 um
  (`3.990 + 0.020 - 4.000` and `8.000 - 7.970 - 0.020`) — a balanced design.
* The rotor attitude EOM closes: `M_gas +22.2` vs `M_bush_curved -21.3 N*m` (96 % cancelled),
  residual 0.019 N*m = 2 % of the reference inertia torque. The bush curved film is the
  structural load path that carries the gas moment.
* Carrying that moment needs `21.28 / 0.025 = 851 N` of curved-film reaction. A direct
  load-capacity sweep of the arc film gives 725 N at `ecc/gap = 0.98` and 1527 N at 0.99, so
  **851 N lands at `ecc/gap ~ 0.985`, i.e. a 0.45 um film** — which is exactly what the
  coupled orbit reports (0.58 um). The thin curved film is a consequence of the load path,
  not a numerical artefact, and it is why the eccentricity sits against the clamp.
* Consequently the curved film runs at `Lambda = h/sigma ~ 1` (boundary regime) by design,
  while the crank-pin journal runs at `Lambda ~ 11` (full film) and needs no contact model.

## 9. Open items

1. The cavitation floor (section 3) — needs an oil/refrigerant solubility model to pin down.
2. `beta`, `eta` for the Greenwood-Tripp contact — need a surface measurement.
3. Whether the curved film really runs at `Lambda ~ 1` in the hardware (wear evidence would
   settle it); the squeeze term is ~16x stronger than the wedge at the same eccentricity, so a
   real machine may ride on squeeze more than this steady-orbit model does.
