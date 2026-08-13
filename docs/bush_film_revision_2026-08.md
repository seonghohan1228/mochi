# Bush-film revision (2026-08) — cavitation, groove gas boundary, mixed lubrication

Record of a review of the swing-bush oil films, triggered by comparing the coupled 9-DOF
model against Pan et al. (2022) and ending in one genuine bug fix, two physics additions, and
several negative results worth not repeating.

**Status of the reference figure**: Pan 2022 is a *different machine* (valveless swing
compressor, swing rod on the rotor, bush in the cylinder, 10.34 cc) at different pressures
(0.867/2.542 MPa vs our 0.820/3.240 MPa) and a different refrigerant model. Its numbers are
used here **only as an order-of-magnitude anchor and a plausibility band** — never as a target
to tune to. Mechanism differences make absolute mechanical efficiency unpredictable across
the two machines, so no parameter in this repo was fitted to reproduce a Pan value.

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
