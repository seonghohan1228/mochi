# Physics and numerical methods

This document is the contract between the physical model and its
implementation. It starts deliberately conservative: unconfirmed compressor
geometry and force laws are listed as open decisions rather than presented as
facts.

## 1. Problem statement and initial scope

The project will become a broader rotary-compressor simulation program. The
current accepted implementation is only a supporting prescribed-motion
visualization for checking mechanism geometry and animation conventions.

Currently in scope:

- a fixed cylinder with an integral fixed-length stepped vane at the top
  center;
- a prescribed eccentric orbit for the rotor center;
- a rotor-fixed circular cutout constrained to the vane centerline;
- a white clearance slot slightly wider than the vane;
- physical crank speed separated from slowed display speed;
- suction/discharge chamber cross-section areas per crank angle from the
  circular-rotor approximation (Section 3.1), a geometric quantity only.

Not implemented by this visualization:

- pressure, flow, force, contact, friction, torque, or deformation;
- initial-value dynamics or numerical time integration;
- bearing, lubrication, thermal, or fluid-structure models;
- the main rotary-compressor solver and its production outputs.

These belong to the broader simulation only after their equations and
validation cases are reviewed.

## 2. Coordinates, signs, and units

All internal quantities use SI units: metre, kilogram, second, newton, pascal,
joule, kelvin, and radian. Input/output adapters may display other units only
through explicit conversion.

The visualization frame is fixed to the cylinder and defined in Section 3;
Section 4 promotes it to the global inertial frame and records the accepted
bodies, frames, degrees of freedom, and sign discipline. No equation of
motion is evaluated yet; the force-on-body convention for each new model
term follows Section 5.

Angles are in radians internally. Each angular input must state its zero
direction and positive sense. Pressure models must state whether pressure is
absolute or gauge and which reference pressure is used.

## 3. Prescribed rotary mechanism

The test GUI uses a fixed cylinder-centered frame. The origin is the cylinder
center, positive $y$ points toward the top vane, and positive $x$ points toward
the inlet side. The crank angle $\theta=0$ places the rotor at the top of its
orbit, $(x_r,y_r)=(0,e)$; increasing $\theta$ is clockwise.

**Top dead center (TDC)** is $\theta=0$: the rotor-cylinder tangency point
coincides with the vane position, the discharge volume reaches zero, and the
suction charge seals over (Section 3.1). **Bottom dead center (BDC)** is
$\theta=\pi$: the tangency point lies diametrically opposite the vane and the
two chambers have equal volume. BDC is a geometric reference only — the
chamber volume change rate is maximal there, and each chamber's volume
extreme occurs at the vane position ($\theta\to2\pi^-$), not at $\theta=\pi$.

The supplied geometry is:

| Symbol | Meaning | Value |
|---|---|---:|
| $D_c$ | cylinder inside diameter | 77.0 mm |
| $D_r$ | main rotor outside diameter | 68.0 mm |
| $e$ | rotor-center eccentricity | 4.5 mm |
| $r_h$ | circular-cutout radius | 8.0 mm |
| $L$ | rotor-center to cutout-center distance | 25.0 mm |
| $w_v$ | vane width | 8.0 mm |
| $a_v$ | rotor-center to vane-tip distance at the top position | 9.0 mm |
| $H_c$ | cylinder axial height | 21.0 mm |
| $f$ | physical rotation frequency | 30 Hz = 1800 rpm |

The eccentricity is the radial difference,

$$
e = \frac{D_c-D_r}{2} = 4.5\ {\rm mm},
$$

so the main rotor remains internally tangent to the cylinder in this ideal
geometry. Its prescribed center is

$$
x_r=e\sin\theta,\qquad y_r=e\cos\theta.
$$

The circular-cutout center $\mathbf{h}=(0,y_h)$ is fixed in the rotor a distance
$L$ from the rotor center and constrained to the stationary vane centerline. Selecting
the upper assembly branch gives

$$
y_h=y_r+\sqrt{L^2-x_r^2}.
$$

The rotor orientation is the direction from the rotor center to the cutout center,

$$
\phi=\operatorname{atan2}(y_h-y_r,-x_r).
$$

The cutout fits inside the 34 mm rotor radius because $L+r_h=33$ mm. A radial
slot opens the circular cutout through the rotor OD and stops at the circular
boundary; no rectangular or square cutout continues toward the rotor center.
Its default displayed width is $w_s=w_v+1$ mm. At the two rotor lips, tangent
fillets with nominal radius $r_f=1.5$ mm round the transition from the rotor OD
to the slot. The displayed fillet radius is reduced automatically only if
edited inputs leave insufficient material.

Let $\mathbf{r}=(x_r,y_r)$ and define the rotor-fixed axial and transverse unit
vectors

$$
\mathbf{u}=\frac{\mathbf{h}-\mathbf{r}}{L},\qquad
\mathbf{n}=(-u_y,u_x).
$$

The stepped vane is integral with the cylinder (Section 3.3), so its tip is
fixed on the vane axis:

$$
y_{tip}=e+a_v=13.5\ {\rm mm},\qquad
\ell_v=\frac{D_c}{2}-y_{tip}=25.0\ {\rm mm},
$$

where the supplied $a_v$ is the rotor-center-to-tip distance at top dead
center. The existing parameter constraints ($a_v<L\le R_r-r_h$ and
$e\le R_c-R_r$) guarantee the fixed tip always lies inside the bore. In the
planar drawing the vane outline may overlap the rotor silhouette below the
circular cutout; the real parts do not interfere because of the axial
stepped structure of Section 3.3.

The light-gray rotor is one parametric polygon built from the rotor OD, two
fillet arcs, the slot sides, and the remaining circular-cutout circumference.
The majority cylinder arc and the right, bottom, and left vane sides are one
continuous line with no fill. Its interior is transparent, so the vane can
overlap the rotor below the circular cutout; the surrounding canvas remains
white. No white shapes are overpainted to manufacture a cutout.

At constant physical frequency,

$$
\theta(t)=\theta_0+2\pi f t.
$$

The display slow factor $s$ changes wall-clock animation only:
$\theta_{display}(t_{wall})=\theta_0+2\pi f s t_{wall}$. It must never alter
physical time, velocity, force, or result metadata.

The cylinder center `C`, rotor center `R`, and dashed circular locus of radius
$e$ are display guides only; they do not add geometry or a physical boundary.

Ports are angular windows on the bore, not single markers; their supplied
angles and the cycle events they define are in Section 3.4. Angles are
measured from positive $y$ toward positive $x$, the same sense as the crank
angle. The earlier symmetric $\pm 30^\circ$ display markers are superseded
and the real windows are not symmetric.

This prototype is prescribed kinematics. The cutout centerline constraint,
ideal tangency, and displayed slot clearance are not contact-force or
lubrication models.

### 3.1 Chamber cross-sections from the circular-rotor approximation

As groundwork for later suction/discharge pressure boundary conditions, the
prescribed mechanism reports the suction and discharge chamber cross-section
areas at each crank angle. This is a purely geometric evaluation; it adds no
pressure, flow, or force model, and the drawn geometry is unchanged.

Approximation: the rotor is replaced by its full outside-diameter disc of
radius $R_r=D_r/2$ centered at $\mathbf{r}=(x_r,y_r)$. The rotor slot and the
circular cutout are ignored. The working region is the crescent between the
cylinder bore (radius $R_c=D_c/2$) and this disc, with total area

$$
A_{cres} = \pi\left(R_c^2-R_r^2\right).
$$

The union of the vane and the rotor disc divides the crescent into two
regions: the suction chamber on the inlet side ($x>0$, marker `IN`) and the
discharge chamber on the outlet side ($x<0$, marker `OUT`). The suction
chamber runs clockwise from the vane's inlet-side face to the rotor-cylinder
tangency point; the discharge chamber continues clockwise back to the vane's
outlet-side face.

Along a ray from the cylinder center at clockwise angle $\varphi$ from
positive $y$, the distance to the rotor-disc boundary is

$$
\rho(\varphi) = e\cos(\varphi-\theta)
  + \sqrt{R_r^2-e^2\sin^2(\varphi-\theta)},
$$

valid because $e<R_r$ keeps the cylinder center inside the disc. The crescent
area swept from the vane centerline to the tangency point has the closed form
$A_{sweep}(\theta)=F(\theta)$ with

$$
F(\psi) = \tfrac12\left(R_c^2-R_r^2\right)\psi
  - \tfrac14 e^2\sin 2\psi
  - \tfrac12\, e\sin\psi\sqrt{R_r^2-e^2\sin^2\psi}
  - \tfrac12 R_r^2\arcsin\!\left(\frac{e\sin\psi}{R_r}\right).
$$

The vane occupies the part of the crescent above the rotor disc within
$|x|\le w_v/2$. With $G(x;R)=\tfrac{x}{2}\sqrt{R^2-x^2}
+\tfrac{R^2}{2}\arcsin(x/R)$, each half-width strip is

$$
A_{strip}^{\pm} = \left|G(\pm\tfrac{w_v}{2};R_c)\right|
  - \tfrac{w_v}{2}\,y_r
  - \left|G(\pm\tfrac{w_v}{2}-x_r;R_r) - G(-x_r;R_r)\right|,
$$

and with $\theta$ normalized to $[0,2\pi)$ the chamber areas are

$$
A_{suc}(\theta) = F(\theta) - A_{strip}^{+},\qquad
A_{dis}(\theta) = A_{cres} - F(\theta) - A_{strip}^{-}.
$$

A chamber volume is the area multiplied by the axial cylinder height,
$V = A\,H_c$ with the supplied $H_c=21$ mm; a different height may be passed
explicitly. The maximum chamber volume is therefore
$A_{cres}H_c \approx 21.5\ {\rm cm^3}$ minus the vane strips.

Validity conditions. The evaluation refuses to produce a value (it raises a
specific error) instead of forcing a split when the two chambers do not
exist:

1. **Tangency seal:** $e = R_c-R_r$ within 1 nm. A larger radial gap connects
   the chambers around the rotor.
2. **Vane seal:** both vane bottom corners $(\pm w_v/2,\,y_{tip})$ must lie
   strictly inside the rotor disc so the vane bridges the cylinder wall to
   the disc. For the supplied geometry the smallest margin over a revolution
   is about 22 mm.
3. **Seal-over window:** for $|\theta| < \theta_s$ mod $2\pi$, with
   $\theta_s=\arcsin\!\left(\frac{w_v}{2R_c}\right)\approx 5.96^\circ$ for the
   supplied geometry, the tangency point lies under the vane width and the
   two chambers form one connected region. This is the physical instant where
   the completed suction chamber seals over and becomes the next compression
   chamber. Inside this window no split is reported.

At $\theta=\theta_s$ the suction area is exactly zero and grows
monotonically to nearly the full crescent at $\theta=2\pi-\theta_s$;
mirrored crank angles swap the suction and discharge areas. The closed forms
are regression-tested against an independent Simpson quadrature of the same
regions and against conservation,
$A_{suc}+A_{dis}+A_{strip}^{+}+A_{strip}^{-}=A_{cres}$.

**Recommendation:** keep the eccentricity locked to the tangency value
$e=(D_c-D_r)/2$ (the GUI `Lock e` checkbox). Any smaller eccentricity opens
a radial rotor-cylinder gap, the chambers stay connected at every crank
angle, and no chamber split, volume, or pressure assignment is defined.

### 3.2 Accepted chamber pressure conditions

All pressures are **absolute**. The working fluid is **R410A**. The supplied
port pressures correspond to standard rating saturation temperatures
(0.85 °C evaporating, 52.5 °C condensing; CoolProp R410A HEOS):

| Symbol | Meaning | Value |
|---|---|---:|
| $p_{suc}$ | suction-port pressure (absolute) | 0.82 MPa |
| $p_{dis}$ | discharge-port pressure (absolute) | 3.24 MPa |
| $n$ | effective R410A polytropic exponent | 1.07 |
| $f_v$ | discharge-valve opening rise, fraction of $p_{dis}$ | 0.05 |
| $p_{open}$ | valve opening pressure $(1+f_v)\,p_{dis}$ | 3.40 MPa |

The exponent is the effective value fitting $p\,v^n=\mathrm{const}$ between
the endpoints of the real-gas isentrope from saturated R410A vapour at
0.82 MPa to 3.24 MPa (CoolProp HEOS: 1.064 for saturated suction, 1.075 at
10 K superheat; note the ratio $c_p/c_v\approx1.38$ near suction is *not*
the correct $p$–$v$ exponent for this real gas). The valve rise
approximates the reed-valve opening loss of small hermetic compressors
(typical 2–6 % of discharge pressure); it awaits confirmation against
actual valve data. Both are case inputs with these defaults.

The accepted per-chamber rule is:

- **Suction chamber:** always at the suction-port pressure $p_{suc}$
  (ideal filling through the inlet, no pressure drop).
- **Compression (discharge) chamber:** seals at $\theta=\theta_s$ with
  volume $V_0=V_{dis}(\theta_s)$ at pressure $p_{suc}$ and follows the
  polytropic relation until it reaches the valve opening pressure
  $p_{open}=(1+f_v)\,p_{dis}$ at the crank angle $\theta_v$; during
  delivery the pressure falls **linearly in crank angle** from $p_{open}$
  to the discharge-port pressure at the seal-over window entry,
  approximating the decaying valve/flow overpressure:

$$
p_{comp}(\theta)=
\begin{cases}
p_{suc}\left(\dfrac{V_0}{V_{dis}(\theta)}\right)^{\!n}, &
\theta_s\le\theta\le\theta_v,\\[1.2em]
p_{open}+\left(p_{dis}-p_{open}\right)
\dfrac{\theta-\theta_v}{(2\pi-\theta_s)-\theta_v}, &
\theta_v\le\theta\le2\pi-\theta_s.
\end{cases}
$$

  Because the axial height is constant, the volume ratio equals the
  cross-section area ratio of Section 3.1. With the accepted values the
  valve opens at $V_0/V_{dis}\approx(p_{open}/p_{suc})^{1/n}\approx3.8$,
  near crank angle $\theta_v\approx221^\circ$.
- **Seal-over window:** in the entering half
  ($2\pi-\theta_s<\theta<2\pi$) the merged region follows a linear mixing
  ramp from $p_{dis}$ down to $p_{suc}$, reaching the suction pressure
  exactly at $360^\circ$; in the leaving half ($0\le\theta<\theta_s$) the
  region stays at $p_{suc}$.

With these ramps the discharge-chamber pressure trace is continuous at
every transition — compression start ($p_{suc}$ at $\theta_s$), valve
opening ($p_{open}$ at $\theta_v$), window entry ($p_{dis}$ at
$2\pi-\theta_s$), and top dead center ($p_{suc}$ at $2\pi$). The delivery
decline and the mixing ramp are prescribed boundary-condition
simplifications chosen for continuity, not valve-dynamics or flow models.

Valve dynamics, discharge-port timing, internal leakage, heat transfer, and
clearance re-expansion are neglected. This is a pressure
boundary-condition rule for later gas-force models, not a full
thermodynamic compression model; pressure forces on the rotor and vane are
not yet implemented.

### 3.3 Axial structure: stepped vane, rotor recess channel, swing bush

The vane is one piece with the cylinder and stepped in the axial direction;
the rotor carries a matching recessed channel. Supplied dimensions
(2026-07-20):

| Item | Value |
|---|---:|
| vane full-thickness segment, depth from the bore | 0 to 15.4 mm (21 mm thick) |
| vane ledge segment, depth from the bore | 15.4 to 25.0 mm |
| ledge thickness (one at the top, one at the bottom) | 2.4 mm each |
| open axial gap between the ledges | 21 − 2(2.4) = 16.2 mm |
| vane tip depth from the bore (global $y_{tip}$) | 25.0 mm ($y=13.5$ mm) |
| vane-root blend into the bore (both flanks) | R2.1 mm (supplied 2026-07-23) |
| vane-tip corner round (both corners) | R1.5 mm (leaves a 5.0 mm tip flat) |
| rotor base thickness | 21.0 mm |
| rotor recess-channel thickness | 14.7 mm |
| recess depth per side | (21 − 14.7)/2 = 3.15 mm |
| recess-channel width, centered on the circular groove | $8\sqrt{2}\approx11.31$ mm |
| recess-channel extent | rotor center to the circular groove, joining the groove arc without interruption |
| ledge-to-recess axial clearance | 3.15 − 2.4 = 0.75 mm per side |
| rotor OD flat (**inlet side only**), at the rotor center | starts at 34.8° on the OD, ends at 13.4° at radius 33.657 mm |
| OD flat length / maximum OD reduction | 12.566 mm / 0.773 mm |
| outlet-side mouth start | lip tangent to the OD at 13° (original design, no flat) |
| mouth radial extent | groove circle to the rotor OD, through the full axial height |
| lip radii (both sides) | R1.5 leaving the OD, R1.0 entering the groove |
| inlet lip (tangent-continuous) | flat end → R1.5 (102.2°) → 0.399 mm straight → R1.0 (52.62°) → groove at 47.74° |
| outlet lip (tangent-continuous) | OD tangent → R1.5 (92.6°) → 0.744 mm straight → R1.0 (52.43°) → groove at 47.97° |
| recess-space pressure (absolute, fixed for now) | 4.0 MPa |

**The mouth is asymmetric.** The inlet side is anchored by the OD flat
(34.8° on the OD to 13.4° at radius 33.657 mm) and a 102.2° first sweep;
the outlet side keeps the original design, with the lip tangent to the OD
at 13° and a 92.6° first sweep. Both use R1.5 leaving the OD and R1.0
entering the groove (the earlier 1.4 / 0.9 ignored the chamfer). Solving
each tangent-continuity (G1) chain then fixes the rest: the inlet straight
is 0.399 mm with a 52.62° blend reaching the groove at 47.74°, and the
outlet straight is 0.744 mm with a 52.43° blend reaching it at 47.97°.
Both blend sweeps sit at the 52° shown on the supplied CAD view, and both
groove positions stay near the 48.5° of the previous design. The
groove-side half-angle still clears the vane at the extreme relative tilt:
the vane flank subtends about 34.2° at the groove center and the rotor
tilts up to 10.4° relative to the vane, totalling 44.6° with a 3.1° margin
on the tighter (inlet) side. The narrow drawn slot of the current GUI
contour ($w_s=w_v+1$ mm with lip fillets, Section 3) is a display
simplification of this mouth and its update is planned.

**Bore clearance:** the profile's maximum radius is exactly the rotor OD —
reached where the inlet flat begins at 34.8° and along the outlet lip's
tangent point — while everything else lies inside that envelope (the inlet
lip peaks at 33.68 mm). A prescribed-motion sweep at 0.2° crank steps
confirms zero bore penetration at every crank angle. This supersedes the
earlier semicircular-lip variants, which would have penetrated the bore by
up to 0.32 mm near top dead center.

**Vane-root blend (supplied 2026-07-23):** each vane flank meets the bore
through an R2.1 mm blend rather than a sharp corner, a stress-relief fillet
at the root of the fixed vane. The blend is tangent to the flank at
$x=\pm w_v/2$ and internally tangent to the bore, so its centre sits one
blend radius outside the flank and one inside the bore; for the supplied
$w_v=8$ mm and $R_c=38.5$ mm the flank tangent is at
$(4.000, 35.885)$ mm and the bore tangent at $(6.452, 37.956)$ mm on the
$+x$ side, mirrored on the $-x$ side. The blend is full axial height, so on
each side it converts a small full-height pocket between the flank, the
bore, and the blend arc from gas into cylinder material. Its only effect on
the model is that volume: it removes about 0.028 cm³ from the chamber gas
(see Section 3.4), lowering the clearance volume from 0.193 to 0.165 cm³,
still inside the 0.1–0.5 cm³ band. It does not touch the seal-over geometry
or the rotor.

**Vane-tip round:** each corner of the fixed vane tip carries an $R=1.5$ mm
round — the flank ($x=\pm w_v/2$) runs straight down to $y=y_{tip}+R=15.0$ mm,
then a quarter-circle to the tip flat at $y_{tip}=13.5$ mm, tangent there at
$x=\pm(w_v/2-R)=\pm 2.5$ mm, leaving a $w_v-2R=5.0$ mm tip flat. The tip sits at
$y=13.5$ mm, deep inside the rotor disc and **below the gas crescent**, so it
removes no chamber-gas volume. Its one modelled effect is on the **swing-bush
flat film** (Section 3.6): near BDC the bush flat overhangs the round, shortening
the flat contact from $11.94$ to $11.47$ mm — but that is where the flat sliding
velocity is near zero (the translation turning point), so the flat friction is
essentially unchanged. The round is carried in `kinematics.vane_tip_round`.

The apparent vane-rotor overlap in the planar drawing is resolved axially:
the full-thickness vane segment always stays within the circular-groove
span, and the thin ledges travel inside the rotor's recessed channel. The
recess spaces above and below the thin rotor section are sealed by the
full-height rotor, the swing bush, and the vane; their pressure is fixed at
4.0 MPa absolute (`RECESS_PRESSURE_PA`, not yet used in any computation)
until the swing-bush pressure boundary conditions are examined. Oil-film
leakage across these seals is acknowledged and neglected until then.

**Swing-bush placement and clearances (supplied 2026-07-20):** the swing
bush consists of two half-moon pieces riding in the $\varnothing16$ groove
with their central slot aligned to the fixed vane. Each piece runs flat face
→ R0.5 fillet → 7.970 mm cylindrical face (106.1°) → R0.5 fillet → flat
face, with the flat 3.990 mm from the piece's own center. The fillet radius
sets the arc extent exactly:

$$
\tfrac12\,\text{arc}=\arccos\frac{3.990+0.5}{7.970-0.5}=53.05^\circ,
$$

reproducing the supplied 106.1°. A piece centered on the groove would
overlap the vane by 0.010 mm; shifting each piece 0.020 mm toward its own
side (inlet or outlet) opens both faces to a uniform 0.010 mm:

| Item | Value |
|---|---:|
| piece outer radius | 7.970 mm (Ø15.94) |
| cylindrical face angular extent | 106.1° (±53.05° about the piece center) |
| flat face from the piece's own center | 3.990 mm |
| fillet joining the flat to the cylindrical face | R0.5 mm, tangent to both |
| fillet center (from the piece center) | (4.490, ±5.970) mm |
| flat contact length | 11.940 mm (\|y\| ≤ 5.970) |
| cylindrical face arc length | 14.76 mm |
| piece-center offset toward its own side | 0.020 mm |
| flat-face-to-vane film (operating position) | 3.990 + 0.020 − 4 = 0.010 mm |
| curved-face-to-groove film (operating position, minimum) | 8 − 7.970 − 0.020 = 0.010 mm |
| in-plane free play per piece (flat contact to groove contact) | 0.020 mm |
| axial clearance to each end plate | 0.0085 mm |
| bush height | 21 − 2(0.0085) = 20.983 mm |

The 0.020 mm offset is exactly the midpoint of the free-play interval
(flat contact at 0.010 mm, groove contact at 0.030 mm), which is why both
films equalize at 0.010 mm — the supplied numbers are mutually consistent.
The rotor's ±10.37° tilt does not change these films
to first order, because the groove is axisymmetric about its center and
the flat face merely slides along the vane. The bush translates along the
vane with the groove center (9.0 mm stroke per revolution, up to 0.86 m/s
at 30 Hz) while the rotor oscillates ±10.37° about it (up to 33.9 rad/s).
Derived geometric conditions:

- **Capture:** a piece is 7.970 − 3.990 = 3.980 mm thick radially, while the
  mouth opening measured on the groove circle — from the vane face to the
  lip's groove blend at 48.5° — is only 2.573 mm, so a piece cannot pass
  through the mouth. Angularly, the piece's cylindrical face reaches 37.07°
  from the groove center (the fillet starts there and the flat ends at
  33.89°) against the groove wall starting at 48.5°, leaving 1.07° of
  overlap even at the extreme 10.37° tilt.
- **Lip clearance:** the mouth lip stays outside the groove circle, so the
  bush-to-lip clearance equals the bush-to-groove radial clearance; no
  separate condition arises.
- **Sealing land:** the groove wall runs from the mouth's blend at 48.5° to
  the channel opening at 135°, so the film land against full-height rotor
  material is 86.5° wide. The piece's cylindrical face (37.07°–142.93°)
  covers that land throughout the cycle except near the extreme positive
  tilt, where the land's far end runs onto the piece's fillet and the
  effective land narrows to 84.07° (at θ = 90°, tilt 10.37°).
- **Sealing height:** the recess channel meets the groove over ±45° of its
  circumference, where the 3.15 mm recess bands (4.0 MPa absolute) end on
  the bush OD. The bush must fill the full 21 mm height (minus its axial
  clearance) or the 4 MPa space would leak past it; this contact presses
  the bush toward the vane.
- **Vane engagement:** the full-thickness vane engages the bush slot over
  14.4 mm at top dead center, shrinking to 5.4 mm at bottom dead center
  (in-plane), with the thin ledges passing through the same slot on the two
  2.4 mm bands.

With the clearances supplied, the remaining bush-related open item is its
pressure boundary conditions (Section 9). The planar chamber model of
Section 3.1 keeps using the rotor OD silhouette, and the vane bottom it
references is the ledge tip; the free gas space inside the mouth opening is
a neglected volume of that approximation. Section 3.4 now measures that
space directly (24.0 mm² at top dead center, of the same order as the 30 mm²
estimate recorded here) and shows it is the entire clearance volume of the
cycle.

### 3.4 Port timing and true-geometry chamber volume

Ports were display markers until this section: one angle each, drawn on the
bore and read by nothing. They are now angular windows, and the four
characteristic angles they define drive the phases of the compression cycle.

**Supplied port timing (2026-07-20).** All measured from the vane
centerline in the crank-angle sense, positive $y$ toward positive $x$:

| Symbol | Meaning | Value |
|---|---|---:|
| $\varphi$ | suction seal angle; the suction port opens here | 10.4° |
| $\beta$ | compression start; the chamber loses the suction port | 27.7° |
| $\gamma$ | discharge port angular span | 7.2° |
| $\delta$ | recompression angle; the discharge port shuts this far before top | 13.2° |

The cycle events follow directly: suction opens at $\varphi$, compression
starts at $\beta$, the discharge port spans
$[2\pi-\delta-\gamma,\; 2\pi-\delta] = [339.6^\circ,\; 346.8^\circ]$, and
recompression runs from $2\pi-\delta$ to the seal-over window entry.

$\gamma$ checks out against the hardware: on the 38.5 mm bore radius it
subtends a 4.84 mm arc, a plausible discharge port diameter for a machine of
this size. $\varphi$ is the suction-port opening angle — the crank position
at which the bore first uncovers the suction port — and it is a **machined
port location, not a vane-derived quantity**, so it is not expected to equal
any vane angle. In particular it is a different kind of thing from the
seal-over half angle $\alpha=5.96^\circ$ (full band $2\arcsin(4/38.5)
=11.93^\circ$) of Section 3.1, which the vane width does set: $\alpha$ is the
crank range where the rotor-cylinder contact hides under the vane and the two
chambers merge, whereas $\varphi$ is where the wall opens to the suction
port. The two measure unrelated features and no design rule forces them
equal, so the fact that $10.4^\circ \ne 11.93^\circ$ is expected, not a
discrepancy. $\varphi$ is carried as an independent supplied angle. The
suction port is the window $[\varphi, \beta] = [10.4^\circ, 27.7^\circ]$ and
the discharge port the window $[2\pi-\delta-\gamma, 2\pi-\delta] =
[339.6^\circ, 346.8^\circ]$; the two sit on opposite sides of top dead
center, suction just after it and discharge just before it.

**Discharge port open area.** The contact point wipes the port closed over
the port's own width rather than switching it shut at one angle. Treating
the port as a circular opening of radius $r_p=\tfrac12 R_c\gamma$ and
writing $u=R_c(\theta-\theta_{port})$ for the arc-length offset of the
contact from the port center,

$$
A_{eff}(\theta)=
\begin{cases}
\pi r_p^2 & u \le -r_p\\[2pt]
r_p^2\arccos(u/r_p)-u\sqrt{r_p^2-u^2} & |u| < r_p\\[2pt]
0 & u \ge r_p .
\end{cases}
$$

The compression chamber runs clockwise from the contact to the vane, so the
port belongs to it while $\theta$ is the smaller angle; $\theta$ and
$\theta_{port}$ are compared after normalizing to $[0,2\pi)$ without
wrapping their difference, which would destroy that ordering. No flow is
computed from $A_{eff}$; it is the geometric factor a later orifice or valve
model needs.

**True-geometry chamber volume.** The circular-rotor approximation of
Section 3.1 cannot support a recompression phase. At $2\pi-\delta$ it gives
a trapped volume of 0.00138 cm³ and drives it to zero by $2\pi$, so
$pV^n=\text{const}$ diverges; it also refuses to split the chambers at all
inside the seal-over window. The volume is therefore integrated on the real
boundary instead: the rotor contour of Section 3.3, both swing-bush pieces,
and the stepped vane.

Occupancy is evaluated on a square grid over the bore. The rotor material is
the outside-diameter disc minus the mouth cavity, where the cavity is the
polygon closed by the short outside-diameter arc across the mouth, the inlet
flat, and the mouth path through both lips and the groove. Each free cell
carries an axial height: the full cylinder height in general, zero under the
full-thickness vane segment, and the 16.2 mm open gap between the ledges
below it. The chamber split follows Section 3.1 — suction runs clockwise
from the vane to the contact, discharge continues back to the vane.

Assumptions:

- The swing bush holds a fixed attitude and only translates with the groove
  center. The rotor's $\pm 10.37^\circ$ relative tilt is not applied to it.
- The 0.010 mm bush films against the vane and the groove are treated as
  filled. The bush therefore seals the 4.0 MPa recess spaces away from the
  chamber gas, and those spaces carry no chamber volume.
- The mouth runs through the full axial height (Section 3.3).

Results at the supplied geometry:

Results at the supplied geometry, with the R2.1 vane-root blend of
Section 3.3 included:

| Quantity | Value |
|---|---:|
| mouth cavity free area at top dead center | 24.0 mm² |
| total gas volume at top dead center | 22.02 cm³ |
| same by the circular-rotor approximation | 21.52 cm³ (−2.3 %) |
| **clearance volume, trapped at $2\pi-\delta$** | **0.165 cm³** |
| same by the circular-rotor approximation | 0.00138 cm³ (120× smaller) |

The clearance volume is a derived quantity, not a supplied one, and it lands
inside the 0.1–0.5 cm³ band the pressure rule set asks for. **The dead
volume of this machine is the rotor mouth itself, not a separate discharge
recess:** at $2\pi-\delta$ the trapped discharge gas lies inside the mouth
cavity. Because the mouth never closes, no artificial volume floor is needed
to keep the recompression finite. Without the vane-root blend the clearance
would be 0.193 cm³; the blend removes about 0.028 cm³ of full-height gas
between the vane flanks and the bore (Section 3.3).

Grid convergence: the clearance volume moves from 0.16569 cm³ at a 0.050 mm
pitch to 0.16517 cm³ at 0.025 mm, a 0.3 % change, which sets the accuracy of
every volume quoted here.

**Port-timed pressure phases.** With a clearance volume available the
pressure rule gains the two phases Section 3.2 could not carry. On the
instantaneous view over $[0,2\pi)$:

| Range | Phase | Rule |
|---|---|---|
| $[0,\varphi)$ | re-expansion | $p=p_0\left(V(0^+)/V\right)^n$, clamped at $p_{suc}$ |
| $[\varphi,2\pi)$ | suction | $p=p_{suc}$ |
| $[\beta,\theta_{vo})$ | compression | $p=p_{suc}\left(V(\beta)/V\right)^n$ |
| $[\theta_{vo},2\pi-\delta)$ | delivery | $p=p_{dis}$ |
| $[2\pi-\delta,2\pi-\alpha)$ | recompression | $p=p_{dis}\left(V(2\pi-\delta)/V\right)^n$ |

$\theta_{vo}$ stays a pressure condition, found by bisection on the
polytropic pressure, so it moves when the operating condition moves; it is
never a constant. Recompression ends at the seal-over window entry, where
Section 3.1 already merges the chambers and no separate compression chamber
exists. Its end pressure is the residual $p_0$ that opens the next
re-expansion, so the cycle closes without hardcoding a residual: with the
supplied ports and the vane-root blend $p_0 \approx 10$ MPa. This value is
grid-sensitive, because it is set by the smallest volumes of the cycle
(R5.2 of the pressure rule set), so it is quoted only to an order.

That recompression peak is a **strict upper bound, not a prediction**:
leakage is neglected, and leakage is exactly what limits an over-compression
spike in a real machine. It is reported because the crank angle of the
over-compression peak is what a later bush-load model needs.

Because the mouth cavity adds volume that the outside-diameter disc
discards, the port-timed compression pressures run about 7 % below the
Section 3.2 rule at mid-stroke (1.73 MPa against 1.86 MPa at 180°). The two
rules are kept side by side: Section 3.2 remains the accepted constant
boundary condition and the regression baseline, and this section is
validated against it before it replaces it.

Position of $\beta$ barely matters. Moving the compression start from the
seal-over exit at 5.96° to $\beta=27.7^\circ$ shifts $\theta_{vo}$ from
220.568° to 220.921°, because the compression-chamber volume is flat near
top dead center. The unsettled reading of the supplied angles (Section 9)
therefore cannot change the compression stroke appreciably.

### 3.5 Indicated work and power

The indicated work of a chamber is the area of its closed pressure-volume loop
over one revolution, and the indicated power is that work times the shaft
frequency:

$$W = -\oint p\,\mathrm{d}V, \qquad P_\text{ind} = W\,f .$$

Sign convention: work done *on* the gas is positive, so the leading minus sign
makes a compressor return positive $W$ and $P_\text{ind}$. The pressure is the
port-timed rule of Section 3.4 (`port_timed_pressures`) and the volume is the
true-geometry trace (`CycleTrace`); the loop is the compression-chamber branch
(maximum to clearance volume) closed by the suction-chamber branch (clearance
to maximum volume).

The integral runs over $[0,\,2\pi-\alpha]$, stopping at the seal-over entry
rather than the full turn. Inside the seal-over window the two chambers merge
and the true-geometry split is degenerate, so the discharge volume balloons
there; excluding it removes that artefact, and the clearance-gas re-expansion
across the merge window carries negligible work.

At the supplied geometry ($0.82\to3.24$ MPa R410A, $n=1.07$, $22$ cm³ swept
volume, $30$ Hz): $W = 24.6$ J per revolution ($43.6$ J into the compression
chamber, $-19.0$ J returned by the filling suction chamber) and
$P_\text{ind} = 738$ W. Convergence with sample count is better than $0.1\%$.
This is a derived quantity on the prescribed pressure/volume boundary
conditions, not a transient thermodynamic solution: leakage, heat transfer, and
mechanical losses are excluded (the bush-film loss is estimated in Section 3.6).
Implemented in `mochi.indicated_work`.

**Validation — thermodynamic cross-check (no experimental P-V diagram).** No
directly comparable measured indicator diagram exists for this machine and point,
so the indicated work is validated by **triangulation** — two independent routes
to the same energy must agree:

* **Route A**, the P-V loop above, $W = 24.6$ J/rev, $P_\text{ind}=738$ W.
* **Route B**, thermodynamics: an ideal compressor's indicated work equals the
  delivered mass times the isentropic enthalpy rise, $W = m\,\Delta h_s$. CoolProp
  HEOS for R410A compressed isentropically from saturated vapour at 0.82 MPa
  ($0.85\,^\circ$C, $h=421.64$ kJ/kg, $\rho=31.4$ kg/m³) to 3.24 MPa
  ($74.1\,^\circ$C, $h=459.31$ kJ/kg) gives $\Delta h_s = 37.67$ kJ/kg (hard-coded
  like $n$, not a runtime dependency).

Because $n=1.07$ is nearly the isentropic pv-index (CoolProp: $1.064$ saturated to
$1.083$ at 20 K superheat), the specific indicated work $W/m$ must equal
$\Delta h_s$. With the delivered mass $m=0.650$ g/rev (Section 3.8),
$W/m = 37.8$ kJ/kg agrees with $\Delta h_s = 37.67$ kJ/kg to **0.4 %**, and the
power routes agree at $738$ vs $735$ W — an independent confirmation of the P-V
number with no measured diagram. Both are the reversible **isentropic-ideal**
lower bound (the full-displacement, $\eta_v=1$ bound is $781$ W); a real indicator
diagram is +15–40 % higher through valve (Section 3.8), throttling, and
heat-transfer losses. Implemented in `mochi.thermo_check`
(`isentropic_cross_check`).

### 3.6 Swing-bush lubrication film

Each swing-bush piece runs on two oil films of uniform gap $h = 10\,\mu$m: a
**curved** film between the bush outside diameter and the rotor groove
(length = the $14.76$ mm cylindrical arc) and a **flat** film between the bush
flat and the fixed vane (length = the $11.94$ mm contact, shortening to
$11.47$ mm near BDC where the bush overhangs the R1.5 vane-tip round, Section
3.3), both spanning the $20.983$ mm bush height $H$. The bush holds a fixed attitude and only
translates with the groove centre (Section 3.4 assumption), so the two films
have different no-slip drivers:

$$U_\text{flat} = \omega\,\frac{\mathrm{d}\,y_\text{groove}}{\mathrm{d}\theta},
\qquad
U_\text{curved} = \omega\,\frac{\mathrm{d}\,\psi}{\mathrm{d}\theta}\,r_\text{bush},$$

where $y_\text{groove}$ is the groove-centre position along the vane and $\psi$
the rotor orientation. The flat film is driven by the bush translation (up to
$0.86$ m/s near $\theta=90^\circ$) and the curved film by the rotor's
$\pm10.37^\circ$ orientation swing about the fixed-attitude bush (up to
$0.27$ m/s) — **not** the shaft speed, because the bush moves with the groove.

Each film is a one-dimensional incompressible thin film. With a uniform gap,
mass conservation makes the flux constant and the pressure linear between the
two end pressures $p_a$ (chamber end) and $p_b$ (recess end):

$$p(x) = p_a + (p_b-p_a)\frac{x}{L}, \qquad
\tau = \frac{\mu U}{h}, \qquad
q = \frac{U h}{2} - \frac{h^3}{12\mu}\frac{\mathrm{d}p}{\mathrm{d}x} .$$

The Couette wall shear $\tau = \mu U/h$ is the traction on the bush; the film
carries the gas load as the linear (Poiseuille) pressure set by the boundary
pressures, and $q$ is the leakage flux per unit width (Couette drag minus
Poiseuille back-flow). **A pure Couette film generates no load-bearing pressure
of its own** — the pressure field comes entirely from the gas boundary
conditions. The boundary conditions are the chamber pressure at the mouth/OD
end (the suction pressure on the IN piece, the compression pressure on the OUT
piece, Section 3.4) and the $4.0$ MPa recess pressure `RECESS_PRESSURE_PA` at
the inner end; this is the first use of that constant. Their exact film-edge
assignment is a first-cut modelling choice (Section 9 open item), so the end
pressures are exposed as parameters.

Lubricant: POE ISO VG68 at $\sim80\,^\circ$C, dynamic viscosity
$\mu = 0.010$ Pa·s (neat oil; refrigerant dilution neglected). At the supplied
geometry the wall shear reaches $\sim860$ Pa (flat) and $\sim270$ Pa (curved),
the cycle-mean film friction loss is about $0.2$ W ($\sim0.03\%$ of the
indicated power), and the leakage is of order $0.1$ mL/s. **Validity:** the
uniform-gap assumption means no hydrodynamic (wedge) load — a Reynolds/
eccentricity model would be needed for that; the fixed-attitude and 10 µm
uniform-film assumptions are the Section 9 open items. Implemented in
`mochi.bush_film`.

### 3.7 Gap leakage and volumetric efficiency

The port-timed recompression branch (§3.4) over-compresses the trapped clearance
gas to a ~10 MPa spike, reported there as a strict upper bound because nothing
bleeds the over-pressurised gas out. This section adds a **single equivalent
leakage orifice** from the compression/recompression chamber to the suction side
and integrates the chamber **mass balance** over one revolution — capping the
spike and yielding the **volumetric efficiency**. It is a *parallel* model; the
accepted `port_timed_pressures` baseline is untouched (as §3.2 and §3.4 sit side
by side). Implemented in `mochi.leakage`.

**Mass-scaled polytropic pressure.** The model has no absolute temperature, so
the mass balance keeps the polytropic form. A chamber sealed at $(p_{ref},
V_{ref}, m_{ref})$ has, with a changed mass,

$$p = p_{ref}\,\frac{m}{m_{ref}}\left(\frac{V_{ref}}{V}\right)^{n}, \qquad
\rho = \rho_{suc}\left(\frac{p}{p_{suc}}\right)^{1/n},$$

the ideal-gas identity $p=mRT/V$ with $T$ on the fixed-mass polytropic path
($n=1.07$); the $m=m_{ref}$ case is the §3.4 rule. Two R410A saturated-vapour
properties are introduced, hard-coded from CoolProp exactly as $n$ was (not a
runtime dependency): suction density $\rho_{suc}=31.4$ kg/m³ and isentropic
(nozzle) exponent $\gamma = a^2\rho/p = 1.10$ at 0.82 MPa.

**Leakage orifice.** The rotor–cylinder tangency is an ideal zero-gap line seal
in the model and no end-face path exists, so an **effective** leakage gap is
assumed (editable; all paths lumped into one equivalent orifice): gap
$\delta = 5\,\mu$m over the cylinder height, $C_d = 0.6$. The mass flow is the
compressible isentropic nozzle, choked-capable,

$$\dot m = C_d A \sqrt{\tfrac{2\gamma}{\gamma-1}\,p_{up}\rho_{up}\,
\bigl(r^{2/\gamma}-r^{(\gamma+1)/\gamma}\bigr)},\quad
r=\tfrac{p_{down}}{p_{up}},$$

with the choked limit $\dot m = C_d A\sqrt{\gamma p_{up}\rho_{up}}\,
(2/(\gamma+1))^{(\gamma+1)/(2(\gamma-1))}$ for $r\le(2/(\gamma+1))^{\gamma/
(\gamma-1)}$.

**Cycle mass balance.** From the seal angle $\beta$ the compression chamber is
marched in crank angle: closed compression leaking to suction until the
valve-opening pressure, delivery at $p_{dis}$, then port-shut recompression of
the trapped clearance gas, $\mathrm{d}m/\mathrm{d}\theta = -\dot m_{leak}/\omega$.

**Results and honest scope.** At the supplied geometry with $\delta=5\,\mu$m the
recompression peak caps only modestly, from $10.0$ to $\approx 9.7$ MPa: the
recompression window is short ($\sim 7^\circ \approx 0.65$ ms at 1800 rpm), so
even choked flow bleeds just a few percent of the trapped mass. The spike is
therefore largely inherent to the **geometric** port-closing model (gas is
trapped and compressed); what fully limits it in a real machine is the discharge
**reed valve reopening** (backflow), a separate valve model. The capping grows
with the assumed gap ($8.8$ MPa at $20\,\mu$m). The solid deliverable is the
**volumetric efficiency** $\eta_v = m_{deliv}/(\rho_{suc}V_{disp})$: clearance
re-expansion alone gives $\eta_v = 0.94$, and the $5\,\mu$m leak lowers it to
$0.92$ (both within the typical R410A rotary $0.85$–$0.95$). **Validity:**
single lumped orifice with an assumed gap; mass-scaled polytropic (no energy
equation); leak to suction only (discharge-side backflow and the reed valve not
modelled). Implemented in `mochi.leakage` (`orifice_mass_flow`, `leaky_cycle`).

### 3.8 Discharge reed valve and overpressure

Sections 3.4 and 3.7 hold the compression chamber at the discharge-port pressure
$p_{dis}$ through delivery. A real discharge is gated by a **reed (check)
valve**: it lifts only when the chamber pressure exceeds the line, and it flows
through the geometric port opening the rotor-cylinder contact closes over the
last few degrees. Adding the valve to the leakage mass balance (Section 3.7) lets
the pressure **float** during delivery, so the finite, shrinking port area
produces a physical **discharge overpressure**. This is a *parallel* model on top
of Section 3.7; the baselines are untouched. Implemented in `mochi.reed_valve`.

**Quasi-static check valve.** The reed opens at $p_{dis}(1+f_v)$ (the 5% valve
rise of Section 3.4 as the spring preload) and vents to the line through the
geometric port area $A(\theta)=$ `ports.port_open_area_m2` — full $\approx
18.4\,$mm² through most of delivery, tapering to zero over the discharge window
$[339.6^\circ, 346.8^\circ]$ as the contact covers the port. The mass flow is the
compressible orifice of Section 3.7, $\dot m_{valve}=\dot m\big(p, p_{dis},
A(\theta)\big)$, with the pressure from the same mass-scaled polytropic; the
tangency-gap leak to suction runs alongside.

**Emergent overpressure.** When the valve cannot pass the swept mass rate through
$A(\theta)$ at a small pressure difference, the chamber pressure rises above
$p_{dis}$. At the supplied geometry the delivery peak is $\approx 3.75$ MPa
($+0.5$ MPa over the $3.24$ MPa line) and the overpressure loss — the extra
indicator-diagram work $\oint (p-p_{dis})\,(-\mathrm{d}V)$ over delivery — is
$1.7$ J/rev, $\approx 51$ W ($\sim 7\%$ of the indicated power). It scales with
the assumed port size and falls toward the reed-opening floor as the area grows
(178 W at half area, 7 W at triple area).

**Valve-aware indicated power.** Adding the overpressure area to the §3.5 baseline
gives the valve-reflected indicated work, $W_{valve} = W_\text{ind} + \oint
(p-p_{dis})(-\mathrm{d}V) = 24.6 + 1.7 = 26.3$ J/rev, $P = 789$ W — the first step
of the loss ladder above the isentropic-ideal $738$ W (§3.5). This is a **separate
performance term, not propagated into the mechanical loads**: the force, torque,
bearing reaction, and shaft power (§4.5–4.7) are all unified on the baseline
(reed-valve-free) indicated $738$ W. It is exposed on `ValvedCycle`
(`valve_indicated_power_w`) and drawn only on the §3.5 P-V diagram, as a
clearly-labelled dashed/shaded overpressure branch. Real machines add suction
throttling and heat transfer on top for the full +15–40 %.

**No back-flow; spike not capped.** Because the chamber pressure stays at or above
$p_{dis}$ throughout the open window, there is **no back-flow** in this
quasi-static limit — back-flow requires the reed to stay open past pressure
reversal, a valve-**dynamics** effect (reed mass/stiffness, Section 9). And after
the geometric port closes ($\theta>346.8^\circ$) the area is zero, so the valve
cannot vent the recompression pocket: that $\approx 9.6$ MPa spike is still
limited only by leakage (Section 3.7), confirming the Section 3.7 finding. The
volumetric efficiency is $\eta_v\approx 0.96$ (clearance-limited; the overpressure
is a **work** loss, not a mass loss). **Validity:** no reed inertia/flutter/back-
flow; geometric port area with an assumed discharge coefficient; mass-scaled
polytropic. Implemented in `mochi.reed_valve` (`valved_cycle`).

## 4. Bodies, frames, and degrees of freedom

This section records the accepted planar rigid-body definition of the
mechanism. No equation of motion is evaluated yet; Section 4.4 keeps the
dynamics interface for later freed degrees of freedom.

### 4.1 Bodies

| Body | Role | Planar DOF before constraints |
|---|---|---:|
| Cylinder assembly (bore, stepped vane, ports, end plates) | ground; carries the global frame | 0 |
| Crankshaft (eccentric journal, axis through $O$) | drive | 1 ($\theta$) |
| Rotor (OD 68 mm, circular groove, recess channel) | compression | 3 ($x_r$, $y_r$, $\phi$) |
| Swing bush (in the rotor groove, wraps the vane) | rotor-vane interface | 3 |
| Gas | massless pressure boundary condition (Section 3.2) | — |

Assumptions: rigid bodies; planar model in the cylinder cross-section
plane (axial motion, tilt, and end-plate effects neglected); gravity
excluded because the shaft is vertical, so gravity acts normal to the
plane; ideal tangency at the rotor-cylinder contact; ideal radial bearing
holding the crank axis at $O$.

### 4.2 Frames and angle conventions

- **Global inertial frame $G$**, fixed to the cylinder: origin $O$ at the
  cylinder center, $\hat{x}$ toward the inlet side, $\hat{y}$ toward the
  vane, $\hat{z}=\hat{x}\times\hat{y}$ out of the drawing (right-handed).
- **Sign discipline:** the crank angle $\theta$ is measured from $+\hat{y}$,
  positive clockwise — a rotation about $-\hat{z}$. Any vector moment or
  angular-velocity balance must be written about $+\hat{z}$
  (counter-clockwise positive) with $\omega_z=-\dot{\theta}$.
- **Crank direction** $\hat{e}_c=(\sin\theta,\cos\theta)$; the rotor center
  is $\mathbf{r}=e\,\hat{e}_c$.
- **Rotor frame:** origin at the rotor center; $\mathbf{u}=(\mathbf{h}-
  \mathbf{r})/L$ toward the groove center, $\mathbf{n}=\hat{z}\times
  \mathbf{u}$; orientation $\phi=\operatorname{atan2}(u_y,u_x)$,
  counter-clockwise from $+\hat{x}$ ($\phi=90^\circ$ at TDC), matching the
  implementation and GUI display.
- **Vane:** part of the ground body; centerline $x=0$, fixed tip at
  $(0,\,e+a_v)$.

### 4.3 Degrees of freedom and constraints

Unconstrained planar count: crank $1$ + rotor $3$ + bush $3$ = $7$.
Constraints:

1. Crank journal bearing: $\mathbf{r}=e(\sin\theta,\cos\theta)$ — 2.
2. Bush-groove pin joint: the bush center coincides with the groove center
   $\mathbf{h}$ — 2.
3. Bush-vane prismatic joint: the bush is centered on the vane axis and
   cannot rotate ($x_b=0$, $\phi_b=0$) — 2.

Net degrees of freedom: $7-6=1$, the crank angle $\theta$. Combining
constraints 2 and 3 yields the groove-center condition $x_h=0$ implemented
in Section 3; rotor position and orientation and the bush position are
functions of $\theta$ alone. The crank angle remains **prescribed**,
$\theta(t)=2\pi f t$, so force analysis on this mechanism is quasi-static:
gas and reaction loads are evaluated along the prescribed motion. Freeing
$\theta$ later requires the drive-torque model and inertia data of
Section 4.4.

### 4.4 Dynamics interface for freed degrees of freedom (not implemented)

The visualization evaluates no equation of motion. A future dynamics module
may introduce generalized position $\mathbf{q}$, velocity
$\mathbf{v}=\dot{\mathbf{q}}$, and mass matrix $\mathbf{M}$:

$$
\mathbf{M}\ddot{\mathbf{q}} = \mathbf{F}_{total}.
$$

The total is assembled from named contributions:

$$
\mathbf{F}_{total} =
\mathbf{F}_{drive} +
\mathbf{F}_{gas} +
\mathbf{F}_{contact} +
\mathbf{F}_{fluid} +
\mathbf{F}_{support} +
\mathbf{F}_{other}.
$$

Only accepted, configured terms are evaluated. An absent model contributes
zero and must be reported in case metadata; it must not be approximated
silently.

If a linear support is later justified, a possible constitutive term is

$$
\mathbf{F}_{support} =
-\mathbf{C}\mathbf{v}
-\mathbf{K}(\mathbf{q}-\mathbf{q}_{ref}),
$$

where $\mathbf{C}$ has units N s/m and $\mathbf{K}$ has units N/m. This is an
interface example, not yet an accepted model of the compressor.

If dynamic rotation is added, define an orientation $\phi$, moment of inertia
$I$, and a moment balance $I\ddot{\phi} = M_{total}$ with the same explicit
sign discipline. Do not infer rotation from planar force data alone.

### 4.5 Gas pressure force and torque

This is the first evaluated term of the $\mathbf{F}_{total}$ assembly: the
$\mathbf{F}_{gas}$ contribution of the chamber gas on the rotor, and the net gas
torque about the crank axis $O$ -- the indicated torque. It is quasi-static
(net DOF $=1$, $\theta$ prescribed): the loads are read off the accepted chamber
pressures along the prescribed motion, with no equation of motion.

**Exerting/receiving bodies.** The chamber gas presses on the rotor OD and on
the two vane flanks. The vane is part of the ground body (Section 4.2), so its
gas load reacts to the cylinder and carries no shaft torque; only the rotor load
drives the crank. The vane force is reported to close the balance.

**Rotor force (closed form).** In the circular-rotor approximation the rotor OD
is a circle of radius $R_r$ about the rotor centre
$\mathbf{r}=\mathbf{O}_r=e(\sin\theta,\cos\theta)$, height $H$. A chamber at
uniform pressure $p$ acts along the outward normal
$\hat{n}(\varphi)=(\sin\varphi,\cos\varphi)$ (angle $\varphi$ from $+\hat{y}$,
clockwise), so the force on an arc is fixed by its endpoints,

$$
\mathbf{F}_{arc} = -\,p\,H\!\int \hat{n}\,\mathrm{d}\ell
= -\,p\,R_r H\,\bigl(\cos\varphi_s-\cos\varphi_e,\;
\sin\varphi_e-\sin\varphi_s\bigr).
$$

The suction chamber wets the OD from the inlet vane foot $\varphi_{in}=
\arcsin\!\big((w/2-e\sin\theta)/R_r\big)$ clockwise to the rotor-cylinder contact
at $\varphi=\theta$; the compression chamber wets it from the contact clockwise
(the long way, through the bottom) to the outlet foot $\varphi_{out}=
\arcsin\!\big((-w/2-e\sin\theta)/R_r\big)$. The vane footprint over the top is
excluded. The rotor force is the sum of the two arc terms with $p_{suc}$ and
$p_{comp}$ (Section 3.4).

**Torque about $O$.** Because every pressure element acts radially through
$\mathbf{O}_r$, the resultant $\mathbf{F}_{rotor}$ passes through $\mathbf{O}_r$
and exerts no couple about it, so the moment about $+\hat{z}$ (Section 4.2) is

$$
T_{gas}(\theta) = \mathbf{O}_r\times\mathbf{F}_{rotor}
= e\sin\theta\,F_y - e\cos\theta\,F_x .
$$

**Independent test (cross-check with Section 3.5).** The gas does work on the
rotor at the rate $\mathbf{F}_{rotor}\cdot\dot{\mathbf{O}}_r=-T_{gas}\dot\theta$,
and that work equals $\sum_c p_c\,\dot V_c$ (only the rotor boundary moves).
Integrating over one revolution gives the identity

$$
\oint T_{gas}\,\mathrm{d}\theta = -\oint p\,\mathrm{d}V = W,
$$

the indicated work of Section 3.5. This holds on the *same* volume basis as the
force: the circular-rotor crescent volumes $V=A\,H$ from `chamber_areas`. At the
supplied geometry $\oint T_{gas}\,\mathrm{d}\theta = 23.8$ J matches the
circular-rotor $-\oint p\,\mathrm{d}V$ to one part in $10^6$, and the implied
power $715$ W is within $3\%$ of the $738$ W Section 3.5 headline; the residual
is the circular-rotor vs true-geometry ($22$ cm³) volume gap.

**True-geometry refinement (mouth-aware).** The $3\%$ gap is closed by integrating
the force and torque over the *same* shape the pressure and volume already use --
the real rotor contour of Section 3.3 (`rotor_profile`), including the **mouth**
cavity the circular disc discards. Two facts change:

* *Force.* The pressure is integrated over the true material boundary (OD arc plus
  the gas-wetted mouth walls, skipping spans a bush piece seals or the vane covers,
  by the same predicates as the Section 3.4 volume). The mouth cavity gas presses
  the lips outward and **lowers the net rotor load $\sim20\%$** below the OD-disc
  closed form (at $\theta=235^\circ$, $2.28$ kN vs $3.01$ kN); a grid-converged
  raster integral of the identical cell classification confirms it to $<1\%$.

* *Torque.* On the true contour the resultant no longer passes through $O_r$, the
  rotor both orbits and swings, and the true indicated work includes work on the
  moving bush -- so $T_{gas}=O_r\times F$ no longer holds. The work-conjugate
  **shaft** gas torque is taken from virtual work,
  $T_{gas}(\theta)=-(p_c\,\mathrm{d}V_c/\mathrm{d}\theta+p_s\,\mathrm{d}V_s/
  \mathrm{d}\theta)$ on the true-geometry trace, whose crank-angle integral equals
  the true indicated work $\oint T_{gas}\,\mathrm{d}\theta=24.6$ J $\to 738$ W
  **exactly** (vs the circular $715$ W). This is the headline for the gas force,
  the crank-pin reaction $\mathbf{R}_j=-\mathbf{F}_{gas}$ (Section 4.6), and the
  indicated gas torque; the circular closed form above is retained as the fast
  analytic check. Implemented in `mochi.true_gas_force` (`true_gas_load`,
  `true_gas_torque_work_j`, `peak_rotor_force_n`).

**Validity.** The circular closed form ignores the rotor slot, cutout, and mouth
(Section 3.1); the true-geometry refinement restores the mouth but still assumes
full cylinder height for the force (the vane's axial ledges, a second-order axial
correction, are not applied). Both ground the vane and exclude the seal-over
window $\lvert\theta\rvert<\alpha$, where the chambers merge and no separated arcs
exist. The vane force is a first cut: the two vertical flanks from the rotor foot
to the bore, with the tip and root fillets neglected. Implemented in
`mochi.gas_force` (`gas_load`, `gas_torque_work_j`) and `mochi.true_gas_force`.

### 4.6 Crank-pin bearing reaction and drive torque

This closes the rotor's quasi-static force balance and gives the two mechanism
loads a designer needs: the crank eccentric (journal) bearing reaction, and the
drive torque the shaft must supply. It resolves the §9 open items "identify all
applied, gas, contact, friction, and fluid forces" and "select the
contact/friction treatment."

**Statical indeterminacy and its resolution.** The crank journal fixes the rotor
centre on the line $O$–$O_r$ (§4.3 constraint 1), and the rotor is ideally
tangent to the cylinder at $C$ with $O,O_r,C$ colinear (§3.1). The
rotor–cylinder contact normal force $N_c$ (radial, through $O_r$) and the radial
component of the journal reaction $\mathbf{R}_j$ (also through $O_r$) lie on the
same line, so translational balance cannot split them. A moment balance about the
rotor centre $O_r$ gives no help either: $N_c$, $\mathbf{R}_j$, and the gas
resultant $\mathbf{F}_{gas}$ (§4.5) all pass through $O_r$ with zero moment arm,
leaving only friction torques. Under the rigid, frictionless, inertialess model
(no mass/inertia data exists, §4.3) $N_c$ is therefore **not a determinate
rigid-body quantity** — a genuine result, not a gap. The accepted treatment
(Yanagisawa & Shimizu; Aw & Ooi review) lumps the radial load onto the journal
bearing and treats the contact as friction-only; a determinate nonzero $N_c$
would need an oil-film/Hertzian **compliance** model (a later rung). Hence

$$\mathbf{R}_j(\theta) = -\mathbf{F}_{gas}(\theta),$$

the journal carries the whole gas force. On the mouth-aware true force (Section
4.5) its peak magnitude is $\approx 2.5$ kN (at $\theta\approx227^\circ$, near the
peak-pressure crank angle) -- the rotor mouth lowers it $\sim20\%$ below the
$\approx 3.3$ kN circular closed form; this refined value is the bearing design
load, and the locus of $\mathbf{R}_j(\theta)$ over a revolution is the polar
bearing-load diagram (circular locus overlaid for reference).

**Drive torque.** In quasi-static motion the shaft supplies the indicated gas
torque plus the mechanical friction:

$$T_{drive}(\theta) = T_{gas}(\theta) + T_{fric}(\theta), \qquad
T_{fric}(\theta) = \frac{P_{bush}(\theta) + P_{journal}(\theta)}{\omega},$$

with $T_{gas}$ from §4.5, $P_{bush}$ the swing-bush film dissipation of §3.6, and
$P_{journal}$ the crank-pin journal bearing friction of §4.7, both referred to
the shaft ($\omega=2\pi f$). Dry rotor-cylinder contact friction is still omitted
(indeterminate contact force, §4.6; no coefficient), so this remains a lower
bound on the true drive torque.

**Shaft power.** The shaft supplies the indicated gas work plus the friction, so
it **exceeds** the indicated power by the friction alone. It is built on the
isentropic-ideal **baseline** indicated power (§3.5, $738$ W):

$$P_{shaft} = P_{ind} + P_{bush} + P_{journal}
 = 738 + 0.2 + 8.9 = 747\ \text{W}.$$

The discharge reed-valve overpressure (§3.8, $\approx 51$ W, giving a valve-aware
indicated $789$ W) is a **separate performance term, not propagated into the
loads** — it is *not* added to the shaft power here; all of the mechanical loads
(force, torque, bearing reaction, shaft power) are unified on the baseline
(reed-valve-free) indicated work. Neither is the shaft power the drive-torque
work: the drive-torque integral
$\oint T_{drive}\,\mathrm{d}\theta = \oint T_{gas}\,\mathrm{d}\theta + \oint
T_{fric}\,\mathrm{d}\theta$ integrates the gas torque on the **circular-rotor**
basis (§4.5), whose power $714.5$ W is ~3% below the true indicated $738$ W — so
$\oint T_{drive}\,\mathrm{d}\theta \cdot f = 724$ W is **not** the shaft power;
using it would let the ~3% basis gap ($23$ W) swamp the smaller friction ($9$ W)
and make the shaft power spuriously fall below the indicated power. The exact
torque-integral split ($\oint T_{drive} = \oint T_{gas} + \oint T_{fric}$) remains
the internal cross-check binding this rung to §3.5, §3.6, §4.5, and §4.7.
**Validity:** journal-lumped radial
load (no independent contact force), bush + journal friction (dry contact
omitted), inertia neglected. Implemented in `mochi.bearing_load`
(`mechanism_load`, `shaft_work_j`).

### 4.7 Journal bearing friction

The crank-pin (eccentric journal) bearing carries the reaction $\mathbf{R}_j$ of
§4.6 (peak $\approx 2.5$ kN, mouth-aware) and is the **dominant mechanical loss** — an order of
magnitude above the swing-bush film. It is modelled with **Petroff's law** for a
concentric, full oil film, using the same POE VG68 lubricant as §3.6
($\mu = 0.010$ Pa·s):

$$T_j(\theta) = \frac{2\pi\,\mu\,\lvert\omega_{rel}\rvert\,r_j^{3}\,L_j}{c_j},
\qquad P_j(\theta) = T_j\,\lvert\omega_{rel}\rvert,$$

with journal radius $r_j$, length $L_j$, radial clearance $c_j$. The oil shears
at the **relative** speed between the crank pin (turning at the shaft speed) and
the rotor bore (turning at the rotor spin $\omega\,\mathrm{d}\phi/\mathrm{d}
\theta$):

$$\omega_{rel}(\theta) = \omega\left(1 - \frac{\mathrm{d}\phi}{\mathrm{d}\theta}
\right).$$

Because the rotor only swings $\pm10.37^\circ$ (§3.3), $\omega_{rel}$ stays near
the shaft speed, $154$–$222$ rad/s. Referred to the shaft the friction adds
$P_j/\omega$ to the drive torque (§4.6).

**Geometry.** The crank-pin journal radius and length are taken from the **CAD
model**: $r_j = 14.2$ mm (bounded by the swing-bush groove it must clear,
$r_j < L - r_{cut} = 25 - 8 = 17$ mm) and $L_j = 21$ mm (the rotor height). The
radial clearance is **not drawn** — the CAD shows a nominal full-contact fit — so
it stays an **assumed** oil-film value $c_j = 15\,\mu$m ($r_j/c_j \approx 947$),
exposed as an editable constant. These give $T_j \approx 0.048$ N·m and a
cycle-mean $P_j \approx 8.9$ W ($\approx 1.2\%$ of the indicated power) — still
dwarfing the $0.2$ W bush film.

**Validity.** Petroff assumes a **concentric** journal, so $T_j$ is
load-independent; the Sommerfeld number $S = (\mu N / P)(r_j/c_j)^2$ with
$N = \lvert\omega_{rel}\rvert/2\pi$ and projected pressure $P = \lvert
\mathbf{R}_j\rvert/(2 r_j L_j)$ is reported as a regime indicator. At the peak
load $S \approx 0.07$ — a heavily loaded bearing running well off centre, so
Petroff **underestimates** the true friction; the short-bearing/Sommerfeld
(eccentricity-coupled) model that uses $\mathbf{R}_j(\theta)$ is **§4.9** (running
$\varepsilon \approx 0.71$, journal loss $\approx 10.8$ W vs this $8.9$ W).
Boundary/dry friction and thermal effects are excluded. Implemented in
`mochi.journal_bearing` (`petroff_friction`, `journal_relative_speed_rad_s`).

### 4.8 Body mass properties (dynamics inertia data)

The multi-body-dynamics rung (D2, prescribed-speed rotordynamics) needs the rotor
and swing-bush **inertia**, which the quasi-static rungs never used. The material is
**not confirmed**, so the density is an **editable assumption** and the mass and
inertia rescale linearly with it:

$$
m = \rho\,H\!\int_{A}\!\mathrm{d}A,\qquad
I_z = \rho\,H\!\int_{A}\!r^2\,\mathrm{d}A,
$$

with the solid cross-section $A$ taken from the *same* regions the pressure/volume
model classifies (Section 3.4): the **rotor** OD disc minus the mouth cavity minus
the central crank-pin bore (radius $r_j = 14.2$ mm), at full height $H = 21$ mm;
the **swing bush** as its two pieces (`SwingBush.occupies`), at height $H -
2\times0.0085$ mm. The integrals run on the same raster as `true_chamber_volumes`
and $I_z$ is the polar second moment about the rotor centre $O_r$ (rotor) / the
groove centre (bush) — the axes the crank-pin journal and the D2 tilt/whirl use.

At the assumed densities (grey cast iron $\rho_r = 7200$ and steel bush $\rho_b =
7850$ kg/m³) this gives **rotor $m_r \approx 0.42$ kg, $I_{r} \approx 2.8\times
10^{-4}$ kg·m²** (within $0.4\%$ of the disc-minus-bore annulus $\tfrac12 m_r(R_r^2
+ r_j^2)$, the mouth removing a little off-axis material) and **bush $m_b \approx
13$ g, $I_b \approx 5.5\times10^{-7}$ kg·m²**. Editing `ROTOR_DENSITY_KG_M3` /
`BUSH_DENSITY_KG_M3` (or passing `density_kg_m3`) rescales these immediately.
**Validity:** the rotor is modelled as a full-height annulus minus the mouth (the
axial recess-channel step and any web/lightening are not resolved); the crank pin
itself belongs to the (prescribed) crankshaft, not the rotor. **This inertia is
dormant data:** it is not yet consumed by any solver -- the quasi-static rungs use
no inertia, and no equation of motion exists -- so it takes effect only when the
future D2 dynamic rung is built. Implemented in `mochi.mass_properties`
(`rotor_mass_properties`, `bush_mass_properties` → `MassProperties`).

### 4.9 Short-bearing (Ocvirk) journal force -- eccentricity-coupled film

Section 4.7 friction is Petroff's **concentric** ($\varepsilon = 0$) law, which is
load-independent and, at the peak reaction ($S \approx 0.07$), *underestimates* the
friction. This rung lifts the concentric assumption to the **short-bearing
(Ocvirk)** hydrodynamic film: a closed-form crank-pin force that depends on the
journal **eccentricity ratio** $\varepsilon = e/c_j$ (0 = centred, 1 = metal
contact), its rate $\dot\varepsilon$, and the film **entrainment speed** -- the
first eccentricity-dependent bearing law and the $\mathbf{F}(\varepsilon,
\dot\varepsilon, \psi)$ element a future D2 rotordynamics EOM would integrate. It is
still closed form (**no lubrication PDE is solved** -- that is the later L2 Reynolds
rung).

**This rung is quasi-static, not dynamic.** It supplies the eccentricity-dependent
*force law*, but here $\varepsilon(\theta)$ is obtained from a **steady
force balance at each crank angle** ($\dot\varepsilon = 0$, no squeeze/whirl
transient) -- **no equation of motion is written and nothing is integrated in
time.** The rotor mass/inertia (§4.8) and the squeeze term $\dot\varepsilon$ are
present in the code but stay dormant; they activate only when the (not-yet-built) D2
rung frees the rotor's motion and marches it in time. The result below is therefore
the quasi-static equilibrium the future dynamic film must reduce to.

**Model.** With $\beta$ measured from the point of maximum film thickness, the film
is $h = c_j(1 + \varepsilon\cos\beta)$. The short-bearing approximation keeps only
the **axial** pressure flow of the Reynolds equation (formally $L/D \lesssim 0.5$;
here $L_j/2r_j = 0.74$, an order-of-magnitude application), so the pressure and its
integral are closed form. Under the **$\pi$-film (Gümbel/half-Sommerfeld)**
cavitation boundary the force resolves into a component along the line of centres
($F_e$, toward maximum film -- the load-supporting direction) and perpendicular
($F_t$):

$$F_e = \frac{\mu r_j L_j^{3}}{c_j^{2}}\left[\dot\varepsilon\,P - \varepsilon\,
\Omega\,Q\right],\qquad
F_t = \frac{\mu r_j L_j^{3}}{c_j^{2}}\left[\dot\varepsilon\,Q - \varepsilon\,
\Omega\,S\right],$$

$$P = \frac{\pi(1+2\varepsilon^2)}{2(1-\varepsilon^2)^{5/2}},\quad
Q = \frac{-2\varepsilon}{(1-\varepsilon^2)^{2}},\quad
S = \frac{\pi}{2(1-\varepsilon^2)^{3/2}}.$$

For **pure rotation** ($\dot\varepsilon = 0$) this is the textbook Ocvirk load
capacity $\lvert W\rvert = \dfrac{\mu\,\Omega\,r_j L_j^{3}}{2c_j^{2}}
\dfrac{\varepsilon\sqrt{16\varepsilon^2 + \pi^2(1-\varepsilon^2)}}{(1-\varepsilon^2
)^{2}}$ and the attitude angle $\tan\psi = \dfrac{\pi\sqrt{1-\varepsilon^2}}
{4\varepsilon}$ (both reproduced to machine precision in the tests). $\lvert W\rvert$
is monotonic in $\varepsilon$, so it **inverts** for the running eccentricity of a
given reaction.

**Entrainment vs shear speed.** The wedge is driven by the **mean** surface speed
$\Omega = \tfrac12(\omega_j + \omega_b) - \dot\psi = \omega - \tfrac12\lvert
\omega_{rel}\rvert$ (crank pin at $\omega$, rotor bore at $\omega\,\mathrm{d}\phi/
\mathrm{d}\theta$; $\Omega = \omega/2$ when the rotor does not spin -- the classic
half-speed), whereas the **shear/friction** uses their *difference* $\omega_{rel}$
of §4.7. Over the cycle $\Omega \approx 77$–$111$ rad/s.

**Running eccentricity and film (headline).** Balancing the film against the true
gas reaction $\lvert\mathbf{R}_j(\theta)\rvert$ (§4.6) at $\Omega(\theta)$ gives the
eccentricity trace $\varepsilon(\theta)$: it ranges $0.13$–$\mathbf{0.71}$ and peaks
**$\varepsilon \approx 0.71$ at $\theta \approx 226^\circ$** -- the peak reaction
($\lvert\mathbf{R}_j\rvert \approx 2.49$ kN), where the minimum film thins to
$h_{min} = c_j(1-\varepsilon) \approx \mathbf{4.4\ \mu m}$ (from the assumed
$c_j = 15\ \mu$m). The bearing runs safely off-centre, not concentric.

**Eccentric friction (reduction check).** The Couette shear over the eccentric film
integrates to

$$T_f = \frac{2\pi\mu\lvert\omega_{rel}\rvert r_j^{3} L_j}{c_j\sqrt{1-\varepsilon^2}}
= \frac{T_{j,\text{Petroff}}}{\sqrt{1-\varepsilon^2}},$$

so as $\varepsilon \to 0$ it **recovers Petroff exactly** (the validation-ladder
reduction), and off-centre it is larger. The cycle-mean journal loss rises from the
Petroff $\approx 9.2$ W to $\approx \mathbf{10.8\ W}$ ($\sim 18\%$), confirming the
§4.7 statement that the concentric model underestimates the friction. (The
pressure-flow shear term is second order and omitted at this L1 fidelity.)

**Validity.** $\pi$-film (superposed rotation + squeeze) cavitation; short-bearing
(axial-flow-only) at a moderate $L/D$; isothermal, iso-viscous oil; the crank-pin
**contact** with the cylinder is a *separate* Hertz/EHL element (not this
hydrodynamic journal). Rotor mass/inertia (§4.8) and the squeeze term
$\dot\varepsilon$ enter the force but are exercised only when D2 frees the rotor
motion; here $\varepsilon(\theta)$ is the quasi-static (steady-eccentricity)
equilibrium. Implemented in `mochi.ocvirk_bearing` (`short_bearing_force` →
`OcvirkForce`, `static_load_capacity_n`, `equilibrium_eccentricity_ratio`,
`eccentric_friction_torque_nm`, `crank_pin_entrainment_speed_rad_s`,
`eccentricity_cycle` → `EccentricityCycle`); figures
`bearing_load/{journal_eccentricity, eccentric_friction_power}.png`.

### 4.10 Crank-pin journal film frame and Reynolds boundary conditions

The coordinate frame, film kinematics, and boundary/cavitation conditions for the
**dynamic (D2) journal film solver and the rotor equation of motion**. It is chosen
so the **1D short-bearing** solver validates against the §4.9 Ocvirk closed form and
extends with no re-derivation to the **2D** film (circumferential + axial + tilt).
This section fixes conventions before any dynamic code is written; none is
implemented yet.

**Three distinct centres.** Do not conflate the macro orbit with the bearing
clearance -- they differ by ~300:1:

- **$O$ — cylinder centre = shaft rotation axis** (fixed global origin, §4.2).
- **$O_j$ — crank-pin centre = the journal** (the shaft eccentric; the inner, fast
  member). It is *prescribed*, orbiting $O$ at the crank-throw radius $e = 4.5$ mm:
  $O_j(\theta)=e(\sin\theta,\cos\theta)$.
- **$O_b = O_r$ — bore centre = rotor-bore centre = rotor centre = the sleeve** (the
  outer member, swings only). It is the **moving body in the EOM**; nominally
  coincident with $O_j$, offset only by the clearance-scale bearing eccentricity.

So neither bearing centre is the shaft axis: both orbit $O$ at $e=4.5$ mm. The film
is drawn in `geometry/journal_film_coordinates.png` (full rotor + journal zoom).

**Film frame $F_j$.** A *local* frame that co-orbits with the crank throw. Origin at
the bore centre $O_b$ (the rotor-bore/sleeve centre, the moving member); axes stay
**parallel to the global frame** §4.2 ($\hat x$ inlet, $\hat y$ vane; non-spinning):

- circumferential angle $\alpha$ from $+\hat x$, counter-clockwise about $+\hat z$,
  $\alpha\in[0,2\pi)$; arc coordinate $\xi = r_j\,\alpha$,
- axial coordinate $\zeta = z\in[-L_j/2,\,+L_j/2]$ (journal axis $\parallel\hat z$).

**Eccentricity state ( = the EOM state).** The **Cartesian** bearing eccentricity is
the rotor-centre offset *from the crank pin*,

$$\mathbf{e}=(e_x,e_y)=O_b-O_j,\qquad |\mathbf{e}|\le c_j\approx15\ \mu\text{m}\ \ (\ll
e=4.5\text{ mm}),$$

with $\varepsilon=|\mathbf{e}|/c_j$, $\psi=\operatorname{atan2}(e_y,e_x)$. The
Cartesian form (not $\varepsilon,\psi$) has **no attitude-angle singularity at
$\varepsilon\!\to\!0$** and is exactly the rotor displacement the equation of motion
integrates: $\mathbf{e}$ is the state and the absolute rotor centre is
$O_b=O_j(\theta)+\mathbf{e}$ (prescribed orbit plus the clearance-scale eccentricity).

**Film thickness.** With $\mathbf{e}=O_b-O_j$ (sleeve centre relative to journal),
$h(\alpha)=c_j+(e_x\cos\alpha+e_y\sin\alpha)$, minimum opposite $\mathbf{e}$. With
journal-axis **misalignment** (tilt slope $\boldsymbol\tau=(\tau_x,\tau_y)$, rad) the
offset varies linearly along $\zeta$, $\mathbf{e}(\zeta)=\mathbf{e}_0+\zeta
\boldsymbol\tau$, so

$$h(\alpha,\zeta)=c_j+\big[(e_x+\zeta\tau_x)\cos\alpha+(e_y+\zeta\tau_y)\sin\alpha\big].$$

Tilt is deferred (the first rungs keep $\boldsymbol\tau=0$); the term is written now
so the 2D extension needs no re-derivation. Axial section:
`geometry/journal_film_axial.png`.

**Reynolds equation.** On $(\xi,\zeta)$, incompressible and iso-viscous,

$$\frac{\partial}{\partial\xi}\!\Big(\frac{h^3}{12\mu}\frac{\partial p}{\partial\xi}
\Big)+\frac{\partial}{\partial\zeta}\!\Big(\frac{h^3}{12\mu}\frac{\partial p}
{\partial\zeta}\Big)=\frac{U}{2}\frac{\partial h}{\partial\xi}+\frac{\partial h}
{\partial t},$$

with entrainment $U=r_j\,\omega\,(1+\mathrm{d}\phi/\mathrm{d}\theta)$ (crank pin at
$\omega$, rotor bore at $\omega\,\mathrm{d}\phi/\mathrm{d}\theta$; consistent with the
§4.9 mean speed $\Omega=U/2r_j$) and squeeze $\partial h/\partial t$ from
$(\dot e_x,\dot e_y)$ [and $\dot{\boldsymbol\tau}$ under tilt].

**Boundary and cavitation conditions.** The crank-pin journal is bathed in the recess
environment at a **uniform $p_{rec}=4.0$ MPa** (§3.3): the suction/discharge split is
**circumferential/radial, not axial**, so the low- and high-pressure chambers do not
straddle the journal along $\zeta$ and the two axial ends carry the **same** pressure,

$$p(\alpha,\pm L_j/2)=p_{rec}\quad\text{(symmetric Dirichlet).}$$

Film rupture in the refrigerant-saturated POE oil is taken at the surrounding
pressure, cavitation floor $p_{cav}=p_{rec}=4.0$ MPa: **Swift–Stieber (Reynolds)** for
the numerical solver, **$\pi$-film (Gümbel)** for the §4.9 closed-form check.

**Offset invariance.** Reynolds is linear in $p$; shifting the axial-end Dirichlet
value and the cavitation floor by the **same** constant shifts the whole field by
that constant, and the net journal force $\oint p\,\mathbf{n}\,\mathrm{d}A$ is
**unchanged** (a uniform pressure integrates to zero over the closed journal
surface). Hence the uniform $4.0$ MPa environment leaves the §4.6/§4.9 load
(peak $\varepsilon\approx0.71$) unchanged; the absolute level only sets the
cavitation location and the axial leakage, **not the load capacity**.

**Force (and 2D moments).** The journal reaction is the pressure integral over the
wetted surface,

$$F_x=-\iint p\cos\alpha\;r_j\,\mathrm{d}\alpha\,\mathrm{d}\zeta,\qquad
F_y=-\iint p\sin\alpha\;r_j\,\mathrm{d}\alpha\,\mathrm{d}\zeta,$$

and under tilt the moments $M_x,M_y=\iint p\,\zeta\,(\cdots)$ close the misaligned
balance. These feed the rotor EOM $\mathbf{M}\ddot{\mathbf q}=\mathbf F_{gas}+\mathbf
F_{film}+\cdots$ (§4.4), with generalized coordinates $(e_x,e_y)$ [+ bush, per the D2
scope].

**1D → 2D ladder.** (i) *1D short bearing:* drop $\partial/\partial\xi$, $h=h(\alpha)$
— the axial solution is parabolic in $\zeta$; integrate the force and **validate
against the §4.9 Ocvirk closed form** (this checks the numerical machinery,
cavitation clamp, and force integration). (ii) *2D:* retain both terms on an
$(\alpha,\zeta)$ grid (finite bearing, Swift–Stieber) whose short-$L/D$ limit must
reproduce Ocvirk; then activate tilt ($\boldsymbol\tau\neq0$, $z$-dependent $h$).
Figures `geometry/{journal_film_coordinates, journal_film_axial}.png`.

### 4.11 Swing-bush film frames and vane-referenced attitudes

The swing bush (§3.6) runs two films per piece -- a **curved** film (bush OD in the
rotor groove) and a **flat** film (bush flat against the fixed vane). For the dynamic
(D2) rung the bush becomes a moving body, so both films need force-carrying frames
like the crank-pin journal (§4.10), and the rotor/bush **attitudes must be referenced
to the vane** (the true ground) rather than to the bush, whose fixed-attitude
assumption (§3.6) no longer holds. This section fixes those conventions; no dynamic
code uses them yet. Figures `geometry/{bush_film_coordinates,
bush_attitude_reference}.png`.

**Curved film $F_c$ (partial journal).** The bush-piece OD (convex, $r_b = 7.97$ mm)
runs in the rotor groove (concave, $r_{cut} = 8$ mm) over the cylindrical arc -- a
partial journal with clearance $c_c\approx10\,\mu$m. Origin at the groove centre
$O_g$; circumferential angle $\beta$ over the arc; axial $\zeta\in[-H/2,H/2]$. The
bearing eccentricity is $\mathbf{e}_c=O_p-O_g$ (bush-piece centre relative to groove
centre) and $h_c(\beta)=c_c-\mathbf{e}_c\!\cdot\!\hat r$. The film shears at the rotor
swing **relative to the bush** (below).

**Flat film $F_f$ (slider).** The bush flat runs on the fixed vane flank -- a
Cartesian slider. Coordinate $s$ along the vane ($0..L_f$, the contact length); axial
$\zeta$. Gap $h_f(s)=c_f-\Delta x+s\,\gamma$ (uniform $c_f\approx10\,\mu$m plus a
wedge $\gamma$ if the bush tilts). The film shears at the bush translation along the
vane. Both films carry the same Reynolds equation and boundary conditions as $F_j$:
chamber pressure at the OD/mouth end, $p_{rec}=4$ MPa at the inner end (§4.10).

**Vane-referenced attitudes.** With the vane as ground ($+\hat y$, §4.2):

- **rotor attitude** $\phi_r=\phi-90^\circ$ (the §4.2 orientation less the vane
  direction) -- *prescribed*, the $\pm10.37^\circ$ swing;
- **bush attitude** $\phi_b$ relative to the vane -- formerly constrained to $0$
  (prismatic joint), now a small **D2 DOF** within the flat-film clearance (an
  **in-plane** rotation about $\hat z$, equal to the flat-film wedge $\gamma$;
  distinct from the *out-of-plane* journal misalignment of §4.10).

The film-driving motions become **differences**: the curved film shears at
$U_c=r_b(\dot\phi_r-\dot\phi_b)$ (rotor relative to the *moving* bush), and the flat
film at the bush translation with $\phi_b$ as its wedge. **Setting $\phi_b=0$ recovers
the §3.6/§4.7 formulas exactly** (rotor swing against a fixed-attitude bush), so this
is a clean generalisation.

**Film pressure -- axial-uniform 1-D reduction.** The two bush films use the
**opposite** 1-D reduction of the Reynolds equation to the crank-pin journal. The
journal (§4.9/§4.12) is short ($L/D\approx0.74<1$), so the pressure escapes
**axially** and the axial (Ocvirk) flow is kept. The bush films are the reverse: the
curved arc ($L/D\approx1.3$) and the flat pad (axial height $H\gg$ sliding length
$L_f$) are **long** in the film-line direction, so the pressure escapes **along the
film line**, not axially. We therefore assume the pressure is **uniform across the
axial width** ($\partial p/\partial\zeta\approx0$), solve the 1-D Reynolds equation
*along the film line*, and multiply by $H$ for the force. Writing the line coordinate
as arc-length $x=r_b\beta$ (curved) or $s$ (flat), the equation
$\frac{d}{dx}\!\left(\frac{h^3}{12\mu}\frac{dp}{dx}\right)=\bar u\,\frac{dh}{dx}+\frac{\partial h}{\partial t}$
becomes, with $p=0$ at both open ends,

$$
\frac{d}{d\beta}\!\Big(h_c^3\frac{dp}{d\beta}\Big)=12\mu r_b^2\,S_c(\beta),\quad
S_c=\Omega\,(e_x\sin\beta-e_y\cos\beta)-(\dot e_x\cos\beta+\dot e_y\sin\beta);
$$
$$
\frac{d}{ds}\!\Big(h_f^3\frac{dp}{ds}\Big)=12\mu\,S_f(s),\quad
S_f=\tfrac{U}{2}\gamma-\dot\Delta+s\,\dot\gamma,
$$

with $\Omega=\bar u/r_b$ the mean surface (entrainment) angular speed. Each is
discretised in conservative flux form and solved as a tridiagonal system
(`mochi.line_reynolds`), with **Gumbel (half-Sommerfeld)** cavitation ($p\ge0$).
The curved-film pressure is radial and passes through the arc's centre of curvature
$O_p$, so it gives a **force but no moment** about $O_p$; the piece attitude $\phi_b$
is driven by the **flat-film moment** $M=H\!\int p\,s\,ds$. Validation of the flat pad
in **both source modes**: the converging wedge reproduces Reynolds' analytic
fixed-incline-slider load, and the parallel squeeze reproduces
$W=\mu\dot\Delta L_f^3 H/h^3$ (both to $\lesssim10^{-4}$); a shallow arc reduces to the
same slider (`mochi.arc_film`, `mochi.slider_film`). Modules `mochi.line_reynolds`,
`mochi.arc_film`, `mochi.slider_film`.

> **Long vs short — the two journal limits.** The curved bush film and the crank-pin
> journal (§4.10) are journals at *opposite* aspect ratios and get opposite reductions.
> The **bush curved film** ($L/D>1$) is **long-bearing** (circumferential escape,
> $dp/dz\approx0$) → benchmarked against the infinitely-long Sommerfeld form above. The
> **crank-pin journal** ($L/D=21/28.4=0.74<1$) is **short-bearing** (axial escape,
> parabolic in $\zeta$) → benchmarked against the Ocvirk closed form (§4.9, §4.12,
> `reynolds_1d` vs Ocvirk to $\sim10^{-5}$ at the crank-pin geometry). Both journals
> are **$\pi$-film (half-Sommerfeld)**: the crank-pin's symmetric $4$ MPa is a full
> *axial-end* Dirichlet, but circumferentially the film still ruptures at the recess
> floor, so it cavitates like the bush — a full-Sommerfeld (non-cavitating) load would
> over-predict it (~2–3×), and the long-bearing form over-predicts a further ~2× on top
> at this $L/D$.

**Long-bearing (Sommerfeld) benchmark for the curved film.** The curved film is a
journal in the *opposite* aspect-ratio limit to the crank-pin short bearing (§4.9,
§4.12): here $L/D>1$, the long-bearing limit. Its analytical counterpart is the
**infinitely-long (Sommerfeld) journal**, whose exact steady pressure with
$h=c(1+\varepsilon\cos\theta)$ is
$p(\theta)=\frac{6\mu U r_b}{c^2}\frac{\varepsilon\sin\theta\,(2+\varepsilon\cos\theta)}{(2+\varepsilon^2)(1+\varepsilon\cos\theta)^2}$
($U=2r_b\Omega$). This pressure is **zero at both the maximum- and minimum-film
points** ($\theta=0,\pi$), so on the converging half $[0,\pi]$ it is exactly the
$p=0$-ended Dirichlet problem `arc_film` solves, with $p>0$ throughout (Gumbel
inactive). `arc_film` therefore **reproduces the closed-form Sommerfeld force to its
discretisation error** ($\sim10^{-7}$ relative over $\varepsilon=0.2$–$0.9$, force
components and peak pressure), and the Gumbel half-load carries the pure-number
attitude invariant $\tan\phi=\pi\sqrt{1-\varepsilon^2}/(2\varepsilon)$ (full-Sommerfeld
gives zero radial load, $\phi=90^\circ$) — an $\mu,U,r_b,L,c$-independent check. This
is the long-bearing analogue of the §4.12 short-bearing (Ocvirk) cross-check. Model
`mochi.long_bearing` (`long_bearing_load`, `sommerfeld_pressure`); tests
`tests/test_long_bearing.py` (22).

**Gas boundary conditions, and cavitation on the *total* field.** Both bush films open
to gas at their ends, so the imposed end pressures add the source-free Poiseuille
field $p_{gas}(\xi)=p_0+(p_L-p_0)\,C(\xi)/C(L)$, $C(\xi)=\int_0^\xi h^{-3}d\xi'$
(`mochi.line_reynolds.poiseuille_bias_pressure`). Reynolds is **linear in $p$**, so the
physical field is the superposition $p=p_{hyd}+p_{gas}$ — but **cavitation is a
constraint on that sum, not on $p_{hyd}$ alone**. An earlier revision disabled the
clamp whenever end pressures were present (on the assumption that gas flooding
suppresses cavitation); evaluated at the real operating state the total reached
$-11.6$ MPa **absolute** on the diverging half-stroke of the reciprocating pad, worth
$-2310$ N of impossible suction. The solvers therefore compute $p_{hyd}$ unclamped,
superpose $p_{gas}$, and clamp the sum at `cavitation_pressure_pa` (expressed in the
same gauge as the end pressures; default $0$ recovers classic Gumbel exactly). The
floor follows one rule for all three films — *a film cannot fall below the lowest
pressure it is connected to* — which for the crank-pin journal (symmetric $4.0$ MPa
ends) degenerates to the existing $p_{cav}=p_{rec}$ of §4.10, and for the bush films
(asymmetric recess $\leftrightarrow$ chamber ends) gives $p_{cav}=\min(0,p_{ch}-p_{rec})$.
The exact floor between the chamber pressure and absolute zero depends on the
oil/refrigerant solubility and is **not** modelled (§9).

**Curved-film sealing land and its ends.** The curved film does not span the whole
piece arc: per §3.3 the groove wall runs from the mouth blend to the channel opening,
an **$86.5^\circ$ land** ($41.5^\circ$ mouth side $+\,45.0^\circ$ channel side about the
piece's arc centre); beyond it the arc faces gas, not oil. Those two ends open to
*different* gases — the rotor **mouth** (the piece's own chamber) and the **recess
channel** at $p_{rec}=4.0$ MPa — so the curved film carries a genuine pressure drop,
mirrored between the two pieces. Modelling it yields $+367$ N per piece **away from the
vane**, largely cancelling the flat film's $-400$ N gas squeeze; omitting it (the earlier
state) left that squeeze unbalanced.

**Mixed lubrication (Greenwood–Tripp).** Where the films run into the roughness
($\Lambda=h_{min}/\sigma\lesssim3$) the asperities share the load with the film. The
contact pressure is the standard elastic Greenwood–Tripp form
$p_{asp}=\tfrac{16\sqrt2}{15}\pi(\eta\beta\sigma)^2E'\sqrt{\sigma/\beta}\,F_{5/2}(h/\sigma)$
(`mochi.asperity_contact`), driven by measurable surface parameters rather than an
arbitrary pressure scale, and applied to **both** bush films by load sharing. Friction
is then boundary ($\mu_b$ on the asperity load) plus viscous shear, reported separately
on `RotorBushOrbit`. Literature-typical ground-steel values ($\beta=10\,\mu$m,
$\eta=10^{10}\,$m$^{-2}$, $E'=115$ GPa, $\eta\beta\sigma\approx0.05$) are placeholders
pending a surface measurement (§9); the flat-film loss varies over 66–91 W across
$\beta=5$–$40\,\mu$m. Tests `tests/test_asperity_contact.py` (6).

**Why the curved film runs thin — a load-path consequence, not an artefact.** The rotor
attitude EOM (§4.14) closes with $M_{gas}=+22.2$ vs $M_{bush,curved}=-21.3$ N·m (96 %
cancelled, residual $0.019$ N·m $=2\%$ of the reference inertia torque): the curved film
**is** the structural path that carries the gas moment. Reacting it at the groove lever
$l_g=25$ mm demands $21.28/0.025=851$ N of curved-film force. A direct load-capacity
sweep of `arc_film` at the swing entrainment ($\Omega\approx34$ rad/s) gives $725$ N at
$\varepsilon=0.98$ and $1527$ N at $0.99$, so $851$ N sits at $\varepsilon\approx0.985$,
i.e. a $\approx0.45\,\mu$m film — which is what the coupled orbit reports ($0.58\,\mu$m).
The thin curved film and its $\varepsilon$ against the safety rail are therefore
**required by the load path**, not a numerical artefact; adding elastic (Winkler)
compliance to the arc moved $\varepsilon/\varepsilon_{max}$ only $1.328\to1.330$ (the
eccentricity is set by geometry — the piece is pinned to the vane by the flat film while
the groove rides the rotor) and was reverted. Full record:
`docs/bush_film_revision_2026-08.md`.

### 4.12 Numerical short-bearing (1-D) Reynolds solver

Step (i) of the §4.10 1D$\to$2D ladder, and the project's **first PDE-style film
solver**: it integrates the **axial** (short-bearing) Reynolds equation numerically
and **validates against the §4.9 Ocvirk closed form**. Dropping the circumferential
pressure-flow term leaves, at each circumferential angle $\beta$, a 1-D axial equation

$$\frac{\partial}{\partial z}\!\Big(\frac{h^3}{12\mu}\frac{\partial p}{\partial z}
\Big)=S(\beta),\quad S=c\,(\dot\varepsilon\cos\beta-\varepsilon\Omega\sin\beta),\quad
p(\pm L/2)=0,$$

with $h(\beta)=c(1+\varepsilon\cos\beta)$ (§4.9 convention). A **tridiagonal finite
difference** solves the axial profile, **Gumbel (half-Sommerfeld)** cavitation clamps
$p\ge0$ to match the §4.9 $\pi$-film, and the pressure is integrated over the wetted
surface to the force $(F_e,F_t)$.

**Validation.** For pure rotation the numeric force reproduces the Ocvirk closed form
to $\sim10^{-4}$ (relative) across $\varepsilon=0.1$–$0.9$ -- the axial solution is an
exact parabola, so the residual is only the circumferential quadrature. The value of
this rung is not a new number but the **numerical machinery** -- tridiagonal axial
solve, cavitation clamp, force integration -- that the later 2-D (circumferential +
axial + tilt) Reynolds solver reuses; here it is checked against a closed form.
**Honest finding:** with a squeeze rate ($\dot\varepsilon\neq0$) the numeric and
closed form differ by $\sim2\%$, because the §4.9 closed form **superposes** two
$\pi$-films (rotation + squeeze) while the solver clamps the **true (shifted)**
cavitation region; neither is exact without a mass-conserving (JFO) cavitation model
(a later refinement). Implemented in `mochi.reynolds_1d` (`solve_short_bearing_1d` →
`ShortBearingSolution`); figure `bearing_load/reynolds_1d_validation.png`; tests
`tests/test_reynolds_1d.py` (15).

### 4.13 Rotor lateral equation of motion (first time integration)

The project's **first real mechanical time integration**. The rotor bore no longer
sits exactly on the crank pin but moves within the crank-pin clearance under gas
load, oil film, and inertia. At constant drive speed the state is the **bearing
eccentricity** $\mathbf{e}=O_b-O_j$ (§4.10); with $O_b=O_j(\theta)+\mathbf{e}$ and
$O_j=e_\text{throw}(\sin\theta,\cos\theta)$, Newton's second law gives

$$m_r\,\ddot{\mathbf{e}} = \mathbf{F}_{gas}(\theta) + \mathbf{F}_{film}(\mathbf{e},
\dot{\mathbf{e}}) + m_r\,\omega^2\,O_j(\theta),$$

the last term the orbital **centrifugal load** $m_r e_\text{throw}\omega^2$ the
quasi-static model dropped (§4.8). $\mathbf{F}_{gas}$ is the mouth-aware gas force
(§4.5); $\mathbf{F}_{film}$ is the Ocvirk reaction (§4.9) at the instantaneous
$(\varepsilon,\dot\varepsilon,\psi)$ -- the whirl rate $\dot\psi$ reduces the
entrainment and the squeeze term $\dot\varepsilon$ is now **driven live**. The gas
force is prescribed on a periodic grid (the small seal-over seam bridged linearly);
the stiff, overdamped system is integrated with an implicit method
(`scipy.integrate.solve_ivp`, BDF) from the quasi-static equilibrium to a steady
periodic orbit.

**Result -- the squeeze film limits the eccentricity.** The centrifugal (inertia)
load $m_r e_\text{throw}\omega^2 \approx 44$ N is only $\sim2\%$ of the $\approx2.5$
kN gas reaction, so inertia is a small correction. The **dominant dynamic effect is
the squeeze-film lag**: the film's squeeze response has a time constant that is a
fixed fraction of the cycle ($\tau/T \approx 0.26$, independent of speed, mass, and
clearance -- a geometric function of $\varepsilon$), so the rotor cannot follow the
fast quasi-static swing. It **attenuates and phase-lags** it -- **peak eccentricity
$\approx0.50$ (vs the quasi-static $0.71$) and a larger minimum film $\approx7.4$
$\mu$m** (vs $4.4$). The whirl orbit is a small bounded loop well inside the
clearance. So the **quasi-static model (§4.9) overestimates the peak eccentricity**;
the journal runs safer than it predicted.

**Journal friction from the dynamic eccentricity.** Feeding the dynamic
$\varepsilon(\theta)$ into the eccentric-friction law (§4.9) gives the realistic
cycle-mean journal loss: **$\approx10.2$ W** — between the concentric Petroff
($\approx9.1$ W, §4.7) and the quasi-static eccentric estimate ($\approx10.7$ W,
§4.9), and ~4 % below the latter. The dynamic $\varepsilon$ swing is *flatter*
(higher trough though lower peak), and the friction $\propto1/\sqrt{1-\varepsilon^2}$
is convex, so the dynamic mean lands just under the wider quasi-static swing. (The
swing-bush film friction, §3.6, is a separate ~0.2 W term; its dynamic value needs
the full bush multi-body rung.)

**Validation.** (i) *Static balance:* at the quasi-static equilibrium with
$\dot{\mathbf{e}}=0$ the film reaction equals $-\mathbf{F}_{gas}$ exactly (to
$10^{-6}$). (ii) *Reduction check:* under a **constant (frozen) load** the orbit
relaxes exactly to the §4.9 quasi-static equilibrium -- so the operating-speed
deviation is the genuine dynamic (squeeze) effect, not numerics. (iii) The steady
orbit is periodic. Confirmed real rotor mass $m_r=0.275$ kg (§4.8 note); editable.
Implemented in `mochi.rotor_dynamics` (`integrate_rotor_orbit` → `RotorOrbit`);
figure `bearing_load/rotor_orbit.png`; tests `tests/test_rotor_dynamics.py` (9).

### 4.14 Nine-DOF rotor + two-piece swing-bush multibody dynamics

The §4.13 rotor rung froze the swing bush on the vane axis. This rung frees it: the
**two bush pieces become independent rigid bodies**, coupled to the rotor through the
curved and flat oil films (§4.11) rebuilt in the axial-uniform 1-D form (§4.11,
`arc_film`/`slider_film`). Nine DOF, all modelled as **deviations about the prescribed
kinematics** (§3.1) -- the rigid-film limit, exactly as $\mathbf{e}_j=O_b-O_j$ already
is in §4.13:

| DOF | meaning | reference (deviation 0) |
|---|---|---|
| $\mathbf{e}_j$ (2) | rotor bore $-$ crank pin, $O_b-O_j$ | $0$ |
| $\delta\phi_r$ (1) | rotor attitude, $\phi_\text{orient}=\phi_\text{orient,ref}+\delta\phi_r$ | $0$ |
| $O_{p,k}$ (2×2) | each piece centre (global) | $O_{p,k,\text{ref}}$ |
| $\phi_k$ (2×1) | each piece attitude (rel. vane) | $0$ |

integrated as an 18-state stiff ODE (BDF), keeping the last of ~3--4 revolutions.
Zeroing all deviations recovers the prescribed pose and the §4.11 $10\,\mu$m films.

**Two frames, deliberately.** The groove rides the rotor, so $O_g=O_j+\mathbf{e}_j+
\ell_g\hat n$ and the curved-film eccentricity $\mathbf{e}_c=O_p-O_g$ use the
*swinging* rotor frame $\hat n=(\cos\phi_\text{orient},\sin\phi_\text{orient})$,
$\hat\tau=(-\sin,\cos)$ ($\ell_g=25$ mm the groove lever). The vane is *fixed* in the
world, so the flat film against it is taken in the **world frame**: the piece slides
along world $y$ (slide speed $=\dot y_k$, matching `flat_slide_velocity` at the
reference) and approaches the vane along world $x$ (approach $\delta_k=\text{shift}-
\text{side}_k\,x_k$). Putting the flat film in the swinging frame would misstate the
slide speed by the $\sim10^\circ$ swing (at $\theta=90^\circ$ the world-$y$ slide is
$0.85$ m/s, the $\hat\tau$-component only $0.15$). The curved-film clearance fed to
`arc_film` is the **concentric** radial gap $r_{cut}-r_b=30\,\mu$m; the reference
$20\,\mu$m piece shift toward the piece's own side then sets the minimum curved film
to $30-20=10\,\mu$m (the §4.11 operating film) -- these are two conventions for the
same geometry.

**The nine equations of motion.** Rotor (3): the §4.13 lateral EOM with the curved
reactions $-\mathbf{F}_{c,k}$ added at $O_g$, plus a rotor-attitude equation

$$m_r\ddot{\mathbf{e}}_j = \mathbf{F}_{gas}+\mathbf{F}_{journal}+\textstyle\sum_k(-\mathbf{F}_{c,k})+\mathbf{F}_{seal}+m_r\omega^2O_j,\qquad
I_r\,\delta\ddot\phi_r = M_{gas}+M_{fric}+\textstyle\sum_k[(O_g-O_b)\times(-\mathbf{F}_{c,k})]_z + M_{seal} - I_r\ddot\phi_\text{orient,ref}.$$

The reference-inertia term $-I_r\ddot\phi_\text{orient,ref}$ ($\lesssim1.4$ N·m) is the
angular analogue of the lateral $m_r\omega^2O_j$ -- $\phi_\text{orient}=\phi_\text{orient,ref}+\delta\phi_r$
and only the deviation is integrated. Each piece (3): $m_p\ddot O_{p,k}=\mathbf{F}_{c,k}+\mathbf{F}_{f,k}$
(curved + flat normal; the flat normal $\text{side}_k N_k\hat x$ pushes off the fixed
vane), $I_p\ddot\phi_k=M_{f,k}$ (the flat slider moment; the curved film is radial
through $O_p$ and adds no moment). $I_p$ is the raster CM inertia (§4.8) rescaled to
the confirmed $m_p=6.341$ g and shifted to the piece centre by parallel axis
($\approx2.7\times10^{-7}$ kg·m²). Curved-film entrainment $\Omega_c=\tfrac12(\dot\phi_\text{orient}+\dot\phi_k)$
(the mean groove + piece spin; the piece barely spins).

**The tenth force — rotor–cylinder sealing contact (default on).** The rotor OD also
seals against the cylinder bore (§4.15). By default (`seal_contact=True`) that contact
is solved **inside this 9-DOF system**, not in the reduced free-rotor model of §4.15: a
compliant Hertz line contact (`hertz_line_contact_force_n`) acts on the rotor along
$\hat u=O_b/|O_b|$ when the OD penetrates the bore ($|O_b|>e$), with boundary friction
tangential at $O_b+R_r\hat u$,
$$\mathbf{F}_{seal}=-N_c\,\hat u-\mu N_c\,\text{sgn}(v_\text{slide})\,\hat t,\qquad
M_{seal}=-R_r\,\mu N_c\,\text{sgn}(v_\text{slide})$$
(the normal is radial through $O_b$ and adds no moment; only the friction does, over the
lever $R_r$). This is the **full methodology**: $N_c$ is driven by the always-outward
centrifugal $m_r\omega^2 e$, so the seal engagement — and the loss — scale with
$\omega^2$; **any operating-speed sweep must keep the coupling on**, or the
(speed-dependent) rotor-cylinder loss is missing entirely. Setting `seal_contact=False`
isolates the bush films (no cylinder reaction), for the reduction cross-checks only.

**Result -- the swing bush is a heavily loaded element.** The rotor gas *moment* about
$O_b$ ($M_{gas}$ up to $\approx22$ N·m) can only be reacted against the fixed vane
**through the bush**: rotor $\to$ curved film $\to$ piece $\to$ flat film $\to$ vane.
Across the $\ell_g=25$ mm lever this is a $\sim$kN force the two films must carry, so
at peak load the **curved film runs very thin** (minimum $\approx0.4\,\mu$m, curved
eccentricity ratio $\approx0.99$) while the flat film thins to $\approx2\,\mu$m. The
rotor bore itself stays well inside its clearance (peak $\varepsilon_j\approx0.64$) and
all nine deviations stay small ($\delta\phi_r\approx1$ mrad, piece centres within
$\sim10\,\mu$m of reference), vindicating the perturbation premise. The **dynamic bush-film
friction $\approx0.6$ W** (Couette over the running films) is **the value of record**
for the swing-bush loss; it is a few times the quasi-static §3.6 estimate
($\approx0.2$ W, kept only as a comparison), reflecting the thinner loaded films -- the
Stage 6 deliverable. **This is a hydrodynamic-only first cut**
(gas-pressure film bias omitted, §9.3): the near-contact curved film is the honest
consequence and flags the bush as the critical wear/lubrication site, to revisit with
the gas bias and an EHL/contact treatment.

**Result — the coupled sealing load and its attribution.** With the seal on (default),
the $\sim6\,\mu$m free-orbit bore penetration collapses to a physical
$\approx2.4\,\mu$m Hertz deflection and the sealing load is **$N_c\approx296$ N mean,
$1455$ N peak** — about **twice** the free-rotor reduced model ($140$ N, §4.15). The
cycle-mean *radial* force balance (outward $+$, N) attributes it cleanly:
gas $-223$, journal $+709$, bush $-234$, seal $-296$, centrifugal $+44$ (sum $\approx0$).
The **bush reaction is inward** — it does *not* push the rotor onto the bore. Rather,
tying the rotor to the **fixed vane loads the crank-pin journal** ($+323\to+709$ N
outward vs the free rotor), and that journal force, net of the inward bush, presses the
rotor $\sim1\,\mu$m deeper into the cylinder. So the free-floating model of §4.15
**under-predicts the seal loss by $2\times$ precisely because it omits the vane
constraint**, and the coupled value of record is $\dot W_{r\text{-}c}\approx11$ W
($\mu=0.1$, Ra 0.3 µm). The seal-friction torque ($R_r\mu N_c\approx1.0$ N·m) is
reacted through the bush without distorting the attitude ($\delta\phi_r$ stays
$\approx1$ mrad, unchanged from the seal-off run), and $N_c$ is converged (identical
over revolutions 4–6).

**Numerics.** The squeeze films are extremely stiff (a fast $\sim3\times10^{-8}$ s
overdamped eigenvalue: $\dot{\mathbf{e}}_c\!=\!0.01$ m/s already gives $\sim$kN). BDF
with a uniformly tight `atol` pins to that timescale forever; a **loose absolute
tolerance on the nine velocity states** (only a means to the positions) lets it step
over the decayed fast mode while the positions -- the films we report -- stay tight.
The pieces start at the reference groove velocity so the squeeze films begin unloaded.

**Validation.** (i) *Film reduction:* at zero deviation the curved eccentricity is the
$20\,\mu$m shift into the $30\,\mu$m concentric gap and the flat approach is zero, so
both films are $10\,\mu$m (§4.11, `film_thicknesses_m`). (ii) *Perturbation:* the
integrated deviations stay small and the orbit is periodic. (iii) *Piece inertia:* the
parallel-axis inertia exceeds the CM value. (iv) *Seal coupling:* the radial force
balance closes, the penetration collapses to sub-3 µm, and $N_c$ is revolution-converged.
(v) *1-D Reynolds vs analytical:* each film's 1-D reduction reproduces its closed-form
reference to $\le10^{-4}$ — bush curved vs the long-bearing Sommerfeld load, bush flat
vs the Reynolds fixed-incline slider, journal vs Ocvirk (figures
`bush_film/reynolds_{curved_vs_long_bearing, flat_vs_incline_slider, journal_vs_ocvirk}.png`,
one per film). Implemented in `mochi.rotor_bush_dynamics`
(`integrate_rotor_bush_orbit` → `RotorBushOrbit`, `seal_contact=True` by default); tests
`tests/test_rotor_bush_dynamics.py` (13); figures
`bush_film/film_clearance_{journal, bush_curved, bush_flat}.png`
(each film's clearance at the most-loaded crank angle),
`bush_film/bush_{curved, flat}_clearance_vs_crank.png` (each bush film's minimum clearance over
the whole cycle for both pieces -- the IN and OUT pieces alternate as the gas moment reverses,
the curved IN film reaching $\approx0.7\,\mu$m),
`assembly/{layout, bush_clearance, journal_clearance}_42deg.png` (the journal, vane, rotor
and two bush pieces at the most-loaded angle -- the full assembly at true scale, and the
bush and journal clearances at one uniform exaggeration factor with the piece/groove centre
offsets, showing the IN piece run into near-contact), and
`bearing_load/friction_dynamic_vs_quasistatic.png` (per-film quasi-static vs dynamic:
bush ×3.3, journal ×1.0, seal ×3.7, total ×1.6). **Gas-pressure film boundary conditions — a finding (`gas_film_boundary=True`, opt-in;
the default keeps the validated hydrodynamic films).** The 1-D
film solvers accept chamber/crank gas pressures at the ends
(`arc_film_force`/`flat_slider_film` `pressure_start_pa`/`pressure_end_pa`), superposing
the Poiseuille gas-bias field (`poiseuille_bias_pressure`; Reynolds is linear in $p$, and
the high gas pressure suppresses cavitation so the hydrodynamic part runs full-Sommerfeld)
and reporting the throughflow leakage. In the coupled EOM the piece is immersed in the
**bore gas $P_\text{bore}=4$ MPa** (above both chambers), so the film pressures are
referenced to it and the uniform immersion cancels by the divergence theorem — the
**flat (vane-sealing) film** then carries the net gas load, running a bore→chamber drop
(outer end at the piece's chamber: IN=suction, OUT=compression; inner end the bore), while
the curved film stays in the bore region (no drop). **Result — a first-order, in fact
dominant, effect:** the bore gas presses each piece onto the vane, driving the flat film
from $\approx1.8\,\mu$m down to the $0.2\,\mu$m contact floor and raising the (Couette)
bush friction from $0.67$ to $\approx5.5$ W (×8); the rotor orbit and the rotor-cylinder
seal are unchanged. Because the flat film reaches metal contact, this pure-hydrodynamic
value is **contact-clamp-limited** — the vane-bush flat interface now needs an EHL/contact
treatment, exactly as the rotor-cylinder seal did (§4.15). Open items §9: that flat-film
EHL/contact, the curved-film gas mapping, curved entrainment sign cross-check.

### 4.15 Rotor-cylinder sealing contact — sliding kinematics and friction

The rotor OD seals against the cylinder bore along a line near the vane. Two facts fix
the friction character. **(1) The geometry is exactly critical:** the radial clearance
$R_{cyl}-R_r=(77-68)/2=4.5$ mm **equals the throw** $e=4.5$ mm, so in the rigid limit
the rotor centre orbits at $e$ about $O$ and the OD just kisses the bore at every crank
angle (the rolling-piston condition, penetration = gap = 0). **(2) The rotor swings,
it does not spin:** tied to the fixed vane through the swing bush, its orientation
oscillates $\pm10.37^\circ$ with **net-zero rotation per revolution**. A body that does
not rotate cannot roll, so unlike a free rolling piston the sealing line **slides**:

$$v_\text{slide}(\theta)=-e\,\omega+\omega\,\frac{d\phi_\text{orient}}{d\theta}\,R_r,$$

the orbital term $-e\omega=-0.85$ m/s (constant) plus the swing modulation
($\pm1.15$ m/s, zero at the $\theta=90^\circ/270^\circ$ orientation extrema, maximal at
TDC/BDC). Over a cycle $|v_\text{slide}|$ averages **0.95 m/s** (peak 2.0) — the seal
genuinely rubs, so this is potentially a *primary* loss, not a small rolling
correction.

**Boundary friction and the missing $N_c$.** In boundary lubrication the loss is
$P=\langle\mu\,N_c\,|v_\text{slide}|\rangle$. Two of three ingredients are in hand:
$v_\text{slide}$ (above) and $\mu$ (steel-on-steel, refrigerant-oil boundary
$\approx0.1$; the mixed-lubrication asperity value from compressor studies is
$0.04$–$0.08$). The **contact normal force $N_c$ is not**: it is statically
indeterminate — the crank-pin bearing, this contact, and the bush/vane over-constrain
the rotor — and is the *same* unknown behind the **§4.13/§4.14 bore-penetration
artifact** (the free journal-clearance motion drives the OD $\sim6\,\mu$m into the bore
/ $\sim7\,\mu$m off it because no contact reaction resists it). A compliant
contact / EHL sealing-film model supplies $N_c$, removes the penetration, and closes
the loss in one step. Implemented: the sliding kinematics, a **quasi-static $N_c(\theta)$** estimate
(`contact_normal_force_n` — the rigid-contact / journal-tangential limit: the contact
reacts the net outward radial load, gas radial + centrifugal, clamped $\ge0$), the
boundary law, a mixed-EHL law, and the **self-consistent contact rung** below
(`mochi.rotor_cylinder`; tests, 15).

**Result — the two closures compared.** The estimated contact load is modest,
$N_c\approx50$ N mean / $190$ N peak (the crank-pin journal takes the gas load; the
contact carries only the $\sim44$ N centrifugal floor plus the outward gas residual,
active ~55 % of the cycle). At $\mu=0.1$, $|v|\approx0.95$ m/s:

| closure | assumption | $\dot W_{r\text{-}c}$ |
|---|---|---|
| **A** pure boundary | $\mu=0.06$ / $0.10$ | **1.9 / 3.1 W** |
| **B** mixed-EHL | $\sigma=0.1/0.2/0.3\,\mu$m | **0.7 / 1.4 / 2.4 W** |

The conformal contact ($R_{eq}\approx291$ mm) forms a partial EHL film — specific
thickness $\Lambda\approx1$–$2$ (mixed regime) — so **B is roughly half of A** at a
typical $\sigma=0.2\,\mu$m ($B/A\approx0.43$), and is strongly roughness-driven (smooth
→ near full-film $\sim0.7$ W, rough → boundary bound). This rigid estimate brackets **$\dot W_{r\text{-}c}\approx1$–3 W** but omits the
inertial loading; the self-consistent rung below is the definitive value.

**Self-consistent closure — reduced free-rotor cross-check (the rung).**
`integrate_sealing_contact` is a **4-state rotor-only** model (lateral bore motion +
journal film + contact + gas + centrifugal, *no bush/vane*); it is the cheap cross-check,
**not** the value of record. Coupling the same contact into the full 9-DOF system (§4.14,
`seal_contact=True`, the default) roughly **doubles** $N_c$ because the bush ties the
rotor to the fixed vane and loads the journal (attribution in §4.14). Take the coupled
$\dot W_{r\text{-}c}\approx11$ W as the value of record; the free-rotor value below
under-predicts by $\sim2\times$. `integrate_sealing_contact` replaces the
arbitrary penalty with the **physical compliant Hertz line contact**
(`hertz_line_contact_force_n`, Palmgren: $\delta[\text{mm}]=3.84\times10^{-5}\,
Q^{0.9}/L^{0.8}$ — $\sim0.22\,\mu$m deflection at $100$ N) added to the §4.13 rotor EOM
as a one-sided radial reaction + sliding friction, and integrates to the steady orbit.
The load then splits between the journal film and the contact **by their true
stiffnesses**, closing $N_c(\theta)$ with no free parameter. Result: the $\sim6\,\mu$m
free-orbit penetration collapses to the physical **$\sim1.4\,\mu$m Hertz deflection**
(and the seal lifts up to $\sim7.9\,\mu$m elsewhere — a leakage window), with **$N_c$
mean $\approx140$ N, peak $\approx780$ N**, contact ~1/3 of the cycle. The friction is
**higher than the rigid estimate** because the always-outward **centrifugal
($\sim44$ N)** and the inertial bore excursions press the rotor onto the bore harder
than the quasi-static gas balance alone:

$$\dot W_{r\text{-}c}\approx 6.7\ \text{W}\ \text{(free-rotor reduced model)}\quad\Rightarrow\quad
\boxed{\dot W_{r\text{-}c}\approx 11\ \text{W}}\ \text{(coupled 9-DOF, §4.14 — value of record)}\quad(\mu=0.1,\ \text{Ra}=0.3\,\mu\text{m}).$$

At the **Ra 0.3 µm design-standard finish** the composite RMS is $\sigma=\sqrt2\cdot
1.25\,\text{Ra}=0.53\,\mu$m, so $\Lambda\approx0.45$ — essentially **boundary**
lubrication, and mixed-EHL coincides with the boundary bound ($6.7$ W). (For a smoother
$\sigma=0.2\,\mu$m the film would be mixed, $\Lambda\approx1.2$, dropping mixed-EHL to
$\approx5.3$ W; roughness is the one free lever.) Insensitive to the contact-onset
smoothing. **So the rotor-cylinder seal is a first-order mechanical loss — in the coupled
model ($\approx11$ W) it is the *largest* single mechanical loss, above the crank-pin
journal ($8.9$ W) and an order above the bush (**$\approx0.6$ W**, the
Section 4.14 *dynamic* value of record; the $0.2$ W quasi-static §3.6 estimate is kept
only for comparison)** — even the free-rotor $6.7$ W already reaches the journal's order,
and the
earlier "tens of W" (pessimistic $N_c$) and "1–3 W" (rigid, no inertia) both bracketed
it loosely. Implemented in
`mochi.rotor_cylinder` (`integrate_sealing_contact` → `SealContactOrbit`;
`hertz_line_contact_force_n`); tests `tests/test_rotor_cylinder.py` (17).

**Full-EHL check — the seal does not float.** Two normal-law idealisations were tested.
Treating the roller–cylinder oil film as a **large-clearance journal bearing** (Ocvirk,
$c=e=4.5$ mm) is *invalid* — the thin-film assumption $c\ll R$ is violated and it
returns spurious $\sim$50 kN film forces. The correct hydrodynamic model is the
**EHL line contact** (conformal, $R_{eq}\approx291$ mm): at the self-consistent loads
the Dowson–Higginson minimum film is **$h_{\min}\approx0.24\,\mu$m**, so at the design
$\sigma=0.53\,\mu$m the film parameter is **$\Lambda\approx0.45$** — essentially
**boundary**, and even at a smooth $\sigma=0.2\,\mu$m only $\Lambda\approx1.2$ (mixed).
Either way the seal runs in asperity contact and does **not** float on a full film. For
this *hard* steel contact the normal load–deflection is therefore elastic (Hertz) — the
oil film ($\ll$ the $\sim1\,\mu$m elastic deflection) rides in the deflected contact and
sets only the friction regime. Hence the "full-EHL normal law" reduces to the standard
**EHL contact element** (Hertz elastic normal + EHL-film mixed friction) already
integrated here, confirming $\dot W_{r\text{-}c}\approx6.7$ W. `ehl_film_thickness_m`
and the orbit's `mean_film_parameter` report the regime. **Roughness convention:**
$\Lambda$ uses the composite RMS ($R_q$), not $R_a$ ($R_q\approx1.25\,R_a$; two surfaces
add in quadrature) — the Ra 0.3 µm spec gives $\sigma=0.53\,\mu$m. Open: coupling the
lift-off leakage window into §3.7.

**Toward the film-parameter (mixed-EHL) closure.** The contact is nearly conformal
(convex $R_r=34$ in concave $R_{cyl}=38.5$ → equivalent radius
$R_r R_{cyl}/(R_{cyl}-R_r)\approx291$ mm), so a hydrodynamic/EHL film can partly carry
the seal. The recommended closure is the **specific film thickness**
$\Lambda=h_{\min}/\sigma$ with $h_{\min}$ from the Dowson–Higginson line-contact EHL
formula (entrainment $\bar u=\tfrac12|v_\text{slide}|$, since the bore is stationary),
composite roughness $\sigma\approx0.1$–$0.3\,\mu$m, POE pressure–viscosity
$\alpha\approx13$–$20$ GPa$^{-1}$, steel $E'\approx226$ GPa. Regimes: $\Lambda>3$
full-film (viscous shear only), $1<\Lambda<3$ mixed, $\Lambda<1$ boundary. In the mixed
band the load and friction split between the film (Patir–Cheng average-flow) and the
asperities (Greenwood–Tripp, carrying $\mu_b\approx0.05$–$0.1$); the EHL film force
then *is* the compliant $N_c$, so this single model resolves the indeterminacy, the
penetration, and the friction together.

**References** (rolling-piston / swing-compressor friction & lubrication):
T. Yanagisawa & T. Shimizu, "Friction losses in rolling piston type rotary
compressors, I–III," *Int. J. Refrigeration* **8** (1985) — the classic contact-by-
contact friction model. "Theoretical model development and mixed-lubrication analyses
of rolling-piston rotary compressors: a review," *Lubricants* **12**(8):273 (2024) —
Greenwood–Tripp + Patir–Cheng methodology, reports $\mu\approx0.04$–$0.08$. For the
**swing** mechanism specifically (Daikin's integrated vane+roller on a swing bush, the
machine modelled here): "Development of swing compressor for alternative refrigerants,"
*Int. Compressor Eng. Conf. at Purdue*; and "Dynamic characteristics of a swing
compressor …," *Int. J. Refrigeration* (2019). Foundational tribology:
Dowson–Higginson (EHL film thickness), Greenwood–Tripp (asperity contact),
Patir–Cheng (average-flow mixed lubrication).

## 5. Kinematics and force interfaces

A compressor case must provide the geometry and prescribed driver state needed
by its force laws. If crank speed is prescribed,

$$
\theta_c(t) = \theta_{c,0} + \int_0^t \omega_c(\tau)\,d\tau.
$$

Constant speed is a special case, not a universal assumption. The crank-angle
origin and positive rotation must be documented with the first case.

Each force contribution should be a deterministic function with an interface
equivalent to

$$
\mathbf{F}_i =
\mathbf{F}_i(t,\mathbf{q},\mathbf{v},\text{case data}).
$$

It must document:

- the component exerting the force and the body receiving it;
- the frame in which inputs and outputs are expressed;
- required geometry and parameter units;
- validity range and singular or discontinuous conditions;
- source, derivation, or reference data;
- an independent test.

For unilateral contact, the eventual model must at least respect the gap and
normal-reaction conditions

$$
g(\mathbf{q}) \ge 0,\qquad \lambda_n \ge 0,\qquad
\lambda_n g(\mathbf{q}) = 0,
$$

before a penalty, complementarity, event-based, or other numerical treatment
is chosen. Friction requires a separate tangential law and sign convention.

## 6. Initial-value formulation

Define the first-order state

$$
\mathbf{y} =
\begin{bmatrix}\mathbf{q} \\ \mathbf{v}\end{bmatrix},\qquad
\dot{\mathbf{y}} =
\begin{bmatrix}
\mathbf{v} \\
\mathbf{F}_{total}(t,\mathbf{q},\mathbf{v})/m_b
\end{bmatrix}.
$$

For smooth forces, the first reference implementation may use SciPy's adaptive
`solve_ivp` interface with an explicit method and reported relative/absolute
tolerances. Output sampling must be independent of internal solver steps.

This choice must be reconsidered when discontinuous contact, stiffness, or
algebraic constraints are introduced. An adaptive solver reporting success is
not by itself evidence that the physical solution is correct.

State variables with different magnitudes need component-appropriate absolute
tolerances or documented nondimensionalization. A case must record the method,
tolerances, step limits, number of accepted/rejected steps when available, and
termination reason.

## 7. Required validation ladder

Implementations should progress through these tests in order:

1. **Zero force:** velocity is constant and position is linear in time.
2. **Constant force:** acceleration is constant and matches the analytic
   quadratic trajectory.
3. **Linear oscillator:** if a support term is implemented, compare frequency
   and damping decay with the analytic solution.
4. **Frame/sign cases:** simple axis-aligned loads produce motion in the stated
   positive direction; action/reaction pairs have opposite signs.
5. **Step/tolerance refinement:** quantities of interest approach a stable
   value as error tolerances tighten or maximum step decreases.
6. **Energy/work balance:** for conservative cases, total energy remains within
   the stated numerical error; for dissipative cases, the loss has the correct
   sign and magnitude.
7. **Compressor reference case:** compare against independently calculated,
   experimental, or trusted published results with provenance and acceptance
   limits.

Regression tolerances must reflect a measured numerical error. They must not be
chosen solely to accept the current output.

## 8. Case and result contract

Every reproducible case should eventually record:

- model version or Git commit;
- geometry and material parameters with units;
- active and inactive force models;
- initial state and simulated interval;
- prescribed speed/load histories and their interpolation rule;
- numerical method, tolerances, and step limits;
- requested outputs and output times.

Every result should include time, position, velocity, total force, each active
force contribution, and solver status. Results from failed or incomplete solves
must remain clearly marked and must not look like successful final states.

## 9. Open model decisions

The following questions block a physical production solver:

- [x] Define the future dynamic bodies and their interfaces: four rigid
  planar bodies — cylinder with integral stepped vane (ground), crankshaft,
  rotor, swing bush — with the gas as a massless pressure boundary condition
  (Section 4.1).
- [x] Define the prototype global frame, angle origin, and rotation direction.
- [ ] Confirm that prototype frame against the compressor CAD convention.
- [x] Decide the required degrees of freedom and constraints: net one degree
  of freedom (crank angle, prescribed) from 7 body DOFs minus 6 joint
  constraints (Section 4.3).
- [x] Document the supplied cylinder, rotor, circular-cutout, and vane dimensions.
- [x] Confirm the vane-tip relation and clearances against CAD: the vane is a
  fixed stepped part of the cylinder with its tip 25 mm from the bore; the
  2.4 mm ledges run inside the rotor's 14.7 mm recess channel with 0.75 mm
  axial clearance per side (Section 3.3). The former invisible reference-line
  rule was removed.
- [x] Obtain the swing-bush clearances: piece outer radius 7.970 mm with
  its flat face 3.990 mm from the piece center, each piece offset 0.020 mm
  toward its own side to give uniform 0.010 mm films; axial clearance
  0.0085 mm per side, bush height 20.983 mm (Section 3.3).
- [~] Revisit the recess/bush pressure boundary conditions; the recess
  spaces stay fixed at 4.0 MPa absolute meanwhile (Section 3.3). **Confirmed for the
  crank-pin journal environment** (2026-07-31): the rotor-centre/crank region is a
  uniform 4.0 MPa (circumferential/radial pressure split, not axial), so the journal
  film carries symmetric $4.0$ MPa axial-end BCs and $p_{cav}=4.0$ MPa (§4.10). The
  bush **film-end** BCs are now resolved from the §3.3 geometry (2026-08, §4.11): the
  flat pad runs recess$\to$chamber, and the curved film lives on the $86.5^\circ$
  sealing land whose ends open to the recess channel and the rotor mouth. Two pieces of
  this remain open: (i) the **cavitation floor** — the rule adopted is "no lower than the
  lowest connected pressure", but the true value lies between the chamber pressure and
  absolute zero and needs an oil/refrigerant **solubility** model (the bush loss is very
  sensitive to it: 48.7 W vs 6.2 W at 1800 rpm across that span); (ii) the
  Greenwood–Tripp surface parameters $\beta,\eta$, which need a **surface measurement**
  (66–91 W across the literature range). Whether the curved film really rides at
  $\Lambda\approx1$ in hardware is a wear question, not a modelling one — note the squeeze
  term is $\sim16\times$ the wedge at equal eccentricity, so the real machine may ride on
  squeeze more than this steady-orbit model does.
- [x] Confirm the crank-pin journal bearing dimensions against CAD. Radius
  $r_j = 14.2$ mm (bounded by the 17 mm bush-groove clearance) and length
  $L_j = 21$ mm are now taken from the CAD model (Section 4.7). The radial
  clearance is not a drawn dimension (nominal full-contact fit), so it stays an
  assumed $15\,\mu$m oil-film value, exposed as an editable constant in
  `mochi.journal_bearing`. The eccentricity-coupled short-bearing model (§4.9,
  `mochi.ocvirk_bearing`) now consumes it: at the assumed $c_j$ the peak reaction
  drives $\varepsilon \approx 0.71$ (minimum film $\approx 4.4\,\mu$m), so the
  clearance value directly sets the running film margin.
- [ ] Decide the oil-film leakage treatment across the recess and bush seals.
- [x] Resolve the mouth-lip/bore interference: the confirmed G1 lip
  construction (OD cut at 13°, R1.4 × 92.6°, 0.97 mm straight,
  R0.9 × 52.3°) touches the rotor OD only at its tangent point, so no bore
  interference exists at any crank angle; earlier semicircular-lip variants
  penetrated up to 0.32 mm and were discarded (Section 3.3).
- [ ] Identify all applied, gas, contact, friction, and fluid forces.
- [x] Confirm the refrigerant and the polytropic exponent: R410A with an
  effective exponent n = 1.07, fitted to the CoolProp HEOS isentrope from
  saturated vapour at 0.82 MPa to 3.24 MPa (1.064 saturated, 1.075 at 10 K
  superheat); the supplied ports match 0.85 degC evaporating and 52.5 degC
  condensing saturation temperatures.
- [~] Confirm the discharge-valve opening rise (currently 5 % of the
  discharge-port pressure) against actual valve or measurement data. Section 3.8
  now adds the reed valve as a parallel model (`mochi.reed_valve`): the 5% rise is
  the reed opening preload, and the finite geometric port area yields an emergent
  discharge overpressure (~+0.5 MPa, ~51 W loss at the supplied port size).
  Remaining: reed dynamics/back-flow and the true port size/discharge coefficient
  vs measurement.
- [x] Obtain the clearance volume the recompression phase starts from: it is
  not a supplied number but a derived one, 0.165 cm³ trapped at
  $2\pi-\delta$ (0.193 cm³ before the vane-root blend), all of it rotor mouth
  cavity gas. The circular-rotor approximation gave 0.00138 cm³ and diverged.
  Cross-checked against the independently recorded 30 mm² mouth free area of
  Section 3.3, measured here as 24.0 mm² (Section 3.4).
- [x] Fix which port $\varphi$ belongs to and how the supplied angles are
  referenced (confirmed 2026-07-23). $\varphi$ is the **suction**-port
  opening angle, measured with $\beta$, $\gamma$, $\delta$ from the vane
  centerline, matching the reference rule set (R2.1 ①) and this
  implementation: suction port $[\varphi, \beta]=[10.4^\circ, 27.7^\circ]$,
  discharge port $[2\pi-\delta-\gamma, 2\pi-\delta]=[339.6^\circ, 346.8^\circ]$.
  The alternative band reading was rejected; in any case it moves
  $\theta_{vo}$ by less than 0.6° (Section 3.4).
- [x] Resolve the earlier note that $\varphi = 10.4^\circ$ does not match the
  vane band ($11.93^\circ$). It is not a discrepancy: $\varphi$ is a machined
  suction-port location and the vane band is the seal-over width, unrelated
  quantities, so there is no reason they should be equal (Section 3.4). The
  only residual check is to confirm the port's angular position on the
  drawing, which does not involve the vane.
- [ ] Confirm the fixed-attitude swing-bush assumption used by the volume
  integration. The rotor tilts $\pm 10.37^\circ$ relative to the bush; the
  groove is axisymmetric about its center so the volume effect is expected
  to be small, but it has not been quantified (Section 3.4).
- [x] Add leakage, or state the recompression peak as an upper bound
  wherever it is used. Added as the parallel Section 3.7 mass-balance model
  (`mochi.leakage`): a single equivalent orifice caps the peak (10.0 → ~9.7 MPa
  at a 5 µm gap) and yields the volumetric efficiency (η_v ≈ 0.92). The
  `port_timed_pressures` baseline still reports the ~10 MPa no-leakage upper
  bound. Remaining: the assumed effective gap (5 µm, editable) and the R410A
  ρ_suc = 31.4 kg/m³ / γ = 1.10 are CoolProp-derived assumptions to confirm
  against measurement; discharge-side reed-valve backflow is still unmodelled.
- [ ] Choose whether time or crank angle is the primary independent variable.
- [ ] Select the contact/friction treatment, if needed.
- [ ] Establish at least one trusted compressor reference case.
- [ ] Set accuracy and conservation acceptance criteria.

Record resolved decisions here with equations and evidence; do not merely check
the box.

## 10. Model-change procedure

Any pull request that changes physics or numerics must update this document in
the same change. State the old and new equation, assumptions, expected effect,
and validation evidence. Update `CHANGES.md`, and add or revise tests that would
fail under the wrong sign, unit, or model.
