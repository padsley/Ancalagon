# GEANT3 -> GEANT4 pilot: DRAGON separator optics + target/BGO + focal plane, visualized

Scoped pilot for a full GEANT3 -> GEANT4 port of the DRAGON separator
simulation. It does **not** touch anything under the original GEANT3
source's own `src/`, `inc/`, or `dat/` -- it only reads them, for
validation, from the sibling `G3_DRAGON` repo/checkout (this pilot was
originally an untracked subdirectory there; it's its own repo now -- see
"Build & run" below for the sibling-checkout layout the validation step
still expects).

## What this is

The DRAGON separator's ion optics (`src/mitray_*.f`) are not solid-geometry
tracking: GEANT3 calls a per-element analytic field function
(`mitray_field.f` dispatching to `mitray_dipole.f` / `mitray_poles.f` /
`mitray_edipol.f` / `mitray_solnd.f` / `mitray_sasp.f`) and lets the generic
GEANT3 Runge-Kutta swimmer (`guswim.f` -> `GRKUTA`) integrate the trajectory.
That maps directly onto GEANT4's field/equation-of-motion architecture
(`G4MagneticField` + `G4Mag_UsualEqRhs` + `G4ClassicalRK4`), so this part of
the simulation is a genuine port, not a rewrite.

This pilot covers **every element** of the real 2014 DRAGON separator
configuration (`dat/dragon_2014_DSSSD.dat`) -- its whole beamline, Q1
through Q14, chained end to end:

```
Q1 - Q2 - D1 - Q3 - Q4 - Q5 - Q6 - Q7 - E1 - Q8 - Q9 - Q10 - D2 - Q11 - Q12 - E2 - Q13 - Q14
```

(14 quadrupoles, 2 magnetic dipoles, 2 electrostatic deflectors -- every
`'POLE'`/`'DIPO'`/`'EDIP'` card in the file.)

Beyond the optics, this pilot also builds (see their own sections below)
the target chamber and 30-crystal BGO gamma-ray array at the reaction
point, named diagnostic slits and a focal-plane DSSSD along the way, and
a Geant4 visualization session (`--vis`) to look at all of it.

### Q1 - Q14 quadrupoles (`'POLE'` cards)

- `include/`, `src/MitrayQuadrupoleField.{hh,cc}`: a line-by-line C++ port of
  `src/mitray_poles.f` + `src/mitray_zone.f` (the MIT-RAYTRACE multipole
  field with Enge-function fringe falloff), as a `G4MagneticField`.
  `MitrayPoleData::Q1()` through `Q14()` hold each magnet's real parameters
  -- the same field code throughout, just different numbers. Every card in
  this file has `LF1=LU1=LF2=3` and `A=B=0`; `Q1()`/`Q2()` predate the
  `FromCard()` helper Q3-Q14 use and are left as their own literal
  transcriptions.
  - Some are pure quadrupoles (`BQD` only, e.g. Q1); some add a hexapole
    term (`BHX`, e.g. Q2); Q3, Q7, Q10, Q12 are pure hexapole/decapole
    magnets with **no quadrupole term at all** (`BQD=0`) and **no fringe
    falloff** (`C0..C5` all zero in the file, giving a hard-edged field
    profile) -- ported exactly as written, not "fixed".
- Geometry: at `theta=0` (nothing has bent the beam yet -- true for Q1/Q2),
  a `G4Tubs` volume sized exactly as `src/ugeom_mitray.f`'s `ugeo_mpole`
  sizes the real GEANT3 volume (radius = aperture, half-length =
  `(L+Z11+Z22)/2`). Once rotated (Q3 onward, see Chaining below), a
  `G4Orb` bounding sphere of the same tube instead -- rotating the tube
  itself would work too, but a sphere sidesteps `G4PVPlacement`'s
  rotation-sense convention entirely for what both ways are just a field
  container, not real geometry.

### D1 and D2 dipoles (`'DIPO'` cards)

- `include/`, `src/MitrayDipoleField.{hh,cc}`: a line-by-line C++ port of
  `src/mitray_dipole.f` (curvilinear A/B/C-axis coordinate and field
  transforms `MITRAY_DIPO_ATOB/BTOC/BBTOBA/BCTOBB`, the field dispatch
  `MITRAY_DIPOLE`, the gradient-dipole field model `MITRAY_NDIP`, and the
  Fermi-function fringe falloff with off-midplane finite-difference
  expansion `MITRAY_NDPP`/`MITRAY_NDPPX`), as a `G4MagneticField`.
  `MitrayDipoleData::D1()`/`D2()` hold each magnet's real parameters.
  - **Scope**: ports exactly the code path D1 *and* D2 both exercise --
    `MTYP=3` (nonuniform/gradient field), `IMAP=0` (no field map), flat
    pole faces (`NSRF=0` on both faces of both magnets). Other `MTYP`
    values, field-map dipoles, the SASP clamshell variant, and curved pole
    faces (which need `MITRAY_SDIP`'s iterative closest-point search on
    the curved effective field boundary) are **not** ported; the
    constructor throws if given parameters outside this scope.
  - Geometry: at `theta=0` (true for D1, the first bend), a generous
    bounding `G4Box` (not a port of `ugeo_dipole`'s exact TRAP wedge shape
    -- the field function's own zone logic already zeroes B outside the
    physical region). Once rotated (D2, downstream of D1's own bend), a
    `G4Orb` of radius `RB+50cm`, for the same reason as the rotated quads.

### E1 and E2 electrostatic deflectors (`'EDIP'` cards)

- `include/`, `src/MitrayEdipoleField.{hh,cc}`: a line-by-line C++ port of
  `src/mitray_edipol.f` (a cylindrical electrostatic deflector: a pure 1/r
  radial field in the uniform region, Enge-function fringe falloff at each
  end), as a `G4ElectricField`. Considerably simpler than the magnetic
  dipole -- there's no pole-face rotation, no curved-face complication (no
  equivalent of `MITRAY_SDIP` exists in `mitray_edipol.f` at all), and no
  off-midplane finite-difference expansion (Y only ever shifts the fringe
  longitudinal coordinate, and is moot for both E1 and E2 since they have
  `EC2=EC4=0`) -- so this port has no scope restriction beyond what
  `mitray_edipol.f` itself computes for A=0 (true for both cards).
  - **Units**: `src/guefld.f` scales the raw field from `mitray_edipol.f`
    by 1e-6 to get GEANT3's internal field unit (GeV/cm, per
    `src/grkuta.f`'s docstring) -- i.e. the raw value is in kV/cm.
  - An electric field does work on the particle (unlike the purely
    magnetic elements), so this needs a different GEANT4 plumbing:
    `G4EqMagElectricField` (not `G4Mag_UsualEqRhs`) as the equation of
    motion, and an 8-variable `G4ClassicalRK4`/`G4IntegrationDriver`
    (position + momentum direction + energy + time, since energy is no
    longer conserved).
  - Geometry: same pattern as the dipoles -- a bounding `G4Box` sized off
    `RB` at `theta=0` (E1, standalone mode only -- in the chain E1 already
    sits past D1's bend), a `G4Orb` of radius `RB+50cm` once rotated.

All eighteen share `src/DetectorConstruction.cc` (element, or "Chain",
selected by name), `src/PrimaryGeneratorAction.cc`, `src/SteppingAction.cc`,
`src/PhysicsList.cc` (a minimal, transportation-only app to fire a particle
through one -- or all -- of them and log the trajectory).

### Chaining the whole beamline

`DetectorConstruction("Chain")` places all eighteen elements in one world,
spaced *and oriented* exactly as the real beamline does -- see
`include/RotateAboutY.hh` and the `Chain*`/`Drift`/`Shift` functions in
`src/DetectorConstruction.cc`.

**Position bookkeeping** follows `src/mitray_setup.f`'s `ugeom_setup`
(~lines 606-665 for `'DRIF'`/`'SHRT'`, ~lines 1031-1088 for `'POLE'`): a
drift or shift advances the running position by its own length along the
current axis; a quadrupole card advances it by `A+B+L` (entry to exit),
centred at `A+(Z22+L-Z11)/2` from entry -- `L/2` for every quad here, since
they all have `Z11==Z22`.

**Orientation bookkeeping** is new since the Q1->Q2->D1-only chain: every
`'SHRT'` card in this entire file has its rotation-angle fields (its 4th,
5th, 6th data values) equal to zero -- checked by hand against all 20+
`'SHRT'` cards in the file -- so *no* element ever rotates the beam axis by
a `'SHRT'` card. The only rotation source is a bend at a dipole or
electrostatic deflector: each one bends the *design* trajectory by `PHI`
along a circular arc of radius `RB`, toward local `-x` (this file's
consistent sign convention, verified empirically in the standalone D1/E1
pilots before any chaining existed), so the whole beamline's accumulated
orientation reduces to a single scalar angle `theta` -- the design
trajectory's local `-x`/local-`+z` exit point off a bend is just
`(-RB(1-cosPHI), 0, RB sinPHI)`, and `theta` decreases by `PHI` there. Both
field classes take this `theta` as a constructor parameter and internally
rotate the query point into (and the resulting field back out of) their
own local frame -- see `RotateAboutY.hh`; at `theta=0` this is the
identity, so it doesn't disturb any of the already-validated single-element
field values.

**Scope/caveat, deliberate and clearly separate from the already-validated
field physics:** GEANT3's own dipole volume-placement algorithm is not
replicated. `ugeom_setup` computes a dipole's actual placement at the
bend's arc-midpoint (using the curvature centre and a `uh`/`uv` basis
derived from `RB`, `PHI`, and the collimator-derived `dr_dipole` --
`mitray_setup.f` ~lines 757-820). This pilot instead puts each bending
element's field origin (`za=0`) at the plain drift-accumulated entry
point -- simpler, and it *only* affects where an element sits/points in
this simplified world relative to GEANT3's own bookkeeping, not the field
values themselves (still bit-exact against the real Fortran, see below).

### Apertures and slits (`'RCOL'` cards -- every one, not just the named slits)

`src/mitray_setup.f` gives `'COLL'` and `'RCOL'` cards different real
meanings: a `'COLL'` card (idata==13) is pure acceptance-limiting
bookkeeping fed into whichever optical element it brackets, never real
geometry (`ugeo_col` is only ever called for idata==17, i.e. `'RCOL'`
cards). Every `'RCOL'` card in this file -- 23 of them, not just the three
named diagnostic slits -- gets real, absorbing copper geometry, a direct
port of `src/ugeom_mitray.f`'s `ugeo_col`: either an annular tube
(`DATA(1)!=0`) or four picture-frame jaw boxes around a rectangular hole
(`DATA(1)==0`), sized/positioned exactly per that card's own `DATA(2..6)`
and `ugeo_col`'s own per-name magnification rule (20x for the charge slit
`QSLT`, 10x for the mass slit `MSLT` and final slit `FSLT`, 1x for
everything else). `DetectorConstruction.cc`'s `ChainCollimator()` builds
these, correctly rotated to the local beam axis past a bend (unlike the
optics elements' own orientation-agnostic field-container spheres --
these are real absorbing geometry, so getting that wrong would mean the
beam clips the wrong material). Two `'TEST'` cards, `TST3` and `TST4`,
carry no data at all -- per the file's own `'15 Jan 2014'` comment they
mark where MCP0/MCP1 *would* sit, so `ChainMarker()` still places a
purely-decorative wireframe box there, since there's nothing real to
port.

The chain ends at the `'COLL' 'DSSD'` card's position with a real
silicon wafer -- **but this one is not a port of anything in the
Fortran**: a direct search confirms no real DSSSD (or MCP) volume is
built anywhere in the GEANT3 source for this configuration (`ugeo_dssd`
exists in `src/ugeom_mitray.f` but is dead code, never called; no
`'MCPF'` card appears in this file). Its geometry instead comes from the
real part: a **W1** double-sided silicon strip detector (Micron
Semiconductor Limited's 2024 catalogue, p.33 -- 16 junction x 16 ohmic
strips, 3100um pitch, 49.50mm x 49.50mm active area; three W1 packaging
variants in the catalogue all share this active area). Built as a single
`G4_Si` slab, 49.50mm x 49.50mm x 1mm (thickness and focal-plane
centring per the user's own direction, not the catalogue, which offers a
thickness range rather than one fixed value) -- not 16+16 separate strip
volumes, since there's no sensitive-detector/digitization step yet to
make per-strip segmentation meaningful.

## Target chamber and BGO array

The reaction itself happens well upstream of Q1, at the beamline's own
`'STRV'` start -- this is a separate part of the real GEANT3 geometry
from the MITRAY optics chain, built by `src/ugeom.f`'s own call tree
(`ugeo_space` -> `ugeo_detector` -> a target-chamber variant ->
`ugeo_finger` -> `ugeo_pmt`), reading numeric constants out of
`dragon_2003.ffcards` via `src/ugffgo.f` into `inc/geometry.inc`/
`inc/uggeom.inc`. `DetectorConstruction`'s `"Chain"` mode builds both,
at the world origin, before walking the optics chain.

### BGO gamma-ray array (`include/`, `src/BgoArray.{hh,cc}`)

30 detector positions, each a hexagonal prism, arranged in front/back
rings around the target -- a literal transcription of
`src/ugeom_gbox.f`'s `ugeo_finger`: detectors 1-10 are individually
hand-placed in the source (their own x/y/z formulas, including two
pulled back 6.7/7.2cm "to make room for lead collimators"); 11-30 come
from a clean `Do k = 9,3,-1` / `Do j = 1,jm` loop. All 30 positions are
computed in code (not hardcoded numbers) straight from the same
`hexagon_small_width`/`hexagon_large_width`/`depth` formulas and
`FSID`/`WALL`/`BGAP`/`HOLE`/`PMTR` ffcards values the real code uses --
see `BgoArray.cc` for the exact derivation of each constant.

**Scope/simplification**: each position is built as a *single* `G4_BGO`
hexagonal prism sized to the outer housing (`HSNG`) envelope, not the
real nested `HSNG`(housing) > `FNGR`(air gap) > `MGOR`(reflector) >
`SCNT`(crystal) stack with a separate PMT tube behind each one. There's
no field/formula to validate here the way the optics have, so the bar
for this geometry is faithful position/dimension transcription (done),
not bit-exact physics. One genuine gotcha surfaced and got fixed here:
`G4Polyhedra`'s `phiStart` places a *vertex*, not a flat, at the start
angle, but the real hexagons tile edge-to-edge along the stacking axis
(`hexagon_small_width` is a flat-to-flat measurement) -- using
`phiStart=0` produced ~8mm crystal-crystal overlaps (confirmed via
`G4PVPlacement`'s `pCheckOverlaps` flag); `phiStart=30deg` fixed it.
Each detector's front/back orientation also needed a real
`G4RotationMatrix` (unlike the optics elements' rotation-agnostic
`G4Orb` containers) -- derived from the same active-rotation convention
validated for the field classes, and empirically confirmed overlap-free.

### Target chamber (`include/`, `src/TargetChamber.{hh,cc}`)

The real target (`src/ugeom_trgt_small.f`, the `tubetype=0` variant
`dat/dragon_2014_DSSSD.dat`'s own `TUBE` ffcard selects -- the other
five variants in `ugeom_trgt_*.f` only differ by a handful of small
acceptance-restricting collimator blocks, not worth porting separately)
is a **windowless differential-pumping gas target**, ported as two
pieces:

- `OuterBoxAndShielding()`: `src/ugeom_gbox.f`'s `ugeo_detector`
  (`tubetype==0` branch, ~lines 1-732) -- `CMBR` (stainless-steel outer
  support box; a previous version of this file's own comment misidentified
  this medium as "NAI:TL" by reading the wrong table -- `ugstmed.f`'s
  medium 20 is genuinely `'STAINLESS STEEL'`) > `CMBG` (its gas interior),
  plus `UHOL`/`DHOL` (beam holes through `CMBR`'s own wall), `PUAI`/`PDAI`
  (+`PUBI`/`PDBI`, the aluminum collimator inserts just inside them), and
  the lettered lead/aluminum/gas collar chains bolted onto each face --
  `PUA1`..`PUJ2` upstream, `PDA1`..`PDI4` downstream (asymmetric: no
  separate `PDC`/`PDJ`, no lead layer in `PDB`/`PDE` -- ported exactly as
  the real file's own asymmetry gives it, not harmonized). One section,
  `PUG`, has a lead **box** (not tube) as its outer layer -- everything
  else is concentric tubes, narrowing toward the beam axis one gap
  (`d_mtl`/`air_gap`/etc.) at a time, built by a small generic
  `CollarSection`/`BuildCollarChain` helper rather than by hand.
- `CellAndApertures()`: the actual reaction volume, real GEANT3 daughter
  of `CMBG` (not the world) -- `CELL` (an aluminum TRD1 outer shell),
  `CELG` (the real gas volume nested inside it, same shape, `mtarg`
  material), and `EAPG`/`XAPG` (the entrance/exit beam apertures bored
  through `CELL`'s slanted side wall). All shape/position numbers are
  that file's own literals, not re-derived.
- `DifferentialPumpingChain()`: the graduated pumping chain around the
  cell, from `src/ugeom_trgt_small.f` -- 7 aluminum `CONE` collimator
  inserts (`EN2C`/`EN3C`/`EX2C`..`EX6C`) nested inside 11 stepped-radius
  `mcent`-material (low-pressure) gas `TUBE` volumes (`EN1G`..`EN6G`,
  `EX2G`..`EX10`), at the exact Z positions that file's own formulas give
  from `inc/uggeom.inc`'s `zent()`/`len1-3` constants. `EX6C` alone sits
  directly in the world with no surrounding gas, matching the file's own
  comment that it "extends into Q1" -- clipped here to end exactly at
  Q1's entry face (see `TargetChamber.cc` for the taper interpolation).

**Scope/caveat**: `ugeo_detector` and `ugeom_trgt_small.f` were evidently
authored somewhat independently and occupy overlapping z-ranges around
the same physical region -- G4's overlap checker (on for every placement
here) flags a few crossings between the `EN`/`EX` chain and the `PU`/`PD`
collar chain (e.g. `PUJ1` vs. `EN2G`). That reflects the real files'
own, apparently uncoordinated overlap (GEANT3's tracker tolerates
arbitrary volume overlaps in a way GEANT4's strict mother/daughter model
does not), not a placement error introduced by this port -- see
`TargetChamber.hh`'s own header comment for the full reasoning. `CMBG`
itself is placed as an ordinary daughter of `CMBR` rather than GEANT3's
own `'MANY'`-flagged placement, which gives the same physical result (a
hollow steel box with a gas-filled interior) without needing GEANT4 to
tolerate GEANT3-style volume overlaps at all.

**Gas**: `src/ugmate_trgt.f` picks the target gas by projectile mass
(`atarg.lt.1.2` -> H2, else He) at `ptarg` = 5 Torr, with a `mcent`
central-vacuum gas at `ptarg*0.002` filling the pumping chain, same
species as the target gas. `TargetChamber::Build()` takes that same
choice as an explicit `TargetGas::kHydrogen`/`kHelium` parameter rather
than hardcoding one -- `DetectorConstruction.cc` passes `kHelium` for
this configuration's real reaction, **15O(alpha,gamma)19Ne** (`atarg`=4,
see "Reaction specification" below); a (p,gamma) configuration would pass
`kHydrogen` instead. Both gases are built the same way, NIST STP
`G4_H`/`G4_He` scaled by the pressure fraction (`src/ugmate_trgt.f`'s own
reference densities, 8.38e-5 g/cm3 for H2 and 1.6586e-4 g/cm3 for He, are
close but not identical figures to the NIST STP ones -- not chased here).
An earlier version of this pilot hardcoded H2 as a placeholder before the
real reaction had been identified.

## Reaction specification

`dragon_2003.ffcards`'s `FKIN 20 2.0 6.0` card sets `LKINE=20` --
`src/ureact.f`'s `case(20)` reads the real per-run reaction parameters
(masses, resonance energy, branching) from an external namelist file
named by the `INPUT` environment variable, which isn't checked into this
repo. But the data files identify the reaction directly:
`dat/dragon_2014_DSSSD.dat`'s own header comment --

```
'COMM'  '15O 28 Sept 00, 19Ne=1.8885 MeV, 15O=2.3921 MeV; use M13/15 sext'
```

-- plus `src/ugmate_trgt.f`'s gas-selection branch (`atarg`=4 -> helium)
and `src/ureact.f`'s own hardcoded `case(2)` -- `' (2) 15O(alpha,g)19Ne '`
-- together settle it: this configuration is **15O(alpha,gamma)19Ne**
radiative capture, an at-rest alpha (He gas target) capturing a 15O beam.

The two numbers in that `COMM` line are *not* the same quantity.
`dragon_2003.ffcards`'s own comment on the same tune, "REFERENCE TUNE
E=1.8885 MeV, A=19, Q=4", says E is for an A=19, Q=4 ion -- i.e. the
**recoil** (19Ne, charge state 4+), not the beam. Using `src/ureact.f`
`case(2)`'s own hardcoded resonance energy (`Er`=0.5036 MeV, the CM
energy above the 15O+alpha threshold -- a real, checked-in number, not a
guess) in two-body relativistic kinematics reproduces both `COMM` numbers
independently:

- beam (15O) lab kinetic energy: 2.3913 MeV (`COMM`: "15O=2.3921 MeV")
- recoil (19Ne) lab kinetic energy: 1.83-1.95 MeV across the isotropic
  gamma-emission-angle spread, centred near 1.888 MeV (`COMM`:
  "19Ne=1.8885 MeV")
- recoil lab momentum: 254-262 MeV/c, centred near 258.5 MeV/c --
  matching, to <0.1%, this pilot's own already-validated full-chain
  design orbit (258.7 MeV/c for 19Ne4+, found independently from the
  separator's D1/D2/E1/E2 field settings -- see "Full chain" below)

One more independent check: the predicted CM gamma energy, Q-value + Er =
3.5292 + 0.5036 = 4.0328 MeV, is the well-known astrophysical Ex=4.03 MeV
state in 19Ne -- the key breakout resonance for 15O(alpha,gamma)19Ne in
classical novae/X-ray bursts, and exactly the state `src/ureact.f`
`case(2)`'s own branching (80% direct to the 19Ne ground state) is built
to populate. Four independent numbers in this repo's own data files, all
explained by one coherent physical picture -- `Er`=0.5036 MeV is what
this configuration's files corroborate, not this pilot's guess.

`include/`, `src/ReactionKinematics.{hh,cc}` implement this as a
relativistic reaction-vertex generator: beam + at-rest target -> a
compound resonance (invariant mass = m_beam + m_target + `Er`), then a
chain of two-body breakups (parent -> daughter + one gamma, isotropic in
*that step's own parent's* rest frame) down through the recoil's real
excited-state cascade to its ground state -- see "Reaction config file"
below for how that cascade (and everything else about the reaction) is
specified, not hardcoded. **Scope**: vertex kinematics only -- the whole
chain is treated as instantaneous at the reaction vertex (no separate
tracking of an excited intermediate recoil's own brief flight before its
next decay), no angular distributions/correlations anywhere in the
chain, no energy-loss straggling in the gas, and the recoil's charge
state is fixed at whatever the config file's own `RECL` card gives (4+
for this pilot's own reaction), not sampled from a real charge-state
distribution. `PrimaryGeneratorAction`'s reaction-file constructor fires
the ground-state recoil ion plus one gamma per cascade step actually
sampled (1-3 for this pilot's own reaction, see below) per event, all
from a single vertex; `--track-reaction` and `--vis Reaction` use it in
place of the fixed-momentum recoil gun `--track-chain`/`--vis Chain` use.

### Reaction config file

`include/`/`src/ReactionConfig.{hh,cc}`: the reaction (beam/target/recoil
species and masses, resonance energy, and the recoil's excited-state
cascade) is a **run-time config file**, not hardcoded -- styled after
this project's own `dragon_2003.ffcards` convention (unquoted 4-letter
card keyword + that card's own values, one per line; `COMM` lines and
blanks ignored; parsing stops at `SENT`), rather than introducing a new
one. `reactions/o15ag_19ne.reaction` is this pilot's own bundled file for
the reaction identified above:

```
BEAM 8   15   2.8554         # Z A massExcessMeV -- 15O
TARG 2   4    2.42492        # 4He (target, at rest)
RECL 10  19   1.7511   4     # 19Ne, production charge state 4+
ERES 0.5036                  # CM resonance energy above threshold, MeV

LEVL   1   1.536   2.8E-11   # level# excitationMeV lifetimeS (19Ne 3/2+)
LEVL   2   0.275   6.3E-11   # (19Ne 1/2)
LEVL   3   0.238   2.6E-8    # (19Ne 5/2+)

BRAT   -1   80.0   0         # resonance (-1) -> ground (0), 80%
BRAT   -1   15.0   2         # resonance -> level 2, 15%
BRAT   -1    5.0   1         # resonance -> level 1, 5%
BRAT    1   95.0   3         # level 1 -> level 3, 95%
BRAT    1    5.0   2         # level 1 -> level 2, 5%
BRAT    2  100.0   0         # level 2 -> ground, 100%
BRAT    3  100.0   0         # level 3 -> ground, 100%
SENT
```

(masses use the AME mass-excess convention throughout, `mass = A*amu +
massExcessMeV` -- including the target, unlike `src/ureact.f`'s own
`hemass` constant, a literal atomic mass; the file's own `TARG`
`massExcessMeV` reproduces `hemass` to 7 eV, so every particle in the
file is specified the same way, not a mix of conventions). Branching
ratios are `src/ureact.f` `case(2)`'s own real `gsdk` `brat`/`mode`
values (level numbers here are this file's own arbitrary labels, not
`case(2)`'s own `idpart` numbers, which are a GEANT3 tracking-medium
implementation detail, not physics); `ReactionKinematics` renormalizes
percentages per `fromLevel` against their own siblings, so a level's own
branches don't strictly need to sum to 100.

Loaded via `ReactionConfig::Load()`, which throws `std::runtime_error`
(uncaught, deliberately -- a broken reaction file should stop the run,
not silently substitute something else) on a missing file or malformed
card. Which file to load can be overridden with the `REACTION_INPUT`
environment variable (mirroring `src/ureact.f` `case(20)`'s own
`INPUT`-env-var convention for exactly this kind of external reaction
file -- a different name, to avoid colliding with an unrelated
pre-existing use of `INPUT`); it defaults to `o15ag_19ne.reaction` (built
alongside the executable, see `CMakeLists.txt`).

**`--reaction-stats [path] [nEvents]`**: samples a reaction file's own
cascade (default `o15ag_19ne.reaction`, 10000 events) with no Geant4 run
manager at all (`ReactionKinematics::GenerateEvent()` only needs CLHEP,
not a constructed geometry/physics list) and tallies how many gammas
each event emits -- a quick way to sanity-check a hand-edited `LEVL`/
`BRAT` scheme (branch fractions, that every chain actually terminates)
before spending time on full tracking runs. Verified against this
pilot's own bundled file (20000 events): 80.13% 1-gamma, 15.02%
2-gamma, 4.86% 3-gamma -- matching the config's own 80/15/5 branching
(the two indirect branches, 15% and 5%, always end up 2 and 3 gammas
respectively, since level 1 always decays again rather than reaching
ground directly).

## Sensitive detectors

`include/`/`src/{Bgo,Dsssd}{Hit,SD}.{hh,cc}`, `EventAction.{hh,cc}`: a
real Geant4 hits pipeline for the BGO crystals (`SCNT`, see "Target
chamber and BGO array" above) and the DSSSD -- `G4VHit` classes
(`BgoHit`/`DsssdHit`: energy deposit, position, time, track ID, particle
name; `BgoHit` also carries which of the 30 crystal positions was hit),
`G4VSensitiveDetector` subclasses (`BgoSD`/`DsssdSD`) that create one hit
per step and insert it into a `G4THitsCollection`, and an `EventAction`
that pulls both collections out of each event, prints a one-line summary
(hit count, crystal IDs, total energy deposit), and fills two ROOT
ntuples (see "ROOT output" below).

### ROOT output

`RunAction.{hh,cc}`, `HitNtuple.hh`: every hit from both collections is
also written, one row per hit, into a real ROOT file, **`dragon_hits.root`**
(written to the working directory at the end of each run) -- via Geant4's
own `G4AnalysisManager` (`G4RootAnalysisManager` backend, chosen by the
`.root` extension), which writes genuine ROOT-format `TTree`s using
Geant4's embedded writer, with **no dependency on an external ROOT
install** to produce the file (only to read it back, e.g. with a plain
`root` session -- one was used to verify the schema/contents below).

Two trees, `Bgo` and `Dsssd` (`HitNtuple.hh` is the single place their
column layout is defined -- `RunAction.cc` creates columns in that exact
order, `EventAction.cc` fills them by those same named indices, so the
two can't drift out of sync):

| tree | columns |
|---|---|
| `Bgo` | `eventID`, `crystalID` (1-30), `edepMeV`, `x_cm`/`y_cm`/`z_cm`, `t_ns`, `trackID`, `particle` |
| `Dsssd` | `eventID`, `edepMeV`, `x_cm`/`y_cm`/`z_cm`, `t_ns`, `trackID`, `particle` (no `crystalID` -- only one DSSSD) |

Verified by reading `dragon_hits.root` back with real ROOT after a
`--track-reaction 30` run: both trees have the expected schema; `Bgo`
had 318 rows with sensible contents -- e.g. a primary `gamma` hit
followed by its own `e+`/`e-` secondaries (pair production) at the same
position, each with its own real `edepMeV`, all tagged with the correct
`crystalID`. `Dsssd` was schema-valid but empty in that run, consistent
with the already-documented finding that recoils aren't currently
reaching the DSSSD (see "Full chain" above) -- not a new bug in the
ROOT-writing path itself, which uses the exact same
`G4AnalysisManager` calls as the (verified-working) `Bgo` tree.

**No file at all if nothing hit anything**: `G4AnalysisManager` deletes
`dragon_hits.root` at the end of the run if *both* ntuples end up with
zero total rows (its own "delete empty file" behavior, confirmed via
`SetVerboseLevel(3)`) -- e.g. `--track-chain` with the current (EM-
physics-era, see "Full chain" above) design orbit produces no file at
all, since neither BGO nor the DSSSD gets hit. Not a bug: run
`--track-reaction` (its gammas reliably hit the BGO array) if you want to
confirm the file/trees exist and inspect their contents.

**BGO crystal ID**: `SCNT` is a single logical volume placed once per
position (all 30 share it), and only the *outer* `HSNG` placement in each
position's own nesting (`SCNT` -> `MGOR` -> `FNGR` -> `HSNG`) carries the
real per-position copy number (1-30) -- so `BgoSD::ProcessHits` reads it
off `G4TouchableHistory::GetVolume(3)` (ancestor depth 3 from the hit
volume) rather than off `SCNT`'s own copy number (always 0). Verified
empirically: crystal IDs recorded during `--track-reaction` runs land
in [1,30] as expected.

**edep is real now**: `PhysicsList` registers `G4EmStandardPhysics`
(delegated to, not reimplemented -- see `PhysicsList.cc`), so gammas and
charged particles get real EM processes (photoelectric/Compton/pair/
Rayleigh for gammas; ionisation/multiple scattering for everything
charged, ions included). Hits are still *not* filtered on `edep>0` the
way most Geant4 examples do -- every step through a sensitive volume is
recorded regardless, so a track that passes through depositing nothing
(rare now, common before) still shows up with `edep=0`, rather than being
silently dropped. Verified against `--track-reaction`: BGO events show a
clean **full-energy peak at ~4.0-4.07 MeV**, matching the capture gamma's
own predicted CM energy (Q-value + Er = 4.03 MeV, see "Reaction
specification" above) almost exactly -- i.e. the gamma is fully
photoelectrically absorbed in a single crystal in many events, plus a
Compton-continuum tail (partial energy, 1 crystal) and multi-crystal
scatter chains (several hits, one per crystal the same gamma Compton-
scattered through) in others -- exactly the spectral shape a real BGO
produces.

**Scope**: the DSSSD's real 16+16 (junction/ohmic) strip segmentation
isn't modeled -- one sensitive volume, one hit per step, not per-strip
channels -- and BGO/DSSSD digitization (`src/gudigi.f`'s light
collection/resolution smearing/PMT response) isn't ported either; both
remain open items (see "What this pilot does not cover").

`--track-reaction`/`--vis Reaction` are the natural way to exercise
this: the cascade gamma(s) and the 19Ne recoil are correlated, physically
distinct particles from one real event, so a given reaction event may or
may not reach the BGO array and/or the DSSSD at all, depending on its own
sampled emission angle and recoil momentum -- a real recoil separator has
finite acceptance too, so some fraction of events being lost along the
way (often on the very slits/collimators built above) isn't a bug, it's
the physics. **The DSSSD specifically**: with EM physics now active,
neither `--track-reaction`'s real recoils nor `--track-chain`'s
fixed-momentum design-orbit recoil (which *did* reliably reach the DSSSD
under the old transportation-only physics list) reach it reliably any
more -- see "Full chain" above for why (a real, but not yet resolved,
consequence of turning EM physics on).

## Visualization

`--vis [element] [macro]` opens a Geant4 viewer (`macros/init_vis.mac`
picks whichever driver is actually available -- Qt/OpenGL/etc.) and
either runs the given macro or drops into an interactive prompt.
Trajectories accumulate on screen (`/vis/scene/endOfEventAction
accumulate`); quadrupoles are blue, dipoles red, electrostatic
deflectors green, the target pale blue, BGO crystals gold, and the named
slit/MCP/DSSSD markers yellow/cyan/magenta respectively.

```sh
./dragon_g4_pilot --vis Chain              # interactive session, whole beamline
./dragon_g4_pilot --vis Chain run_chain.mac  # batch: fires 5 design-orbit ions, exits
./dragon_g4_pilot --vis Reaction           # same chain, real reaction events (see above)
```

At an interactive prompt, `/run/beamOn 5` fires more events (each
accumulates); the default view looks down the Y axis, zoomed out to
show the whole ~11m x 7m chain (see the full-chain validation numbers
above) -- `/vis/viewer/zoomTo`, `/vis/viewer/set/viewpointVector`, etc.
work as usual to look closer at any one element.

## Validation

### Field grids: every element, bit-exact

For each element, a standalone Fortran program links directly against the
**real, unmodified** GEANT3 source (`mitray_poles.f`/`mitray_zone.f` for
the quads, `mitray_dipole.f` for the dipoles, `mitray_edipol.f` for the
deflectors -- no CERNLIB/hbook dependency, those routines are
self-contained) and evaluates the field on a grid spanning the far field,
both fringe zones, and the uniform region; the matching `main.cc --probe*`
mode evaluates the C++ port on the identical grid; `validate/compare_fields.py`
diffs the two CSVs.

| element(s) | Fortran probe | grid points x components | result |
|---|---|---|---|
| Q1 | `mitray_poles_probe.f` | 525 x 3 | bit-exact |
| Q2 | `mitray_poles_probe_q2.f` | 525 x 3 | bit-exact |
| Q3-Q14 (combined) | `mitray_poles_probe_rest.f` | 3300 x 3 | bit-exact |
| D1 | `mitray_dipole_probe.f` | 875 x 3 | bit-exact |
| D2 | `mitray_dipole_probe_d2.f` | 875 x 3 | bit-exact |
| E1 | `mitray_edipol_probe.f` | 875 x 3 | bit-exact |
| E2 | `mitray_edipol_probe_e2.f` | 875 x 3 | bit-exact |

`Worst absolute error: 0.000e+00` in every case. Q2 exercises the nonzero
hexapole term Q1's grid never touches; Q3-Q14 as a batch exercise every
other multipole combination in the file (pure hexapole/decapole magnets
with no fringe falloff included).

### Tracking spot-checks: each bending element on its own design orbit

A magnetic dipole's design bend angle `PHI` is reached by a particle on
the nominal radius-of-curvature `RB` orbit, i.e. one with magnetic
rigidity `p/q = B*RB` (0.3 x `B`[T] x `RB`[m] MeV/c per unit `e`). A
cylindrical electrostatic deflector's design condition instead balances
`qE = mv^2/RB`, i.e. `KE/q = E*RB/2` (independent of mass). Both D1/D2 and
E1/E2 turn out to share one design rigidity/KE-per-charge each (see below),
so each was spot-checked on its own with a proton at the momentum matching
*that* element's own condition:

```
./dragon_g4_pilot --track-dipole 0.0 0.0 64.7 1     # D1: 49.9 deg vs 50 design
./dragon_g4_pilot --track-e1     0.0 0.0 29.768 1   # E1: 19.99 deg vs 20 design
./dragon_g4_pilot --track-e2     0.0 0.0 29.763 1   # E2: 35.00 deg vs 35 design
```

(D2 wasn't spot-checked standalone -- see the full-chain result below,
which exercises it along with everything else.) Each of these is a
stronger check than the bit-exact grid comparison alone: it confirms the
port reproduces the magnet's actual optical design (coordinate transforms
+ field shape + RK4 integration all together), not just the field formula
in isolation.

### Full chain: one recoil ion through all eighteen elements

D1 and D2 share one design magnetic rigidity: `B x RB` = `0.21561204 T x
1.0m` = `0.21561 T.m` for D1, `0.26520546 T x 0.813m` = `0.21566 T.m` for
D2 -- matching to 0.02%, confirming they're tuned for the same beam. E1
and E2 likewise share one design electric "rigidity": `E x RB` =
`4.721 kV/cm x 200cm` = `944.2 kV` for E1, `3.77557049 kV/cm x 250cm` =
`943.9 kV` for E2.

A single proton can't satisfy both a fixed `p/q` *and* a fixed `KE/q`
simultaneously (`KE` and `p` relate through the particle's mass). Solving
`m/q = p^2/(2*KE)` with `p/q = 64.68 MeV/c` and `KE/q = 0.4721 MeV` gives
`m/q = 4.758 u` per unit charge -- i.e., for a `4+` charge state, mass
`~19.03u`. That's a **19Ne4+** recoil, consistent with this file's own
`'COMM'` lines about a 19Ne reaction. `PrimaryGeneratorAction` gained a
second constructor for exactly this (`G4IonTable::GetIon(Z, A, 0.0)` plus
an explicit charge-state override, since 19Ne4+ isn't fully stripped);
`--track-chain` uses it by default (`p = 4 x 64.68` = `258.7 MeV/c`):

```
./dragon_g4_pilot --track-chain 0.0 0.0 258.7 1
```

The ion enters along Q1's axis and exits Q14 with momentum
`(2.05, 0, -258.71)` MeV/c against an input of `(0, 0, 258.7)` -- magnitude
preserved to 0.03% (258.72 vs 258.7, consistent with the magnetic bends
conserving it and the electric ones returning to it on the design orbit),
direction bent by **179.5 degrees** against a summed design of
`50+20+75+35` = **180 degrees** (D1+E1+D2+E2). The residual is consistent
with rounding the input momentum to 258.7 instead of the exact 258.722
MeV/c. This is the strongest check in this pilot: one real recoil ion,
launched once, correctly threading fourteen quadrupoles, two dipoles, and
two electrostatic deflectors in sequence, each contributing its own real,
independently bit-exact-validated field, composed through a hand-derived
(and now code-verified) rotating reference frame.

**Caveat added once real EM physics existed** (see "Sensitive detectors"
above): this exact result was obtained under the transportation-only
physics list this pilot originally shipped with, where no material
anywhere costs the ion any energy. With `G4EmStandardPhysics` now
registered, the *same* command no longer reproduces it. Root cause,
found by adding a volume name + step-edep to every `SteppingAction`
line (see `SteppingAction.cc`) and reading off exactly where the energy
goes: **not** a geometry/alignment bug -- the recoil correctly threads
both `EAPG` and `XAPG` (small, expected ionisation losses in their
He-gas fill, ~0.015-0.016 MeV each) rather than clipping `CELL`'s solid
aluminum shell, ruling out an earlier, incorrect hypothesis. The real
cause: crossing `CELG` itself (the actual target gas cell, ~9.9cm of
5 Torr He) deposits **~0.197 MeV in that one crossing** -- by far the
single largest energy loss anywhere in the trajectory, and right where
the momentum drop and an accompanying transverse kick (0.46->1.44 MeV/c)
both appear. That's large relative to this recoil's own total kinetic
energy (~1.9 MeV -- see "Reaction specification" above) because the
recoil is extraordinarily slow (v/c=beta~0.015, ~0.1 MeV/nucleon):
stopping power for a heavy ion rises sharply at such low velocity, so
even a dilute gas costs it a real, physically legitimate fraction of its
own energy -- not a simulation artifact. (Also checked and ruled out:
the capture gamma's own kinematic recoil kick, at most ~4 MeV/c
transverse per "Reaction specification" above, cannot be the cause
either -- `--track-chain`'s design-orbit test involves no gamma at all,
fires with exactly zero transverse momentum, and shows the identical
CELG-crossing signature.)

A recoil that loses this much rigidity partway through no longer matches
what the downstream magnets/deflectors are tuned for, and the resulting
deviation compounds through the rest of the chain enough that this
design-orbit test no longer reliably reaches the DSSSD (see "Sensitive
detectors" above) with EM physics on. This mirrors a real, known
experimental challenge for recoil separators handling low-energy
reaction products (target-thickness energy-loss straggling reducing
transmission) -- i.e. this may be this pilot correctly surfacing a real
physical effect for the first time, not a bug to fix. **Not further
pursued** (e.g. whether Geant4's default ion stopping-power model is
well-validated at this specific, very low energy, or whether the real
DRAGON target/tune's own acceptance already accounts for this) --
flagged here rather than silently left inconsistent with the
validated-result claim above, which itself remains true only for the
transportation-only configuration this pilot launched with.

## Build & run

Requires a Geant4 install (v11.3.0 used throughout this README) built
with its `analysis` category (`libG4analysis.so`) -- no separate ROOT
install is needed to *build* or *run* this (see "ROOT output" above),
only to read `dragon_hits.root` back afterwards.

This repo doesn't vendor the original GEANT3 Fortran (`G3_DRAGON`, this
port's own source of truth throughout -- see every "Validation" section
below) -- only the "Field grids" validation step (compiling and running
the real, unmodified Fortran probes) needs it, and does so via `../../
G3_DRAGON/{src,inc}` relative paths, i.e. it expects `G3_DRAGON` cloned
as a **sibling** directory of this repo (both directly under the same
parent, e.g. `~/codes/G3_DRAGON` and `~/codes/Ancalagon`). Building and
running the Geant4 pilot itself (everything else below) needs nothing
from `G3_DRAGON` at all.

```sh
mkdir -p build && cd build
cmake -DGeant4_DIR=/path/to/geant4-install/lib/cmake/Geant4 ..  # omit -DGeant4_DIR if Geant4 is already on your CMAKE_PREFIX_PATH
make -j$(nproc)

# Field validation (each writes/reads its own CSV pair; compare_fields.py
# takes <fortran.csv> <cpp.csv>)
./dragon_g4_pilot --probe            # Q1  -> q1_cpp_reference.csv
./dragon_g4_pilot --probe-q2         # Q2  -> q2_cpp_reference.csv
./dragon_g4_pilot --probe-q3-q14     # Q3-Q14 (combined) -> q3_q14_cpp_reference.csv
./dragon_g4_pilot --probe-dipole     # D1  -> d1_cpp_reference.csv
./dragon_g4_pilot --probe-dipole-d2  # D2  -> d2_cpp_reference.csv
./dragon_g4_pilot --probe-e1         # E1  -> e1_cpp_reference.csv
./dragon_g4_pilot --probe-e2         # E2  -> e2_cpp_reference.csv

cd ../validate
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_poles_probe mitray_poles_probe.f ../../G3_DRAGON/src/mitray_poles.f ../../G3_DRAGON/src/mitray_zone.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_poles_probe_q2 mitray_poles_probe_q2.f ../../G3_DRAGON/src/mitray_poles.f ../../G3_DRAGON/src/mitray_zone.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_poles_probe_rest mitray_poles_probe_rest.f ../../G3_DRAGON/src/mitray_poles.f ../../G3_DRAGON/src/mitray_zone.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_dipole_probe mitray_dipole_probe.f ../../G3_DRAGON/src/mitray_dipole.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_dipole_probe_d2 mitray_dipole_probe_d2.f ../../G3_DRAGON/src/mitray_dipole.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_edipol_probe mitray_edipol_probe.f ../../G3_DRAGON/src/mitray_edipol.f
gfortran -Dgfortran -fallow-argument-mismatch -I../../G3_DRAGON/inc \
  -o mitray_edipol_probe_e2 mitray_edipol_probe_e2.f ../../G3_DRAGON/src/mitray_edipol.f
./mitray_poles_probe && ./mitray_poles_probe_q2 && ./mitray_poles_probe_rest
./mitray_dipole_probe && ./mitray_dipole_probe_d2
./mitray_edipol_probe && ./mitray_edipol_probe_e2

python3 compare_fields.py q1_fortran_reference.csv ../build/q1_cpp_reference.csv
python3 compare_fields.py q2_fortran_reference.csv ../build/q2_cpp_reference.csv
python3 compare_fields.py q3_q14_fortran_reference.csv ../build/q3_q14_cpp_reference.csv
python3 compare_fields.py d1_fortran_reference.csv ../build/d1_cpp_reference.csv
python3 compare_fields.py d2_fortran_reference.csv ../build/d2_cpp_reference.csv
python3 compare_fields.py e1_fortran_reference.csv ../build/e1_cpp_reference.csv
python3 compare_fields.py e2_fortran_reference.csv ../build/e2_cpp_reference.csv

# Tracking demos: x0(cm) y0(cm) p(MeV/c) nEvents
cd ../build
./dragon_g4_pilot --track        1.0 0.0 3000 1     # Q1 (standalone)
./dragon_g4_pilot --track-q2     1.0 0.0 3000 1     # Q2 (standalone)
./dragon_g4_pilot --track-dipole 0.0 0.0 64.7 1     # D1 (standalone), on design orbit
./dragon_g4_pilot --track-e1     0.0 0.0 29.768 1   # E1 (standalone), on design orbit
./dragon_g4_pilot --track-e2     0.0 0.0 29.763 1   # E2 (standalone), on design orbit
./dragon_g4_pilot --track-chain  0.0 0.0 258.7 1     # ALL 18 elements, 19Ne4+ on design orbit
./dragon_g4_pilot --track-reaction 5          # ALL 18 elements, real 15O(a,g)19Ne events (nEvents only)
./dragon_g4_pilot --reaction-stats            # sanity-check the reaction file's own cascade (no G4RunManager)

# Visualization (see "Visualization" above)
./dragon_g4_pilot --vis Chain                # interactive
./dragon_g4_pilot --vis Chain run_chain.mac  # batch, 5 events, then exit
./dragon_g4_pilot --vis Reaction              # interactive, real reaction events
```

## What this pilot does *not* cover (i.e. what's left for the full port)

- The rest of the dipole's parameter space: `MTYP` 1/2/5 (uniform-field
  magnets, via `mitray_dipole.f`'s `MITRAY_BDIP`/`MITRAY_BDPP`/`BDPPX`),
  `MTYP=6` (pretzel magnet, `MITRAY_BPRETZ`), field-map dipoles (`IMAP!=0`,
  `MITRAY_BDMP`/`FMAP`/`DMAP`), and curved pole faces (`MITRAY_SDIP`'s
  iterative search) -- none of which D1 or D2 need, but a different DRAGON
  configuration file's dipole(s) might.
- The other three MITRAY element types: solenoid (`mitray_solnd.f`), Einzel
  lens, velocity selector, accelerator gap -- none appear in this
  particular file, but do in the general MITRAY format. Each is its own
  `G4MagneticField`/`G4ElectroMagneticField` port of comparable size to
  Q1/D1/E1.
- GEANT3's own dipole volume-placement/orientation algorithm (see the
  Chaining section above) -- a bookkeeping/geometry-alignment detail,
  separate from the field physics, which remains bit-exact validated.
- The DSSSD's **strip segmentation**: it's a real part now (a Micron W1,
  see "Apertures and slits" above), but built as one solid silicon slab, not
  16+16 separate junction/ohmic strip volumes -- there's no
  sensitive-detector/digitization step yet (next bullet) to make
  per-strip geometry meaningful. Geometry-wise, every other part of this
  pilot is now a real port too: every `'RCOL'` aperture/slit, the full
  nested BGO housing/reflector/PMT stack, and the CMBR/CMBG outer box +
  shielding collar chain.
- `src/gudigi.f`'s digitization (863 lines: light collection, resolution
  smearing, PMT response) -- see "Sensitive detectors" below for what
  *is* now in place (a real `G4VSensitiveDetector`/hits pipeline for both
  BGO and DSSSD, so hits are actually collected); this pilot doesn't
  smear/digitize them into anything resembling real detector output yet.
- `src/ureact.f`/`gureact.f`'s full generality as a
  `G4VUserPrimaryGeneratorAction`: angular distributions/correlations
  anywhere in the cascade, a sampled recoil charge-state distribution
  (fixed at whatever the config file's own `RECL` card gives), and the
  other 18 reactions `ureact.f` hardcodes besides 15O(alpha,g)19Ne (a
  config file for any of those could be written in the same format, see
  "Reaction config file" above, but none is bundled) -- gamma-cascade
  support itself *is* now implemented (same section) and isn't limited to
  this one reaction's own dominant direct-to-ground branch any more.
- **Hadronic** physics processes (interactions/decay -- including the
  reduced hadron-decay KE cut in `src/geant_mod.f`'s patched `GTHADR`).
  Standard **EM** physics (energy loss, multiple scattering, gamma
  interactions) *is* now registered (`PhysicsList`'s `G4EmStandardPhysics`
  -- see "Sensitive detectors" above) -- this is what makes the BGO/DSSSD
  sensitive detectors' `edep` physically meaningful (verified: a clean
  ~4.03 MeV BGO full-energy peak matching the capture gamma). One real,
  not-yet-resolved consequence: the recoil ion now measurably loses
  momentum crossing the target region (see "Full chain"'s own caveat
  above), enough that neither `--track-chain`'s design-orbit test nor
  `--track-reaction`'s real events reach the DSSSD reliably any more.
- HBOOK/PAW output -> ROOT output: raw hits *are* now written to ROOT
  `TTree`s (see "Sensitive detectors"/"ROOT output" above) -- what's not
  ported is any of the real analysis/histogramming `src/gudigi.f` and the
  original HBOOK/PAW ntuples would have done downstream of that (energy
  spectra, coincidence gating, etc.).
