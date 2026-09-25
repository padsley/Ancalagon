# Element-by-element fan: where the MSLT focus was lost

Matched single-ray 19Ne<sup>4+</sup> fans in both codes, fired from z = 90 cm
(inside Q1, past the target gas; x0 = 90 cm · tan θ) at each code's own
tuned rigidity, compared at every field-free plane from Q1 to FSLT.

```
./run_fans.sh <workdir>      # GEANT3 (KINE rays) -> g3/, Ancalagon -> ancalagon/
python3 compare.py [dir]     # dir: ancalagon (current), ancalagon_origin_fix, ancalagon_before
./run_section.sh <workdir>   # rays launched at QSLT with GEANT3's state there
```

Two Ancalagon bugs, found in this order:

1. **Quad field origin** (`ChainQuad`). `MitrayQuadrupoleField` measures local
   z from the element's entrance reference plane (`zb = A - za`), as GEANT3
   does (`ugeom_mitray.f` offsets each POLE by `(Z11 - Z22 - L)/2` "to position
   the A-frame at the edge of the entrance EFB"). Ancalagon passed the
   container centre instead, shifting every quad's field downstream by
   A + (L + Z22 - Z11)/2 (12.6 cm for Q1). Present since the initial commit;
   the bit-exact field-value validation compared local coordinates, so it
   could not see a placement error. Fixing it makes Q1–Q2 match GEANT3 at
   RC9 and QSLT (`ancalagon_origin_fix`).

2. **No field superposition in overlaps.** GEANT3's `gufld.f` adds the fields
   of every overlapping MANY volume. Geant4 gives a point to one volume, so
   wherever neighbouring quad containers overlap (Q3/Q4, Q4/Q5, ...) one
   fringe was dropped. The Q3–Q7 section still disagreed after fix 1, even
   for rays launched at QSLT with GEANT3's own state (`run_section.sh`,
   `*_qslt_origin_fix`). `ChainMagneticFieldSum` now sums every chain
   magnetic element whose container holds the point
   (`FIELD_SUPERPOSITION=0` restores the old behaviour).

3. **Dipole fringe truncation** (follow-up, same day). D1's and D2's
   containers are capped to keep clear of neighbouring collimators, which
   cut ~14 cm off D1's 30 cm fringes (its on-axis ray left D1 1.6 mrad off
   GEANT3's). Each dipole and e-dipole now registers an unplaced
   natural-reach field region with the sum, and the World carries the
   summed field (a `G4ElectroMagneticField`, so B and E add), so fringes act
   in the gaps without any container swallowing a slit. The target chamber
   and BGO volumes are reset to no field.

4. **Corrections tuned on the old optics.** `kD2ResidualTrim` (2.6 %)
   compensated D2's truncated fringe and now over-bends (on-axis ray 3.6 cm
   off at RC40), so it applies only with `FIELD_SUPERPOSITION=0`. The Q14
   steering default (0.21) moves the beam away from GEANT3 once the rest is
   fixed (k39pg_40ca mean angle at the DSSSD −9.6 mrad with it, −3.3 without,
   GEANT3 −1.2) and no single value matches both reactions, so its default
   is now 0 (0.21 in legacy mode; `Q14_STEER_COEFF` still overrides).

| (x\|θ) mm/mrad | GEANT3 | before | origin fix | + superposition | all fixes (final) |
|---|---|---|---|---|---|
| RC9 | +2.691 | +3.213 | +2.718 | +2.691 | +2.691 |
| QSLT | −0.004 | −0.105 | +0.063 | +0.015 | −0.003 |
| RC23 | −2.746 | −3.148 | −3.339 | −2.693 | −2.744 |
| MSLT | +0.007 | +0.735 | −0.769 | +0.024 | +0.006 |
| RC40 | −0.100 | −0.838 | +0.765 | +0.058 | −0.102 |
| FSLT | −0.021 | +2.592 | −1.574 | −0.192 | +0.052 |

The final column uses fine-stepped Ancalagon rays (`FINE_STEP_CM=0.5`, now
in `run_fans.sh`); the earlier columns were coarse-stepped, which is exact
only where the World carried no field. (a|θ) also agrees to <1 % through
RC40, 5–8 % after E2.

## E1/E2: the remaining difference is GEANT3's electrostatic transport

The "6 cm earlier E2 bend" was an artifact of measuring path length from
RC40: measured from each code's own crossing of E2's entrance EFB, the
curvature profiles along the on-axis ray agree (full field 1/250 cm⁻¹ in
both; Ancalagon's entrance fringe up to ~10 % stronger at the EFB). GEANT3's
e-dipole A-frame (`ugeo_edipol`) works out exactly at the deck's entrance
point and direction when Z11 = Z22, as Ancalagon places it.

Section tests (`run_section_plane.sh`, `section_matrix.py`) launch identical
rays at a field-free plane in both codes -- GEANT3's own on-axis state plus
±0.2 cm and ±2 mrad -- and measure the element's (x, a) transfer matrix, with
angles from two field-free planes (Ancalagon prints z to 1e-4 cm, too coarse
for single-step directions). Phase-space conservation requires det ≈ 1:

| section (launch → planes) | det GEANT3 | det Ancalagon | largest difference |
|---|---|---|---|
| Q8–Q10 + D2 (MSLT → RC40, RC42) | 0.980 | 0.999 | none (≤ 0.5 %) |
| E1 (RC25 → MSLT, RC28) | 1.075 | 0.999 | (a\|a) 0.887 vs 0.811 |
| E2 (RC49 → RC51, RC55) | 1.305 | 0.999 | (a\|a) 0.837 vs 0.539 |

Magnetic transport agrees; GEANT3's electrostatic transport does not conserve
phase space. Its energy bookkeeping is also off: with a temporary momentum
column in FOCUSTEST (reverted), a ray entering E2 3.19 mm off-orbit should
lose qE0·u ≈ 4.5 keV; Ancalagon gives −4.2 keV, GEANT3 −1.4 keV. GEANT3's
gradient across the gap is right (13.8 keV/cm) but its zero sits ~2 mm off
the design orbit. The likely source is `grkuta.f`'s 1997 electric-field
extension (a Nyström step treating the E impulse as a direction kick, with
the momentum magnitude averaged from two stage estimates) together with
`gthion.f` applying the momentum change only when it exceeds 1 keV/c per step.
Not traced further. The E1/E2 differences, the on-axis offset after E2
(5.7 mm at RC55) and the residual at FSLT are therefore not Ancalagon bugs.
