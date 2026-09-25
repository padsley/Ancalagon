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

Still open: the on-axis ray is 0.6 mm apart at RC40 and 5.7 mm apart by
RC55. Through E2 its direction runs a constant 24 mrad ahead of GEANT3's,
i.e. Ancalagon's E2 bend starts ~6 cm earlier along the path, although both
codes put E2's entrance EFB at the same chain position and bend at the same
radius. GEANT3 builds its e-dipole volume and A-frame from a trapezoid
construction (`ugeo_edipol` in `ugeom_mitray.f`) that may place the A-frame
differently from the deck; not yet traced. E1 (smaller bend) agrees.
