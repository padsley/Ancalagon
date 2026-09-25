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

| (x\|θ) mm/mrad | GEANT3 | before | origin fix | both fixes |
|---|---|---|---|---|
| RC9 | +2.691 | +3.213 | +2.718 | +2.691 |
| QSLT | −0.004 | −0.105 | +0.063 | +0.015 |
| RC23 | −2.746 | −3.148 | −3.339 | −2.693 |
| MSLT | +0.007 | +0.735 | −0.769 | +0.024 |
| FSLT | −0.021 | +2.592 | −1.574 | −0.192 |

Still open: after D2 a smaller residual remains ((a|θ) 9 % high at RC40,
−0.19 mm/mrad at FSLT), and the on-axis ray leaves D1 about 1.6 mrad
different from GEANT3's (1.9 mm off at MSLT). The Q14 steering correction
(0.21) and the D2 residual trim were calibrated against the buggy optics and
need re-deriving.
