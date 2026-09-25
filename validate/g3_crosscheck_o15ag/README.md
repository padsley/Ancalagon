# GEANT3 cross-check: 15O(α,γ)19Ne, 2014 tune (2026-09-25)

Same comparison as `../g3_crosscheck_k39pg`, for the bundled
`reactions/o15ag_19ne.reaction` (Er = 503.6 keV, 19Ne<sup>4+</sup>, 5 Torr He),
at the real 2014 hardware tune (MSLT 0.75 × 1.25 cm) on both sides.

```
./run_g3.sh <workdir>          # 8 x 1250 GEANT3 events (+800 global-frame); copy outputs into data/g3*
./run_g3_planes.sh <workdir>   # 4 x 400 GEANT3 recoil tracks near target, MSLT, FSLT
./run_ancalagon.sh <workdir>   # 10k reaction events, track captures, scatter-free angle fan
python3 analyze.py             # where recoils end up; final-drift optics
python3 plane_fit.py           # population (x|theta) at MSLT and FSLT, both codes
python3 fan.py                 # Ancalagon single-ray fan through Q1/QSLT/MSLT/FSLT
```

## GEANT3 input

GEANT3's built-in reaction 2 (15O(α,γ)19Ne) cannot be used: in `ureact.f`
its 0.275 and 0.238 MeV states decay to particle 95, which is never
defined (the 19Ne ground state is particle 85), so about 20 % of events
never produce a proper recoil. `o15ag.dat` is a file-driven (`FKIN 20 4.0 4.0`)
copy with the same masses, levels, lifetimes and branching ratios as both
the built-in reaction and Ancalagon's reaction file, and Ancalagon's total
width (50.6 eV). Both codes tune the separator to the recoil momentum after
the target gas; GEANT3's magnetic scale factor is 0.9736, Ancalagon's RTUN
0.9641 (1 % apart, against 0.07 % for k39pg).

## Results

| Where recoils end (% of reacting events) | GEANT3 | Ancalagon |
|---|---|---|
| DSSSD | 82.3 | 51.2 |
| MSLT | 0.9 | 39.8 |
| FSLT | 0.1 | 4.6 |
| QSLT | 8.2 | 1.0 |
| Target and pumping | 7.1 | 2.8 |
| Other beamline | 1.5 | 0.7 |

| (x\|θ) from each code's own recoils | GEANT3 | Ancalagon |
|---|---|---|
| MSLT slope | +0.03 mm/mrad | +0.79 mm/mrad |
| MSLT x rms | 0.09 cm | 0.43 cm |
| FSLT slope | +0.04 mm/mrad | +1.91 mm/mrad |
| FSLT x rms | 0.25 cm | 1.32 cm |

GEANT3 forms angle-independent foci at MSLT and FSLT. Ancalagon does not:
its MSLT losses rise with recoil angle (0.6 % below 4 mrad, 55 % above
12 mrad), and the 19Ne cone reaches about 16 mrad, against about 7 mrad for
40Ca in k39pg. That is the whole transmission gap. GEANT3's extra QSLT and
target-region losses are consistent with its finite beam spot and
divergence, which Ancalagon does not model.

The scatter-free Ancalagon fan (`fan.py`) is linear at Q1 and shows
(x|θ) ≈ +0.8 mm/mrad at MSLT, with the on-axis ray 1.9 mm off centre.

## Method notes

- Single-ray fans must start past the dense target gas. Fired from the
  target centre, each ray picks up a random 2–3 mrad multiple-scattering
  kick in its first step, which can make a fan look discontinuous.
- GEANT3 prints recoil steps (FOCUSTEST) only from Q1 onward (z ≈ 88 cm), so
  its angle is taken from its first steps there; Ancalagon's over z < 60 cm.
  Both codes take long straight steps through the field-free drifts, so
  slit-plane positions are interpolated along those straight segments.
- Raw track captures (`data/*/planes.txt`, up to 143 MB) are gitignored;
  `plane_fit.py` reduces them to `plane_points.csv` on first read.
