# GEANT3 cross-check: 39K(p,g)40Ca, 2014 tune (2026-09-25)

Full-distribution comparison of Ancalagon against the GEANT3 DRAGON
simulation (`G3_DRAGON_Claude`, `dat/k39pg.dat`), at the real 2014 hardware
tune (MSLT 0.75 x 1.25 cm), on Ancalagon `main` after the quad-container fix
and the Q14 steering correction.

```
./run_g3.sh <workdir>          # 8 x 1250 GEANT3 events + a 400-event global-frame job
./run_ancalagon.sh <workdir>   # 10k reaction events, trajectory extracts, angle fans
python3 analyze.py             # -> results/summary.txt, compare.png, chart_data.json
```

`run_g3.sh` leaves its output in `<workdir>`; copy `j*/dragon1.root` to
`data/g3/dragon1_j*.root` and `jg/{endv,tail}.txt` to `data/g3_global/`.
`run_ancalagon.sh` writes into `data/` directly. `*.root` is gitignored, so the
ROOT files in `data/` stay local; everything else here is small text.

## Results

| | GEANT3 | Ancalagon |
|---|---|---|
| Recoil reaches DSSSD | 99.39 ± 0.08 % | 99.97 ± 0.02 % |
| BGO any crystal hit | 71.1 % | 71.9 % |
| BGO highest crystal > 2 MeV | 61.1 % | 62.1 % |
| BGO full-energy sum | 30.3 % | 31.9 % |
| x rms at FSLT | 0.125 cm | 1.069 cm |
| x rms at DSSSD | 0.318 cm | 0.890 cm |
| mean x' at DSSSD | −1.2 mrad | −11.3 mrad |
| x waist | at FSLT | none before the DSSSD |

Transmission and the BGO array agree. The final horizontal focus does not:
Ancalagon has (x|θ) = −2.5 mm/mrad at FSLT, confirmed by a fit to 3k real
recoils (residual 0.11 cm) and by single-ray fans of 40Ca8+ and 19Ne4+.

## Conventions and gotchas

- GEANT3 efficiencies are per reacting event (HISTORY `react == 1`). About 5 %
  of GEANT3 events never react; counting them gives a false 94.5 %.
- BGO: GEANT3's `gudigi.f` digitisation is applied to Ancalagon's hits (sum
  per crystal, 0.1 MeV per-crystal threshold, highest crystal first).
- Compare DSSSD/FSLT positions in global coordinates. GEANT3's ENDV
  histograms (h11–h14) are centred on its own ENDV volume, about 3.6 mm from
  Ancalagon's DSSSD centre (x = −1024.717 cm). h15 (recoil energy) spans
  0–20 MeV and misses these ~23 MeV recoils.
- Ancalagon fires recoils from a point with no beam emittance, so its y spot
  is smaller than GEANT3's (which has a 1.16 mm rms beam spot, 1.1 mrad
  divergence and a 1.1 cm rms spread of interaction points in the gas). This
  is an input difference, not optics.
- GEANT3 run notes: FFCARD/INPUT/MITRAY paths are truncated at 80
  characters; the FOCUSTEST print in `src/gustep_mitray.f` fires on every
  recoil step; `dsbatch` runs about 1.15 events/s per process.
