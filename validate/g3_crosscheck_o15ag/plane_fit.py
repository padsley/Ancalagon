"""Population (x|theta) at MSLT and FSLT for both codes, from recoil tracks
captured near the target and near each slit (data/*/planes.txt, written by
run_g3_planes.sh and run_ancalagon.sh, reduced on first read to
data/*/plane_points.csv, which is what the repo keeps). Same method for both codes:

  theta : horizontal angle of the recoil's track just after the target,
          from a straight fit: Ancalagon over z < 60 cm (after the target
          gas), GEANT3 over z = 88-92 cm (its first printed steps, at the
          upstream edge of Q1's fringe field)
  x     : transverse position where the track crosses the slit plane,
          in the slit's local frame (same sign convention for both codes)

Only recoils that actually cross a plane contribute, so the fit covers the
angles each code transmits to that point.
"""
import json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLANES = {"MSLT": (-498.206223, 657.134274, -70.0), "FSLT": (-1024.717008, -6.075286, -180.0)}

def load(path):
    evs, cur = [], []
    for l in open(path):
        if l.startswith("EVT"):
            if cur: evs.append(np.array(cur))
            cur = []
        else:
            p = l.split(); cur.append([float(p[2]), float(p[3]), float(p[4])])  # x y z
    if cur: evs.append(np.array(cur))
    return evs

def reduce(evs):
    """Per recoil: theta and, where crossed, x at each slit plane."""
    th = []
    for e in evs:
        t = e[(e[:, 2] < 92) & (np.abs(e[:, 0]) < 5)]
        if t[0, 2] > 80: t = t[t[:, 2] < t[0, 2] + 4]      # GEANT3: first ~4 cm from z = 88
        if len(t) < 3 or np.ptp(t[:, 2]) < 2: continue
        a = 1000 * np.polyfit(t[:, 2], t[:, 0], 1)[0]
        row = {"th": a}
        for k, (xm, zm, deg) in PLANES.items():
            r = math.radians(deg); d = np.array([math.sin(r), math.cos(r)]); tr = np.array([math.cos(r), -math.sin(r)])
            rel = np.c_[e[:, 0] - xm, e[:, 2] - zm]; s = rel @ d; u = rel @ tr
            near = np.hypot(rel[:, 0], rel[:, 1]) < 450   # both ends near this slit (straight drift)
            idx = np.where((s[:-1] < 0) & (s[1:] >= 0) & near[:-1] & near[1:])[0]
            if len(idx):
                i = idx[0]; f = -s[i] / (s[i + 1] - s[i]); row[k] = u[i] + f * (u[i + 1] - u[i])
        th.append(row)
    return th

def analyse(th):
    out = {}
    for k in PLANES:
        pts = np.array([[r["th"], r[k]] for r in th if k in r])
        if len(pts) < 20: out[k] = None; continue
        m = np.abs(pts[:, 0]) < 8
        c1 = np.polyfit(pts[m, 0], pts[m, 1], 1)
        c2 = np.polyfit(pts[:, 0], pts[:, 1], 2)
        out[k] = dict(n=len(pts), slope_mm_per_mrad=10 * c1[0], offset_cm=c1[1], quad_cm_per_mrad2=c2[0],
                      rms_cm=float(pts[:, 1].std()), resid_after_linear_cm=float((pts[m, 1] - np.polyval(c1, pts[m, 0])).std()),
                      theta_rms_mrad=float(pts[:, 0].std()),
                      binned=[(lo, float(pts[(pts[:, 0] >= lo) & (pts[:, 0] < lo + 4), 1].mean()) if ((pts[:, 0] >= lo) & (pts[:, 0] < lo + 4)).sum() > 5 else None)
                              for lo in range(-16, 16, 4)])
    out["n_events"] = len(th)
    return out

def rows_for(d):
    """Reduce data/<d>/planes.txt (large; regenerable) to plane_points.csv, or read the CSV."""
    raw, csv = f"{HERE}/data/{d}/planes.txt", f"{HERE}/data/{d}/plane_points.csv"
    if os.path.exists(raw):
        th = reduce(load(raw))
        with open(csv, "w") as f:
            f.write("theta_mrad,x_mslt_cm,x_fslt_cm\n")
            for r in th: f.write(f"{r['th']:.4f},{r.get('MSLT', float('nan')):.5f},{r.get('FSLT', float('nan')):.5f}\n")
        return th
    a = np.genfromtxt(csv, delimiter=",", skip_header=1)
    return [dict(th=t, **({"MSLT": m} if m == m else {}), **({"FSLT": f} if f == f else {})) for t, m, f in a]

res = {code: analyse(rows_for(d)) for code, d in (("GEANT3", "g3_global"), ("Ancalagon", "ancalagon"))}
for k in PLANES:
    print(f"== {k} ==")
    for code in res:
        r = res[code][k]
        if r is None: print(f"  {code}: too few crossings"); continue
        print(f"  {code:9s} n={r['n']:5d}  (x|th) {r['slope_mm_per_mrad']:+6.2f} mm/mrad (|th|<8)  offset {r['offset_cm']:+.3f} cm  "
              f"quad {r['quad_cm_per_mrad2']:+.5f} cm/mrad^2  rms {r['rms_cm']:.3f} cm  resid {r['resid_after_linear_cm']:.3f} cm  th rms {r['theta_rms_mrad']:.1f}")
        print("            mean x per 4-mrad bin from -16: " + " ".join("   -  " if v is None else f"{v:+.3f}" for _, v in r["binned"]))
json.dump(res, open(f"{HERE}/results/plane_fit.json", "w"), indent=1)
