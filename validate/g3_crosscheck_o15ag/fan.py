"""Ancalagon single-ray 19Ne4+ angle fan: x where each ray crosses Q1,
QSLT, MSLT and FSLT (slit-local frame), from data/fan/f.<mrad>.txt.

The rays start at z0 = 74 cm, x0 = 74 cm * tan(theta), i.e. as if from the
target centre, but past the dense target gas: fired from inside the gas,
multiple scattering gives each single ray a random 2-3 mrad kick, which
is enough to make a fan look discontinuous. Slits opened to +-20 cm so no
ray is stopped before FSLT. See run_ancalagon.sh for the commands.
"""
import glob, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PLANES = {"Q1": (0.0, 106.885, 0), "QSLT": (-59.376, 351.542243, -50),
          "MSLT": (-498.206223, 657.134274, -70), "FSLT": (-1024.717008, -6.075286, -180)}
angs = sorted(int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{HERE}/data/fan/f.*.txt"))
print("mrad   " + "  ".join(f"{k:>8s}" for k in PLANES) + "   last volume")
res = {}
for a in angs:
    rows = [l.split() for l in open(f"{HERE}/data/fan/f.{a}.txt")]
    P = np.array([[float(r[1]), float(r[2])] for r in rows])  # z, x
    out = {}
    for k, (xm, zm, th) in PLANES.items():
        t = math.radians(th); d = np.array([math.sin(t), math.cos(t)]); tr = np.array([math.cos(t), -math.sin(t)])
        rel = np.c_[P[:, 1] - xm, P[:, 0] - zm]; s = rel @ d; u = rel @ tr
        idx = np.where((s[:-1] < 0) & (s[1:] >= 0) & (np.abs(u[:-1]) < 30))[0]
        if len(idx): i = idx[0]; f = -s[i] / (s[i + 1] - s[i]); out[k] = float(u[i] + f * (u[i + 1] - u[i]))
    res[a] = out
    print(f"{a:4d}   " + "  ".join(f"{out[k]:8.3f}" if k in out else "       -" for k in PLANES) + "   " + rows[-1][7])
json_out = {str(a): v for a, v in res.items()}
import json; json.dump(json_out, open(f"{HERE}/results/fan.json", "w"), indent=1)
