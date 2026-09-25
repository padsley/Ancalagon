"""Element-by-element comparison of matched single-ray 19Ne4+ fans
(run_fans.sh) between Q1 and FSLT.

At each reference plane (an element position from Ancalagon's chain
geometry, local frame), finds where each ray crosses it and its local
angle there, then fits against the initial angle theta:
    u  = u0  + (x|theta) * theta
    u' = u0' + (a|theta) * theta
A plane is only trusted when the crossing segment lies in field-free
space in both codes: Ancalagon prints one point per step, and inside a
field volume one step can span a whole element.
"""
import glob, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G4DIR = sys.argv[1] if len(sys.argv) > 1 else "ancalagon"   # or ancalagon_before, ancalagon_origin_fix
GEO = {}
for l in open(f"{HERE}/../data/ancalagon/chain_geometry.txt"):
    p = l.split(); GEO[p[1]] = (float(p[2].split("=")[1]), float(p[3].split("=")[1]), float(p[4].split("=")[1]))
FIELD = ("Q", "D1", "D2", "E1", "E2")   # volume-name prefixes that carry a field

def load(code, a):
    pts = []
    for l in open(f"{HERE}/{G4DIR if code == 'ancalagon' else code}/ray.{a}.txt"):
        p = l.split()
        if code == "g3": pts.append((float(p[2]), float(p[4]), p[1]))            # x, z, volume
        else: pts.append((float(p[2]), float(p[1]), p[7]))                        # TRAJ: z x y ...
    return pts

def crossing(pts, name):
    xm, zm, th = GEO[name]
    r = math.radians(th); d = np.array([math.sin(r), math.cos(r)]); tr = np.array([math.cos(r), -math.sin(r)])
    P = np.array([[x, z] for x, z, _ in pts]); rel = P - [xm, zm]
    s, u = rel @ d, rel @ tr
    for i in range(len(pts) - 1):
        if s[i] < 0 <= s[i + 1] and abs(u[i]) < 50 and np.hypot(*rel[i]) < 500:
            f = -s[i] / (s[i + 1] - s[i]); seg = P[i + 1] - P[i]
            ang = 1000 * math.atan2(seg @ tr, seg @ d)
            vols = (pts[i][2], pts[i + 1][2])
            return u[i] + f * (u[i + 1] - u[i]), ang, vols, float(np.hypot(*seg))
    return None

angles = sorted(int(f.rsplit(".", 2)[1]) for f in glob.glob(f"{HERE}/g3/ray.*.txt"))
planes = ["Q1", "RC9", "D1", "QSLT", "RC12", "Q3", "Q4", "Q5", "Q6", "Q7", "RC23", "RC25", "FC1", "E1", "E1_exit",
          "FC2", "RC27", "MSLT", "RC28", "RC30", "Q8", "D2_exit", "RC40", "E2_exit", "RC55", "FSLT"]
print(f"{'plane':8s} | {'GEANT3 (x|th)':>13s} {'(a|th)':>8s} {'u0':>7s} | {'Ancalagon (x|th)':>16s} {'(a|th)':>8s} {'u0':>7s} | segment volumes (G3 / G4), G4 step cm")
rows = []
for name in planes:
    if name not in GEO: continue
    res = {}
    for code in ("g3", "ancalagon"):
        c = [(a, crossing(load(code, a), name)) for a in angles]
        c = [(a, v) for a, v in c if v is not None]
        if len(c) < 3: res[code] = None; continue
        A = np.array([a for a, _ in c]); U = np.array([v[0] for _, v in c]); T = np.array([v[1] for _, v in c])
        res[code] = dict(xt=10 * np.polyfit(A, U, 1)[0], at=np.polyfit(A, T, 1)[0], u0=float(np.interp(0, A, U)),
                         vols=c[len(c) // 2][1][2], step=c[len(c) // 2][1][3])
    g, m = res["g3"], res["ancalagon"]
    fmt = lambda r: "            -        -       -" if r is None else f"{r['xt']:+13.3f} {r['at']:+8.3f} {r['u0']:+7.3f}"
    fmt4 = lambda r: "               -        -       -" if r is None else f"{r['xt']:+16.3f} {r['at']:+8.3f} {r['u0']:+7.3f}"
    vol = (("/".join(g["vols"]) if g else "-") + "  " + ("/".join(m["vols"]) if m else "-") + (f"  {m['step']:.0f}" if m else ""))
    print(f"{name:8s} | {fmt(g)} | {fmt4(m)} | {vol}")
    rows.append(dict(plane=name, g3=g, g4=m))
print("\n(x|th) in mm/mrad, (a|th) in mrad/mrad, u0 = on-axis ray offset in cm (slit-local frame).")
