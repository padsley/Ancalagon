"""Ancalagon vs GEANT3 cross-check, 15O(a,g)19Ne Er=503.6 keV, 2014 tune.

Reads data/ (see run_g3.sh / run_ancalagon.sh), writes results/summary.txt
and results/chart_data.json. Same conventions as ../g3_crosscheck_k39pg:
GEANT3 efficiencies per reacting event, positions in the global frame
relative to the final-drift axis x = -1024.717 cm.

Loss locations: Ancalagon reports the volume each recoil stopped in.
GEANT3 writes each stopped recoil's position to fort.4; each stop is
assigned to the nearest beamline element in Ancalagon's own chain
geometry (same global coordinates), then both are grouped the same way.
"""
import glob, json, os, re
import numpy as np
import uproot

HERE = os.path.dirname(os.path.abspath(__file__))
D, R = f"{HERE}/data", f"{HERE}/results"
os.makedirs(R, exist_ok=True)
XC, Z_FSLT, Z_DSSSD = -1024.717008, -6.075286, -72.37
RECOIL_M = "m=  17698"

def group(name):
    for g in ("QSLT", "MSLT", "FSLT"):
        if name.startswith(g): return g
    if name == "DSSSD": return "DSSSD"
    if re.match(r"(CEL|EX|PD|EN|GAS|TARG|CENT|HOLE)", name): return "Target and pumping"
    return "Other beamline"

# ---------------- geometry for GEANT3 stop assignment ----------------
geo = []
for l in open(f"{D}/ancalagon/chain_geometry.txt"):
    p = l.split(); geo.append((p[1], float(p[2].split("=")[1]), float(p[3].split("=")[1])))
gx = np.array([g[1] for g in geo]); gz = np.array([g[2] for g in geo])
def nearest_element(x, z):
    d = (gx - x) ** 2 + (gz - z) ** 2; i = int(np.argmin(d))
    return geo[i][0], float(np.sqrt(d[i]))

# ---------------- GEANT3 ----------------
H = [uproot.open(f)["h1000"].arrays(["react", "recdet"], library="np") for f in sorted(glob.glob(f"{D}/g3/dragon1_j*.root"))]
react = np.concatenate([h["react"] for h in H]) == 1
recdet = np.concatenate([h["recdet"] for h in H]) == 1
N3all, N3 = len(react), int(react.sum()); k3 = int((recdet & react).sum())
stops = []
for f in sorted(glob.glob(f"{D}/g3/fort4_j*.txt")):
    for l in open(f):
        v = l.split()
        if len(v) >= 3:
            try: stops.append((float(v[0]), float(v[1]), float(v[2])))
            except ValueError: pass
g3_groups = {}
g3_far = 0
for x, y, z in stops:
    n, d = nearest_element(x, z)
    if z < 80 and abs(x) < 30:  # target chamber and pumping tubes, before Q1 (z = 106.9)
        n = "EX"
    g3_groups[group(n)] = g3_groups.get(group(n), 0) + 1
    g3_far += d > 60

# ---------------- Ancalagon ----------------
fates, pre = [], []
for l in open(f"{D}/ancalagon/recoil_fate.txt"):
    p = l.split()
    if p[0] == "END": fates.append(p[8])
    elif p[0] == "PRE": pre.append([float(t) for t in p[2:8]])
N4 = len(fates); k4 = fates.count("DSSSD")
g4_groups = {}
for n in fates: g4_groups[group(n)] = g4_groups.get(group(n), 0) + 1
g4_side = {n: fates.count(n) for n in set(fates) if n[:4] in ("MSLT", "FSLT", "QSLT")}
P4 = np.array(pre); z4, x4, y4, px4, py4, pz4 = P4.T
def g4_x(z): return x4 + px4 / pz4 * (z - z4) - XC

# ---------------- GEANT3 global frame (final drift) ----------------
evs, cur = [], []
for l in open(f"{D}/g3_global/tail.txt"):
    if l.startswith("EVT"):
        if cur: evs.append(np.array(cur))
        cur = []
    else:
        p = l.split(); cur.append([float(p[2]), float(p[3]), float(p[4])])
if cur: evs.append(np.array(cur))
A3 = np.array([np.r_[e[-1], e[-1] - e[0]] for e in evs if len(e) > 1 and abs(e[-1, 2] - e[0, 2]) > 0.1])
def g3_x(z): return A3[:, 0] + A3[:, 3] / A3[:, 5] * (z - A3[:, 2]) - XC
g3_xp = 1000 * A3[:, 3] / -A3[:, 5]; g4_xp = 1000 * px4 / -pz4
zs = np.linspace(110, Z_DSSSD, 183)
rms3 = np.array([g3_x(z).std() for z in zs]); rms4 = np.array([g4_x(z).std() for z in zs])

# ---------------- summary ----------------
def pct(k, n): p = k / n; return 100 * p, 100 * np.sqrt(p * (1 - p) / n)
order = ["DSSSD", "MSLT", "FSLT", "QSLT", "Target and pumping", "Other beamline"]
out = [f"GEANT3: {N3all} generated, {N3} reacting; {len(stops)} recoil stop points; global-frame sample {len(A3)} recoils at the DSSSD.",
       f"Ancalagon: {N4} events.", "", "== Where recoils end up (% of reacting events) =="]
g3_groups["DSSSD"] = k3
tab = []
for g in order:
    a, b = pct(g3_groups.get(g, 0), N3), pct(g4_groups.get(g, 0), N4)
    tab.append(dict(label=g, g3=a, g4=b)); out.append(f"  {g:20s} GEANT3 {a[0]:6.2f} +- {a[1]:4.2f}   Ancalagon {b[0]:6.2f} +- {b[1]:4.2f}")
out.append(f"  (GEANT3 stops unassigned/far from any element: {g3_far}; GEANT3 stops total {len(stops)} vs reacting-but-lost {N3-k3})")
out.append("  Ancalagon slit jaws: " + ", ".join(f"{k} {v}" for k, v in sorted(g4_side.items())))
optics = [("x at FSLT, rms (cm)", g3_x(Z_FSLT).std(), g4_x(Z_FSLT).std()),
          ("x at FSLT, mean (cm)", g3_x(Z_FSLT).mean(), g4_x(Z_FSLT).mean()),
          ("x at DSSSD, rms (cm)", g3_x(Z_DSSSD).std(), g4_x(Z_DSSSD).std()),
          ("x' at DSSSD, mean (mrad)", g3_xp.mean(), g4_xp.mean()),
          ("x waist position z (cm)", zs[np.argmin(rms3)], zs[np.argmin(rms4)])]
out += ["", "== Final-focus optics, transmitted recoils (global frame) =="] + [f"  {a:28s} GEANT3 {b:8.3f}   Ancalagon {c:8.3f}" for a, b, c in optics]
open(f"{R}/summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
json.dump(dict(fates=tab, optics=[dict(label=a, g3=round(float(b), 4), g4=round(float(c), 4)) for a, b, c in optics],
               waist=dict(dist=(Z_FSLT - zs).round(2).tolist(), g3=rms3.round(4).tolist(), g4=rms4.round(4).tolist()),
               n=dict(g3_all=N3all, g3=N3, g3_global=len(A3), g4=N4), g4_jaws=g4_side),
          open(f"{R}/chart_data.json", "w"))
