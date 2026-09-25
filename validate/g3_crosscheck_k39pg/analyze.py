"""Ancalagon (GEANT4) vs GEANT3 DRAGON cross-check, 39K(p,g)40Ca Er=606 keV, 2014 tune.

Reads data/ (see run_g3.sh / run_ancalagon.sh) and writes results/summary.txt,
results/compare.png and results/chart_data.json.

Conventions:
  * GEANT3 efficiencies are per *reacting* event (HISTORY react==1). ~5% of
    GEANT3 events never react (beam crosses the target off-resonance) and
    must not count as separator losses.
  * BGO: GEANT3's gudigi.f digitisation is replicated for Ancalagon (sum
    per crystal, drop crystals < 0.1 MeV, highest crystal = "first").
  * Positions/angles near the DSSSD are compared in GLOBAL coordinates,
    relative to Ancalagon's DSSSD centre x = -1024.717 cm. GEANT3's own
    ENDV histograms (h11-h14) are centred on its ENDV volume, which sits
    ~3.6 mm away, so they are not used for centroids.
  * Everything between z = +110 cm (end of Q14's fringe field) and the
    DSSSD (z = -72.37 cm) is field-free drift, so straight-line
    extrapolation is exact there. The beam travels toward -z.
"""
import glob, json, os
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D, R = f"{HERE}/data", f"{HERE}/results"
os.makedirs(R, exist_ok=True)
XC = -1024.717008          # Ancalagon DSSSD / FSLT axis x (cm)
Z_FSLT, Z_DSSSD = -6.075286, -72.37
EG = 8.934                 # single gamma, Ex(40Ca) -> g.s.

# ---------------- GEANT3 (10k, ntuples) ----------------
H, G = [], []
for rf in sorted(glob.glob(f"{D}/g3/dragon1_j*.root")):
    f = uproot.open(rf)
    H.append(f["h1000"].arrays(["react", "recdet"], library="np"))
    G.append(f["h1001"].arrays(["e_bgos_total", "e_bgo_first", "num_bgos_hit"], library="np"))
g3 = {k: np.concatenate([h[k] for h in H]) for k in H[0]}
g3g = {k: np.concatenate([g[k] for g in G]) for k in G[0]}
react = g3["react"] == 1
N3all, N3 = len(react), int(react.sum())
g3_det = (g3["recdet"] == 1) & react
g3_first, g3_tot, g3_n = g3g["e_bgo_first"][react], g3g["e_bgos_total"][react], g3g["num_bgos_hit"][react]

# ---------------- GEANT3 (400, global coordinates) ----------------
evs, cur = [], []
for l in open(f"{D}/g3_global/tail.txt"):
    if l.startswith("EVT"):
        if cur: evs.append(np.array(cur))
        cur = []
    else:
        p = l.split(); cur.append([float(p[2]), float(p[3]), float(p[4])])
if cur: evs.append(np.array(cur))
# last two recorded points (TST4 marker, 0.2 cm apart, field-free) give the line
A3 = np.array([np.r_[e[-1], e[-1] - e[0]] for e in evs if len(e) > 1 and abs(e[-1, 2] - e[0, 2]) > 0.1])
def g3_x(z): return A3[:, 0] + A3[:, 3] / A3[:, 5] * (z - A3[:, 2]) - XC
def g3_y(z): return A3[:, 1] + A3[:, 4] / A3[:, 5] * (z - A3[:, 2])
g3_xp = 1000 * A3[:, 3] / -A3[:, 5]; g3_yp = 1000 * A3[:, 4] / -A3[:, 5]

# ---------------- Ancalagon (10k) ----------------
f4 = uproot.open(f"{D}/ancalagon/dragon_hits.root")
bgo = f4["Bgo"].arrays(["eventID", "crystalID", "edepMeV"], library="np")
N4 = 10000
n4_det = len(np.unique(f4["Dsssd"].arrays(["eventID"], library="np")["eventID"]))
e = np.zeros((N4, 31))
np.add.at(e, (bgo["eventID"].astype(int), bgo["crystalID"].astype(int)), bgo["edepMeV"])
e[e < 0.1] = 0
g4_tot, g4_first, g4_n = e.sum(1), e.max(1), (e > 0).sum(1)
P4 = np.array([[float(v) for v in l.split()[2:8]] for l in open(f"{D}/ancalagon/dsssd_traj.txt") if l.startswith("PRE")])
z4, x4, y4, px4, py4, pz4 = P4.T
def g4_x(z): return x4 + px4 / pz4 * (z - z4) - XC
def g4_y(z): return y4 + py4 / pz4 * (z - z4)
g4_xp = 1000 * px4 / -pz4; g4_yp = 1000 * py4 / -pz4

# ---------------- Ancalagon (x|theta) fit at FSLT ----------------
L = open(f"{D}/ancalagon/vp.txt").read().split("\n")
V = np.array([[float(t) for t in l.split()[2:8]] for l in L if l.startswith("V")])
Pp = np.array([[float(t) for t in l.split()[2:8]] for l in L if l.startswith("P")])
th = 1000 * V[:, 3] / V[:, 5]
p0 = np.linalg.norm(V[:, 3:6], axis=1); dp = 100 * (p0 - p0.mean()) / p0.mean()
xF = Pp[:, 1] + Pp[:, 3] / Pp[:, 5] * (Z_FSLT - Pp[:, 0]) - XC
M = np.c_[np.ones_like(th), th, dp, th**2, th * dp, dp**2]
coef, *_ = np.linalg.lstsq(M, xF, rcond=None)
fit_resid = (xF - M @ coef).std()

# ---------------- single-ray angle fans ----------------
def fan(prefix):
    out = []
    for fn in glob.glob(f"{D}/fan/{prefix}.*.txt"):
        a = int(fn.rsplit(".", 2)[1])
        rows = [l.split() for l in open(fn)]
        last = None
        for r in rows:
            z, x = float(r[1]), float(r[2])
            if last and last[0] > Z_FSLT >= z and x < -1000:
                t = (Z_FSLT - last[0]) / (z - last[0]); out.append((a, last[1] + t * (x - last[1]) - XC)); break
            last = (z, x)
    return sorted(out)
fan_ca, fan_ne = fan("ca"), fan("ne")
ca = np.array(fan_ca); ca_slope = np.polyfit(ca[:, 0], ca[:, 1], 1)[0]

# ---------------- waist scan over the field-free drift ----------------
zs = np.linspace(110, Z_DSSSD, 183)
rms3 = np.array([g3_x(z).std() for z in zs]); rms4 = np.array([g4_x(z).std() for z in zs])

# ---------------- summary ----------------
def eff(k, n): p = k / n; return 100 * p, 100 * np.sqrt(p * (1 - p) / n)
rows = []
def row(label, k3, k4):
    a, b = eff(k3, N3), eff(k4, N4); rows.append(dict(label=label, g3=a, g4=b))
row("Recoil reaches DSSSD", int(g3_det.sum()), n4_det)
row("BGO: any crystal hit", int((g3_first > 0.1).sum()), int((g4_first > 0.1).sum()))
row("BGO: highest crystal > 2 MeV", int((g3_first > 2).sum()), int((g4_first > 2).sum()))
row("BGO: highest crystal > 5 MeV", int((g3_first > 5).sum()), int((g4_first > 5).sum()))
row("BGO: full-energy sum > 8.63 MeV", int((g3_tot > EG - 0.3).sum()), int((g4_tot > EG - 0.3).sum()))
det4 = np.zeros(N4, bool); det4[np.unique(f4["Dsssd"].arrays(["eventID"], library="np")["eventID"]).astype(int)] = True
row("DSSSD and highest BGO > 2 MeV", int((g3_det[react] & (g3_first > 2)).sum()), int((det4 & (g4_first > 2)).sum()))

optics = [
    ("x at FSLT, rms (cm)", g3_x(Z_FSLT).std(), g4_x(Z_FSLT).std()),
    ("x at FSLT, mean (cm)", g3_x(Z_FSLT).mean(), g4_x(Z_FSLT).mean()),
    ("x at DSSSD, rms (cm)", g3_x(Z_DSSSD).std(), g4_x(Z_DSSSD).std()),
    ("x at DSSSD, mean (cm)", g3_x(Z_DSSSD).mean(), g4_x(Z_DSSSD).mean()),
    ("x' at DSSSD, mean (mrad)", g3_xp.mean(), g4_xp.mean()),
    ("x' at DSSSD, rms (mrad)", g3_xp.std(), g4_xp.std()),
    ("y at DSSSD, rms (cm)", g3_y(Z_DSSSD).std(), g4_y(Z_DSSSD).std()),
    ("y' at DSSSD, rms (mrad)", g3_yp.std(), g4_yp.std()),
    ("x waist position z (cm)", zs[np.argmin(rms3)], zs[np.argmin(rms4)]),
]
txt = [f"GEANT3: {N3all} events generated, {N3} reacting (denominator). GEANT3 global-frame sample: {len(A3)} recoils.",
       f"Ancalagon: {N4} events.", "", "== Efficiencies (%) =="]
txt += [f"  {r['label']:34s} GEANT3 {r['g3'][0]:6.2f} +- {r['g3'][1]:4.2f}   Ancalagon {r['g4'][0]:6.2f} +- {r['g4'][1]:4.2f}" for r in rows]
txt += [f"  BGO crystal multiplicity (hit events): GEANT3 {g3_n[g3_n>0].mean():.3f}  Ancalagon {g4_n[g4_n>0].mean():.3f}", "", "== Final-focus optics (global frame) =="]
txt += [f"  {a:28s} GEANT3 {b:8.3f}   Ancalagon {c:8.3f}" for a, b, c in optics]
txt += ["", "== Ancalagon x at FSLT vs initial recoil angle/momentum (3k events) ==",
        "  x = %.3f %+.4f*th[mrad] %+.3f*dp[%%] %+.5f*th^2 %+.4f*th*dp %+.3f*dp^2   residual rms %.3f cm" % (*coef, fit_resid),
        f"  single-ray 40Ca8+ fan slope at FSLT: {10*ca_slope:.2f} mm/mrad",
        "  40Ca8+ fan (mrad, x@FSLT cm): " + ", ".join(f"({a:+d}, {x:+.3f})" for a, x in fan_ca),
        "  19Ne4+ fan (mrad, x@FSLT cm): " + ", ".join(f"({a:+d}, {x:+.3f})" for a, x in fan_ne)]
open(f"{R}/summary.txt", "w").write("\n".join(txt) + "\n"); print("\n".join(txt))

# ---------------- chart data for the report page ----------------
def hist(a, bins): h, _ = np.histogram(a, bins); return (h / max(len(a), 1)).round(5).tolist()
xb = np.linspace(-3, 3, 61); eb = np.linspace(0, 10, 101)
json.dump(dict(
    rows=rows, optics=[dict(label=a, g3=round(float(b), 4), g4=round(float(c), 4)) for a, b, c in optics],
    waist=dict(dist=(Z_FSLT - zs).round(2).tolist(), g3=rms3.round(4).tolist(), g4=rms4.round(4).tolist()),
    xfslt=dict(edges=xb.round(3).tolist(), g3=hist(g3_x(Z_FSLT), xb), g4=hist(g4_x(Z_FSLT), xb)),
    fan=dict(ca=fan_ca, ne=fan_ne, slope_mm_per_mrad=round(float(10 * ca_slope), 3),
             fit_slope_mm_per_mrad=round(float(10 * coef[1]), 3), fit_resid_cm=round(float(fit_resid), 3)),
    bgo=dict(edges=eb.round(2).tolist(), g3=hist(g3_first[g3_first > 0.1], eb), g4=hist(g4_first[g4_first > 0.1], eb),
             g3_frac=round(float((g3_first > 0.1).mean()), 4), g4_frac=round(float((g4_first > 0.1).mean()), 4)),
    n=dict(g3_all=N3all, g3=N3, g3_global=len(A3), g4=N4, fit=len(xF)),
), open(f"{R}/chart_data.json", "w"))

# ---------------- static figure ----------------
C3, C4 = "#2a78d6", "#eb6834"
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
a = ax[0, 0]; dist = Z_FSLT - zs
a.plot(dist, rms3, color=C3, lw=2, label="GEANT3"); a.plot(dist, rms4, color=C4, lw=2, label="Ancalagon")
for d, n in ((0, "FSLT"), (Z_FSLT - Z_DSSSD, "DSSSD")): a.axvline(d, color="0.6", lw=1, ls=":"); a.text(d, a.get_ylim()[1] * 0.95, n, ha="center", fontsize=9)
a.set_xlabel("distance downstream of FSLT (cm)"); a.set_ylabel("x rms (cm)"); a.set_title("Horizontal beam size along the final drift"); a.legend()
a = ax[0, 1]
a.stairs(np.array(hist(g3_x(Z_FSLT), xb)), xb, color=C3, lw=2, label=f"GEANT3 (n={len(A3)})")
a.stairs(np.array(hist(g4_x(Z_FSLT), xb)), xb, color=C4, lw=2, label=f"Ancalagon (n={N4})")
a.set_xlabel("x at FSLT (cm, rel. to axis)"); a.set_ylabel("fraction per 0.1 cm"); a.set_title("Horizontal profile at FSLT"); a.legend()
a = ax[1, 0]; ne = np.array(fan_ne)
a.plot(ca[:, 0], ca[:, 1], "o-", color=C4, lw=2, label=f"40Ca8+ (k39pg), {10*ca_slope:.2f} mm/mrad")
a.plot(ne[:, 0], ne[:, 1], "s--", color=C4, lw=1.5, mfc="white", label="19Ne4+ (o15ag)")
a.axhline(0, color="0.6", lw=1); a.set_xlabel("initial angle (mrad)"); a.set_ylabel("x at FSLT (cm)"); a.set_title("Ancalagon single-ray angle fan"); a.legend()
a = ax[1, 1]
a.stairs(np.array(hist(g3_first[g3_first > 0.1], eb)), eb, color=C3, lw=1.5, label="GEANT3")
a.stairs(np.array(hist(g4_first[g4_first > 0.1], eb)), eb, color=C4, lw=1.5, label="Ancalagon")
a.set_yscale("log"); a.set_xlabel("highest-crystal energy (MeV)"); a.set_ylabel("fraction of hit events per 0.1 MeV"); a.set_title("BGO response"); a.legend()
fig.suptitle("39K(p,g)40Ca Er=606 keV, 2014 tune: Ancalagon vs GEANT3")
fig.tight_layout(); fig.savefig(f"{R}/compare.png", dpi=110)
