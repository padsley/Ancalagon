"""Idealized check of GEANT3's E-dipole problem, independent of both codes:
rays through a hard-edge cylindrical electrostatic deflector with E2's
geometry (R = 250 cm, 35 deg, field E0*R/r) plus the RC49->EFB and
exit->RC51 drifts from the chain geometry, integrated two ways --
  exact : kinetic energy follows the potential (as Geant4 does)
  frozen: speed held fixed inside the field (what GEANT3's gthion.f does
          when each step's momentum change is below its 1 keV/c cut)
-- then the (x, a) transfer matrix and its determinant, as in section_matrix.py.
Non-relativistic 19Ne4+ at 1.79 MeV (beta = 0.0142)."""
import math
import numpy as np

R, PHI = 250.0, math.radians(35.0)
HERE = __file__.rsplit("/", 1)[0]
geo = {}
for l in open(f"{HERE}/../data/ancalagon/chain_geometry.txt"):
    p = l.split(); geo[p[1]] = (float(p[2].split("=")[1]), float(p[3].split("=")[1]))
dist = lambda a, b: math.hypot(geo[a][0] - geo[b][0], geo[a][1] - geo[b][1])
L1, L2 = dist("RC49", "E2"), dist("E2_exit", "RC51")
M, KE0, q = 17698.0981, 1.7924, 4.0            # MeV, MeV, e
E0 = 2 * KE0 / (q * R)                          # MV/cm so the design ray stays at r = R

def track(x0, a0, mode, h=0.05):
    # drift L1
    x = x0 + L1 * math.tan(a0 * 1e-3)
    # polar coords about the deflector centre; start at angle 0, radius R + x
    r, th = R + x, 0.0
    ke = KE0
    p = math.sqrt(ke * (ke + 2 * M))
    vr, vt = math.sin(a0 * 1e-3), math.cos(a0 * 1e-3)   # direction (radial, tangential), outward +
    # integrate in cartesian: centre at origin, start at (r, 0) moving +y (tangential)
    X, Y = r, 0.0; ux, uy = vr, vt
    while math.atan2(Y, X) < PHI:
        rr = math.hypot(X, Y); Er = E0 * R / rr        # radial field magnitude (points inward: force inward)
        fx, fy = -q * Er * X / rr, -q * Er * Y / rr     # force (MeV/cm), toward the centre
        pvec = np.array([ux, uy]) * p
        # momentum change over ds: dp = F * dt = F * ds / v ; v = p c^2 / E -> dp = F * ds * Etot / (p) (c=1)
        etot = ke + M
        dp = np.array([fx, fy]) * h * etot / p
        if mode == "frozen":
            pn = pvec + dp; pn = pn / np.linalg.norm(pn) * p                  # direction only
        else:
            pn = pvec + dp
        X += h * pvec[0] / p; Y += h * pvec[1] / p
        if mode == "exact":
            # kinetic energy from the potential V(r) = q*E0*R*ln(r/R)
            ke = KE0 - q * E0 * R * math.log(math.hypot(X, Y) / R)
            p = math.sqrt(ke * (ke + 2 * M))
        ux, uy = pn / np.linalg.norm(pn)
    # exit: transverse offset and angle relative to the exit axis (radial direction at theta=PHI)
    rr = math.hypot(X, Y); ang_pos = math.atan2(Y, X)
    er = np.array([math.cos(PHI), math.sin(PHI)]); et = np.array([-math.sin(PHI), math.cos(PHI)])
    x_out = rr - R
    a_out = math.atan2(np.dot([ux, uy], er), np.dot([ux, uy], et))
    x_out += L2 * math.tan(a_out)
    return x_out, a_out * 1e3

for mode in ("exact", "frozen"):
    xp, ap = track(0.2, 0, mode); xm, am = track(-0.2, 0, mode)
    xa1, aa1 = track(0, 2, mode); xa2, aa2 = track(0, -2, mode)
    xx = (xp - xm) / 0.4; ax = (ap - am) / 0.4; xa = (xa1 - xa2) / 4 * 10; aa = (aa1 - aa2) / 4
    print(f"{mode:6s}: (x|x)={xx:+.4f} (x|a)={xa:+.3f} mm/mrad (a|x)={ax:+.4f} mrad/cm (a|a)={aa:+.4f}  det={xx*aa-(xa/10)*ax:.4f}")
print("measured: GEANT3 (a|a)=+0.837 det=1.305 ; Ancalagon (a|a)=+0.537 det=0.998 (RC49->RC51, with fringes)")
