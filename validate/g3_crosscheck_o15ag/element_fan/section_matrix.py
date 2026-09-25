"""Transfer matrix of one element in each code, from run_section_plane.sh
rays (central, +-0.2 cm in u, +-2 mrad in angle), at a downstream plane.
Angles come from positions at two field-free planes (Ancalagon prints z to
only 1e-4 cm, so a single 0.5 cm step's direction can be off by ~0.2 mrad
on diagonal lines). Also prints det = (x|x)(a|a) - (x|a)(a|x), which
phase-space conservation requires to be ~1 (p_in/p_out; the energy
changes here are < 0.5 %).
Usage: python3 section_matrix.py <tag> <plane1> <plane2>   (both field-free, same axis)"""
import sys
exec(open(__file__.replace("section_matrix.py", "compare.py")).read().split("angles = sorted")[0])
tag = sys.argv[1]
def ld(code, n):
    pts = []
    for l in open(f"{HERE}/{code}_{tag}/ray.{n}.txt"):
        p = l.split()
        pts.append((float(p[2]), float(p[4]), p[1]) if code == "g3" else (float(p[2]), float(p[1]), p[7]))
    return pts
p1, p2 = sys.argv[2], sys.argv[3]
(x1, z1, _), (x2, z2, _) = GEO[p1], GEO[p2]
L = math.hypot(x2 - x1, z2 - z1)   # cm along the axis
print(f"== {p1} (x) and {p2}, {L:.2f} cm apart (angle = du/L)")
for code in ("g3", "ancalagon"):
    u1 = {n: crossing(ld(code, n), p1)[0] for n in ("c", "up", "um", "ap", "am")}
    u2 = {n: crossing(ld(code, n), p2)[0] for n in ("c", "up", "um", "ap", "am")}
    a = {n: 1000 * (u2[n] - u1[n]) / L for n in u1}                                  # mrad
    xx = (u1["up"] - u1["um"]) / 0.4; ax = (a["up"] - a["um"]) / 0.4
    xa = (u1["ap"] - u1["am"]) / 4 * 10; aa = (a["ap"] - a["am"]) / 4
    det = xx * aa - (xa / 10) * ax
    print(f"  {code:9s} central u={u1['c']:+.4f} cm a={a['c']:+.4f} mrad | (x|x)={xx:+.4f} (x|a)={xa:+.3f} mm/mrad "
          f"(a|x)={ax:+.4f} mrad/cm (a|a)={aa:+.4f} | det={det:.4f}")
