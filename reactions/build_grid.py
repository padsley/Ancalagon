#!/usr/bin/env python3
"""Drives generate_reactions.py to populate the dense grid: full 1g/2g
extensions (cheap, ungated) plus a Sobol-sampled 3g supplement and new
4g coverage, sized from Phase 3's measured per-run disk/time cost."""
import os
from scipy.stats import qmc

import generate_reactions as gr

REACTIONS_DIR = os.path.dirname(os.path.abspath(__file__))


def existing_calib1g_energies_mev():
    out = set()
    for i in range(1, 26):
        out.add(round(i * 0.1, 4))
    return out


def build_calib1g_dense():
    existing = existing_calib1g_energies_mev()
    targets = [round(0.05 + 0.05 * i, 4) for i in range(176) if round(0.05 + 0.05 * i, 4) <= 8.80 + 1e-9]
    new = [e for e in targets if round(e, 4) not in existing]
    n = 0
    for e in new:
        path = os.path.join(REACTIONS_DIR, f"k39pg_40ca_calib1g_d{gr._fmt_kev(e)}keV.reaction")
        gr.write_calib1g_file(path, e)
        n += 1
    return n, len(new)


def build_2g_dense(terminus, max_mev, prefix):
    # existing files sit at multiples of 0.1 in (0, max_mev]; add the
    # interstitial half-steps (0.05, 0.15, ..., max-0.05).
    targets = []
    e = 0.05
    while e < max_mev - 1e-9:
        targets.append(round(e, 4))
        e += 0.10
    n = 0
    for L in targets:
        path = os.path.join(REACTIONS_DIR, f"{prefix}_d{gr._fmt_kev(L)}keV.reaction")
        gr.write_reaction_file(path, [L], terminus, header_extra=f"2-gamma dense interstitial point, L1={L:.4f} MeV.")
        n += 1
    return n, len(targets)


def build_sobol_sample(terminus, e_budget, out_dir, n_levels, n_points, seed, tag):
    """Quasi-random (Sobol) sample of the ordered n_levels-simplex
    0 < L_1 < ... < L_{n_levels} < e_budget -- order-statistics-of-uniform
    method: draw n_levels independent uniforms per point, sort ascending."""
    os.makedirs(out_dir, exist_ok=True)
    m = max(1, (n_points - 1).bit_length())  # Sobol wants a power-of-2 sample size
    sampler = qmc.Sobol(d=n_levels, scramble=True, seed=seed)
    pts = sampler.random_base2(m=m)[:n_points] * e_budget
    n = 0
    for row in pts:
        levels = sorted(round(x, 4) for x in row)
        if len(set(levels)) != n_levels or levels[0] <= 0 or levels[-1] >= e_budget:
            continue
        stem = "_".join(f"L{i+1}-{gr._fmt_kev(L)}keV" for i, L in enumerate(levels))
        path = os.path.join(out_dir, f"{tag}_{stem}.reaction")
        gr.write_reaction_file(path, levels, terminus,
                                header_extra=f"{n_levels + 1}-gamma Sobol grid-supplement point.")
        n += 1
    return n


if __name__ == "__main__":
    n_calib, n_calib_target = build_calib1g_dense()
    print(f"calib1g dense: wrote {n_calib} new files (target {n_calib_target})")

    n_2g, n_2g_target = build_2g_dense("ground", 8.8, "k39pg_40ca_cascade")
    print(f"2g ground dense: wrote {n_2g} new files (target {n_2g_target})")

    n_2g0p2, n_2g0p2_target = build_2g_dense("0p2", 5.4, "k39pg_40ca_cascade0p2")
    print(f"2g 0+_2 dense: wrote {n_2g0p2} new files (target {n_2g0p2_target})")

    # Final campaign, sized from Phase 3's measured per-run cost (~9.5MB,
    # ~23s at 8000 events regardless of multiplicity -- cost tracks total
    # gamma energy carried, not gamma count) and padsley's chosen budget
    # (~5000 new runs total, ~47GB, ~2hr wall-clock parallelized across 16
    # cores). 3g/0+_2 split proportional to the existing 3828:1431 grid.
    n_3g = build_sobol_sample("ground", gr.E_BUDGET_GROUND,
                               os.path.join(REACTIONS_DIR, "k39pg_40ca_cascade3g"),
                               n_levels=2, n_points=2184, seed=10, tag="sobol")
    print(f"3g ground supplement: wrote {n_3g} files")

    n_3g0p2 = build_sobol_sample("0p2", gr.E_BUDGET_0P2,
                                  os.path.join(REACTIONS_DIR, "k39pg_40ca_cascade3g0p2"),
                                  n_levels=2, n_points=816, seed=11, tag="sobol")
    print(f"3g 0+_2 supplement: wrote {n_3g0p2} files")

    n_4g = build_sobol_sample("ground", gr.E_BUDGET_GROUND,
                               os.path.join(REACTIONS_DIR, "k39pg_40ca_cascade4g"),
                               n_levels=3, n_points=1237, seed=12, tag="sobol")
    print(f"4g ground: wrote {n_4g} files")

    n_4g0p2 = build_sobol_sample("0p2", gr.E_BUDGET_0P2,
                                  os.path.join(REACTIONS_DIR, "k39pg_40ca_cascade4g0p2"),
                                  n_levels=3, n_points=463, seed=13, tag="sobol")
    print(f"4g 0+_2: wrote {n_4g0p2} files")
