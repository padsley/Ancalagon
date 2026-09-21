#!/usr/bin/env python3
"""Generate k39(p,gamma)40Ca gamma-cascade .reaction files on an
energy-combination grid, for the Ancalagon_k39_emulator dense-grid
simulation campaign (see that repo's README, "Dense combinatorial
simulation grid" section, for why).

Every intermediate "level" this pilot's cascades pass through is a
FICTIONAL placeholder (LEVL cards carry unused 1e-12 s lifetimes; see
ReactionConfig.hh's own LEVL comment -- this pilot's cascades are
"instantaneous at vertex", not a real modeled flight time), so gamma
energies can be placed anywhere on a grid without violating any real
nuclear-structure constraint.

Unifying RECL formula (verified against every existing hand-written
reaction file in this directory before trusting it here): for a cascade
terminating at a state with energy `e_terminus_excitation` above true
ground, with a total energy budget `e_budget` = Ex - e_terminus_excitation
available to be shared across the cascade's gammas,

    RECL(e_budget) = BEAM_MASS_EXCESS + TARG_MASS_EXCESS - (e_budget - ERES)
                    = -25.91223 - e_budget   (this reaction, specifically)

True ground: e_budget = Ex = 8.93377 -> RECL = -34.8460 (matches every
existing k39pg_40ca_cascade*/cascade3g* file). Real 0+_2 state
(3.3526 MeV): e_budget = Ex - 3.3526 = 5.58117 -> RECL = -31.4934
(matches every existing *0p2* file). A calib1g-style single-gamma
source at an arbitrary chosen energy E: e_budget = E (matches every
existing calib1g file).

Cascade energy bookkeeping, N gammas via N-1 ordered levels
L_1 < L_2 < ... < L_{N-1} (all as excitation above the terminus, in
(0, e_budget)):
    gamma_top    = e_budget - L_{N-1}   (resonance -> highest level)
    gamma_i      = L_{i+1} - L_i        (each intermediate step)
    gamma_bottom = L_1                  (lowest level -> terminus)
(before recoil-kinematics corrections, applied downstream by
ReactionKinematics -- unchanged from every existing file).
"""
import argparse
import os

BEAM_MASS_EXCESS = -33.8072
TARG_MASS_EXCESS = 7.28897
ERES = 0.606
EX_TOTAL = 8.93377  # Q (8.32777) + Er (0.606)
E_0P2 = 3.3526  # real 40Ca 0+_2 state excitation above true ground
E_BUDGET_GROUND = EX_TOTAL
E_BUDGET_0P2 = EX_TOTAL - E_0P2  # 5.58117

RECL_BASE = -25.91223  # = BEAM_MASS_EXCESS + TARG_MASS_EXCESS + ERES


def recl_for_budget(e_budget: float) -> float:
    return RECL_BASE - e_budget


def _fmt_kev(e_mev: float) -> str:
    return f"{round(e_mev * 1000):04d}"


def generate_reaction_text(
    level_energies_mev: list[float],
    e_budget: float,
    terminus_label: str,
    recl_massexcess: float,
    header_extra: str,
) -> str:
    """`level_energies_mev` ascending, in (0, e_budget); empty for a
    single-gamma (calib1g-style) source."""
    n_levels = len(level_energies_mev)
    n_gammas = n_levels + 1

    gammas = []
    prev = 0.0
    for L in level_energies_mev:
        gammas.append(L - prev)
        prev = L
    gammas.append(e_budget - prev)
    gammas = gammas[::-1]  # top-down emission order (highest first)

    lines = []
    lines.append(f"COMM 39K(p,gamma)40Ca {n_gammas}-gamma cascade, script-generated dense-grid")
    lines.append(f"COMM point (generate_reactions.py) -- terminus: {terminus_label}, e_budget =")
    lines.append(f"COMM {e_budget:.5f} MeV. Every level is a FICTIONAL placeholder, unused by")
    lines.append("COMM kinematics beyond fixing gamma energies (see ReactionConfig.hh's own LEVL")
    lines.append("COMM comment) -- same convention as every hand-written file in this directory.")
    if header_extra:
        lines.append(f"COMM {header_extra}")
    lvl_str = ", ".join(f"L{i+1}={e:.4f}" for i, e in enumerate(level_energies_mev))
    lines.append(f"COMM Levels (excitation above terminus, MeV): {lvl_str or '(none -- single gamma)'}")
    g_str = ", ".join(f"gamma{i+1}={g:.4f}" for i, g in enumerate(gammas))
    lines.append(f"COMM Gamma energies before recoil-kinematics corrections (MeV): {g_str}")
    lines.append("COMM Verify with `--reaction-stats <this file>`.")
    lines.append("")
    lines.append(f"BEAM 19  39   {BEAM_MASS_EXCESS:<10}   # Z A massExcessMeV -- 39K")
    lines.append(f"TARG 1   1    {TARG_MASS_EXCESS:<10}   # 1H (target, at rest)")
    lines.append(f"RECL 20  40   {recl_massexcess:<10}  8  # 40Ca -- {terminus_label}")
    lines.append(f"ERES {ERES}                  # CM resonance energy above threshold, MeV")
    lines.append("")
    for i, L in enumerate(level_energies_mev, start=1):
        lines.append(f"LEVL   {i}   {L:.4f}   1.0E-12   # fictional level {i}")
    lines.append("")
    if n_levels == 0:
        lines.append("BRAT   -1  100.0   0        # resonance -> terminus directly, 100%, single gamma")
    else:
        lines.append(f"BRAT   -1  100.0   {n_levels}        # resonance -> highest level, 100%")
        for i in range(n_levels, 1, -1):
            lines.append(f"BRAT    {i}  100.0   {i - 1}        # level {i} -> level {i - 1}, 100%")
        lines.append(f"BRAT    1  100.0   0        # level 1 -> terminus, 100%")
    lines.append("SENT")
    return "\n".join(lines) + "\n"


def write_reaction_file(out_path, level_energies_mev, terminus, header_extra=""):
    if terminus == "ground":
        e_budget = E_BUDGET_GROUND
        recl = -34.8460
        terminus_label = "true ground state"
    elif terminus == "0p2":
        e_budget = E_BUDGET_0P2
        recl = -31.4934
        terminus_label = "real 0+_2 state (3.3526 MeV, E0-conversion terminus)"
    else:
        raise ValueError(f"unknown terminus {terminus!r}")
    text = generate_reaction_text(level_energies_mev, e_budget, terminus_label, recl, header_extra)
    with open(out_path, "w") as f:
        f.write(text)


def write_calib1g_file(out_path, e_gamma_mev):
    e_budget = e_gamma_mev
    recl = recl_for_budget(e_budget)
    text = generate_reaction_text(
        [], e_budget,
        f"fictional single-gamma terminus at {e_gamma_mev:.4f} MeV excitation",
        recl,
        f"Single-gamma calibration source, RECL(e_budget={e_budget:.5f})={recl:.4f}.",
    )
    with open(out_path, "w") as f:
        f.write(text)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reactions-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    print(f"RECL(Ex_total={EX_TOTAL})={recl_for_budget(EX_TOTAL):.4f} (expect -34.8460)")
    print(f"RECL(Ex-E0p2={E_BUDGET_0P2:.5f})={recl_for_budget(E_BUDGET_0P2):.4f} (expect -31.4934)")
