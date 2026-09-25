#!/usr/bin/env bash
# Ancalagon side of the 39K(p,g)40Ca cross-check, real 2014 tune.
# Usage: run_ancalagon.sh <workdir>   (writes straight into this directory's data/)
#
#  1. 10k --track-reaction events -> dragon_hits.root (BGO + DSSSD hit trees)
#  2. the same 10k events (same RANDOM_SEED) with the TRAJ dump on, keeping
#     only the last step before each DSSSD hit (position + momentum there)
#  3. 3k events (seed 777) keeping each recoil's first step and its last
#     step before the DSSSD, for the (x|theta) fit at FSLT
#  4. single-ray angle fans (needs --track-chain firing the reaction's own
#     recoil species): 40Ca8+ at 1320 MeV/c and 19Ne4+ at 258 MeV/c, both
#     from the target centre (z0 = 0)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
W=$(realpath -m "$1"); mkdir -p "$W"; cd "$W"
A=$REPO/build/Ancalagon; OUT=$HERE/data
K39=$REPO/reactions/k39pg_40ca.reaction; O15=$REPO/reactions/o15ag_19ne.reaction
export MSLT_HALFGAP_X_CM=0.75 MSLT_HALFGAP_Y_CM=1.25

REACTION_INPUT=$K39 TRAJ_QUIET=1 RANDOM_SEED=12345 "$A" --track-reaction 10000 > g4.log 2>&1
cp dragon_hits.root "$OUT/ancalagon/"

REACTION_INPUT=$K39 RANDOM_SEED=12345 "$A" --track-reaction 10000 2>&1 | awk '
  /^TRAJ/ { if ($8=="DSSSD" && prev!="") print "PRE", prev; prev=$0; next }
  /^HITS DSSSD/ {print}' > "$OUT/ancalagon/dsssd_traj.txt"

REACTION_INPUT=$K39 RANDOM_SEED=777 "$A" --track-reaction 3000 2>&1 | awk '
  /^TRAJ/ && /m=  37220/ { if (!seen) {v=$0; seen=1}
                           if ($8=="DSSSD" && prev!="") print "V", v, "\nP", prev; prev=$0 }
  /^HITS DSSSD/ {seen=0; prev=""}' > "$OUT/ancalagon/vp.txt"

deg() { python3 -c "print($1*0.001*57.29578)"; }
for a in -4 -2 -1 0 1 2 4; do
  REACTION_INPUT=$K39 TRACK_CHAIN_ANGLE_DEG=$(deg $a) "$A" --track-chain 0.0 0.0 1320.0 1 0.0 2>&1 \
    | grep "^TRAJ" > "$OUT/fan/ca.$a.txt"
done
unset MSLT_HALFGAP_X_CM MSLT_HALFGAP_Y_CM
for a in -2 -1 0 1 2; do
  REACTION_INPUT=$O15 TRACK_CHAIN_ANGLE_DEG=$(deg $a) "$A" --track-chain 0.0 0.0 258.0 1 0.0 2>&1 \
    | grep "^TRAJ" > "$OUT/fan/ne.$a.txt"
done
echo "done"
