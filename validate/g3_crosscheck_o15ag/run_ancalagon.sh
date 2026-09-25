#!/usr/bin/env bash
# Ancalagon side of the 15O(a,g)19Ne cross-check, real 2014 tune.
# Usage: run_ancalagon.sh <workdir>   (writes straight into this directory's data/)
#
# One 10k-event --track-reaction run with the TRAJ dump on, reduced per
# event to: the recoil's last volume and position (where it stopped, or
# DSSSD), and its last step before the DSSSD (position + momentum).
# The same run writes dragon_hits.root (BGO + DSSSD hit trees).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
W=$(realpath -m "$1"); mkdir -p "$W"; cd "$W"
export REACTION_INPUT=$REPO/reactions/o15ag_19ne.reaction MSLT_HALFGAP_X_CM=0.75 MSLT_HALFGAP_Y_CM=1.25 RANDOM_SEED=12345
"$REPO/build/Ancalagon" --track-reaction 10000 2>&1 | awk '
  function flush() { if (have) print "END", last; have=0 }
  /^TRAJ/ && /m=  17698/ { if ($8=="DSSSD" && prev!="") print "PRE", prev; prev=$0; last=$0; have=1; next }
  /^HITS BGO/ { flush(); prev="" }' > "$HERE/data/ancalagon/recoil_fate.txt"
cp dragon_hits.root "$HERE/data/ancalagon/"
echo done

# Recoil tracks near the target, MSLT and FSLT (same regions as
# run_g3_planes.sh), 2000 events, for the population (x|theta) fit.
RANDOM_SEED=2024 "$REPO/build/Ancalagon" --track-reaction 2000 2>&1 | awk '
  /^HITS BGO/ {print "EVT"}
  /^TRAJ/ && /m=  17698/ { z=$2; x=$3
    if (z < 60 && x > -3 && x < 3) print "FOCUSTEST", $8, $3, $4, $2
    else if ((x+498.206)^2 + (z-657.134)^2 < 10000) print "FOCUSTEST", $8, $3, $4, $2
    else if ((x+1024.717)^2 + (z+6.075)^2 < 10000) print "FOCUSTEST", $8, $3, $4, $2 }' > "$HERE/data/ancalagon/planes.txt"
echo done planes

# Scatter-free single-ray 19Ne4+ angle fan (see fan.py): fired from
# z0 = 74 cm, past the dense target gas, slits opened wide.
unset RANDOM_SEED
for a in -16 -12 -8 -4 -2 -1 0 1 2 4 8 12 16; do
  deg=$(python3 -c "print($a*0.001*57.29578)"); x0=$(python3 -c "import math; print(74*math.tan($a*0.001))")
  QSLT_HALFGAP_X_CM=20 QSLT_HALFGAP_Y_CM=20 MSLT_HALFGAP_X_CM=20 MSLT_HALFGAP_Y_CM=20 TRACK_CHAIN_ANGLE_DEG=$deg \
    "$REPO/build/Ancalagon" --track-chain "$x0" 0.0 248.9 1 74.0 2>&1 | grep "^TRAJ" > "$HERE/data/fan/f.$a.txt"
done
echo done fan
