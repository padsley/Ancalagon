#!/usr/bin/env bash
# Matched single-ray 19Ne4+ fans in both codes, for an element-by-element
# comparison between Q1 and MSLT (see compare.py). Each ray starts at
# z0 = 90 cm (inside Q1, past the target gas, where GEANT3's tracking
# starts), x0 = 90 cm * tan(theta), i.e. as if from the target centre,
# at the code's own tuned rigidity (GEANT3 0.9736, Ancalagon RTUN 0.9641,
# times 258.72 MeV/c). Ancalagon's slits are opened so no ray stops early.
# Usage: run_fans.sh <workdir> [G3 checkout]
set -euo pipefail
W=$(realpath -m "$1"); G=$(realpath "${2:-$HOME/codes/G3_DRAGON_Claude}")
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../../.." && pwd)
ANGLES="-8 -4 -2 0 2 4 8"
P3=251.89; P4=249.43; Z0=90.0

mkdir -p "$W/g3"; cd "$W/g3"
for f in "$G"/*.dat "$G"/*.txt "$G"/*.root; do ln -sf "$f" .; done; ln -sfn "$G/dat" dat; cp "$HERE/../o15ag.dat" .
export DSROOT=$G DSDAT=$G/dat DSBIN=$G/bin MITRAY=dat/dragon_2014_DSSSD.dat FFCARD=ff.cards INPUT=o15ag.dat
for a in $ANGLES; do
  x0=$(python3 -c "import math; print(f'{$Z0*math.tan($a*1e-3):.6f}')")
  af=$(python3 -c "print(float($a))")   # FFREAD misreads a bare integer in a REAL card
  sed -e "s/^FKIN .*/FKIN 20 4.0 4.0/" -e "s/^KINE .*/KINE 85 $x0 0.0 $Z0 $P3 $af 0.0 0.0 0.0 0.0/" \
      "$G/dragon_focustest.ffcards" > ff.cards
  "$G/bin/dsbatch" > "run.$a.log" 2>&1
  grep "^FOCUSTEST" "run.$a.log" > "$HERE/g3/ray.$a.txt"
  grep -E "Magnetic element scale|Electric element scale" "run.$a.log" | head -2 > "$HERE/g3/scale.txt"
done

mkdir -p "$W/g4"; cd "$W/g4"
for a in $ANGLES; do
  x0=$(python3 -c "import math; print(f'{$Z0*math.tan($a*1e-3):.6f}')")
  deg=$(python3 -c "print($a*1e-3*180/3.141592653589793)")
  REACTION_INPUT=$REPO/reactions/o15ag_19ne.reaction TRACK_CHAIN_ANGLE_DEG=$deg \
  QSLT_HALFGAP_X_CM=20 QSLT_HALFGAP_Y_CM=20 MSLT_HALFGAP_X_CM=20 MSLT_HALFGAP_Y_CM=20 \
    "$REPO/build/Ancalagon" --track-chain "$x0" 0.0 $P4 1 $Z0 2>&1 | grep "^TRAJ" > "$HERE/ancalagon/ray.$a.txt"
done
echo done
