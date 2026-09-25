#!/usr/bin/env bash
# Section test: identical rays launched AT a field-free plane in both codes
# -- GEANT3's own on-axis state there plus offsets in u (+-0.2 cm) and in
# angle (+-2 mrad) -- to measure one element's transfer matrix in each code
# (section_matrix.py). Usage:
#   run_section_plane.sh <workdir> <tag> <x_cm> <z_cm> <theta_deg> <u0_cm> <a0_mrad> [G3]
# e.g. E2: ... e2 -964.045378 510.33554 -145 0.3190 -0.6667   (RC49)
#      E1: ... e1 -320.428625 570.591405 -50 0.0221 -0.4393   (RC25)
set -euo pipefail
W=$(realpath -m "$1"); TAG=$2; PX=$3; PZ=$4; PT=$5; U0=$6; A0=$7
G=$(realpath "${8:-$HOME/codes/G3_DRAGON_Claude}")
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../../.." && pwd)
mkdir -p "$HERE/g3_$TAG" "$HERE/ancalagon_$TAG" "$W/g3"
RAYS="c:0:0 up:0.2:0 um:-0.2:0 ap:0:2 am:0:-2"
state() {  # name du da -> x0 z0 alpha_mrad alpha_deg
  python3 - "$1" "$2" "$PX" "$PZ" "$PT" "$U0" "$A0" <<'PY'
import math, sys
du, da, xm, zm, deg, u0, a0 = map(float, sys.argv[1:8])
r = math.radians(deg)
u = u0 + du; a = a0 + da
x0 = xm + u * math.cos(r); z0 = zm - u * math.sin(r); al = r * 1000 + a
print(f"{x0:.6f} {z0:.6f} {al:.4f} {al * 1e-3 * 180 / math.pi:.8f}")
PY
}
cd "$W/g3"
for f in "$G"/*.dat "$G"/*.txt "$G"/*.root; do ln -sf "$f" .; done; ln -sfn "$G/dat" dat; cp "$HERE/../o15ag.dat" .
export DSROOT=$G DSDAT=$G/dat DSBIN=$G/bin MITRAY=dat/dragon_2014_DSSSD.dat FFCARD=ff.cards INPUT=o15ag.dat
for r in $RAYS; do
  IFS=: read name du da <<< "$r"; read x0 z0 al deg < <(state $du $da)
  sed -e "s/^FKIN .*/FKIN 20 4.0 4.0/" -e "s/^KINE .*/KINE 85 $x0 0.0 $z0 251.89 $al 0.0 0.0 0.0 0.0/" "$G/dragon_focustest.ffcards" > ff.cards
  "$G/bin/dsbatch" > "$TAG.$name.log" 2>&1; grep "^FOCUSTEST" "$TAG.$name.log" > "$HERE/g3_$TAG/ray.$name.txt"
  FINE_STEP_CM=0.5 REACTION_INPUT=$REPO/reactions/o15ag_19ne.reaction TRACK_CHAIN_ANGLE_DEG=$deg \
    "$REPO/build/Ancalagon" --track-chain "$x0" 0.0 249.43 1 "$z0" 2>&1 | grep "^TRAJ" > "$HERE/ancalagon_$TAG/ray.$name.txt"
done
echo done
