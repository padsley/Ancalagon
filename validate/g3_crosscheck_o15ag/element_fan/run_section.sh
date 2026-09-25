#!/usr/bin/env bash
# Section test: identical rays launched AT the QSLT plane in both codes,
# using GEANT3's own on-axis state there (u = 0.0050 - 0.00039*theta cm,
# a = 0.4763 - 2.2720*theta mrad, theta = the original angle at the
# target), so D1's different central-ray bend is taken out and only
# Q3..E1..MSLT is compared. Usage: run_section.sh <workdir> [G3 checkout]
set -euo pipefail
W=$(realpath -m "$1"); G=$(realpath "${2:-$HOME/codes/G3_DRAGON_Claude}")
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../../.." && pwd)
mkdir -p "$HERE/g3_qslt" "$HERE/ancalagon_qslt" "$W/g3"
state() {  # prints x0 z0 alpha_mrad alpha_deg for original angle $1
  python3 - "$1" <<'PY'
import math, sys
t = float(sys.argv[1]); xm, zm, r = -59.376, 351.542243, math.radians(-50.0)
u = 0.0050 - 0.00039 * t; a = 0.4763 - 2.2720 * t
x0 = xm + u * math.cos(r); z0 = zm - u * math.sin(r)
al = r * 1000 + a
print(f"{x0:.6f} {z0:.6f} {al:.4f} {al * 1e-3 * 180 / math.pi:.8f}")
PY
}
cd "$W/g3"
for f in "$G"/*.dat "$G"/*.txt "$G"/*.root; do ln -sf "$f" .; done; ln -sfn "$G/dat" dat; cp "$HERE/../o15ag.dat" .
export DSROOT=$G DSDAT=$G/dat DSBIN=$G/bin MITRAY=dat/dragon_2014_DSSSD.dat FFCARD=ff.cards INPUT=o15ag.dat
for t in -8 -4 -2 0 2 4 8; do
  read x0 z0 al deg < <(state $t)
  sed -e "s/^FKIN .*/FKIN 20 4.0 4.0/" -e "s/^KINE .*/KINE 85 $x0 0.0 $z0 251.89 $al 0.0 0.0 0.0 0.0/" "$G/dragon_focustest.ffcards" > ff.cards
  "$G/bin/dsbatch" > "sec.$t.log" 2>&1; grep "^FOCUSTEST" "sec.$t.log" > "$HERE/g3_qslt/ray.$t.txt"
  REACTION_INPUT=$REPO/reactions/o15ag_19ne.reaction TRACK_CHAIN_ANGLE_DEG=$deg \
  QSLT_HALFGAP_X_CM=20 QSLT_HALFGAP_Y_CM=20 MSLT_HALFGAP_X_CM=20 MSLT_HALFGAP_Y_CM=20 \
    "$REPO/build/Ancalagon" --track-chain "$x0" 0.0 249.43 1 "$z0" 2>&1 | grep "^TRAJ" > "$HERE/ancalagon_qslt/ray.$t.txt"
done
echo done
