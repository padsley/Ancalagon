#!/usr/bin/env bash
# GEANT3 side of the 39K(p,g)40Ca cross-check: 8 parallel dsbatch jobs of
# 1250 events each (distinct RNDM seeds), plus one 400-event job that keeps
# the recoil's global coordinates near the DSSSD (FOCUSTEST lines).
#
# Usage: run_g3.sh <workdir> [G3_DRAGON checkout, default ~/codes/G3_DRAGON_Claude]
# Outputs land in <workdir>/j1..j8/dragon1.root and <workdir>/jg/{endv,tail}.txt;
# copy them into data/g3 and data/g3_global to rerun analyze.py.
#
# Gotchas this works around:
#  - dsbatch truncates FFCARD/INPUT/MITRAY paths at 80 chars, so each job
#    runs from its own directory with relative paths.
#  - src/gustep_mitray.f's FOCUSTEST print fires on every recoil step, not
#    only in single-ray runs (~4 MB of log per 20 s), so it is filtered out.
#  - ~1.15 events/s per process.
set -euo pipefail
W=$(realpath -m "$1"); G=$(realpath "${2:-$HOME/codes/G3_DRAGON_Claude}")
HERE=$(cd "$(dirname "$0")" && pwd)
BASE=$HERE/data/g3/dragon_k39pg_10k.ffcards   # dragon_k39pg.ffcards with TRIG 10000

setup() {  # setup <dir> <trig> <seed1> <seed2>
  mkdir -p "$1"; cd "$1"
  for f in "$G"/*.dat "$G"/*.txt "$G"/*.root; do ln -sf "$f" .; done
  ln -sfn "$G/dat" dat
  sed "s/^TRIG .*/TRIG $2/; s/^RNDM .*/RNDM $3 $4/" "$BASE" > ff.cards
}
export DSROOT=$G DSDAT=$G/dat DSBIN=$G/bin MITRAY=dat/dragon_2014_DSSSD.dat FFCARD=ff.cards INPUT=dat/k39pg.dat

for i in 1 2 3 4 5 6 7 8; do
  ( setup "$W/j$i" 1250 $((1000*i+7)) $((31*i+11))
    "$G/bin/dsbatch" 2>&1 | grep -v -E "^FOCUSTEST|Whats stopping|istop:" > g3.log
    h2root dragon1.hbook dragon1.root > h2r.log 2>&1 ) &
done

( setup "$W/jg" 400 99007 7771
  "$G/bin/dsbatch" 2>&1 | awk '
    /GTRIGI/ {print "EVT" > "tail.txt"; prevvol=""}
    /^FOCUSTEST/ {
      if ($2=="ENDV" && prevvol!="ENDV") {print "PRE", prev > "endv.txt"; print "IN", $0 > "endv.txt"}
      if ($5 < -5.0 && $5 > -72.49) print > "tail.txt"
      prevvol=$2; prev=$0 }' ) &
wait
echo "done: $W"
