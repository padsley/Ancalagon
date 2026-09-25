#!/usr/bin/env bash
# GEANT3 recoil tracks near the target, MSLT and FSLT (global FOCUSTEST
# coordinates), for a population (x|theta) fit comparable with Ancalagon's.
# GEANT3 only prints from Q1 onward (z ~ 88 cm), and takes long straight
# steps through the field-free drifts, hence the wide capture radii.
# Usage: run_g3_planes.sh <workdir> [G3 checkout]; 4 jobs x 400 events.
set -euo pipefail
W=$(realpath -m "$1"); G=$(realpath "${2:-$HOME/codes/G3_DRAGON_Claude}")
HERE=$(cd "$(dirname "$0")" && pwd); BASE=$HERE/data/g3/dragon_o15ag_10k.ffcards
export DSROOT=$G DSDAT=$G/dat DSBIN=$G/bin MITRAY=dat/dragon_2014_DSSSD.dat FFCARD=ff.cards INPUT=o15ag.dat
for i in 1 2 3 4; do
  ( d=$W/p$i; mkdir -p "$d"; cd "$d"
    for f in "$G"/*.dat "$G"/*.txt "$G"/*.root; do ln -sf "$f" .; done; ln -sfn "$G/dat" dat; cp "$HERE/o15ag.dat" .
    sed "s/^TRIG .*/TRIG 400/; s/^RNDM .*/RNDM $((5000+37*i)) $((77*i+3))/" "$BASE" > ff.cards
    "$G/bin/dsbatch" 2>&1 | awk '
      /GTRIGI/ {print "EVT"}
      /^FOCUSTEST/ { x=$3; z=$5
        if (z < 92 && x > -5 && x < 5) print
        else if ((x+498.206)^2 + (z-657.134)^2 < 160000) print
        else if ((x+1024.717)^2 + (z+6.075)^2 < 40000) print }' > planes.txt ) &
done
wait; cat "$W"/p*/planes.txt > "$HERE/data/g3_global/planes.txt"; echo done
