#!/usr/bin/env bash
# export.sh <runs_dir> <dest_dir>: the small files of every run (curve, damage
# history, status, convergence, provenance, log tail) - what a presentation
# needs and what fits in git.  The .frd fields stay in the container.
# Renamed to .txt because the repository's .gitignore drops *.sta, *.cvg,
# *.damage and *.log.
R=$(readlink -m "$1"); D=$(readlink -m "$2"); mkdir -p "$D"
python3 /home/user/fn-series/tools/curve.py "$R"/fn* > "$D/series.txt" 2>&1
cp -f "$R/summary.txt" "$D/summary.txt" 2>/dev/null
for r in "$R"/fn*/; do
  n=$(basename "$r"); mkdir -p "$D/$n"
  cp -f "$r/curve.csv" "$D/$n/" 2>/dev/null
  cp -f "$r/provenance.txt" "$D/$n/" 2>/dev/null
  cp -f "$r/m.sta" "$D/$n/sta.txt" 2>/dev/null
  cp -f "$r/m.cvg" "$D/$n/cvg.txt" 2>/dev/null
  cp -f "$r/m.damage" "$D/$n/damage.txt" 2>/dev/null
  tail -n 400 "$r/run.log" > "$D/$n/run_tail.txt" 2>/dev/null
  grep -a -E 'FRACTURE COMPLETE|DAMAGE COMMIT|DE1.2 COMMIT|HINGE CLUSTER' "$r/run.log" > "$D/$n/events.txt" 2>/dev/null
done
