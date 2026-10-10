#!/usr/bin/env bash
# run_series.sh <out_dir> [deck ...]: the Fn series, JOBS runs at a time
# (default 4, one thread each), most informative decks first.
OUT=$(readlink -m "${1:-/home/user/fn-series/runs_v1}"); shift
DECKS=${*:-"fn0 fn100 fn25 fn40 fn20 fn33 fn60 fn80"}
mkdir -p "$OUT"
for d in $DECKS; do echo "$d"; done |
  xargs -P "${JOBS:-4}" -I{} bash /home/user/fn-series/tools/run_one.sh \
    /home/user/fn-series/v1/{}.inp "$OUT/{}"
echo "SERIES DONE" >> "$OUT/summary.txt"
