#!/usr/bin/env bash
# saver.sh <runs_dir> <repo_dir> <subdir> [period_s]: every period, export the
# runs' small files into <repo_dir>/<subdir>, commit if anything changed and
# push the current branch; stops after the series has finished and the last
# state is pushed.  Push retries 4 times with exponential backoff.
R=$1; REPO=$2; SUB=$3; P=${4:-900}
BR=$(git -C "$REPO" rev-parse --abbrev-ref HEAD)
push() {
  for d in 0 2 4 8 16; do
    sleep $d
    git -C "$REPO" push -q -u origin "$BR" && return 0
  done
  echo "$(date -u +%FT%TZ) push FAILED" ; return 1
}
while true; do
  done_flag=0; grep -q 'SERIES DONE' "$R/summary.txt" 2>/dev/null && done_flag=1
  bash /home/user/fn-series/tools/export.sh "$R" "$REPO/$SUB"
  git -C "$REPO" add -A "$SUB"
  if ! git -C "$REPO" diff --cached --quiet; then
    s=$(cat "$REPO/$SUB/summary.txt" 2>/dev/null | grep -c rc=)
    git -C "$REPO" commit -q -m "Fn series v1: results snapshot ($s of 8 runs finished)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Rw7uNBF1wKvvyKCt49DLBF" && push && echo "$(date -u +%FT%TZ) pushed ($s finished)"
  fi
  [ $done_flag = 1 ] && { echo "$(date -u +%FT%TZ) series done, saver exits"; exit 0; }
  sleep "$P"
done
