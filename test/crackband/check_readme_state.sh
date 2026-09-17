#!/usr/bin/env bash
# Do this directory's README claims about the TREE still hold?
#
#   test/crackband/check_readme_state.sh
#
# WHY THIS EXISTS, and it is not tidiness.  In one day this directory's
# documentation was caught four times asserting a state that had moved
# underneath it:
#
#   * a printed column labelled "theta=" that was the attempt counter, and the
#     label travelled into the other agent's notes before anyone read the
#     column header;
#   * a script's closing paragraph saying an experiment "is not established
#     here" and naming what would settle it, printed on the one path a reader
#     reaches, after that experiment had been run;
#   * two predictions printed against the wrong side of a switch, one of them
#     concluding "the lengths are not reconciled" under the switch that exists
#     to reconcile them;
#   * the README telling readers the width unit test was "a live proposal and
#     not mine to land" while it was a gate preflight.
#
# All four were claims about STATE, and none was a number.  That is the
# pattern: a number is reprinted by every run and so cannot go stale, while a
# claim about state is reprinted by nobody.  Each of the four was caught by
# luck - by running a tool rather than reading it, which is not a method.
#
# So the mechanically checkable claims are checked here, and this script goes
# red when one drifts.  It deliberately does NOT try to check prose: a claim
# about what a measurement means is not machine-checkable, and pretending
# otherwise would put a green tick on the half that actually goes wrong.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
R=$ROOT/test/crackband/README.md
bad=0
say(){ printf "  %-58s %s\n" "$1" "$2"; }
fail(){ say "$1" "FAIL - $2"; bad=$((bad+1)); }

# 1. the file count in the opening line
n=$(ls "$ROOT/test/crackband" | grep -v -e '^README\.md$' -e '__pycache__' | wc -l)
words="zero one two three four five six seven eight nine ten eleven twelve
thirteen fourteen fifteen sixteen seventeen eighteen nineteen twenty"
w=$(echo $words | cut -d' ' -f$((n+1)))
if head -1 "$R" | grep -qi "^# " && grep -qiE "^${w} files," "$R"; then
  say "opening line counts $n files ($w)" "ok"
else
  fail "opening line counts $n files ($w)" "$(sed -n '3p' "$R" | cut -c1-30)..."
fi

# 2. every file in the directory appears in the Files table
miss=0
for f in $(ls "$ROOT/test/crackband" | grep -v -e '^README\.md$' -e '__pycache__'); do
  if grep -q "| \`$f\`" "$R"; then :; else
    fail "$f listed in the Files table" "missing"; miss=$((miss+1)); fi
done
# Only claim this once nothing was missing.  The first version printed it
# unconditionally, right under its own FAIL lines - a green line next to the
# red one it contradicts, which is the shape of the defect this script exists
# to catch, reproduced inside the script on its first run.
[ "$miss" -eq 0 ] && say "every file appears in the Files table" "ok"

# 3. the claim that the gate runs the width unit test
if grep -q '\*\*The gate runs it\*\*' "$R"; then
  if grep -q 'test/crackband/run_cbwidth_test.sh' "$ROOT/test/regress/run.py"; then
    say "README says the gate runs the unit test, and it does" "ok"
  else
    fail "README says the gate runs the unit test" "run.py does not"
  fi
else
  if grep -q 'test/crackband/run_cbwidth_test.sh' "$ROOT/test/regress/run.py"; then
    fail "README no longer claims the gate runs it" "but run.py does"
  else
    say "README does not claim it, and the gate does not run it" "ok"
  fi
fi

# 4. every commit a measurement is stamped with still exists
for h in $(grep -o 'Measured on `[0-9a-f]\{7,40\}`' "$R" | tr -d '`' | awk '{print $3}'); do
  if git -C "$ROOT" cat-file -e "$h^{commit}" 2>/dev/null; then
    say "measurement stamp $h is a real commit" "ok"
  else
    fail "measurement stamp $h is a real commit" "no such commit"
  fi
done

echo
if [ "$bad" -eq 0 ]; then
  echo "  check_readme_state: all state claims still hold"
  exit 0
fi
echo "  check_readme_state: $bad state claim(s) have drifted"
exit 1
