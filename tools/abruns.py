#!/usr/bin/env python3
"""Compare two whole gate runs byte for byte, and say what is allowed to differ.

    tools/abruns.py A_RUNDIR B_RUNDIR

WHY THIS EXISTS
---------------
The gate already proves a lot: sixteen cases, their scalars against the
record, and `equal_to` between cases WITHIN one run.  What it could not say
is the thing an extraction actually claims - that this binary produces the
same bytes as the one before it.  That comparison was being done by hand,
and a comparison done by hand is one that gets skipped on the day it
matters.

WHAT IS ALLOWED TO DIFFER, AND WHY
----------------------------------
Exactly two things, both of them clocks rather than arithmetic:

  * the UDATE and UTIME records of a .frd file, which are the calendar date
    and wall-clock time the run started.  frdheader.c writes them from
    time(); two runs of the identical binary differ there too, so treating
    them as differences would make the check cry wolf on every pair.
    UTIME was stripped from the first version and UDATE was not, and the
    omission showed itself the way these always do: a comparison run either
    side of midnight reported sixteen differing .frd files, one per case,
    all of them identical in length.  A check that goes red once a night is
    a check people learn to ignore;
  * provenance.txt, which records the binary's sha256 on purpose - it is
    SUPPOSED to change when the binary changes.

Everything else must match to the byte.  A FILE present on one side and
missing on the other is a failure, not a skip: a report that stopped being
written is a regression the scalars cannot see.  A CASE present only in the
reference run is the same failure - it stopped running.  A case present
only in the new run is a case added since, and is reported as a note: there
is nothing for it to differ from.

Exit status is the number of differing files, so it is usable from a hook.
"""
import argparse,pathlib,re,sys

# .log carries timings, thread ids and progress chatter by design.
SKIP_SUFFIX={'.log'}
SKIP_NAME={'provenance.txt'}
CLOCK=re.compile(rb'^\s*1U(TIME|DATE)\b')

def canon(path):
    """File bytes with the clock records removed, or None if unreadable."""
    try: b=path.read_bytes()
    except OSError: return None
    if path.suffix=='.frd':
        return b"".join(l for l in b.splitlines(keepends=True)
                        if not CLOCK.match(l))
    return b

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('-v',action='store_true',help='list every file compared')
    o=ap.parse_args()
    A=pathlib.Path(o.a).resolve(); B=pathlib.Path(o.b).resolve()
    cases=sorted({d.name for d in A.iterdir() if d.is_dir()} &
                 {d.name for d in B.iterdir() if d.is_dir()})
    onlyA=sorted({d.name for d in A.iterdir() if d.is_dir()}-set(cases))
    onlyB=sorted({d.name for d in B.iterdir() if d.is_dir()}-set(cases))
    nfile=0;nbad=0;nskip=0
    for c in cases:
        fa={f.name for f in (A/c).iterdir() if f.is_file()}
        fb={f.name for f in (B/c).iterdir() if f.is_file()}
        for name in sorted(fa|fb):
            p=pathlib.Path(name)
            if p.suffix in SKIP_SUFFIX or name in SKIP_NAME: nskip+=1; continue
            if name not in fa or name not in fb:
                print("MISSING  %-24s %s  (only in %s)"
                      %(c,name,"A" if name in fa else "B")); nbad+=1; continue
            nfile+=1
            x=canon(A/c/name); y=canon(B/c/name)
            if x!=y:
                nbad+=1
                extra=""
                if x is not None and y is not None:
                    extra=" (%d vs %d bytes)"%(len(x),len(y))
                print("DIFFERS  %-24s %s%s"%(c,name,extra))
            elif o.v:
                print("same     %-24s %s"%(c,name))
    # A case present on one side only is not symmetric.  Gone from B is a
    # case that STOPPED RUNNING, which is a regression the scalars cannot
    # see.  Present only in B is a case ADDED since the reference run, which
    # cannot be a difference in output because there is nothing to differ
    # from - counting it made every comparison against an older baseline
    # report one phantom failure for ever.
    for c in onlyA: print("MISSING  case %s ran only in A - it stopped running"%c); nbad+=1
    for c in onlyB: print("note     case %s is new since the reference run"%c)
    print("\n[ABRUNS] %d case(s), %d file(s) compared byte for byte, "
          "%d skipped as clocks, %d differ"%(len(cases),nfile,nskip,nbad))
    if nbad==0:
        print("[ABRUNS] the two binaries produce identical output")
    return nbad

if __name__=='__main__':
    sys.exit(min(main(),255))
