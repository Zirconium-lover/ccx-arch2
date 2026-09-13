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

  * the UTIME record of a .frd file, which is the wall-clock time the run
    started.  frdheader.c writes it from time(); two runs of the identical
    binary differ there too, so treating it as a difference would make the
    check cry wolf on every pair;
  * provenance.txt, which records the binary's sha256 on purpose - it is
    SUPPOSED to change when the binary changes.

Everything else must match to the byte.  A file present on one side and
missing on the other is a failure, not a skip: a report that stopped being
written is a regression the scalars cannot see.

Exit status is the number of differing files, so it is usable from a hook.
"""
import argparse,pathlib,re,sys

# .log carries timings, thread ids and progress chatter by design.
SKIP_SUFFIX={'.log'}
SKIP_NAME={'provenance.txt'}
UTIME=re.compile(rb'^\s*1UTIME\b')

def canon(path):
    """File bytes with the clock records removed, or None if unreadable."""
    try: b=path.read_bytes()
    except OSError: return None
    if path.suffix=='.frd':
        return b"".join(l for l in b.splitlines(keepends=True)
                        if not UTIME.match(l))
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
    for c in onlyA: print("MISSING  case %s ran only in A"%c); nbad+=1
    for c in onlyB: print("MISSING  case %s ran only in B"%c); nbad+=1
    print("\n[ABRUNS] %d case(s), %d file(s) compared byte for byte, "
          "%d skipped as clocks, %d differ"%(len(cases),nfile,nskip,nbad))
    if nbad==0:
        print("[ABRUNS] the two binaries produce identical output")
    return nbad

if __name__=='__main__':
    sys.exit(min(main(),255))
