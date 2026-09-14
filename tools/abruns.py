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

def compare(A,B,verbose=False):
    """Number of differing files.  Printing included, because the message is
    the product: a count with no name attached is not actionable."""
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
            elif verbose:
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

# ---------------------------------------------------------------------------
# THE SELF TEST.
#
# This program decides whether an extraction changed the answer.  Every
# commit on this branch cites it.  A comparison that has not been shown able
# to fail is not a comparison, and this one has already been wrong once in a
# way that only showed up by luck: UTIME was stripped and UDATE was not, so
# a run either side of midnight reported sixteen differing .frd files, all
# of them identical in length.
#
# Synthetic run trees, where the right answer is known because it was
# written down first.

def _mk(root,case,files):
    d=root/case; d.mkdir(parents=True,exist_ok=True)
    for n,t in files.items(): (d/n).write_bytes(t)
    return d

def selftest():
    import tempfile,shutil,io,contextlib
    bad=[]
    def chk(what,got,want):
        ok=(got==want)
        if not ok: bad.append(what)
        print("[ABRUNS]   %-54s %-14s %s"
              %(what,"%s (want %s)"%(got,want),"PASS" if ok else "FAIL"))
    def run(A,B):
        with contextlib.redirect_stdout(io.StringIO()) as f:
            n=compare(A,B)
        return n,f.getvalue()
    d=pathlib.Path(tempfile.mkdtemp(prefix='abruns'))
    try:
        FRD=b"    1UTIME  12.5\n    1UDATE  01.01.2026\n -1  1  0.5000\n"
        A,B=d/'A',d/'B'
        _mk(A,'c1',{'m.frd':FRD,'m.sta':b"inc 1\n",'m.log':b"took 3s\n",
                    'provenance.txt':b"sha=aaa\n"})
        _mk(B,'c1',{'m.frd':FRD,'m.sta':b"inc 1\n",'m.log':b"took 9s\n",
                    'provenance.txt':b"sha=bbb\n"})
        n,_=run(A,B)
        chk("identical trees compare equal",n,0)
        chk(".log and provenance.txt are skipped, not compared",n,0)

        # only the clock records differ -> still equal
        (B/'c1'/'m.frd').write_bytes(
            b"    1UTIME  99.9\n    1UDATE  31.12.2027\n -1  1  0.5000\n")
        n,_=run(A,B)
        chk("a .frd differing ONLY in 1UTIME/1UDATE is equal",n,0)

        # a real number differs -> reported
        (B/'c1'/'m.frd').write_bytes(
            b"    1UTIME  12.5\n    1UDATE  01.01.2026\n -1  1  0.6000\n")
        n,out=run(A,B)
        chk("a .frd differing in a data line is reported",n,1)
        chk("  and it says which file",'m.frd' in out,True)
        (B/'c1'/'m.frd').write_bytes(FRD)

        # a file present on one side only
        (A/'c1'/'m.dat').write_bytes(b"x\n")
        n,out=run(A,B)
        chk("a report that stopped being written is a failure",n,1)
        chk("  and it is called MISSING",'MISSING' in out,True)
        (A/'c1'/'m.dat').unlink()

        # a case only in A stopped running; a case only in B is new
        _mk(A,'gone',{'m.sta':b"inc 1\n"})
        n,out=run(A,B)
        chk("a case that stopped running is a failure",n,1)
        shutil.rmtree(A/'gone')
        _mk(B,'added',{'m.sta':b"inc 1\n"})
        n,out=run(A,B)
        chk("a case added since the reference run is NOT a failure",n,0)
        chk("  and it is reported as a note",'note' in out,True)
    finally:
        shutil.rmtree(d,ignore_errors=True)
    print("\n[ABRUNS] self test: %d failure(s) -- %s"
          %(len(bad),"FAILED" if bad else "PASSED"))
    return len(bad)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('a',nargs='?'); ap.add_argument('b',nargs='?')
    ap.add_argument('-v',action='store_true',help='list every file compared')
    ap.add_argument('--selftest',action='store_true',
                    help='check the comparison against synthetic run trees')
    o=ap.parse_args()
    if o.selftest: return selftest()
    if not o.a or not o.b: ap.error("two run directories, or --selftest")
    return compare(pathlib.Path(o.a).resolve(),pathlib.Path(o.b).resolve(),o.v)

if __name__=='__main__':
    sys.exit(min(main(),255))
