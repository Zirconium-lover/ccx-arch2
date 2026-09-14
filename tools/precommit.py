#!/usr/bin/env python3
"""Refuse the commits this repository has already made by accident.

    tools/precommit.py            # check the staged index
    tools/precommit.py --selftest # check the checks

WHY THIS EXISTS
---------------
Twice is chance.  Three times is a missing check.

A statically linked MKL binary is 113 to 131 MB, and `git add -A' has put
one into a commit three times on this branch's history - ccx_stepB, then
ccx_glob, then ccx_selftest.  Every time the remote's pre-receive hook
caught it, the commit had to be rewritten, and .gitignore gained one more
NAME.  A rule that lists the binary which already got through does not stop
the next one; a size limit does, whatever it is called.

There is no good generic pattern to ignore instead: `ccx_*' matches
ccx_2.23.c and would ignore source.

The second check is staleness.  src/ccxopt_list.h is generated from the
sources by tools/mkswitches.py, and a commit that edits a switch without
regenerating it leaves the binary reading a registry that does not describe
it.  The gate catches this in preflight, but the gate runs after the commit
and the commit is what gets pushed.
"""
import argparse,pathlib,subprocess,sys

ROOT=pathlib.Path(__file__).resolve().parent.parent
LIMIT=10*1024*1024          # no source file in this tree is near this

def staged():
    """(path, size) for every file staged for commit, deletions excluded."""
    out=subprocess.run(['git','-C',str(ROOT),'diff','--cached','--name-only',
                        '--diff-filter=d'],capture_output=True,text=True)
    for n in out.stdout.split('\n'):
        if not n.strip(): continue
        p=ROOT/n
        if p.is_file(): yield n,p.stat().st_size

def check_sizes(files):
    return ["%s is %.1f MB"%(n,s/1048576.) for n,s in files if s>LIMIT]

def check_switches():
    r=subprocess.run([sys.executable,str(ROOT/'tools'/'mkswitches.py'),'--check'],
                     capture_output=True,text=True)
    if r.returncode!=0:
        return ["src/ccxopt_list.h is stale; run tools/mkswitches.py"]
    return []

def selftest():
    bad=[]
    def chk(what,got,want):
        ok=(got==want)
        if not ok: bad.append(what)
        print("[PRECOMMIT]   %-50s %-16s %s"
              %(what,"%s (want %s)"%(got,want),"PASS" if ok else "FAIL"))
    chk("a file over the limit is refused",
        len(check_sizes([('src/ccx_selftest',119099384)])),1)
    chk("  and the message names it and its size",
        'ccx_selftest' in check_sizes([('src/ccx_selftest',119099384)])[0],True)
    chk("an ordinary source file is not refused",
        len(check_sizes([('src/nonlingeo.c',520000)])),0)
    chk("a file exactly at the limit is not refused",
        len(check_sizes([('x',LIMIT)])),0)
    chk("one byte over is refused",len(check_sizes([('x',LIMIT+1)])),1)
    print("\n[PRECOMMIT] self test: %d failure(s) -- %s"
          %(len(bad),"FAILED" if bad else "PASSED"))
    return len(bad)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--selftest',action='store_true')
    a=ap.parse_args()
    if a.selftest: return selftest()
    files=list(staged())
    bad=check_sizes(files)+check_switches()
    for b in bad: print("[PRECOMMIT] REFUSED: %s"%b)
    print("[PRECOMMIT] %d staged file(s), %d problem(s)"%(len(files),len(bad)))
    return 1 if bad else 0

if __name__=='__main__':
    sys.exit(main())
