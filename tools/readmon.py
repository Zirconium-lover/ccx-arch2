#!/usr/bin/env python3
"""Read a run's structured records and apply the playbook to them.

    tools/readmon.py <rundir-or-run.log>

handover/02-DIAGNOSTICS.md is eleven readings a human performs by eye on
prose output.  It cannot be anything else while the diagnostics are printf.
Now that src/monitor.c emits one JSON object per record, section 3 - the
stiffness census - can be a program instead, and this is it.

What it reads (per 02-DIAGNOSTICS.md section 3):

  worst ~ 3.7e-07 = gmin * (Kn*A / K_bulk)   the node has lost ALL its bulk
                                             and hangs on FAILED facets
  worst ~ 4e-02                              lost its bulk, facets intact
  worst 1e-1 .. 3e-1                         still has live elements
  nonpositive > 0                            a diagonal that is not positive

It is deliberately small.  The point is not this reader; it is that a
reading which used to require a person can now be one, for every record in
the run, without anyone remembering to look.
"""
import json,pathlib,re,sys

def records(path):
    p=pathlib.Path(path)
    if p.is_dir(): p=p/'run.log'
    for line in open(p,errors='replace'):
        if line.startswith('[MON] '):
            try: yield json.loads(line[6:])
            except ValueError: pass

# The bands 02-DIAGNOSTICS.md section 3 names, and NOTHING BETWEEN THEM.
# The document's 3.7e-07 is gmin*(Kn*A/K_bulk) computed for s3rad, so it is
# a deck constant and not a universal one; a reader that hard-coded it as a
# universal threshold would be inventing precision.  Where a value falls
# between two named bands the honest answer is that it does.
BANDS=[(3.e-7,"has lost ALL its bulk and hangs on FAILED facets "
               "(s3rad's gmin*Kn*A/K_bulk; check it for this deck)"),
       (6.e-2,"has lost its bulk but its facets are still intact"),
       (3.e-1,"still has live elements")]

def read_stiffness(r):
    """Returns (the band, whether the value is near its upper edge)."""
    w=r['worst_ratio']
    if r['nonpositive']>0:
        return "NONPOSITIVE diagonal on %d node(s)"%r['nonpositive'],False
    for i,(hi,text) in enumerate(BANDS):
        if w<hi: return "worst node "+text,(i>0 and w>0.3*hi)
    return "nothing has lost significant support",False

def main():
    if len(sys.argv)!=2: sys.exit(__doc__)
    rs=[r for r in records(sys.argv[1]) if r.get('kind')=='stiffness']
    if not rs:
        sys.exit("no [MON] stiffness records in %s.  Is CCX_DAMAGE_AUTOSPC or "
                 "CCX_DAMAGE_STIFF_PROBE set on that run?"%sys.argv[1])
    print("%d stiffness census record(s), increments %d to %d"
          %(len(rs),rs[0]['inc'],rs[-1]['inc']))
    print("\n  %-6s %-14s %-10s %-8s %s"
          %("inc","worst ratio","node","below1e-3","reading"))
    last=None
    for r in rs:
        reading,edge=read_stiffness(r)
        # print a line only when the BAND changes, plus the last record: a run
        # that says the same thing for 400 increments should say it once
        if (reading!=last) or (r is rs[-1]):
            print("  %-6d %-14.4e %-10d %-8d %s%s"
                  %(r['inc'],r['worst_ratio'],r['worst_node'],
                    r['below_1e-3'],reading,
                    "  [near the band edge - a trend, not a class]" if edge else ""))
            last=reading
    worst=min(rs,key=lambda r:r['worst_ratio'])
    print("\n  extreme over the run: node %d at %.4e (increment %d) - %s"
          %(worst['worst_node'],worst['worst_ratio'],worst['inc'],
            read_stiffness(worst)[0]))
    return 0

if __name__=='__main__': sys.exit(main())
