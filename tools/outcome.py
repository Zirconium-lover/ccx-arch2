#!/usr/bin/env python3
"""Run a whole case under named arms and compare what the run DID.

    CCX_EXE=... tools/outcome.py --variant wrapped \
        --arm stock CCX_DAMAGE_TANGENT= \
        --arm unsym CCX_DAMAGE_TANGENT=UNSYM

opcheck.py answers "is the operator the differential".  This answers the
question that decides whether that matters: does the run come out different.
It reports the physics (last increment, theta, deletion count and set) and
the cost (attempts, Newton iterations) side by side.

It exists because the same mistake was made twice in this project: an A/B
whose two arms were configured identically, because test/fast/run_fast.sh
delegates to run_s3rad.sh, which `export`s twelve CCX_* names unconditionally
and only then applies its positional overrides.  An empty override is how you
UNSET one; leaving it out of the environment does nothing.  So this script
passes overrides positionally, reads back each run's own [SWITCHES] banner,
and REFUSES to report a comparison whose arms ended up identically
configured.
"""
import argparse,json,os,pathlib,re,shutil,subprocess,sys,time

ROOT=pathlib.Path(__file__).resolve().parent.parent

def switches(log):
    txt=open(log,errors='replace').read()
    return {m.group(1):m.group(2).strip() for m in
            re.finditer(r'^\[SWITCHES\]   (CCX_[A-Z0-9_]+) = (.*?)(?:   \[|$)',txt,re.M)}

def sta(p):
    rows=[l.split() for l in open(p,errors='replace') if len(l.split())==7
          and l.split()[0].isdigit()]
    return rows

def outcome(d):
    rows=sta(d/'m.sta') if (d/'m.sta').exists() else []
    o={'attempts':len(rows),
       'last_inc':int(rows[-1][1]) if rows else None,
       'theta':rows[-1][4] if rows else None}
    try: o['iters']=sum(1 for l in open(d/'m.cvg',errors='replace')
                        if len(l.split())>=7 and l.split()[0].isdigit())
    except OSError: o['iters']=None
    try:
        dl=[l.split()[0] for l in open(d/'m.damage',errors='replace')
            if not l.startswith('#') and l.strip()]
        o['deletions']=len(dl); o['deletion_set']=set(dl)
    except OSError: o['deletions']=None; o['deletion_set']=set()
    return o

def main():
    ap=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--variant',default='wrapped',choices=['wrapped','plain'])
    ap.add_argument('--threads',type=int,default=1)
    ap.add_argument('-o',default=None)
    ap.add_argument('--json',default=None)
    ap.add_argument('--arm',action='append',nargs='+',default=[],metavar='NAME K=V')
    a=ap.parse_args()
    exe=os.environ.get('CCX_EXE')
    if not exe or not os.access(exe,os.X_OK):
        sys.exit("set CCX_EXE to a PARDISO-enabled ccx_2.23 binary")
    outroot=pathlib.Path(a.o or ('/tmp/outcome-'+time.strftime('%H%M%S'))).resolve()
    outroot.mkdir(parents=True,exist_ok=True)
    res={}; sw={}
    print("variant %s, %d thread(s)\n"%(a.variant,a.threads))
    print("  %-16s %-5s %-14s %-6s %-9s %-7s %s"
          %("arm","rc","theta","inc","attempts","iters","deletions"))
    for arm in a.arm:
        name,over=arm[0],arm[1:]
        d=outroot/name; shutil.rmtree(d,ignore_errors=True); d.mkdir(parents=True)
        env=dict(os.environ)
        env['OMP_NUM_THREADS']=str(a.threads); env['MKL_NUM_THREADS']=str(a.threads)
        env['FAST_VARIANT']=a.variant; env['CCX_EXE']=exe
        r=subprocess.run('%s/test/fast/run_fast.sh %s %s'%(ROOT,d," ".join(over)),
                         shell=True,env=env,stdout=subprocess.DEVNULL,
                         stderr=subprocess.STDOUT)
        got=switches(d/'run.log') if (d/'run.log').exists() else {}
        bad=[]
        for kv in over:
            k,_,v=kv.partition('=')
            if v=='':
                if k in got: bad.append("%s asked UNSET, run has %s"%(k,got[k]))
            elif got.get(k)!=v: bad.append("%s asked %r, run has %r"%(k,v,got.get(k)))
        if bad:
            print("  %-16s CONFIGURATION DID NOT REACH THE BINARY: %s"%(name,"; ".join(bad)))
            continue
        o=outcome(d); o['rc']=r.returncode; res[name]=o
        sw[name]=tuple(sorted(got.items()))
        print("  %-16s %-5d %-14s %-6s %-9d %-7s %s"
              %(name,o['rc'],o['theta'],o['last_inc'],o['attempts'],
                o['iters'],o['deletions']))
    names=list(sw)
    dup=False
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            if sw[names[i]]==sw[names[j]]:
                print("\n*** %s and %s ran with IDENTICAL switches: any difference "
                      "between them is noise, and any AGREEMENT is a tautology"
                      %(names[i],names[j])); dup=True
    if len(names)>1 and not dup:
        base=names[0]
        print("\n  deletion set against %s:"%base)
        for n in names[1:]:
            d0,d1=res[base]['deletion_set'],res[n]['deletion_set']
            print("    %-16s %d only in %s, %d only in %s"
                  %(n,len(d0-d1),base,len(d1-d0),n))
    if a.json:
        j={k:{kk:vv for kk,vv in v.items() if kk!='deletion_set'}
           for k,v in res.items()}
        pathlib.Path(a.json).write_text(json.dumps(j,indent=2,sort_keys=True)+"\n")
    print("\nruns in %s"%outroot)
    return 2 if dup else 0

if __name__=='__main__': sys.exit(main())
