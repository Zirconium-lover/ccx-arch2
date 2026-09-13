#!/usr/bin/env python3
"""Run the operator check under a named configuration and report the verdict.

    CCX_EXE=... tools/opcheck.py --inc 50 --elem 1238 \
        --arm shipped CCX_DAMAGE_TANGENT=UNSYM CCX_DAMAGE_VISCOSITY=1.e-4 \
        --arm no-visc CCX_DAMAGE_TANGENT=UNSYM CCX_DAMAGE_VISCOSITY=0

Why this exists rather than another shell loop: the probe is run dozens of
times while a tangent is being chased, every run is an A/B, and this project
has already paid once for an A/B whose two arms were configured identically
(handover/04-REFUTED.md).  So this script

  - passes overrides POSITIONALLY, because test/fast/run_fast.sh delegates to
    run_s3rad.sh, which exports twelve CCX_* names unconditionally and would
    clobber anything set in the environment;
  - VERIFIES against each run's own [SWITCHES] banner that what was asked for
    is what the binary got, and fails the arm if not;
  - refuses a comparison whose arms end up with identical switch sets.

It reads the machine-readable [MON] records, not the human table.
"""
import argparse,json,os,pathlib,re,shutil,subprocess,sys,time

ROOT=pathlib.Path(__file__).resolve().parent.parent

def gen_deck(dst,variant):
    """The fast deck, generated rather than committed, as its README requires."""
    dst.parent.mkdir(parents=True,exist_ok=True)
    args={'plain':'','wrapped':'--incl 6,3,0,2,3,3'}.get(variant)
    if args is None: sys.exit("unknown variant %s"%variant)
    r=subprocess.run('python3 %s/test/fast/mkfast.py -o %s --nx 10 --ny 6 --nz 6 '
                     '--from-deck %s/test/s3rad/m12_s3rad_gc24_w.inp %s'
                     %(ROOT,dst,ROOT,args),shell=True,capture_output=True,text=True)
    if r.returncode!=0: sys.exit("mkfast failed:\n"+r.stdout+r.stderr)

def switches(log):
    txt=open(log,errors='replace').read()
    return {m.group(1):m.group(2).strip() for m in
            re.finditer(r'^\[SWITCHES\]   (CCX_[A-Z0-9_]+) = (.*?)(?:   \[|$)',txt,re.M)}

def records(log,kind):
    txt=open(log,errors='replace').read()
    out=[]
    for m in re.finditer(r'^\[MON\] (\{.*\})$',txt,re.M):
        try: r=json.loads(m.group(1))
        except ValueError: continue
        if r.get('kind')==kind: out.append(r)
    return out

def run_arm(name,deck,env_over,a,outroot):
    d=outroot/name; shutil.rmtree(d,ignore_errors=True); d.mkdir(parents=True)
    shutil.copy(deck,d/'m.inp')
    env=dict(os.environ)
    env['OMP_NUM_THREADS']='1'; env['MKL_NUM_THREADS']='1'
    fixed=['CCX_STRUCT_FD_INC=%d'%a.inc,'CCX_STRUCT_FD_ITER=%d'%a.iter,
           'CCX_STRUCT_FD_H=%s'%a.h,'CCX_STRUCT_FD_BASE=%s'%a.base]
    if a.elem: fixed.append('CCX_STRUCT_FD_ELEM=%d'%a.elem)
    if a.step: fixed.append('CCX_STRUCT_FD_STEP=%d'%a.step)
    want={}
    for kv in fixed+env_over:
        k,_,v=kv.partition('='); env[k]=v; want[k]=v
    subprocess.run('%s -i m > run.log 2>&1'%a.exe,shell=True,env=env,cwd=d)
    log=d/'run.log'
    got=switches(log)
    bad=[]
    for k,v in want.items():
        if v=='':
            if k in got: bad.append("%s asked UNSET, run has %s"%(k,got[k]))
        elif got.get(k)!=v:
            bad.append("%s asked %r, run has %r"%(k,v,got.get(k)))
    return log,got,bad

def main():
    ap=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--deck',default=None,help='deck file; default: generate the fast variant')
    ap.add_argument('--variant',default='wrapped')
    ap.add_argument('--inc',type=int,required=True)
    ap.add_argument('--iter',type=int,default=2)
    ap.add_argument('--step',type=int,default=0)
    ap.add_argument('--elem',type=int,default=0)
    ap.add_argument('--h',default='1e-7')
    ap.add_argument('--base',default='VOLD',choices=['V','VOLD'])
    ap.add_argument('-o',default=None)
    ap.add_argument('--json',default=None)
    ap.add_argument('--arm',action='append',nargs='+',default=[],metavar='NAME K=V',
                    help='an arm: its name, then its overrides')
    a=ap.parse_args()
    a.exe=os.environ.get('CCX_EXE')
    if not a.exe or not os.access(a.exe,os.X_OK):
        sys.exit("set CCX_EXE to a PARDISO-enabled ccx_2.23 binary")
    outroot=pathlib.Path(a.o or ('/tmp/opcheck-'+time.strftime('%H%M%S'))).resolve()
    outroot.mkdir(parents=True,exist_ok=True)
    deck=pathlib.Path(a.deck).resolve() if a.deck else outroot/'gen.inp'
    if not a.deck: gen_deck(deck,a.variant)
    if not a.arm: a.arm=[['stock']]

    print("deck %s   increment %d iteration %d   base %s   h %s%s"
          %(deck.name,a.inc,a.iter,a.base,a.h,
            "   element %d"%a.elem if a.elem else "   element: most damaged"))
    print("\n  %-26s %-9s %-11s %-11s %-11s %6s %6s %6s %6s"
          %("arm","dam","|ctr-asm|","|fwd-bwd|","min|sd-asm|","ok","kink","wrong","both"))
    out={}; swsets={}
    for arm in a.arm:
        name,over=arm[0],arm[1:]
        log,got,bad=run_arm(name,deck,over,a,outroot)
        if bad:
            print("  %-26s CONFIGURATION DID NOT REACH THE BINARY: %s"%(name,"; ".join(bad)))
            out[name]={'error':bad}; continue
        tot=records(log,'opcheck'); col=records(log,'opcheck_column')
        if not tot:
            why=[l for l in open(log,errors='replace') if l.startswith('[OPCHECK]')]
            print("  %-26s no verdict%s"%(name," - "+why[0].strip() if why else ""))
            out[name]={'error':'no verdict'}; continue
        b=[t for t in tot if t.get('element')] or tot
        b=b[0]
        cc=[c for c in col] or [{}]
        worst=max(cc,key=lambda c:c.get('rel_ctr',0))
        print("  %-26s %-9.4f %-11.4e %-11.4e %-11.4e %6d %6d %6d %6d"
              %(name,b.get('dam',float('nan')),worst.get('rel_ctr',0),
                worst.get('rel_side',0),worst.get('rel_best',0),
                b['ok'],b['kink'],b['wrong'],b['both']))
        out[name]={'total':b,'worst_column':worst,'switches':got}
        swsets[name]=tuple(sorted(got.items()))
    names=[n for n in swsets]
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            if swsets[names[i]]==swsets[names[j]]:
                print("\n*** %s and %s ran with IDENTICAL switches: any difference "
                      "between them is noise"%(names[i],names[j]))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
        print("\nwrote %s"%a.json)
    print("\nruns in %s"%outroot)
    return 0

if __name__=='__main__': sys.exit(main())
