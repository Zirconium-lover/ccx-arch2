#!/usr/bin/env python3
"""Where does the run's time go?  Run one gate case under the event timer.

    CCX_EXE=/path/to/ccx_2.23_pardiso tools/profile.py fast-wrapped

Reuses test/regress/cases.json so a profile is taken of a case whose purpose
somebody has written down, and pins OMP_NUM_THREADS/MKL_NUM_THREADS to 1 by
default: a timing comparison across thread counts is meaningless here, and so
is a timing comparison against a machine that was busy doing something else.

Cases are run ONE AT A TIME on purpose.  test/regress/run.py runs them
concurrently, which is right for a pass/fail gate and wrong for a stopwatch.

The solver's own instrument (src/logview.c) does the attribution; this script
only arranges the run, refuses to believe a perturbed one, and prints the
result.  --json writes the record so two arms can be diffed.
"""
import argparse,json,os,pathlib,re,shutil,subprocess,sys,time

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'test/regress'))

def load_cases():
    return {c['name']:c for c in json.load(open(ROOT/'test/regress/cases.json'))['cases']}

def run(case,rundir,exe,threads,extra_env):
    env=dict(os.environ)
    env['OMP_NUM_THREADS']=str(threads); env['MKL_NUM_THREADS']=str(threads)
    env['CCX_LOG_VIEW']='1'; env['CCX_EXE']=exe
    rundir.mkdir(parents=True,exist_ok=True)
    if case['kind']=='fast':
        # POSITIONAL overrides, not environment: run_fast.sh delegates to
        # run_s3rad.sh, which exports twelve CCX_* names unconditionally and
        # would silently clobber anything set the other way.  Getting this
        # wrong is how an A/B of CCX_PARDISO_REUSE_SYMBOLIC came out with the
        # switch ON in BOTH arms and produced a published, wrong refutation.
        env['FAST_VARIANT']=case['variant']
        cmd='%s/test/fast/run_fast.sh %s %s'%(ROOT,rundir,
                                              " ".join(case['env']+extra_env))
        cwd=None
    else:
        # close and mixed run the binary directly; nothing clobbers them
        for kv in case['env']+extra_env:
            k,_,v=kv.partition('='); env[k]=v
        env.setdefault('CCX_DAMAGE_AUTOSPC','1.e-3')
        if case['kind']=='close':
            subprocess.run('python3 %s/test/pathfollow/mkclose.py -o %s'
                           %(ROOT,rundir/'close.inp'),shell=True,check=True)
            cmd='%s -i close > run.log 2>&1'%exe
        else:
            shutil.copy(ROOT/'test/pathfollow/mixed.inp',rundir/'mixed.inp')
            cmd='%s -i mixed > run.log 2>&1'%exe
        cwd=rundir
    t0=time.time()
    subprocess.run(cmd,shell=True,env=env,cwd=cwd,stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)
    return time.time()-t0,rundir/'run.log'

def switches(log):
    """Which CCX_* names were in force, from the run's own [SWITCHES] banner.

    This exists because of a measurement that went wrong here and produced a
    confident, published, WRONG conclusion.  test/fast/run_fast.sh delegates
    to run_s3rad.sh, which exports CCX_PARDISO_REUSE_SYMBOLIC=1, so an A/B
    run as "stock against the switch" had the switch on in BOTH arms.  The
    5% difference measured was noise; the real effect is 29%.

    handover/02-DIAGNOSTICS.md section 9 says to check this block before
    believing any comparison.  A human who has to remember will not.  So the
    tool reads it and refuses to be quiet when two arms are the same arm."""
    txt=open(log,errors='replace').read()
    return {m.group(1):m.group(2).strip() for m in
            re.finditer(r'^\[SWITCHES\]   (CCX_[A-Z0-9_]+) = (.*?)(?:   \[|$)',
                        txt,re.M)}

def parse(log):
    txt=open(log,errors='replace').read()
    m=re.search(r'^\[LOGVIEW_JSON\] (.*)$',txt,re.M)
    if not m:
        bad=re.findall(r'^\[LOGVIEW\].*$',txt,re.M)
        raise SystemExit("no profile in %s%s"%(log,
            "\n  "+"\n  ".join(bad) if bad else
            "\n  the run printed no [LOGVIEW] block at all: is CCX_LOG_VIEW read "
            "by this binary?  check the [SWITCHES] banner."))
    return json.loads(m.group(1))

def show(name,prof,wall):
    ev=sorted(prof['events'],key=lambda e:-e['self'])
    tot=prof['total']
    acc=sum(e['self'] for e in ev)
    print("\n%s   %.2f s measured by the solver, %.2f s of wall clock"%(name,tot,wall))
    print("  %-28s %8s %10s %10s %7s %10s"%("event","calls","incl (s)","self (s)","%run","ms/call"))
    for e in ev:
        print("  %-28s %8d %10.3f %10.3f %6.1f%% %10.2f"
              %(e['name'],e['calls'],e['incl'],e['self'],100.*e['self']/tot,
                1000.*e['incl']/e['calls']))
    print("  %-28s %8s %10s %10.3f %6.1f%%"
          %("not instrumented","","",tot-acc,100.*(tot-acc)/tot))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('case',nargs='+',help='case name(s) from test/regress/cases.json')
    ap.add_argument('-t','--threads',type=int,default=1)
    ap.add_argument('-e','--env',action='append',default=[],metavar='K=V',
                    help='extra environment, e.g. CCX_PARDISO_REUSE_SYMBOLIC=1')
    ap.add_argument('-o',default=None,help='where to put the runs')
    ap.add_argument('--json',default=None,help='write the profiles here')
    ap.add_argument('--compare',default=None,metavar='CASE',
                    help='treat CASE as the control arm: print which switches '
                         'the other arms differ in, and FAIL if any of them is '
                         'configured identically to it')
    a=ap.parse_args()
    exe=os.environ.get('CCX_EXE')
    if not exe or not os.access(exe,os.X_OK):
        sys.exit("set CCX_EXE to a PARDISO-enabled ccx_2.23 binary")
    cases=load_cases()
    outroot=pathlib.Path(a.o or (ROOT/'test/regress/_runs'/
             ('profile-'+time.strftime('%Y%m%d-%H%M%S')))).resolve()
    out={}
    for name in a.case:
        if name not in cases: sys.exit("no such case: %s"%name)
        wall,log=run(cases[name],outroot/name,exe,a.threads,a.env)
        prof=parse(log)
        prof['wall']=wall; prof['threads']=a.threads; prof['env']=a.env
        prof['switches']=switches(log)
        out[name]=prof
        show(name,prof,wall)
    # An A/B whose two arms are configured identically is not wrong, it is
    # UNINFORMATIVE, and that is the expensive kind: it costs a whole run to
    # notice, and it can cost a published conclusion if nobody does.
    if a.compare is not None:
        base=out.get(a.compare)
        if base is None:
            sys.exit("--compare names a case that was not run: %s"%a.compare)
        for name,prof in out.items():
            if name==a.compare: continue
            if prof['switches']==base['switches']:
                print("\n*** %s and %s were run with IDENTICAL switches.  Any "
                      "difference between them is noise, not an effect.\n"
                      "    In force in both: %s"
                      %(name,a.compare,
                        ", ".join("%s=%s"%kv for kv in sorted(base['switches'].items()))
                        or "(none)"))
                return 2
            diff=sorted(set(prof['switches'].items())^set(base['switches'].items()))
            print("\n%s vs %s differ only in: %s"
                  %(name,a.compare,", ".join("%s=%s"%kv for kv in diff)))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
        print("\nwrote %s"%a.json)

if __name__=='__main__': sys.exit(main() or 0)
