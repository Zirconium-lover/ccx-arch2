#!/usr/bin/env python3
"""Run every regression that finishes in minutes and compare it to the record.

    CCX_EXE=/path/to/ccx_2.23_pardiso test/regress/run.py [-j N] [-k NAME]

Exit status is the number of failed cases, so it is usable from a hook or a
CI step.  Cases are defined in cases.json, which also carries what each one
is FOR - a case whose purpose nobody can state is a case nobody will fix.

Why this exists: the architecture audit's finding was that validation cost
2.5 hours, so nothing could be checked cheaply and every fix was a point fix
made blind.  These cases cover the load-path judgement, the crack-face kink,
bulk damage with deletion and cutbacks, and the analytical mixed-mode branch,
and they run on one core in the time it takes to read a diff.
"""
import argparse,concurrent.futures as cf,io,contextlib,json,os,pathlib,re,shutil,subprocess,sys,time

HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parent.parent
sys.path.insert(0,str(ROOT/'tools'))
import ccxdiff

def sh(cmd,env,cwd=None,timeout=3600):
    return subprocess.run(cmd,shell=True,env=env,cwd=cwd,timeout=timeout,
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)

def base_env(extra):
    e=dict(os.environ)
    e.setdefault('OMP_NUM_THREADS','1'); e.setdefault('MKL_NUM_THREADS','1')
    for kv in extra:
        k,_,v=kv.partition('='); e[k]=v
    return e

def last_sta(p):
    try: lines=[l for l in open(p) if l.strip()]
    except OSError: return None,None
    if not lines: return None,None
    f=lines[-1].split()
    return int(f[1]),f[4]

def ndel(p):
    try: return sum(1 for l in open(p) if not l.startswith('#') and l.strip())
    except OSError: return None

def delset(p):
    try: return {l.split()[0] for l in open(p) if not l.startswith('#') and l.strip()}
    except OSError: return set()

def census(log):
    """max nodes masked at the standard threshold, and the worst ratio."""
    mx=0; worst=None
    try: txt=open(log,errors='replace').read()
    except OSError: return 0,None
    for m in re.finditer(r'below_1e-3=(\d+)',txt): mx=max(mx,int(m.group(1)))
    for m in re.finditer(r'worst=node_\d+_at_([0-9.eE+-]+)',txt):
        v=float(m.group(1))
        if worst is None or v<worst: worst=v
    return mx,worst

def check_switches(log,want):
    """Did the switches the case asked for actually reach the binary?

    Reads the run's own [SWITCHES] banner - 02-DIAGNOSTICS.md section 9 says
    to check it before believing any comparison, and a check a human has to
    remember is a check that does not happen."""
    bad=[]
    try: txt=open(log,errors='replace').read()
    except OSError: return ["no log"]
    got={m.group(1):m.group(2).strip() for m in
         re.finditer(r'^\[SWITCHES\]   (CCX_[A-Z0-9_]+) = (.*?)(?:   \[|$)',txt,re.M)}
    for kv in want:
        k,_,v=kv.partition('=')
        if v=='':
            if k in got:
                bad.append("%s was asked to be UNSET and the run has it = %s"
                           %(k,got[k]))
        elif got.get(k)!=v:
            bad.append("%s was asked for as %r and the run has it as %r"
                       %(k,v,got.get(k)))
    return bad

def selftests(log,required,lines):
    """every named self test must report PASSED, no unit may report a failure
    or an error, and every required line must be present - a report that
    silently stopped being emitted is a regression too."""
    try: txt=open(log,errors='replace').read()
    except OSError: return ["no log"]
    bad=[]
    for name in required:
        if not re.search(re.escape(name)+r'.*(PASSED|0 failure)',txt):
            if name in txt: bad.append("%s did not report PASSED"%name)
    for m in re.finditer(r'\[([A-Z0-9 _]+)\][^\n]*?([1-9]\d*) failure',txt):
        bad.append("%s reported %s failure(s)"%(m.group(1),m.group(2)))
    for m in re.finditer(r'\[([A-Z0-9 _]+)\][^\n]*?\*ERROR([^\n]*)',txt):
        bad.append("%s reported an error:%s"%(m.group(1),m.group(2)[:80]))
    for want in lines:
        if want not in txt: bad.append("log does not contain %r"%want)
    return bad

def run_fast(case,rundir,exe):
    """A case's env goes through as POSITIONAL overrides, not as environment.

    run_fast.sh delegates to run_s3rad.sh, which `export`s twelve CCX_* names
    unconditionally and only then applies the NAME=VALUE arguments.  So a
    case that set one of those twelve in the environment was silently
    overwritten, and fast-wrapped-nospc - whose whole stated purpose is to
    run WITHOUT the load-path mask - ran with CCX_DAMAGE_AUTOSPC=1.e-3, i.e.
    as a second copy of fast-wrapped.  It passed every time, and proved
    nothing.  Not wrong: UNINFORMATIVE, which is the expensive kind.

    The `[SWITCHES]` banner said so in every one of those logs.  Nobody read
    it, which is why check_switches() below now does."""
    env=base_env([]); env['FAST_VARIANT']=case['variant']; env['CCX_EXE']=exe
    args=" ".join(case['env'])
    r=sh('%s/test/fast/run_fast.sh %s %s'%(ROOT,rundir,args),env)
    return r.returncode,rundir/'run.log',rundir/'m.sta',rundir/'m.damage'

def run_close(case,rundir,exe):
    """the single-facet closing benchmark: verifies the NORMAL law itself."""
    rundir.mkdir(parents=True,exist_ok=True)
    r=sh('python3 %s/test/pathfollow/mkclose.py -o %s'%(ROOT,rundir/'close.inp'),
         base_env([]))
    if r.returncode!=0: return r.returncode,rundir/'run.log',None,None
    env=base_env(case['env']); env.setdefault('CCX_DAMAGE_AUTOSPC','1.e-3')
    r=sh('%s -i close > run.log 2>&1'%exe,env,cwd=rundir)
    return r.returncode,rundir/'run.log',rundir/'close.sta',None

def run_mixed(case,rundir,exe):
    rundir.mkdir(parents=True,exist_ok=True)
    shutil.copy(ROOT/'test/pathfollow/mixed.inp',rundir/'mixed.inp')
    env=base_env(case['env']); env.setdefault('CCX_DAMAGE_AUTOSPC','1.e-3')
    r=sh('%s -i mixed > run.log 2>&1'%exe,env,cwd=rundir)
    return r.returncode,rundir/'run.log',rundir/'mixed.sta',None

def one(case,outroot,exe,required,lines):
    rundir=outroot/case['name']
    if rundir.exists(): shutil.rmtree(rundir)
    t0=time.time()
    if   case['kind']=='fast':  rc,log,sta,dam=run_fast(case,rundir,exe)
    elif case['kind']=='close': rc,log,sta,dam=run_close(case,rundir,exe)
    else:                       rc,log,sta,dam=run_mixed(case,rundir,exe)
    got={'rc':rc}
    inc,theta=last_sta(sta)
    got['last_inc'],got['theta']=inc,theta
    if dam is not None: got['deletions']=ndel(dam)
    exp=case['expect']
    if 'masked_max' in exp or 'worst_ratio' in exp:
        got['masked_max'],got['worst_ratio']=census(log)
    if 'check_close' in exp:
        r=sh('python3 %s/test/pathfollow/check_close.py %s --zeta %s'
             %(ROOT,rundir,case.get('zeta','0')),base_env([]))
        got['check_close']='PASSED' if 'PASSED' in r.stdout else 'FAILED'
        m=re.search(r'ratio ([0-9.e+-]+) \(1/g=([0-9.e+-]+)\)',r.stdout)
        if m: got['tangent_ratio']=float(m.group(1))
        m=re.search(r'inside the blend band: (\d+)',r.stdout)
        if m: got['in_band']=int(m.group(1))
        m=re.search(r'worst relative error against the law: ([0-9.e+-]+)',r.stdout)
        if m: got['law_error']=float(m.group(1))
    if 'check_mixed' in exp:
        r=sh('python3 %s/test/pathfollow/check_mixed.py %s'%(ROOT,rundir),base_env([]))
        got['check_mixed']='PASSED' if 'PASSED' in r.stdout else 'FAILED'
        m=re.search(r'(\d+) accepted increments, (\d+) of them post-peak',r.stdout)
        if m: got['accepted'],got['postpeak']=int(m.group(1)),int(m.group(2))
    fails=[]
    for k,want in exp.items():
        have=got.get(k)
        if k in ('worst_ratio','tangent_ratio'):
            ok = have is not None and abs(have-want)<=1e-2*abs(want)
        elif k=='law_error':
            ok = have is not None and have<=want
        else:
            ok = have==want
        if not ok: fails.append("%s: want %r, got %r"%(k,want,have))
    fails+=selftests(log,required,lines)
    fails+=check_switches(log,case['env'])
    return {'name':case['name'],'what':case['what'],'seconds':round(time.time()-t0,1),
            'got':got,'fails':fails,'rundir':str(rundir)}

def provenance(exe):
    """What a later reader needs in order to believe a recorded number.

    A baseline that does not state its binary and its thread count is not a
    baseline: runs are NOT reproducible across thread counts (measured -
    handover/05-DEBT.md item 5), so a scalar without that context cannot be
    compared against anything."""
    import hashlib,platform
    def cmd(c):
        try: return subprocess.run(c,shell=True,capture_output=True,text=True).stdout.strip().split('\n')[0]
        except Exception: return None
    h=hashlib.sha256()
    with open(exe,'rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''): h.update(chunk)
    e=base_env([])
    return {'binary':os.path.realpath(exe),
            'binary_sha256':h.hexdigest(),
            'omp_num_threads':e.get('OMP_NUM_THREADS'),
            'mkl_num_threads':e.get('MKL_NUM_THREADS'),
            'mkl_cbwr':e.get('MKL_CBWR'),
            'git_commit':cmd('git -C %s rev-parse HEAD'%ROOT),
            'git_dirty':bool(cmd('git -C %s status --porcelain'%ROOT)),
            'cc':cmd('gcc --version'),'fc':cmd('gfortran --version'),
            'machine':platform.machine(),'python':platform.python_version()}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('-j',type=int,default=2,help='cases to run at once')
    ap.add_argument('-k',default=None,help='run only cases whose name contains this')
    ap.add_argument('-o',default=None,help='where to put the runs')
    ap.add_argument('--json',default=None,metavar='PATH',
                    help='also write the whole run - provenance, every '
                         'expectation and every measured scalar - as JSON, '
                         'so "it still passes" can be checked rather than '
                         'asserted')
    ap.add_argument('--record-coverage',action='store_true',
                    help='rewrite test/regress/covered.txt from this run, so '
                         'docs/SWITCHES.md can say which switches any test '
                         'actually exercises')
    a=ap.parse_args()
    exe=os.environ.get('CCX_EXE')
    if not exe or not os.access(exe,os.X_OK):
        sys.exit("set CCX_EXE to a PARDISO-enabled ccx_2.23 binary")
    # Preflight: the switch registry is generated from the sources, so a
    # new switch that nobody regenerated would make every run.log understate
    # its own configuration.  Cheap, and it fails before any solver runs.
    pre=sh('python3 %s/tools/mkswitches.py --check'%ROOT,base_env([]))
    print("preflight  switch registry: %s"%pre.stdout.strip().replace('\n','; '))
    preflight_bad = 1 if pre.returncode!=0 else 0
    # A comparison that has not been shown able to fail is not a comparison,
    # and it is about to decide whether cases pass.
    cd=sh('python3 %s/tools/ccxdiff.py --selftest'%ROOT,base_env([]))
    print("preflight  comparison layer: %s"%
          (cd.stdout.strip().splitlines() or ['no output'])[-1])
    if cd.returncode!=0: preflight_bad+=1

    # The path-follower diagnostic judges a mechanism nobody could judge
    # before; a diagnostic that has itself gone wrong is worse than none,
    # so its own self test runs here with the others.
    pf=sh('python3 %s/tools/pathfollow.py --selftest'%ROOT,base_env([]))
    print("preflight  path-follower diagnostic: %s"%
          (pf.stdout.strip().splitlines()[-1] if pf.stdout.strip() else "no output"))
    if pf.returncode!=0: preflight_bad+=1
    spec=json.load(open(HERE/'cases.json'))
    cases=[c for c in spec['cases'] if not a.k or a.k in c['name']]
    outroot=pathlib.Path(a.o or (HERE/'_runs'/time.strftime('%Y%m%d-%H%M%S'))).resolve()
    outroot.mkdir(parents=True,exist_ok=True)
    print("binary %s\nruns   %s\ncases  %d, %d at a time\n"%(exe,outroot,len(cases),a.j))
    res={}
    with cf.ThreadPoolExecutor(max_workers=a.j) as ex:
        futs={ex.submit(one,c,outroot,exe,spec['selftests_required'],
                        spec.get('required_lines',[])):c for c in cases}
        for f in cf.as_completed(futs):
            r=f.result(); res[r['name']]=r
            print("  %-24s %6.1fs  %s"%(r['name'],r['seconds'],
                  "ok" if not r['fails'] else "FAILED"))
    # cross-case relations, checked after everything has run
    for c in cases:
        r=res.get(c['name'])
        if r is None: continue
        other=c.get('same_deletion_set_as')
        if other and other in res:
            A=delset(pathlib.Path(res[other]['rundir'])/'m.damage')
            B=delset(pathlib.Path(r['rundir'])/'m.damage')
            if A!=B: r['fails'].append("deletion set differs from %s by %d element(s)"
                                       %(other,len(A^B)))
        other=c.get('identical_sta_as')
        if other and other in res:
            name='mixed.sta' if c['kind']=='mixed' else 'm.sta'
            A=(pathlib.Path(res[other]['rundir'])/name).read_bytes()
            B=(pathlib.Path(r['rundir'])/name).read_bytes()
            if A!=B: r['fails'].append("%s is not identical to %s"%(name,other))
        # equal_to: the same relation, stated as a TOLERANCE instead of as
        # byte identity, over whichever files the case names.  Byte identity
        # can prove a change is a no-op and cannot prove one is correct to
        # 1e-10, so an extraction that necessarily reorders a summation had
        # no way to be asserted at all.  tools/ccxdiff.py --selftest proves
        # this comparison goes red, and goes red for the right reason.
        specs=c.get('equal_to')
        if isinstance(specs,dict): specs=[specs]
        for spec in (specs or []):
          if spec['case'] in res:
              refdir=pathlib.Path(res[spec['case']]['rundir'])
              newdir=pathlib.Path(r['rundir'])
              rtol=spec.get('rtol',1.e-10); atol=spec.get('atol',1.e-12)
              exact=spec.get('exact',False)
              buf=io.StringIO()
              with contextlib.redirect_stdout(buf):
                  nbadf=sum(ccxdiff.compare_file(refdir,newdir,f,rtol,atol,exact)
                            for f in spec['files'])
              r['got']['equal_to']="%s at rtol=%g atol=%g%s"%(
                  spec['case'],rtol,atol," (exact)" if exact else "")
              if nbadf:
                  lines=[l.strip() for l in buf.getvalue().splitlines() if l.strip()]
                  head=[l for l in lines if ('DIFFERS' in l) or ('MISSING' in l)]
                  detail=[l for l in lines if l not in head]
                  r['fails'].append("not equal to %s: %s"%(spec['case'],"; ".join(head)))
                  # the first few offenders and a count.  A gate that prints
                  # four hundred lines of difference is a gate nobody reads.
                  for l in detail[:5]: r['fails'].append("   %s"%l)
                  if len(detail)>5:
                      r['fails'].append("   ... and %d more differing value(s)"
                                        %(len(detail)-5))
    # Which switches did any case actually put in force?  The [SWITCHES]
    # banner makes this measurable instead of assumed, and the number is
    # worth knowing: a switch no test ever sets is a switch whose behaviour
    # nobody is checking.
    cov=set()
    for r in res.values():
        try: txt=open(pathlib.Path(r['rundir'])/'run.log',errors='replace').read()
        except OSError: continue
        for m in re.finditer(r'^\[SWITCHES\]   (CCX_[A-Z0-9_]+) =',txt,re.M):
            cov.add(m.group(1))
    (outroot/'covered.txt').write_text("".join(x+"\n" for x in sorted(cov)))
    if a.record_coverage and not a.k:
        (HERE/'covered.txt').write_text("".join(x+"\n" for x in sorted(cov)))
        print("\nrecorded %d exercised switches in test/regress/covered.txt"%len(cov))
    print()
    nbad=0
    for c in cases:
        r=res[c['name']]
        print("%-24s %s"%(r['name'],r['what']))
        print("   %s"%"  ".join("%s=%s"%(k,v) for k,v in r['got'].items()))
        if r['fails']:
            nbad+=1
            for f in r['fails']: print("   FAIL %s"%f)
        else:
            print("   ok")
    if preflight_bad:
        print("\npreflight FAILED: the switch registry is stale")
    print("\nswitches put in force by these cases: %d"%len(cov))
    print("\n%d of %d case(s) failed%s"%(nbad,len(cases),
          ", plus the preflight" if preflight_bad else ""))
    if a.json:
        rec={'provenance':provenance(exe),
             'when':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
             'switches_exercised':sorted(cov),
             'failed':nbad,'preflight_failed':preflight_bad,
             'cases':[{'name':c['name'],'what':c['what'],
                       'expect':c['expect'],'got':res[c['name']]['got'],
                       'seconds':res[c['name']]['seconds'],
                       'fails':res[c['name']]['fails']} for c in cases]}
        pathlib.Path(a.json).write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
        print("wrote %s"%a.json)
    return nbad+preflight_bad

if __name__=='__main__':
    sys.exit(main())
