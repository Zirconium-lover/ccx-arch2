#!/usr/bin/env python3
"""Is the dissipation path follower doing its job?

    CCX_EXE=... tools/pathfollow.py --variant plain --tau 1e-4
    tools/pathfollow.py --selftest

CCX_PATHFOLLOW implements the Gutierrez (2004) / Verhoosel et al. (2009)
dissipation constraint, and it prints a line per accepted increment.  What
it has never had is a VERDICT: the lines say what happened, nobody says
whether what happened is the method working.

Five properties, each checked separately, because they fail separately:

  ARMED     the switch reached the binary and the domain gate let it arm.
            "not armed: outside the verified domain" is a silent no-op
            otherwise - the run simply proceeds under ordinary control.

  ENGAGED   the measured dissipation reached 0.2*tau and the constraint
            took over.  A run that arms and never engages has paid for
            nothing.

  ENFORCED  once engaged, the achieved dG equals the requested tau.  This
            is the one that decides whether the METHOD works: the whole
            point of a constraint is that it is satisfied.  A Newton
            iteration whose Jacobian row is not the derivative of its
            residual row converges to something else, and the symptom is
            exactly dG != tau with the run otherwise looking healthy.
            src/pathfollow.c says the PREVIOUS implementation in this tree
            had that defect; this checks the current one.

  FREE      lambda decreased at least once.  A dissipation follower exists
            to walk a branch where the load parameter falls.  If lambda is
            monotone the follower has done nothing a step-time control
            could not, whatever else it reports.

  STABLE    refusals do not grow without bound.  Every refusal is an
            increment where the constraint could not be applied and the
            ordinary path was taken instead.
"""
import argparse,os,pathlib,re,shutil,subprocess,sys

ROOT=pathlib.Path(__file__).resolve().parent.parent

ARM=re.compile(r'^\[PATHFOLLOW\] armed: tau=([-\d.eE+]+)')
NOARM=re.compile(r'^\[PATHFOLLOW\] not armed: (.*)')
ERR=re.compile(r'^\[PATHFOLLOW\] \*ERROR: (.*)')
ENG=re.compile(r'^\[PATHFOLLOW\] engaged at inc=(\d+) lambda=([-\d.eE+]+)')
ACC=re.compile(r'^\[PATHFOLLOW\] inc=(\d+) ACCEPTED lambda=([-\d.eE+]+) '
               r'dG=([-\d.eE+]+) P=([-\d.eE+]+) ff=([-\d.eE+]+) '
               r'engaged=(\d+) refusals=(\d+)')

def parse(text):
    r={'tau':None,'noarm':None,'error':None,'engaged_at':None,'rows':[]}
    for l in text.splitlines():
        m=ARM.match(l)
        if m: r['tau']=float(m.group(1)); continue
        m=NOARM.match(l)
        if m: r['noarm']=m.group(1); continue
        m=ERR.match(l)
        if m: r['error']=m.group(1); continue
        m=ENG.match(l)
        if m: r['engaged_at']=int(m.group(1)); continue
        m=ACC.match(l)
        if m:
            r['rows'].append({'inc':int(m.group(1)),'lam':float(m.group(2)),
                              'dG':float(m.group(3)),'P':float(m.group(4)),
                              'ff':float(m.group(5)),'eng':int(m.group(6)),
                              'ref':int(m.group(7))})
    return r

def verdict(r,rtol=0.05,minfrac=0.8):
    """Returns (list of (name, ok, detail), nbad)."""
    out=[]
    armed=(r['tau'] is not None) and (r['error'] is None)
    out.append(("ARMED",armed,
                r['error'] or r['noarm'] or
                ("tau=%.6e"%r['tau'] if r['tau'] is not None else
                 "no [PATHFOLLOW] arming line at all")))
    eng=r['engaged_at'] is not None
    out.append(("ENGAGED",eng,
                ("at increment %d"%r['engaged_at']) if eng else
                "the measured dissipation never reached 0.2*tau"))
    rows=[x for x in r['rows'] if x['eng']==1]
    if rows and r['tau']:
        err=[abs(x['dG']-r['tau'])/r['tau'] for x in rows]
        good=sum(1 for e in err if e<=rtol)
        frac=good/len(err)
        worst=max(err)
        out.append(("ENFORCED",frac>=minfrac,
                    "%d of %d engaged increments within %.0f%% of tau "
                    "(worst %.1f%%)"%(good,len(err),100*rtol,100*worst)))
    else:
        out.append(("ENFORCED",False,"no engaged increment to judge"))
    if len(rows)>1:
        drops=sum(1 for a,b in zip(rows,rows[1:]) if b['lam']<a['lam']-1e-12)
        out.append(("FREE",drops>0,
                    "lambda decreased on %d of %d engaged increments"
                    %(drops,len(rows)-1)))
    else:
        out.append(("FREE",False,"fewer than two engaged increments"))
    if r['rows']:
        ref=[x['ref'] for x in r['rows']]
        grew=ref[-1]-ref[0]
        out.append(("STABLE",grew<=max(5,0.1*len(r['rows'])),
                    "refusals went %d -> %d over %d accepted increments"
                    %(ref[0],ref[-1],len(r['rows']))))
    else:
        out.append(("STABLE",False,"no accepted increment reported"))
    return out,sum(1 for _,ok,_ in out if not ok)

def report(name,r,rtol,minfrac):
    print("\n%s"%name)
    rows,nbad=verdict(r,rtol,minfrac)
    for n,ok,detail in rows:
        print("  %-4s %-9s %s"%("ok" if ok else "FAIL",n,detail))
    if r['rows']:
        print("  %d accepted increment(s); lambda %.6f -> %.6f"
              %(len(r['rows']),r['rows'][0]['lam'],r['rows'][-1]['lam']))
    return nbad

def _fake(tau=1.e-4,n=12,enforce=True,free=True,engage=True,refusals=0):
    L=[]
    if tau is not None: L.append("[PATHFOLLOW] armed: tau=%.6e per increment, x"%tau)
    if engage: L.append("[PATHFOLLOW] engaged at inc=3 lambda=0.50000000: y")
    lam=0.5
    for i in range(n):
        lam=lam-0.001 if free else lam+0.001
        dG=tau if enforce else tau*3.0
        L.append("[PATHFOLLOW] inc=%d ACCEPTED lambda=%.8f dG=%.6e P=1.0 "
                 "ff=1.0 engaged=%d refusals=%d"
                 %(10+i,lam,dG,1 if engage else 0,refusals*i))
    return "\n".join(L)

def selftest():
    bad=0
    def chk(what,text,want_bad):
        nonlocal bad
        rows,nb=verdict(parse(text))
        ok=(nb>0)==want_bad
        if not ok: bad+=1
        names=[n for n,o,_ in rows if not o]
        print("  %-4s %-52s %s"%("ok" if ok else "FAIL",what,
              ("red: "+",".join(names)) if nb else "green"))
    chk("a healthy follower passes",_fake(),False)
    chk("never arming is caught","[PATHFOLLOW] not armed: outside the "
        "verified domain (ncont=3)",True)
    chk("arming and never engaging is caught",_fake(engage=False),True)
    chk("dG three times tau is caught - the constraint is NOT enforced",
        _fake(enforce=False),True)
    chk("a monotone lambda is caught - the follower did nothing",
        _fake(free=False),True)
    chk("runaway refusals are caught",_fake(refusals=3),True)
    chk("an empty log is caught","",True)
    print("\n[PATHFOLLOW DIAG SELFTEST] %s"
          %("PASSED" if bad==0 else "%d check(s) FAILED"%bad))
    return bad

def main():
    ap=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--variant',default='plain')
    ap.add_argument('--tau',default='1.e-4')
    ap.add_argument('--threads',type=int,default=2)
    ap.add_argument('--rtol',type=float,default=0.05)
    ap.add_argument('--minfrac',type=float,default=0.8)
    ap.add_argument('-o',default=None)
    ap.add_argument('--log',default=None,help="judge an existing run.log")
    a=ap.parse_args()
    if a.selftest: return selftest()
    if a.log:
        return report(a.log,parse(open(a.log,errors='replace').read()),
                      a.rtol,a.minfrac)
    if selftest()!=0:
        print("the diagnostic's own self test failed; refusing to report")
        return 2
    exe=os.environ.get('CCX_EXE')
    if not exe or not os.access(exe,os.X_OK):
        sys.exit("set CCX_EXE to a PARDISO-enabled ccx_2.23 binary")
    d=pathlib.Path(a.o or ('/tmp/pf-'+a.variant)).resolve()
    shutil.rmtree(d,ignore_errors=True); d.mkdir(parents=True)
    env=dict(os.environ)
    env['OMP_NUM_THREADS']=str(a.threads); env['MKL_NUM_THREADS']=str(a.threads)
    env['FAST_VARIANT']=a.variant; env['CCX_EXE']=exe
    r=subprocess.run("%s/test/fast/run_fast.sh %s CCX_PATHFOLLOW=%s"
                     %(ROOT,d,a.tau),shell=True,env=env,
                     stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    print("rc=%d, run in %s"%(r.returncode,d))
    return report("variant %s, tau=%s"%(a.variant,a.tau),
                  parse((d/'run.log').read_text(errors='replace')),
                  a.rtol,a.minfrac)

if __name__=='__main__': sys.exit(main())
