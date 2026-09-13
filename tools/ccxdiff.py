#!/usr/bin/env python3
"""Compare two runs field by field, to a stated tolerance.

    tools/ccxdiff.py REF NEW [-f m.sta] [--rtol 1e-10] [--atol 1e-12] [--exact]

The gate this replaces the strictest half of (test/regress/run.py) pins a
case by scalars and by BYTE identity of m.sta and m.damage.  Byte identity is
a strong check and a terrible foundation for refactoring: it can prove a
change is a no-op and it cannot express "the answer is the same to 1e-10", so
the first decomposition that reorders a summation fails it for the wrong
reason.  That constraint is part of what produced the code this project
exists to take apart - see NEXT_TASK.md item 2.

So: byte identity stays available and stays the strictest setting (--exact),
and below it there is a real comparison.

WHAT IS COMPARED, and why these files

  m.sta      the accepted increment trajectory: step, increment, attempt,
             iterations, and the three times
  m.damage   the committed deletion history.  The element SET is compared
             element by element (a stronger statement than a count, which is
             why the gate already did it); the deletion TIMES are compared
             numerically
  m.cvg      residual and correction per iteration, per attempt
  m.dat      whatever the deck asked to be printed - here the grip reaction,
             which is the physically meaningful scalar
  m.frd      the nodal fields, block by block (DISP, STRESS, ...), per node
             and per component

TOLERANCE SEMANTICS, following Code_Aster's CRITERE='RELATIF'/'ABSOLU' pair
and MOOSE's CSVDiff rel_err/abs_zero: a pair (a,b) agrees when

    |a-b| <= atol      OR      |a-b| <= rtol * max(|a|,|b|)

atol is what keeps a field whose true value is machine zero from failing on
relative error alone.  Labels - node and element numbers, step and increment
counters, material ids - are NEVER compared with a tolerance; they must match
exactly, because a tolerance on an identity is a way of not noticing that two
runs did different things.

The tolerance actually used is printed with every comparison, passing or
failing, so a green line cannot quietly mean "compared with 1e-3".
"""
import argparse,json,math,os,pathlib,re,sys

# ---------------------------------------------------------------- readers

def read_sta(path):
    """SUMMARY OF JOB INFORMATION: one row per accepted increment."""
    rows=[]
    for line in open(path,errors='replace'):
        f=line.split()
        if len(f)!=7: continue
        try: rows.append(((int(f[0]),int(f[1]),int(f[2])),
                          {'ITRS':int(f[3])},
                          {'TOT TIME':float(f[4]),'STEP TIME':float(f[5]),
                           'INC TIME':float(f[6])}))
        except ValueError: continue
    return {'kind':'rows','rows':rows,'keyname':'(step,inc,att)'}

def read_cvg(path):
    """SUMMARY OF CONVERGENCE INFORMATION: one row per Newton iteration."""
    rows=[]
    for line in open(path,errors='replace'):
        f=line.split()
        if len(f)<5: continue
        try: step,inc,att,it=int(f[0]),int(f[1]),int(f[2]),int(f[3])
        except ValueError: continue
        vals={}
        for i,name in enumerate(('CONT.EL.','RESID.FORCE','CORR.DISP',
                                 'RESID.FLUX','CORR.TEMP')):
            if 4+i<len(f):
                try: vals[name]=float(f[4+i])
                except ValueError: pass
        rows.append(((step,inc,att,it),{},vals))
    return {'kind':'rows','rows':rows,'keyname':'(step,inc,att,iter)'}

def read_damage(path):
    """element step increment step_time total_time material damage ip batch"""
    rows=[]
    for line in open(path,errors='replace'):
        if line.startswith('#') or not line.strip(): continue
        f=line.split()
        if len(f)<9: continue
        rows.append(((int(f[0]),),
                     {'step':int(f[1]),'increment':int(f[2]),
                      'material':int(f[5]),'critical_ip':int(f[7]),
                      'batch':int(f[8])},
                     {'step_time':float(f[3]),'total_time':float(f[4]),
                      'damage':float(f[6])}))
    return {'kind':'rows','rows':rows,'keyname':'(element)'}

_NUM=re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$')
_INT=re.compile(r'^\d+$')

def read_dat(path):
    """Whatever the deck printed.  A section is a heading line plus the rows
    under it; the heading carries the time, so two runs that printed at
    different times are a structural difference, not a numeric one.

    Two shapes of row occur and they must not be confused.  A nodal row
    begins with an integer node number and the rest are values.  A SUMMED
    row - "total force (fx,fy,fz) for set ... and time ..." - has no label
    at all, only the components, and the first of them is the grip reaction,
    which is the single most interesting number the fast decks print.  Until
    2026-09-11 both shapes went through the nodal branch, so the reaction
    became part of the KEY: a run whose grip force changed was reported as
    one record vanishing and another appearing, never as a value differing,
    and no tolerance was ever applied to it.  Headless rows are therefore
    keyed by their position within the section and compared in full."""
    rows=[]; head=None; ordinal={}
    for line in open(path,errors='replace'):
        s=line.rstrip()
        if not s.strip(): continue
        f=s.split()
        if _NUM.match(f[0]) and all(_NUM.match(x) for x in f[1:]):
            if head is None: continue
            if _INT.match(f[0]) and len(f)>1:
                key=(head,f[0]); vals=f[1:]
            else:
                n=ordinal.get(head,0); ordinal[head]=n+1
                key=(head,'#%d'%n); vals=f
            rows.append((key,{},
                         {'c%d'%i:float(x.replace('D','E').replace('d','e'))
                          for i,x in enumerate(vals)}))
        else:
            head=' '.join(f)
    return {'kind':'rows','rows':rows,'keyname':'(section,row)'}

def read_frd(path):
    """The nodal result blocks.  100CL gives the block's time and step; -4
    names the block and -5 its components; -1 records carry one node each."""
    rows=[]; block=None; comps=[]; time=None; step=None
    for line in open(path,errors='replace'):
        if line.startswith('  100CL'):
            f=line.split()
            try: time,step=float(f[2]),int(f[3])
            except (ValueError,IndexError): time,step=None,None
            block=None; comps=[]
        elif line.startswith(' -4'):
            block=line.split()[1]; comps=[]
        elif line.startswith(' -5'):
            if block is not None: comps.append(line.split()[1])
        elif line.startswith(' -1') and block is not None:
            # fixed width: ' -1' then a 10-char node field then 12-char values
            node=line[3:13].strip()
            vals=[]; i=13
            while i+12<=len(line.rstrip('\n')):
                t=line[i:i+12].strip()
                if not t: break
                try: vals.append(float(t))
                except ValueError: break
                i+=12
            names=comps if len(comps)==len(vals) else ['c%d'%k for k in range(len(vals))]
            rows.append(((block,step,node),{},dict(zip(names,vals))))
        elif line.startswith(' -3'):
            block=None
    return {'kind':'rows','rows':rows,'keyname':'(block,step,node)'}

READERS={'m.sta':read_sta,'mixed.sta':read_sta,'close.sta':read_sta,
         'm.cvg':read_cvg,'mixed.cvg':read_cvg,'close.cvg':read_cvg,
         'm.damage':read_damage,
         'm.dat':read_dat,'mixed.dat':read_dat,'close.dat':read_dat,
         'm.frd':read_frd,'mixed.frd':read_frd,'close.frd':read_frd}

def reader_for(name):
    if name in READERS: return READERS[name]
    for suffix,fn in (('.sta',read_sta),('.cvg',read_cvg),('.damage',read_damage),
                      ('.dat',read_dat),('.frd',read_frd)):
        if name.endswith(suffix): return fn
    return None

# ------------------------------------------------------------- comparison

def agree(a,b,rtol,atol):
    d=abs(a-b)
    if d<=atol: return True,d,0.
    m=max(abs(a),abs(b))
    r=d/m if m>0. else float('inf')
    return (r<=rtol),d,r

def compare(ref,new,rtol,atol,maxreport=5):
    """Returns (findings, worst) where findings is a list of strings and worst
    is the largest achieved (abs,rel) error over everything that DID agree -
    the number a reader needs in order to choose the next tolerance."""
    out=[]; wabs=0.; wrel=0.; wwhere=None
    # A record's natural key is not always unique - an frd block repeats its
    # step number at every output time, and a re-attempted increment repeats
    # its counters - so the key carries the occurrence index too.  That keeps
    # a missing record LOCAL: without it, one dropped row shifts every
    # comparison after it and the report becomes unreadable.
    def index(rows):
        d={}; seen={}
        for k,lab,val in rows:
            n=seen.get(k,0); seen[k]=n+1
            d[(k,n)]=(lab,val)
        return d
    A=index(ref['rows']); B=index(new['rows'])
    show=lambda kk: kk[0] if kk[1]==0 else "%s#%d"%(kk[0],kk[1]+1)
    only_ref=sorted(set(A)-set(B),key=repr); only_new=sorted(set(B)-set(A),key=repr)
    if only_ref:
        out.append("%d record(s) only in the reference, first %s"
                   %(len(only_ref),[show(k) for k in only_ref[:maxreport]]))
    if only_new:
        out.append("%d record(s) only in the new run, first %s"
                   %(len(only_new),[show(k) for k in only_new[:maxreport]]))
    nnum=0
    for k in sorted(set(A)&set(B),key=repr):
        la,va=A[k]; lb,vb=B[k]
        for name in sorted(set(la)|set(lb)):
            if la.get(name)!=lb.get(name):
                out.append("%s %s: label %s differs exactly: %r vs %r"
                           %(ref['keyname'],show(k),name,la.get(name),lb.get(name)))
        for name in sorted(set(va)|set(vb)):
            if name not in va or name not in vb:
                out.append("%s %s: %s present in only one run"
                           %(ref['keyname'],show(k),name))
                continue
            nnum+=1
            ok,d,r=agree(va[name],vb[name],rtol,atol)
            if not ok:
                if len(out)<maxreport*4:
                    out.append("%s %s: %s %.17g vs %.17g  (abs %.3g, rel %.3g)"
                               %(ref['keyname'],show(k),name,va[name],vb[name],d,r))
            elif (r>wrel) or (r==wrel and d>wabs):
                wabs,wrel,wwhere=d,r,"%s %s"%(show(k),name)
    return out,(wabs,wrel,wwhere,nnum)

_TOK=re.compile(r'[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?')

def resolution(path,limit=400000):
    """How finely can this file's own digits resolve a relative difference?

    A tolerance tighter than the file's printed precision is not a strict
    comparison, it is a comparison of rounding.  Getting "equal" back would
    be the "not wrong but uninformative" answer that costs a run.

    Computed per token, not from a digit count: what matters is the QUANTUM
    of the last printed digit relative to the value itself, and that depends
    on the format AND on the magnitude.  `0.9639E+00` quantises at 1e-4 on a
    value of 0.96, so it resolves 1e-4; `0.100000E+01` quantises at 1e-5 on
    a value of 1.0, so it resolves 1e-5.  A digit count alone gets both of
    those wrong, in opposite directions.

    Integers are skipped - a node number is exact, not quantised - and the
    median is reported rather than an extreme, because one badly scaled
    value should not speak for the file."""
    rels=[]
    txt=open(path,errors='replace').read(limit)
    for m in _TOK.finditer(txt):
        t=m.group(0).replace('D','E').replace('d','e')
        if ('.' not in t) and ('e' not in t) and ('E' not in t): continue
        mant,_,exp=t.replace('E','e').partition('e')
        if '.' not in mant: continue
        frac=len(mant.split('.')[1])
        try:
            v=abs(float(t)); e=int(exp) if exp else 0
        except ValueError:
            continue
        if v==0.: continue
        rels.append((10.**(e-frac))/v)
    if not rels: return 0.
    rels.sort()
    return rels[len(rels)//2]

_CLOCK=re.compile(rb'^(    1UTIME +)\d\d:\d\d:\d\d *$',re.M)
_DATE=re.compile(rb'^(    1UDATE +)\d{1,2}\.[A-Za-z]+\.\d{4} *$',re.M)

def _strip_clock(data):
    """The ONE exception to byte identity, and it is not a tolerance.

    CalculiX stamps the wall clock into the .frd header as a 1UTIME
    record and the calendar date as a 1UDATE one.  The first version of
    this excused only the clock, and the omission surfaced the honest way:
    two gate runs either side of midnight reported all twelve cases as
    differing, and the difference was 11.september against 12.september.  Two runs of the same binary on the same deck therefore
    never compare byte-identical, which quietly made --exact useless on
    the file that carries the nodal results - the extraction of
    damrank1.f on 2026-09-11 was reported as changing m.frd in all nine
    cases, and the change was the clock.

    Only that record is blanked, only in its exact fixed-width form, and
    the report says so whenever it fired.  Anything else in the header -
    the version, the node count, a block name - is still compared byte
    for byte.
    """
    out,n=_CLOCK.subn(rb'\g<1>HH:MM:SS',data)
    out,m=_DATE.subn(rb'\g<1>DD.MONTH.YYYY',out)
    return out,n+m

def compare_file(refdir,newdir,fname,rtol,atol,exact,quiet=False):
    a=pathlib.Path(refdir)/fname; b=pathlib.Path(newdir)/fname
    if not a.exists() and not b.exists():
        print("  %-12s absent from both runs"%fname); return 0
    for p in (a,b):
        if not p.exists():
            print("  %-12s MISSING: %s"%(fname,p)); return 1
    if exact:
        ra,na=_strip_clock(a.read_bytes()); rb,nb=_strip_clock(b.read_bytes())
        same=ra==rb
        print("  %-12s %s  (exact: byte for byte%s)"
              %(fname,"identical" if same else "DIFFERS",
                "" if not (na or nb) else ", wall clock excepted"))
        return 0 if same else 1
    rd=reader_for(fname)
    if rd is None:
        same=a.read_bytes()==b.read_bytes()
        print("  %-12s %s  (no reader for this name, fell back to byte identity)"
              %(fname,"identical" if same else "DIFFERS"))
        return 0 if same else 1
    res=resolution(a)
    note=("  [this file's digits resolve about %.0e relative, so rtol=%g is "
          "below what they can answer]"%(res,rtol)) if rtol<res else ""
    findings,(wabs,wrel,wwhere,nnum)=compare(rd(a),rd(b),rtol,atol)
    if findings:
        print("  %-12s DIFFERS  (rtol=%g atol=%g, %d numbers compared)%s"
              %(fname,rtol,atol,nnum,note))
        for f in findings[:20]: print("       %s"%f)
        if len(findings)>20: print("       ... and %d more"%(len(findings)-20))
        return 1
    if not quiet:
        print("  %-12s equal to rtol=%g atol=%g  (%d numbers; worst rel %.3g%s)%s"
              %(fname,rtol,atol,nnum,wrel,(" at "+wwhere) if wwhere else "",note))
    return 0

# ---------------------------------------------------------------- self test

_STA_HEAD="SUMMARY OF JOB INFORMATION\n  STEP      INC     ATT  ITRS     TOT TIME     STEP TIME      INC TIME\n"

def _sta(rows):
    t=_STA_HEAD
    for st,inc,att,it,a,b,c in rows:
        t+="%6d %10d %5d %5d  %12.6E  %12.6E  %12.6E\n"%(st,inc,att,it,a,b,c)
    return t

def _dam(rows):
    t="# CalculiX committed damage deletion history v2\n"
    for el,st,inc,a,b,mat,dmg,ip,batch in rows:
        t+="%d %d %d %.15e %.15e %d %.15e %d %d\n"%(el,st,inc,a,b,mat,dmg,ip,batch)
    return t

def _mk(d,name,text):
    d.mkdir(parents=True,exist_ok=True); (d/name).write_text(text); return d

def selftest():
    """Demonstrate the comparison going red, and going red for the right
    reason.  A suite that cannot fail is decoration; so is a tolerance whose
    dial nobody has turned in both directions."""
    import io,contextlib,tempfile
    base=[(1,i,1,3,1.e-2*i,1.e-2*i,1.e-2) for i in range(1,21)]
    dbase=[(100+i,1,i,1.e-2*i,1.e-2*i,2,1.0,1,i) for i in range(1,6)]
    bad=0
    def run(what,ref,new,expect_bad,**kw):
        nonlocal bad
        buf=io.StringIO()
        with contextlib.redirect_stdout(buf):
            n=sum(compare_file(ref,new,f,kw.get('rtol',1.e-10),kw.get('atol',1.e-12),
                               kw.get('exact',False)) for f in kw.get('files',['m.sta']))
        ok=(n>0)==expect_bad
        if not ok: bad+=1
        print("  %-4s %-62s %s"%("ok" if ok else "FAIL",what,
              "red" if n>0 else "green"))
        if not ok:
            for l in buf.getvalue().splitlines(): print("         %s"%l)
    with tempfile.TemporaryDirectory() as tmp:
        T=pathlib.Path(tmp)
        R=_mk(T/'ref','m.sta',_sta(base)); _mk(T/'ref','m.damage',_dam(dbase))

        _mk(T/'same','m.sta',_sta(base)); _mk(T/'same','m.damage',_dam(dbase))
        run("identical runs agree",R,T/'same',False,files=['m.sta','m.damage'])

        # A relative perturbation of 1e-8 on one number: the dial must turn
        # both ways.  It is done on m.damage and not on m.sta because m.sta
        # carries seven significant digits and CANNOT represent a 1e-8
        # change - which is itself worth knowing, and is why the tool reports
        # each file resolution alongside the tolerance.
        pert=[r if i!=2 else (r[0],r[1],r[2],r[3]*(1+1.e-8),r[4],r[5],r[6],r[7],r[8])
              for i,r in enumerate(dbase)]
        _mk(T/'p8','m.damage',_dam(pert))
        run("a 1e-8 relative change is caught at rtol=1e-10",R,T/'p8',True,
            files=['m.damage'])
        run("the same change passes at rtol=1e-6",R,T/'p8',False,
            files=['m.damage'],rtol=1.e-6)
        _mk(T/'p8s','m.sta',_sta(base))
        run("a 1e-8 change m.sta cannot represent is reported as agreement",
            R,T/'p8s',False)

        # a near-zero absolute change: atol is what keeps machine zero honest
        z=[(1,1,1,3,0.0,0.0,0.0)]; zp=[(1,1,1,3,1.e-13,0.0,0.0)]
        _mk(T/'z','m.sta',_sta(z)); _mk(T/'zp','m.sta',_sta(zp))
        run("1e-13 against zero passes at atol=1e-12",T/'z',T/'zp',False)
        run("the same passes nothing at atol=1e-14",T/'z',T/'zp',True,atol=1.e-14)

        # a LABEL must never be forgiven by a tolerance
        dl=[(r if i!=2 else (r[0],r[1],r[2],r[3],r[4],99,r[6],r[7],r[8]))
            for i,r in enumerate(dbase)]
        _mk(T/'lab','m.damage',_dam(dl))
        run("a changed material id is caught at any tolerance",R,T/'lab',True,
            files=['m.damage'],rtol=1.,atol=1.e30)

        # The grip reaction: a SUMMED m.dat row has no node number.  If the
        # reader keys it by its first field the reaction is a label, so the
        # dial does not turn on it at all - which is the failure this pair
        # of checks demonstrates.  Red at a tight tolerance AND green at a
        # loose one together prove the number is compared as a number.
        def _dat(fx):
            return (" forces (fx,fy,fz) for set G and time  0.1000000E+01\n\n"
                    "       519 -1.463154E+00  0.000000E+00  0.000000E+00\n\n"
                    " total force (fx,fy,fz) for set G and time  0.1000000E+01\n\n"
                    "        %.6E  0.000000E+00  0.000000E+00\n"%fx)
        _mk(T/'dat','m.dat',_dat(61.22186)); _mk(T/'datp','m.dat',_dat(61.22265))
        run("a changed grip reaction is caught at rtol=1e-6",T/'dat',T/'datp',
            True,files=['m.dat'],rtol=1.e-6)
        run("the same is a VALUE difference, forgiven at rtol=1e-4",
            T/'dat',T/'datp',False,files=['m.dat'],rtol=1.e-4)

        # a different deletion SET, which is the check the gate already valued
        _mk(T/'set','m.damage',_dam(dbase[:-1]))
        run("a missing deletion record is caught",R,T/'set',True,files=['m.damage'])
        _mk(T/'set2','m.damage',_dam(dbase[:-1]+[(999,1,5,5.e-2,5.e-2,2,1.0,1,5)]))
        run("a substituted element id is caught",R,T/'set2',True,files=['m.damage'])

        # the wall clock in the .frd header, and ONLY the wall clock
        frd=(" "*4+"1UTIME              17:38:47"+" "*24+"\n"
             +"    1UVERSION  CalculiX 2.23\n")
        _mk(T/'clk','m.frd',frd)
        _mk(T/'clk2','m.frd',frd.replace("17:38:47","18:15:39"))
        run("a differing .frd wall clock is not a difference",
            T/'clk',T/'clk2',False,files=['m.frd'],exact=True)
        _mk(T/'clk3','m.frd',frd.replace("2.23","2.24"))
        run("anything else in the same header still is",
            T/'clk',T/'clk3',True,files=['m.frd'],exact=True)
        dated=("    1UDATE              11.september.2026"+" "*31+"\n"
               +"    1UVERSION  CalculiX 2.23\n")
        _mk(T/'dt','m.frd',dated)
        _mk(T/'dt2','m.frd',dated.replace("11.september","12.september"))
        run("a run either side of midnight is not a difference",
            T/'dt',T/'dt2',False,files=['m.frd'],exact=True)
        _mk(T/'dt3','m.frd',dated.replace("2026","2027").replace("1UDATE","1UDATX"))
        run("and a mangled date record IS still a difference",
            T/'dt',T/'dt3',True,files=['m.frd'],exact=True)

        # and --exact must be strictly stronger than any tolerance
        _mk(T/'ws','m.sta',_sta(base).replace('\n',' \n',1))
        run("trailing whitespace is invisible to the tolerant comparison",
            R,T/'ws',False)
        run("trailing whitespace IS visible to --exact",R,T/'ws',True,exact=True)
    print("\n[CCXDIFF SELFTEST] %s"%("PASSED" if bad==0 else
          "%d check(s) FAILED"%bad))
    return bad

def main():
    ap=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ref',nargs='?'); ap.add_argument('new',nargs='?')
    ap.add_argument('--selftest',action='store_true',
                    help='prove the comparison can go red, and for the right reason')
    ap.add_argument('-f','--file',action='append',default=[],
                    help='compare only these names (default: everything both runs have)')
    ap.add_argument('--rtol',type=float,default=1.e-10)
    ap.add_argument('--atol',type=float,default=1.e-12)
    ap.add_argument('--exact',action='store_true',
                    help='byte identity - the strictest setting, still available')
    a=ap.parse_args()
    if a.selftest: return selftest()
    if not a.ref or not a.new: ap.error("give two run directories, or --selftest")
    names=a.file or sorted({p.name for d in (a.ref,a.new)
                            for p in pathlib.Path(d).iterdir()
                            if reader_for(p.name)})
    print("reference %s\nnew       %s"%(a.ref,a.new))
    bad=sum(compare_file(a.ref,a.new,n,a.rtol,a.atol,a.exact) for n in names)
    print("\n%d of %d file(s) differ"%(bad,len(names)))
    return bad

if __name__=='__main__': sys.exit(main())
