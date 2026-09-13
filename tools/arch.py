#!/usr/bin/env python3
"""Measure the architecture, so "it is better now" can be checked instead of said.

    tools/arch.py                 # the tables
    tools/arch.py --check         # fail if anything got worse than the budget
    tools/arch.py --record        # rewrite the budget from this tree

WHAT IS BEING MEASURED AND WHY
------------------------------
The complaint this repository was handed is one sentence: to find out where a
decision is taken you have to read nonlingeo() whole.  That is a statement
about LOCALITY, and locality is measurable.  Three numbers say it:

  1. the length of the longest function, because a 16,000-line body is the
     shape the complaint describes;

  2. the number of LOCALS that function declares, because every local is a
     piece of state with no owner - it is reachable, and writable, from
     every one of those lines.  This is the number that actually decides
     whether a block can be moved out: a phase that touches 200 locals
     cannot become a function, whatever anybody wants;

  3. for each named subsystem, its SPAN - how many separate places in
     nonlingeo.c its code sits in, and how far apart the first and last are.
     A subsystem living in one file has a span of 0 sites and an answer to
     "where is this decided".  One with eleven sites 9,000 lines apart has
     no answer, which is the complaint restated as an integer.

The budget file makes them a ratchet.  A number that is allowed to drift
back up is not a measurement, it is a mood.
"""
import argparse,json,pathlib,re,sys

ROOT=pathlib.Path(__file__).resolve().parent.parent
SRC=ROOT/'src'
BUDGET=pathlib.Path(__file__).resolve().parent/'arch_budget.json'

# The subsystems whose code the complaint is about.  A tag is what the source
# itself prints, so this list is not a private taxonomy: grep finds exactly
# what the run's own log calls the thing.
TAGS=['DAMAGE CT','PATHFOLLOW','WALLDIAG','CRACKCTL','OPCHECK','LOADCUT',
      'TOPOLOGY','TOPODIAG','GLOBALIZE','DAMSTATE','CONVERGE','DAMAGE RAY',
      'DAMAGE TR','DAMAGE BT','DAMAGE CORR','DAMAGE ABA','DAMAGE EVT',
      'DAMAGE RELEASE','DAMAGE DE1.3','SWITCHES','DISSIPATION','CENSUS',
      'BATCHTRACE','LSLADDER','DAMAGE REG','FRACTURE']

def strip(line):
    """Enough comment/string removal to count braces without being fooled."""
    s=re.sub(r'/\*.*?\*/','',line)
    s=re.sub(r'//.*','',s)
    s=re.sub(r'"(\\.|[^"\\])*"','""',s)
    s=re.sub(r"'(\\.|[^'\\])*'","''",s)
    return s

def functions(path):
    """(name, first line, last line) for every function defined at file scope.

    A definition is a line starting in column 1 that carries a '(' and ends
    the parameter list with '{' within a few lines.  Crude, and it does not
    have to be better: the point is to find the long ones."""
    txt=path.read_text(errors='replace').split('\n')
    out=[];i=0;n=len(txt)
    while i<n:
        l=txt[i]
        m=re.match(r'^(?:static\s+|const\s+|unsigned\s+)*[A-Za-z_][\w \t*]*?([A-Za-z_]\w*)\s*\(',l)
        if m and not l.startswith('#') and 'typedef' not in l:
            # walk forward to the brace that opens the body
            j=i;depth=0;opened=False
            while j<min(n,i+120):
                s=strip(txt[j])
                if '{' in s: opened=True;break
                if ';' in s: break
                j+=1
            if opened:
                # now count to the matching close
                depth=0;k=j;started=False
                while k<n:
                    s=strip(txt[k])
                    for ch in s:
                        if ch=='{': depth+=1;started=True
                        elif ch=='}': depth-=1
                    if started and depth<=0: break
                    k+=1
                if k<n and k>j:
                    out.append((m.group(1),i+1,j+1,k+1))
                    i=k
        i+=1
    return out

DECLTYPES=(r'char|double|ITG|FILE|int|float|long|unsigned|size_t|'
           r'topo_txn|glob_census|converge|lsladder|crackcontrol_census|'
           r'topodiag_report|damstate|damcfg|damct|damdiag_ray')

def locals_of(path,first,last):
    """Every name the function declares.

    Statement-wise, not line-wise: declarations here run over several lines
    and carry brace initialisers, so a line-by-line rule mistakes `mass[2]=
    {0,0}` for the end of the prologue and stops counting after 138 of
    them.  The prologue ends at the first statement that does not begin
    with a type.  It counts NAMES, not bytes: a name is a piece of state
    with no owner."""
    txt="\n".join(path.read_text(errors='replace').split('\n')[first:last])
    txt=re.sub(r'/\*.*?\*/','',txt,flags=re.S)
    txt=re.sub(r'"(\\.|[^"\\])*"','""',txt)
    names=set();depth=0;cur=''
    for ch in txt:
        if ch in '([{': depth+=1
        elif ch in ')]}': depth-=1
        if ch==';' and depth==0:
            s=cur.strip();cur=''
            m=re.match(r'^(?:const\s+|static\s+|unsigned\s+|long\s+)*(?:%s)\b(.*)$'%DECLTYPES,
                       s,flags=re.S)
            if not m: break
            d=0;piece='';parts=[]
            for c in m.group(1):
                if c in '([{': d+=1
                elif c in ')]}': d-=1
                if c==',' and d==0: parts.append(piece);piece=''
                else: piece+=c
            parts.append(piece)
            for p in parts:
                p=p.split('=')[0].replace('*',' ')
                mm=re.search(r'([A-Za-z_]\w*)\s*(?:\[[^\]]*\])*\s*$',p)
                if mm: names.add(mm.group(1))
        else:
            cur+=ch
    return names

def sites(text,tag,gap=40):
    """Contiguous regions of one file that mention a tag.

    Two mentions closer than `gap` lines are one site: the question is how
    many PLACES a reader has to find, not how many lines mention it."""
    lines=text.split('\n')
    hit=[i for i,l in enumerate(lines,1) if '['+tag+']' in l or '['+tag+' ' in l]
    if not hit: return [],0
    reg=[[hit[0],hit[0]]]
    for h in hit[1:]:
        if h-reg[-1][1]<=gap: reg[-1][1]=h
        else: reg.append([h,h])
    return reg,hit[-1]-hit[0]

def ext_modules():
    """The files this fork added or rewrote around the nonlinear solve.

    Named explicitly rather than detected: a list that guesses would quietly
    start counting stock CalculiX and the number would stop meaning anything."""
    names=['nonlingeo.c','converge.c','topology.c','topodiag.c','globalize.c',
           'lsladder.c','crackcontrol.c','pathfollow.c','ccxopt.c',
           'damstate.c','stiffcensus.c','loadcut.c','fracture.c','mincut.c',
           'damdiag.c','erosion.c','dogleg.c','damcont.c','damstats.c',
           'slownewton.c','damcfg.c']
    return [SRC/n for n in names if (SRC/n).exists()]

def measure():
    ng=SRC/'nonlingeo.c'
    txt=ng.read_text(errors='replace')
    fns=functions(ng)
    big=max(fns,key=lambda f:f[3]-f[1])
    loc=locals_of(ng,big[2],big[3])
    pre=[('damage_ct_','the bounded continuation'),
         ('pf_','dissipation path following'),
         ('damage_diss_','dissipation control'),
         ('damage_corr_','the recovery corridor'),
         ('damage_slow_','the Newton iteration budget'),
         ('damage_release_','the release probe'),
         ('damage_bt_','transactional backtracking'),
         ('damage_ray_','the residual ray'),
         ('td_','topology diagnostics'),
         ('damage_de13_','terminal deletion'),
         ('damage_wall_','the wall probes'),
         ('damage_evt_','the event census'),
         ('damage_fracture_','the termination test'),
         ('damage_stab_','stabilisation'),
         ('damage_','the rest of the damage extension')]
    clusters=[]
    left=set(loc)
    for p,what in pre:
        got={n for n in left if n.startswith(p)}
        left-=got
        clusters.append({'prefix':p,'what':what,'locals':len(got)})
    m={'nonlingeo_lines':len(txt.split('\n')),
       'longest_function':big[0],
       'longest_function_lines':big[3]-big[1]+1,
       'locals_in_longest':len(loc),
       'fork_locals':sum(c['locals'] for c in clusters),
       'stock_locals':len(left),
       'clusters':clusters,
       'file_scope_functions':len(fns),
       'subsystems':{}}
    for t in TAGS:
        reg,span=sites(txt,t)
        owner=[p.name for p in ext_modules()
               if p.name!='nonlingeo.c' and ('['+t+']' in p.read_text(errors='replace')
                                             or '['+t+' ' in p.read_text(errors='replace'))]
        m['subsystems'][t]={'sites_in_nonlingeo':len(reg),'span_lines':span,
                            'files':owner}
    m['extension_files']=len(ext_modules())
    return m

def table(m):
    o=[]
    o.append("THE FUNCTION")
    o.append("  %-46s %s"%("longest function in src/nonlingeo.c",m['longest_function']+'()'))
    o.append("  %-46s %d"%("its length, lines",m['longest_function_lines']))
    o.append("  %-46s %d"%("locals it declares",m['locals_in_longest']))
    o.append("  %-46s %d"%("  of them stock CalculiX",m['stock_locals']))
    o.append("  %-46s %d"%("  of them this fork's",m['fork_locals']))
    o.append("  %-46s %d"%("functions at file scope in that file",m['file_scope_functions']))
    o.append("")
    o.append("STATE WITH NO OWNER  (locals of the longest function, by cluster)")
    for c in m['clusters']:
        if c['locals']:
            o.append("  %-22s %4d   %s"%(c['prefix']+'*',c['locals'],c['what']))
    o.append("")
    o.append("WHERE IS IT DECIDED  (sites = separate places in nonlingeo.c to read)")
    o.append("  %-18s %6s %8s  %s"%("subsystem","sites","span","file that owns it"))
    for t,s in sorted(m['subsystems'].items(),
                      key=lambda kv:(-kv[1]['sites_in_nonlingeo'],kv[0])):
        own=",".join(s['files']) if s['files'] else "-"
        o.append("  %-18s %6d %8d  %s"%(t,s['sites_in_nonlingeo'],s['span_lines'],own))
    return "\n".join(o)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--check',action='store_true',
                    help='exit nonzero if any measured number is worse than the budget')
    ap.add_argument('--record',action='store_true',help='rewrite the budget')
    ap.add_argument('--json',default=None)
    a=ap.parse_args()
    m=measure()
    if a.json: pathlib.Path(a.json).write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
    if a.record:
        BUDGET.write_text(json.dumps(
            {'nonlingeo_lines':m['nonlingeo_lines'],
             'longest_function_lines':m['longest_function_lines'],
             'locals_in_longest':m['locals_in_longest'],
             'fork_locals':m['fork_locals'],
             'sites':{t:s['sites_in_nonlingeo'] for t,s in m['subsystems'].items()}},
            indent=2,sort_keys=True)+"\n")
        print("recorded %s"%BUDGET)
        return 0
    print(table(m))
    if a.check:
        if not BUDGET.exists():
            print("\nno budget recorded; run --record"); return 1
        b=json.loads(BUDGET.read_text());bad=[]
        for k in ('nonlingeo_lines','longest_function_lines','locals_in_longest',
                  'fork_locals'):
            if m[k]>b[k]: bad.append("%s: %d, budget %d"%(k,m[k],b[k]))
        for t,n in b.get('sites',{}).items():
            got=m['subsystems'].get(t,{}).get('sites_in_nonlingeo',0)
            if got>n: bad.append("%s: %d sites, budget %d"%(t,got,n))
        print()
        if bad:
            for x in bad: print("OVER BUDGET  %s"%x)
            print("\n%d number(s) worse than tools/arch_budget.json"%len(bad))
            return 1
        print("[ARCH] every measured number is at or under tools/arch_budget.json")
    return 0

if __name__=='__main__':
    sys.exit(main())
