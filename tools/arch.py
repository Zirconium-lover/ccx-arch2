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
import argparse,json,pathlib,re,signal,sys

try: signal.signal(signal.SIGPIPE,signal.SIG_DFL)   # `| head` is not an error
except Exception: pass

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
                # k>=j, not k>j: a definition written entirely on ONE line -
                # `void damstats_init(damstats *d){ memset(d,0,sizeof(*d)); }'
                # - opens and closes on the same line, so k==j.  Requiring
                # k>j dropped every one of them, which mattered far beyond
                # the length table: module_symbols() is built from this, so
                # the layering check and the public-API count were blind to
                # them too.  Found by a partitioner that could not work out
                # who owned damstats_init().
                if k<n and k>=j:
                    out.append((m.group(1),i+1,j+1,k+1))
                    i=k
        i+=1
    return out

def typedef_names():
    """Every struct typedef CalculiX.h declares.

    Read rather than listed: a hand-kept list silently stops recognising the
    next module's type, the declaration scan then stops at its first use,
    and the locals count drops by ninety for no reason at all.  Measured the
    hard way."""
    h=(SRC/'CalculiX.h').read_text(errors='replace')
    return sorted(set(re.findall(r'^\}\s*([A-Za-z_]\w*)\s*;',h,re.M)))

DECLTYPES=(r'char|double|ITG|FILE|int|float|long|unsigned|size_t|'
           +"|".join(typedef_names()))

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
           'lsladder.c','crackcontrol.c','pathfollow.c','ccxopt.c','logview.c',
           'damstate.c','stiffcensus.c','loadcut.c','opcheck.c',
           'damdiag.c','erosion.c','dogleg.c','damcont.c','damstats.c',
           'slownewton.c','rescue.c','trial.c','loadctl.c','dammat.c',
           'census.c','monitor.c']
    return [SRC/n for n in names if (SRC/n).exists()]

# ---------------------------------------------------------------------------
# The second half of the measurement: THE INTERFACE.
#
# Everything above measures locality - how far apart the pieces of one
# decision sit.  That was the complaint as it was handed over, and it is only
# half of the thing.  A tree can have every decision in its own file and
# still be miserable to work in, and the numbers that say so are different
# ones:
#
#   4. the WIDEST PUBLIC SIGNATURE.  A function taking twenty-two arguments
#      does not have an interface; it has a copy of its caller's locals, and
#      every change to the caller is a change to it.  Measured here: the
#      worst signature, and how many are over the budget.  The three reasons
#      the lists got long are worth naming because they have three different
#      fixes: derived scalars passed by hand (mi[0], mi[1]+1 - thirty-three
#      hand-offs), iteration state that no object owns (ram, uam, qam, iit),
#      and model state the evaluator context never learned about (dambase,
#      dmcon);
#
#   5. LAYERS - which module may call which.  Stated and checked, a cycle
#      cannot appear by accident.  Unstated, the first one appears the day
#      somebody needs a number and reaches for it, and after that no piece
#      can be understood or tested alone.  The table below is not an
#      aspiration: it was derived from the call graph as it already is, and
#      it has zero violations on the tree that introduced it.  That is the
#      point - a rule adopted while it is free stays cheap; one adopted after
#      the cycles exist never gets adopted;
#
#   6. the REBUILD FAN-OUT: how many translation units recompile when one
#      extension interface changes.  This is the cost of an experiment, paid
#      every time anybody tries anything, and it is why a seven-thousand-line
#      shared header is not a question of taste;
#
#   7. how many self tests can run WITHOUT a finite-element deck.  A test
#      that needs a full analysis to run is a test that runs once a day, and
#      a test that runs once a day does not catch the mistake while the hand
#      is still on it.
#
# WHY THESE AND NOT OTHERS.  Each is a cost somebody pays on every change,
# and each has a direction that is not arguable: nobody wants a wider
# signature, a cycle, a bigger rebuild or a slower test.  Numbers that only
# go one way can be a ratchet; numbers that trade off cannot.

# Which layer each extension file sits in.  A module may call strictly DOWN
# and never sideways or up.
#
# Read the layers as answers to "what does this need in order to be true":
#   L0  the switch registry - depends on nothing, everything depends on it
#   L1  things that need nothing but the platform: instrumentation, and
#       plain statements about what the deck declared
#   L2  the model and its topology: what IS, before anybody has an opinion
#   L3  measurement of the model - observers, and the shared question of who
#       is driving the load factor.  These take their inputs const and cannot
#       change an answer
#   L4  mechanisms: single decisions taken about the model
#   L5  strategies that COMPOSE mechanisms - rescue ladders, continuation
#
# loadctl sits at L3 rather than among the mechanisms deliberately.  It does
# not decide anything; it answers "who owns the load factor this increment",
# which is the question every mechanism's arming block asks first.  Shared
# state read by many is not a mechanism, and putting it at L4 is what would
# force the sideways edges.
LAYERS={'ccxopt':0,
        'census':1,'dammat':1,'logview':1,'monitor':1,
        'damstate':2,'topology':2,'trial':2,
        'damdiag':3,'damstats':3,'loadctl':3,'opcheck':3,'topodiag':3,
        'converge':4,'crackcontrol':4,'damcont':4,'erosion':4,'globalize':4,
        'loadcut':4,'lsladder':4,
        'dogleg':5,'pathfollow':5,'rescue':5,'slownewton':5}

# The budget a public signature has to fit in.  Six is not a round number
# chosen to flatter: with the evaluator context, the iteration state and the
# derived-scalar accessors in place, the widest honest signature in this
# solver is an object, a context, a destination and a couple of scalars.
# Anything past that is a parameter list standing in for a missing object.
PARAM_BUDGET=6

def module_symbols():
    """name -> module, for every function an extension file DEFINES at file
    scope and does not mark static.

    Derived rather than listed, for the same reason typedef_names() is: a
    prefix table has to be maintained, and on the day a module gains a
    function that does not begin with the module's own name, every edge that
    function creates silently stops being counted.  Statics are excluded
    because they cannot be called across a file boundary - counting them
    would invent an edge out of two files that happen to name a helper the
    same thing."""
    own={}
    for p in ext_modules():
        if p.name=='nonlingeo.c': continue
        txt=p.read_text(errors='replace').split('\n')
        for name,first,_,_ in functions(p):
            if txt[first-1].lstrip().startswith('static'): continue
            own.setdefault(name,p.stem)
    return own

def module_deps(own):
    """module -> set of other modules it calls into."""
    deps={}
    for p in ext_modules():
        if p.name=='nonlingeo.c': continue
        src=p.read_text(errors='replace')
        src=re.sub(r'/\*.*?\*/','',src,flags=re.S)
        src=re.sub(r'"(\\.|[^"\\])*"','""',src)
        d=set()
        for call in set(re.findall(r'\b([A-Za-z_]\w*)\s*\(',src)):
            o=own.get(call)
            if o and o!=p.stem: d.add(o)
        deps[p.stem]=d
    return deps

def layer_violations(deps):
    """Every call that does not go strictly down.  Unknown modules are
    skipped rather than guessed: a file nobody has placed yet is a gap in
    the table, and reporting it as a violation would teach people to ignore
    the check."""
    bad=[]
    for m,ds in sorted(deps.items()):
        if m not in LAYERS: continue
        for d in sorted(ds):
            if d in LAYERS and LAYERS[d]>=LAYERS[m]:
                bad.append((m,LAYERS[m],d,LAYERS[d]))
    return bad

def public_api(own):
    """(name, params) for every extension function DECLARED in the shared
    header - i.e. every one whose signature the whole tree can see, and
    therefore every one that costs a full rebuild to change."""
    h=(SRC/'CalculiX.h').read_text(errors='replace')
    h=re.sub(r'/\*.*?\*/','',h,flags=re.S)
    out=[]
    for name,args in re.findall(
            r'\n(?:[A-Za-z_][\w \t*]*?)\b([A-Za-z_]\w*)\s*\(([^;{}]*?)\)\s*;',h,re.S):
        if name not in own: continue
        a=args.strip()
        out.append((name,0 if a in ('','void') else a.count(',')+1))
    return sorted(set(out),key=lambda x:(-x[1],x[0]))

def header_fanout():
    """How many .c files recompile when the shared header changes."""
    n=0
    for p in sorted(SRC.glob('*.c')):
        if re.search(r'#\s*include\s*"CalculiX\.h"',p.read_text(errors='replace')):
            n+=1
    return n,len(list(SRC.glob('*.c')))

def selftests():
    """(how many self tests exist, how many can be run without a deck).

    The second number is the one that matters.  A test reachable only from
    inside nonlingeo() costs a full analysis to run, and a test that costs
    minutes is a test that runs once a day.

    Counted from the TABLE in selftest.c rather than by matching function
    names.  Name matching was the first version and it was wrong twice over:
    it counted selftest_run_all() and selftest_gate() - the machinery, not
    tests - and it would have missed any test whose name did not contain the
    word.  The table is what both callers actually iterate, so it is what
    there are."""
    src=(SRC/'selftest.c')
    if not src.exists(): return 0,0
    txt=src.read_text(errors='replace')
    m=re.search(r'SELFTESTS\[\]\s*=\s*\{(.*?)\n\};',txt,re.S)
    total=len(re.findall(r'\{\s*"[^"]+"\s*,',m.group(1))) if m else 0
    # Deck-free means a runner exists that can reach them without a solve.
    runner=(SRC/'selftest_main.c').exists()
    return total,(total if runner else 0)

def measure():
    ng=SRC/'nonlingeo.c'
    txt=ng.read_text(errors='replace')
    fns=functions(ng)
    big=max(fns,key=lambda f:f[3]-f[1])
    loc=locals_of(ng,big[2],big[3])
    # The clusters are DERIVED, not listed: a hand-kept taxonomy flatters
    # itself by leaving out what it forgot.  Any prefix that four or more
    # locals share is a piece of state somebody already named; the map below
    # only supplies English for the ones that have been looked at.
    what={'damage_ct_':'the bounded continuation',
          'pf_':'dissipation path following',
          'damage_dl_':'the trust-region dogleg',
          'damage_diss_':'dissipation control',
          'damage_corr_':'the recovery corridor',
          'damage_fd_':'the operator finite-difference probe',
          'damage_slow_':'the Newton iteration budget',
          'damage_release_':'the release probe',
          'damage_bt_':'transactional backtracking',
          'damage_ray_':'the residual ray',
          'td_':'topology diagnostics',
          'damage_de13_':'terminal deletion',
          'damage_wall_':'the wall probes',
          'damage_evt_':'the event census',
          'damage_fracture_':'the termination test',
          'damage_stab_':'stabilisation',
          'damage_aba_':'the A-B-A purity test',
          'damage_linesearch_':'the damage line search',
          'damage_unsym_':'the asymmetric tangent',
          'damage_path_':'path control'}
    import collections as _c
    cnt=_c.Counter()
    for n in loc:
        m2=re.match(r'^(damage_[a-z0-9]+_|pf_|td_)',n)
        if m2: cnt[m2.group(1)]+=1
    clusters=[{'prefix':p,'what':what.get(p,'-'),'locals':k}
              for p,k in cnt.most_common() if k>=4]
    named=sum(c['locals'] for c in clusters)
    fork=sum(1 for n in loc if n.startswith(('damage_','pf_','td_')))
    clusters.append({'prefix':'(smaller clusters)','what':'fewer than four locals each',
                     'locals':fork-named})
    m={'nonlingeo_lines':len(txt.split('\n')),
       'longest_function':big[0],
       'longest_function_lines':big[3]-big[1]+1,
       'locals_in_longest':len(loc),
       'fork_locals':fork,
       'stock_locals':len(loc)-fork,
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

    own=module_symbols()
    deps=module_deps(own)
    api=public_api(own)
    fan,ncfile=header_fanout()
    ntest,ndeckless=selftests()
    m['layers']={k:sorted(v) for k,v in deps.items()}
    m['layer_violations']=[{'from':a,'from_layer':b,'to':c,'to_layer':d}
                           for a,b,c,d in layer_violations(deps)]
    m['unplaced_modules']=sorted(k for k in deps if k not in LAYERS)
    m['public_api']=len(api)
    m['widest_signature']=api[0][1] if api else 0
    m['widest_signature_name']=api[0][0] if api else '-'
    m['over_param_budget']=sum(1 for _,n in api if n>PARAM_BUDGET)
    m['worst_signatures']=[{'name':n,'params':k} for n,k in api[:10] if k>PARAM_BUDGET]
    m['header_fanout']=fan
    m['c_files']=ncfile
    m['selftests']=ntest
    m['selftests_deckless']=ndeckless
    return m

def table(m):
    o=[]
    o.append("THE FILE")
    o.append("  %-46s %d"%("src/nonlingeo.c, lines",m['nonlingeo_lines']))
    o.append("  %-46s %d"%("files the extension is spread over",m['extension_files']))
    o.append("")
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
    o.append("")
    o.append("THE INTERFACE  (what a change costs, rather than where it sits)")
    o.append("  %-46s %d"%("extension functions in the shared header",m['public_api']))
    o.append("  %-46s %d  (%s)"%("widest public signature, arguments",
                                 m['widest_signature'],m['widest_signature_name']))
    o.append("  %-46s %d"%("signatures over the %d-argument budget"%PARAM_BUDGET,
                           m['over_param_budget']))
    o.append("  %-46s %d of %d"%("files recompiled by a header change",
                                 m['header_fanout'],m['c_files']))
    o.append("  %-46s %d of %d"%("self tests runnable without a deck",
                                 m['selftests_deckless'],m['selftests']))
    if m['worst_signatures']:
        o.append("")
        o.append("  the widest, which are the missing objects named:")
        for w in m['worst_signatures']:
            o.append("    %-26s %2d"%(w['name'],w['params']))
    o.append("")
    o.append("LAYERS  (a module may call strictly DOWN, never sideways or up)")
    for L in sorted(set(LAYERS.values())):
        o.append("  L%d  %s"%(L,"  ".join(sorted(k for k,v in LAYERS.items() if v==L))))
    if m['unplaced_modules']:
        o.append("  not placed in the table: %s"%", ".join(m['unplaced_modules']))
    if m['layer_violations']:
        o.append("")
        for v in m['layer_violations']:
            o.append("  VIOLATION  %s(L%d) -> %s(L%d)"
                     %(v['from'],v['from_layer'],v['to'],v['to_layer']))
    else:
        o.append("  violations: 0")
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
             'public_api':m['public_api'],
             'widest_signature':m['widest_signature'],
             'over_param_budget':m['over_param_budget'],
             'header_fanout':m['header_fanout'],
             'layer_violations':len(m['layer_violations']),
             'selftests_deckless':-m['selftests_deckless'],
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
                  'fork_locals','public_api','widest_signature',
                  'over_param_budget','header_fanout'):
            if k in b and m[k]>b[k]: bad.append("%s: %d, budget %d"%(k,m[k],b[k]))
        # Layering is not a ratchet, it is a floor: the table was adopted on a
        # tree with zero violations, so any violation at all is new.
        if len(m['layer_violations'])>b.get('layer_violations',0):
            for v in m['layer_violations']:
                bad.append("layer violation %s(L%d) -> %s(L%d)"
                           %(v['from'],v['from_layer'],v['to'],v['to_layer']))
        # Stored negated, so the same "must not go up" rule makes deck-free
        # tests a number that must not go DOWN.
        if -m['selftests_deckless']>b.get('selftests_deckless',0):
            bad.append("self tests runnable without a deck: %d, budget %d"
                       %(m['selftests_deckless'],-b['selftests_deckless']))
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
