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

# A struct is written two ways in these headers: the anonymous
# `typedef struct{...}name;' and, for the contexts that have to be
# forward-declared, the tagged `struct name{...};'.  EVERY scan for type
# names has to know both.  Knowing only the first is a bug that has now
# appeared TWICE - once in struct_fields(), where it made the field set
# empty, and once in duplicated_params(), where tagging erosion_batch
# dropped it out of the known-type set and its parameter began counting as
# duplication.  Written once so there is one place left to be wrong.
def typedefs_in(text):
    return (set(re.findall(r'^\}\s*([A-Za-z_]\w*)\s*;',text,re.M))
           |set(re.findall(r'^struct\s+([A-Za-z_]\w*)\s*\{',text,re.M)))

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
    """Every struct typedef the tree declares, from EVERY header.

    It read CalculiX.h alone, which was right until the extension's
    declarations moved into ccxfork.h.  After that it returned two names,
    DECLTYPES stopped recognising `trialctx nlgt;' as a declaration, and the
    prologue scan in locals_of() - which stops at the first statement that
    does not begin with a type - would stop early at whichever module handle
    came first.  A hand-kept list was rejected here for the same reason; a
    hand-kept FILE is no better."""
    out=set()
    for h in sorted(SRC.glob('*.h')):
        out|=typedefs_in(h.read_text(errors='replace'))
    return sorted(out)


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
           'slownewton.c','rescue.c','trial.c','loadctl.c','dammat.c','release.c',
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
        'damstate':2,'nlstate':2,'topology':2,'trial':2,
        'damdiag':3,'damstats':3,'loadctl':3,'opcheck':3,'release':3,'topodiag':3,
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

def fork_headers():
    """The headers that declare this extension's interface.

    Found by asking which headers declare a function some extension file
    defines, rather than by keeping a list: the block is being split module
    by module, and a list would go stale on the commit after the one that
    wrote it."""
    own=module_symbols()
    out=[]
    for h in sorted(SRC.glob('*.h')):
        if h.name=='CalculiX.h': continue
        txt=re.sub(r'/\*.*?\*/','',h.read_text(errors='replace'),flags=re.S)
        if any(re.search(r'\b'+re.escape(n)+r'\s*\(',txt) for n in own):
            out.append(h)
    return out

def public_api(own):
    """(name, params) for every extension function DECLARED in a header -
    i.e. every one whose signature is visible outside its own file, and
    therefore every one that costs a rebuild to change."""
    out=[]
    for h in [SRC/'CalculiX.h']+fork_headers():
        txt=re.sub(r'/\*.*?\*/','',h.read_text(errors='replace'),flags=re.S)
        for name,args in re.findall(
                r'\n(?:[A-Za-z_][\w \t*]*?)\b([A-Za-z_]\w*)\s*\(([^;{}]*?)\)\s*;',txt,re.S):
            if name not in own: continue
            a=args.strip()
            out.append((name,0 if a in ('','void') else a.count(',')+1))
    return sorted(set(out),key=lambda x:(-x[1],x[0]))

# The scalars every caller recomputes from the model by hand.  mt is mi[1]+1,
# mi0 is mi[0], mi2 is mi[2], nstate is *nstate_.  Derived, so they are not
# state anybody owns - they are an accessor somebody has not written.
DERIVED={'mt','mi0','mi2','nstate','ndmat_','ntmat_','ndmat','ntmat'}

# The assembled operator and how to solve it.  These travel together through
# four functions and are not held by any object; they are a context nobody
# has built yet.
LINSYS={'ad','au','adb','aub','icol','irow','jq','nzs','isolver',
        'symmetryflag','inputformat','sigma','nrhs','nzs3','neq','ndof'}

def struct_fields(name):
    """The field names of one typedef'd struct in the extension headers.

    The closing `}name;' is found first and the opening searched BACKWARDS
    from it.  Searching forwards from the first `typedef struct' matches the
    wrong struct entirely - a non-greedy span from the first one runs through
    every struct in between, which reported trialctx as holding 367 fields
    and owning half of nlstate's.  Measured the hard way, on a number that
    then went into a plan."""
    for h in [SRC/'CalculiX.h']+fork_headers():
        t=h.read_text(errors='replace')
        # Two spellings, because the extension uses both: the anonymous
        # `typedef struct{...}name;' and, for the contexts that have to be
        # forward-declared, the tagged `struct name{...};'.  Handling only
        # the first was WRONG AND SILENT: tagging trialctx made this return
        # the empty set, the context-field set collapsed to the derived
        # scalars, and the duplicated-parameter count fell from 263 to 73
        # without a line of solver code changing.  A metric that improves
        # when you rename a struct is not measuring anything.
        start=t.find('struct %s{'%name)
        if start>=0:
            end=t.find('\n};',start)
        else:
            end=t.find('}%s;'%name)
            if end<0: continue
            start=t.rfind('typedef struct',0,end)
        if start<0 or end<0: continue
        body=re.sub(r'/\*.*?\*/','',t[start:end],flags=re.S)
        # One name per DECLARATOR, not per statement.  `double *ram,*ram1,
        # *ram2;' declares three fields and a regex ending at the semicolon
        # sees one - which read nlstate as holding seven of its nine.
        # trialctx happens to put one field per line, so the bug was
        # invisible there and would have stayed invisible.
        out=set()
        for stmt in body.split(';'):
            stmt=stmt.split('{')[-1]
            for part in stmt.split(','):
                m=re.search(r'([A-Za-z_]\w*)\s*(?:\[[^\]]*\])*\s*$',part)
                if m: out.add(m.group(1))
        return out
    return set()

def context_fields():
    """Everything some context already holds, or trivially could."""
    out=set()
    for n in ('trialctx','nlstate'): out|=struct_fields(n)
    return out|DERIVED|LINSYS

def selftest_reachable(own):
    """Public functions a self test calls directly.

    These are the ones that must KEEP taking plain arrays.  A self test
    builds a synthetic four-element mesh and three materials; it cannot
    build a 180-field trialctx bound to the locals of a running solver, and
    it should not have to.  Converting such a function to take the context
    would not tidy it - it would delete the only check on it.

    It is not an accident that the widest three functions in this tree are
    all in this set.  A function that takes plain arrays is a function you
    can test; a function that takes the context is one you can only run."""
    out=set()
    for p in sorted(SRC.glob('*.c')):
        txt=p.read_text(errors='replace')
        for name,first,body,last in functions(p):
            if 'selftest' not in name and not name.endswith('_legacycheck'):
                continue
            seg="\n".join(txt.split('\n')[first-1:last])
            for c in set(re.findall(r'\b([a-z_]\w*)\s*\(',seg)):
                if c in own and c!=name: out.add(c)
    return out

def duplicated_params(own):
    """How many parameters name something a context already holds.

    THE number this work is judged by, and the reason the obvious metric -
    how wide is the widest signature - is the wrong one.  A function that
    composes seven mechanisms and names all seven is honest at seven
    arguments; forcing it under six would mean inventing an eighth object
    whose only purpose is to hide the seven.  Width is a symptom.  The
    disease is a caller taking apart a context that already exists and
    passing the pieces in one at a time - which is exactly what the call to
    slownewton_allow() was doing with its own slownewton object.

    A parameter whose TYPE is one of the extension's structs is not counted:
    passing an object is the cure, not the disease."""
    ctx=context_fields()
    keep=selftest_reachable(own)
    # Two different diseases hide in one number, and they have two different
    # cures.  A parameter naming a trialctx or nlstate field is a caller
    # taking apart an object that EXISTS; a LINSYS name is a parameter of the
    # context this tree decided not to build (see the plan: a separate object
    # would split the matrix between two contexts, so `how to solve it' is to
    # become a parameter of the pre_solve extension point instead).  Only the
    # first is fixable today, and the second is the evidence for when to fix
    # it, so the report says which is which.
    split={'context':0,'derived':0,'linsys':0}
    own_ctx=struct_fields('trialctx')|struct_fields('nlstate')
    types=set()
    for h in [SRC/'CalculiX.h']+fork_headers():
        types|=typedefs_in(h.read_text(errors='replace'))
    n=0; worst=[]
    for h in [SRC/'CalculiX.h']+fork_headers():
        txt=re.sub(r'/\*.*?\*/','',h.read_text(errors='replace'),flags=re.S)
        for name,args in re.findall(
                r'\n(?:[A-Za-z_][\w \t*]*?)\b([A-Za-z_]\w*)\s*\(([^;{}]*?)\)\s*;',txt,re.S):
            if name not in own: continue
            if name in keep: continue     # must stay callable from a test
            k=0
            for a in args.split(','):
                a=a.strip()
                ty=re.match(r'(?:const\s+)?(\w+)',a)
                if ty and ty.group(1) in types: continue
                g=re.findall(r'(\w+)\s*(?:\[\s*\])?$',a)
                if g and g[0] in ctx:
                    k+=1
                    split['context' if g[0] in own_ctx else
                          'derived' if g[0] in DERIVED else 'linsys']+=1
            if k: worst.append((k,name))
            n+=k
    worst.sort(reverse=True)
    return n,worst,len(keep),split

def includers():
    """header name -> the .c files that reach it, following headers.

    Transitive, because that is what make rebuilds: ccxfork.h includes
    ccxopt.h, so touching ccxopt.h costs everything that reaches ccxfork.h
    too.  Counting only direct includes would have flattered the split that
    introduced this."""
    inc={}
    for f in list(SRC.glob('*.c'))+list(SRC.glob('*.h')):
        inc[f.name]=set(re.findall(r'#\s*include\s*"([^"]+)"',
                                   f.read_text(errors='replace')))
    def reaches(start,target,seen=None):
        if seen is None: seen=set()
        if start in seen: return False
        seen.add(start)
        for h in inc.get(start,()):
            if h==target or reaches(h,target,seen): return True
        return False
    out={}
    for h in [SRC/'CalculiX.h']+fork_headers():
        out[h.name]=sum(1 for c in SRC.glob('*.c') if reaches(c.name,h.name))
    return out

def header_fanout():
    """The cost of changing an extension INTERFACE: how many .c files
    recompile when the worst-placed extension header is touched.

    Not CalculiX.h's fan-out any more.  That number is now about stock
    CalculiX, which this work is not changing and whose 5,500 lines are not
    edited in the course of it; what matters is what an extension change
    costs, and that is the widest of the headers the extension owns."""
    inc=includers()
    fork={k:v for k,v in inc.items() if k!='CalculiX.h'}
    worst=max(fork.values()) if fork else 0
    return worst,len(list(SRC.glob('*.c')))

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
    m['header_cost']=includers()
    m['c_files']=ncfile
    m['selftests']=ntest
    m['selftests_deckless']=ndeckless
    ndup,worstdup,nkeep,dupsplit=duplicated_params(own)
    m['selftest_reachable']=nkeep
    m['duplicated_params']=ndup
    m['dup_split']=dupsplit
    m['worst_duplicators']=[{'name':x[1],'params':x[0]} for x in worstdup[:10]]
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
    o.append("  %-46s %d"%("parameters a context already holds",
                           m['duplicated_params']))
    s=m['dup_split']
    o.append("  %-46s %d"%("  a context that EXISTS: trialctx, nlstate",
                           s['context']))
    o.append("  %-46s %d"%("  a derived scalar: mt, mi0, nstate",s['derived']))
    o.append("  %-46s %d"%("  the linear system, which has no object yet",
                           s['linsys']))
    o.append("  %-46s %d"%("  functions exempt: a self test calls them",
                           m['selftest_reachable']))
    o.append("  %-46s %d of %d"%("files recompiled by an interface change",
                                 m['header_fanout'],m['c_files']))
    for h,n in sorted(m['header_cost'].items(),key=lambda kv:-kv[1]):
        o.append("      %-42s %d"%(h,n))
    o.append("  %-46s %d of %d"%("self tests runnable without a deck",
                                 m['selftests_deckless'],m['selftests']))
    if m['worst_signatures']:
        o.append("")
        o.append("  the widest, which are the missing objects named:")
        for w in m['worst_signatures']:
            o.append("    %-26s %2d"%(w['name'],w['params']))
    if m['worst_duplicators']:
        o.append("")
        o.append("  taking apart a context that already exists:")
        for w in m['worst_duplicators']:
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

# ---------------------------------------------------------------------------
# THE SELF TEST, and why this file of all files needed one.
#
# arch.py decides whether a commit is acceptable: --check is the ratchet.  It
# had no self test, and in one session it broke SILENTLY twice.
#
#   * struct_fields() knew only the anonymous `typedef struct{...}name;'
#     spelling.  Giving trialctx a tag made it return the empty set, the
#     context-field set collapsed to the derived scalars, and the headline
#     count fell from 263 to 73 WITHOUT A LINE OF SOLVER CODE CHANGING.  A
#     metric that improves when you rename a struct is measuring nothing;
#
#   * it counted one name per STATEMENT, so `double *ram,*ram1,*ram2;' read
#     as one field and nlstate came out as 7 of its 9.  trialctx happens to
#     put one field per line, so the bug was invisible there and would have
#     stayed invisible.
#
# Both were caught by disbelieving a number, which is not a method.  This is
# the method: synthetic sources in a temporary directory, where the right
# answer is known because it was written down first.
#
# The fixtures use REAL module names (trial.c, erosion.c, ...) because
# ext_modules() is an explicit list - a synthetic `foo.c' would be invisible
# to the very code under test.

def _chk(what,got,want,bad):
    ok=(got==want)
    if not ok: bad.append(what)
    print("[ARCH]   %-56s %-18s %s"
          %(what,"%s (want %s)"%(got,want),"PASS" if ok else "FAIL"))
    return ok

def selftest():
    """Every measurement, against sources whose answer is known."""
    import tempfile,shutil
    global SRC
    real=SRC
    bad=[]
    d=pathlib.Path(tempfile.mkdtemp(prefix='archtest'))
    try:
        SRC=d
        # ---- struct_fields: both spellings, multi-declarator, nesting ----
        (d/'CalculiX.h').write_text("#define ITG int\n")
        (d/'trial.h').write_text("""
typedef struct{
  double **co;
  ITG    **nk,**ne;
  double *a,*b,*c;
}trialctx;
struct nlstate{
  ITG    *iit;
  double *ram,*ram1,*ram2;
};
typedef struct{ ITG x; }other;
void trial_results(const trialctx *mdl);
void trial_check(const trialctx *mdl);
""")
        (d/'trial.c').write_text(
          "void trial_results(const trialctx *mdl){}\n"
          "void trial_check(const trialctx *mdl){}\n")
        f=struct_fields('trialctx')
        _chk("struct_fields: anonymous typedef spelling",
             sorted(f),['a','b','c','co','ne','nk'],bad)
        _chk("struct_fields: multi-declarator counts each name",
             len(struct_fields('trialctx')&{'a','b','c'}),3,bad)
        _chk("struct_fields: tagged `struct name{...};' spelling",
             sorted(struct_fields('nlstate')),
             ['iit','ram','ram1','ram2'],bad)
        _chk("struct_fields: a different struct does not leak in",
             'x' in struct_fields('trialctx'),False,bad)

        # ---- typedef_names and locals_of --------------------------------
        # This case exists because the bug happened.  typedef_names() read
        # CalculiX.h alone; when the extension's declarations moved to
        # ccxfork.h it returned two names, DECLTYPES stopped recognising
        # `trialctx nlgt;' as a declaration, and locals_of() - which ends the
        # prologue at the first statement that is not a declaration - stopped
        # at whichever module handle came first.  nonlingeo()'s locals read
        # 579 when they were 618, for four commits, and nothing said so.
        _chk("typedef_names: a type declared outside CalculiX.h is found",
             'trialctx' in typedef_names(),True,bad)
        (d/'nonlingeo.c').write_text(
          "void nonlingeo(void){\n"
          "  double aa=0.,bb=0.;\n"
          "  trialctx handle;\n"          # only found if the typedef is known
          "  ITG cc=0,dd=0;\n"
          "  aa=1.;\n"
          "}\n")
        import importlib
        global DECLTYPES
        DECLTYPES=(r'char|double|ITG|FILE|int|float|long|unsigned|size_t|'
                   +"|".join(typedef_names()))
        fns=functions(d/'nonlingeo.c')
        big=[f for f in fns if f[0]=='nonlingeo'][0]
        _chk("locals_of: the prologue does not stop at a module handle",
             sorted(locals_of(d/'nonlingeo.c',big[2],big[3])),
             ['aa','bb','cc','dd','handle'],bad)

        # ---- layering: one edge upward, and only one ---------------------
        (d/'erosion.c').write_text(
          '#include "topology.h"\nvoid erosion_mark(void){ topo_selftest(); }\n')
        (d/'topology.h').write_text("ITG topo_selftest(void);\nvoid topo_up(void);\n")
        (d/'topology.c').write_text(
          'ITG topo_selftest(void){ return 0; }\n'
          'void topo_up(void){ erosion_mark(); }\n')     # L2 -> L4, upward
        own=module_symbols()
        v=layer_violations(module_deps(own))
        _chk("layer_violations: an upward edge is reported",
             [(a,c) for a,_,c,_ in v],[('topology','erosion')],bad)

        # ---- duplicated params: field yes, object no, selftest exempt ----
        (d/'damdiag.h').write_text("""
void damdiag_a(const double *co,ITG nk,double z);
void damdiag_b(const trialctx *co,double z);
void damdiag_c(const double *co,ITG nk);
void damdiag_d(const nlstate *co,double z);
void damdiag_e(ITG nk,ITG mt,double *ad);
ITG  damdiag_selftest(void);
""")
        (d/'damdiag.c').write_text(
          "void damdiag_a(const double *co,ITG nk,double z){}\n"
          "void damdiag_b(const trialctx *co,double z){}\n"
          "void damdiag_c(const double *co,ITG nk){}\n"
          "void damdiag_d(const nlstate *co,double z){}\n"
          "void damdiag_e(ITG nk,ITG mt,double *ad){}\n"
          "ITG damdiag_selftest(void){ damdiag_c(0,0); return 0; }\n")
        own=module_symbols()
        n,worst,nkeep,sp=duplicated_params(own)
        _chk("duplicated_params: a context field counts",
             dict(( (w,k) for k,w in worst )).get('damdiag_a'),2,bad)
        # damdiag_b's parameter is deliberately NAMED `co', which IS a
        # trialctx field.  Only its TYPE distinguishes it, so a version that
        # forgot the type check would count it - with any other name the
        # case could not fail and would be a test in appearance only.
        _chk("duplicated_params: an object-typed parameter does not",
             'damdiag_b' in [w for _,w in worst],False,bad)
        # ...and the same when the type uses the TAGGED spelling.  Added
        # because it was MISSING and the ratchet, not this test, is what
        # caught the consequence: tagging erosion_batch dropped it out of
        # the known-type set, and release_arm's `b' - an erosion_batch, not
        # trialctx's right-hand side - began counting as duplication.
        _chk("duplicated_params: a TAGGED object type is recognised too",
             'damdiag_d' in [w for _,w in worst],False,bad)
        _chk("duplicated_params: a self test's callee is exempt",
             'damdiag_c' in [w for _,w in worst],False,bad)
        # The split, because one number was hiding two diseases with two
        # different cures.  damdiag_e takes one of each: `nk' is a field of
        # a context that EXISTS, `mt' is a derived scalar, `ad' belongs to
        # the linear system, which has no object.  A classifier that put
        # everything in the first bucket - the easy mistake, since that is
        # the bucket the total used to be - leaves the other two at zero.
        _chk("dup_split: a field of an existing context",
             sp['context'],3,bad)
        _chk("dup_split: a derived scalar is not that",
             sp['derived'],1,bad)
        _chk("dup_split: a linear-system name is not that either",
             sp['linsys'],1,bad)

        # ---- includers: transitive, because that is what make rebuilds ---
        (d/'logview.h').write_text("void logview_report(double s);\n")
        (d/'logview.c').write_text('#include "logview.h"\n'
                                   "void logview_report(double s){}\n")
        (d/'damstats.c').write_text('#include "damdiag.h"\n')
        (d/'damdiag.h').write_text('#include "logview.h"\n'
                                   +(d/'damdiag.h').read_text())
        inc=includers()
        # logview.c reaches it directly; damstats.c only through damdiag.h.
        # Counting direct includes alone would say 1, and would have
        # flattered exactly the umbrella-header arrangement that the real
        # measurement rejected.
        _chk("includers: counts a file that reaches a header indirectly",
             inc.get('logview.h'),2,bad)

        # ---- selftests(): from the table, not from function names --------
        (d/'selftest.c').write_text("""
const selftest_entry SELFTESTS[]={
  {"one",one_selftest},
  {"two",two_selftest},
};
ITG selftest_run_all(ITG v){return 0;}
void selftest_gate(void){}
/* A helper that a by-name count would wrongly call a third test.  Without
   it the fixture could not tell the two counting rules apart, which is the
   shape of a test that only looks like one. */
static ITG two_coeff_selftest(void){ return 0; }
""")
        (d/'selftest_main.c').write_text("int main(void){return 0;}\n")
        tot,deck=selftests()
        _chk("selftests: counted from the table, not by name",tot,2,bad)
        _chk("selftests: a helper matching the name is not counted",
             tot,2,bad)
    finally:
        SRC=real
        shutil.rmtree(d,ignore_errors=True)
    print("\n[ARCH] self test: %d failure(s) -- %s"
          %(len(bad),"FAILED" if bad else "PASSED"))
    return len(bad)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--check',action='store_true',
                    help='exit nonzero if any measured number is worse than the budget')
    ap.add_argument('--record',action='store_true',help='rewrite the budget')
    ap.add_argument('--selftest',action='store_true',
                    help='check every measurement against synthetic sources')
    ap.add_argument('--json',default=None)
    a=ap.parse_args()
    if a.selftest: return selftest()
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
             'duplicated_params':m['duplicated_params'],
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
                  'over_param_budget','header_fanout','duplicated_params'):
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
