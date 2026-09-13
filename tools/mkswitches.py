#!/usr/bin/env python3
"""Regenerate the option registry from the sources AND from the declarations.

    tools/mkswitches.py            # rewrite src/ccxopt_list.h and docs/SWITCHES.md
    tools/mkswitches.py --check    # exit 1 if either is out of date

Two sources, and the difference between them is the point.

  SCRAPED: every CCX_* name the binary reads, found by scanning for getenv in
  the C and Fortran sources.  This is what ccxopt.c reports from, so a run
  states its own configuration and a misspelt name is caught rather than
  silently doing nothing - the failure that makes an A/B uninformative rather
  than wrong, and therefore the expensive kind.

  DECLARED: src/ccxopt_decl.h, hand written, giving each option a type, a
  default, a range, its legal spellings and one line of prose.  The generated
  documentation now comes FROM the declarations for anything that has one,
  instead of being scraped from the line that reads it.

An option that is scraped and not declared is reported as UNDECLARED, here
and at run time.  That list is the retirement queue: a switch with no
declaration, no test and no prose has no defenders.
"""
import argparse,os,re,sys,collections

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC=os.path.join(ROOT,'src')
HDR=os.path.join(SRC,'ccxopt_list.h')
DOC=os.path.join(ROOT,'docs','SWITCHES.md')
COV=os.path.join(ROOT,'test','regress','covered.txt')
GETENV=re.compile(r'getenv\s*\(\s*[\'"](CCX_[A-Z0-9_]+)[\'"]')

BANNER=re.compile(r'printf\s*\(\s*"((?:[^"\\\\]|\\\\.)*)"((?:\s*"(?:[^"\\\\]|\\\\.)*")*)')

def banner(txt,pos):
    """The banner this switch prints, taken verbatim from the source.

    Only a printf that appears BEFORE the next getenv of any CCX_ name is
    accepted, so a switch cannot borrow its neighbour's description - which
    is what a naive nearest-printf search does, and it produces confident
    nonsense.  A switch that shares a parse block with another therefore gets
    no description here, which is the honest answer."""
    nxt=GETENV.search(txt,pos)
    end=nxt.start() if nxt else min(len(txt),pos+2500)
    m=BANNER.search(txt,pos,end)
    if not m: return None
    s=m.group(1)+m.group(2)
    s=re.sub(r'"\s*"','',s)
    s=s.replace('\\n',' ').replace('\\"',"'")
    s=re.sub(r'%[-0-9.*]*(?:"\s*ITGFORMAT\s*"|[a-zA-Z])','<value>',s)
    s=s.replace('"','')
    s=re.sub(r'\s+',' ',s).strip()
    s=re.sub(r'(\s*<value>)+$','',s).strip()
    if s.startswith('*ERROR') or len(s)<30: return None
    return s

COMMENT=re.compile(r'/\*(.*?)\*/',re.S)

def preceding_comment(txt,pos):
    """The C comment block that ends just above this getenv, if it does.

    Falls back to this only when the switch prints no banner, and only when
    the comment is genuinely adjacent - at most two lines above - so a
    comment belonging to something else higher up is not attributed here."""
    line_start=txt.rfind('\n',0,pos)
    if line_start<0: return None
    head=txt[:line_start]
    best=None
    for m in COMMENT.finditer(head): best=m
    if best is None: return None
    between=head[best.end():]
    if between.count('\n')>2 or between.strip(): return None
    body=re.sub(r'^\s*\*+','',best.group(1),flags=re.M)
    body=re.sub(r'\s+',' ',body).strip()
    if len(body)<40: return None
    return body

def scan():
    where=collections.defaultdict(set); desc={}
    for dirpath,_,files in os.walk(SRC):
        for fn in sorted(files):
            if not fn.endswith(('.c','.f','.h')): continue
            if fn in ('ccxopt_list.h','ccxopt_decl.h'): continue
            txt=open(os.path.join(dirpath,fn),errors='replace').read()
            for m in GETENV.finditer(txt):
                where[m.group(1)].add(fn)
                if m.group(1) not in desc:
                    b=banner(txt,m.end())
                    if b is None and fn.endswith(('.c','.h')):
                        c=preceding_comment(txt,m.start())
                        if c: b='(from the comment above it) '+c
                    if b: desc[m.group(1)]=b
    return {k:sorted(v) for k,v in sorted(where.items())},desc

def docs_text():
    out=[]
    for dirpath,_,files in os.walk(ROOT):
        if '.git' in dirpath.split(os.sep): continue
        for fn in files:
            if fn.endswith('.md') and os.path.join(dirpath,fn)!=DOC:
                out.append(open(os.path.join(dirpath,fn),errors='replace').read())
    return '\n'.join(out)

def src_text():
    out=[]
    for dirpath,_,files in os.walk(SRC):
        for fn in files:
            if fn.endswith(('.c','.f','.h')) and fn not in ('ccxopt_list.h','ccxopt_decl.h'):
                out.append(open(os.path.join(dirpath,fn),errors='replace').read())
    return '\n'.join(out)

DECL=os.path.join(SRC,'ccxopt_decl.h')

def _c_strings(text):
    """Adjacent C string literals, concatenated as C does."""
    return "".join(m.group(1).replace('\\"','"')
                   for m in re.finditer(r'"((?:[^"\\]|\\.)*)"',text))

def declarations():
    """Read src/ccxopt_decl.h - the hand-written table - into a dict.

    A small parser and not a real one: the table's shape is fixed and lives
    in this repository, so the cost of getting it wrong is a build that fails
    loudly rather than documentation that is quietly false."""
    if not os.path.exists(DECL): return {}
    txt=open(DECL,errors='replace').read()
    txt=re.sub(r'/\*.*?\*/','',txt,flags=re.S)
    i=txt.find('ccxopt_decl_table[]={')
    if i<0: return {}
    body=txt[i+len('ccxopt_decl_table[]={'):]
    out={}
    for m in re.finditer(r'\{\s*("(?:[^"\\]|\\.)*"[\s\S]*?)\},',body):
        e=m.group(1)
        # split on commas that are outside string literals
        parts=[];cur='';instr=False;esc=False
        for ch in e:
            if esc: cur+=ch; esc=False; continue
            if ch=='\\': cur+=ch; esc=True; continue
            if ch=='"': instr=not instr; cur+=ch; continue
            if ch==',' and not instr: parts.append(cur); cur=''; continue
            cur+=ch
        parts.append(cur)
        # CCXOPT_UNBOUNDED is one macro standing for two fields; expand it so
        # the positions line up whether an entry spells the range out or not
        for k,q in enumerate(parts):
            if q.strip()=='CCXOPT_UNBOUNDED':
                parts[k:k+1]=['1.','0.']; break
        if len(parts)<7: continue
        name=_c_strings(parts[0])
        if not name.startswith('CCX_'): continue
        typ=parts[1].strip()
        dflt=_c_strings(parts[2])
        lo,hi=parts[3].strip(),parts[4].strip()
        if lo=='CCXOPT_UNBOUNDED': lo,hi='1.','0.'
        choices=_c_strings(parts[5]) if '"' in parts[5] else ''
        doc=_c_strings(parts[6])
        dep=_c_strings(parts[7]) if len(parts)>7 and '"' in parts[7] else ''
        try: rng='' if float(lo)>float(hi) else "[%s, %s]"%(lo.rstrip('.'),hi.rstrip('.'))
        except ValueError: rng=''
        out[name]={'type':typ.replace('CCXOPT_','').lower(),'default':dflt,
                   'range':rng,'choices':choices,'doc':doc,'deprecated':dep}
    return out

def header(sw):
    fort={k:any(f.endswith('.f') for f in v) for k,v in sw.items()}
    L=["/* GENERATED by tools/mkswitches.py - do not edit.",
       "",
       "   Every CCX_* name this binary reads.  ccxopt.c reports the ones",
       "   that are set, validates the ones that are declared, names anything",
       "   else in the environment that looks like one of ours, and says at",
       "   the end which were set and never read.",
       "",
       "   ccxopt_known_fortran marks a name whose only reads are in Fortran.",
       "   Those sites do not route through ccxopt_getenv yet, so the",
       "   set-but-never-read report says nothing about them rather than",
       "   claiming they were unread.",
       "",
       "   Regenerate with tools/mkswitches.py; tools/mkswitches.py --check",
       "   fails if this file is stale. */",
       "",
       "#define CCXOPT_KNOWN_COUNT %d"%len(sw),
       "",
       "static const char *const ccxopt_known_name[CCXOPT_KNOWN_COUNT]={"]
    for k in sw: L.append('  "%s",'%k)
    L+=["};",
        "",
        "static const char ccxopt_known_fortran[CCXOPT_KNOWN_COUNT]={"]
    for k in sw: L.append('  %d,'%(1 if fort[k] else 0))
    L+=["};",""]
    return "\n".join(L)

def coverage():
    """Switches that the three-minute gate actually puts in force.

    Recorded by test/regress/run.py --record-coverage from the [SWITCHES]
    banner of each case, so it is measured rather than declared."""
    if not os.path.exists(COV): return set()
    return {l.strip() for l in open(COV) if l.strip()}

def document(sw,dtxt,stxt,desc):
    cov=coverage(); dec=declarations()
    nd=[k for k in sw if k not in dtxt and k not in desc and k not in dec]
    undeclared=[k for k in sw if k not in dec]
    orphan=[k for k in undeclared if k not in cov and k not in dtxt and k not in desc]
    L=["# Options",
       "",
       "GENERATED by `tools/mkswitches.py` - do not edit by hand.",
       "",
       "Two tables, and the difference between them is the point.",
       "",
       "| | |",
       "|---|---|",
       "| names the binary reads | **%d** |"%len(sw),
       "| **declared** in `src/ccxopt_decl.h` - type, default, range, spellings, prose | **%d** |"%len(dec),
       "| undeclared, documented only by whatever the source says of them | %d |"%len(undeclared),
       "| **explained nowhere at all** | **%d** |"%len(nd),
       "| exercised by `test/regress/run.py` | %d |"%len([k for k in sw if k in cov]),
       "| **never set by any test in this tree** | **%d** |"%len([k for k in sw if k not in cov]),
       "| **no declaration, no test and no prose - the retirement queue** | **%d** |"%len(orphan),
       "",
       "Every run prints the ones that are set (`[SWITCHES]` at the top of any",
       "`run.log`), validates the declared ones against their type and range,",
       "names anything in the environment starting with `CCX_` that is not in",
       "this list, and reports at the end which options were set and never",
       "read (`[SWITCHES LEFT]`).",
       "",
       "## Declared",
       "",
       "These have an owner. The type, default, range and description below",
       "are read from the declaration, not from the line that reads it.",
       "",
       "| option | type | default | range / values | in the gate | what it is |",
       "|---|---|---|---|---|---|"]
    for k in sorted(dec):
        d=dec[k]
        rng=d['choices'].replace('|','\\|') if d['choices'] else d['range']
        L.append("| `%s` | %s | %s | %s | %s | %s |"
                 %(k,d['type'],d['default'] or '-',rng or '-',
                   "yes" if k in cov else "-",
                   (d['doc']+(" **DEPRECATED: "+d['deprecated']+"**"
                              if d['deprecated'] else "")).replace("|","\\|")))
    L+=["",
        "## Undeclared",
        "",
        "No type, no stated default, no range, and a description that is",
        "whatever the source happens to say: the banner printed by the block",
        "that reads it, or failing that the comment immediately above the",
        "line that reads it.  Nothing here is written by hand or inferred,",
        "which is the point - and also the limit.  Several are parsed in a",
        "shared block with one banner between them, so what you see may",
        "describe the mechanism rather than that one knob.  A blank means the",
        "code says nothing there at all.",
        "",
        "| option | read by | elsewhere | in the gate | what it says of itself |",
        "|---|---|---|---|---|"]
    for k in undeclared:
        L.append("| `%s` | %s | %s | %s | %s |"
                 %(k,", ".join("`%s`"%f for f in sw[k]),
                   "yes" if k in dtxt else "**no**",
                   "yes" if k in cov else "-",
                   desc.get(k,"").replace("|","\\|")))
    L+=["",
        "\"elsewhere\" means the name appears in some other markdown file in",
        "this tree.  It is a presence check, not a quality one.",
        "",
        "## The retirement queue",
        "",
        "%d option(s) have no declaration, no test that sets them and no"%len(orphan),
        "prose anywhere but the line that reads them.  A switch in this state",
        "has no defenders: retiring it means making its behaviour the default",
        "or deleting it, and either is progress where leaving it is not.",
        ""]
    for k in orphan: L.append("- `%s` (`%s`)"%(k,"`, `".join(sw[k])))
    return "\n".join(L)+"\n"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--check',action='store_true')
    a=ap.parse_args()
    sw,desc=scan(); h=header(sw); d=document(sw,docs_text(),src_text(),desc)
    if a.check:
        bad=[]
        for path,want in ((HDR,h),(DOC,d)):
            have=open(path).read() if os.path.exists(path) else None
            if have!=want: bad.append(path)
        if bad:
            print("stale, regenerate with tools/mkswitches.py:")
            for p in bad: print("  "+os.path.relpath(p,ROOT))
            return 1
        print("switch registry is up to date (%d switches)"%len(sw))
        return 0
    os.makedirs(os.path.dirname(DOC),exist_ok=True)
    open(HDR,'w').write(h); open(DOC,'w').write(d)
    print("wrote %s and %s: %d switches"%(os.path.relpath(HDR,ROOT),
                                          os.path.relpath(DOC,ROOT),len(sw)))
    return 0

if __name__=='__main__': sys.exit(main())
