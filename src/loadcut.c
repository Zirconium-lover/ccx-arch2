/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2025 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);             */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

/*
  How WIDE is the remaining load path.

  WHY THIS EXISTS
  ---------------
  damconnect.f answers whether a chain of surviving elements still links
  the two grips, and its own comment states the limit honestly:

      "Elements that are merely degraded still conduct: a nonzero
       stiffness is a load path, and calling it broken would be a
       modelling decision rather than a topological fact."

  That is correct, and it is why the question can never be answered there.
  Measured on the target deck at the end of a complete run
  (research/09-SEVERANCE.md): every one of the six topological rules the
  tree implements says CONNECTED, and the two grips are in fact joined by
  ONE triangular element face of area 0.0109 against a nominal section of
  5.7 - five ten-thousandths of a cross-section, carrying 2.1% of peak
  load.  The specimen is two pieces on a hinge and the boolean is right to
  say otherwise.

  The history is the same argument one level down.  CCX_FRACTURE_LINK
  moved from NODE to FACE conduction after a run finished in two pieces of
  12500 and 13300 elements joined through a single VERTEX (E-75).  A
  vertex is not a load path; the fix was a stricter boolean.  A single
  face is not a load path either - three shared nodes transmit force and
  no bending, so the joint has near-zero-energy rotational modes - and a
  stricter boolean would not have caught it, because the connection
  genuinely IS a whole face.  The question was never boolean.

  WHAT IS COMPUTED
  ----------------
  The minimum cut between the two node sets, over the graph whose vertices
  are surviving elements and whose edges are shared faces, each weighted

      capacity = (area of the shared triangle) * min(g_i, g_j)

  with g = 1 - D floored at gmin, the same residual stiffness the assembly
  uses.  By max-flow/min-cut that is the smallest total weighted area one
  would have to break to separate the grips: the width of the load path,
  in the units of an area.

  It subsumes the boolean.  Cut = 0 is exactly "disconnected"; a cut of
  0.05% of a section is a specimen that is topologically whole and
  mechanically finished, which is the state this routine exists to name.

  COST, AND WHY IT IS NOT A PROBLEM
  ---------------------------------
  Max-flow is run with an EARLY EXIT at a caller-supplied target.  A
  healthy specimen stops after a handful of augmenting paths, because the
  target is reached immediately; a nearly severed one stops because there
  is almost nothing to push.  The expensive case - a cut just at the
  target - is the only one that runs to completion, and it is the case
  worth paying for.  It is evaluated where damconnectsets is evaluated
  today: once per committed deletion batch, not per iteration.

  RESTRICTIONS, CHECKED RATHER THAN ASSUMED
  -----------------------------------------
  C3D4 bulk and 6-node cohesive facets, which is what the progressive
  damage path supports anyway.  Any other live element type makes the
  routine REFUSE and report, rather than guess: an absent measurement is
  recoverable and a wrong one is not.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalculiX.h"

#define LC_INF 1.e30

/* ------------------------------------------------------------------ */
/* Pure geometry                                                       */
/* ------------------------------------------------------------------ */

static double lc_tri(const double *co,ITG a,ITG b,ITG c)
{
  double u[3],v[3],w[3];
  ITG k;
  for(k=0;k<3;k++){
    u[k]=co[3*(b-1)+k]-co[3*(a-1)+k];
    v[k]=co[3*(c-1)+k]-co[3*(a-1)+k];
  }
  w[0]=u[1]*v[2]-u[2]*v[1];
  w[1]=u[2]*v[0]-u[0]*v[2];
  w[2]=u[0]*v[1]-u[1]*v[0];
  return 0.5*sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
}

/* ------------------------------------------------------------------ */
/* A small max-flow, written so the self test can drive it directly    */
/* ------------------------------------------------------------------ */

typedef struct {
  ITG n;            /* vertices                                        */
  ITG ne,ncap;      /* directed edges used / allocated                 */
  ITG *head,*next;  /* adjacency as linked lists                       */
  ITG *to;
  double *cap;
} lcgraph;

static ITG lc_init(lcgraph *G,ITG n,ITG ecap)
{
  ITG i;
  G->n=n; G->ne=0; G->ncap=2*ecap+4;
  G->head=(ITG*)malloc(sizeof(ITG)*(size_t)(n>0?n:1));
  G->next=(ITG*)malloc(sizeof(ITG)*(size_t)G->ncap);
  G->to  =(ITG*)malloc(sizeof(ITG)*(size_t)G->ncap);
  G->cap =(double*)malloc(sizeof(double)*(size_t)G->ncap);
  if((G->head==NULL)||(G->next==NULL)||(G->to==NULL)||(G->cap==NULL)) return 0;
  for(i=0;i<n;i++) G->head[i]=-1;
  return 1;
}

static void lc_free(lcgraph *G)
{
  free(G->head); free(G->next); free(G->to); free(G->cap);
  G->head=NULL; G->next=NULL; G->to=NULL; G->cap=NULL;
}

/* One UNDIRECTED edge of capacity c: two directed arcs that are each
   other's residual, which is the standard encoding and is what makes an
   undirected min cut come out of a directed max flow. */

static void lc_edge(lcgraph *G,ITG u,ITG v,double c)
{
  if(G->ne+2>G->ncap) return;
  G->to[G->ne]=v; G->cap[G->ne]=c; G->next[G->ne]=G->head[u]; G->head[u]=G->ne++;
  G->to[G->ne]=u; G->cap[G->ne]=c; G->next[G->ne]=G->head[v]; G->head[v]=G->ne++;
}

/* Max flow, by Dinic with CAPACITY SCALING, and with a work budget.

   CORRECTION.  An earlier version of this comment claimed plain
   Edmonds-Karp HUNG the target deck at increment 173, and a second version
   claimed the run was killed externally.  Both were wrong: the run was
   alive and progressing throughout, past increment 430 and still going.

   Two bad liveness checks produced them, and the habit is worth recording
   because neither looks wrong when you read it:

     - `pgrep -c ccx_2.23_pardiso` without -f returns ZERO for a running
       process.  Linux truncates a process's comm field to 15 characters
       and this name is 16.  pgrep prints a warning saying exactly that,
       and the check had stderr redirected to /dev/null.
     - a file's mtime read as "last write" is only "it stopped there" if
       you compare it against the CURRENT time.  Read on its own it is
       just now.

   Use `ps -eo pid,etimes,comm`, or `pgrep -af` with the full path, and
   check that a counter is CHANGING rather than that a timestamp exists.

   What IS true, and is why the algorithm changed anyway: the exact mode
   is expensive.  Same deck, same machine, 4 threads - 599 increments in
   69 minutes without it, 274 in 66 minutes with it, so roughly half the
   throughput.  And the reason is real: the capacities span four decades,
   from a healthy face at ~1e-2 down to area*gmin ~1e-6 against an element
   pinned at the residual-stiffness floor.  Each augmenting path moves
   only its own bottleneck, so plain augmentation degrades badly as damage
   spreads.  With integer capacities Edmonds-Karp is bounded; with real
   ones spanning decades it is not.

   Scaling fixes exactly that.  Work at a threshold D, admitting only arcs
   whose residual capacity is at least D, and halve D when no more flow can
   be found at that threshold.  Every augmentation then moves at least D,
   so the augmentations per phase are bounded by flow/D and the phases are
   logarithmic in the capacity range.

   The budget is the second half of the fix.  If the flow still cannot be
   completed within it the routine REFUSES - returns a negative value -
   rather than spending a run's wall clock inside a diagnostic.  A
   measurement that is absent can be recovered; a run that never finishes
   cannot. */

static double lc_maxflow(lcgraph *G,ITG s,ITG t,double target,
                         double cmax,ITG *iwork)
{
  ITG *lvl,*it,*st,*pe,sp,u,v,e,i,naug=0;
  double flow=0.,push,D;
  const ITG budget=4000000;

  if(iwork!=NULL) *iwork=0;
  if((s<0)||(t<0)||(s>=G->n)||(t>=G->n)) return 0.;
  lvl=(ITG*)malloc(sizeof(ITG)*(size_t)G->n);
  it =(ITG*)malloc(sizeof(ITG)*(size_t)G->n);
  st =(ITG*)malloc(sizeof(ITG)*(size_t)(G->n+2));
  pe =(ITG*)malloc(sizeof(ITG)*(size_t)(G->n+2));
  if((lvl==NULL)||(it==NULL)||(st==NULL)||(pe==NULL)){
    free(lvl);free(it);free(st);free(pe);return 0.;}

  for(D=(cmax>0.)?cmax:1.;D>1.e-14;D*=0.5){
    while(1){
      if((target>0.)&&(flow>=target)) break;
      if(naug>budget) break;
      {
        ITG qh=0,qt=0,*q=st;            /* st doubles as the BFS queue */
        for(i=0;i<G->n;i++) lvl[i]=-1;
        lvl[s]=0; q[qt++]=s;
        while(qh<qt){
          u=q[qh++];
          for(e=G->head[u];e>=0;e=G->next[e]){
            v=G->to[e];
            if((lvl[v]<0)&&(G->cap[e]>=D)){lvl[v]=lvl[u]+1;q[qt++]=v;}
          }
        }
      }
      if(lvl[t]<0) break;
      for(i=0;i<G->n;i++) it[i]=G->head[i];
      while(1){                         /* one blocking flow */
        sp=0; st[0]=s; u=s;
        while(1){
          if(u==t) break;
          for(e=it[u];e>=0;e=G->next[e]){
            v=G->to[e];
            if((G->cap[e]>=D)&&(lvl[v]==lvl[u]+1)) break;
          }
          it[u]=e;
          if(e<0){                      /* dead end: retreat */
            lvl[u]=-1;
            if(sp==0) break;
            sp--; u=st[sp];
            if(it[u]>=0) it[u]=G->next[it[u]];
            continue;
          }
          pe[sp]=e; sp++; st[sp]=G->to[e]; u=st[sp];
        }
        if(u!=t) break;
        push=LC_INF;
        for(i=0;i<sp;i++) if(G->cap[pe[i]]<push) push=G->cap[pe[i]];
        if(!(push>0.)) break;
        for(i=0;i<sp;i++){G->cap[pe[i]]-=push;G->cap[pe[i]^1]+=push;}
        flow+=push; naug++;
        if((target>0.)&&(flow>=target)) break;
        if(naug>budget) break;
      }
    }
    if((target>0.)&&(flow>=target)) break;
  }
  free(lvl);free(it);free(st);free(pe);
  if(iwork!=NULL) *iwork=naug;
  if(naug>budget) return -2.;
  return flow;
}

/* ------------------------------------------------------------------ */
/* The element graph                                                   */
/* ------------------------------------------------------------------ */

static ITG lc_nope(const char *lak)
{
  if(lak[0]=='U'){
    if((lak[1]=='C')&&(lak[2]=='6')) return 6;
    return 0;
  }
  if((lak[0]=='C')&&(lak[3]=='4')) return 4;
  return -1;               /* live and not understood: refuse */
}

double loadcut_width(double *co,ITG *ipkon,ITG *kon,char *lakon,ITG ne,
                     ITG nk,const ITG *nodesa,ITG na,
                     const ITG *nodesb,ITG nb,
                     const double *dam,const ITG *mi,double gmin,
                     const double *xstate,ITG nstate,
                     const ITG *ifacdead,double target,
                     ITG *nfaces,ITG *nelem,ITG *ibelow,ITG *iexact,
                     ITG *nwork)
{
  ITG *cnt=NULL,*noel=NULL,*ipnoel=NULL,*mark=NULL,*stamp=NULL;
  ITG *shared=NULL,*isa=NULL,*isb=NULL,*eid=NULL,*rid=NULL;
  ITG i,j,k,n,e,o,nn,nlive=0,nedge=0,S,T,bad=0;
  double *g=NULL,cut=0.,a,cmax=0.;
  lcgraph G;

  if(nfaces!=NULL) *nfaces=0;
  if(nelem!=NULL)  *nelem=0;
  if(ibelow!=NULL) *ibelow=0;
  if(iexact!=NULL) *iexact=1;
  if(nwork!=NULL) *nwork=0;
  if((ne<=0)||(nk<=0)||(na<=0)||(nb<=0)) return -1.;

  eid=(ITG*)calloc((size_t)ne,sizeof(ITG));
  rid=(ITG*)malloc(sizeof(ITG)*(size_t)ne);
  g  =(double*)malloc(sizeof(double)*(size_t)ne);
  if((eid==NULL)||(rid==NULL)||(g==NULL)){free(eid);free(rid);free(g);return -1.;}

  /* which elements are live, and what residual stiffness each carries */

  for(i=0;i<ne;i++){
    eid[i]=-1;
    if(ipkon[i]<0) continue;
    if((ifacdead!=NULL)&&(ifacdead[i]!=0)) continue;
    n=lc_nope(lakon+8*i);
    if(n<0){bad++;continue;}
    if(n==0) continue;
    g[i]=1.;
    if((dam!=NULL)&&(mi!=NULL)&&(lakon[8*i]=='C')){
      double d=dam[mi[0]*i]-1.;
      if(d<0.) d=0.;
      if(d>1.) d=1.;
      g[i]=1.-d;
      if(g[i]<gmin) g[i]=gmin;
    }

    /* A cohesive facet carries its own damage and must be weighted by it.
       The first version of this gave EVERY surviving facet g=1, and the
       error is the one this whole file exists to correct, wearing
       different clothes: a separated interface treated as a rigid link.
       It was caught the only way it could be - two independent
       implementations of the same measure, this one and the offline
       analysis in research/09-SEVERANCE.md, disagreeing by a factor of
       380 on the same final state (cut 1.184 against 0.0031) because the
       offline one excluded facets entirely and this one let 5316 of them
       conduct at full strength.

       xstate index 1 is the facet damage over its integration points and
       index 3 is its failed flag, read exactly as damstate_facet_dead
       reads it.  The WORST point governs, because a facet transmits
       across its whole area and is as weak as its weakest place. */

    if((xstate!=NULL)&&(mi!=NULL)&&(nstate>=4)&&(lakon[8*i]=='U')){
      ITG ip,nip=(mi[0]<3)?mi[0]:3,ndead=0;
      double dw=0.;
      for(ip=0;ip<nip;ip++){
        double d=xstate[nstate*(ip+mi[0]*i)+1];
        if(d<0.) d=0.;
        if(d>1.) d=1.;
        if(d>dw) dw=d;
        if(xstate[nstate*(ip+mi[0]*i)+3]>=0.5) ndead++;
      }
      g[i]=1.-dw;
      if((nip>0)&&(ndead==nip)) g[i]=0.;
      if(g[i]<gmin) g[i]=gmin;
    }
    eid[i]=nlive; rid[nlive]=i; nlive++;
  }
  if(bad>0){
    printf("[LOADCUT] refusing: %" ITGFORMAT " live element(s) are neither "
           "C3D4 nor UC6.  A measurement that is absent can be recovered; "
           "one that is wrong cannot.\n",bad);
    free(eid);free(rid);free(g); return -1.;
  }
  if(nlive<=0){free(eid);free(rid);free(g); return 0.;}
  if(nelem!=NULL) *nelem=nlive;

  /* node -> elements, as a CSR built with one counting pass */

  ipnoel=(ITG*)calloc((size_t)(nk+2),sizeof(ITG));
  if(ipnoel==NULL){free(eid);free(rid);free(g);return -1.;}
  for(k=0;k<nlive;k++){
    i=rid[k]; n=lc_nope(lakon+8*i);
    for(j=0;j<n;j++){
      ITG nd=kon[ipkon[i]+j];
      if((nd>=1)&&(nd<=nk)) ipnoel[nd]++;
    }
  }
  for(i=1;i<=nk;i++) ipnoel[i]+=ipnoel[i-1];
  noel=(ITG*)malloc(sizeof(ITG)*(size_t)(ipnoel[nk]>0?ipnoel[nk]:1));
  cnt =(ITG*)calloc((size_t)(nk+2),sizeof(ITG));
  if((noel==NULL)||(cnt==NULL)){
    free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);return -1.;}
  for(k=0;k<nlive;k++){
    i=rid[k]; n=lc_nope(lakon+8*i);
    for(j=0;j<n;j++){
      ITG nd=kon[ipkon[i]+j];
      if((nd>=1)&&(nd<=nk)) noel[ipnoel[nd-1]+cnt[nd]++]=k;
    }
  }

  /* the graph: an edge wherever two live elements share three nodes */

  mark =(ITG*)calloc((size_t)nlive,sizeof(ITG));
  stamp=(ITG*)calloc((size_t)nlive,sizeof(ITG));
  shared=(ITG*)malloc(sizeof(ITG)*8);
  if((mark==NULL)||(stamp==NULL)||(shared==NULL)){
    free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
    free(mark);free(stamp);free(shared);return -1.;}
  for(k=0;k<nlive;k++) stamp[k]=-1;

  if(lc_init(&G,nlive+2,6*nlive+2*(na+nb)+8)==0){
    free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
    free(mark);free(stamp);free(shared);return -1.;}
  S=nlive; T=nlive+1;

  for(k=0;k<nlive;k++){
    i=rid[k]; n=lc_nope(lakon+8*i);
    for(j=0;j<n;j++){
      ITG nd=kon[ipkon[i]+j];
      if((nd<1)||(nd>nk)) continue;
      for(e=ipnoel[nd-1];e<ipnoel[nd];e++){
        o=noel[e];
        if(o<=k) continue;                       /* each pair once */
        if(stamp[o]!=k){stamp[o]=k;mark[o]=0;}
        mark[o]++;
      }
    }
    for(j=0;j<n;j++){
      ITG nd=kon[ipkon[i]+j];
      if((nd<1)||(nd>nk)) continue;
      for(e=ipnoel[nd-1];e<ipnoel[nd];e++){
        o=noel[e];
        if(o<=k) continue;
        if((stamp[o]!=k)||(mark[o]<3)) continue;
        mark[o]=-1;                              /* emit once */
        /* the shared nodes, taken from k's connectivity */
        nn=0;
        for(ITG j2=0;j2<n;j2++){
          ITG n2=kon[ipkon[i]+j2],m2,no=lc_nope(lakon+8*rid[o]);
          for(m2=0;m2<no;m2++)
            if(kon[ipkon[rid[o]]+m2]==n2){ if(nn<8) shared[nn++]=n2; break; }
        }
        if(nn<3) continue;
        a=lc_tri(co,shared[0],shared[1],shared[2]);
        if(!(a>0.)) continue;
        {
          double cc=a*(g[i]<g[rid[o]]?g[i]:g[rid[o]]);
          lc_edge(&G,k,o,cc);
          if(cc>cmax) cmax=cc;
        }
        nedge++;
      }
    }
  }
  if(nfaces!=NULL) *nfaces=nedge;

  /* the two grips */

  isa=(ITG*)calloc((size_t)nlive,sizeof(ITG));
  isb=(ITG*)calloc((size_t)nlive,sizeof(ITG));
  if((isa==NULL)||(isb==NULL)){
    lc_free(&G);free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
    free(mark);free(stamp);free(shared);free(isa);free(isb);return -1.;}
  for(j=0;j<na;j++){
    ITG nd=nodesa[j];
    if((nd<1)||(nd>nk)) continue;
    for(e=ipnoel[nd-1];e<ipnoel[nd];e++) isa[noel[e]]=1;
  }
  for(j=0;j<nb;j++){
    ITG nd=nodesb[j];
    if((nd<1)||(nd>nk)) continue;
    for(e=ipnoel[nd-1];e<ipnoel[nd];e++) isb[noel[e]]=1;
  }
  for(k=0;k<nlive;k++){
    if(isa[k]) lc_edge(&G,S,k,LC_INF);
    if(isb[k]) lc_edge(&G,k,T,LC_INF);
  }

  {
    ITG naug=0;
    cut=lc_maxflow(&G,S,T,target,cmax,&naug);
    if(cut>=1.e29){
      printf("[LOADCUT] refusing: at least one element touches BOTH "
             "termination sets, so no finite cut separates them.  The sets "
             "have to be disjoint in the mesh, not only in the deck.\n");
      lc_free(&G);
      free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
      free(mark);free(stamp);free(shared);free(isa);free(isb);
      return -1.;
    }
    if(cut<-1.5){
      printf("[LOADCUT] refusing: the flow did not complete within %"
             ITGFORMAT " augmentations.  Reporting nothing rather than "
             "spending the run inside a diagnostic.\n",naug);
      lc_free(&G);
      free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
      free(mark);free(stamp);free(shared);free(isa);free(isb);
      return -1.;
    }
    if(nwork!=NULL) *nwork=naug;
  }
  if((ibelow!=NULL)&&(target>0.)) *ibelow=(cut<target)?1:0;

  /* The early exit stops as soon as the flow reaches the target, so in
     that case the returned number is a LOWER BOUND on the cut and not
     the cut.  Saying which is the difference between a measurement and a
     number that looks like one. */

  if(iexact!=NULL) *iexact=((target>0.)&&(cut>=target))?0:1;

  lc_free(&G);
  free(eid);free(rid);free(g);free(ipnoel);free(noel);free(cnt);
  free(mark);free(stamp);free(shared);free(isa);free(isb);
  return cut;
}

/* ------------------------------------------------------------------ */
/* Resolving a named node set, in C                                    */
/*                                                                     */
/* CalculiX stores node sets with a trailing 'N', so FACE_X0_NSET in a  */
/* deck is FACE_X0_NSETN internally, and a negative entry in ialset     */
/* encodes a generated range a, b, -c meaning a, a+c, ... up to b.      */
/* damconnectsets.f does this in Fortran; it is repeated here rather    */
/* than called so that this file has no Fortran dependency and so that  */
/* the expansion is covered by this file's own self test.               */
/* ------------------------------------------------------------------ */

ITG loadcut_sets(char *set,ITG nset,ITG *istartset,ITG *iendset,ITG *ialset,
                 const char *name,ITG **nodes,ITG *nn)
{
  ITG i,j,id=-1,n=0,k,a,b,c,cap;
  char want[82];
  size_t L;

  *nodes=NULL; *nn=0;
  if((set==NULL)||(nset<=0)||(name==NULL)) return 0;
  L=strlen(name);
  while((L>0)&&(name[L-1]==' ')) L--;
  if((L==0)||(L>79)) return 0;
  for(i=0;i<(ITG)L;i++)
    want[i]=(char)((name[i]>='a'&&name[i]<='z')?name[i]-32:name[i]);
  want[L]='N';
  for(i=(ITG)L+1;i<81;i++) want[i]=' ';

  for(i=0;i<nset;i++){
    for(j=0;j<81;j++){
      char ch=set[81*i+j];
      if((ch>='a')&&(ch<='z')) ch=(char)(ch-32);
      if(ch!=want[j]) break;
    }
    if(j==81){id=i;break;}
  }
  if(id<0) return 0;

  /* Two passes - count, then fill - with the expansion written
     identically in both so they cannot disagree about the length.  This
     mirrors damsetcount/damsetfill exactly: on a negative entry the two
     preceding values a and b have ALREADY been emitted, so the counter
     backs up by two and the whole range a, a+c, ... b replaces them. */

  for(k=0;k<2;k++){
    n=0;
    for(j=istartset[id]-1;j<iendset[id];j++){
      if(ialset[j]>0){
        if(k==1) (*nodes)[n]=ialset[j];
        n++;
      }else{
        if(j<2) continue;
        a=ialset[j-2]; b=ialset[j-1]; c=-ialset[j];
        if(c<=0) continue;
        n-=2;
        if(n<0) n=0;
        for(i=a;i<=b;i+=c){
          if(k==1) (*nodes)[n]=i;
          n++;
        }
      }
    }
    if(k==0){
      cap=(n>0)?n:1;
      *nodes=(ITG*)malloc(sizeof(ITG)*(size_t)cap);
      if(*nodes==NULL) return 0;
    }
  }
  *nn=n;
  return 1;
}

/* ------------------------------------------------------------------ */
/* Self test                                                           */
/*                                                                     */
/* Hand-built meshes whose answer is known on paper.  The last two are  */
/* the point: one shows the measure catching a specimen that every      */
/* boolean rule calls intact, and one breaks the measure deliberately   */
/* and shows it going red.                                             */
/* ------------------------------------------------------------------ */

static ITG lc_chk(const char *what,double got,double want,double tol,ITG *nbad)
{
  ITG ok=(fabs(got-want)<=tol*(fabs(want)>1.?fabs(want):1.));
  printf("  %-4s %-54s got %.9f want %.9f\n",ok?"ok":"FAIL",what,got,want);
  if(!ok) (*nbad)++;
  return ok;
}

ITG loadcut_selftest(void)
{
  ITG nbad=0;
  ITG nk,ne,i;
  double co[3*12];
  ITG kon[4*8],ipkon[8];
  char lakon[8*8];
  double dam[8];
  ITG mi[3];
  ITG nodesa[3],nodesb[3],nf,nel,below,exact;
  double cut;

  mi[0]=1;mi[1]=3;mi[2]=1;
  memset(lakon,' ',sizeof(lakon));

  /* Two tetrahedra sharing the face (1,2,3), which is the triangle
     (0,0,0)-(1,0,0)-(0,1,0) of area 0.5.  Node 4 is above it and node 5
     below, so the two tets are distinct.  A = node 4, B = node 5. */

  nk=5; ne=2;
  {
    double c[15]={0.,0.,0.,  1.,0.,0.,  0.,1.,0.,  0.,0.,1.,  0.,0.,-1.};
    for(i=0;i<15;i++) co[i]=c[i];
  }
  kon[0]=1;kon[1]=2;kon[2]=3;kon[3]=4;   ipkon[0]=0;
  kon[4]=1;kon[5]=2;kon[6]=3;kon[7]=5;   ipkon[1]=4;
  memcpy(lakon+0,"C3D4    ",8);
  memcpy(lakon+8,"C3D4    ",8);
  dam[0]=1.;dam[1]=1.;                    /* undamaged: dam = 1 + D, D = 0 */
  nodesa[0]=4; nodesb[0]=5;

  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("two tets on one face: the cut is that face's area",
         cut,0.5,1.e-9,&nbad);

  /* Damage one of them to D = 0.6.  g = 0.4 and the cut follows the
     WEAKER of the two, because that is what the assembly would carry. */

  dam[1]=1.6;
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("a neighbour at D=0.6 scales the cut by g=0.4",
         cut,0.2,1.e-9,&nbad);

  /* The residual-stiffness floor is where a fully damaged element lands,
     and it is NOT zero: the cut has to show gmin, not vanish, because the
     element is still in the operator. */

  dam[1]=2.0;
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("a fully damaged neighbour leaves gmin, not zero",
         cut,0.5e-4,1.e-9,&nbad);
  dam[1]=1.;

  /* Delete the second tet: no path at all. */

  ipkon[1]=-1;
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("with the far element deleted the cut is zero",cut,0.,1.e-12,&nbad);
  ipkon[1]=4;

  /* THE CASE THIS FILE EXISTS FOR.

     A wide bridge and a hair-thin one, in parallel.  Every boolean rule -
     node conduction, face conduction, dead-facet exclusion - says
     CONNECTED for both, and says the same thing after the wide bridge is
     deleted.  The cut says the specimen lost 99.98% of its load path.

     Geometry: two "grip" tets joined twice.  Bridge 1 shares the big
     triangle (1,2,3); bridge 2 shares the sliver (6,7,8).  */

  nk=8; ne=4;
  {
    double c[24]={0.,0.,0.,   1.,0.,0.,   0.,1.,0.,   0.,0.,1.,
                  0.,0.,-1.,  5.,0.,0.,   5.01,0.,0., 5.,0.1,0.};
    for(i=0;i<24;i++) co[i]=c[i];
  }
  kon[0]=1;kon[1]=2;kon[2]=3;kon[3]=4;   ipkon[0]=0;   /* grip A side  */
  kon[4]=1;kon[5]=2;kon[6]=3;kon[7]=5;   ipkon[1]=4;   /* grip B side  */
  kon[8]=6;kon[9]=7;kon[10]=8;kon[11]=4; ipkon[2]=8;   /* sliver, A    */
  kon[12]=6;kon[13]=7;kon[14]=8;kon[15]=5;ipkon[3]=12; /* sliver, B    */
  memcpy(lakon+16,"C3D4    ",8);
  memcpy(lakon+24,"C3D4    ",8);
  for(i=0;i<4;i++) dam[i]=1.;
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("two bridges in parallel: the cut is their sum",
         cut,0.5+5.e-4,1.e-9,&nbad);

  ipkon[0]=-1;ipkon[1]=-1;              /* delete the wide bridge only */
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("with only the sliver left the cut collapses 99.9 percent",
         cut,5.e-4,1.e-9,&nbad);
  if(cut>0.) {
    printf("  ok   and a boolean sweep still calls this CONNECTED - which "
           "is\n       the whole reason this file exists\n");
  }

  /* The early exit must not change the verdict it is asked for. */

  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,1.e-2,&nf,&nel,&below,&exact,NULL);
  if(below!=1){ printf("  FAIL %-54s below=%" ITGFORMAT "\n",
                       "a target above the true cut reports BELOW",below);
                nbad++; }
  else printf("  ok   %-54s below=1\n","a target above the true cut reports BELOW");
  ipkon[0]=0;ipkon[1]=4;
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,1.e-2,&nf,&nel,&below,&exact,NULL);
  if(below!=0){ printf("  FAIL %-54s below=%" ITGFORMAT "\n",
                       "a target below the true cut reports NOT below",below);
                nbad++; }
  else printf("  ok   %-54s below=0\n","a target below the true cut reports NOT below");

  /* An element type the routine does not understand must make it refuse,
     not guess.  Stopping a run too early destroys a result; failing to
     report merely costs a measurement. */

  printf("  ..   the next line is the refusal this check provokes\n");
  memcpy(lakon+24,"C3D8    ",8);
  cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                    dam,mi,1.e-4,NULL,0,NULL,-1.,&nf,&nel,&below,&exact,NULL);
  lc_chk("an unsupported live element makes it refuse (-1)",cut,-1.,1.e-12,&nbad);
  memcpy(lakon+24,"C3D4    ",8);

  /* A cohesive facet must be weighted by its OWN damage.  Giving every
     surviving facet g=1 is what made this routine disagree with the
     offline analysis by a factor of 380 on the same state, and it is the
     boolean error in weighted clothing: a separated interface read as a
     rigid link.

     The mesh mirrors the real one - a zero-thickness facet between two
     tetrahedra, with the two sides carrying DUPLICATED nodes at the same
     coordinates, which is exactly how the s3rad deck builds its ZrH
     interface.  So the only path from grip to grip is tet - facet - tet,
     and the facet's own stiffness is the whole answer. */
  {
    double xs[4*3*3];
    ITG m;
    nk=8; ne=3;
    {
      double c[24]={0.,0.,0.,  1.,0.,0.,  0.,1.,0.,  0.,0.,1.,
                    0.,0.,0.,  1.,0.,0.,  0.,1.,0.,  0.,0.,-1.};
      for(i=0;i<24;i++) co[i]=c[i];
    }
    kon[0]=1;kon[1]=2;kon[2]=3;kon[3]=4;               ipkon[0]=0;
    kon[4]=1;kon[5]=2;kon[6]=3;kon[7]=5;kon[8]=6;kon[9]=7; ipkon[1]=4;
    kon[10]=5;kon[11]=6;kon[12]=7;kon[13]=8;           ipkon[2]=10;
    memcpy(lakon+0 ,"C3D4    ",8);
    memcpy(lakon+8 ,"UC6     ",8);
    memcpy(lakon+16,"C3D4    ",8);
    dam[0]=1.;dam[1]=1.;dam[2]=1.;
    nodesa[0]=4; nodesb[0]=8;
    mi[0]=3;
    for(m=0;m<4*3*3;m++) xs[m]=0.;
    cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                      dam,mi,1.e-4,xs,4,NULL,-1.,&nf,&nel,&below,&exact,NULL);
    lc_chk("an undamaged cohesive facet conducts at full area",
           cut,0.5,1.e-9,&nbad);
    for(m=0;m<3;m++) xs[4*(m+3*1)+1]=0.75;
    cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                      dam,mi,1.e-4,xs,4,NULL,-1.,&nf,&nel,&below,&exact,NULL);
    lc_chk("a facet at D=0.75 conducts at a quarter",cut,0.125,1.e-9,&nbad);
    for(m=0;m<3;m++) xs[4*(m+3*1)+3]=1.0;
    cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                      dam,mi,1.e-4,xs,4,NULL,-1.,&nf,&nel,&below,&exact,NULL);
    lc_chk("a fully failed facet drops to gmin, not to full strength",
           cut,0.5e-4,1.e-9,&nbad);

    /* and the degenerate deck: both grips on the same element */
    printf("  ..   the next line is the refusal this check provokes\n");
    nodesb[0]=1;
    cut=loadcut_width(co,ipkon,kon,lakon,ne,nk,nodesa,1,nodesb,1,
                      dam,mi,1.e-4,xs,4,NULL,-1.,&nf,&nel,&below,&exact,NULL);
    lc_chk("grips sharing an element make it refuse, not report infinity",
           cut,-1.,1.e-12,&nbad);
    mi[0]=1;
  }

  /* Capacity scaling, which is what makes the exact mode usable at all.
     A fan of many thin paths in parallel with one thick one: with plain
     augmentation the thin ones are found first and the count explodes
     with their number, which is how the target deck hung at increment
     173.  With scaling the thick path is saturated in the first phase.
     The check is on the augmentation COUNT, because that is the thing
     that would regress if somebody replaced the algorithm. */
  {
    lcgraph G2; ITG m,naug=0; double f;
    const ITG NP=400;
    if(lc_init(&G2,NP+4,2*NP+8)){
      /* 0 = source, 1 = sink; NP thin parallel arcs plus one thick */
      for(m=0;m<NP;m++){
        lc_edge(&G2,0,2+m,1.e-6);
        lc_edge(&G2,2+m,1,1.e-6);
      }
      lc_edge(&G2,0,2+NP,1.0);
      lc_edge(&G2,2+NP,1,1.0);
      f=lc_maxflow(&G2,0,1,-1.,1.0,&naug);
      if((fabs(f-(1.0+NP*1.e-6))<1.e-9)&&(naug<=NP+4)){
        printf("  ok   %-54s %d augmentations\n",
               "400 thin paths beside one thick: scaling keeps the count low",
               (int)naug);
      }else{
        printf("  FAIL %-54s f=%.9f naug=%d\n",
               "400 thin paths beside one thick: scaling keeps the count low",
               f,(int)naug);
        nbad++;
      }
      lc_free(&G2);
    }
  }

  /* The set resolver: the trailing N, the case folding, and the
     generated-range expansion a, b, -c.  It is repeated from Fortran, so
     it is exactly the kind of code that drifts silently. */
  {
    char set[81*2];
    ITG istartset[2]={1,6},iendset[2]={5,6},ialset[6];
    ITG *nds=NULL,nn2=0,ok;
    memset(set,' ',sizeof(set));
    memcpy(set+0,"GRIPN",5);
    memcpy(set+81,"OTHERN",6);
    /* set 1: 7, then the generated range 10,20,-5 -> 10,15,20 */
    ialset[0]=7; ialset[1]=10; ialset[2]=20; ialset[3]=-5; ialset[4]=99;
    ialset[5]=3;
    ok=loadcut_sets(set,2,istartset,iendset,ialset,"grip",&nds,&nn2);
    if((ok==1)&&(nn2==5)&&nds[0]==7&&nds[1]==10&&nds[2]==15&&nds[3]==20
       &&nds[4]==99)
      printf("  ok   %-54s %d nodes\n","a set name resolves, lower case, with the range expanded",(int)nn2);
    else{
      printf("  FAIL %-54s ok=%d n=%d\n","a set name resolves, lower case, with the range expanded",(int)ok,(int)nn2);
      nbad++;
    }
    free(nds); nds=NULL;
    ok=loadcut_sets(set,2,istartset,iendset,ialset,"nosuchset",&nds,&nn2);
    if(ok==0) printf("  ok   %-54s\n","a missing set is reported, not invented");
    else{ printf("  FAIL %-54s\n","a missing set is reported, not invented"); nbad++; }
    free(nds);
  }

  printf("[LOADCUT SELFTEST] %s\n",
         nbad?"FAILED":"PASSED");
  return nbad;
}
