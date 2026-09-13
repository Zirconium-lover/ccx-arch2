/*     Where the run actually spends its time.
 *
 *     The precondition nobody had met: this code has never been profiled.
 *     The target deck takes 2.3 hours and the candidate optimisations in
 *     handover/08-OBJECT-MODEL.md section 4 - a sparsity pattern held fixed
 *     under erosion, an interpolating line search instead of a fixed ladder
 *     of full residual evaluations, a tangent reuse policy - are each a
 *     hypothesis about where that time goes.  A hypothesis about a number
 *     nobody has measured is a guess.
 *
 *     This is PETSc's -log_view idea at the scale this tree needs it:
 *     named events, entered and left explicitly, each carrying a call count,
 *     an INCLUSIVE time (the event and everything under it) and a SELF time
 *     (the event with its children removed), plus who called whom.  It is a
 *     measuring instrument, not a decision: nothing here is on a solution
 *     path and nothing here reads a physical quantity.
 *
 *     Why explicit events rather than a sampling profiler.  Two reasons that
 *     are specific to this code and not a matter of taste:
 *
 *       - the questions are about CALLS, not about lines.  "How many full
 *         residual evaluations does the line-search ladder consume per
 *         increment" is a count, and a sampler cannot answer it.
 *       - the expensive work is inside MKL, which is a static library with
 *         no symbols worth sampling.  An event around the call attributes it
 *         correctly; a PC histogram does not.
 *
 *     Contract:
 *
 *       - OFF unless CCX_LOG_VIEW is set, and off means two predictable
 *         branches per instrumented call and nothing else.  No arithmetic on
 *         a physical quantity is touched in either state, which is why the
 *         feature-off arm is bit-identical by construction and is checked
 *         anyway.
 *       - the self test runs before anything is reported, and the report is
 *         SUPPRESSED if it fails: a timing table that may be wrong is worse
 *         than none, because it will be believed.
 *       - one line of the report is machine readable, so a profile can be
 *         diffed between two arms instead of read by eye.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "CalculiX.h"

#define LOGVIEW_MAXEVENT 48
#define LOGVIEW_MAXDEPTH 32
#define LOGVIEW_ROOT     LOGVIEW_MAXEVENT   /* the caller of the outermost */

static ITG    logview_armed=-1;     /* -1 not yet decided, 0 off, 1 on     */
static ITG    logview_reported=0;   /* the FINAL report is printed once     */
static double logview_t0=0.;        /* when the instrument was armed        */
static double logview_every=600.;   /* seconds between interim reports      */
static double logview_last=0.;      /* when the last interim report went out*/
static ITG    logview_ninterim=0;
static ITG    logview_nevent=0;
static ITG    logview_broken=0;     /* an unbalanced end or a full stack   */
static const char *logview_name[LOGVIEW_MAXEVENT];
static ITG    logview_calls[LOGVIEW_MAXEVENT];
static double logview_incl[LOGVIEW_MAXEVENT];
static double logview_self[LOGVIEW_MAXEVENT];
/* who called whom: [parent][child], parent LOGVIEW_ROOT means top level */
static ITG    logview_pcalls[LOGVIEW_MAXEVENT+1][LOGVIEW_MAXEVENT];
static double logview_pincl[LOGVIEW_MAXEVENT+1][LOGVIEW_MAXEVENT];
/* the open events, innermost last */
static ITG    logview_depth=0;
static ITG    logview_stack[LOGVIEW_MAXDEPTH];
static double logview_start[LOGVIEW_MAXDEPTH];
static double logview_child[LOGVIEW_MAXDEPTH];  /* time spent in children */

static double logview_now(void){
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1.e-9*(double)t.tv_nsec;
}

static void logview_atexit(void);
static void logview_table(double totalseconds,ITG interim);

ITG logview_enabled(void){
  if(logview_armed<0){
    const char *e=ccxopt_getenv("CCX_LOG_VIEW");
    logview_armed=((e!=NULL)&&(e[0]!='\0')&&(strcmp(e,"0")!=0))?1:0;
    if(logview_armed){
      const char *w=ccxopt_getenv("CCX_LOG_VIEW_EVERY");
      if(w!=NULL) logview_every=atof(w);
      logview_t0=logview_now();
      logview_last=logview_t0;
      /* Most of the runs worth profiling do not reach the end of main: a
         deck that walls leaves through the *ERROR path, and the wrapped gate
         case is exactly such a run.  A profile you can only collect from a
         successful run is a profile of the wrong runs. */
      atexit(logview_atexit);
    }
  }
  return logview_armed;
}

/* Idempotent: the same name always gets the same id, so a caller can hold a
   static id and pay the lookup once.  Returns -1 when the table is full,
   which every entry point then treats as "not instrumented" rather than as
   an error - a missing row is visible in the report, a crash is not. */
ITG logview_event(const char *name){
  ITG i;
  for(i=0;i<logview_nevent;i++){
    if(strcmp(logview_name[i],name)==0) return i;
  }
  if(logview_nevent>=LOGVIEW_MAXEVENT) return -1;
  logview_name[logview_nevent]=name;
  return logview_nevent++;
}

void logview_begin(ITG id){
  if(!logview_enabled()) return;
  if((id<0)||(id>=logview_nevent)) return;
  if(logview_depth>=LOGVIEW_MAXDEPTH){logview_broken++;return;}
  logview_stack[logview_depth]=id;
  logview_child[logview_depth]=0.;
  logview_start[logview_depth]=logview_now();
  logview_depth++;
}

void logview_end(ITG id){
  double t,dt;
  ITG parent;
  if(!logview_enabled()) return;
  if((id<0)||(id>=logview_nevent)) return;
  if(logview_depth<=0){logview_broken++;return;}
  if(logview_stack[logview_depth-1]!=id){
    /* an event was left in the wrong order, so every time below this point
       is suspect; say so rather than reporting a plausible table */
    logview_broken++;return;
  }
  t=logview_now();
  logview_depth--;
  dt=t-logview_start[logview_depth];
  logview_calls[id]++;
  logview_incl[id]+=dt;
  logview_self[id]+=dt-logview_child[logview_depth];
  parent=(logview_depth>0)?logview_stack[logview_depth-1]:LOGVIEW_ROOT;
  logview_pcalls[parent][id]++;
  logview_pincl[parent][id]+=dt;
  if(logview_depth>0) logview_child[logview_depth-1]+=dt;

  /* An interim report, at the top level only.
     Two 2.3-hour runs of the target deck were killed part way through and
     produced NO profile at all, because the table was printed from atexit.
     A profiler that reports only at the end is useless on exactly the runs
     it exists for: the long ones, which are also the ones most likely to be
     interrupted.  The cost is one table per CCX_LOG_VIEW_EVERY seconds. */
  if((logview_depth==0)&&(logview_every>0.)&&
     (t-logview_last>=logview_every)){
    logview_last=t;
    logview_ninterim++;
    logview_table(t-logview_t0,1);
  }
}

/* Self test.  What it has to establish is that the two numbers the report
   is read for - inclusive and self - mean what they say when events nest,
   and that a mis-nesting is DETECTED rather than absorbed.  Everything else
   in the module is bookkeeping around those two facts. */
ITG logview_selftest(void){
  ITG bad=0,outer,inner,other,i,save_armed,save_nevent,save_depth,save_broken;
  double t0,x=0.;
  ITG s_calls[LOGVIEW_MAXEVENT];
  double s_incl[LOGVIEW_MAXEVENT],s_self[LOGVIEW_MAXEVENT];

  /* the test measures with the real machinery, so the run's own totals are
     saved and restored: a self test that pollutes what it certifies is not
     a self test */
  save_armed=logview_armed; save_nevent=logview_nevent;
  save_depth=logview_depth; save_broken=logview_broken;
  for(i=0;i<LOGVIEW_MAXEVENT;i++){
    s_calls[i]=logview_calls[i];s_incl[i]=logview_incl[i];
    s_self[i]=logview_self[i];
  }
  logview_armed=1; logview_depth=0; logview_broken=0;

  outer=logview_event("selftest outer");
  inner=logview_event("selftest inner");
  other=logview_event("selftest other");
  if((outer<0)||(inner<0)||(other<0)){
    printf("[LOGVIEW] *ERROR: the event table is too small for the self test\n");
    bad++;
  }else{
    if(logview_event("selftest outer")!=outer){
      printf("[LOGVIEW] *ERROR: registration is not idempotent\n");bad++;}
    if(outer==inner){
      printf("[LOGVIEW] *ERROR: two names share one id\n");bad++;}

    logview_calls[outer]=0;logview_incl[outer]=0.;logview_self[outer]=0.;
    logview_calls[inner]=0;logview_incl[inner]=0.;logview_self[inner]=0.;

    /* the clock has to advance over an interval this machine cannot
       optimise away, or nothing below means anything */
    t0=logview_now();
    for(i=0;i<2000000;i++) x+=1./(1.+(double)i);
    if(logview_now()-t0<=0.){
      printf("[LOGVIEW] *ERROR: the monotonic clock did not advance\n");bad++;}

    /* outer( idle ; inner( busy ) ; idle ) : inner must appear in outer's
       inclusive time and must NOT appear in outer's self time */
    logview_begin(outer);
    for(i=0;i<1000000;i++) x+=1./(1.+(double)i);
    logview_begin(inner);
    for(i=0;i<4000000;i++) x+=1./(1.+(double)i);
    logview_end(inner);
    for(i=0;i<1000000;i++) x+=1./(1.+(double)i);
    logview_end(outer);

    if(x<=0.){printf("[LOGVIEW] *ERROR: the busy loop was elided\n");bad++;}
    if(logview_calls[outer]!=1){
      printf("[LOGVIEW] *ERROR: outer counted %" ITGFORMAT " calls, not 1\n",
             logview_calls[outer]);bad++;}
    if(logview_incl[inner]<=0.){
      printf("[LOGVIEW] *ERROR: the inner event measured no time\n");bad++;}
    if(logview_incl[outer]<=logview_incl[inner]){
      printf("[LOGVIEW] *ERROR: inclusive time does not contain the child\n");
      bad++;}
    /* the identity the whole report rests on, to the clock's own accuracy */
    if(fabs((logview_self[outer]+logview_incl[inner])-logview_incl[outer])
       >1.e-6*logview_incl[outer]+1.e-7){
      printf("[LOGVIEW] *ERROR: self+child is not inclusive: %e + %e != %e\n",
             logview_self[outer],logview_incl[inner],logview_incl[outer]);bad++;}
    /* and the child must dominate: it ran four times the work */
    if(logview_incl[inner]<=logview_self[outer]){
      printf("[LOGVIEW] *ERROR: the child did not dominate its parent\n");
      bad++;}
    if(logview_pcalls[outer][inner]!=1){
      printf("[LOGVIEW] *ERROR: the caller of the inner event was not recorded\n");
      bad++;}

    /* a mis-nesting must be caught, not absorbed */
    logview_broken=0;
    logview_begin(outer);
    logview_end(other);          /* wrong id */
    logview_end(outer);
    if(logview_broken==0){
      printf("[LOGVIEW] *ERROR: leaving an event out of order went unnoticed\n");
      bad++;}
    logview_broken=0;
    logview_end(other);          /* nothing is open */
    if(logview_broken==0){
      printf("[LOGVIEW] *ERROR: closing an event that is not open went unnoticed\n");
      bad++;}
    if(logview_depth!=0){
      printf("[LOGVIEW] *ERROR: the stack did not unwind: depth %" ITGFORMAT "\n",
             logview_depth);bad++;}
  }

  /* undo everything the test did, including the three events it registered */
  for(i=0;i<LOGVIEW_MAXEVENT;i++){
    logview_calls[i]=s_calls[i];logview_incl[i]=s_incl[i];
    logview_self[i]=s_self[i];
  }
  for(i=save_nevent;i<logview_nevent;i++){
    ITG p;
    for(p=0;p<=LOGVIEW_MAXEVENT;p++){logview_pcalls[p][i]=0;logview_pincl[p][i]=0.;}
    logview_name[i]=NULL;
  }
  for(i=0;i<=LOGVIEW_MAXEVENT;i++){
    ITG c;
    for(c=save_nevent;c<logview_nevent;c++){logview_pcalls[i][c]=0;logview_pincl[i][c]=0.;}
  }
  logview_nevent=save_nevent; logview_armed=save_armed;
  logview_depth=save_depth;   logview_broken=save_broken;
  return bad;
}

/* The two calls a caller actually writes.  The name lookup is a strcmp over
   a table of at most 48 entries and it happens only when the instrument is
   armed, which is a few microseconds against events that cost milliseconds;
   paying it buys call sites that carry their own name and cannot get out of
   step with an enum. */
void logview_begin_named(const char *name){
  if(!logview_enabled()) return;
  logview_begin(logview_event(name));
}

void logview_end_named(const char *name){
  if(!logview_enabled()) return;
  logview_end(logview_event(name));
}

/* Called from the end of main with the total that run printed, and from
   atexit with nothing, whichever happens first. */
static void logview_atexit(void){
  logview_report(logview_now()-logview_t0);
}

void logview_report(double totalseconds){
  if(!logview_enabled()) return;
  if(logview_reported) return;
  logview_reported=1;
  if(logview_ninterim>0)
    printf("[LOGVIEW] final report; %" ITGFORMAT " interim report(s) preceded "
           "it in this log\n",logview_ninterim);
  logview_table(totalseconds,0);
}

static void logview_table(double totalseconds,ITG interim){
  ITG i,j,order[LOGVIEW_MAXEVENT],n,a,b,tmp;
  if(logview_nevent<=0){
    printf("[LOGVIEW] armed, but nothing was instrumented\n");fflush(stdout);return;}
  if(logview_selftest()!=0){
    printf("[LOGVIEW] the self test failed; reporting nothing rather than "
           "reporting a timing table that may be wrong\n");
    fflush(stdout);return;
  }
  if(logview_depth!=0){
    printf("[LOGVIEW] *ERROR: %" ITGFORMAT " event(s) still open at the report\n",
           logview_depth);
    logview_broken++;
  }
  if(logview_broken){
    printf("[LOGVIEW] the event nesting was violated %" ITGFORMAT " time(s); "
           "reporting nothing rather than reporting a table that may be wrong\n",
           logview_broken);
    fflush(stdout);return;
  }
  n=logview_nevent;
  for(i=0;i<n;i++) order[i]=i;
  for(i=0;i<n;i++){                      /* by self time, descending */
    for(j=i+1;j<n;j++){
      a=order[i];b=order[j];
      if(logview_self[b]>logview_self[a]){tmp=order[i];order[i]=order[j];order[j]=tmp;}
    }
  }
  printf("\n[LOGVIEW]%s %.3f s of run, %" ITGFORMAT " instrumented event(s).  "
         "self = the event with its children removed.\n",
         interim?" INTERIM":"",totalseconds,n);
  printf("[LOGVIEW] %-28s %10s %12s %12s %7s %12s\n",
         "event","calls","incl (s)","self (s)","%run","us/call");
  for(i=0;i<n;i++){
    ITG e=order[i];
    if(logview_calls[e]==0) continue;
    printf("[LOGVIEW] %-28s %10" ITGFORMAT " %12.3f %12.3f %6.2f%% %12.1f\n",
           logview_name[e],logview_calls[e],logview_incl[e],logview_self[e],
           (totalseconds>0.)?100.*logview_self[e]/totalseconds:0.,
           1.e6*logview_incl[e]/(double)logview_calls[e]);
  }
  printf("[LOGVIEW] who called whom (calls, inclusive seconds):\n");
  for(i=0;i<=LOGVIEW_MAXEVENT;i++){
    for(j=0;j<n;j++){
      if(logview_pcalls[i][j]==0) continue;
      printf("[LOGVIEW]   %-24s -> %-24s %10" ITGFORMAT " %10.3f\n",
             (i==LOGVIEW_ROOT)?"(top level)":logview_name[i],
             logview_name[j],logview_pcalls[i][j],logview_pincl[i][j]);
    }
  }
  /* one machine-readable line, so two arms can be diffed rather than read */
  printf("[LOGVIEW_JSON] {\"interim\":%d,\"total\":%.6f,\"events\":[",
         interim?1:0,totalseconds);
  for(i=0,j=0;i<n;i++){
    if(logview_calls[i]==0) continue;
    printf("%s{\"name\":\"%s\",\"calls\":%" ITGFORMAT
           ",\"incl\":%.6f,\"self\":%.6f}",(j++)?",":"",
           logview_name[i],logview_calls[i],logview_incl[i],logview_self[i]);
  }
  printf("]}\n");
  fflush(stdout);
}
