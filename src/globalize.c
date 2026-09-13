/*     CalculiX - damage/fracture extension                              */
/*     globalize.c: which globalization mechanism actually did anything. */

/* Why this module exists
   ----------------------
   handover/12-GLOBALIZATION.md has the argument.  The short form:
   nonlingeo.c stacks six globalization mechanisms in a fixed order - path
   following, transactional backtracking, two rescue levels, a dogleg trust
   region and the adaptive line-search ladder - and 05-DEBT.md section 3
   records that at one wall THREE CONSECUTIVE ATTEMPTS produced
   bit-identical residual sequences.  Two rescue levels ran and changed
   nothing.

   That is an anecdote.  This module turns it into a counter.

   The brief's test for keeping a mechanism is to name the failure it
   addresses and the gate case that would go red without it.  Nobody can do
   that for these six because nobody has the data, so the first step here is
   an instrument and not an extraction.  A mechanism that never changes an
   iterate has no defender, and then deleting it is an argument from data.

   How it measures
   ---------------
   Per ATTEMPT, not per mechanism, because that is the level at which the
   observation was made.  Each mechanism says "I acted" at the site where it
   already announces itself; at the end of the attempt the correction vector
   is hashed and compared with the previous attempt's.  Same hash, after a
   mechanism claimed to act, means the mechanism changed nothing that
   reached the solve.

   The hash is topodiag's, which is already in the tree and already has a
   self test - measuring with an instrument that has never been checked is
   the failure this project keeps finding in other people's code.

   Reason codes follow PETSc SNESLineSearchReason (petscsnes.h:918-924):
   SUCCEEDED, and the failures named separately rather than collapsed into
   "it did not help".

   Contract
   --------
     - it counts and it reports; it decides nothing and changes no
       arithmetic, so a run with it is bit-identical to a run without it;
     - a mechanism that does not call glob_fired is invisible here, and the
       report says how many attempts had no mechanism at all so that
       missing instrumentation cannot be mistaken for a quiet solver.   */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

static const glob_census *glob_active=NULL;
static ITG glob_reported=0;

static const char *glob_mech_name_tab[GLOB_NMECH]={
  "line-search ladder",
  "transactional backtrack",
  "rescue level 1",
  "rescue level 2",
  "trust region (dogleg)",
  "path following"
};

const char *glob_mech_name(glob_mech m)
{
  if((m<0)||(m>=GLOB_NMECH)) return "unknown";
  return glob_mech_name_tab[m];
}

void glob_census_init(glob_census *g)
{
  ITG k;
  for(k=0;k<GLOB_NMECH;k++){
    g->fired[k]=0; g->effective[k]=0; g->inert[k]=0; g->unknown[k]=0;
  }
  g->attempts=0; g->bare_attempts=0; g->repeats=0;
  g->prev_hash=0ULL; g->cur_hash=0ULL; g->have_prev=0;
  g->pending=0; g->niter=0;
}

/* Mechanism m acted in the current attempt.  Called at the site where the
   mechanism already prints that it fired, so instrumenting it is one line
   and cannot drift from what the log says. */
void glob_fired(glob_census *g,glob_mech m)
{
  if((m<0)||(m>=GLOB_NMECH)) return;
  g->fired[m]++;
  g->pending|=(1u<<m);
}

/* One Newton iteration of the current attempt; b is the correction the
   solve produced.  Folded into a running hash rather than compared here,
   because the observation this exists to count was made over a SEQUENCE:
   three consecutive attempts whose residual sequences were bit-identical.
   Comparing consecutive iterations would answer a different question and
   would almost always say "different". */
void glob_iterate(glob_census *g,const double *b,ITG n)
{
  if(g->niter==0) g->cur_hash=topodiag_hash_seed();
  g->cur_hash=topodiag_hash_d(b,n,g->cur_hash);
  g->niter++;
}

/* The attempt is over.  Its whole iteration sequence hashing to what the
   previous attempt's did means every mechanism that claimed to act in
   between changed nothing that reached the solve. */
void glob_attempt_end(glob_census *g)
{
  unsigned long long h;
  ITG k;
  unsigned same;

  if(g->niter==0) return;     /* nothing was solved: not an attempt */

  g->attempts++;
  h=g->cur_hash;

  same=(g->have_prev&&(h==g->prev_hash));
  if(same) g->repeats++;

  if(g->pending==0u){
    g->bare_attempts++;
  }else{
    for(k=0;k<GLOB_NMECH;k++){
      if(g->pending&(1u<<k)){
        /* The FIRST attempt has no predecessor, so nothing about this
           mechanism has been demonstrated.  Counting it effective would
           inflate exactly the number that protects a useless mechanism
           from deletion - the instrument erring in the direction that
           makes its own finding unreachable.  The self test caught this
           and it is why "unknown" exists. */
        if(!g->have_prev)   g->unknown[k]++;
        else if(same)       g->inert[k]++;
        else                g->effective[k]++;
      }
    }
  }

  g->prev_hash=h; g->have_prev=1; g->pending=0u; g->niter=0;
}

void glob_census_report(const glob_census *g)
{
  ITG k;

  glob_reported=1;

  printf("\n[GLOBALIZE CENSUS] %" ITGFORMAT " attempt(s), %" ITGFORMAT
         " with no mechanism, %" ITGFORMAT " that reproduced the previous "
         "correction exactly\n",g->attempts,g->bare_attempts,g->repeats);
  /* Two different units in one table is how a table gets misread:
     firings counts calls, everything after it counts ATTEMPTS.  A
     mechanism that loops trials inside one attempt shows a large firing
     count and a small attempt count, and that is not a contradiction. */
  printf("[GLOBALIZE CENSUS] %-26s %8s %8s %8s %6s %8s\n",
         "mechanism","firings","attempts","changed","inert","unknown");
  for(k=0;k<GLOB_NMECH;k++){
    printf("[GLOBALIZE CENSUS] %-26s %8" ITGFORMAT " %8" ITGFORMAT
           " %8" ITGFORMAT " %6" ITGFORMAT " %8" ITGFORMAT "%s\n",
           glob_mech_name((glob_mech)k),g->fired[k],
           g->effective[k]+g->inert[k]+g->unknown[k],
           g->effective[k],g->inert[k],g->unknown[k],
           (g->fired[k]==0)?"   (never fired)":
           ((g->effective[k]==0)?"   (NEVER CHANGED ANYTHING)":""));
  }
  fflush(stdout);
}

/* Reporting where the run actually ends.

   The first version printed only at the end of nonlingeo(), and the gate
   proved immediately why that is wrong: the four cases that hit a wall
   exit through checkconvergence's FORTRAN(stop) and produced NO census at
   all.  The runs an instrument exists for are exactly the runs that do not
   reach the bottom of the function.

   This tree already learned that once - CCX_LOG_VIEW's defender in
   docs/SWITCHES.md is two 2.3-hour runs killed part way through that
   produced no profile because the table was printed from atexit.  The
   remedy there was interim reports; the remedy here is atexit, and the
   difference matters: FORTRAN(stop) is an ORDERLY exit, so atexit runs.  A
   run killed with a signal still loses the census, and if that becomes a
   problem the answer is LOG_VIEW's, not another exit hook.

   One module-level pointer, because there is one solver and one census per
   process, and the alternative is threading a pointer through
   checkconvergence to reach an exit it does not own. */

static void glob_atexit_report(void)
{
  if((glob_active!=NULL)&&(glob_reported==0)){
    glob_reported=1;
    glob_census_report(glob_active);
  }
}

void glob_census_arm(const glob_census *g)
{
  glob_active=g;
  glob_reported=0;
  atexit(glob_atexit_report);
}

/* ---------------------------------------------------------------- tests */

static ITG glob_chki(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-38s got=%-8" ITGFORMAT " want=%-8" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG glob_selftest(void)
{
  ITG nbad=0,k;
  glob_census g;
  double b[4];

  printf("[GLOBALIZE] self test%s","\n");

  /* every mechanism has a name and no two share one - the complaint
     against a stack of six is that none of them has one */
  {
    ITG distinct=1;
    for(k=0;k<GLOB_NMECH;k++){
      ITG j;
      if(strcmp(glob_mech_name((glob_mech)k),"unknown")==0) distinct=0;
      for(j=k+1;j<GLOB_NMECH;j++)
        if(strcmp(glob_mech_name((glob_mech)k),
                  glob_mech_name((glob_mech)j))==0) distinct=0;
    }
    glob_chki("all six mechanisms named, no duplicates",distinct,1,&nbad);
  }

  glob_census_init(&g);
  glob_chki("a fresh census has counted nothing",
            (g.attempts==0)&&(g.repeats==0),1,&nbad);

  /* attempt 1: the ladder acts, the correction is something */
  for(k=0;k<4;k++) b[k]=(double)(k+1);
  glob_fired(&g,GLOB_LADDER);
  glob_iterate(&g,b,4);
  glob_attempt_end(&g);
  glob_chki("first attempt cannot repeat anything",g.repeats,0,&nbad);
  glob_chki("  the ladder is NOT credited on attempt 1",
            g.effective[GLOB_LADDER],0,&nbad);
  glob_chki("  it is recorded unknown instead",
            g.unknown[GLOB_LADDER],1,&nbad);

  /* attempt 2: rescue 1 acts and the correction comes back IDENTICAL -
     this is the debt entry's observation, as a counter */
  glob_fired(&g,GLOB_RESCUE1);
  glob_iterate(&g,b,4);
  glob_attempt_end(&g);
  glob_chki("an identical correction is a repeat",g.repeats,1,&nbad);
  glob_chki("  and the mechanism is recorded inert",
            g.inert[GLOB_RESCUE1],1,&nbad);
  glob_chki("  and NOT credited as effective",
            g.effective[GLOB_RESCUE1],0,&nbad);

  /* attempt 3: two mechanisms act together and the correction changes -
     both are credited, because at this resolution neither can be
     separated from the other and claiming otherwise would be a guess */
  b[2]=99.;
  glob_fired(&g,GLOB_RESCUE1);
  glob_fired(&g,GLOB_RESCUE2);
  glob_iterate(&g,b,4);
  glob_attempt_end(&g);
  glob_chki("two mechanisms in one attempt are both credited",
            (g.effective[GLOB_RESCUE1]==1)&&(g.effective[GLOB_RESCUE2]==1),
            1,&nbad);

  /* attempt 4: nobody acts.  An attempt with no mechanism must be
     counted separately, or missing instrumentation reads as a quiet
     solver and the census lies by omission. */
  b[0]=-5.;
  glob_iterate(&g,b,4);
  glob_attempt_end(&g);
  glob_chki("an attempt with no mechanism is counted apart",
            g.bare_attempts,1,&nbad);
  glob_chki("total attempts",g.attempts,4,&nbad);

  /* the pending set must not leak across attempts */
  b[1]=-7.;
  glob_iterate(&g,b,4);
  glob_attempt_end(&g);
  glob_chki("a mechanism is not credited twice",
            g.effective[GLOB_RESCUE2],1,&nbad);

  /* THE SEMANTIC THIS EXISTS FOR: a whole iteration SEQUENCE repeating,
     not one correction.  The debt entry's observation was three
     consecutive attempts with bit-identical residual sequences, and an
     instrument that compared single iterations would answer a different
     question and almost always say "different". */
  {
    glob_census q;
    double s1[3]={1.,2.,3.},s2[3]={1.,2.,4.};
    ITG i;

    glob_census_init(&q);

    /* attempt 1: three iterations */
    glob_fired(&q,GLOB_BACKTRACK);
    for(i=0;i<3;i++){ s1[0]=(double)i; glob_iterate(&q,s1,3); }
    glob_attempt_end(&q);

    /* attempt 2: the same three iterations, after a mechanism acted */
    glob_fired(&q,GLOB_BACKTRACK);
    for(i=0;i<3;i++){ s1[0]=(double)i; glob_iterate(&q,s1,3); }
    glob_attempt_end(&q);
    glob_chki("an identical SEQUENCE is a repeat",q.repeats,1,&nbad);
    glob_chki("  and the mechanism is inert",q.inert[GLOB_BACKTRACK],1,&nbad);

    /* attempt 3: the same but for the last iteration only */
    glob_fired(&q,GLOB_BACKTRACK);
    for(i=0;i<2;i++){ s1[0]=(double)i; glob_iterate(&q,s1,3); }
    glob_iterate(&q,s2,3);
    glob_attempt_end(&q);
    glob_chki("a sequence differing in ONE iteration is not a repeat",
              q.repeats,1,&nbad);
    glob_chki("  and the mechanism is effective",
              q.effective[GLOB_BACKTRACK],1,&nbad);
    glob_chki("  with attempt 1 still unknown, not effective",
              q.unknown[GLOB_BACKTRACK],1,&nbad);

    /* an attempt in which nothing was solved is not an attempt: a cutback
       that never reaches a solve must not dilute the denominator */
    glob_attempt_end(&q);
    glob_chki("an attempt with no iteration is not counted",q.attempts,3,
              &nbad);
  }

  /* firings and attempts are different units and the table says so; a
     mechanism that loops inside one attempt must show many firings and
     one attempt, not a contradiction */
  {
    glob_census u;
    double v[2]={1.,2.};
    glob_census_init(&u);
    glob_iterate(&u,v,2); glob_attempt_end(&u);      /* a predecessor */
    glob_fired(&u,GLOB_TRUSTREGION);
    glob_fired(&u,GLOB_TRUSTREGION);
    glob_fired(&u,GLOB_TRUSTREGION);
    v[0]=9.; glob_iterate(&u,v,2); glob_attempt_end(&u);
    glob_chki("three firings in one attempt count as three",
              u.fired[GLOB_TRUSTREGION],3,&nbad);
    glob_chki("  but as ONE attempt",
              u.effective[GLOB_TRUSTREGION]+u.inert[GLOB_TRUSTREGION]+
              u.unknown[GLOB_TRUSTREGION],1,&nbad);
  }

  printf("[GLOBALIZE] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
