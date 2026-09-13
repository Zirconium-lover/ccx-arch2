/*     CalculiX - damage/fracture extension                              */
/*     trial.c: evaluate the residual at a trial state.                  */

/* Why this module exists
   ----------------------
   Every globalisation mechanism in this solver - the backtracking ladder,
   the transactional backtracking, the dogleg, the residual ray, the
   line-search probe, the continuation's finite differences - is built out
   of ONE operation: put a step in b, evaluate the model there, and read
   the residual.  That operation was written out by hand SEVENTEEN times
   inside nonlingeo(), each time as the same thirty-two lines:

       SFREE(v);SFREE(stx);SFREE(fn);
       ... the step goes into b ...
       MNEW(v,...); cpypardou(v,vold,...); NNEW(stx,...); MNEW(fn,...);
       if(ne1d2d==1)NNEW(inum,...);
       results( <130 arguments, twenty-three lines> );
       if(ne1d2d==1)SFREE(inum);
       calcresidual( <thirty-five arguments, seven lines> );

   Fifteen of those seventeen are character-identical apart from the array
   the residual lands in.  Measured on the parent commit: 28 calls to
   results() of which 25 have a character-identical argument list, and 20
   calls to calcresidual() of which 19 differ only in the destination.
   That is roughly 700 lines of one function that are one operation with
   no name.

   And it is the reason nothing could be moved out.  A block that needs to
   evaluate a residual needs those 130 locals; as long as the only way to
   name them was to be inside nonlingeo(), every mechanism had to be
   written inside nonlingeo() too.  The trust-region dogleg lives there
   for that reason and no other.

   What trialctx IS, and what it is not
   ------------------------------------
   It is NOT an abstraction of the model.  It is the argument list of
   results() and calcresidual(), written down once, and it is honest about
   that: 180 fields, one per argument, in the order the call makes them.
   The model this describes is one where results() takes 130 parameters;
   this file does not fix that and does not pretend to.  What it buys is
   that the list is written ONCE instead of twenty-five times, and that the
   operation has a name a caller in another file can use.

   Every field holds the ADDRESS of the caller's local, never its value.
   That is the whole trick and it is why the binding needs no maintenance:
   nonlingeo() reallocates v, stx, fn, inum, f, b and a dozen more on every
   iteration, and remastruct() swaps several of them wholesale.  A context
   holding values would be stale after the first NNEW.  One holding
   addresses cannot be.  The rule is uniform and mechanical - the field is
   `&x` for every argument the call writes as `x` or `&x` - which is why
   TRIAL_BIND and the calls below were GENERATED from the call text and the
   prototypes rather than typed out.

   The same rule is what makes trial_check() possible: every field is the
   address of a local, so every field is non-NULL, so a field the bind
   macro forgot is caught by one loop at arming instead of by a wrong
   answer at increment 300.

   Contract
   --------
     - trial_results() is the results() call and nothing else.  No
       allocation, no residual: for the sites that want the model evaluated
       in place;
     - trial_residual() is the whole atom.  It OWNS the scratch arrays v,
       stx, fn and inum for its duration: it frees what was there,
       allocates fresh, evaluates, and leaves v/stx/fn allocated exactly as
       the hand-written copies did.  The caller fills b first; that is the
       one thing it does not own, because the step is the caller's decision
       and this file takes no view on it;
     - it evaluates.  It does not decide, does not print and does not
       change b.

   One deliberate cosmetic difference: NNEW/SFREE record the TEXT of their
   argument for the CCX_LOG_ALLOC trace, so allocations made here are
   logged as `*(t->v)' rather than `v', with this file and line.  That
   trace is opt-in and the pointer is the same one; no output the gate
   compares is affected, and the gate is 144 files byte for byte.       */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void trial_results(const trialctx *t)
{
  results(*(t->co),*(t->nk),*(t->kon),*(t->ipkon),*(t->lakon),*(t->ne),
          *(t->v),*(t->stn),*(t->inum),*(t->stx),*(t->elcon),*(t->nelcon),
          *(t->rhcon),*(t->nrhcon),*(t->alcon),*(t->nalcon),*(t->alzero),
          *(t->ielmat),*(t->ielorien),*(t->norien),*(t->orab),
          *(t->ntmat_),*(t->t0),*(t->t1act),*(t->ithermal),*(t->prestr),
          *(t->iprestr),*(t->filab),*(t->eme),*(t->emn),*(t->een),
          *(t->iperturb),*(t->f),*(t->fn),*(t->nactdof),t->iout,t->qa,
          *(t->vold),*(t->b),*(t->nodeboun),*(t->ndirboun),*(t->xbounact),
          *(t->nboun),*(t->ipompc),*(t->nodempc),*(t->coefmpc),
          *(t->labmpc),*(t->nmpc),*(t->nmethod),t->cam,t->neq1,
          *(t->veold),*(t->accold),t->bet,t->gam,t->dtime,t->time,
          *(t->ttime),*(t->plicon),*(t->nplicon),*(t->plkcon),
          *(t->nplkcon),*(t->xstateini),*(t->xstiff),*(t->xstate),
          *(t->npmat_),*(t->epn),*(t->matname),*(t->mi),t->ielas,t->icmd,
          *(t->ncmat_),*(t->nstate_),*(t->stiini),*(t->vini),*(t->ikboun),
          *(t->ilboun),*(t->ener),*(t->enern),*(t->emeini),*(t->xstaten),
          *(t->eei),*(t->enerini),*(t->cocon),*(t->ncocon),*(t->set),
          *(t->nset),*(t->istartset),*(t->iendset),*(t->ialset),
          *(t->nprint),*(t->prlab),*(t->prset),*(t->qfx),*(t->qfn),
          *(t->trab),*(t->inotr),*(t->ntrans),*(t->fmpc),*(t->nelemload),
          *(t->nload),*(t->ikmpc),*(t->ilmpc),*(t->istep),t->iinc,
          *(t->springarea),t->reltime,t->ne0,*(t->thicke),*(t->shcon),
          *(t->nshcon),*(t->sideload),*(t->xloadact),*(t->xloadold),
          t->icfd,*(t->inomat),*(t->pslavsurf),*(t->pmastsurf),
          *(t->mortar),*(t->islavact),*(t->cdn),*(t->islavnode),
          *(t->nslavnode),*(t->ntie),*(t->clearini),*(t->islavsurf),
          *(t->ielprop),*(t->prop),t->energyini,*(t->energy),t->kscale,
          *(t->iponoeln),*(t->inoeln),*(t->nener),*(t->orname),
          *(t->network),*(t->ipobody),*(t->xbodyact),*(t->ibody),
          *(t->typeboun),*(t->itiefac),*(t->tieset),*(t->smscale),
          t->mscalmethod,*(t->nbody),*(t->t0g),*(t->t1g),
          *(t->islavquadel),*(t->aut),*(t->irowt),*(t->jqt),
          t->mortartrafoflag,t->intscheme,*(t->physcon),*(t->dam),
          *(t->damn),*(t->iponoel));
}

/* The reduction half: turn the model state results() just built into the
   residual the solver works in.  Five sites called this with the identical
   thirty-five arguments and differed only in where the answer went. */
void trial_reduce(const trialctx *t,double *dst)
{
  calcresidual(*(t->nmethod),*(t->neq),dst,*(t->fext),*(t->f),*(t->iexpl),
          *(t->nactdof),*(t->aux2),*(t->vold),*(t->vini),t->dtime,
          *(t->accold),*(t->nk),*(t->adb),*(t->aub),*(t->jq),*(t->irow),
          *(t->nzl),*(t->alpha),*(t->fextini),*(t->fini),*(t->islavnode),
          *(t->nslavnode),*(t->mortar),*(t->ntie),*(t->mi),*(t->nzs),
          t->nasym,t->idamping,*(t->veold),*(t->adc),*(t->auc),
          *(t->cvini),*(t->cv),t->alpham,t->num_cpus);
}

/* Scratch and evaluate, without the reduce.  Two sites want the model
   evaluated at a trial state and then read the STATE - the stored
   tractions, the damage - rather than the residual, and both had written
   these eight lines out by hand. */
void trial_evaluate(const trialctx *t)
{
  ITG *nk=*(t->nk),*ne=*(t->ne),*mi=*(t->mi);
  ITG mt=mi[1]+1,isiz;

  /* The scratch arrays belong to this operation.  The hand-written copies
     freed them BEFORE filling b and allocated after; freeing here, i.e.
     after the caller has filled b, is the only ordering difference and it
     is the safe direction - the arrays are alive for longer, never for
     less. */
  SFREE(*(t->v)); SFREE(*(t->stx)); SFREE(*(t->fn));
  MNEW(*(t->v),double,mt**nk);
  isiz=mt**nk; cpypardou(*(t->v),*(t->vold),&isiz,t->num_cpus);
  NNEW(*(t->stx),double,6*mi[0]**ne);
  MNEW(*(t->fn),double,mt**nk);
  if(*(t->ne1d2d)==1) NNEW(*(t->inum),ITG,*nk);

  trial_results(t);

  if(*(t->ne1d2d)==1) SFREE(*(t->inum));
}

/* The whole atom: scratch, evaluate, reduce.  The caller fills b first -
   the step is the caller's decision and this file takes no view on it. */
void trial_residual(const trialctx *t,double *dst)
{
  trial_evaluate(t);
  trial_reduce(t,dst);
}

/* Every field is the address of one of the caller's locals, so every field
   is non-NULL.  A field TRIAL_BIND forgot is therefore a NULL, and this
   loop is the only thing between that and a wrong answer three hundred
   increments later.  It walks the struct as an array of pointers, which is
   legal here because every member is a pointer and C guarantees no padding
   between members of identical type. */
ITG trial_check(const trialctx *t)
{
  const void *const *p=(const void *const *)t;
  ITG n=(ITG)(sizeof(trialctx)/sizeof(void *)),i,nbad=0;

  for(i=0;i<n;i++) if(p[i]==NULL) nbad++;
  printf("[TRIAL] context bound: %" ITGFORMAT " field(s), %" ITGFORMAT
         " unbound -- %s\n",n,nbad,nbad?"FAILED":"PASSED");
  fflush(stdout);
  return nbad;
}
