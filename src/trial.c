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

   Where this shape comes from
   ---------------------------
   It is the standard answer in this field and it is worth naming the
   places it is already load-bearing, because none of them is a nonlinear
   FE solver that got away with 130 loose arguments either.

     - PETSc's SNES passes application state to the residual evaluation
       through one opaque `ctx' pointer (SNESSetFunction, petscsnes.h).
       Their own documentation states the purpose in as many words: a
       user-defined context is a structure in which objects are stashed,
       and it is how the library avoids global variables while the solver
       never sees application data.  trialctx is that context, with the
       opacity dropped because here the caller and the callee are the same
       program;
     - deal.II's WorkStream carries per-cell state in a ScratchData object
       for exactly this reason - so the worker does not take many
       individual arguments - and its CopyData is the output half.  The
       split here is the same: the context is what the evaluation reads,
       `dst' is where the answer goes.

   This file borrows the shape, not the ambition.  Neither of those
   libraries would write a 180-field struct; neither of them has a
   results() with 130 parameters to write it about.

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
   logged as `*(mdl->v)' rather than `v', with this file and line.  That
   trace is opt-in and the pointer is the same one; no output the gate
   compares is affected, and the gate is 144 files byte for byte.       */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#include "ccxfork.h"
void trial_results(const trialctx *mdl)
{
  results(*(mdl->co),*(mdl->nk),*(mdl->kon),*(mdl->ipkon),*(mdl->lakon),*(mdl->ne),
          *(mdl->v),*(mdl->stn),*(mdl->inum),*(mdl->stx),*(mdl->elcon),*(mdl->nelcon),
          *(mdl->rhcon),*(mdl->nrhcon),*(mdl->alcon),*(mdl->nalcon),*(mdl->alzero),
          *(mdl->ielmat),*(mdl->ielorien),*(mdl->norien),*(mdl->orab),
          *(mdl->ntmat_),*(mdl->t0),*(mdl->t1act),*(mdl->ithermal),*(mdl->prestr),
          *(mdl->iprestr),*(mdl->filab),*(mdl->eme),*(mdl->emn),*(mdl->een),
          *(mdl->iperturb),*(mdl->f),*(mdl->fn),*(mdl->nactdof),mdl->iout,mdl->qa,
          *(mdl->vold),*(mdl->b),*(mdl->nodeboun),*(mdl->ndirboun),*(mdl->xbounact),
          *(mdl->nboun),*(mdl->ipompc),*(mdl->nodempc),*(mdl->coefmpc),
          *(mdl->labmpc),*(mdl->nmpc),*(mdl->nmethod),mdl->cam,mdl->neq1,
          *(mdl->veold),*(mdl->accold),mdl->bet,mdl->gam,mdl->dtime,mdl->time,
          *(mdl->ttime),*(mdl->plicon),*(mdl->nplicon),*(mdl->plkcon),
          *(mdl->nplkcon),*(mdl->xstateini),*(mdl->xstiff),*(mdl->xstate),
          *(mdl->npmat_),*(mdl->epn),*(mdl->matname),*(mdl->mi),mdl->ielas,mdl->icmd,
          *(mdl->ncmat_),*(mdl->nstate_),*(mdl->stiini),*(mdl->vini),*(mdl->ikboun),
          *(mdl->ilboun),*(mdl->ener),*(mdl->enern),*(mdl->emeini),*(mdl->xstaten),
          *(mdl->eei),*(mdl->enerini),*(mdl->cocon),*(mdl->ncocon),*(mdl->set),
          *(mdl->nset),*(mdl->istartset),*(mdl->iendset),*(mdl->ialset),
          *(mdl->nprint),*(mdl->prlab),*(mdl->prset),*(mdl->qfx),*(mdl->qfn),
          *(mdl->trab),*(mdl->inotr),*(mdl->ntrans),*(mdl->fmpc),*(mdl->nelemload),
          *(mdl->nload),*(mdl->ikmpc),*(mdl->ilmpc),*(mdl->istep),mdl->iinc,
          *(mdl->springarea),mdl->reltime,mdl->ne0,*(mdl->thicke),*(mdl->shcon),
          *(mdl->nshcon),*(mdl->sideload),*(mdl->xloadact),*(mdl->xloadold),
          mdl->icfd,*(mdl->inomat),*(mdl->pslavsurf),*(mdl->pmastsurf),
          *(mdl->mortar),*(mdl->islavact),*(mdl->cdn),*(mdl->islavnode),
          *(mdl->nslavnode),*(mdl->ntie),*(mdl->clearini),*(mdl->islavsurf),
          *(mdl->ielprop),*(mdl->prop),mdl->energyini,*(mdl->energy),mdl->kscale,
          *(mdl->iponoeln),*(mdl->inoeln),*(mdl->nener),*(mdl->orname),
          *(mdl->network),*(mdl->ipobody),*(mdl->xbodyact),*(mdl->ibody),
          *(mdl->typeboun),*(mdl->itiefac),*(mdl->tieset),*(mdl->smscale),
          mdl->mscalmethod,*(mdl->nbody),*(mdl->t0g),*(mdl->t1g),
          *(mdl->islavquadel),*(mdl->aut),*(mdl->irowt),*(mdl->jqt),
          mdl->mortartrafoflag,mdl->intscheme,*(mdl->physcon),*(mdl->dam),
          *(mdl->damn),*(mdl->iponoel));
}

/* The reduction half: turn the model state results() just built into the
   residual the solver works in.  Five sites called this with the identical
   thirty-five arguments and differed only in where the answer went. */
void trial_reduce(const trialctx *mdl,double *dst)
{
  calcresidual(*(mdl->nmethod),*(mdl->neq),dst,*(mdl->fext),*(mdl->f),*(mdl->iexpl),
          *(mdl->nactdof),*(mdl->aux2),*(mdl->vold),*(mdl->vini),mdl->dtime,
          *(mdl->accold),*(mdl->nk),*(mdl->adb),*(mdl->aub),*(mdl->jq),*(mdl->irow),
          *(mdl->nzl),*(mdl->alpha),*(mdl->fextini),*(mdl->fini),*(mdl->islavnode),
          *(mdl->nslavnode),*(mdl->mortar),*(mdl->ntie),*(mdl->mi),*(mdl->nzs),
          mdl->nasym,mdl->idamping,*(mdl->veold),*(mdl->adc),*(mdl->auc),
          *(mdl->cvini),*(mdl->cv),mdl->alpham,mdl->num_cpus);
}

/* Scratch and evaluate, without the reduce.  Two sites want the model
   evaluated at a trial state and then read the STATE - the stored
   tractions, the damage - rather than the residual, and both had written
   these eight lines out by hand. */
void trial_evaluate(const trialctx *mdl)
{
  ITG *nk=*(mdl->nk),*ne=*(mdl->ne),*mi=*(mdl->mi);
  ITG mt=mi[1]+1,isiz;

  /* The scratch arrays belong to this operation.  The hand-written copies
     freed them BEFORE filling b and allocated after; freeing here, i.e.
     after the caller has filled b, is the only ordering difference and it
     is the safe direction - the arrays are alive for longer, never for
     less. */
  SFREE(*(mdl->v)); SFREE(*(mdl->stx)); SFREE(*(mdl->fn));
  MNEW(*(mdl->v),double,mt**nk);
  isiz=mt**nk; cpypardou(*(mdl->v),*(mdl->vold),&isiz,mdl->num_cpus);
  NNEW(*(mdl->stx),double,6*mi[0]**ne);
  MNEW(*(mdl->fn),double,mt**nk);
  if(*(mdl->ne1d2d)==1) NNEW(*(mdl->inum),ITG,*nk);

  trial_results(mdl);

  if(*(mdl->ne1d2d)==1) SFREE(*(mdl->inum));
}

/* The whole atom: scratch, evaluate, reduce.  The caller fills b first -
   the step is the caller's decision and this file takes no view on it. */
void trial_residual(const trialctx *mdl,double *dst)
{
  trial_evaluate(mdl);
  trial_reduce(mdl,dst);
}

/* Every field is the address of one of the caller's locals, so every field
   is non-NULL.  A field TRIAL_BIND forgot is therefore a NULL, and this
   loop is the only thing between that and a wrong answer three hundred
   increments later.  It walks the struct as an array of pointers, which is
   legal here because every member is a pointer and C guarantees no padding
   between members of identical type. */
ITG trial_check(const trialctx *mdl)
{
  const void *const *p=(const void *const *)mdl;
  ITG n=(ITG)(sizeof(trialctx)/sizeof(void *)),i,nbad=0;

  for(i=0;i<n;i++) if(p[i]==NULL) nbad++;
  printf("[TRIAL] context bound: %" ITGFORMAT " field(s), %" ITGFORMAT
         " unbound -- %s\n",n,nbad,nbad?"FAILED":"PASSED");
  fflush(stdout);
  return nbad;
}
