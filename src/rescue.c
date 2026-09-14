/*     CalculiX - damage/fracture extension                              */
/*     rescue.c: the state of what happens when an increment fails.      */

/* Why this module exists
   ----------------------
   Six clusters of locals in nonlingeo() - sixty-seven of them - are one
   mechanism seen from six angles: the damage line search, its probe,
   transactional backtracking, the same-load re-equilibration that follows
   a deletion, the bounded recovery corridor, and the levels that order
   them.  Together they are the answer to "what does this solver do when an
   increment will not converge", and that answer had no place to live.

   The fields keep their cluster's name on purpose.  Four of the six had a
   `mode'; flattening the names would have merged four distinct variables
   into one field, compiled cleanly, and changed the answer.

   Contract
   --------
     - rescue is the state: what each level is configured to do, what it
       has spent, and the snapshots it restores from.  One object, one
       lifetime, initialised once;
     - the LADDER - which level fires, in what order, and what a failure at
       the top does - is still nonlingeo's, because it cuts increments,
       rolls back topology and ends steps.  Moving the state is what makes
       moving the ladder possible later; it is not the same thing and this
       file does not claim to be it.

   The measured reason any of this exists is in globalize.c: six
   mechanisms stacked in a fixed order with no interface, three
   consecutive attempts producing bit-identical residual sequences, and two
   rescue levels that had run and changed nothing.                      */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#include "ccxopt.h"
#include "ccxfork.h"
/* The handshake between the rescue ladder and the convergence verdict.
   checkconvergence.c SETS ccx_rescue_req when it is about to stop a run
   that a rescue could still save, and reads ccx_rescue_arm to know whether
   one is available; nonlingeo() arms and consumes.  Three globals, because
   the verdict is reached inside a routine that this object cannot be
   passed into - that is the honest description, not a design.  They are
   defined HERE rather than in nonlingeo.c so that the state and the flags
   that publish it are in one file. */
ITG ccx_rescue_active=0,ccx_rescue_arm=0,ccx_rescue_req=0;

/* The values the sixty-seven locals carried at their declarations. */
void rescue_init(rescue *r)
{
  memset(r,0,sizeof(*r));
  r->bt_dam=NULL; r->bt_visc=NULL; r->bt_xs=NULL;
  r->linesearch_env=NULL; r->linesearch_step=NULL; r->reeq_scale_env=NULL;
  r->bt_floor=0.015625;
  r->bt_growth=1.;
  r->bt_window=1;
  r->corr_exit=5;
  r->corr_grace=10;
  r->corr_maxesc=3;
  r->corr_maxinc=60;
  r->corr_maxwall=30;
  r->corr_stableneed=3;
  r->corr_tryevery=5;
  r->corr_minfrac=0.05;
  r->ls_trials=DAMAGE_LINESEARCH_MAX_TRIALS;
  r->ls_min=DAMAGE_LINESEARCH_MIN;
  r->rec_maxunrec=3;
  r->rec_window=5;
  r->rescue_maxlevel=1;
}

/* ---- transactional backtracking -----------------------------------------

   The second loop to move, and it moved for the same two reasons as the
   first: its state has an owner here, and the residual evaluations it
   makes go through trial.c.  148 lines naming twenty-seven things outside
   themselves, of which nineteen now come through trialctx.

   What it does: on a re-equilibration step that makes the residual worse,
   walk alpha down 1, 1/2 ... 1/64 with a committed baseline restored
   before every probe, accept on Armijo against a non-monotone reference,
   and restore the full step if nothing is acceptable.  THIS CHANGES THE
   ANSWER, which is why it is off by default and why the guard stays with
   the caller.

   stx is read back after trial_residual() has reallocated it, so it is
   dereferenced from the context at each use rather than aliased.  v, fn
   and inum are not named here at all any more: the evaluation that used to
   open-code the scratch dance is trial_evaluate().                     */

void rescue_backtrack(rescue *r,const trialctx *mdl,glob_census *g,
                      probedrv *p,double *damvisc,ITG iit)
{
  double *b=*(mdl->b),*dam=*(mdl->dam),*xstate=*(mdl->xstate);
  ITG *neq=*(mdl->neq),*mi=*(mdl->mi),*ne=*(mdl->ne);
  ITG *nstate_=*(mdl->nstate_);
  ITG ne0=*(mdl->ne0),iinc=*(mdl->iinc),num_cpus=*(mdl->num_cpus),isiz;

  /* gcc warns that bref `may be used uninitialized' here and did NOT warn
     about the identical lines inside nonlingeo() - checked, by compiling
     the pre-move file.  It is a false positive that only became reachable
     when the function got small enough to analyse: bii runs from 0, the
     only early break is guarded on bii>=2, and the bii==0 pass assigns
     bref before any pass can read it.  The code is unchanged. */

  static const double btA[]={0.,1.,0.5,0.25,0.125,0.0625,0.03125,
                             0.015625};
  ITG bii,bjj,bacc,bnst,bevt;
  double ba,br,bref,br2;

  bnst=*nstate_;
  if(p->ray_p==NULL){
    NNEW(p->ray_p,double,neq[1]);
    NNEW(p->ray_res,double,neq[1]);
  }
  if(r->bt_dam==NULL){
    NNEW(r->bt_dam,double,mi[0]**ne);
    NNEW(r->bt_visc,double,mi[0]**ne);
    if(bnst>0) NNEW(r->bt_xs,double,bnst*mi[0]**ne);
  }
  isiz=neq[1];cpypardou(p->ray_p,b,&isiz,&num_cpus);
  isiz=mi[0]**ne;cpypardou(r->bt_dam,dam,&isiz,&num_cpus);
  if(damvisc!=NULL){
    isiz=mi[0]**ne;
    cpypardou(r->bt_visc,damvisc,&isiz,&num_cpus);
  }
  if((bnst>0)&&(r->bt_xs!=NULL)){
    isiz=bnst*mi[0]**ne;cpypardou(r->bt_xs,xstate,&isiz,&num_cpus);
  }

  bacc=-1; r->bt_ntrial=0;
  /* non-monotone reference: the largest residual over the window.  The
     ring is reset at the first iteration of each same-load solve. */
  if(iit<=1){ r->bt_nring=0; }
  for(bii=0;bii<8;bii++){
    ba=btA[bii];
    if((bii>=2)&&(ba<r->bt_floor-1.e-12)) break;
    isiz=mi[0]**ne;cpypardou(dam,r->bt_dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,r->bt_visc,&isiz,&num_cpus);
    }
    if((bnst>0)&&(r->bt_xs!=NULL)){
      isiz=bnst*mi[0]**ne;
      cpypardou(xstate,r->bt_xs,&isiz,&num_cpus);
    }
    for(bjj=0;bjj<neq[1];bjj++) b[bjj]=ba*p->ray_p[bjj];
    trial_residual(mdl,p->ray_res);
    br=0.;br2=0.;
    for(bjj=0;bjj<neq[0];bjj++){
      if(fabs(p->ray_res[bjj])>br) br=fabs(p->ray_res[bjj]);
      br2+=p->ray_res[bjj]*p->ray_res[bjj];
    }
    br2=sqrt(br2);

    if((p->evt_on==1)&&(bii==0)){
      if(p->evt_sgn==NULL) NNEW(p->evt_sgn,ITG,mi[0]*ne0);
      damage_evt_sign(mdl,p->evt_sgn);
      for(bjj=0;bjj<8;bjj++){
        p->evt_nsw[bjj]=0;p->evt_fe[bjj]=0;p->evt_fp[bjj]=0;
      }
    }else if((p->evt_on==1)&&(bii<8)){
      p->evt_nsw[bii]=damage_evt_flips(mdl,p->evt_sgn,
                                           &p->evt_fe[bii],
                                           &p->evt_fp[bii]);
    }
    if(bii==0){
      r->bt_r0=br;
      bref=br;
      for(bjj=0;bjj<r->bt_nring;bjj++)
        if(r->bt_ring[bjj]>bref) bref=r->bt_ring[bjj];
      continue;
    }
    if(bii==1){
      /* engagement gate: a full step that is not much worse than the
         current point is taken unchanged.  Newton is allowed to be
         temporarily worse; that is how the undamped run crosses. */
      if(br<=r->bt_growth*r->bt_r0){ bacc=1; break; }
    }
    r->bt_ntrial++;
    printf("[DAMAGE BT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
           " trial=%" ITGFORMAT " alpha=%.6f R=%.6e R0=%.6e"
           " Rref=%.6e R2=%.6e%s",
           iinc,iit,r->bt_ntrial,ba,br,r->bt_r0,bref,br2,"\n");
    if(br<=(1.-1.e-4*ba)*bref){ bacc=bii; break; }
  }
  /* push the residual actually kept into the non-monotone ring */
  if(r->bt_window>1){
    if(r->bt_nring<r->bt_window){
      r->bt_ring[r->bt_nring++]=r->bt_r0;
    }else{
      for(bjj=1;bjj<r->bt_window;bjj++)
        r->bt_ring[bjj-1]=r->bt_ring[bjj];
      r->bt_ring[r->bt_window-1]=r->bt_r0;
    }
  }

  /* [DAMAGE EVT] Nothing satisfied Armijo.  Level one restores the
     full step here and hands the increment to the standard cutback.
     Level two first asks whether the obstruction is a UC6 crossing: if
     some ladder alpha changes the compression set, take the SMALLEST
     such alpha.  That leaves the facet definitely on its new branch, so
     the next assembly - built from vold, e_c3d_uc6.f:45 - carries
     ctan(1,1)=kn for it, which is the whole point.  If no alpha changes
     the set, nothing happens and level one's behaviour stands. */
  bevt=-1;
  if((bacc<0)&&(p->evt_on==1)){
    for(bii=7;bii>=1;bii--){
      if(p->evt_nsw[bii]>0){ bevt=bii; break; }
    }
    if(bevt>=0){
      bacc=bevt;
      p->evt_nstep++;
      printf("[DAMAGE EVT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
             " no alpha satisfied Armijo; the UC6 compression set first "
             "changes at alpha=%.6f (%" ITGFORMAT " point(s) cross, "
             "first el=%" ITGFORMAT " ip=%" ITGFORMAT
             "); taking that EVENT STEP so the next tangent is built on "
             "the new branch.  event steps so far: %" ITGFORMAT "%s",
             iinc,iit,btA[bevt],p->evt_nsw[bevt],
             p->evt_fe[bevt],p->evt_fp[bevt],
             p->evt_nstep,"\n");
    }else{
      printf("[DAMAGE EVT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
             " no alpha satisfied Armijo and NO ladder alpha changes the "
             "UC6 compression set; the obstruction is not a crossing, so "
             "level one behaviour stands (full step restored)%s",
             iinc,iit,"\n");
    }
    fflush(stdout);
  }

  ba=(bacc>=0)?btA[bacc]:1.;
  isiz=mi[0]**ne;cpypardou(dam,r->bt_dam,&isiz,&num_cpus);
  if(damvisc!=NULL){
    isiz=mi[0]**ne;
    cpypardou(damvisc,r->bt_visc,&isiz,&num_cpus);
  }
  if((bnst>0)&&(r->bt_xs!=NULL)){
    isiz=bnst*mi[0]**ne;cpypardou(xstate,r->bt_xs,&isiz,&num_cpus);
  }
  for(bjj=0;bjj<neq[1];bjj++) b[bjj]=ba*p->ray_p[bjj];
  trial_evaluate(mdl);
  glob_fired(&*g,GLOB_BACKTRACK);
  printf("[DAMAGE BT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
         " R0=%.6e trials=%" ITGFORMAT " %s%.6f%s",
         iinc,iit,r->bt_r0,r->bt_ntrial,
         (bevt>=0)?"EVENT-STEP alpha=":
         ((bacc>=0)?"ACCEPTED alpha=":
                    "NONE-ACCEPTED-full-step-restored alpha="),
         ba,"\n");
  fflush(stdout);
}

/* ---- arming ------------------------------------------------------------

   Two blocks, kept as two calls at the two points they occupied, because
   the levels block reads what the backtracking block wrote: always-on
   backtracking and the rescue levels are mutually exclusive and the refusal
   is in the second. */

void rescue_configure_backtrack(rescue *r)
{
  const char *e;

  if(ccxopt_getenv("CCX_DAMAGE_REEQ_BACKTRACK")!=NULL){
    r->bt_mode=1;
    /* Three tunables, each aimed at a MEASURED failure of the
       first version (J-15 -> bandrad regressed 25%).
       _GROWTH : engage only when the full step makes the residual
                 worse by more than this factor.  Damping a step
                 that merely fails Armijo is what made the method
                 more aggressive than BK3 (which needs 1.10) and
                 is what stalled bandrad.  1.0 = old behaviour.
       _WINDOW : non-monotone reference (Grippo-Lampariello-
                 Lucidi).  Acceptance compares against the MAX of
                 the last WINDOW residuals, not the current one,
                 so Newton may worsen the residual briefly and
                 cross the kink - which is exactly what the
                 undamped control does.  1 = monotone = old.
       _FLOOR  : refuse to accept a step shorter than this.  The
                 measured death mode was a chain of accepts at
                 alpha=0.031 and 0.016 buying 1-3% each while the
                 iteration budget drained.  0.015625 = old. */
    if((e=ccxopt_getenv("CCX_DAMAGE_BT_GROWTH"))!=NULL){
      r->bt_growth=atof(e);
      if(r->bt_growth<1.) r->bt_growth=1.;
    }
    if((e=ccxopt_getenv("CCX_DAMAGE_BT_WINDOW"))!=NULL){
      r->bt_window=atoi(e);
      if(r->bt_window<1) r->bt_window=1;
      if(r->bt_window>8) r->bt_window=8;
    }
    if((e=ccxopt_getenv("CCX_DAMAGE_BT_FLOOR"))!=NULL){
      r->bt_floor=atof(e);
      if(r->bt_floor<0.015625) r->bt_floor=0.015625;
      if(r->bt_floor>1.) r->bt_floor=1.;
    }
    printf("[DAMAGE BT] transactional backtracking ENABLED in "
           "idamagereeq: alpha 1, 1/2 ... 1/64, Armijo on |R|inf with "
           "c1=1e-4, committed baseline restored before every probe, "
           "full step restored and the increment handed to the standard "
           "cutback if no probe is acceptable.  THIS CHANGES THE "
           "ANSWER.  growth=%.3f window=%" ITGFORMAT
           " floor=%.6f%s",r->bt_growth,r->bt_window,
           r->bt_floor,"\n");
    fflush(stdout);
  }
}

void rescue_configure_levels(rescue *r,loadctl *c,probedrv *p)
{
  const char *e;
  ITG i;                       /* the banner walks the lambda ladder */

  /* ---- CCX_DAMAGE_REEQ_RESCUE ------------------------------------
     Emergency-only backtracking.  Always-on BT is EXPERIMENTAL and was
     measured to shorten solver survival on three placements of four and
     to destroy the bandrad severance the control reaches (J-17), so the
     two must never run together. */

  if((ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE")!=NULL)||
     (ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL)||
     (ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)||
     (ccxopt_getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)){
    r->rescue_mode=1;
    ccx_rescue_active=1;
    if(ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL){
      r->rescue_maxlevel=2;
      p->evt_nstep=0;
    }
    if((ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)||
       (ccxopt_getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)){
      p->evt_nstep=0;
      c->reg_nlam=5;
      r->rescue_maxlevel=2+c->reg_nlam;
    }
    if(ccxopt_getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL){
      r->corr_mode=1;
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_MAXINC"))!=NULL)
        r->corr_maxinc=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_EXIT"))!=NULL)
        r->corr_exit=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_TRY"))!=NULL)
        r->corr_tryevery=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_GRACE"))!=NULL)
        r->corr_grace=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_MAXWALL"))!=NULL)
        r->corr_maxwall=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_MAXESC"))!=NULL)
        r->corr_maxesc=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_STABLE"))!=NULL)
        r->corr_stableneed=atoi(e);
      if((e=ccxopt_getenv("CCX_DAMAGE_CORR_MINFRAC"))!=NULL)
        r->corr_minfrac=atof(e);
      if(r->corr_maxinc<1) r->corr_maxinc=1;
      if(r->corr_exit<1) r->corr_exit=1;
      if(r->corr_tryevery<1) r->corr_tryevery=1;
      if(r->corr_maxwall<1) r->corr_maxwall=1;
      if(r->corr_maxesc<1) r->corr_maxesc=1;
      if(r->corr_stableneed<1) r->corr_stableneed=1;
      printf("[DAMAGE CORR] bounded recovery CORRIDOR enabled.  On a "
             "wall levels 1 and 2 cannot touch (idamagereeq=0) the "
             "regularization that made the increment converge is HELD, "
             "so the next increments start already regularized.  dtime "
             "stays with the stock controller.  Three separate states "
             "are kept: the lambda in use, the PROVEN lambda (one that "
             "survived %" ITGFORMAT " converged increments) and a probe "
             "flag.  Every %" ITGFORMAT " increments the help is probed "
             "downwards (lambda/4, then 0).  A wall on a PROBE returns "
             "to the proven lambda; a wall on the HELD lambda would be "
             "an identical repeat, so lambda is escalated one ladder "
             "step instead, at most %" ITGFORMAT " times, and the "
             "corridor closes if the ladder runs out.  Wall counters "
             "are evaluated AT THE WALL, so a chain of walls cannot "
             "run unbounded.  Exit after %" ITGFORMAT " increments with "
             "NO help.  Breakers: <=%" ITGFORMAT " walls, <=%" ITGFORMAT
             " increments, and after a grace of %" ITGFORMAT " the mean "
             "dtime inside must stay above %.3f of the dtime at the "
             "last clean increment before entry.  On failure the "
             "guarantee is t_end NOT LOWER than rescue-2; a byte-exact "
             "rescue-2 result is impossible once corridor increments "
             "have been accepted, since no entry snapshot is taken.%s",
             r->corr_stableneed,r->corr_tryevery,
             r->corr_maxesc,r->corr_exit,r->corr_maxwall,
             r->corr_maxinc,r->corr_grace,
             r->corr_minfrac,"\n");
      fflush(stdout);
    }
    if((e=ccxopt_getenv("CCX_DAMAGE_RESCUE_WINDOW"))!=NULL){
      r->rec_window=atoi(e);
      if(r->rec_window<1) r->rec_window=1;
    }
    if((e=ccxopt_getenv("CCX_DAMAGE_RESCUE_MAXUNREC"))!=NULL){
      r->rec_maxunrec=atoi(e);
      if(r->rec_maxunrec<1) r->rec_maxunrec=1;
    }
    printf("[DAMAGE RESCUE] bounded recovery window: a rescue counts "
           "as RECOVERED only after %" ITGFORMAT " consecutive "
           "increments converge with no intervention; after %"
           ITGFORMAT " consecutive un-recovered rescues the mechanism "
           "DISARMS itself and the wall goes to the original stock "
           "stop.  This exists because a run that needs rescuing at "
           "nearly every increment is crawling, not passing a wall: "
           "measured, 92 regularized rescues bought 2.2e-4 of step "
           "time on s3rad.%s",r->rec_window,r->rec_maxunrec,"\n");
    if(r->bt_mode==1){
      printf("[DAMAGE RESCUE] CCX_DAMAGE_REEQ_BACKTRACK (always-on, "
             "experimental) must not run together with rescue; it is "
             "switched OFF for this run.%s","\n");
      r->bt_mode=0;
    }
    printf("[DAMAGE RESCUE] emergency rescue backtracking ENABLED.  The "
           "trajectory, the stock Newton and every stock cutback are "
           "unchanged.  Only where the next stock cutback would put "
           "dtheta below tmin and the run would stop, the increment is "
           "rolled back by the STANDARD cutback path and retried ONCE at "
           "the last admissible dtheta with transactional BT active for "
           "that attempt alone.  One attempt per wall; re-armed after any "
           "increment that converges.%s","\n");
    if(r->rescue_maxlevel==2){
      printf("[DAMAGE RESCUE2] second level ARMED.  A wall now gets two "
             "attempts.  The first is the accepted level-one behaviour, "
             "unchanged.  Only if it fails does the second run, and there "
             "the single branch \"nothing accepted -> restore the full "
             "step\" is replaced by an EVENT STEP: the smallest ladder "
             "alpha at which the set of UC6 points in compression differs "
             "from alpha=0, read from sign(stx(1)) over every live UC6 "
             "point.  e_c3d_uc6.f:45 assembles the stiffness from vold, so "
             "the next assembly picks up ctan(1,1)=kn on the crossed facet "
             "by itself.  No constitutive law, no kn, no g and no material "
             "parameter is touched, and no element, ip or increment is "
             "named.  If no ladder alpha changes the set, this level does "
             "nothing.%s","\n");
    }
    if(c->reg_nlam>0){
      printf("[DAMAGE RESCUE3] third level ARMED with %" ITGFORMAT
             " attempt(s): positive diagonal regularization K+lambda*D."
             "  A wall whose failing solve has idamagereeq=0 goes "
             "straight here - it never enters a same-load solve, so "
             "levels 1 and 2 are gated out and would only repeat the "
             "identical attempt.  ad[k] += lambda*D[k] with "
             "D[k]=max(|ad[k]|,1e-6*mean|ad|) > 0, immediately before "
             "the solver dispatch, where the existing stabiliser "
             "already edits the same diagonal.  NOT ad*=(1+lambda), "
             "which shifts only where the diagonal is positive; and D "
             "is NOT plain |ad|, which is zero where the diagonal is "
             "zero and, for ad<0, gives |ad|*(lambda-1) so lambda=1 "
             "lands exactly on zero.  The per-attempt sign census is a "
             "diagnostic of the diagonal, NOT evidence about the "
             "definiteness of K.  Acceptance stays on the UNMODIFIED "
             "residual: only ad is shifted.  If mean|ad| is zero, NaN "
             "or infinite the shift is skipped and the attempt runs "
             "stock.  lambda ladder:",c->reg_nlam);
      for(i=0;i<c->reg_nlam;i++) printf(" %.3e",c->reg_lam[i]);
      printf(".  When it is exhausted the wall is left to the original "
             "stock stop, so the run ends exactly where rescue-2 ends "
             "it.%s","\n");
    }
    fflush(stdout);
  }
}

/* ---- the rescue attempt ------------------------------------------------

   checkconvergence() has just deferred a stop: it left the increment
   exactly as an ordinary cutback leaves it, so the standard rollback
   further down nonlingeo() restores the start of the increment from the
   baselines it already saved.  What this does is undo the step reduction
   and decide which level gets the one extra attempt - corridor, event
   step, dogleg or continuation - counting walls and disarming the ladder
   when it is buying nothing.

   THIRTEEN PARAMETERS, and the count is the point.  Six are the objects
   this reads and writes; the other seven are nonlingeo's own control flow -
   the step time it restores, the iteration and the re-equilibration flag it
   judges on.  Bundling those into a struct would make the signature shorter
   and the coupling identical, and it would put stock CalculiX locals
   (theta, dtheta, iit) into a fork-specific type, which costs something
   real on the next upstream merge.  So they are spelled out: the length of
   this list is the honest measure of how much of the increment this reaches
   into, and it should be read as a number to be reduced, not hidden.

   Only dtheta and dthetaref are written, which is why only those two are
   pointers.                                                            */

void rescue_attempt(rescue *r,dogleg *d,damcont *k,
                    loadctl *c,probedrv *p,glob_census *g,
                    const trialctx *mdl,
                    double *dtheta,double *dthetaref,
                    double theta,const double *tper,ITG iit,ITG idamagereeq)
{
  ITG iinc=*(mdl->iinc),*ipkon=*(mdl->ipkon),*ielprop=*(mdl->ielprop);
  ITG *mi=*(mdl->mi),ne0=*(mdl->ne0);
  double time=*(mdl->time),dtime=*(mdl->dtime),*prop=*(mdl->prop);
  char *lakon=*(mdl->lakon);

if(ccx_rescue_req==1){
  ccx_rescue_req=0;
  ccx_rescue_arm=0;
  if(r->corr_on==1){
    /* Wall inside the corridor.  Counted and judged HERE, not
       after some later converged increment, so a chain of walls
       is bounded even if nothing ever converges again. */
    ITG cwclose=0;
    r->corr_nint++;
    r->corr_nwalltot++;
    if(r->corr_nwalltot>r->corr_maxwall){
      printf("[DAMAGE CORR] CIRCUIT BREAKER (walls): %" ITGFORMAT
             " walls inside the corridor%s",
             r->corr_nwalltot,"\n");
      cwclose=1;
    }else if(r->corr_trial==1){
      /* the probe failed - fall back to the proven lambda.  This
         is not a repeat: lambda changes. */
      r->corr_trial=0;
      r->corr_lam=r->corr_lamstable;
      r->corr_nsince=0;
      printf("[DAMAGE CORR] the probe did not hold at inc=%"
             ITGFORMAT "; back to the proven lambda=%.3e "
             "(intervention %" ITGFORMAT ")%s",
             iinc,r->corr_lam,r->corr_nint,"\n");
    }else{
      /* Wall on the HELD lambda.  Repeating it would recompute
         the identical attempt, so escalate one ladder step -
         bounded - or close. */
      ITG ei;double lnext=0.;
      r->corr_nwallstab++;
      for(ei=0;ei<c->reg_nlam;ei++){
        if(c->reg_lam[ei]>r->corr_lam*1.0000001){
          lnext=c->reg_lam[ei];break;
        }
      }
      if((r->corr_nwallstab>r->corr_maxesc)||(lnext<=0.)){
        printf("[DAMAGE CORR] CIRCUIT BREAKER (escalation): wall "
               "on the held lambda=%.3e at inc=%" ITGFORMAT ", %"
               ITGFORMAT " consecutive, ladder %s%s",
               r->corr_lam,iinc,r->corr_nwallstab,
               (lnext<=0.)?"exhausted":"still open","\n");
        cwclose=1;
      }else{
        printf("[DAMAGE CORR] wall on the held lambda at inc=%"
               ITGFORMAT ": an identical retry is refused; lambda "
               "escalated %.3e -> %.3e (escalation %" ITGFORMAT
               " of %" ITGFORMAT ")%s",
               iinc,r->corr_lam,lnext,r->corr_nwallstab,
               r->corr_maxesc,"\n");
        r->corr_lam=lnext;
        r->corr_lamstable=lnext;
        r->corr_nsince=0;
      }
    }
    if(cwclose==1){
      printf("[DAMAGE CORR] corridor CLOSES: %" ITGFORMAT
             " increments, %" ITGFORMAT " intervention(s), %"
             ITGFORMAT " regularized factorisation(s), step time "
             "gained %.6e.  The mechanism DISARMS; the next wall "
             "goes to the original stock stop%s",
             r->corr_ninc,r->corr_nint,r->corr_nfact,
             (theta-r->corr_theta0)**tper,"\n");
      r->corr_on=0;r->corr_lam=0.;
      r->rec_disarmed=1;
    }
    c->reg_lambda=r->corr_lam;
    c->reg_on=(r->corr_lam>0.)?1:0;
    /* Bank BEFORE zeroing.  Dropping this counts only the
       factorisations of attempts that succeeded, so the price
       of the corridor would be understated by exactly the cost
       of its failures. */
    r->corr_nfact+=c->reg_napply;
    c->reg_napply=0;
    r->rescue_used=0;
    r->rescue_bt_on=0;
    p->evt_on=0;
    r->rec_used_in_inc=1;
    fflush(stdout);
  }else{
  /* Recovery accounting.  A rescue that fires before WINDOW
     untouched increments have gone by has not recovered
     anything; MAXUNREC of those in a row and the mechanism
     stops pretending and disarms. */
  if(r->rescue_used==0){
    if(r->rec_healthy<r->rec_window){
      r->rec_unrec++;
    }else{
      r->rec_unrec=0;
    }
    r->rec_healthy=0;
    if(r->rec_unrec>r->rec_maxunrec){
      r->rec_disarmed=1;
      printf("[DAMAGE RESCUE] recovery window BLOWN: %" ITGFORMAT
             " consecutive rescues without %" ITGFORMAT
             " clean increments in between.  The solver is being "
             "carried, not recovering, so the mechanism DISARMS "
             "and this wall goes to the original stock stop%s",
             r->rec_unrec,r->rec_window,"\n");
      fflush(stdout);
      ccx_rescue_arm=0;
      ccx_rescue_req=0;
    }
  }
  if(r->rec_disarmed==1){
    r->rescue_used=r->rescue_maxlevel;
  }else{
  r->rescue_used++;
  r->rec_used_in_inc=1;
  r->rescue_bt_on=1;
  /* A wall reached with idamagereeq==0 never enters a same-load
     solve, so BT and the event step - both gated on
     idamagereeq==1 - cannot act on it, and levels 1 and 2 would
     recompute the identical attempt.  Measured: s3rad inc=569
     carries four [DAMAGE RESCUE] lines and not one [DAMAGE BT],
     and both retries there failed identically. */
  if((idamagereeq==0)&&(c->reg_nlam>0)&&
     (r->rescue_used<2)){
    printf("[DAMAGE RESCUE] this wall has idamagereeq=0: no "
           "same-load solve, so levels 1 and 2 cannot act on it; "
           "skipping straight to the regularized level%s","\n");
    r->rescue_used=2;
  }
  /* [DAMAGE CT] level 4: Rescue2 levels 1-2 and the dogleg have
     all failed on this wall.  The standard cutback rollback has
     already restored the increment-start state, so everything below
     is measured on the last COMMITTED state.  Any refusal leaves
     that state untouched and hands the wall to the stock stop. */
  if((k->mode==1)&&(r->rescue_used>=4)&&
     (k->on==0)&&(k->refused==0)){
    ITG cse,csi,csn,csr;
    double csm[3],csk,csd,cst;
    if((k->alloc==1)&&(k->nring>=6)&&
       (damcont_select(k->ring,k->fl,k->dt,
                         ipkon,lakon,ielprop,prop,ne0,mi[0],
                         k->head,&cse,&csi,csm,&csk,&csd,
                         &cst,&csn,&csr,k->kaptol,k->bl,k->nbl)==1)){
      k->arm=1;   /* already gated above; this only reports */
      printf("[DAMAGE CT] level 4 pre-check PASSED at inc=%" ITGFORMAT
             ": %" ITGFORMAT " candidate(s) with five admissible "
             "committed intervals; best element %" ITGFORMAT " ip %"
             ITGFORMAT ", kappa=%.6e (max/min=%.4f), ds0=%.6e.  The "
             "increment gets one continuation attempt.%s",
             iinc,csn,cse+1,csi+1,csk,cst,csd,"\n");
    }else{
      printf("[DAMAGE CT] level 4 REFUSED at inc=%" ITGFORMAT
             " (%s; candidates=%" ITGFORMAT ", ring=%" ITGFORMAT
             "/6).  No continuation attempt is granted, so the stock "
             "path is untouched and this wall goes to the ORIGINAL "
             "stock stop.%s",iinc,
             (k->alloc==0)?"ring not allocated":
             ((k->nring<6)?"ring incomplete":
              ((csr==2)?"kappa unstable over the five intervals":
               "no candidate with five admissible intervals")),
             (k->alloc==1)?csn:0,k->nring,"\n");
      k->refused=1;
      r->rescue_used=r->rescue_maxlevel;
      ccx_rescue_arm=0;
    }
    fflush(stdout);
  }
  p->evt_on=(r->rescue_used==2)?1:0;
  c->reg_on=0;
  /* [DAMAGE TR] guard: without c->reg_nlam>0 this block indexes
     c->reg_lam[-1] and switches the diagonal shift on with a
     garbage lambda.  Unreachable while only RESCUE3/CORRIDOR could
     raise the level count; reachable the moment any other mechanism
     claims level 3, which the dogleg does. */
  if((r->rescue_used>=3)&&(c->reg_nlam>0)){
    c->reg_level=r->rescue_used-3;
    if(c->reg_level>=c->reg_nlam)
      c->reg_level=c->reg_nlam-1;
    c->reg_lambda=c->reg_lam[c->reg_level];
    c->reg_on=1;
    c->reg_napply=0;
  }
  }
  }
  /* [DAMAGE TR] level 3.  Levels 1 and 2 have already run and
     failed on THIS wall - unchanged Rescue2 - so the trajectory up
     to here is the Rescue2 trajectory.  Only now is the increment
     given one more attempt, with the trust region choosing every
     correction inside it. */
  d->lasthelp=iinc;
  d->selfrec=0;
  d->on=0;
  if((d->mode==1)&&(r->rescue_used>=3)&&
     (r->rec_disarmed==0)){
    if(d->have<0){
      printf("[DAMAGE TR] not re-arming: the transpose check has "
             "already failed once in this run.  The wall goes to "
             "the original stock stop.%s","\n");
      fflush(stdout);
      r->rescue_used=r->rescue_maxlevel;
    }else if(d->narm>=d->maxarm){
      printf("[DAMAGE TR] ATTEMPT BUDGET EXHAUSTED: %" ITGFORMAT
             " armed attempts of %" ITGFORMAT " used.  The mechanism "
             "DISARMS, the state is left as the standard cutback "
             "rollback restored it, and this wall goes to the "
             "ORIGINAL stock stop.%s",
             d->narm,d->maxarm,"\n");
      fflush(stdout);
      r->rec_disarmed=1;
      r->rescue_used=r->rescue_maxlevel;
    }else{
      if(d->narm>0){
        printf("[DAMAGE TR] inc=%" ITGFORMAT " running total "
               "BEFORE this attempt: "
               "armed %" ITGFORMAT ", accepted steps %" ITGFORMAT
               " (Newton %" ITGFORMAT ", Cauchy %" ITGFORMAT
               ", dogleg %" ITGFORMAT "), iterations that accepted "
               "nothing %" ITGFORMAT ", rejected trials %" ITGFORMAT
               ", residual evaluations %" ITGFORMAT
               ", armed factorisations %" ITGFORMAT "%s",
               iinc,d->narm,d->nacc,d->nnewt,
               d->ncau,d->ndog,d->nfail,
               d->nrej,d->neval,d->nfact,"\n");
        fflush(stdout);
      }
      d->narm++;
      d->on=1;
      d->have=0;
      d->delta=0.;
      d->banner=0;
      d->incarm=iinc;
      printf("[DAMAGE TR] LEVEL 3 armed for inc=%" ITGFORMAT
             " (attempt %" ITGFORMAT " of %" ITGFORMAT ").  Levels 1 "
             "and 2 both failed on this wall; idamagereeq=%" ITGFORMAT
             ".  The increment is retried once more at the last "
             "admissible (*dtheta), and inside it BK3 is replaced by a "
             "dogleg trust region on phi=1/2|R|^2.  Nothing else "
             "changes: no matrix entry, no material constant, no "
             "deletion rule, and convergence is still decided by "
             "checkconvergence on the unmodified residual.%s",
             iinc,d->narm,d->maxarm,idamagereeq,"\n");
      fflush(stdout);
    }
  }
  r->rescue_nfired++;
  (*dtheta)=r->rescue_dtheta_last;
  (*dthetaref)=r->rescue_dthetaref_last;
  glob_fired(&*g,(r->rescue_used>=2)?
             GLOB_RESCUE2:GLOB_RESCUE1);
  printf("[DAMAGE RESCUE] FIRED #%" ITGFORMAT " inc=%" ITGFORMAT
         " iter=%" ITGFORMAT " time=%.12e dtime=%.12e; (*dtheta) restored "
         "to the last admissible %.12e; transactional BT armed for "
         "THIS attempt only (level %" ITGFORMAT ")%s",
         r->rescue_nfired,iinc,iit,time,dtime,(*dtheta),
         r->rescue_used,"\n");
  fflush(stdout);
}
}
