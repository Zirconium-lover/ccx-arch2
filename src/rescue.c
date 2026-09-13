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

void rescue_backtrack(rescue *r,const trialctx *t,glob_census *g,
                      probedrv *p,double *damvisc,ITG iit)
{
  double *b=*(t->b),*dam=*(t->dam),*xstate=*(t->xstate);
  ITG *neq=*(t->neq),*mi=*(t->mi),*ne=*(t->ne);
  ITG *nstate_=*(t->nstate_),*ipkon=*(t->ipkon);
  char *lakon=*(t->lakon);
  ITG ne0=*(t->ne0),iinc=*(t->iinc),num_cpus=*(t->num_cpus),isiz;

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
    trial_residual(t,p->ray_res);
    br=0.;br2=0.;
    for(bjj=0;bjj<neq[0];bjj++){
      if(fabs(p->ray_res[bjj])>br) br=fabs(p->ray_res[bjj]);
      br2+=p->ray_res[bjj]*p->ray_res[bjj];
    }
    br2=sqrt(br2);

    if((p->evt_on==1)&&(bii==0)){
      if(p->evt_sgn==NULL) NNEW(p->evt_sgn,ITG,mi[0]*ne0);
      damage_evt_sign((*(t->stx)),ipkon,lakon,ne0,mi[0],p->evt_sgn);
      for(bjj=0;bjj<8;bjj++){
        p->evt_nsw[bjj]=0;p->evt_fe[bjj]=0;p->evt_fp[bjj]=0;
      }
    }else if((p->evt_on==1)&&(bii<8)){
      p->evt_nsw[bii]=damage_evt_flips((*(t->stx)),ipkon,lakon,ne0,mi[0],
                                           p->evt_sgn,
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
  trial_evaluate(t);
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
