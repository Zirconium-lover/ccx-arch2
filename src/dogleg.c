/*     CalculiX - damage/fracture extension                              */
/*     dogleg.c: the trust region, and the state it works in.            */

/* Why this module exists
   ----------------------
   [DAMAGE TR] was the worst row in tools/arch.py: twenty separate places
   in nonlingeo.c, spread over 15,498 lines, and forty-nine locals - the
   trust radius, the five scalars the step is chosen from, five workspace
   vectors, three snapshot arrays, eight configuration knobs and thirteen
   counters - every one of them reachable and writable from the whole
   function.

   The forty-nine are one object.  They are created together, they mean
   nothing apart, and every one of them is named `damage_dl_' because
   somebody already knew that.

   Contract
   --------
     - dogleg is the state: policy, workspace and census, in that order,
       and the three are separated by comment and not by accident;
     - dogleg_pick() is the step CHOICE and nothing else.  Given a radius
       and five scalars it returns which of NEWTON / CAUCHY / DOGLEG the
       step is and the two coefficients that build it.  It touches no
       vector, allocates nothing and prints nothing, which is what makes
       it testable against a closed form;
     - dogleg_pred() is the model's predicted reduction for that step -
       the denominator of rho.  Same properties;
     - dogleg_selftest() runs both against a J whose answer is known
       exactly, and every refusal branch against a NaN.  It is required by
       the gate.

   What this is NOT
   ----------------
   The trust-region LOOP - accept or reject, grow or shrink, and what to do
   when it runs out - is still in nonlingeo(), because it also commits
   topology and hands control to the next rescue level.  What has changed
   is that it can now be moved: its state has an owner here, and the
   residual evaluations it makes go through trial.c.  Before those two it
   could not be, and that is the whole reason it was where it was.     */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

/* The initial values the forty-nine locals used to carry at their
   declarations.  Written once, here, where the struct is. */
void dogleg_init(dogleg *d)
{
  memset(d,0,sizeof(*d));
  /* memset zeroes the counters and the five scalars; the pointers are set
     explicitly because all-bits-zero is not NULL by the standard, and the
     eight non-zero defaults are the values the declarations carried. */
  d->d=NULL;   d->pn=NULL;  d->pm=NULL; d->r0=NULL; d->res=NULL;
  d->w=NULL;   d->wm=NULL;  d->xs=NULL; d->dam=NULL; d->visc=NULL;
  d->d0fac=1.;
  d->maxarm=12;
  d->maxeval=600;
  d->maxfact=250;
  d->maxtrial=6;
  d->lc_it=1;
  d->lc_nit=1;
  d->recdone=-1;
}

/* ---- [DAMAGE TR] the dogleg step, isolated ------------------------------

   Given the trust radius dl and the five scalars the trust region works in
   (nd2=|d|^2, nw2=|Jd|^2, npn2=|p_N|^2, dtpn=dot(d,p_N), and tc=nd2/nw2 for
   the Cauchy point), return the step as p = (*pa)*d + (*pb)*p_N and its norm.

   Return value: 1 NEWTON (the full step is inside the radius), 2 CAUCHY (the
   Cauchy point is already outside it), 3 DOGLEG (the blend that meets the
   boundary), 0 REFUSE - a degenerate or non-finite input, which the caller
   must treat as "no step", never as a step of zero length.

   Every test is written as a POSITIVE comparison (!(x>0.) rather than x<=0.)
   so that a NaN falls into the refusal branch instead of silently passing. */

ITG dogleg_pick(double dl,double nd2,double nw2,double npn2,
                          double dtpn,double *pa,double *pb,double *nrm)
{
  double tc,pcn,pnn,uu,uv,vv,disc,tau,q;

  *pa=0.;*pb=0.;*nrm=0.;
  if(!(dl>0.)) return 0;
  if(!(nw2>0.)) return 0;
  if(!(nd2>0.)) return 0;
  if(!(npn2>0.)) return 0;
  if(!(dtpn==dtpn)) return 0;

  tc=nd2/nw2;
  if(!(tc>0.)) return 0;
  pnn=sqrt(npn2);
  pcn=tc*sqrt(nd2);
  if(!(pnn>0.)||!(pcn>0.)) return 0;

  if(pnn<=dl){
    *pa=0.;*pb=1.;
  }else if(pcn>=dl){
    *pa=dl/sqrt(nd2);*pb=0.;
  }else{
    uu=pcn*pcn;
    uv=tc*dtpn-uu;
    vv=npn2-2.*tc*dtpn+uu;
    if(!(vv>0.)){
      *pa=0.;*pb=1.;
      q=npn2;
      *nrm=sqrt(q);
      return 1;
    }
    disc=uv*uv-vv*(uu-dl*dl);
    if(!(disc>=0.)) disc=0.;
    tau=(-uv+sqrt(disc))/vv;
    if(!(tau==tau)) return 0;
    if(tau<0.) tau=0.;
    if(tau>1.) tau=1.;
    *pa=(1.-tau)*tc;*pb=tau;
    q=(*pa)*(*pa)*nd2+2.*(*pa)*(*pb)*dtpn+(*pb)*(*pb)*npn2;
    *nrm=(q>0.)?sqrt(q):0.;
    if(!((*nrm)>0.)) return 0;
    return 3;
  }
  q=(*pa)*(*pa)*nd2+2.*(*pa)*(*pb)*dtpn+(*pb)*(*pb)*npn2;
  *nrm=(q>0.)?sqrt(q):0.;
  if(!((*nrm)>0.)) return 0;
  return (*pb>0.)?1:2;
}

/* predicted reduction of phi=1/2|R|^2 for p = pa*d + pb*p_N, using
   J p = pa*w + pb*r0 and R = -r0:  phi(u) - 1/2|R+Jp|^2 */

double dogleg_pred(double pa,double pb,double nb2,double nd2,
                             double nw2)
{
  return 0.5*nb2-0.5*((pb-1.)*(pb-1.)*nb2+2.*pa*(pb-1.)*nd2+pa*pa*nw2);
}

/* ---- [DAMAGE TR] geometry self-test ------------------------------------

   Runs on a J whose answer is known in closed form, so a wrong branch is a
   failure of arithmetic and not of the model:  J=diag(1,10), r0=(1,1).
       p_N = (1, 0.1)      |p_N|^2 = 1.01
       d   = J^T r0 = (1,10)   |d|^2   = 101
       w   = J d    = (1,100)  |w|^2   = 10001
       |r0|^2 = 2  and  dot(d,p_N) = 2   - the same identity the run checks.
   Prints PASS/FAIL per case and returns the number of failures. */

ITG dogleg_selftest(void)
{
  const double nb2=2.,nd2=101.,nw2=10001.,npn2=1.01,dtpn=2.;
  double pa,pb,nr,pr,pnn,pcn,tc;
  ITG k,nbad=0;

  tc=nd2/nw2;pnn=sqrt(npn2);pcn=tc*sqrt(nd2);
  printf("[DAMAGE TR SELFTEST] closed-form case: |p_N|=%.12e |p_C|=%.12e\n",
         pnn,pcn);

  /* 1. the full Newton step lies inside the radius */
  k=dogleg_pick(2.,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=dogleg_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=2 -> kind=%" ITGFORMAT
         " (want 1 NEWTON) pa=%.6e pb=%.6e |p|=%.12e (want %.12e) "
         "pred=%.12e %s\n",k,pa,pb,nr,pnn,pr,
         ((k==1)&&(fabs(nr-pnn)<=1.e-12*pnn)&&(pr>0.))?"PASS":"FAIL");
  if(!((k==1)&&(fabs(nr-pnn)<=1.e-12*pnn)&&(pr>0.))) nbad++;

  /* 2. the Cauchy point is already outside the radius */
  k=dogleg_pick(0.05,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=dogleg_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=0.05 -> kind=%" ITGFORMAT
         " (want 2 CAUCHY) pa=%.6e pb=%.6e |p|=%.12e (want 5.0e-02) "
         "pred=%.12e %s\n",k,pa,pb,nr,pr,
         ((k==2)&&(pb==0.)&&(fabs(nr-0.05)<=1.e-12*0.05)&&(pr>0.))?
         "PASS":"FAIL");
  if(!((k==2)&&(pb==0.)&&(fabs(nr-0.05)<=1.e-12*0.05)&&(pr>0.))) nbad++;

  /* 3. the intermediate step must sit exactly ON the boundary */
  k=dogleg_pick(0.5,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=dogleg_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=0.5 -> kind=%" ITGFORMAT
         " (want 3 DOGLEG) pa=%.6e pb=%.6e |p|=%.12e (want 5.0e-01) "
         "pred=%.12e %s\n",k,pa,pb,nr,pr,
         ((k==3)&&(fabs(nr-0.5)<=1.e-10*0.5)&&(pr>0.)&&(pa>0.)&&(pb>0.))?
         "PASS":"FAIL");
  if(!((k==3)&&(fabs(nr-0.5)<=1.e-10*0.5)&&(pr>0.)&&(pa>0.)&&(pb>0.))) nbad++;

  /* 4. the model prediction is exact for the full Newton step */
  pr=dogleg_pred(0.,1.,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] pred(full Newton)=%.12e (want %.12e, i.e. the "
         "linear model predicts phi=0) %s\n",pr,0.5*nb2,
         (fabs(pr-0.5*nb2)<=1.e-14*nb2)?"PASS":"FAIL");
  if(!(fabs(pr-0.5*nb2)<=1.e-14*nb2)) nbad++;

  /* 5. the Cauchy point is the exact minimiser along d */
  pr=dogleg_pred(tc,0.,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] pred(Cauchy point)=%.12e (want %.12e = "
         "|d|^4/(2|Jd|^2)) %s\n",pr,0.5*nd2*nd2/nw2,
         (fabs(pr-0.5*nd2*nd2/nw2)<=1.e-12*fabs(pr))?"PASS":"FAIL");
  if(!(fabs(pr-0.5*nd2*nd2/nw2)<=1.e-12*fabs(pr))) nbad++;

  /* 6-9. degenerate and non-finite input must REFUSE, not return a step */
  {
    double dnan=0.,dzero=0.;
    ITG k6,k7,k8,k9;
    dnan=dzero/dzero;                       /* NaN without a literal */
    k6=dogleg_pick(0.,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
    k7=dogleg_pick(0.5,nd2,0.,npn2,dtpn,&pa,&pb,&nr);
    k8=dogleg_pick(dnan,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
    k9=dogleg_pick(0.5,nd2,nw2,npn2,dnan,&pa,&pb,&nr);
    printf("[DAMAGE TR SELFTEST] refusals: Delta=0 -> %" ITGFORMAT
           ", |Jd|^2=0 -> %" ITGFORMAT ", Delta=NaN -> %" ITGFORMAT
           ", dot(d,p_N)=NaN -> %" ITGFORMAT " (all want 0) %s\n",
           k6,k7,k8,k9,
           ((k6==0)&&(k7==0)&&(k8==0)&&(k9==0))?"PASS":"FAIL");
    if(!((k6==0)&&(k7==0)&&(k8==0)&&(k9==0))) nbad++;
  }

  printf("[DAMAGE TR SELFTEST] %" ITGFORMAT " failure(s)\n",nbad);
  fflush(stdout);
  return nbad;
}

/* ---- the trust-region loop ---------------------------------------------

   This is the block that could not be moved.  It is 239 lines and it names
   twenty things outside itself; before trial.c and the dogleg object,
   thirteen of those twenty were raw locals of nonlingeo() and the other
   seven were the thirty-line results()/calcresidual() pair written out by
   hand.  There was no signature to give it.

   There is now.  Everything the model side needs comes through trialctx,
   everything the mechanism owns is in dogleg, and what is left is seven
   arguments: the census it reports to, the viscous damage it snapshots,
   and the three numbers it prints.

   The caller keeps the GUARD.  Whether the trust region is allowed to fire
   at all - not thermal, not dynamic, no contact, no continuation running -
   is a decision about the increment, and it belongs where the increment
   is.  What belongs here is what the region then does.

   The local aliases below are exactly the names the block used inside
   nonlingeo(), bound once from the context, so the body that moved is the
   body that was there.  None of the arrays they name is reallocated while
   this runs: trial_residual() replaces v, stx and fn, and nothing else. */

void dogleg_rescue(dogleg *d,const trialctx *t,glob_census *g,
                   double *damvisc,ITG iit,ITG icutb,double *uam)
{
  double *b=*(t->b),*xstate=*(t->xstate),*dam=*(t->dam);
  double *qa=t->qa,*cam=t->cam;
  ITG *neq=*(t->neq),*mi=*(t->mi),*ne=*(t->ne),*nstate_=*(t->nstate_);
  ITG num_cpus=*(t->num_cpus),iinc=*(t->iinc);
  ITG isiz;

  ITG tnst,tii,tjj,tacc,tkind,tbnd,tkkind;
  double tpa,tpb,tnrm,tphi,tl2,tinf,tpred,tared,trho;
  double tpcn,tpnn,tdold;
  double tkpa,tkpb,tkphi,tknrm,tkrho;
  double tqas[4],tuams[2];
  const char *tname;

  tnst=*nstate_;
  if(d->res==NULL){
    NNEW(d->res,double,neq[1]);
    NNEW(d->dam,double,mi[0]**ne);
    NNEW(d->visc,double,mi[0]**ne);
    if(tnst>0) NNEW(d->xs,double,tnst*mi[0]**ne);
  }

  tpnn=sqrt(d->npn2);
  tpcn=d->tc*sqrt(d->nd2);
  d->phi0=0.5*d->nb2;
  tname="NEWTON";tkind=1;tii=0;
  if(d->delta<=0.){
    d->delta=d->d0fac*tpnn;
    d->dmax=1.e3*tpnn;
    if(d->delta<=0.){d->delta=1.;d->dmax=1.e3;}
    printf("[DAMAGE TR] ARMED inc=%" ITGFORMAT " iter=%" ITGFORMAT
           " attempt=%" ITGFORMAT ": |p_N|=%.6e |p_C|=%.6e |R|2=%.6e "
           "phi=%.6e Delta0=%.6e.  BK3 is bypassed for this attempt; "
           "the correction is chosen by the trust region.%s",
           iinc,iit,icutb+1,tpnn,tpcn,sqrt(d->nb2),
           d->phi0,d->delta,"\n");
    fflush(stdout);
  }

  /* snapshot the trial-derived state and the convergence bookkeeping */
  isiz=mi[0]**ne;cpypardou(d->dam,dam,&isiz,&num_cpus);
  if(damvisc!=NULL){
    isiz=mi[0]**ne;
    cpypardou(d->visc,damvisc,&isiz,&num_cpus);
  }
  if((tnst>0)&&(d->xs!=NULL)){
    isiz=tnst*mi[0]**ne;
    cpypardou(d->xs,xstate,&isiz,&num_cpus);
  }
  /* uam is the running maximum correction over the WHOLE step and is
     updated only after this block, so the value captured here is the
     value BEFORE the rejected full Newton step.  Restoring it before
     the final evaluation is what keeps a rejected step out of the
     displacement criterion for the rest of the step. */
  for(tjj=0;tjj<4;tjj++) tqas[tjj]=qa[tjj];
  for(tjj=0;tjj<2;tjj++) tuams[tjj]=uam[tjj];

  tacc=-1;tkpa=0.;tkpb=1.;tkphi=0.;tknrm=tpnn;tkrho=0.;tkkind=1;
  for(tii=0;tii<d->maxtrial;tii++){
    if(d->neval>=d->maxeval) break;
    tdold=d->delta;

    /* ---- the dogleg step for the current radius.  Same function
       the geometry self-test exercised before the run started; a
       degenerate or non-finite input returns 0 and is a refusal, never
       a step of zero length. ---- */
    tkind=dogleg_pick(d->delta,d->nd2,d->nw2,
                         d->npn2,d->dtpn,&tpa,&tpb,&tnrm);
    /* a REFUSAL is not a firing: tkind==0 means no step was
       constructed, and counting it would credit the mechanism for
       declining to act */
    if(tkind!=0) glob_fired(&*g,GLOB_TRUSTREGION);
    if(tkind==0){
      printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
             " trial=%" ITGFORMAT ": the step construction REFUSED a "
             "degenerate or non-finite input (Delta=%.6e |d|^2=%.6e "
             "|Jd|^2=%.6e |p_N|^2=%.6e dot=%.6e).  Nothing is accepted; "
             "the full Newton step is restored and the stock divergence "
             "and cutback machinery takes over unchanged.%s",
             iinc,iit,tii+1,d->delta,d->nd2,
             d->nw2,d->npn2,d->dtpn,"\n");
      fflush(stdout);
      break;
    }
    tname=(tkind==1)?"NEWTON":((tkind==2)?"CAUCHY":"DOGLEG");
    tbnd=(tnrm>=0.99*d->delta)?1:0;
    tpred=dogleg_pred(tpa,tpb,d->nb2,d->nd2,
                         d->nw2);
    if(tpred<=0.){
      printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
             " trial=%" ITGFORMAT " %s Delta=%.6e |p|=%.6e pred=%.6e "
             "<= 0 - the model itself promises nothing, shrinking "
             "without evaluating%s",
             iinc,iit,tii+1,tname,d->delta,tnrm,tpred,"\n");
      fflush(stdout);
      d->delta=0.25*((tnrm>0.)?tnrm:d->delta);
      if(d->delta<1.e-12*d->dmax) break;
      continue;
    }

    /* ---- transactional trial evaluation of the ORIGINAL residual --- */
    /* [DAMAGE TR] cam starts from the STOCK CLEAN state, never from a
       snapshot: cam[0] is a running maximum, and any snapshot taken
       after the full-Newton results() already carries that step's
       correction.  Same initialisation the iteration top uses. */
    for(tjj=0;tjj<3;tjj++) cam[tjj]=0.;
    for(tjj=3;tjj<5;tjj++) cam[tjj]=0.5;
    isiz=mi[0]**ne;cpypardou(dam,d->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,d->visc,&isiz,&num_cpus);
    }
    if((tnst>0)&&(d->xs!=NULL)){
      isiz=tnst*mi[0]**ne;
      cpypardou(xstate,d->xs,&isiz,&num_cpus);
    }
    for(tjj=0;tjj<neq[1];tjj++)
      b[tjj]=tpa*d->d[tjj]+tpb*d->pn[tjj];
    trial_residual(t,d->res);
    d->neval++;
    tphi=0.;tinf=0.;
    for(tjj=0;tjj<neq[1];tjj++){
      tphi+=d->res[tjj]*d->res[tjj];
      if(fabs(d->res[tjj])>tinf) tinf=fabs(d->res[tjj]);
    }
    tl2=sqrt(tphi);tphi*=0.5;

    tared=d->phi0-tphi;
    trho=tared/tpred;
    if(trho>1.e-4){
      tacc=tii;tkpa=tpa;tkpb=tpb;tkphi=tphi;tknrm=tnrm;
      tkrho=trho;tkkind=tkind;
    }
    if(trho<0.25){
      d->delta=0.25*tnrm;
    }else if((trho>0.75)&&(tbnd==1)){
      d->delta=2.*tnrm;
      if(d->delta>d->dmax)
        d->delta=d->dmax;
    }
    printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT " trial=%"
           ITGFORMAT " %s Delta=%.6e |p|=%.6e boundary=%" ITGFORMAT
           " |R|2=%.6e |R|inf=%.6e phi=%.6e pred=%.6e ared=%.6e "
           "rho=%.6e %s Delta %.6e -> %.6e%s",
           iinc,iit,tii+1,tname,tdold,tnrm,tbnd,tl2,tinf,tphi,tpred,
           tared,trho,(trho>1.e-4)?"ACCEPT":"REJECT",tdold,
           d->delta,"\n");
    fflush(stdout);
    if(tacc>=0) break;
    if(d->delta<1.e-12*d->dmax){
      printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
             ": the radius COLLAPSED to %.6e - no step of any length "
             "along either leg reduces phi%s",
             iinc,iit,d->delta,"\n");
      fflush(stdout);
      break;
    }
  }

  /* ---- leave the point that is kept, and prove the transaction ----
     cam/qa/uam go back first, so that only the kept step own
     evaluation contributes to the convergence bookkeeping. */
  for(tjj=0;tjj<4;tjj++) qa[tjj]=tqas[tjj];
  for(tjj=0;tjj<2;tjj++) uam[tjj]=tuams[tjj];
  if(tacc<0){
    tkpa=0.;tkpb=1.;tkkind=0;
    d->nfail++;
    d->nrej+=tii;
  }else{
    d->nacc++;
    d->nrej+=tacc;
    if(tkkind==1) d->nnewt++;
    else if(tkkind==2) d->ncau++;
    else d->ndog++;
  }
  tpa=tkpa;tpb=tkpb;
    /* [DAMAGE TR] cam starts from the STOCK CLEAN state, never from a
       snapshot: cam[0] is a running maximum, and any snapshot taken
       after the full-Newton results() already carries that step's
       correction.  Same initialisation the iteration top uses. */
    for(tjj=0;tjj<3;tjj++) cam[tjj]=0.;
    for(tjj=3;tjj<5;tjj++) cam[tjj]=0.5;
    isiz=mi[0]**ne;cpypardou(dam,d->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,d->visc,&isiz,&num_cpus);
    }
    if((tnst>0)&&(d->xs!=NULL)){
      isiz=tnst*mi[0]**ne;
      cpypardou(xstate,d->xs,&isiz,&num_cpus);
    }
    for(tjj=0;tjj<neq[1];tjj++)
      b[tjj]=tpa*d->d[tjj]+tpb*d->pn[tjj];
    trial_residual(t,d->res);
    d->neval++;
    tphi=0.;tinf=0.;
    for(tjj=0;tjj<neq[1];tjj++){
      tphi+=d->res[tjj]*d->res[tjj];
      if(fabs(d->res[tjj])>tinf) tinf=fabs(d->res[tjj]);
    }
    tl2=sqrt(tphi);tphi*=0.5;

  printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT " KEEP %s "
         "|p|=%.6e phi_kept=%.6e phi_start=%.6e rho=%.6e "
         "Delta_next=%.6e evals=%" ITGFORMAT " armed_factorisations=%"
         ITGFORMAT " rollback_purity=%.3e%s",
         iinc,iit,
         (tkkind==0)?"FULL-NEWTON (nothing accepted; the stock "
         "divergence and cutback machinery takes over unchanged)":
         ((tkkind==1)?"NEWTON":((tkkind==2)?"CAUCHY":"DOGLEG")),
         tknrm,tphi,d->phi0,tkrho,d->delta,
         d->neval,d->nfact,
         ((tacc>=0)&&(tkphi>0.))?fabs(tphi-tkphi)/tkphi:0.,"\n");
  fflush(stdout);
  if((tacc>=0)&&(tkphi>0.)&&(fabs(tphi-tkphi)>1.e-12*tkphi)){
    printf("[DAMAGE TR] *WARNING IMPURITY: re-evaluating the kept step "
           "gave phi=%.17e against %.17e in the trial pass.  The trial "
           "is then NOT a pure function of the step, and every rho "
           "above compared different physical states.%s",
           tphi,tkphi,"\n");
    fflush(stdout);
  }
  if(cam[0]>tknrm*(1.+1.e-9)+1.e-30){
    printf("*ERROR [DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
           ": cam[0]=%.12e exceeds |p_kept|_2=%.12e.  cam[0] is the "
           "largest component of the correction actually applied and "
           "cannot exceed its 2-norm, so the state left behind does "
           "NOT belong to the step that was kept and every "
           "displacement criterion from here on is meaningless.  "
           "Stopping the experiment.%s",iinc,iit,cam[0],tknrm,"\n");
    fflush(stdout);FORTRAN(stop,());
  }
  d->used=1;

  if((d->neval>=d->maxeval)||
     (d->nfact>=d->maxfact)){
    printf("[DAMAGE TR] BUDGET EXHAUSTED (evaluations %" ITGFORMAT
           "/%" ITGFORMAT ", armed factorisations %" ITGFORMAT "/%"
           ITGFORMAT ").  The trust region switches OFF; the rest of "
           "this attempt runs the stock path and the wall, if it "
           "returns, goes to the original stock stop.%s",
           d->neval,d->maxeval,d->nfact,
           d->maxfact,"\n");
    fflush(stdout);
    d->on=0;
  }
}

/* ---- arming ------------------------------------------------------------

   The block that reads this mechanism's own switches, moved beside the
   mechanism.  It stays a separate call at the exact point it used to sit,
   because these arming blocks REFUSE TO ARM ON EACH OTHER'S STATE - the
   dogleg needs Rescue2 already read, the continuation needs the dogleg -
   so the order in which they run is part of the behaviour, not an
   accident.  PETSc has the same shape and the same name for it:
   XXXSetFromOptions, one per object, called in the order the objects are
   created.

   It refuses rather than degrades, and it stops the run rather than
   printing a warning nobody reads: a method whose step construction is
   wrong, or which is sharing a rescue level with a mechanism measured
   negative, produces numbers that look like an answer.                */

void dogleg_configure(dogleg *d,rescue *r,const loadctl *c)
{
  const char *e;

  /* ---- CCX_DAMAGE_TR_DOGLEG ---------------------------------------

     A root-finding TRUST REGION with a dogleg step, on the ORIGINAL
     equilibrium residual.  It is NOT another regularisation: no matrix
     is modified, no diagonal is shifted, no constitutive law, no Kn, no
     g and no deletion criterion is touched.  The only thing that changes
     is HOW LONG and IN WHICH DIRECTION the correction is, inside one
     armed increment attempt.

     Why here and not from the start: J-19 rejected the corridor because
     a held regularisation carried the solver instead of returning it to
     itself.  So this arms ONLY as rescue LEVEL 3, i.e. only after the
     Rescue2 levels 1 and 2 have both failed terminally on the same wall,
     and only on a wall with idamagereeq==0 - where the measured failure
     is a LINE SEARCH failure: at s3rad inc=569 BK3 reports lambda pinned
     at its floor 0.100000 with res_damped 2.001658e-03 ABOVE res_old
     1.942365e-03 on every one of the four identical attempts.  A floor
     of 0.1 on the Newton DIRECTION is exactly what a trust region does
     not have: it may go shorter, and it may leave that direction.

     Model, following PETSc SNESNEWTONTRDC:
         phi(u) = 1/2 |R(u)|^2 ,  R = f_int - f_ext = -b
         J p_N  = -R = b        (p_N is what PARDISO returns)
         g      = J^T R = -d    with d := J^T b
         p_C    = (|d|^2/|Jd|^2) d          (Cauchy point)
         p      = dogleg(p_C,p_N,Delta)
         rho    = [phi(u)-phi(u+p)] / [phi(u) - 1/2|R+Jp|^2]
     Acceptance of the STEP is rho; acceptance of the INCREMENT stays
     with checkconvergence on the unmodified residual.

     J^T IS FORMED AS J^T, and the run proves it: dot(d,p_N) must equal
     |b|^2 exactly, because dot(J^T b, J^-1 b) = b^T b.  That identity is
     printed as TRANSPOSE-CHECK on the first armed iteration and the
     mechanism REFUSES TO ARM if it is not 1 to 1e-8.  The asymmetry of
     the operator is measured at the same point, not assumed. */

  if(ccxopt_getenv("CCX_DAMAGE_TR_DOGLEG")!=NULL){
    if(r->rescue_mode==0){
      printf("*ERROR: CCX_DAMAGE_TR_DOGLEG requires "
             "CCX_DAMAGE_REEQ_RESCUE2; it is a level ON TOP of "
             "Rescue2, not a replacement.  Stopping.%s","\n");
      fflush(stdout);FORTRAN(stop,());
    }
    if((r->corr_mode==1)||(c->reg_nlam>0)){
      printf("*ERROR: CCX_DAMAGE_TR_DOGLEG must not run together with "
             "CCX_DAMAGE_REEQ_RESCUE3 or CCX_DAMAGE_RESCUE_CORRIDOR - "
             "both were measured NEGATIVE (J-19) and both would occupy "
             "the same rescue levels.  Stopping.%s","\n");
      fflush(stdout);FORTRAN(stop,());
    }
    if(r->bt_mode==1){
      printf("*ERROR: CCX_DAMAGE_TR_DOGLEG must not run together with "
             "always-on CCX_DAMAGE_REEQ_BACKTRACK (rejected, J-17).  "
             "Stopping.%s","\n");
      fflush(stdout);FORTRAN(stop,());
    }
    d->mode=1;
    r->rescue_maxlevel=3;
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_MAXTRIAL"))!=NULL)
      d->maxtrial=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_MAXEVAL"))!=NULL)
      d->maxeval=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_MAXFACT"))!=NULL)
      d->maxfact=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_MAXARM"))!=NULL)
      d->maxarm=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_D0"))!=NULL)
      d->d0fac=atof(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_LINCHECK"))!=NULL)
      d->lincheck=atoi(e);
    if(d->lincheck<0) d->lincheck=0;
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_LINCHECK_IT"))!=NULL)
      d->lc_it=atoi(e);
    if(d->lc_it<1) d->lc_it=1;
    if((e=ccxopt_getenv("CCX_DAMAGE_TR_LINCHECK_NIT"))!=NULL)
      d->lc_nit=atoi(e);
    if(d->lc_nit<1) d->lc_nit=1;
    /* the geometry of the step is proved before the first increment, on
       a J whose dogleg is known in closed form.  A failure here is
       arithmetic, so the run must not start. */
    if(dogleg_selftest()!=0){
      printf("*ERROR: the trust-region geometry self-test FAILED.  "
             "Stopping rather than running a method whose step "
             "construction is wrong.%s","\n");
      fflush(stdout);FORTRAN(stop,());
    }
    if(d->maxtrial<1) d->maxtrial=1;
    if(d->maxtrial>12) d->maxtrial=12;
    if(d->maxeval<1) d->maxeval=1;
    if(d->maxfact<1) d->maxfact=1;
    if(d->maxarm<1) d->maxarm=1;
    if(d->d0fac<=0.) d->d0fac=1.;
    printf("[DAMAGE TR] trust-region DOGLEG armed as rescue LEVEL 3.  "
           "Levels 1 and 2 are the unchanged Rescue2 behaviour and run "
           "first; only when BOTH have failed on the same wall does the "
           "increment get one more attempt, and in that attempt every "
           "Newton correction is chosen by a dogleg trust region on "
           "phi=1/2|R|^2 instead of by BK3.  No matrix entry, no "
           "material constant and no deletion rule is touched, and "
           "convergence is still judged by checkconvergence on the "
           "UNMODIFIED residual.  Budget: <=%" ITGFORMAT " trial steps "
           "per iteration, <=%" ITGFORMAT " residual evaluations, <=%"
           ITGFORMAT " armed factorisations, <=%" ITGFORMAT " armed "
           "attempts in the whole run; on exhaustion the state is "
           "restored and the ORIGINAL stock stop runs.  Initial radius "
           "= %.3f * |p_Newton| at the first armed iteration.%s",
           d->maxtrial,d->maxeval,d->maxfact,
           d->maxarm,d->d0fac,"\n");
    fflush(stdout);
  }
}
