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
