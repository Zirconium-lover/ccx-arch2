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
