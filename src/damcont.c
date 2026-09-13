/*     CalculiX - damage/fracture extension                              */
/*     damcont.c: the bounded local continuation, and the state it uses.  */

/* Why this module exists
   ----------------------
   [DAMAGE CT] was the second worst row in tools/arch.py: seventeen
   separate places in nonlingeo.c over 13,661 lines, and EIGHTY-EIGHT
   locals - more than any other cluster in the file, more than a sixth of
   everything this fork added to a function that already had five hundred.

   They are one object.  The kinematic snapshot, the bordered system, the
   candidate ring, the step and retry budgets, the counters: created
   together at one arming point, meaningless apart, and every one of them
   already carried the same prefix.

   Contract
   --------
     - damcont is the state: what was frozen, what the bordered system is
       solving, what the budgets are, and what has been spent.  One
       object with one lifetime;
     - damcont_kin() and damcont_snap() READ the mesh.  They mirror
       cohesive_uc6.f exactly and are built from the REFERENCE geometry,
       which is what makes them a mesh constant rather than a state;
     - damcont_bordered() and damcont_rhoden() are the algebra alone - two
       scalars in, one out - which is what lets damcont_selftest() check
       them against a closed form instead of against a run;
     - damcont_select() chooses the control point.  It refuses, with a
       reason code, rather than choosing badly: reason 1 no candidate,
       reason 2 the curvature tolerance not met.

   What this is NOT
   ----------------
   Not a production continuation, and the block comment this file inherits
   says so at length: no terminal landing at lambda=1, no return to stock
   control, no completed step, no restart.  Every ending is PARTIAL.  That
   was true before this file existed and is unchanged by it.

   The DRIVER - arming, the corrector loop, the commit lifecycle and the
   partial exit - is still in nonlingeo(), for the same reason the dogleg
   loop is: it commits topology and hands control on.  What moved is the
   state and the arithmetic.                                            */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

/* The values the eighty-eight locals carried at their declarations. */
void damcont_init(damcont *c)
{
  memset(c,0,sizeof(*c));
  c->beps=NULL; c->dam=NULL; c->jac=NULL; c->qh=NULL; c->r0=NULL;
  c->ring=NULL; c->visc=NULL; c->w=NULL;  c->xs=NULL; c->y=NULL;
  c->z=NULL;    c->bl=NULL;   c->fl=NULL; c->sgn=NULL;
  c->clim=20.;
  c->ulim=20.;
  c->eps=1.e-6;
  c->kaptol=1.5;
  c->rhomin=1.e-4;
  c->maxcorr=15;
  c->maxeval=1600;
  c->maxfact=700;
  c->maxstep=60;
  c->newstep=1;
  c->elem=-1;
  c->ip=-1;
}

/* ==================================================================
   [DAMAGE CT] bounded experimental coupled local continuation (SPEC
   FREEZE v1).  Opt-in, arms only as rescue LEVEL 4 after Rescue2 and the
   dogleg have both failed on one wall.  Everything below is inert unless
   CCX_DAMAGE_CONTINUATION is set.
   ================================================================== */

/* Exact UC6 kinematics, mirroring cohesive_uc6.f:96-160.  rmat rows map a
   global vector to (normal, shear-1, shear-2) and are built from the
   REFERENCE geometry co alone, so they are a mesh constant.  shape is the
   integration-point specific weight set: 1/6 off-point, 2/3 on-point. */

void damcont_kin(const double *co,const ITG *kon,ITG indexe,
                          const double *v,ITG mt,ITG mint,
                          double *dl,double *rmat,double *shape)
{
  double e1[3],e2[3],cv[3],jump[3],n1,nc;
  ITG i,k,nm,np;

  for(k=0;k<3;k++){
    e1[k]=co[3*(kon[indexe+1]-1)+k]-co[3*(kon[indexe]-1)+k];
    e2[k]=co[3*(kon[indexe+2]-1)+k]-co[3*(kon[indexe]-1)+k];
  }
  cv[0]=e1[1]*e2[2]-e1[2]*e2[1];
  cv[1]=e1[2]*e2[0]-e1[0]*e2[2];
  cv[2]=e1[0]*e2[1]-e1[1]*e2[0];
  n1=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
  nc=sqrt(cv[0]*cv[0]+cv[1]*cv[1]+cv[2]*cv[2]);
  if((n1<=1.e-30)||(nc<=1.e-30)){
    for(i=0;i<3;i++){dl[i]=0.;for(k=0;k<3;k++) rmat[3*i+k]=0.;}
    return;
  }
  for(k=0;k<3;k++){rmat[k]=cv[k]/nc;rmat[3+k]=e1[k]/n1;}
  rmat[6]=rmat[1]*rmat[5]-rmat[2]*rmat[4];
  rmat[7]=rmat[2]*rmat[3]-rmat[0]*rmat[5];
  rmat[8]=rmat[0]*rmat[4]-rmat[1]*rmat[3];

  for(i=0;i<3;i++) shape[i]=1./6.;
  shape[mint]=2./3.;
  for(k=0;k<3;k++){
    jump[k]=0.;
    for(i=0;i<3;i++){
      nm=kon[indexe+i]-1;
      np=kon[indexe+i+3]-1;
      jump[k]+=shape[i]*(v[mt*np+k+1]-v[mt*nm+k+1]);
    }
  }
  for(i=0;i<3;i++){
    dl[i]=0.;
    for(k=0;k<3;k++) dl[i]+=rmat[3*i+k]*jump[k];
  }
}

/* One committed endpoint of the ring: the three local separations of every
   live UC6 integration point, plus its failed flag and its compression
   category.  Both flags are needed because the historical admissibility
   guards test them at BOTH endpoints of an interval. */

void damcont_snap(const double *co,const ITG *kon,const ITG *ipkon,
                           const char *lakon,const double *v,
                           const double *stx,const double *xstate,
                           ITG ne0,ITG mi0,ITG nstate,ITG mt,
                           double *ring,ITG *fl)
{
  double dl[3],rmat[9],shape[3];
  ITG i,j,np,idx;

  np=(mi0<3)?mi0:3;
  for(i=0;i<ne0;i++){
    for(j=0;j<mi0;j++){
      idx=mi0*i+j;
      ring[3*idx]=0.;ring[3*idx+1]=0.;ring[3*idx+2]=0.;
      fl[idx]=-1;                                  /* -1: not a live UC6 ip */
    }
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++){
      idx=mi0*i+j;
      damcont_kin(co,kon,ipkon[i],v,mt,j,dl,rmat,shape);
      ring[3*idx]=dl[0];ring[3*idx+1]=dl[1];ring[3*idx+2]=dl[2];
      fl[idx]=((xstate[nstate*idx+3]>=0.5)?1:0)
             |((stx[6*idx]<0.)?2:0);              /* bit0 failed, bit1 comp */
    }
  }
}

/* Predicted reduction is not used here; the continuation needs only the
   bordered algebra.  Isolated so the self-test exercises the same code.

   den = c_lambda + c_u^T y ,  dlambda = (-c - c_u^T z)/den
   Returns 0 on a degenerate or non-finite input - a REFUSAL, never a step
   of zero length.  Every test is written positively so NaN falls through
   to the refusal branch. */

ITG damcont_bordered(double clam,double cuz,double cuy,double c,
                              double *den,double *dlam)
{
  *den=clam+cuy;
  if(!(*den==*den)) return 0;
  if(!(fabs(*den)>0.)) return 0;
  *dlam=(-c-cuz)/(*den);
  if(!(*dlam==*dlam)) return 0;
  return 1;
}

/* Conditioning of the constraint row, measured on the constraint SUPPORT
   only: the raw |den| is not comparable between candidates. */

double damcont_rhoden(double clam,double cunorm,double ysupp,
                               double den)
{
  double d=fabs(clam)+cunorm*ysupp;
  if(!(d>0.)) return 0.;
  return fabs(den)/d;
}

/* [DAMAGE CT] candidate scan.  Uses ONLY committed ring data, so it can be
   evaluated BEFORE the level-4 retry is granted: a refusal then costs the
   stock path nothing at all and the wall goes to the ORIGINAL stock stop
   byte-identically.  Returns 1 and fills the frozen quantities on success,
   0 on refusal.  reason: 1 no candidate, 2 kappa unstable. */

ITG damcont_select(const double *ring,const ITG *fl,
                            const double *dtr,const ITG *ipkon,
                            const char *lakon,const ITG *ielprop,
                            const double *prop,ITG ne0,ITG mi0,ITG head,
                            ITG *belem,ITG *bip,double *m,double *kappa,
                            double *ds0,double *tau,ITG *ncand,ITG *reason,
                            double kaptol,const ITG *bl,ITG nbl)
{
  ITG i,j,k,n,s0,s1,f0,f1,ok,cnt=0;
  /* d0[3] was declared here and never read; gcc says so now that the
     function is compiled on its own.  Removed. */
  double tn0,ts0,gc,beta,df,atau,d1[3],mm[3],deff,adv,best=-1.;
  double q[5],kp[5],sq[5],sk[5],t;

  n=mi0*ne0;*ncand=0;*reason=1;*belem=-1;*bip=-1;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    if(ielprop[i]<0) continue;
    tn0=prop[ielprop[i]+1];ts0=prop[ielprop[i]+2];gc=prop[ielprop[i]+3];
    if((tn0<=0.)||(ts0<=0.)||(gc<=0.)) continue;
    beta=(ts0/tn0)*(ts0/tn0);
    df=2.*gc/tn0;
    atau=1.e-9*df;
    for(j=0;j<3;j++){
      if(j>=mi0) break;
      f1=fl[n*head+mi0*i+j];
      if(f1<0) continue;
      if(f1&1) continue;
      {ITG bk,skip=0;
       for(bk=0;bk<nbl;bk++) if(bl[bk]==4*i+j){skip=1;break;}
       if(skip) continue;}
      for(k=0;k<3;k++) d1[k]=ring[3*(n*head+mi0*i+j)+k];
      deff=(d1[0]>0.?d1[0]*d1[0]:0.)+beta*(d1[1]*d1[1]+d1[2]*d1[2]);
      deff=(deff>0.)?sqrt(deff):0.;
      if(!(deff>0.)) continue;
      if(fabs(d1[0])<0.05*deff) continue;
      mm[0]=(d1[0]>0.?d1[0]:0.)/deff;
      mm[1]=beta*d1[1]/deff;
      mm[2]=beta*d1[2]/deff;
      ok=1;
      for(k=0;k<5;k++){
        s1=(head-k+6)%6;s0=(head-k-1+6)%6;
        f1=fl[n*s1+mi0*i+j];f0=fl[n*s0+mi0*i+j];
        if((f1<0)||(f0<0)){ok=0;break;}
        if((f1&1)||(f0&1)){ok=0;break;}
        if((f1&2)!=(f0&2)){ok=0;break;}
        if((ring[3*(n*s1+mi0*i+j)]>0.)!=(ring[3*(n*s0+mi0*i+j)]>0.)){
          ok=0;break;}
        adv=mm[0]*(ring[3*(n*s1+mi0*i+j)]-ring[3*(n*s0+mi0*i+j)])
           +mm[1]*(ring[3*(n*s1+mi0*i+j)+1]-ring[3*(n*s0+mi0*i+j)+1])
           +mm[2]*(ring[3*(n*s1+mi0*i+j)+2]-ring[3*(n*s0+mi0*i+j)+2]);
        if(!(adv>atau)){ok=0;break;}
        if(!(dtr[s1]>0.)){ok=0;break;}
        q[k]=adv;kp[k]=dtr[s1]/adv;
        if(!(kp[k]>0.)){ok=0;break;}
      }
      if(ok==0) continue;
      cnt++;
      if(q[0]>best){
        ITG b1,b2;
        best=q[0];*belem=i;*bip=j;
        for(k=0;k<3;k++) m[k]=mm[k];
        for(k=0;k<5;k++){sq[k]=q[k];sk[k]=kp[k];}
        for(b1=1;b1<5;b1++){
          t=sq[b1];for(b2=b1;(b2>0)&&(sq[b2-1]>t);b2--) sq[b2]=sq[b2-1];
          sq[b2]=t;
          t=sk[b1];for(b2=b1;(b2>0)&&(sk[b2-1]>t);b2--) sk[b2]=sk[b2-1];
          sk[b2]=t;
        }
        *ds0=sq[2];*kappa=sk[2];
        *tau=(sk[0]>0.)?sk[4]/sk[0]:1.e30;
      }
    }
  }
  *ncand=cnt;
  if((cnt==0)||(*belem<0)){*reason=1;return 0;}
  if(!(*tau<=kaptol)){*reason=2;return 0;}
  *reason=0;
  return 1;
}

/* [DAMAGE CT SELFTEST] bordered algebra on a case with a closed-form
   answer.  Returns the number of failures; a non-zero result stops the run
   before the first increment. */

ITG damcont_selftest(void)
{
  double den,dlam,r;
  ITG k,nbad=0;

  /* c_lambda=0, c_u^T y = 4, c_u^T z = 1, c = 2  ->  den=4, dlam=-3/4 */
  k=damcont_bordered(0.,1.,4.,2.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] den=%.12e (want 4) dlam=%.12e (want -0.75) "
         "rc=%" ITGFORMAT " %s\n",den,dlam,k,
         ((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam+0.75)<1.e-14))?
         "PASS":"FAIL");
  if(!((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam+0.75)<1.e-14))) nbad++;

  /* c_lambda non-zero must enter den */
  k=damcont_bordered(3.,0.,1.,-8.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] den=%.12e (want 4) dlam=%.12e (want 2) "
         "rc=%" ITGFORMAT " %s\n",den,dlam,k,
         ((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam-2.)<1.e-14))?
         "PASS":"FAIL");
  if(!((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam-2.)<1.e-14))) nbad++;

  /* c=0 at a converged constraint must give dlam = -c_u^T z / den */
  k=damcont_bordered(0.,2.,2.,0.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] c=0 -> dlam=%.12e (want -1) %s\n",dlam,
         ((k==1)&&(fabs(dlam+1.)<1.e-14))?"PASS":"FAIL");
  if(!((k==1)&&(fabs(dlam+1.)<1.e-14))) nbad++;

  /* degenerate and non-finite input must REFUSE */
  {
    double z=0.,dn=z/z;
    ITG k1,k2,k3;
    k1=damcont_bordered(0.,1.,0.,1.,&den,&dlam);      /* den=0     */
    k2=damcont_bordered(dn,1.,1.,1.,&den,&dlam);      /* clam=NaN  */
    k3=damcont_bordered(0.,dn,1.,1.,&den,&dlam);      /* cuz=NaN   */
    printf("[DAMAGE CT SELFTEST] refusals: den=0 -> %" ITGFORMAT
           ", clam=NaN -> %" ITGFORMAT ", cuz=NaN -> %" ITGFORMAT
           " (all want 0) %s\n",k1,k2,k3,
           ((k1==0)&&(k2==0)&&(k3==0))?"PASS":"FAIL");
    if(!((k1==0)&&(k2==0)&&(k3==0))) nbad++;
  }

  /* rho_den is a cancellation measure in [0,1] and 0 on a null row */
  r=damcont_rhoden(0.,1.,4.,4.);
  printf("[DAMAGE CT SELFTEST] rho_den=%.12e (want 1) ; null row -> %.12e "
         "(want 0) %s\n",r,damcont_rhoden(0.,0.,0.,0.),
         ((fabs(r-1.)<1.e-14)&&(damcont_rhoden(0.,0.,0.,0.)==0.))?
         "PASS":"FAIL");
  if(!((fabs(r-1.)<1.e-14)&&(damcont_rhoden(0.,0.,0.,0.)==0.))) nbad++;

  printf("[DAMAGE CT SELFTEST] %" ITGFORMAT " failure(s)\n",nbad);
  fflush(stdout);
  return nbad;
}
