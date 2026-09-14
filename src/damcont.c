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

#include "ccxopt.h"
#include "ccxfork.h"
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

void damcont_snap(const trialctx *mdl,const double *v,const double *stx,
                  double *ring,ITG *fl)
{
  /* Unpacked once, so the body below is the body that was there. */
  const double *co=*(mdl->co);
  const ITG *kon=*(mdl->kon),*ipkon=*(mdl->ipkon);
  const char *lakon=*(mdl->lakon);
  const double *xstate=*(mdl->xstate);
  const ITG ne0=*(mdl->ne0),mi0=(*(mdl->mi))[0];
  const ITG nstate=**(mdl->nstate_),mt=(*(mdl->mi))[1]+1;

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

/* ---- arming ------------------------------------------------------------

   Level 4 arms only when levels 1-3 are armed and NOTHING is driving the
   load parameter, because this mechanism takes the boundary over for the
   rest of the step.  That second condition used to be four loose integers
   in this block; it is loadctl_driving() now, and adding a fifth driver is
   one edit in loadctl.c instead of one here and one in dogleg.c.       */

void damcont_configure(damcont *k,dogleg *d,rescue *r,const loadctl *c)
{
  const char *e;

  /* ---- CCX_DAMAGE_CONTINUATION (SPEC FREEZE v1) ------------------

     Bounded EXPERIMENTAL coupled local continuation.  Arms only as
     rescue LEVEL 4, i.e. only after Rescue2 levels 1 and 2 AND the
     dogleg have all failed on one wall.  Its single purpose is to find
     out whether coupled local continuation crosses the s3rad wall near
     inc=589 with physical front advance.

         R(u,lambda) = f(u,xbounact(lambda)) - fext = 0
         c(u,lambda) = m.(delta - delta_c) - ds    = 0

     lambda is a genuine unknown of a bordered system, not a corrected
     theta.  It owns the boundary for the rest of the step once armed.
     This is NOT a production continuation: there is no terminal landing
     at lambda=1, no return to stock control, no completed step and no
     restart.  Every ending is PARTIAL. */

  if(ccxopt_getenv("CCX_DAMAGE_CONTINUATION")!=NULL){
    if((r->rescue_mode==0)||(d->mode==0)){
      printf("*ERROR: CCX_DAMAGE_CONTINUATION requires BOTH "
             "CCX_DAMAGE_REEQ_RESCUE2 and CCX_DAMAGE_TR_DOGLEG; it is a "
             "level ON TOP of them, never a replacement.  Stopping.%s",
             "\n");
      fflush(stdout);FORTRAN(stop,());
    }
    /* the four load-parameter terms are loadctl_driving(); the other two
       are rescue levels, which are rescue.c's to answer for. */
    if((r->corr_mode==1)||(r->bt_mode==1)||loadctl_driving(c)){
      printf("*ERROR: CCX_DAMAGE_CONTINUATION conflicts with "
             "CCX_DAMAGE_ARCLENGTH, CCX_DISSIPATION_CONTROL, "
             "CCX_DAMAGE_PATH, CCX_DAMAGE_REEQ_RESCUE3, "
             "CCX_DAMAGE_RESCUE_CORRIDOR and always-on "
             "CCX_DAMAGE_REEQ_BACKTRACK.  None of them is used as a "
             "foundation and simultaneous operation is refused.  "
             "Stopping.%s","\n");
      fflush(stdout);FORTRAN(stop,());
    }
    k->mode=1;
    r->rescue_maxlevel=4;
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_RHOMIN"))!=NULL)
      k->rhomin=atof(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_CLIM"))!=NULL)
      k->clim=atof(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_ULIM"))!=NULL)
      k->ulim=atof(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_EPS"))!=NULL)
      k->eps=atof(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_KAPTOL"))!=NULL)
      k->kaptol=atof(e);
    if(k->kaptol<1.) k->kaptol=1.5;
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_MAXSTEP"))!=NULL)
      k->maxstep=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_MAXCORR"))!=NULL)
      k->maxcorr=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_MAXFACT"))!=NULL)
      k->maxfact=atoi(e);
    if((e=ccxopt_getenv("CCX_DAMAGE_CT_MAXEVAL"))!=NULL)
      k->maxeval=atoi(e);
    if(k->rhomin<=0.) k->rhomin=1.e-4;
    if(k->clim<=0.) k->clim=20.;
    if(k->ulim<=0.) k->ulim=20.;
    if(k->eps<=0.) k->eps=1.e-6;
    if(k->maxstep<1) k->maxstep=1;
    if(k->maxcorr<1) k->maxcorr=1;
    if(k->maxfact<1) k->maxfact=1;
    if(k->maxeval<1) k->maxeval=1;
    printf("[DAMAGE CT] bounded EXPERIMENTAL continuation armed as rescue "
           "LEVEL 4.  Levels 1-3 (Rescue2 and the dogleg) are unchanged "
           "and run first.  On a wall none of them takes, lambda becomes "
           "a genuine unknown of a bordered system with one FROZEN "
           "local mixed-mode UC6 constraint, and owns the boundary for "
           "the rest of the step.  There is NO terminal landing, NO "
           "return to stock control, NO completed step and NO restart: "
           "every ending is PARTIAL.  Parameters (all opt-in, printed as "
           "actually used): rho_den_min=%.3e C_lambda=%.1f C_u=%.1f "
           "eps_FD=%.3e kappa_tol=%.3f; budget <=%" ITGFORMAT " steps, <=%" ITGFORMAT
           " corrector iterations, <=%" ITGFORMAT " factorisations, <=%"
           ITGFORMAT " residual evaluations.%s",
           k->rhomin,k->clim,k->ulim,k->eps,
           k->kaptol,
           k->maxstep,k->maxcorr,k->maxfact,
           k->maxeval,"\n");
    fflush(stdout);
  }
}

/* ---- the corrector -----------------------------------------------------

   NOTHING IN THIS TREE EXECUTES THIS FUNCTION, and that is worth knowing
   before trusting it.  Level 4 arms on every deck here and is then REFUSED
   by damcont_select() - on the wrapped fast deck the curvature is unstable
   over the five committed intervals, and raising CCX_DAMAGE_CT_KAPTOL only
   moves the refusal to `no candidate with five admissible intervals'.  The
   mechanism was built for the s3rad target and its own banner says every
   ending is PARTIAL.

   So when this body moved out of nonlingeo(), the gate's byte identity
   said nothing about it: the gate does not run it.  What WAS checked is
   the arming and selection path, which does run - the whole [DAMAGE CT]
   output of the wrapped deck, 28 lines, is identical across the move.  The
   body itself is covered by the compiler and by reading, and by nothing
   else.  Whoever next runs a continuation on s3rad is the first person to
   execute this code since it was moved; if it is wrong, it is wrong here
   and not in a thirteen-thousand-line function, which is the only thing
   the move bought.

   One Newton iteration of the bordered system, once level 4 owns the
   boundary: solve the frozen constraint for dlambda, re-solve the
   equilibrium row with the existing factorisation, and report rho and the
   constraint residual rather than assuming either.

   The GUARD stays with the caller - whether the continuation may run at
   all this iteration is a decision about the increment.  The mesh side
   arrives through trialctx; xbounold, the two damage arrays and the three
   solver flags are passed because they are neither results() arguments nor
   state of this object.                                                */

void damcont_corrector(damcont *k,const trialctx *mdl,
                       const double *xboun,const double *xbounold,
                       double *uam,
                       const double *damjac,const double *damvisc,
                       ITG inputformat,ITG nrhs,ITG symmetryflag)
{
  double *b=*(mdl->b),*co=*(mdl->co),*vold=*(mdl->vold),*stx=*(mdl->stx);
  double *qa=mdl->qa,*cam=mdl->cam;
  double *xbounact=*(mdl->xbounact),*dam=*(mdl->dam),*xstate=*(mdl->xstate);
  ITG *neq=*(mdl->neq),*mi=*(mdl->mi),*ne=*(mdl->ne),*nk=*(mdl->nk);
  ITG *kon=*(mdl->kon),*ipkon=*(mdl->ipkon),*nactdof=*(mdl->nactdof);
  ITG *nboun=*(mdl->nboun),*nstate_=*(mdl->nstate_);
  char *lakon=*(mdl->lakon);
  ITG num_cpus=*(mdl->num_cpus),iinc=*(mdl->iinc),ne0=*(mdl->ne0);
  ITG mt=mi[1]+1,isiz;

  ITG ctj,ctk,ctnst,ctbad=0,ctnsw=0;
  k->used=1;
  double ctlam,ctc,ctcuz,ctcuy,ctden,ctdlam,ctrho,ctyinf;
  double ctdl[3],ctrm[9],ctsh[3],ctnrm,ctb0;

  ctnst=*nstate_;
  k->it++;k->ncorr++;
  /* z is what the solver returned in b.  The base snapshot was taken
     BEFORE the solve, at this iterate; it is only restored here. */
  isiz=neq[1];cpypardou(k->z,b,&isiz,&num_cpus);

  /* ---- q by transactional finite difference at the CURRENT u ----
     arc_r0 and the whole read-before-write set were captured BEFORE
     the solve, at this iterate, so every probe starts from ONE base
     state.  The ladder runs on a step's first corrector iteration and
     whenever a probe reports a category switch; the ACCEPTED q is the
     refined q(eps_final). */
  {
    ITG cthalve,ctacc=0,cte,ctp;
    double cteps,ctqn,ctdn,ctbase;
    cteps=k->eps;
    if(k->qh==NULL) NNEW(k->qh,double,neq[1]);

    /* unperturbed baseline through the probe path: it fixes the
       category map AND proves the base state is reproduced */
    ctlam=k->lam;
    for(ctj=0;ctj<neq[1];ctj++) b[ctj]=0.;
    for(ctj=0;ctj<4;ctj++) qa[ctj]=k->qa[ctj];
    for(ctj=0;ctj<5;ctj++) cam[ctj]=k->cam[ctj];
    for(ctj=0;ctj<2;ctj++) uam[ctj]=k->uam[ctj];
    isiz=mi[0]**ne;cpypardou(dam,k->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,k->visc,&isiz,&num_cpus);
    }
    if((ctnst>0)&&(k->xs!=NULL)){
      isiz=ctnst*mi[0]**ne;
      cpypardou(xstate,k->xs,&isiz,&num_cpus);
    }
    if((damjac!=NULL)&&(k->jac!=NULL)){
      isiz=12*mi[0]**ne;
      cpypardou(damjac,k->jac,&isiz,&num_cpus);
    }
    for(ctj=0;ctj<*nboun;ctj++)
      xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
    trial_residual(mdl,k->beps);
    k->neval++;
    damage_evt_sign(stx,ipkon,lakon,ne0,mi[0],k->sgn);
    ctbase=0.;
    for(ctj=0;ctj<neq[1];ctj++)
      if(fabs(k->beps[ctj]-k->r0[ctj])>ctbase)
        ctbase=fabs(k->beps[ctj]-k->r0[ctj]);
    if(k->it<=1){
      printf("[DAMAGE CT FD] inc=%" ITGFORMAT " baseline reproduced: "
             "max|R(base)-arc_r0|=%.6e (must be 0; a non-zero value "
             "means the probe path does not start from the iterate the "
             "residual was taken at)%s",iinc,ctbase,"\n");
      fflush(stdout);
    }

    for(cthalve=0;cthalve<5;cthalve++){
      ctlam=k->lam+cteps;
      for(ctj=0;ctj<neq[1];ctj++) b[ctj]=0.;
    for(ctj=0;ctj<4;ctj++) qa[ctj]=k->qa[ctj];
    for(ctj=0;ctj<5;ctj++) cam[ctj]=k->cam[ctj];
    for(ctj=0;ctj<2;ctj++) uam[ctj]=k->uam[ctj];
    isiz=mi[0]**ne;cpypardou(dam,k->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,k->visc,&isiz,&num_cpus);
    }
    if((ctnst>0)&&(k->xs!=NULL)){
      isiz=ctnst*mi[0]**ne;
      cpypardou(xstate,k->xs,&isiz,&num_cpus);
    }
    if((damjac!=NULL)&&(k->jac!=NULL)){
      isiz=12*mi[0]**ne;
      cpypardou(damjac,k->jac,&isiz,&num_cpus);
    }
    for(ctj=0;ctj<*nboun;ctj++)
      xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
    trial_residual(mdl,k->beps);
    k->neval++;

      ctnsw=damage_evt_flips(stx,ipkon,lakon,ne0,mi[0],k->sgn,
                             &cte,&ctp);
      for(ctj=0;ctj<neq[1];ctj++)
        k->y[ctj]=(k->beps[ctj]-k->r0[ctj])/cteps;
      if(ctnsw!=0){
        printf("[DAMAGE CT FD] inc=%" ITGFORMAT " it=%" ITGFORMAT
               " eps=%.3e REJECTED: %" ITGFORMAT " UC6 category "
               "switch(es) between R(lambda) and R(lambda+eps)%s",
               iinc,k->it,cteps,ctnsw,"\n");
        cteps*=0.5;continue;
      }
      if((k->epsok==1)&&(cthalve==0)){ctacc=1;break;}
      if(cthalve>0){
        ctqn=0.;ctdn=0.;
        for(ctj=0;ctj<neq[1];ctj++){
          ctqn+=k->y[ctj]*k->y[ctj];
          ctdn+=(k->y[ctj]-k->qh[ctj])
               *(k->y[ctj]-k->qh[ctj]);
        }
        ctqn=sqrt(ctqn);ctdn=sqrt(ctdn);
        printf("[DAMAGE CT FD] inc=%" ITGFORMAT " it=%" ITGFORMAT
               " eps=%.3e |q|=%.6e rel.change vs 2eps=%.6e %s%s",
               iinc,k->it,cteps,ctqn,(ctqn>0.)?ctdn/ctqn:0.,
               ((ctqn>0.)&&(ctdn/ctqn<0.1))?
               "ACCEPT (refined q used)":"halve again","\n");
        if((ctqn>0.)&&(ctdn/ctqn<0.1)){
          ctacc=1;k->eps=cteps;k->epsok=1;break;
        }
      }
      isiz=neq[1];cpypardou(k->qh,k->y,&isiz,&num_cpus);
      cteps*=0.5;
    }
    fflush(stdout);
    if(ctacc==0) ctbad=9;
  }
#ifdef PARDISO
  pardiso_solve(k->y,&neq[0],&symmetryflag,&inputformat,&nrhs);
#endif

  /* ---- the constraint at the current iterate ----
     delta_c is anchored HERE, at the first corrector iteration of each
     continuation step, from the same vold the corrector itself reads.
     Anchoring it at arming instead left c off by 8e-10 against
     ds=9.5e-12 on the short deck even though the ring endpoint and
     vold agreed exactly at arming: the state moves between the load
     build and the first corrector call.  With this anchor c = -ds
     holds to machine precision at the start of every step, by
     construction rather than by argument. */
  damcont_kin(co,kon,ipkon[k->elem],vold,mt,k->ip,
                ctdl,ctrm,ctsh);
  if(k->newstep==1){
    k->newstep=0;
    for(ctj=0;ctj<3;ctj++) k->dc[ctj]=ctdl[ctj];
    k->lamc=k->lam;
    k->epsok=0;
  }
  ctc=k->m[0]*(ctdl[0]-k->dc[0])
     +k->m[1]*(ctdl[1]-k->dc[1])
     +k->m[2]*(ctdl[2]-k->dc[2])-k->ds;
  ctcuz=0.;ctcuy=0.;ctnrm=0.;ctyinf=0.;
  for(ctj=0;ctj<k->nsupp;ctj++){
    ctk=k->supp[ctj];
    if(ctk<0) continue;
    ctcuz+=k->w[ctj]*k->z[ctk];
    ctcuy+=k->w[ctj]*k->y[ctk];
    ctnrm+=k->w[ctj]*k->w[ctj];
  }
  ctnrm=sqrt(ctnrm);
  ctb0=0.;
  for(ctj=0;ctj<k->nsupp;ctj++){
    ctk=k->supp[ctj];
    if(ctk<0) continue;
    ctb0+=k->y[ctk]*k->y[ctk];
  }
  ctb0=sqrt(ctb0);
  for(ctj=0;ctj<neq[1];ctj++)
    if(fabs(k->y[ctj])>ctyinf) ctyinf=fabs(k->y[ctj]);

  if(damcont_bordered(k->clam,ctcuz,ctcuy,ctc,
                        &ctden,&ctdlam)==0) ctbad=1;
  ctrho=damcont_rhoden(k->clam,ctnrm,ctb0,ctden);
  k->den=ctden;k->rho=ctrho;k->cprev=ctc;

  if((ctbad==0)&&(ctrho<k->rhomin)) ctbad=2;
  if((ctbad==0)&&
     (fabs(ctdlam)>k->clim*k->lamref)) ctbad=3;
  if((ctbad==0)&&
     (fabs(ctdlam)*ctyinf>k->ulim*k->duref)) ctbad=4;
  if((ctbad==0)&&(ctnsw!=0)) ctbad=5;
  if((ctbad==0)&&(k->neval>=k->maxeval)) ctbad=6;
  if((ctbad==0)&&(k->nfact>=k->maxfact)) ctbad=7;
  if((ctbad==0)&&(k->it>k->maxcorr)) ctbad=8;

  printf("[DAMAGE CT] inc=%" ITGFORMAT " it=%" ITGFORMAT " step=%"
         ITGFORMAT " lambda=%.12e ds=%.6e c=%.6e den=%.6e rho_den=%.4e "
         "dlambda=%.6e |y|inf=%.4e cuz=%.4e cuy=%.4e switches=%"
         ITGFORMAT " evals=%" ITGFORMAT " fact=%" ITGFORMAT " %s%s",
         iinc,k->it,k->step,k->lam,k->ds,
         ctc,ctden,ctrho,ctdlam,ctyinf,ctcuz,ctcuy,ctnsw,
         k->neval,k->nfact,
         (ctbad==0)?"OK":
         ((ctbad==1)?"REFUSE:degenerate-den":
         ((ctbad==2)?"REFUSE:rho_den":
         ((ctbad==3)?"REFUSE:dlambda-bound":
         ((ctbad==4)?"REFUSE:du-bound":
         ((ctbad==5)?"REFUSE:category-switch-in-FD":
         ((ctbad==6)?"REFUSE:eval-budget":
         ((ctbad==7)?"REFUSE:fact-budget":
         ((ctbad==9)?"REFUSE:FD-eps-ladder":
                     "REFUSE:corrector-limit")))))))),
         "\n");
  fflush(stdout);

  if(ctbad!=0){
    /* A refusal attributable to THIS control point blacklists it for
       the epoch and allows the next wall to try the next candidate;
       a budget or FD refusal is not the point's fault and disarms the
       mechanism outright.  The epoch clears on proven progress. */
    if((ctbad==2)||(ctbad==3)||(ctbad==4)){
      if(k->bl==NULL) NNEW(k->bl,ITG,64);
      if(k->nbl<64){
        k->bl[k->nbl++]=4*k->elem+k->ip;
        printf("[DAMAGE CT] control point element %" ITGFORMAT " ip %"
               ITGFORMAT " BLACKLISTED for this epoch (%" ITGFORMAT
               " blacklisted); the next wall may try another "
               "candidate.%s",k->elem+1,k->ip+1,
               k->nbl,"\n");
        k->refused=0;    /* a re-arm is allowed */
      }else{
        k->refused=1;
      }
    }else{
      k->refused=1;
    }
    k->on=0;
    if(k->refused==1) k->partial=1;
    k->lam=k->lamc;
    ctlam=k->lamc;
    for(ctj=0;ctj<neq[1];ctj++) b[ctj]=k->z[ctj];
    for(ctj=0;ctj<4;ctj++) qa[ctj]=k->qa[ctj];
    for(ctj=0;ctj<5;ctj++) cam[ctj]=k->cam[ctj];
    for(ctj=0;ctj<2;ctj++) uam[ctj]=k->uam[ctj];
    isiz=mi[0]**ne;cpypardou(dam,k->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,k->visc,&isiz,&num_cpus);
    }
    if((ctnst>0)&&(k->xs!=NULL)){
      isiz=ctnst*mi[0]**ne;
      cpypardou(xstate,k->xs,&isiz,&num_cpus);
    }
    if((damjac!=NULL)&&(k->jac!=NULL)){
      isiz=12*mi[0]**ne;
      cpypardou(damjac,k->jac,&isiz,&num_cpus);
    }
    for(ctj=0;ctj<*nboun;ctj++)
      xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
    trial_residual(mdl,k->beps);
    k->neval++;

  }else{
    k->lam+=ctdlam;
    ctlam=k->lam;
    for(ctj=0;ctj<neq[1];ctj++)
      b[ctj]=k->z[ctj]+k->y[ctj]*ctdlam;
    for(ctj=0;ctj<4;ctj++) qa[ctj]=k->qa[ctj];
    for(ctj=0;ctj<5;ctj++) cam[ctj]=k->cam[ctj];
    for(ctj=0;ctj<2;ctj++) uam[ctj]=k->uam[ctj];
    isiz=mi[0]**ne;cpypardou(dam,k->dam,&isiz,&num_cpus);
    if(damvisc!=NULL){
      isiz=mi[0]**ne;
      cpypardou(damvisc,k->visc,&isiz,&num_cpus);
    }
    if((ctnst>0)&&(k->xs!=NULL)){
      isiz=ctnst*mi[0]**ne;
      cpypardou(xstate,k->xs,&isiz,&num_cpus);
    }
    if((damjac!=NULL)&&(k->jac!=NULL)){
      isiz=12*mi[0]**ne;
      cpypardou(damjac,k->jac,&isiz,&num_cpus);
    }
    for(ctj=0;ctj<*nboun;ctj++)
      xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
    trial_residual(mdl,k->beps);
    k->neval++;

  }
}
