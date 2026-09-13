/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2025 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);             */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

/*
  Dissipation-based path following (bordered / arc-length family).

  METHOD AND PROVENANCE
  ---------------------
  The constraint is the energy-release ("dissipation") constraint of

      M.G.D. Gutierrez, "Energy release control for numerical simulations
      of failure in quasi-brittle solids", Communications in Numerical
      Methods in Engineering 20 (2004) 19-29,

  in the incremental form used by

      C.V. Verhoosel, J.J.C. Remmers, M.A. Gutierrez, "A dissipation-based
      arc-length method for robust simulation of brittle and ductile
      failure", International Journal for Numerical Methods in Engineering
      77 (2009) 1290-1321.

  Only the published equations are used.  No third-party source code was
  copied: the one open implementation found during the survey
  (github.com/jfriedlein/energy-based_arc-length_method-dealii_public)
  carries no license at all and therefore may not be reused.  Everything
  below is written against the equations and verified numerically by
  pathfollow_selftest().

  WHY DISSIPATION AND NOT A GEOMETRIC ARC LENGTH
  ---------------------------------------------
  A Crisfield/Riks spherical constraint |du|^2 + psi^2 dlambda^2 = dl^2
  measures the whole displacement vector.  Once damage localises into a
  band of a few elements the norm is dominated by the elastic bulk that is
  UNLOADING, so the constraint stops seeing the mechanism that actually
  drives the branch and the method degenerates exactly where it is needed.
  The dissipation constraint measures the energy the mechanism releases, so
  it stays informative through localisation.  That is the argument in
  Gutierrez (2004) and it is why the dissipation form became the standard
  for cohesive/damage failure.

  FORMULATION
  -----------
  lambda multiplies the prescribed load pattern.  f_hat is the reference
  load vector conjugate to lambda in equation (active-dof) space,

      f_hat := -dR/dlambda ,      R := f_int - f_ext ,

  held FIXED over one increment.  Holding it fixed is what makes the
  constraint an exactly differentiable function of the unknowns, which is
  the property the previous implementation lacked.

  With P(u) := f_hat^T u and the committed state (u_n, lambda_n),
  P_n := f_hat^T u_n, the dissipation released over the increment is

      dG(u,lambda) = 1/2 ( lambda_n * P(u) - lambda * P_n )              (1)

  and the constraint is

      g(u,lambda) = dG(u,lambda) - tau = 0                              (2)

  with tau > 0 the prescribed dissipation increment.  Because f_hat is
  fixed, (1) is bilinear and its derivatives are EXACT:

      dg/du     = 1/2 * lambda_n * f_hat                    =: a        (3)
      dg/dlambda= -1/2 * P_n                                =: bb       (4)

  Sanity of the sign: in the elastic range u ~ lambda*u_1 so P ~ lambda*P_1
  with P_1 = f_hat^T K^-1 f_hat > 0.  On a softening branch u keeps growing
  while lambda falls, hence lambda_n*P > lambda*P_n and dG > 0.  The
  constraint therefore drives the solution FORWARD along the dissipating
  branch whether lambda rises or falls, which is precisely what a snap-back
  needs and what monotone step-time control cannot do.

  EXTENDED SYSTEM
  ---------------
  Newton on (u, lambda) gives the bordered system

      [ K     -f_hat ] [ du     ]   [ -R ]
      [ a^T    bb    ] [ dlambda ] = [ -g ]                             (5)

  solved by the standard two-solve (Riks) decomposition

      du_R = K^-1 (-R)          (the ordinary CalculiX right-hand side)
      du_F = K^-1 f_hat
      du   = du_R + dlambda * du_F

  and the second row then gives the scalar

      dlambda = -( g + a^T du_R ) / ( a^T du_F + bb )
              = -( g + 1/2*lambda_n*fr ) / ( 1/2*lambda_n*ff - 1/2*P_n ) (6)

  with fr := f_hat^T du_R and ff := f_hat^T du_F.

  RELATION TO THE PREVIOUS IMPLEMENTATION IN nonlingeo.c
  ------------------------------------------------------
  The pre-existing CCX_DISSIPATION_CONTROL=2 path measures dG with one
  functional and differentiates a different one:

    * it measures P as the reaction work sum(fn*xboun) over the PRESCRIBED
      dofs, but builds the constraint gradient from f_hat on the FREE dofs.
      d/du of the former is not f_hat, so its row 2 is not the derivative
      of its own row-2 residual;
    * independently of that, its denominator carries +1/2*P_n where the
      linearisation of its own dG expression requires the opposite sign.

  A Newton iteration whose Jacobian row is not the derivative of its
  residual row does not converge to the constraint; it stalls with a small
  correction and a stagnant residual, which is the behaviour recorded in
  the comments there.  Rather than patch signs into that path, this file
  provides one consistent constraint whose derivatives are verified
  numerically, and pathfollow_selftest() fails the build-time check if they
  ever stop matching.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CalculiX.h"

/* ------------------------------------------------------------------ */
/* Pure core.  No globals, no I/O, no solver: everything here is a     */
/* function of its arguments so that it can be checked in isolation.   */
/* ------------------------------------------------------------------ */

/* Dissipation increment, equation (1). */

double pathfollow_dg(double Pn,double lamn,double P,double lam){

  return 0.5*(lamn*P-lam*Pn);
}

/* Explicit constraint gradient, equations (3) and (4), expressed as the
   two scalar contractions the bordered solve needs.  Kept separate from
   pathfollow_dlam so the self test can differentiate (1) numerically and
   compare against these without going through the division. */

void pathfollow_dgrad(double lamn,double Pn,double fr,double ff,
                      double *adur,double *aduf,double *bb){

  *adur=0.5*lamn*fr;      /* a^T du_R */
  *aduf=0.5*lamn*ff;      /* a^T du_F */
  *bb=-0.5*Pn;            /* dg/dlambda */
}

/* Scalar row of the bordered system, equation (6).

   Returns 1 and writes *dlam on success.  Returns 0 and sets *reason on
   refusal, leaving *dlam untouched, so a caller can always fall back to
   the ordinary Newton step without having corrupted anything:

     reason 1 : a non-finite input
     reason 2 : the bordered denominator is singular to working precision,
                i.e. the constraint is locally blind to lambda
     reason 3 : the computed step is not finite
     reason 4 : the step was clipped to dlmax (still a success, reported
                so the caller can log it)

   dlmax<=0 disables clipping. */

ITG pathfollow_dlam(double g,double lamn,double Pn,double fr,double ff,
                    double dlmax,double *dlam,ITG *reason){

  double adur,aduf,bb,den,num,d;

  *reason=0;

  if(!(g==g)||!(lamn==lamn)||!(Pn==Pn)||!(fr==fr)||!(ff==ff)){
    *reason=1;return 0;
  }

  pathfollow_dgrad(lamn,Pn,fr,ff,&adur,&aduf,&bb);

  den=aduf+bb;
  num=-(g+adur);

  /* The scale against which "singular" is judged has to be the size of the
     terms that built den, not an absolute number: den is an energy and its
     natural magnitude changes by orders over a run. */

  d=fabs(aduf)+fabs(bb);
  if(!(fabs(den)>1.e-12*d)||!(d>0.)){
    *reason=2;return 0;
  }

  d=num/den;
  if(!(d==d)){
    *reason=3;return 0;
  }

  if(dlmax>0.){
    if(d>dlmax){d=dlmax;*reason=4;}
    if(d<-dlmax){d=-dlmax;*reason=4;}
  }

  *dlam=d;
  return 1;
}

/* ------------------------------------------------------------------ */
/* Self test.                                                          */
/*                                                                     */
/* Checks, in order:                                                   */
/*   A  dG evaluates to its definition                                 */
/*   B  ROW 2 by finite differences: the analytic gradient (3)/(4)      */
/*      reproduces the directional derivative of (1)                   */
/*   C  ROW 1 and ROW 2 of the full bordered system (5) are both        */
/*      satisfied by the two-solve decomposition, on a dense SPD K      */
/*   D  the same on an INDEFINITE K, i.e. past a limit point where an   */
/*      ordinary Newton step is exactly what fails                      */
/*   E  degenerate inputs are refused rather than propagated            */
/*                                                                     */
/* Returns the number of failures; 0 means every check passed.          */
/* ------------------------------------------------------------------ */

/* Dense LU with partial pivoting, local to the test so that the check
   does not depend on any external solver. */

static ITG pf_lusolve(double *A,double *rhs,ITG n,double *x){

  ITG i,j,k,p;
  double mx,t,*M=NULL,*b=NULL;

  M=(double *)malloc(n*n*sizeof(double));
  b=(double *)malloc(n*sizeof(double));
  if((M==NULL)||(b==NULL)){if(M!=NULL)free(M);if(b!=NULL)free(b);return 0;}

  for(i=0;i<n*n;i++) M[i]=A[i];
  for(i=0;i<n;i++) b[i]=rhs[i];

  for(k=0;k<n;k++){
    p=k;mx=fabs(M[k*n+k]);
    for(i=k+1;i<n;i++){if(fabs(M[i*n+k])>mx){mx=fabs(M[i*n+k]);p=i;}}
    if(!(mx>1.e-300)){free(M);free(b);return 0;}
    if(p!=k){
      for(j=0;j<n;j++){t=M[k*n+j];M[k*n+j]=M[p*n+j];M[p*n+j]=t;}
      t=b[k];b[k]=b[p];b[p]=t;
    }
    for(i=k+1;i<n;i++){
      t=M[i*n+k]/M[k*n+k];
      if(t==0.) continue;
      for(j=k;j<n;j++) M[i*n+j]-=t*M[k*n+j];
      b[i]-=t*b[k];
    }
  }
  for(i=n-1;i>=0;i--){
    t=b[i];
    for(j=i+1;j<n;j++) t-=M[i*n+j]*x[j];
    x[i]=t/M[i*n+i];
  }
  free(M);free(b);
  return 1;
}

/* One bordered-system case: build K, f_hat, R, run the two-solve, and
   verify BOTH rows of (5) directly.  Returns the number of failures. */

static ITG pf_bordered_case(const char *name,double *K,double *fhat,
                            double *R,ITG n,double lamn,double Pn,
                            double tau,double P,double lam){

  ITG i,j,ok,reason,nbad=0;
  double *duR=NULL,*duF=NULL,*du=NULL,*mrhs=NULL;
  double fr=0.,ff=0.,g,dlam=0.,adur,aduf,bb,r1=0.,r2,t,nrm=0.;

  duR=(double *)malloc(n*sizeof(double));
  duF=(double *)malloc(n*sizeof(double));
  du =(double *)malloc(n*sizeof(double));
  mrhs=(double *)malloc(n*sizeof(double));
  if((duR==NULL)||(duF==NULL)||(du==NULL)||(mrhs==NULL)){
    printf("   %-22s ALLOCATION FAILED\n",name);
    if(duR!=NULL) free(duR);
    if(duF!=NULL) free(duF);
    if(du!=NULL) free(du);
    if(mrhs!=NULL) free(mrhs);
    return 1;
  }

  for(i=0;i<n;i++) mrhs[i]=-R[i];
  if(pf_lusolve(K,mrhs,n,duR)==0){nbad++;goto done;}
  if(pf_lusolve(K,fhat,n,duF)==0){nbad++;goto done;}

  for(i=0;i<n;i++){fr+=fhat[i]*duR[i];ff+=fhat[i]*duF[i];}

  g=pathfollow_dg(Pn,lamn,P,lam)-tau;

  ok=pathfollow_dlam(g,lamn,Pn,fr,ff,0.,&dlam,&reason);
  if(ok==0){
    printf("   %-22s REFUSED (reason %" ITGFORMAT ")\n",name,reason);
    nbad++;goto done;
  }

  for(i=0;i<n;i++) du[i]=duR[i]+dlam*duF[i];

  /* row 1:  K*du - f_hat*dlambda + R  ==  0 */

  for(i=0;i<n;i++){
    t=R[i]-fhat[i]*dlam;
    for(j=0;j<n;j++) t+=K[i*n+j]*du[j];
    r1+=t*t;
    nrm+=R[i]*R[i];
  }
  r1=sqrt(r1);nrm=sqrt(nrm)+1.;

  /* row 2:  a^T du + bb*dlambda + g  ==  0 */

  pathfollow_dgrad(lamn,Pn,0.,0.,&adur,&aduf,&bb);
  t=0.;
  for(i=0;i<n;i++) t+=0.5*lamn*fhat[i]*du[i];
  r2=t+bb*dlam+g;

  printf("   %-22s dlambda=%+.6e  row1=%.3e  row2=%.3e  %s\n",
         name,dlam,r1/nrm,fabs(r2),
         ((r1/nrm<1.e-10)&&(fabs(r2)<1.e-10*(fabs(g)+fabs(t)+1.)))?
         "ok":"FAIL");

  if(!(r1/nrm<1.e-10)) nbad++;
  if(!(fabs(r2)<1.e-10*(fabs(g)+fabs(t)+1.))) nbad++;

 done:
  free(duR);free(duF);free(du);free(mrhs);
  return nbad;
}

ITG pathfollow_selftest(void){

  ITG i,j,n=6,nbad=0,ok,reason;
  double K[36],fhat[6],R[6],u[6],du[6];
  double lamn=0.7,Pn=1.3,tau=2.5e-3,lam=0.65,P,eps,dlam=0.;
  double g0,g1,ana,fd,adur,aduf,bb,dnan;

  printf("[PATHFOLLOW] self test\n");

  /* ---- A: dG matches its definition -------------------------------- */

  P=2.1;
  ana=pathfollow_dg(Pn,lamn,P,lam);
  fd=0.5*(lamn*P-lam*Pn);
  printf("   %-22s dG=%+.6e  %s\n","A definition",ana,
         (fabs(ana-fd)<1.e-15)?"ok":"FAIL");
  if(!(fabs(ana-fd)<1.e-15)) nbad++;

  /* ---- B: row 2 gradient by finite differences ---------------------- */
  /* g depends on u only through P=f_hat^T u, so a directional derivative
     along du changes P by f_hat^T du.  Perturb BOTH u and lambda at once,
     which is the combination the bordered row actually applies. */

  for(i=0;i<n;i++){
    fhat[i]=0.3+0.11*(double)i;
    u[i]=0.05*(double)(i+1);
    du[i]=(i%2==0)?0.7:-0.4;
  }
  P=0.;for(i=0;i<n;i++) P+=fhat[i]*u[i];
  dlam=-0.031;                       /* lambda is allowed to DECREASE */

  pathfollow_dgrad(lamn,Pn,0.,0.,&adur,&aduf,&bb);
  ana=0.;for(i=0;i<n;i++) ana+=0.5*lamn*fhat[i]*du[i];
  ana+=bb*dlam;

  eps=1.e-7;
  {
    double Pp=0.;
    for(i=0;i<n;i++) Pp+=fhat[i]*(u[i]+eps*du[i]);
    g0=pathfollow_dg(Pn,lamn,P,lam);
    g1=pathfollow_dg(Pn,lamn,Pp,lam+eps*dlam);
    fd=(g1-g0)/eps;
  }
  printf("   %-22s analytic=%+.10e  fd=%+.10e  rel=%.2e  %s\n",
         "B row2 gradient",ana,fd,fabs(ana-fd)/(fabs(ana)+1.e-30),
         (fabs(ana-fd)<1.e-8*(fabs(ana)+1.))?"ok":"FAIL");
  if(!(fabs(ana-fd)<1.e-8*(fabs(ana)+1.))) nbad++;

  /* ---- C: both rows on an SPD tangent ------------------------------- */

  for(i=0;i<n;i++){
    for(j=0;j<n;j++) K[i*n+j]=(i==j)?(4.+0.5*(double)i):
                       ((abs((int)(i-j))==1)?-1.:0.);
    R[i]=0.02*(double)(i+1)-0.05;
  }
  P=0.;for(i=0;i<n;i++) P+=fhat[i]*u[i];
  nbad+=pf_bordered_case("C bordered SPD",K,fhat,R,n,lamn,Pn,tau,P,lam);

  /* ---- D: both rows on an INDEFINITE tangent ------------------------ */
  /* Past a limit point the tangent has a negative eigenvalue.  An
     ordinary Newton step is meaningless there, but the bordered system
     is still solvable and must satisfy both rows exactly.  This is the
     case the whole method exists for, so it is checked explicitly. */

  K[0]=-3.0;
  nbad+=pf_bordered_case("D bordered indefinite",K,fhat,R,n,lamn,Pn,tau,
                         P,lam);

  /* ---- E: degenerate inputs are refused ----------------------------- */

  dnan=0.;dnan=dnan/((dnan==0.)?dnan:1.);       /* quiet NaN, no literal */

  ok=pathfollow_dlam(dnan,lamn,Pn,0.1,0.2,0.,&dlam,&reason);
  if(!((ok==0)&&(reason==1))){printf("   E nan input NOT refused\n");nbad++;}

  /* lamn=0 and Pn=0 make the whole second row vanish */
  ok=pathfollow_dlam(1.e-3,0.,0.,0.1,0.2,0.,&dlam,&reason);
  if(!((ok==0)&&(reason==2))){
    printf("   E singular row NOT refused (ok=%" ITGFORMAT " reason=%"
           ITGFORMAT ")\n",ok,reason);nbad++;
  }

  /* clipping reports itself but still succeeds */
  dlam=0.;
  ok=pathfollow_dlam(1.,0.5,1.,0.,0.,1.e-3,&dlam,&reason);
  if(!((ok==1)&&(reason==4)&&(fabs(fabs(dlam)-1.e-3)<1.e-15))){
    printf("   E clipping wrong (ok=%" ITGFORMAT " reason=%" ITGFORMAT
           " dlam=%.3e)\n",ok,reason,dlam);nbad++;
  }
  printf("   %-22s %s\n","E degenerate inputs",(nbad==0)?"ok":"see above");

  printf("[PATHFOLLOW] self test %s (%" ITGFORMAT " failure(s))\n",
         (nbad==0)?"PASSED":"FAILED",nbad);
  return nbad;
}

/* ------------------------------------------------------------------ */
/* Diagnostic: show, numerically, that the legacy CCX_DISSIPATION_     */
/* CONTROL=2 scalar row is not the derivative of its own constraint.   */
/*                                                                     */
/* This is not part of the method.  It exists so that the claim "the   */
/* old path stalls because its Jacobian row is inconsistent" is a      */
/* measurement rather than an assertion, and so that the measurement   */
/* can be repeated by anyone.                                          */
/*                                                                     */
/* Legacy row (nonlingeo.c, CCX_DISSIPATION_CONTROL=2):                */
/*     dlam = -(g + 0.5*lamn*fr) / (0.5*lamn*ff + 0.5*Pn)              */
/* Consistent row for the same dG, equation (6):                       */
/*     dlam = -(g + 0.5*lamn*fr) / (0.5*lamn*ff - 0.5*Pn)              */
/*                                                                     */
/* Returns the number of cases in which the legacy step leaves a       */
/* non-zero row-2 residual, i.e. fails to enforce the constraint it    */
/* is supposed to enforce.                                             */
/* ------------------------------------------------------------------ */

ITG pathfollow_legacycheck(void){

  ITG i,n=6,nbad=0,ok,reason,c;
  double fhat[6],duR[6],duF[6];
  double lamn,Pn,tau,lam,P,g,fr,ff,dlam_new=0.,dlam_old,den_old;
  double r2new,r2old,t;

  printf("[PATHFOLLOW] legacy row-2 consistency check\n");
  printf("   case   dlam(consistent)   dlam(legacy)   row2(consistent)"
         "   row2(legacy)\n");

  for(c=0;c<3;c++){

    lamn=0.8-0.25*(double)c;
    Pn=1.4+0.6*(double)c;
    tau=2.0e-3;
    lam=lamn-0.02;

    for(i=0;i<n;i++){
      fhat[i]=0.3+0.11*(double)i;
      duR[i]=0.004*(double)(i+1)*(1.+0.3*(double)c);
      duF[i]=0.09-0.011*(double)i;
    }
    P=0.;fr=0.;ff=0.;
    for(i=0;i<n;i++){
      P+=fhat[i]*(0.05*(double)(i+1));
      fr+=fhat[i]*duR[i];
      ff+=fhat[i]*duF[i];
    }
    g=pathfollow_dg(Pn,lamn,P,lam)-tau;

    ok=pathfollow_dlam(g,lamn,Pn,fr,ff,0.,&dlam_new,&reason);
    if(ok==0){printf("   consistent row refused, reason %" ITGFORMAT "\n",
                     reason);nbad++;continue;}

    den_old=0.5*lamn*ff+0.5*Pn;
    if(!(fabs(den_old)>0.)) continue;
    dlam_old=-(g+0.5*lamn*fr)/den_old;

    /* row 2 residual  a^T du + bb*dlam + g  with du = duR + dlam*duF,
       a = 0.5*lamn*f_hat, bb = -0.5*Pn.  Zero means the step actually
       enforces the constraint. */

    t=0.;for(i=0;i<n;i++) t+=0.5*lamn*fhat[i]*(duR[i]+dlam_new*duF[i]);
    r2new=t-0.5*Pn*dlam_new+g;

    t=0.;for(i=0;i<n;i++) t+=0.5*lamn*fhat[i]*(duR[i]+dlam_old*duF[i]);
    r2old=t-0.5*Pn*dlam_old+g;

    printf("   %-6" ITGFORMAT " %+.8e   %+.8e   %.3e        %.3e\n",
           c,dlam_new,dlam_old,fabs(r2new),fabs(r2old));

    if(!(fabs(r2new)<1.e-12*(fabs(g)+1.))) nbad++;
    if(fabs(r2old)>1.e-12*(fabs(g)+1.)){
      /* expected: the legacy row does NOT satisfy the constraint */
    }else{
      printf("   case %" ITGFORMAT ": legacy row unexpectedly consistent\n",c);
    }
  }
  printf("[PATHFOLLOW] consistent row leaves row2=0; the legacy row does "
         "not.  %" ITGFORMAT " failure(s) in the consistent row.\n",nbad);
  return nbad;
}

/* ------------------------------------------------------------------ */
/* Stateful driver.                                                    */
/*                                                                     */
/* DESIGN NOTE - why nothing here shadows the model.                   */
/*                                                                     */
/* The first version of this driver accumulated the corrections handed  */
/* to results() and used that sum as the trial displacement.  That is   */
/* one assumption too many: it is only right if nothing else in a       */
/* 14000-line Newton loop moves the state.  Measured with               */
/* CCX_PATHFOLLOW_ACCUMCHECK, it is not right - at the first iteration  */
/* of an engaged attempt the model already differed from the sum by     */
/* f_hat^T(vold-vini) = 1.607e-01 while the sum was still exactly zero, */
/* and the rollback probe showed |vold-vini| = 0 at the top of the same */
/* attempt, so the state moves between the increment setup and the      */
/* first iteration.                                                     */
/*                                                                     */
/* Rather than hunt that down and depend on the answer staying true,    */
/* the driver now READS the model: the caller projects f_hat onto       */
/* (vold - vini) and onto (vini - u_ref) and passes the two scalars in. */
/* Whatever else moves the state, the constraint is then evaluated at   */
/* the state that actually exists.                                      */
/*                                                                     */
/* Everything is transactional: pathfollow_commit() is the only writer  */
/* of the committed pair (lambda_n, P_n), and a rejected attempt simply */
/* never reaches it.                                                    */
/* ------------------------------------------------------------------ */

static ITG     pf_on=0;          /* armed                                 */
static ITG     pf_neq=0;
static double  pf_tau=0.;        /* prescribed dissipation increment      */
static double *pf_fh=NULL;       /* f_hat, frozen over one increment      */
static double  pf_lamn=0.;       /* committed load factor                 */
static double  pf_Pn=0.;         /* f_hat^T u at the committed state      */
static ITG     pf_have=0;        /* f_hat captured for this increment     */
static ITG     pf_nref=0;        /* constraint refusals                   */
static double  pf_fflast=0.;     /* f_hat^T K^-1 f_hat, last evaluation   */
static ITG     pf_frozen=0;      /* f_hat is the fixed reference load      */

ITG pathfollow_arm(double tau,ITG neq){

  if(!(tau>0.)||(neq<=0)) return 0;
  pf_neq=neq;
  pf_tau=tau;
  pf_fh=(double *)calloc((size_t)neq,sizeof(double));
  if(pf_fh==NULL) return 0;
  pf_lamn=0.;pf_Pn=0.;pf_have=0;pf_on=1;pf_nref=0;pf_fflast=0.;pf_frozen=0;
  return 1;
}

void pathfollow_disarm(void){

  if(pf_fh!=NULL){free(pf_fh);pf_fh=NULL;}
  pf_on=0;pf_neq=0;pf_have=0;
}

ITG pathfollow_armed(void){return pf_on;}
double pathfollow_lamn(void){return pf_lamn;}
double pathfollow_Pn(void){return pf_Pn;}
ITG pathfollow_refusals(void){return pf_nref;}
ITG pathfollow_have(void){return pf_have;}
const double *pathfollow_fhat(void){return pf_fh;}
double pathfollow_ff(void){return pf_fflast;}

/* The dissipation increment is the step size of this method.  When an
   attempt fails, cutting the STEP TIME achieves nothing - lambda is
   decoupled from it, so the load step is unchanged and the retry is the
   same problem.  What has to shrink is tau. */

void pathfollow_settau(double tau){ if(tau>0.) pf_tau=tau; }
double pathfollow_gettau(void){return pf_tau;}

/* Equation-space dimension changed (remastruct after element deletion).
   f_hat no longer refers to the same dofs, so the constraint origin is
   restarted at the current state: lambda_n is kept, P_n goes to zero.
   The constraint is incremental, so restarting its origin is legitimate;
   carrying a stale vector would not be. */

ITG pathfollow_resize(ITG neq){

  double *a;

  if(pf_on==0) return 0;
  if(neq<=0){pathfollow_disarm();return 0;}
  a=(double *)calloc((size_t)neq,sizeof(double));
  if(a==NULL){pathfollow_disarm();return 0;}
  if(pf_fh!=NULL) free(pf_fh);
  pf_fh=a;pf_neq=neq;pf_Pn=0.;pf_have=0;pf_fflast=0.;pf_frozen=0;
  return 1;
}

/* Start of an attempt.

   f_hat is deliberately NOT discarded here.  It is obtained by dividing
   the first residual of an attempt by the predictor jump, which is a
   finite difference whose step is the predictor itself - and the
   predictor is proportional to tau and shrinks as the branch turns.
   Refreshing f_hat every increment therefore divides a residual that is
   going to zero by a jump that is going to zero, and the scale runs away:
   measured, ff went 1.0 -> 1.6e5 -> 9.8e8 -> 1.8e12 -> 2.8e18 over four
   increments while lambda froze, because each inflated f_hat inflated the
   predictor denominator, which shrank the next jump, which inflated
   f_hat again.

   The scale of f_hat is not free: dG is linear in it, so tau would mean
   something different every increment.  Verhoosel et al. use a FIXED
   reference load vector, and that is what pathfollow_freeze installs -
   f_hat is refreshed only while the constraint is not yet engaged, where
   the jump is the healthy step-time increment, and is frozen from then
   on. */

void pathfollow_incstart(void){ (void)0; }

/* Stop refreshing f_hat; the vector in hand becomes the fixed reference
   load for the rest of the run. */

void pathfollow_freeze(void){ pf_frozen=1; }
ITG pathfollow_frozen(void){ return pf_frozen; }

/* Allow the NEXT capture to refresh f_hat, keeping the current vector as
   the fallback if the jump turns out to be too small to divide by.

   This is safe only where the capture is exact, i.e. where the caller has
   suppressed the displacement extrapolation so that the first residual of
   an attempt differs from the committed one by the prescribed pattern
   alone.  It is not a licence to refresh under the stock predictor: that
   is what produced the runaway recorded in pathfollow_incstart. */

void pathfollow_unfreeze(void){ pf_frozen=0; }

/* First Newton iteration, BEFORE the solve.  b holds fext-f = -R, and the
   only thing that moved since the committed state is the prescribed
   pattern, by dlampred, so -R = f_hat*dlampred with f_hat := -dR/dlambda.
   f_hat is then frozen for the rest of the increment, which is what makes
   the constraint exactly differentiable. */

void pathfollow_capture(const double *b,double dlampred){

  ITG k;

  if(pf_on==0) return;
  if((pf_frozen!=0)&&(pf_have!=0)) return;

  /* The jump has to be big enough for the difference to mean anything.
     1.e-30 was the original guard and it is useless: it admits exactly
     the vanishing jumps that produced the runaway above. */

  if(!(fabs(dlampred)>1.e-8)) return;
  for(k=0;k<pf_neq;k++) pf_fh[k]=b[k]/dlampred;
  pf_have=1;

  /* Freeze on the FIRST capture, which happens in the elastic range.

     There the response to the prescribed jump is linear, so b/dlambda is
     not a secant at all - it is exactly -dR/dlambda, and the reference
     load vector is obtained without approximation.  A capture taken later,
     near the limit point, is a secant of a strongly nonlinear response
     over a finite jump, and it is both badly scaled and tangent-dependent:
     taken at engagement it gave ff = 1.003 and P_n = 0.0222 where the
     elastic relation P_n = lambda_n*ff requires 0.737, i.e. a reference
     vector wrong by a factor of 33.

     This is also what Verhoosel et al. prescribe: one FIXED reference load
     vector for the whole analysis. */

  pf_frozen=1;
}

/* P_n must be re-projected with the f_hat just captured, from the
   committed displacement.  The caller computes it against the model. */

void pathfollow_setPn(double Pn){ if(pf_on!=0) pf_Pn=Pn; }

/* Record ff without touching anything.  Called on every iterate once
   f_hat exists, engaged or not, so the tangent predictor has a stiffness
   the moment the constraint takes over. */

void pathfollow_measure(const double *uf){

  ITG k;
  double ff=0.;

  if((pf_on==0)||(pf_have==0)) return;
  for(k=0;k<pf_neq;k++) ff+=pf_fh[k]*uf[k];
  pf_fflast=ff;
}

/* Tangent predictor.

   THIS is what lets the method turn a limit point, and getting it wrong
   is why the first integration stalled: a predictor that always steps
   lambda FORWARD walks into the region where no equilibrium exists, and
   shrinking it does not help, because the direction is wrong rather than
   the length.

   Apply the constraint to the tangent step itself.  With no correction
   yet the trial displacement is du = dlambda*du_F, so

       dG = 1/2 * dlambda * ( lambda_n*ff - P_n )

   and dG = tau gives

       dlambda = 2*tau / ( lambda_n*ff - P_n ) .                       (7)

   The sign is not imposed, it FALLS OUT: lambda_n*ff-P_n compares the
   tangent compliance with the secant one, so it is positive while the
   structure still stiffens and changes sign exactly at the limit point.
   That is the whole mechanism of a snap-back, in one line, and it is
   available only because f_hat is frozen. */

/* Put lambda exactly on the constraint for the CURRENT displacement.

   CalculiX starts a Newton iteration from an extrapolated state, not from
   the committed one: measured with CCX_PATHFOLLOW_ACCUMCHECK, |vold-vini|
   is 0 at the top of an attempt and 2.7e-03 one statement before the
   Newton loop, every increment from the second on.  The constraint is
   defined incrementally from the COMMITTED state, so it sees that jump as
   dissipation that already happened and opens with a large g - measured
   dG = 2.46e-01 against tau = 2e-03, a factor of 120.

   Correcting that by clipped Newton steps costs an iteration per clip and
   fights the equilibrium iteration.  It is unnecessary: dG is LINEAR in
   lambda at fixed u,

       dG = 1/2 ( lambda_n*P - lambda*P_n ) = tau
   =>  lambda = ( lambda_n*P - 2*tau ) / P_n                           (8)

   so lambda can be placed on the constraint exactly, in closed form,
   before the coupled iteration starts.  Newton then only has to close the
   equilibrium row.

   Returns 0 and leaves *lam alone when P_n is too small to divide by,
   which is the elastic regime where the constraint is degenerate anyway. */

ITG pathfollow_project_lambda(double pdu,double *lam){

  double P,d;

  if((pf_on==0)||(pf_have==0)) return 0;
  P=pf_Pn+pdu;
  if(!(fabs(pf_Pn)>1.e-30*(fabs(P)+1.))) return 0;
  d=(pf_lamn*P-2.*pf_tau)/pf_Pn;
  if(!(d==d)) return 0;
  *lam=d;
  return 1;
}

ITG pathfollow_predictor(double *dlam){

  double den,d;

  if(pf_on==0) return 0;
  den=pf_lamn*pf_fflast-pf_Pn;
  if(!(den==den)) return 0;
  if(!(fabs(den)>1.e-30*(fabs(pf_lamn*pf_fflast)+fabs(pf_Pn)+1.))) return 0;
  d=2.*pf_tau/den;
  if(!(d==d)) return 0;
  *dlam=d;
  return 1;
}

/* After the solve.  b holds du_R, uf holds K^-1 f_hat, and pdu is
   f_hat^T(u_current - u_committed) READ FROM THE MODEL.

   Evaluates the constraint at the current iterate, closes the bordered
   scalar row, applies dlambda*du_F to b and advances lambda.  b is left
   as the full correction the caller hands to results().

   Returns 1 when the constraint step was applied; on refusal b and lam
   are untouched and the caller performs an ordinary Newton step. */

ITG pathfollow_step(double *b,const double *uf,double pdu,double *lam,
                    double dlmax,double *dgout,double *gout,double *dlamout,
                    ITG *reason){

  ITG k,ok;
  double fr=0.,ff=0.,P,dg,g,dlam=0.;

  *reason=0;
  if((pf_on==0)||(pf_have==0)) return 0;

  for(k=0;k<pf_neq;k++){
    fr+=pf_fh[k]*b[k];
    ff+=pf_fh[k]*uf[k];
  }
  pf_fflast=ff;

  P=pf_Pn+pdu;
  dg=pathfollow_dg(pf_Pn,pf_lamn,P,*lam);
  g=dg-pf_tau;

  if(dgout!=NULL) *dgout=dg;
  if(gout!=NULL) *gout=g;

  ok=pathfollow_dlam(g,pf_lamn,pf_Pn,fr,ff,dlmax,&dlam,reason);
  if(ok==0){pf_nref++;return 0;}

  for(k=0;k<pf_neq;k++) b[k]+=dlam*uf[k];
  *lam+=dlam;
  if(dlamout!=NULL) *dlamout=dlam;
  return 1;
}

/* Commit an accepted increment.  pdu is again read from the model.  This
   is the only writer of (lambda_n, P_n). */

void pathfollow_commit(double lam,double pdu,double *dgcommit){

  double P;

  if(pf_on==0) return;
  P=pf_Pn+pdu;
  if(dgcommit!=NULL) *dgcommit=pathfollow_dg(pf_Pn,pf_lamn,P,lam);
  pf_lamn=lam;
  pf_Pn=P;

  /* pf_have is NOT cleared.  Clearing it here is what defeated the freeze
     on the first attempt: capture skips only when f_hat is both frozen
     and present, so a commit that dropped it let the next increment
     re-derive f_hat from a vanishing jump and the scale ran away again
     (ff 1.0 -> 110 -> 3856 -> 8.5e5).  While the constraint is not yet
     engaged, capture refreshes f_hat every increment on its own. */
}

/* ==================================================================== */
/* Crack-opening control.                                               */
/*                                                                      */
/* WHY THIS AND NOT THE DISSIPATION CONSTRAINT.                         */
/*                                                                      */
/* On the clean cohesive benchmark the Gutierrez constraint is          */
/* identically zero.  Measured: P = f_hat^T u tracked lambda*ff to every */
/* printed digit through the whole run (lambda=0.81486456, ff=4.999972e+02,
   P=4.074302e+02 = lambda*ff) and dG stayed at 1e-12 while the interface
   was failing.  The reason is structural, not numerical: under displacement
   control f_hat = -dR/dlambda is supported on the dofs next to the LOADED
   FACE, and once the interface fails the whole block on that side moves
   rigidly with the prescribed face, so f_hat^T u stays proportional to
   lambda no matter what the crack does.  A dissipation constraint built
   from that vector cannot see a crack that is not adjacent to the loading
   boundary.
                                                                        
   Crack-opening control has no such blind spot, and it is the natural
   control for cohesive failure: the opening increases monotonically along
   the whole equilibrium branch, including the part where BOTH the load and
   the end displacement run backwards.  That is exactly a snap-back.
                                                                        
   The control functional is linear in u,
                                                                        
       phi(u) = c^T u  =  mean normal separation over the cohesive facets
                                                                        
   with c built once from the reference geometry, so
                                                                        
       g(u,lambda) = c^T u - phi_target ,  dg/du = c ,  dg/dlambda = 0
                                                                        
   and the second row of the bordered system is
                                                                        
       dlambda = -( g + c^T du_R ) / ( c^T du_F ) .                  (9)
                                                                        
   Both derivatives are exact by construction - c is a constant vector -
   so this constraint cannot suffer the inconsistency that the
   linearisation probe was built to detect.                             */

static double *pf_c=NULL;        /* control functional, equation space   */
static double  pf_target=0.;     /* prescribed opening                   */
static ITG     pf_cod=0;         /* crack-opening control armed          */

ITG pathfollow_cod_arm(const double *c,ITG neq){

  ITG k;

  if((c==NULL)||(neq<=0)) return 0;
  if(pf_c!=NULL) free(pf_c);
  pf_c=(double *)calloc((size_t)neq,sizeof(double));
  if(pf_c==NULL) return 0;
  for(k=0;k<neq;k++) pf_c[k]=c[k];
  pf_cod=1;
  return 1;
}

ITG pathfollow_cod(void){return pf_cod;}
double pathfollow_cod_target(void){return pf_target;}
void pathfollow_cod_settarget(double t){pf_target=t;}
const double *pathfollow_cod_c(void){return pf_c;}

/* Scalar row (9).  cu is c^T u at the current iterate, READ FROM THE
   MODEL by the caller, exactly as for the dissipation constraint. */

ITG pathfollow_cod_step(double *b,const double *uf,double cu,double *lam,
                        double dlmax,double *gout,double *dlamout,
                        ITG *reason){

  ITG k;
  double cdr=0.,cdf=0.,g,dlam=0.,adur,aduf,bb;

  *reason=0;
  if((pf_on==0)||(pf_cod==0)) return 0;

  for(k=0;k<pf_neq;k++){
    cdr+=pf_c[k]*b[k];
    cdf+=pf_c[k]*uf[k];
  }
  g=cu-pf_target;
  if(gout!=NULL) *gout=g;
  if(ccxopt_getenv("CCX_PATHFOLLOW_PROBE")!=NULL){
    printf("[COD] phi=%.6e target=%.6e g=%.4e  c.duR=%.6e c.duF=%.6e\n",
           cu,pf_target,g,cdr,cdf);
    fflush(stdout);
  }

  /* Reuse the same guarded division as the dissipation row: a^T du_R
     plays the role of c^T du_R, a^T du_F that of c^T du_F, and the
     explicit lambda derivative is zero. */

  adur=cdr;aduf=cdf;bb=0.;
  {
    double den=aduf+bb,num=-(g+adur),d,sc;
    if(!(g==g)||!(cdr==cdr)||!(cdf==cdf)){*reason=1;return 0;}
    sc=fabs(aduf)+fabs(bb);
    if(!(fabs(den)>1.e-12*sc)||!(sc>0.)){*reason=2;return 0;}
    d=num/den;
    if(!(d==d)){*reason=3;return 0;}
    if(dlmax>0.){
      if(d>dlmax){d=dlmax;*reason=4;}
      if(d<-dlmax){d=-dlmax;*reason=4;}
    }
    dlam=d;
  }

  for(k=0;k<pf_neq;k++) b[k]+=dlam*uf[k];
  *lam+=dlam;
  if(dlamout!=NULL) *dlamout=dlam;
  return 1;
}
