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

#include "ccxopt.h"
#include "ccxfork.h"
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

/* ---- two helpers the driver in nonlingeo() used ------------------------
   They were file-statics there, and both are pure: a dot product over the
   equation space, and the projection of an equation-space vector onto a
   nodal displacement difference.  The projection is how the method reads
   the model instead of shadowing it, which is a property of the method and
   so belongs here.                                               */

/* Project f_hat (equation space) onto a nodal displacement difference.
   Only free dofs contribute, and j starts at 1 because j=0 is the thermal
   dof, exactly as resultsini.c applies the correction.  This is how the
   path following reads the model instead of shadowing it. */

double pathfollow_dot(const double *a,const double *b,ITG n){
  ITG i;
  double s=0.;
  for(i=0;i<n;i++) s+=a[i]*b[i];
  return s;
}

double pathfollow_project(const double *fh,const double *a,const double *c,
                         const ITG *nactdof,ITG nk,ITG mt){
  ITG i,j,k;
  double p=0.;
  for(i=0;i<nk;i++){
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if(k>0) p+=fh[k-1]*(a[mt*i+j]-c[mt*i+j]);
    }
  }
  return p;
}

/* The values the sixty-three driver locals carried at their declarations. */
void pathdrv_init(pathdrv *p)
{
  memset(p,0,sizeof(*p));
  p->env=NULL;
  p->cvec=NULL; p->lhs=NULL; p->rhs=NULL; p->rhs0=NULL; p->p=NULL;
  p->q=NULL;    p->r0=NULL;  p->r1=NULL;  p->uf=NULL;   p->uref=NULL;
  p->y=NULL;    p->sv=NULL;  p->sxs=NULL; p->sxst=NULL; p->sf=NULL;
  p->sfn=NULL;  p->sstx=NULL;p->sdam=NULL;p->sxb=NULL;
  p->clip=0.05;
  p->eps=1.e-5;
  p->dtheta_eng=1.e-3;
  p->ccgrow=1.1;
}

/* ---- arming ------------------------------------------------------------

   The driver's own switches, read beside the method they drive.  It
   refuses rather than degrades, and every refusal here is a REASON rather
   than a fallback: a positive tau, nothing else already driving the load
   factor, and a domain in which the bordered system is exactly the one
   this file verifies.  A path follower that silently ran the stock control
   instead would be worse than none, because the run would look like a
   measurement of path following.

   Crack control rides on the same machinery and arms in the same block, so
   it is configured here too; the two share pathdrv and are mutually
   exclusive with each other as well.

   Everything the mesh side needs comes through trialctx.  isolver and
   ncont are passed because they are neither results() arguments nor state
   of this object - they are facts about the run this mechanism refuses
   on. */

void pathdrv_configure(pathdrv *p,const loadctl *c,const trialctx *mdl,
                       const ITG *isolver,ITG ncont)
{
  const char *e;
    double *co=*(mdl->co),*vold=*(mdl->vold);
    ITG *ipkon=*(mdl->ipkon),*kon=*(mdl->kon),*nactdof=*(mdl->nactdof);
    ITG *ne=*(mdl->ne),*nk=*(mdl->nk),*neq=*(mdl->neq),*nboun=*(mdl->nboun);
    ITG *nmethod=*(mdl->nmethod),*ithermal=*(mdl->ithermal),*mortar=*(mdl->mortar);
    ITG *iexpl=*(mdl->iexpl),*mi=*(mdl->mi);
    char *lakon=*(mdl->lakon);
    ITG num_cpus=*(mdl->num_cpus),mt=mi[1]+1,isiz;

  p->env=ccxopt_getenv("CCX_PATHFOLLOW");
  if(p->env!=NULL){
    p->tauv=atof(p->env);
    if(!(p->tauv>0.)){
      printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW must be a positive "
             "dissipation increment; got \"%s\".  Not armed.\n",p->env);
    }else if(c->diss_ctrl>=1){
      printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW and "
             "CCX_DISSIPATION_CONTROL both drive the load factor; "
             "set only one.  Not armed.\n");
    }else if((*nmethod!=1)||(*ithermal>=2)||(*mortar>1)||(ncont!=0)||
             (*iexpl>1)||(*nboun<=0)||
             ((*isolver!=0)&&(*isolver!=7))){
      printf("[PATHFOLLOW] not armed: outside the verified domain "
             "(nmethod=%" ITGFORMAT " ithermal=%" ITGFORMAT " mortar=%"
             ITGFORMAT " ncont=%" ITGFORMAT " iexpl=%" ITGFORMAT
             " nboun=%" ITGFORMAT " isolver=%" ITGFORMAT
             "; needs static, ithermal<2, no contact, implicit, "
             "prescribed dofs, SPOOLES or PARDISO)\n",
             *nmethod,*ithermal,*mortar,ncont,*iexpl,*nboun,*isolver);
    }else if(pathfollow_arm(p->tauv,neq[1])==0){
      printf("[PATHFOLLOW] *ERROR: could not allocate; not armed.\n");
    }else{
      NNEW(p->uf,double,neq[1]);
      NNEW(p->uref,double,mt**nk);
      isiz=mt**nk;cpypardou(p->uref,vold,&isiz,&num_cpus);
      p->neqarm=neq[1];
      p->taucur=p->tauv;
      p->on=1;
      if(ccxopt_getenv("CCX_PATHFOLLOW_CLIP")!=NULL)
        p->clip=atof(ccxopt_getenv("CCX_PATHFOLLOW_CLIP"));
      if(!(p->clip>0.)) p->clip=0.05;
      if(ccxopt_getenv("CCX_PATHFOLLOW_LINCHECK")!=NULL)
        p->lincheck=atoi(ccxopt_getenv("CCX_PATHFOLLOW_LINCHECK"));
      if(ccxopt_getenv("CCX_PATHFOLLOW_DTHETA")!=NULL)
        p->dtheta_eng=atof(ccxopt_getenv("CCX_PATHFOLLOW_DTHETA"));
      if(!(p->dtheta_eng>0.)) p->dtheta_eng=1.e-3;
      printf("[PATHFOLLOW] armed: tau=%.6e per increment, "
             "|dlambda| clipped at %.3e per iteration.\n",
             p->tauv,p->clip);
      printf("[PATHFOLLOW] lambda is decoupled from the step time and "
             "MAY DECREASE; theta stays monotone so dtime>0.\n");
      printf("[PATHFOLLOW] ordinary control until the measured "
             "dissipation of an accepted increment reaches 0.2*tau, "
             "then the constraint takes over.\n");

      /* ---- crack-opening control ---------------------------------
         Build the control functional c once from the REFERENCE geometry
         of the UC6 facets: phi(u)=c^T u is the mean normal separation, a
         linear functional, so dg/du=c and dg/dlambda=0 are exact by
         construction.  See the block comment in pathfollow.c for why this
         replaces the dissipation constraint on a localised cohesive
         crack. */

      /* ---- mixed-mode crack control ------------------------------
         CCX_CRACK_CONTROL=<dphi> arms the generalisation of the above:
         the control coordinate is the effective separation the UC6 law
         itself advances along, deff^2 = max(dn,0)^2+beta*|ds|^2, frozen
         into an affine functional once per attempt.  The two are
         mutually exclusive; CCX_PATHFOLLOW_COD is kept unchanged so that
         the Mode-I result stays a regression test. */

      if((ccxopt_getenv("CCX_CRACK_CONTROL")!=NULL)&&
         (ccxopt_getenv("CCX_PATHFOLLOW_COD")!=NULL)){
        printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL and "
               "CCX_PATHFOLLOW_COD both define the control coordinate; "
               "set only one.  Not armed.\n");
      }else if(ccxopt_getenv("CCX_CRACK_CONTROL")!=NULL){
        char *cce;
        ITG ce,ncoh=0;
        for(ce=0;ce<*ne;ce++){
          if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
             (lakon[8*ce+2]!='6')) continue;
          ncoh++;
        }
        if(ncoh==0){
          printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL needs UC6 "
                 "cohesive elements; none found.  Not armed.\n");
        }else{
          p->dphi=atof(ccxopt_getenv("CCX_CRACK_CONTROL"));
          if(!(p->dphi>0.)){
            printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL must be a "
                   "positive control increment.  Not armed.\n");
          }else{
            /* Default DISS: the process zone restricted to where it is
               LOADING.  Measured on the target, the unrestricted zone
               mean runs backwards while the loading mean advances
               monotonically - see crackcontrol.c. */
            p->ccmode=2;
            cce=ccxopt_getenv("CCX_CRACK_CONTROL_MODE");
            if(cce!=NULL){
              if((strcmp(cce,"MEAN")==0)||(strcmp(cce,"0")==0)) p->ccmode=0;
              else if((strcmp(cce,"ZONE")==0)||(strcmp(cce,"1")==0)) p->ccmode=1;
              else if((strcmp(cce,"DISS")==0)||(strcmp(cce,"2")==0)) p->ccmode=2;
              else printf("[CRACKCTL] unknown CCX_CRACK_CONTROL_MODE "
                          "\"%s\"; keeping DISS\n",cce);
            }
            cce=ccxopt_getenv("CCX_CRACK_CONTROL_ENGAGE");
            if(cce!=NULL) p->ccengage=atoi(cce);
            cce=ccxopt_getenv("CCX_CRACK_CONTROL_EPS");
            if(cce!=NULL) p->eps=atof(cce);
            if(!(p->eps>0.)) p->eps=1.e-5;
            cce=ccxopt_getenv("CCX_CRACK_CONTROL_GROW");
            if(cce!=NULL) p->ccgrow=atof(cce);
            if(!(p->ccgrow>=1.)) p->ccgrow=1.1;
            NNEW(p->cvec,double,neq[1]);
            if(pathfollow_cod_arm(p->cvec,neq[1])==1){
              p->codmode=2;
              p->ccarmed=1;
              p->dphicur=p->dphi;
              p->engaged=0;      /* ordinary control until the crack has
                                    a process zone to control          */
              printf("[CRACKCTL] armed on %" ITGFORMAT " UC6 facet(s): "
                     "mode=%s, control increment dphi=%.6e per "
                     "increment.\n",ncoh,
                     (p->ccmode==0)?"MEAN (all facets)":
                     ((p->ccmode==1)?"ZONE (process zone)":
                      "DISS (cohesive dissipation)"),p->dphi);
              printf("[CRACKCTL] the functional is refrozen from the "
                     "committed state at every attempt and the "
                     "constraint is incremental: g = c^T(u-u_n) - "
                     "dphi.\n");
              if(p->ccengage>0)
                printf("[CRACKCTL] engagement deferred to increment %"
                       ITGFORMAT ".\n",p->ccengage);
              else
                printf("[CRACKCTL] engages as soon as the process zone "
                       "is non-empty.\n");
            }else{
              printf("[CRACKCTL] *ERROR: could not arm.\n");
            }
          }
        }
      }

      if(ccxopt_getenv("CCX_PATHFOLLOW_COD")!=NULL){
        ITG ce,ci,ck,cip,cn,cnp,cdof,ncoh=0;
        double ca[3],cb[3],cnv[3],cnorm,csh[3],cw;
        for(ce=0;ce<*ne;ce++){
          if(ipkon[ce]<0) continue;
          if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
             (lakon[8*ce+2]!='6')) continue;
          ncoh++;
        }
        if(ncoh==0){
          printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW_COD needs UC6 "
                 "cohesive elements; none found.  Not armed.\n");
        }else{
          NNEW(p->cvec,double,neq[1]);
          cw=1./(3.*(double)ncoh);
          for(ce=0;ce<*ne;ce++){
            if(ipkon[ce]<0) continue;
            if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
               (lakon[8*ce+2]!='6')) continue;
            for(ck=0;ck<3;ck++){
              ca[ck]=co[3*(kon[ipkon[ce]+1]-1)+ck]
                    -co[3*(kon[ipkon[ce]+0]-1)+ck];
              cb[ck]=co[3*(kon[ipkon[ce]+2]-1)+ck]
                    -co[3*(kon[ipkon[ce]+0]-1)+ck];
            }
            cnv[0]=ca[1]*cb[2]-ca[2]*cb[1];
            cnv[1]=ca[2]*cb[0]-ca[0]*cb[2];
            cnv[2]=ca[0]*cb[1]-ca[1]*cb[0];
            cnorm=sqrt(cnv[0]*cnv[0]+cnv[1]*cnv[1]+cnv[2]*cnv[2]);
            if(!(cnorm>1.e-30)) continue;
            for(ck=0;ck<3;ck++) cnv[ck]/=cnorm;
            for(cip=0;cip<3;cip++){
              for(ci=0;ci<3;ci++) csh[ci]=1./6.;
              csh[cip]=2./3.;
              for(ci=0;ci<3;ci++){
                cn =kon[ipkon[ce]+ci];      /* minus side */
                cnp=kon[ipkon[ce]+ci+3];    /* plus side  */
                for(ck=1;ck<mt;ck++){
                  cdof=nactdof[mt*(cnp-1)+ck];
                  if(cdof>0) p->cvec[cdof-1]+=cw*csh[ci]*cnv[ck-1];
                  cdof=nactdof[mt*(cn-1)+ck];
                  if(cdof>0) p->cvec[cdof-1]-=cw*csh[ci]*cnv[ck-1];
                }
              }
            }
          }
          p->dphi=atof(ccxopt_getenv("CCX_PATHFOLLOW_COD"));
          if(!(p->dphi>0.)) p->dphi=5.e-5;
          if(pathfollow_cod_arm(p->cvec,neq[1])==1){
            p->codmode=1;
            p->engaged=1;     /* no warm-up: phi is meaningful at once */
            printf("[PATHFOLLOW] crack-opening control armed on %"
                   ITGFORMAT " UC6 facet(s); mean normal separation "
                   "advances by %.6e per increment.\n",ncoh,p->dphi);
            {ITG cnz=0; double cnn=0.;
             for(ce=0;ce<neq[1];ce++){
               cnn+=p->cvec[ce]*p->cvec[ce];
               if(p->cvec[ce]!=0.) cnz++;}
             printf("[PATHFOLLOW] control functional: |c|=%.6e, %"
                    ITGFORMAT " non-zero of %" ITGFORMAT " equations\n",
                    sqrt(cnn),cnz,neq[1]);}
            printf("[PATHFOLLOW] the control is monotone through a "
                   "snap-back by construction; lambda is free.\n");
          }
        }
      }
    }
  }
}

/* ---- the predictor -----------------------------------------------------

   The largest block that was inline: 343 lines that capture f_hat, solve
   for the reference direction, and place the load factor for this attempt
   before the corrector runs.

   It calls the LINEAR SOLVER, which is why ad, au, icol, isolver, sigma
   and the three format flags are parameters: they are the factorisation's
   arguments, not this object's state and not results() arguments either.
   Everything else arrives through trialctx.

   The guard stays with the caller, as everywhere else here.            */

void pathdrv_predictor(pathdrv *p,glob_census *g,const trialctx *mdl,
                       const double *xboun,const double *xbounold,
                       double *ad,double *au,ITG *icol,const ITG *isolver,
                       double sigma,ITG inputformat,ITG nrhs,
                       ITG symmetryflag,ITG iit)
{
  double *b=*(mdl->b),*vold=*(mdl->vold),*vini=*(mdl->vini);
  double *xbounact=*(mdl->xbounact),*adb=*(mdl->adb),*aub=*(mdl->aub);
  double *f=*(mdl->f),*fext=*(mdl->fext);
  ITG *neq=*(mdl->neq),*nzs=*(mdl->nzs),*irow=*(mdl->irow),*jq=*(mdl->jq);
  ITG *nactdof=*(mdl->nactdof),*mi=*(mdl->mi),*nk=*(mdl->nk),*nboun=*(mdl->nboun);
  double *dam=*(mdl->dam),*xstate=*(mdl->xstate),*fn=*(mdl->fn),*stx=*(mdl->stx);
  double *xstiff=*(mdl->xstiff),*qa=mdl->qa,*cam=mdl->cam;
  ITG *nstate_=*(mdl->nstate_),*ne=*(mdl->ne);
  ITG num_cpus=*(mdl->num_cpus),iinc=*(mdl->iinc),nasym=*(mdl->nasym);
  ITG mt=mi[1]+1,isiz,k;

	  const double *fhat=pathfollow_fhat();
	  for(k=0;k<neq[1];k++) p->uf[k]=fhat[k];
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles(ad,au,adb,aub,&sigma,p->uf,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,&nzs[2]);
#endif
	  }else if(*isolver==7){
#ifdef PARDISO
	    pardiso_main(ad,au,adb,aub,&sigma,p->uf,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#endif
	  }
	  p->nstep++;

	  /* ff is recorded on every iterate, engaged or not, so that the
	     tangent predictor has a stiffness the moment the constraint
	     takes over.  Only the APPLICATION is gated on engagement. */

	  if((p->rhs0!=NULL)&&(ccxopt_getenv("CCX_PATHFOLLOW_SOLVECHECK")!=NULL)){
	    double pfr1=0.,pfr2=0.,pfn1=0.,pfn2=0.,pft;
	    ITG pfi,pfone=1;
	    const double *pffh2=pathfollow_fhat();
	    for(pfi=0;pfi<neq[1];pfi++) p->y[pfi]=0.;
	    FORTRAN(op,(b,p->y,ad,au,jq,irow,&pfone,&neq[0]));
	    for(pfi=0;pfi<neq[0];pfi++){
	      pft=p->y[pfi]-p->rhs0[pfi];
	      pfr1+=pft*pft;pfn1+=p->rhs0[pfi]*p->rhs0[pfi];
	    }
	    for(pfi=0;pfi<neq[1];pfi++) p->y[pfi]=0.;
	    FORTRAN(op,(p->uf,p->y,ad,au,jq,irow,&pfone,&neq[0]));
	    for(pfi=0;pfi<neq[0];pfi++){
	      pft=p->y[pfi]-pffh2[pfi];
	      pfr2+=pft*pft;pfn2+=pffh2[pfi]*pffh2[pfi];
	    }
	    printf("[PF-SOLVE] it=%" ITGFORMAT " |K*duR-rhs|/|rhs|=%.3e "
	           "|K*duF-fhat|/|fhat|=%.3e\n",iit,
	           sqrt(pfr1)/(sqrt(pfn1)+1.e-300),
	           sqrt(pfr2)/(sqrt(pfn2)+1.e-300));
	    fflush(stdout);
	  }

/* Evaluate CalculiX's own residual at (vold, LAMV) and copy fext-f into
   DST.  b is zeroed first so results() leaves the displacement where it
   is: the probe must read the residual AT a state, never move it.

   iout=-1, not 0.  resultsini gates the correction on iout>-1 and the
   prescribed-dof update on abs(iout)<2, so -1 sets the prescribed dofs
   from xbounact (which IS "evaluating R at this lambda") while leaving the
   free dofs alone - and, crucially, it does not request the output fields.
   iout=0 does request them, and results() then writes inum/een/emn/epn,
   which are not allocated on this path: measured as a SIGSEGV inside
   memset called from nonlingeo. */

#define PF_LIN_RESID(LAMV,DST) do{                                    \
  ITG _k,_io=*(mdl->iout);                                              \
  /* Every evaluation starts from the identical base state.  Without  \
     this the probe measures its own leakage: idempotency (evaluating \
     R at lambda0 twice) came out at 1.7e+05 relative near the limit  \
     point, which is exactly the spurious "discontinuity" the sweep   \
     then reported.  vold is NOT restored here - the caller perturbs it\
     on purpose. */                                                   \
  isiz=*nstate_*mi[0]**ne;cpypardou(xstate,p->sxs,&isiz,&num_cpus);   \
  isiz=27*mi[0]**ne;cpypardou(xstiff,p->sxst,&isiz,&num_cpus);        \
  isiz=neq[1];cpypardou(f,p->sf,&isiz,&num_cpus);                     \
  isiz=mt**nk;cpypardou(fn,p->sfn,&isiz,&num_cpus);                   \
  isiz=6*mi[0]**ne;cpypardou(stx,p->sstx,&isiz,&num_cpus);            \
  if((dam!=NULL)&&(p->sdam!=NULL)){                                   \
    isiz=mi[0]**ne;cpypardou(dam,p->sdam,&isiz,&num_cpus);}           \
  for(_k=0;_k<4;_k++) qa[_k]=p->sqa[_k];                              \
  for(_k=0;_k<5;_k++) cam[_k]=p->scam[_k];                            \
  for(_k=0;_k<*nboun;_k++)                                            \
    xbounact[_k]=xbounold[_k]+(xboun[_k]-xbounold[_k])*(LAMV);        \
  for(_k=0;_k<neq[1];_k++) b[_k]=0.;                                  \
  *(mdl->iout)=-1;                                                      \
  isiz=mt**nk;cpypardou(*(mdl->v),vold,&isiz,&num_cpus);                \
  trial_results(mdl);                                                   \
  trial_reduce(mdl,b);                                                  \
  isiz=neq[1];cpypardou((DST),b,&isiz,&num_cpus);                     \
  *(mdl->iout)=_io;                                                     \
    }while(0)

	  /* ================= CCX_PATHFOLLOW_LINCHECK =====================
	     The decisive gate.  A synthetic self test proves the formulas of
	     pathfollow.c; it proves nothing about how u, lambda, xbounact, b
	     and K are actually wired together inside CalculiX.  This measures
	     that wiring, at the REAL base point - after prediction() has
	     extrapolated and after the operator has been assembled - by
	     finite differences of the residual CalculiX itself computes.

	     Checked, from one base point, restoring the full trial state
	     between evaluations, for a sweep of decreasing eps:

	       q_FD  = [R(u,lambda+eps) - R(u,lambda)]/eps
	       row 1 : [R(u+eps*p, lambda+eps*dl) - R(u,lambda)]/eps
	               against  K*p + q_FD*dl

	     and the frozen reference vector f_hat against -q_FD, which is
	     the claim "a constant reference load reproduces the actual
	     dR/dlambda" that must not be assumed. */

	  if((p->lincheck>0)&&(iinc==p->lincheck)&&(iit==1)){
	    ITG pq,pfone=1,pj,pnode,pi;
	    double peps,plam0,pdl,pn1,pn2,pnq,pt,pfhq,pnf,pkp,pql,pqc;
	    double pepsv[5]={1.e-3,1.e-4,1.e-5,1.e-6,1.e-7};

	    if(p->p==NULL){
	      NNEW(p->p,double,neq[1]);NNEW(p->r0,double,neq[1]);
	      NNEW(p->r1,double,neq[1]);NNEW(p->q,double,neq[1]);
	      NNEW(p->sv,double,mt**nk);
	      NNEW(p->sxs,double,*nstate_*mi[0]**ne);
	      if(dam!=NULL) NNEW(p->sdam,double,mi[0]**ne);
	      NNEW(p->sf,double,neq[1]);NNEW(p->sfn,double,mt**nk);
	      NNEW(p->sstx,double,6*mi[0]**ne);
	      NNEW(p->sxb,double,*nboun);
	      /* p->y is the mat-vec scratch; it is otherwise allocated only
	         by SOLVECHECK, and LINCHECK dereferenced it as NULL. */
	      if(p->y==NULL) NNEW(p->y,double,neq[1]);
	      NNEW(p->sxst,double,27*mi[0]**ne);
	      NNEW(p->lhs,double,neq[1]);NNEW(p->rhs,double,neq[1]);
	    }

	    /* --- snapshot everything the probe is about to disturb --- */
	    isiz=mt**nk;cpypardou(p->sv,vold,&isiz,&num_cpus);
	    isiz=*nstate_*mi[0]**ne;cpypardou(p->sxs,xstate,&isiz,&num_cpus);
	    if((dam!=NULL)&&(p->sdam!=NULL)){
	      isiz=mi[0]**ne;cpypardou(p->sdam,dam,&isiz,&num_cpus);}
	    isiz=neq[1];cpypardou(p->sf,f,&isiz,&num_cpus);
	    isiz=mt**nk;cpypardou(p->sfn,fn,&isiz,&num_cpus);
	    isiz=6*mi[0]**ne;cpypardou(p->sstx,stx,&isiz,&num_cpus);
	    isiz=*nboun;cpypardou(p->sxb,xbounact,&isiz,&num_cpus);
	    isiz=27*mi[0]**ne;cpypardou(p->sxst,xstiff,&isiz,&num_cpus);
	    for(k=0;k<4;k++) p->sqa[k]=qa[k];
	    for(k=0;k<5;k++) p->scam[k]=cam[k];
	    isiz=neq[1];cpypardou(p->p,b,&isiz,&num_cpus);   /* direction p */
	    plam0=p->lam;
	    pdl=(fabs(p->dlampred)>0.)?p->dlampred:1.e-3;

	    /* p is left as du_R, the Newton correction itself, and is NOT
	       normalised.  Rescaling it to unit norm makes eps*p a
	       perturbation of order 1e-3 in displacement, which is larger
	       than the cohesive d0 = 4e-4 and therefore no longer probes the
	       branch the solver is on; the sweep then stopped converging for
	       a reason that has nothing to do with the tangent.  du_R is the
	       direction the method actually takes, which is the one worth
	       verifying. */
	    pt=sqrt(pathfollow_dot(p->p,p->p,neq[0]));

	    printf("\n[PF-LIN] increment %" ITGFORMAT ", base point AFTER "
	           "prediction(); lambda=%.10f  |p|=%.4e  dlambda_dir=%.4e\n",
	           iinc,plam0,sqrt(pathfollow_dot(p->p,p->p,neq[1])),pdl);

	    /* R at the base point */
	    PF_LIN_RESID(plam0,p->r0);
	    pn1=sqrt(pathfollow_dot(p->r0,p->r0,neq[0]));
	    printf("[PF-LIN]   |R0| = %.6e   neq0=%" ITGFORMAT " neq1=%"
	           ITGFORMAT "\n",pn1,neq[0],neq[1]);

	    /* bisection: is ONE evaluation self-contaminating, or is it the
	       perturbed ones that leave something behind? */
	    PF_LIN_RESID(plam0,p->r1);
	    pn2=0.;
	    for(k=0;k<neq[0];k++){pt=p->r1[k]-p->r0[k];pn2+=pt*pt;}
	    printf("[PF-LIN]   immediate re-evaluation at lambda0: "
	           "|dR|/|R0| = %.6e\n",(pn1>0.)?sqrt(pn2)/pn1:sqrt(pn2));
	    fflush(stdout);

	    printf("[PF-LIN]   eps        |q_FD|      "
	           "cos(f_hat,-q_FD)  |f_hat|/|q_FD|   row1 rel.err\n");
	    for(pq=0;pq<10;pq++){
	      if(pq==5){
	        /* Second sweep with dlambda = 0.  This removes the load term
	           entirely and tests ONLY whether the assembled operator is
	           the derivative of the residual CalculiX computes, which is
	           a statement about the element tangent and has nothing to do
	           with path following. */
	        pdl=0.;
	        printf("[PF-LIN]   ---- dlambda = 0: pure tangent test, "
	               "[R(u+eps*p)-R(u)]/eps  against  K*p ----\n");
	      }
	      peps=pepsv[pq%5];

	      /* Re-anchor the base residual for every eps.  A single
	         evaluation is idempotent (measured: |dR|/|R0| = 0 exactly),
	         but the perturbed ones must not be allowed to bias the next
	         difference. */
	      PF_LIN_RESID(plam0,p->r0);

	      /* ORDER MATTERS.  p->q was originally computed first and read
	         last, and it did not survive the calls in between: |q| taken
	         at the norm loop disagreed with the |q_FD| printed from the
	         same buffer moments earlier (5.3e+03 against 8.1e-11).  The
	         perturbed evaluation and the mat-vec are therefore done
	         FIRST, and q is rebuilt immediately before it is used, so
	         nothing can run between its definition and its use. */

	      /* (a) full perturbation, in u and lambda together */
	      for(pi=0;pi<*nk;pi++){
	        for(pj=1;pj<mt;pj++){
	          pnode=nactdof[mt*pi+pj];
	          if(pnode>0) vold[mt*pi+pj]+=peps*p->p[pnode-1];
	        }
	      }
	      PF_LIN_RESID(plam0+peps*pdl,p->r1);
	      isiz=mt**nk;cpypardou(vold,p->sv,&isiz,&num_cpus);
	      for(k=0;k<neq[0];k++) p->lhs[k]=(p->r1[k]-p->r0[k])/peps;

	      /* (b) K*p */
	      for(k=0;k<neq[1];k++) p->y[k]=0.;
	      FORTRAN(op,(p->p,p->y,ad,au,jq,irow,&pfone,&neq[0]));

	      /* (c) q_FD, built last */
	      PF_LIN_RESID(plam0+peps,p->r1);
	      for(k=0;k<neq[0];k++) p->q[k]=(p->r1[k]-p->r0[k])/peps;
	      pnq=sqrt(pathfollow_dot(p->q,p->q,neq[0]));

	      pfhq=0.;pnf=0.;
	      if(pathfollow_have()==1){
	        const double *pfh=pathfollow_fhat();
	        pfhq=pathfollow_dot(pfh,p->q,neq[0]);
	        pnf=sqrt(pathfollow_dot(pfh,pfh,neq[0]));
	      }

	      pn2=0.;pn1=0.;pkp=0.;pql=0.;pqc=0.;
	      for(k=0;k<neq[0];k++) p->rhs[k]=-p->y[k]+pdl*p->q[k];
	      for(k=0;k<neq[0];k++){
	        pt=p->lhs[k]-p->rhs[k];  pn2+=pt*pt;
	        pn1+=p->lhs[k]*p->lhs[k];
	        pkp+=p->y[k]*p->y[k];
	        pql+=p->rhs[k]*p->rhs[k];
	        pqc+=p->q[k]*p->q[k];
	      }
	      printf("[PF-LIN]   %.1e  %.6e  %+.10f   %12.6e   %.6e"
	             "   |lhs|=%.4e |Kp|=%.4e |rhs|=%.4e  |q|inloop=%.6e\n",
	             peps,pnq,(pnf*pnq>0.)?pfhq/(pnf*pnq):0.,
	             (pnq>0.)?pnf/pnq:0.,
	             (pn1>0.)?sqrt(pn2)/sqrt(pn1):sqrt(pn2),
	             sqrt(pn1),sqrt(pkp),sqrt(pql),sqrt(pqc));
	      fflush(stdout);
	    }

	    /* Idempotency of the probe itself.  Everything above is only
	       meaningful if evaluating R at the SAME point twice returns the
	       same vector; if it does not, the sweep is measuring the probe's
	       own leakage rather than dR/dlambda. */

	    PF_LIN_RESID(plam0,p->r1);
	    pn1=0.;pn2=0.;
	    for(k=0;k<neq[0];k++){
	      pt=p->r1[k]-p->r0[k];pn2+=pt*pt;pn1+=p->r0[k]*p->r0[k];
	    }
	    printf("[PF-LIN]   idempotency: |R(lam0) again - R0|/|R0| = %.6e"
	           "   %s\n",(pn1>0.)?sqrt(pn2)/sqrt(pn1):sqrt(pn2),
	           (sqrt(pn2)<=1.e-10*sqrt(pn1))?"clean":"PROBE LEAKS STATE");

	    /* --- restore --- */
	    isiz=mt**nk;cpypardou(vold,p->sv,&isiz,&num_cpus);
	    isiz=*nstate_*mi[0]**ne;cpypardou(xstate,p->sxs,&isiz,&num_cpus);
	    if((dam!=NULL)&&(p->sdam!=NULL)){
	      isiz=mi[0]**ne;cpypardou(dam,p->sdam,&isiz,&num_cpus);}
	    isiz=neq[1];cpypardou(f,p->sf,&isiz,&num_cpus);
	    isiz=mt**nk;cpypardou(fn,p->sfn,&isiz,&num_cpus);
	    isiz=6*mi[0]**ne;cpypardou(stx,p->sstx,&isiz,&num_cpus);
	    isiz=*nboun;cpypardou(xbounact,p->sxb,&isiz,&num_cpus);
	    isiz=27*mi[0]**ne;cpypardou(xstiff,p->sxst,&isiz,&num_cpus);
	    for(k=0;k<4;k++) qa[k]=p->sqa[k];
	    for(k=0;k<5;k++) cam[k]=p->scam[k];
	    isiz=neq[1];cpypardou(b,p->p,&isiz,&num_cpus);
	    p->lam=plam0;
	    printf("[PF-LIN] state restored\n\n");
	    fflush(stdout);
	  }
#undef PF_LIN_RESID

	  pathfollow_measure(p->uf);
	  p->applied=0;
	  if(p->engaged==1){

	    /* the trial displacement is READ FROM THE MODEL, not
	       accumulated; see the design note in pathfollow.c */

	    if(p->codmode>=1){

	      /* Mode 1 (CCX_PATHFOLLOW_COD): phi is ABSOLUTE, the mean normal
	         separation measured from the reference configuration.
	         Mode 2 (CCX_CRACK_CONTROL): phi is INCREMENTAL, c^T(u-u_n),
	         because the functional itself is refrozen every attempt and
	         only its increment is meaningful across a refreeze.

	         c only couples plus/minus node pairs, so in both cases the
	         projection is exactly the controlled opening. */

	      p->lam0it=p->lam;
	      p->dlamit=0.;
	      p->cu=pathfollow_project(pathfollow_cod_c(),vold,
	                       (p->codmode==2)?vini:p->uref,nactdof,*nk,mt);
	      p->applied=pathfollow_cod_step(b,p->uf,p->cu,&p->lam,p->clip,
	                                     &p->g,&p->dlam,&p->reason);
	      p->dg=p->cu;
	      if(p->applied==1) glob_fired(g,GLOB_PATHFOLLOW);
	      if(p->applied==1){
	        p->dlamit=p->lam-p->lam0it;
	        for(k=0;k<*nboun;k++){
	          xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*p->lam;
	        }
	      }
	      goto pf_after_step;
	    }
	    p->pdu=pathfollow_project(pathfollow_fhat(),vold,vini,nactdof,*nk,mt);

	    /* At the first iteration put lambda exactly on the constraint,
	       in closed form, so that the coupled iteration does not open
	       with a residual the stock extrapolation put there. */

	    /* Measured: this does NOT rescue the descending branch (1 accepted
	       increment past engagement against 2 without it), so it is
	       opt-in and off by default.  It is kept because it isolates one
	       hypothesis cleanly - the opening constraint residual is not
	       what stops the run. */

	    if((iit==1)&&(ccxopt_getenv("CCX_PATHFOLLOW_PROJECT")!=NULL)){
	      if(pathfollow_project_lambda(p->pdu,&p->lam)==1){
	        for(k=0;k<*nboun;k++){
	          xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*p->lam;
	        }
	      }
	    }
	    p->applied=pathfollow_step(b,p->uf,p->pdu,&p->lam,p->clip,
	                               &p->dg,&p->g,&p->dlam,&p->reason);
	  }
	  if(p->applied==1){
	    for(k=0;k<*nboun;k++){
	      xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*p->lam;
	    }
	  }
	 pf_after_step:
	  if(ccxopt_getenv("CCX_PATHFOLLOW_PROBE")!=NULL){
	    printf("[PATHFOLLOW] it=%" ITGFORMAT " lambda=%.8f dlambda=%+.4e "
	           "dG=%.6e g=%.4e applied=%" ITGFORMAT " reason=%" ITGFORMAT
	           "\n",iit,p->lam,p->applied?p->dlam:0.,p->dg,p->g,
	           p->applied,p->reason);
	    fflush(stdout);
	  }
}
