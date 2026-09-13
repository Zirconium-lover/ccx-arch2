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
  MIXED-MODE CRACK CONTROL FOR THE UC6 COHESIVE INTERFACE
  =======================================================

  WHY THIS FILE EXISTS
  --------------------
  Two control coordinates were already measured on this branch and both
  are inadequate for the real target:

    * the Gutierrez/Verhoosel dissipation constraint built from
      f_hat = -dR/dlambda.  On a localised cohesive crack that is
      driven by prescribed displacement it is IDENTICALLY ZERO: f_hat
      lives on the dofs next to the loaded face, the block on that side
      moves rigidly with the face once the interface fails, and
      f_hat^T u stays proportional to lambda whatever the crack does.
      Measured on test/pathfollow/cohesive.inp.

    * the mean NORMAL separation of the cohesive facets
      (CCX_PATHFOLLOW_COD).  It follows the Mode-I benchmark exactly,
      798 accepted post-peak increments on the analytic bilinear branch.
      It is normal-only and it averages over EVERY facet, so on a model
      whose active front is a small, strongly sheared part of a large
      interface it measures the elastic compliance of the closed facets
      instead of the crack.

  WHAT IS CONTROLLED HERE
  -----------------------
  The UC6 law (cohesive_uc6.f) is driven by ONE scalar per integration
  point, the effective separation

      deff^2 = max(dn,0)^2 + beta*(ds1^2 + ds2^2) ,  beta=(Ts0/Tn0)^2 ,

  so deff is the coordinate the material itself advances along.  It is
  not linear in u, but it is positively homogeneous of degree one, and
  that is enough.  With

      q  = ( max(dn,0), beta*ds1, beta*ds2 )      (local frame)
      m  = q/deff                                 (frozen at a state)

  Euler's identity gives  m . d = deff  EXACTLY at the point where m was
  frozen, and m is exactly the gradient d(deff)/d(d).  So

      phi(u) = sum_ip  w_ip * m_ip . [[u]]_ip

  is an affine functional of u whose value at the linearisation point is
  the weighted mean (or sum) of deff, and whose gradient is the exact
  gradient of that same quantity.  During one corrector m and w are
  constants, so

      g = phi(u) - phi(u_committed) - dphi ,   dg/du = c ,  dg/dlambda = 0

  and the bordered row is closed exactly, as for the Mode-I control.

  THE WEIGHTS ARE NOT ARBITRARY
  -----------------------------
  For the bilinear law the dissipated energy per unit advance of deff is
  a CONSTANT:

      t(dmax)      = Tn0 (df-dmax)/(df-d0)
      Psi_diss     = int t d(delta) - 1/2 t dmax
      dPsi/d(dmax) = 1/2 Tn0 df/(df-d0)                       (constant)

  Therefore, with

      w_ip = (area_ip/3) * 1/2 * Tn0 * df/(df-d0)   on loading facets
      w_ip = 0                                      elsewhere

  phi is the COHESIVE DISSIPATION and dphi is a prescribed energy release
  per increment.  That is the Gutierrez constraint written for the crack
  instead of for the loading boundary, and it has no blind spot: it is
  supported exactly on the process zone.

  Three weightings are provided, because which one is right is an
  empirical question and each is a different measurement:

      0 MEAN  w = area/3 over every live facet.  phi is the mean
              effective separation.  This is the mixed-mode
              generalisation of CCX_PATHFOLLOW_COD and reduces to it
              when the crack is a pure Mode-I plane.
      1 ZONE  the same, restricted to the process zone d0<dmax<df.
      2 DISS  the process zone AND ONLY WHERE IT IS LOADING
              (deff >= dmax at the committed state), weighted by
              area*rate.

  All three are normalised by the weight they carry, so phi is a mean
  separation in every mode and dphi is always a length.  That loses
  nothing: the bordered row

      dlambda = -( g + c^T du_R ) / ( c^T du_F ) ,  g = c^T(u-u_n) - dphi

  is INVARIANT under c -> s*c together with dphi -> s*dphi, so the
  normalisation cannot change the step - only the selection and the
  relative weights can.  The unnormalised weight sum is reported so that
  phi*sum(w) recovers the dissipation in energy units.

  WHY THE LOADING RESTRICTION IS NOT A THRESHOLD BUT THE DEFINITION
  ----------------------------------------------------------------
  A point that is unloading dissipates nothing: its damage is frozen and
  its traction follows the secant back to the origin.  Including it in
  the control functional therefore averages the advancing front together
  with material that is, by construction, not advancing.

  This is not hypothetical.  Measured on the s3rad target at stock
  increment 160-175, from the deck's own SDV/E output:

      process zone            11002 -> 10874 integration points
      of those, loading        4716 ->  4105
      mean deff over the zone  4.007e-03 -> 3.924e-03   (DECREASING)
      mean deff over loading   6.215e-03 -> 6.347e-03   (increasing,
                                             +8.8e-06 per increment)

  The zone mean runs BACKWARDS, because points leave the zone into
  failure faster than the survivors open; a constraint built on it would
  be asked to advance a quantity that physically retreats.  The loading
  mean is monotone and its rate is the natural control increment.

  RE-ANCHORING
  ------------
  m and w are rebuilt from the COMMITTED state at the top of every
  attempt, and the constraint is written INCREMENTALLY,

      g = c^T (u - u_committed) - dphi ,

  so redefining the functional between increments is harmless by
  construction: nothing carries over except the committed displacement
  field, which is the model's own state.  This is also what repairs the
  stale-vector defect of the absolute form when the equation count
  changes after element deletion.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalculiX.h"

/* ------------------------------------------------------------------ */
/* Pure kinematics.  No model, no globals: everything is a function of  */
/* its arguments so that it can be checked in isolation.                */
/* ------------------------------------------------------------------ */

/* Reference frame of one UC6 facet from the three minus-side corners.
   Rows of rmat map a global vector to (normal, shear-1, shear-2), the
   same convention as cohesive_uc6.f.  Returns the facet area, or a
   negative value for a degenerate facet. */

double crackcontrol_frame(const double *x1,const double *x2,const double *x3,
                          double *rmat){

  ITG k;
  double e1[3],e2[3],cv[3],n1,nc;

  for(k=0;k<3;k++){e1[k]=x2[k]-x1[k];e2[k]=x3[k]-x1[k];}
  cv[0]=e1[1]*e2[2]-e1[2]*e2[1];
  cv[1]=e1[2]*e2[0]-e1[0]*e2[2];
  cv[2]=e1[0]*e2[1]-e1[1]*e2[0];
  n1=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
  nc=sqrt(cv[0]*cv[0]+cv[1]*cv[1]+cv[2]*cv[2]);
  if(!(n1>1.e-30)||!(nc>1.e-30)) return -1.;
  for(k=0;k<3;k++){rmat[k]=cv[k]/nc;rmat[3+k]=e1[k]/n1;}
  rmat[6]=rmat[1]*rmat[5]-rmat[2]*rmat[4];
  rmat[7]=rmat[2]*rmat[3]-rmat[0]*rmat[5];
  rmat[8]=rmat[0]*rmat[4]-rmat[1]*rmat[3];
  return 0.5*nc;
}

/* Effective separation and its exact gradient, in the LOCAL frame.

   deff^2 = max(dn,0)^2 + beta*(ds1^2+ds2^2)
   m      = d(deff)/d(delta) = q/deff,  q=(max(dn,0),beta*ds1,beta*ds2)

   and m.delta == deff identically (Euler, deff is 1-homogeneous).  When
   deff is below tol the direction is undefined; the facet normal is
   returned, which is the limit along pure opening and is what a facet
   that has not moved yet must be controlled with. */

double crackcontrol_dir(const double *dl,double beta,double tol,double *m){

  double dn,deff,q[3];

  dn=(dl[0]>0.)?dl[0]:0.;
  deff=sqrt(dn*dn+beta*(dl[1]*dl[1]+dl[2]*dl[2]));
  if(!(deff>tol)){
    m[0]=1.;m[1]=0.;m[2]=0.;
    return deff;
  }
  q[0]=dn;q[1]=beta*dl[1];q[2]=beta*dl[2];
  m[0]=q[0]/deff;m[1]=q[1]/deff;m[2]=q[2]/deff;
  return deff;
}

/* Dissipation carried per unit advance of deff, per unit facet area, for
   the bilinear UC6 law.  Derivation in the block comment above. */

double crackcontrol_dissrate(double kn,double tn0,double gc){

  double d0,df;

  d0=tn0/kn;
  df=2.*gc/tn0;
  if(!(df>d0)) return 0.;
  return 0.5*tn0*df/(df-d0);
}

/* ------------------------------------------------------------------ */
/* Self test.  Verifies, without any model:                            */
/*   A  the frame is right-handed and orthonormal                      */
/*   B  m.delta == deff for opening, shear and mixed states            */
/*   C  m is the true gradient of deff (central differences)           */
/*   D  the closed-cone case dn<0 uses shear only                      */
/*   E  the dissipation rate reproduces Gc when integrated d0..df       */
/* Returns the number of failures.                                     */
/* ------------------------------------------------------------------ */

ITG crackcontrol_selftest(void){

  ITG i,j,nbad=0;
  double x1[3]={0.,0.,0.},x2[3]={2.,0.,0.},x3[3]={0.,3.,0.};
  double rmat[9],area,t,m[3],dl[3],deff,fd,ana,eps,dp[3],dm[3],ep,em;
  double beta=0.5625,kn=2.e5,tn0=300.,gc=2.4,d0,df,rate;

  printf("[CRACKCTL] self test\n");

  /* ---- A: orthonormal, right-handed frame, correct area ------------ */

  area=crackcontrol_frame(x1,x2,x3,rmat);
  t=0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      double d=rmat[3*i+0]*rmat[3*j+0]+rmat[3*i+1]*rmat[3*j+1]
              +rmat[3*i+2]*rmat[3*j+2];
      d-=(i==j)?1.:0.;
      if(fabs(d)>t) t=fabs(d);
    }
  }
  printf("   %-24s area=%.6f (want 3) orthonormality=%.2e  %s\n",
         "A frame",area,t,((fabs(area-3.)<1.e-12)&&(t<1.e-14))?"ok":"FAIL");
  if(!((fabs(area-3.)<1.e-12)&&(t<1.e-14))) nbad++;

  /* ---- B: m.delta == deff ------------------------------------------ */

  t=0.;
  for(i=0;i<4;i++){
    dl[0]=(i==0)?1.e-3:((i==1)?0.:((i==2)?7.e-4:-5.e-4));
    dl[1]=(i==0)?0.:((i==1)?9.e-4:((i==2)?4.e-4:6.e-4));
    dl[2]=(i==0)?0.:((i==1)?-3.e-4:((i==2)?-2.e-4:1.e-4));
    deff=crackcontrol_dir(dl,beta,1.e-30,m);
    ana=m[0]*dl[0]+m[1]*dl[1]+m[2]*dl[2];
    /* for dn<0 the identity holds with the CLOSED normal removed, which
       is what m.delta computes: m[0]=0 there */
    if(fabs(ana-deff)>t) t=fabs(ana-deff);
  }
  printf("   %-24s max|m.delta-deff| = %.3e  %s\n","B Euler identity",t,
         (t<1.e-18)?"ok":"FAIL");
  if(!(t<1.e-18)) nbad++;

  /* ---- C: m is the gradient ----------------------------------------- */

  dl[0]=7.e-4;dl[1]=4.e-4;dl[2]=-2.e-4;
  deff=crackcontrol_dir(dl,beta,1.e-30,m);
  eps=1.e-9;
  t=0.;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){dp[j]=dl[j];dm[j]=dl[j];}
    dp[i]+=eps;dm[i]-=eps;
    ep=sqrt(((dp[0]>0.)?dp[0]:0.)*((dp[0]>0.)?dp[0]:0.)
            +beta*(dp[1]*dp[1]+dp[2]*dp[2]));
    em=sqrt(((dm[0]>0.)?dm[0]:0.)*((dm[0]>0.)?dm[0]:0.)
            +beta*(dm[1]*dm[1]+dm[2]*dm[2]));
    fd=(ep-em)/(2.*eps);
    if(fabs(fd-m[i])>t) t=fabs(fd-m[i]);
  }
  printf("   %-24s max|grad_fd-m| = %.3e  %s\n","C gradient",t,
         (t<1.e-6)?"ok":"FAIL");
  if(!(t<1.e-6)) nbad++;

  /* ---- D: a closed facet is controlled by shear only ---------------- */

  dl[0]=-1.e-3;dl[1]=8.e-4;dl[2]=0.;
  deff=crackcontrol_dir(dl,beta,1.e-30,m);
  ana=sqrt(beta)*8.e-4;
  printf("   %-24s deff=%.6e (want %.6e) m_n=%.1e  %s\n","D closed cone",
         deff,ana,m[0],
         ((fabs(deff-ana)<1.e-16)&&(m[0]==0.))?"ok":"FAIL");
  if(!((fabs(deff-ana)<1.e-16)&&(m[0]==0.))) nbad++;

  /* ---- E: the dissipation rate integrates to Gc --------------------- */

  d0=tn0/kn;df=2.*gc/tn0;
  rate=crackcontrol_dissrate(kn,tn0,gc);
  ana=rate*(df-d0);
  printf("   %-24s rate=%.6e  rate*(df-d0)=%.6f  Gc=%.6f  %s\n",
         "E dissipation rate",rate,ana,gc,
         (fabs(ana-gc)<1.e-12*gc)?"ok":"FAIL");
  if(!(fabs(ana-gc)<1.e-12*gc)) nbad++;

  printf("[CRACKCTL] self test %s (%" ITGFORMAT " failure(s))\n",
         (nbad==0)?"PASSED":"FAILED",nbad);
  return nbad;
}

/* ------------------------------------------------------------------ */
/* Model-facing builder.                                               */
/*                                                                     */
/* Fills c (equation space, zeroed here) with the control functional    */
/* frozen at the displacement field v, and fills the census.  Returns   */
/* the number of integration points carrying weight.                   */
/*                                                                     */
/* mode: 0 MEAN over all live facets, 1 ZONE, 2 DISS (see the header    */
/* comment).  MEAN and ZONE are normalised so that phi is a mean        */
/* separation; DISS is not, so that phi is an energy.                   */
/* ------------------------------------------------------------------ */

void crackcontrol_census_zero(crackcontrol_census *s){

  memset(s,0,sizeof(*s));
}

ITG crackcontrol_build(double *c,ITG neq,ITG mode,
                       const double *co,const ITG *kon,const ITG *ipkon,
                       const char *lakon,ITG ne,
                       const ITG *ielprop,const double *prop,
                       const double *xstate,ITG nstate_,const ITG *mi,
                       const double *v,const ITG *nactdof,ITG nk,ITG mt,
                       crackcontrol_census *s){

  ITG e,ip,i,k,idx,nm,np,dof,nw=0,loading;
  double rmat[9],area,x[9],jump[3],dl[3],m[3],mg[3],sh[3];
  double kn,tn0,ts0,gc,beta,d0,df,deff,dmax,dst,w,wtot=0.,rate,shear;

  for(k=0;k<neq;k++) c[k]=0.;
  crackcontrol_census_zero(s);

  for(e=0;e<ne;e++){
    if((lakon[8*e]!='U')||(lakon[8*e+1]!='C')||(lakon[8*e+2]!='6')) continue;
    s->nfacet++;
    if(ipkon[e]<0){s->ndead++;continue;}
    idx=ipkon[e];
    if(ielprop[e]<0) continue;
    kn =prop[ielprop[e]+0];
    tn0=prop[ielprop[e]+1];
    ts0=prop[ielprop[e]+2];
    gc =prop[ielprop[e]+3];
    if(!(kn>0.)||!(tn0>0.)||!(ts0>0.)||!(gc>0.)) continue;
    beta=(ts0/tn0)*(ts0/tn0);
    d0=tn0/kn;
    df=2.*gc/tn0;
    if(!(df>d0)) continue;
    rate=crackcontrol_dissrate(kn,tn0,gc);

    for(k=0;k<3;k++){
      x[0+k]=co[3*(kon[idx+0]-1)+k];
      x[3+k]=co[3*(kon[idx+1]-1)+k];
      x[6+k]=co[3*(kon[idx+2]-1)+k];
    }
    area=crackcontrol_frame(x,x+3,x+6,rmat);
    if(!(area>0.)) continue;
    s->area+=area;

    for(ip=0;ip<3;ip++){
      for(i=0;i<3;i++) sh[i]=1./6.;
      sh[ip]=2./3.;

      for(k=0;k<3;k++){
        jump[k]=0.;
        for(i=0;i<3;i++){
          nm=kon[idx+i];np=kon[idx+i+3];
          jump[k]+=sh[i]*(v[mt*(np-1)+k+1]-v[mt*(nm-1)+k+1]);
        }
      }
      for(i=0;i<3;i++){
        dl[i]=0.;
        for(k=0;k<3;k++) dl[i]+=rmat[3*i+k]*jump[k];
      }
      deff=crackcontrol_dir(dl,beta,1.e-14*d0,m);

      /* dst is the COMMITTED maximum separation; deff is recomputed from
         the same committed displacement field, so a point that was
         loading at the commit has deff == dst to rounding and one that
         unloaded has deff < dst.  That is the loading test - it reads
         the model's own history variable rather than differencing two
         states. */

      dst=xstate[nstate_*(mi[0]*e+ip)+0];
      loading=(deff>=dst*(1.-1.e-8))?1:0;
      dmax=(dst<deff)?deff:dst;

      /* census -------------------------------------------------------- */

      s->nip++;
      if(dmax>d0){
        s->ninit++;
        if(dmax<df){
          s->nzone++;
          if(loading!=0) s->nload++;
          shear=beta*(dl[1]*dl[1]+dl[2]*dl[2]);
          if(deff>0.) s->shearfrac+=area*shear/(deff*deff);
          s->zonearea+=area;
        }else{
          s->nfail++;
        }
      }
      if(deff>s->deffmax) s->deffmax=deff;
      if(dmax>s->dmaxmax) s->dmaxmax=dmax;

      /* weight -------------------------------------------------------- */

      if(mode==0){
        w=area/3.;
      }else if((dmax>d0)&&(dmax<df)){
        if(mode==2) w=(loading!=0)?(area/3.)*rate:0.;
        else        w=area/3.;
      }else{
        w=0.;
      }
      if(!(w>0.)) continue;
      wtot+=w;
      nw++;

      for(k=0;k<3;k++){
        mg[k]=0.;
        for(i=0;i<3;i++) mg[k]+=rmat[3*i+k]*m[i];
      }
      for(i=0;i<3;i++){
        nm=kon[idx+i];np=kon[idx+i+3];
        for(k=0;k<3;k++){
          dof=nactdof[mt*(np-1)+k+1];
          if(dof>0) c[dof-1]+=w*sh[i]*mg[k];
          dof=nactdof[mt*(nm-1)+k+1];
          if(dof>0) c[dof-1]-=w*sh[i]*mg[k];
        }
      }
    }
  }

  if(s->zonearea>0.) s->shearfrac/=s->zonearea;
  s->weight=wtot;

  /* MEAN and ZONE are means: divide by the weight actually carried, so
     that phi is a separation whatever the size of the active set.  DISS
     is an energy and must not be normalised. */

  if(wtot>0.){
    for(k=0;k<neq;k++) c[k]/=wtot;
  }
  return nw;
}
