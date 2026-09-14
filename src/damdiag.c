/*     CalculiX - damage/fracture extension                              */
/*     damdiag.c: the probes.  They look; they never touch.              */

/* Why this module exists
   ----------------------
   Ten functions that OBSERVE a run - the discrete-branch census, the
   tension/compression event map, the two wall probes, the set difference
   and the A-B-A comparison - were file-statics of nonlingeo.c.  Every one
   of them is a pure function of state handed to it explicitly: not one
   reads a local of the solver, and not one writes anything the solver
   reads back.  They lived inside a 16,000-line function for exactly one
   reason, which is that there was nowhere else to put them.

   What that cost shows up in the measurement.  Before this file existed,
   tools/arch.py reported WALLDIAG with twelve separate sites in
   nonlingeo.c spread over 9,923 lines, and DAMAGE RAY with four over
   11,279.  "What does this probe actually measure" was answerable only by
   reading all of them.

   Contract
   --------
     - everything here is a pure function of its arguments.  It prints and
       it returns numbers; it changes no solver state.  That is the whole
       boundary, and it is why this extraction is free: a probe that cannot
       change an answer cannot change one by being moved;
     - damage_ray_catof() - what an integration point's category IS - is
       file-static, because all three functions that ask are here;
     - the category bits are the vocabulary the CALLER names (nonlingeo
       reads DAMCAT_USOFT and DAMCAT_UADV off a census it was handed), so
       they are declared in CalculiX.h beside the functions that set them.

   What this is NOT
   ----------------
   Not a diagnostic framework, and no attempt at one.  These are the ten
   probes that exist, moved.  Which of them fires, when, and what is
   concluded from it stays where it was.                                */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

/* Number of active damage integration points for the standard 3-D
   continuum elements used by calcdamage.  For uncommon/composite
   formulations fall back to mi[0], i.e. the allocated damage stride. */
void damage_aba_cmp(const char *name,const double *a,
                    const double *b,ITG n,ITG *nbad)
{
  ITG k,first=-1,ndiff=0;
  double mx=0.,d;
  if((a==NULL)||(b==NULL)||(n<=0)){
    printf("   %-14s SKIPPED (null or empty)%s",name,"\n");
    return;
  }
  for(k=0;k<n;k++){
    if(a[k]!=b[k]){
      ndiff++;
      if(first<0) first=k;
      d=fabs(a[k]-b[k]);
      if(d>mx) mx=d;
    }
  }
  if(ndiff==0){
    printf("   %-14s BITWISE IDENTICAL  (n=%" ITGFORMAT ")%s",
           name,n,"\n");
  }else{
    (*nbad)++;
    printf("   %-14s *** DIFFERS *** %" ITGFORMAT " of %" ITGFORMAT
           " entries, first k=%" ITGFORMAT " (%.17e vs %.17e) maxabs=%.6e%s",
           name,ndiff,n,first,a[first],b[first],mx,"\n");
  }
}

static ITG damage_ray_catof(const double *xstate,const double *xstateini,
                            const double *dam,const double *dambase,
                            const double *visc,const double *stx,
                            const char *lakonel,
                            ITG i,ITG j,ITG mi0,ITG nstate)
{
  ITG c=0,ix,is;
  double d,db;
  ix=nstate*mi0*i+nstate*j;
  is=6*mi0*i+6*j;
  if(lakonel[0]=='C'){
    if(nstate>0){
      if(xstate[ix]-xstateini[ix]>1.e-14) c|=DAMCAT_PLAST;
    }
    /* NO real bulk branch flag.  Exporting mattyp into a spare xstate slot
       was attempted and REVERTED: incplas_lin.f:217 owns slots 2..7 for the
       plastic strain tensor, so slot 2 is epl(1), and writing there destroys
       the return map ("no convergence in incplas" at the first increment).
       With *Depvar 4 there is no spare slot at all.  PLAST above stays as
       the observable: peeq grows if and only if the return map ran, so it is
       equivalent by construction, but it IS an inference, not the flag. */
    d=dam[mi0*i+j];
    if(d>=1.) c|=DAMCAT_DINIT;
    if(dambase!=NULL){
      db=dambase[mi0*i+j];
      if(d-db>1.e-14) c|=DAMCAT_DGROW;
    }
  }else if(lakonel[0]=='U'){
    if(nstate>0){ if(xstate[ix]>1.e-14)   c|=DAMCAT_USOFT; }
    if(nstate>1){ if(xstate[ix+1]>1.e-14) c|=DAMCAT_UVISC; }
    if(nstate>3){ if(xstate[ix+3]>0.5)    c|=DAMCAT_UFAIL; }
    /* LOADING/UNLOADING.  cohesive_uc6.f stores dmax=max(dmax0,deff) in slot
       1 and takes the softening tangent only while deff>=dmax0, so
       xstate>xstateini in that slot IS "this point advanced its maximum on
       this trial" - the one active-set indicator of the cohesive law that
       the census did not carry.  A point sitting exactly on the kink has
       xstate==xstateini and is counted as NOT advancing, which is what makes
       a branch crossing along a search direction visible as a transition. */
    if(nstate>0){ if(xstate[ix]>xstateini[ix]) c|=DAMCAT_UADV; }
    /* tension/compression: resultsmech_uc6.f:55 stores traction(1), and both
       branches multiply deltal(1) by a positive factor (kn or g*kn), so the
       sign of the stored traction IS the sign of the normal opening.  This
       switch changes the facet stiffness by 1/g ~ 1e4 and was entirely absent
       from the earlier census. */
    if(stx!=NULL){ if(stx[is]<0.) c|=DAMCAT_UCOMP; }
  }
  return c;
}

void damage_ray_census(ITG *cat,const double *xstate,
                       const double *xstateini,const double *dam,
                       const double *dambase,const double *visc,
                       const double *stx,
                       const ITG *ipkon,const char *lakon,
                       ITG ne0,ITG mi0,ITG nstate)
{
  ITG i,j,nip;
  for(i=0;i<ne0;i++){
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<mi0;j++) cat[mi0*i+j]=0;
    if(ipkon[i]<0) continue;
    for(j=0;j<nip;j++)
      cat[mi0*i+j]=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                                    &lakon[8*i],i,j,mi0,nstate);
  }
}

/* [DAMAGE RAY] the two branch indicators, counted over the live mesh at
   whatever state the caller has just built.  Same iteration as the census
   above, so the two always speak about the same integration points. */

void damage_ray_tally(const double *xstate,const double *xstateini,
                      const double *dam,const double *dambase,
                      const double *visc,const double *stx,
                      const ITG *ipkon,const char *lakon,
                      ITG ne0,ITG mi0,ITG nstate,
                      ITG *nplast,ITG *nucomp)
{
  ITG i,j,nip,c;
  *nplast=0;*nucomp=0;
  for(i=0;i<ne0;i++){
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    if(ipkon[i]<0) continue;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      if(c&DAMCAT_PLAST) (*nplast)++;
      if(c&DAMCAT_UCOMP) (*nucomp)++;
    }
  }
}

/* [DAMAGE EVT] The UC6 tension/compression set.  traction(1) is kn*deltal(1)
   in compression and g*kn*deltal(1) in tension, and both kn and g are
   positive, so sign(stx(1)) IS sign(deltal(1)) and the branch is readable
   without storing deltal anywhere.  UC6 carries exactly three integration
   points. */

void damage_evt_sign(const double *stx,const ITG *ipkon,
                     const char *lakon,ITG ne0,ITG mi0,ITG *sgn)
{
  ITG i,j,np;
  np=(mi0<3)?mi0:3;
  for(i=0;i<ne0;i++){
    for(j=0;j<mi0;j++) sgn[mi0*i+j]=0;
    if(ipkon[i]<0) continue;
    /* UC6 is labelled 'U'; 'C' is the C3D4 bulk (damage_ray_catof splits on
       exactly this).  Getting the letter wrong built the map over the bulk
       and tracked nothing - measured, s3rad inc=231 reported "no ladder alpha
       changes the set" while the ray saw the crossing. */
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++)
      sgn[mi0*i+j]=(stx[6*mi0*i+6*j]<0.)?-1:1;
  }
}

ITG damage_evt_flips(const double *stx,const ITG *ipkon,
                     const char *lakon,ITG ne0,ITG mi0,const ITG *sgn,
                     ITG *firste,ITG *firstip)
{
  ITG i,j,np,n=0,sg;
  np=(mi0<3)?mi0:3;
  *firste=0;*firstip=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++){
      if(sgn[mi0*i+j]==0) continue;
      sg=(stx[6*mi0*i+6*j]<0.)?-1:1;
      if(sg!=sgn[mi0*i+j]){
        if(n==0){ *firste=i+1; *firstip=j+1; }
        n++;
      }
    }
  }
  return n;
}

ITG damage_ray_census_diff(const ITG *cat,const double *xstate,
                           const double *xstateini,const double *dam,
                           const double *dambase,const double *visc,
                           const double *stx,
                           const ITG *ipkon,const char *lakon,
                           ITG ne0,ITG mi0,ITG nstate,
                           ITG *firste,ITG *firstip,
                           ITG *firsta,ITG *firstb)
{
  ITG i,j,nip,c,n=0;
  *firste=-1;*firstip=-1;*firsta=0;*firstb=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      if(c!=cat[mi0*i+j]){
        n++;
        if(*firste<0){
          *firste=i+1;*firstip=j+1;*firsta=cat[mi0*i+j];*firstb=c;
        }
      }
    }
  }
  return n;
}

/* [WALLDIAG] The same census difference, but resolved per category bit and
   split between the two element families, because "the active set moved" and
   "8000 cohesive points crossed their loading/unloading kink" are different
   findings and the aggregate count cannot tell them apart.

   nb[k] counts the integration points whose bit k differs from the reference
   census cat.  Bits are the DAMCAT_ masks; the caller names them. */

/* [WALLDIAG] WHERE the residual and the Newton correction live.

   A norm says how big they are and nothing about what they touch, and the
   two competing explanations of the second wall differ precisely in what
   they touch: a globalisation defect spreads the correction over the mesh,
   a branch-switching one concentrates it on the fracture front.  So the
   probe reports the largest entries by NODE, and for each such node the
   local state of the front - how much of its support is left, how soft its
   cohesive facets are and how many of them are carrying compression.

   nactdof is the post-SPC numbering (1-based, mt-strided) that b, ad and au
   share, so inverting it is what turns an equation index back into a node
   and a direction.  It is inverted here rather than assumed. */

void damage_wall_where(const char *tag,const double *x,ITG neq1,
                       const ITG *nactdof,ITG mt,ITG nk,ITG ntop,
                       const ITG *ipkon,const ITG *kon,
                       const char *lakon,const double *xstate,
                       const double *stx,const double *dam,
                       ITG ne,ITG ne0,ITG mi0,ITG nstate)
{
  ITG *inode=NULL,*idir=NULL,i,j,k,t,ip,np,nb,nu,ncomp,ibest,idx;
  double gmn,gmx,dv,a;

  if((ntop<1)||(neq1<1)) return;
  NNEW(inode,ITG,neq1);
  NNEW(idir,ITG,neq1);
  for(i=0;i<neq1;i++){inode[i]=-1;idir[i]=0;}
  for(i=0;i<nk;i++){
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if((k>0)&&(k<=neq1)){inode[k-1]=i;idir[k-1]=j;}
    }
  }

  for(t=0;t<ntop;t++){
    ibest=-1;a=-1.;
    for(i=0;i<neq1;i++){
      if(inode[i]<0) continue;          /* already reported, or unmapped */
      if(fabs(x[i])>a){a=fabs(x[i]);ibest=i;}
    }
    if(ibest<0) break;
    i=inode[ibest];
    nb=0;nu=0;ncomp=0;gmn=2.;gmx=-1.;
    for(j=0;j<ne;j++){
      if(ipkon[j]<0) continue;
      if(lakon[8*j]=='U') np=6; else if(lakon[8*j]=='C') np=4; else continue;
      ip=0;
      for(k=0;k<np;k++) if(kon[ipkon[j]+k]-1==i) ip=1;
      if(ip==0) continue;
      if(lakon[8*j]=='C'){nb++;continue;}
      nu++;
      if(j>=ne0) continue;
      for(k=0;k<3;k++){
        idx=mi0*j+k;
        if(nstate>1){
          dv=1.-xstate[nstate*idx+1];
          if(dv<gmn) gmn=dv;
          if(dv>gmx) gmx=dv;
        }
        if(stx[6*idx]<0.) ncomp++;
      }
    }
    printf("[WALLDIAG]   %s #%" ITGFORMAT ": node %" ITGFORMAT " dir %"
           ITGFORMAT " value %.6e   live bulk %" ITGFORMAT " live UC6 %"
           ITGFORMAT " facet g in [%.3e,%.3e] UC6 ips in compression %"
           ITGFORMAT "%s",tag,t+1,i+1,idir[ibest],x[ibest],nb,nu,
           (gmn<=1.)?gmn:0.,(gmx>=0.)?gmx:0.,ncomp,"\n");
    inode[ibest]=-1;                       /* do not report it twice */
  }
  SFREE(inode);SFREE(idir);
}

/* [WALLDIAG] WHERE a vector lives, split by what its node touches.

   The linear-model defect at small eps IS eps*(J - dR/du)p, so localising it
   localises the missing term: a defect that sits on nodes with cohesive
   facets accuses the interface tangent, one that sits on nodes whose
   elements are softening accuses the damage rank-1 term, and one spread over
   ordinary plastic bulk accuses the return map.  Shares of |x|_2^2, so they
   add to 1. */

void damage_wall_split(const char *tag,const double *x,ITG neq1,
                       const ITG *nactdof,ITG mt,ITG nk,
                       const ITG *ipkon,const ITG *kon,
                       const char *lakon,const double *dam,
                       const double *xstate,ITG ne,ITG ne0,ITG mi0,
                       ITG nstate)
{
  ITG *touch=NULL,i,j,k,np,idx,nu=0,ns=0,np2=0;
  double tot=0.,su=0.,ss=0.,sp=0.,v;

  NNEW(touch,ITG,nk);
  for(i=0;i<nk;i++) touch[i]=0;
  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]=='U') np=6;
    else if(lakon[8*i]=='C') np=4;
    else continue;
    if(lakon[8*i]=='U'){
      for(k=0;k<np;k++){
        j=kon[ipkon[i]+k]-1;
        if((j>=0)&&(j<nk)) touch[j]|=1;            /* bit0: a live facet   */
      }
    }else{
      idx=mi0*i;
      if((i<ne0)&&(dam[idx]>1.+1.e-12)){
        for(k=0;k<np;k++){
          j=kon[ipkon[i]+k]-1;
          if((j>=0)&&(j<nk)) touch[j]|=2;          /* bit1: softening bulk */
        }
      }
      if((i<ne0)&&(nstate>0)&&(xstate[nstate*idx]>1.e-14)){
        for(k=0;k<np;k++){
          j=kon[ipkon[i]+k]-1;
          if((j>=0)&&(j<nk)) touch[j]|=4;          /* bit2: plastified     */
        }
      }
    }
  }
  for(i=0;i<nk;i++){
    if(touch[i]&1) nu++;
    if(touch[i]&2) ns++;
    if(touch[i]&4) np2++;
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if((k<=0)||(k>neq1)) continue;
      v=x[k-1]*x[k-1];
      tot+=v;
      if(touch[i]&1) su+=v;
      if(touch[i]&2) ss+=v;
      if(touch[i]&4) sp+=v;
    }
  }
  if(tot<=0.) tot=1.;
  printf("[WALLDIAG]   %s lives on: nodes with a live UC6 facet %.6f (%"
         ITGFORMAT " nodes), nodes on softening bulk %.6f (%" ITGFORMAT
         "), nodes on plastified bulk %.6f (%" ITGFORMAT
         "); |x|2=%.6e%s",tag,su/tot,nu,ss/tot,ns,sp/tot,np2,sqrt(tot),"\n");
  SFREE(touch);
}

ITG damage_wall_setdiff(const ITG *cat,const double *xstate,
                        const double *xstateini,const double *dam,
                        const double *dambase,const double *visc,
                        const double *stx,
                        const ITG *ipkon,const char *lakon,
                        ITG ne0,ITG mi0,ITG nstate,ITG *nb)
{
  ITG i,j,nip,c,d,k,n=0;

  for(k=0;k<8;k++) nb[k]=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      d=c^cat[mi0*i+j];
      if(d==0) continue;
      n++;
      for(k=0;k<8;k++) if(d&(1<<k)) nb[k]++;
    }
  }
  return n;
}

/* The values the fifty driver locals carried at their declarations in
   nonlingeo().  Every probe here is opt-in; none of them is armed unless a
   switch says so, which is why the defaults are almost all zero. */
void probedrv_init(probedrv *p)
{
  memset(p,0,sizeof(*p));
  p->ray_cat=NULL;
  p->ray_growth=1.10;
  p->ray_incok=1;
  p->ray_max=8;
  p->ray_p=NULL;
  p->ray_r0=NULL;
  p->ray_res=NULL;
  p->wall_cat=NULL;
  p->wall_def=NULL;
  p->wall_theta=-1.;
  p->aba_a=0.25;
  p->evt_sgn=NULL;
  p->dump_inc=1;
  p->null_nit=4;
  p->null_seed=987654321;
  p->null_x=NULL;
}

/* ---- arming ------------------------------------------------------------

   The A-B-A purity test reads its own switches.  It is a diagnostic that
   fires once and then lets the run continue from the full step, so nothing
   downstream depends on when it is configured - but it is read here, at
   the point it was read before, because that is the rule for all of these
   and an exception would have to be argued rather than assumed. */

void probedrv_configure_aba(probedrv *p)
{
  const char *e;

  /* CCX_DAMAGE_RELEASE_PROBE - PURE DIAGNOSTIC, reads only.

     How much internal force does a topology event actually release?
     The run log cannot answer it: "largest residual force" is printed
     AFTER Newton has already made its first correction, so it is what
     survived the release, not the release.

     The probe differences f_int(u*) across the event with the
     displacement state held fixed, and - the reason it exists - splits
     the result by what happened to the equation:

       surv     DOF active BEFORE and AFTER.  The equation still exists,
                so this is the only place a perturbation of the system
                Newton solves can live.  dF_surv_max/qam is the number.
       removed  DOF active BEFORE, gone AFTER.  The equation does not
                exist any more; Newton neither resolves it nor owes it
                anything.  Diagnosis only, never a criterion.

     Reporting one number for both is what made the raw residual
     unreadable in the first place.

     =1 per batch, =2 adds one line per deleted element.  Unset = off and
     nothing is allocated.  This flag CHANGES NO BIT OF THE ANSWER, and
     that is gated both ways: s0_coarse_ts must give m.damage md5
     74212e957d7cd649 with the probe off AND with it on. */
  /* CCX_DAMAGE_RESIDUAL_RAY - PURE DIAGNOSTIC, reads only.

     J-13 established that the force RELEASED by a topology event does not
     order fatal against non-fatal: the fatal event ranked 487th of 563 on
     bandrad, and 486 larger releases were survived.  So the perturbation
     is not the right-hand side.  What is left is the STEP: the residual
     contracts three times (0.265 -> 0.019 -> 0.0022) and then flies up by
     265x on the next full Newton step, at frozen load and frozen topology.
     That is the classic signature of a full step leaving the basin, and
     the standard instrument for it is a scan of the residual along the
     Newton direction.

     There is no such scan in this tree, and there cannot be a line search
     either: BK3 is excluded from re-equilibration by TWO independent
     gates - its own conjunction (idamagereeq==0, below) and the resold
     store, which is guarded by the same conjunction, so resold is stale
     throughout an idamagereeq pass.  This probe needs neither: it
     evaluates alpha=0 itself and uses that as the reference.

     What it does: at a re-equilibration iteration it evaluates the
     residual at alpha=0 and alpha=1.  If the full step grew the norm by
     more than the growth factor - the fatal signature - it walks the
     whole ray, repeats one alpha to prove the evaluation is a pure
     function of the step length, and ends at alpha=1, which is exactly
     the state the unprobed code would have had.

     PURITY IS THE GATE, NOT AN ASIDE.  The repeated alpha must reproduce
     bitwise.  It can only do so for the residual vector, at fixed dtime,
     with CCX_DAMAGE_NONLOCAL unset and no contact: damjac/xstiff are not
     rebuilt from a baseline, and the nonlocal field (dpsave/ebar) is
     trial-derived, never snapshotted and CG-warm-started to 1e-10.  A
     mismatch under those conditions is a real impurity, not a bug in the
     probe, and it would void any line search built on top.

     Value = max number of rays to walk (default 8).  Unset = off,
     nothing allocated.  Changes no bit of the answer: b is saved and
     restored exactly, and the last evaluation is the alpha=1 state. */
  /* CCX_DAMAGE_REEQ_BACKTRACK - SOLVER CHANGE, not a diagnostic.
     Damps the Newton step during same-load re-equilibration, restoring
     the committed baseline before every probe and restoring the full
     step when nothing is acceptable.  THIS CHANGES THE ANSWER: a run
     that goes further with it is not thereby a success, and adopting it
     needs the full verify + ladder gate.  Default off = bit-identical. */
  /* CCX_DAMAGE_ABA=<alpha> - PURE DIAGNOSTIC.  Proves, or refutes, that a
     trial evaluation is a pure function of the step length.  The earlier
     A-B-A compared ONE scalar (|R|inf); one scalar agreeing proves
     nothing about the rest of the state, and the backtracking snapshot
     only held dam/damvisc/xstate while results() writes more than that.
     This evaluates A, snapshots EVERY array results()/calcresidual
     touch, evaluates B, evaluates A again, and compares byte for byte.
     Fires once, then the run continues from the full step. */
  if((e=ccxopt_getenv("CCX_DAMAGE_ABA"))!=NULL){
    p->aba_mode=1;
    p->aba_a=atof(e);
    if((p->aba_a<=0.)||(p->aba_a>=1.)) p->aba_a=0.25;
    /* CCX_DAMAGE_ABA_INC="143,209" - fire at those increments instead of
       at the first opportunity.  Purity proved on one activated path does
       not prove it on another: a different increment reaches the same
       code through a different constitutive state, and that is exactly
       what has to be shown before an ensemble rests on it. */
    p->aba_ninc=0;
    if((e=ccxopt_getenv("CCX_DAMAGE_ABA_INC"))!=NULL){
      char *acp=e;
      while((*acp!=0)&&(p->aba_ninc<4)){
        while((*acp==' ')||(*acp==',')) acp++;
        if(*acp==0) break;
        p->aba_inc[p->aba_ninc++]=atoi(acp);
        while((*acp!=0)&&(*acp!=',')) acp++;
      }
    }
    printf("[DAMAGE ABA] DIAGNOSTIC: full-state A-B-A at alpha=%.6f;\n"
           "   every array written by results()/calcresidual is compared\n"
           "   byte for byte between two evaluations at the same alpha.%s",
           p->aba_a,"\n");
    fflush(stdout);
  }
}
