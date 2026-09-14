/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2025 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "CalculiX.h"
#include "mortar.h"
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Damage Patch A3 adaptive-fast controller.
   The values below only affect trial load placement.  No trial damage is
   committed unless the post-deletion Newton/active-set solve converges. */
#define DAMAGE_FAST_BATCH_MAX 24
#define DAMAGE_FAST_MIN_FRACTION 0.25
#define DAMAGE_FAST_DIRECT_ALPHA 0.35

/* DE1.1 fixed-point controller.
   DE1 damage is coupled to equilibrium by same-load staggered iterations.
   Requiring bitwise/no-change closure is unnecessarily strict and can make
   one physical increment repeat indefinitely.  Use an engineering tolerance
   on the largest element damage correction, allow a slightly relaxed closure
   at the pass limit, and otherwise cut the physical increment back. */
/* the four now live in CalculiX.h, beside damstats. */

/* DE1.3 terminal-failure controller.
   Continuous DE1.2 degradation remains fully Newton-integrated.  Topology
   changes are allowed only very near complete material failure.  Limit one
   topology update to a modest batch; additional terminal elements are picked
   up by the existing same-load active-set closure after redistribution. */
/* bounds on how fast the dissipation-based step control may change
   the increment: shrinking hard is safe, growing back has to be
   gentle or the controller oscillates across the event */
#define DAMAGE_DISS_GROW 1.25
#define DAMAGE_DISS_SHRINK 0.20
/* largest lambda move one iteration may make, in units of the
   nominal increment: a bad secant must not throw the load factor
   across the whole step */
#define DAMAGE_DISS_DLAM 0.50
/* fraction of the target dissipation an accepted increment has to
   reach before the load factor is handed to the constraint */
#define DAMAGE_DISS_ENGAGE 0.10

#define DAMAGE_DE13_DELETE_D 0.999
#define DAMAGE_DE13_BATCH_MAX 64

/* a node whose entire live support is below this fraction of its original
   stiffness is treated as free for stabilisation purposes (E-61) */
#define DAMAGE_STAB_GDEAD 1.e-2

/* Damage-aware slow-Newton controller.
   The stock convergence controller starts estimating convergence at ir and
   cuts the increment at ic (normally 16), even when the residual and the
   displacement correction are both contracting monotonically.  DE1.2 uses
   an approximate damaged tangent, so linear convergence is expected while D
   evolves.  Permit a bounded extension only for that state.  The ordinary
   checkconvergence() residual/correction criteria and all of its divergence
   paths remain unchanged. */
/* the three numbers now live in CalculiX.h, next to the decision they
   bound; the fourth, D_TOL, went with erosion_softening(). */

/* Damage-specific Newton globalization.
   A full Newton correction is retained whenever it contracts the force
   residual.  Only a genuine residual increase in an active DE1.2/DM2.0
   softening trial triggers one safeguarded secant line-search correction.
   This adds a constitutive/residual evaluation, but no PARDISO
   factorization.  Final convergence tolerances remain the stock ones. */
/* the five now live in CalculiX.h, beside the rescue state. */

/* AUTOSPC.  Nodes whose assembled diagonal has collapsed relative to their
   own intact value are excluded from the DISPLACEMENT convergence norm.
   Read by resultsini.c, which is where the Newton correction is still in
   hand - by the time control returns to nonlingeo, b holds residual forces
   and cam[0] cannot be recomputed.
   NULL, or a zero length, means the feature is off.

   WHY THIS IS SAFE.  The excluded DOF is still solved, still updated, and
   still measured by the FORCE residual ram[0], which is untouched.  Only its
   veto over cam[0] is removed.  A node that still carries load therefore
   still blocks convergence through the force criterion; what is dropped is a
   correction whose size is set by conditioning rather than by physics.
   Nothing is deleted, which is what separates this from E-22, E-57 and E-72.

   WHY IT IS NEEDED.  checkconvergence.c gates on cam[0] <= c2*uam[0] with
   c2 = 0.01, and on m12_epsf50 the run died with a force residual 338x
   inside tolerance and a correction/increment ratio of 0.29, both reported
   at a node with no stiffness left (E-67).  62-68% of the blocking reports
   in the final phase name a topologically defective node. */
ITG *damage_spc_mask=NULL;
ITG damage_spc_nk=0;
ITG damage_spc_count=0;

/* [DAMSTATE] the single owner of the load-path judgement; the three
   globals above are views onto it once it has been updated. */
damstate damage_dstate={0,NULL,NULL,NULL,NULL,0,0.,0};





/* the three rescue flags are defined in rescue.c now, beside the state
   they belong to; CalculiX.h declares them. */


































void nonlingeo(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	       ITG *ne,
	       ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	       ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
	       ITG *nmpc,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	       ITG **nelemloadp,char **sideloadp,double *xload,ITG *nload,
	       ITG *nactdof,
	       ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	       ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,double *t1old,
	       ITG *ithermal,double *prestr,ITG *iprestr,
	       double **voldp,ITG *iperturb,double *sti,ITG *nzs, 
	       ITG *kode,char *filab,
	       ITG *idrct,ITG *jmax,ITG *jout,double *timepar,
	       double *eme,
	       double *xbounold,double *xforcold,double *xloadold,
	       double *veold,double *accold,
	       char *amname,double *amta,ITG *namta,ITG *nam,
	       ITG *iamforc,ITG **iamloadp,
	       ITG *iamt1,double *alpha,ITG *iexpl,
	       ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
	       ITG *nplkcon,
	       double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
	       char *matname,double *qaold,ITG *mi,
	       ITG *isolver,ITG *ncmat_,ITG *nstate_,
	       double *cs,ITG *mcs,ITG *nkon,double **enerp,ITG *mpcinfo,
	       char *output,
	       double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
	       double *physcon,ITG *nflow,double *ctrl,
	       char *set,ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
	       ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
	       ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
	       ITG *ielprop,double *prop,ITG *ntie,char *tieset,
	       ITG *itpamp,ITG *iviewfile,char *jobnamec,double *tietol,
	       ITG *nslavs,double *thicke,ITG *ics,
	       ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	       ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	       ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext,ITG *nprop,
	       ITG *network,char *orname,double *vel,ITG *nef,
	       double *velo,double *veloo,double *energy,ITG *itempuser,
	       ITG *ipobody,ITG *inewton,double *t0g,double *t1g,
	       ITG *ifreebody,ITG *nlabel,ITG *ndmat_,ITG *ndmcon,
	       double *dmcon,double *dam){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
    *sideface=NULL,*labmpc=NULL,*lakonf=NULL,*env,*envsys,fneig[132]="",
    *sideloadref=NULL,*sideload=NULL,stiffmatrix[132]="",
    damagefilename[160]="",*sideloadf=NULL,cflag[1]=" ",
    *damage_tangent_env=NULL,*damage_topology_env=NULL,
    *damage_de13_env=NULL,
    *damage_visc_env=NULL,
    *damage_fracture_env=NULL,
    *damage_fracture_seta=NULL,*damage_fracture_setb=NULL;
  char damage_fracture_a[81],damage_fracture_b[81]; 
 
  ITG *inum=NULL,k,l,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
    icutb=0,istab=0,uncoupled,n1,n2,itruecontact=1,iclean=0,
    iperturb_sav[2],iforbou,*icol=NULL,*irow=NULL,ielas=0,icmd=0,
    memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
    *nodempcref=NULL,memmpcref_,mpcfreeref,*itg=NULL,*ineighe=NULL,
    *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,idamping=0,
    *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,im,
    ntrit,*inocs=NULL,*nacteq=NULL,*ipface=NULL,masslesslinear=0,
    *nactdog=NULL,nteq,*itietri=NULL,*koncont=NULL,istrainfree=0,
    ncont,ne0,nkon0,*ipkon=NULL,*kon=NULL,*ielorien=NULL,
    *ielmat=NULL,itp=0,symmetryflag=0,inputformat=0,kscale=1,
    *iruc=NULL,iitterm=0,iturbulent,ngraph=1,ismallsliding=0,
    *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
    *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,*imastnode=NULL,
    *nslavnode=NULL,*nmastnode=NULL,*imastop=NULL,imat,
    *iponoels=NULL,*inoels=NULL,*islavsurfold=NULL,maxlenmpcref,
    *islavact=NULL,mt=mi[1]+1,*nactdofinv=NULL,*ipe=NULL, 
    *ime=NULL,*ikactmech=NULL,nactmech,inode,idir,neold,neini,
    iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
    *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,negpres=0,
    mass[2]={0,0},stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
    coriolis=0,*ipneigh=NULL,*neigh=NULL,maxprevcontel,nslavs_prev_step,
    *nelemface=NULL,*ipoface=NULL,*nodface=NULL,*ifreestream=NULL,
    *isolidsurf=NULL,*neighsolidsurf=NULL,*iponoeln=NULL,*inoeln=NULL,
    nface,nfreestream,nsolidsurf,i,icfd=0,id,nslavquadel=0,
    node,networknode,iflagact=0,*nodorig=NULL,*ipivr=NULL,iglob=0,
    *inomat=NULL,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,nforcrhs,nloadrhs,
    idampingwithoutcontact=0,*nactdoh=NULL,*nactdohinv=NULL,*ipkonf=NULL,
    *ielmatf=NULL,*ielorienf=NULL,ialeatoric=0,nloadref,isym,
    *nelemloadref=NULL,*iamloadref=NULL,*idefload=NULL,nload_,
    *nelemload=NULL,*iamload=NULL,ncontacts=0,inccontact=0,nrhs=1,
    j=0,inoelnsize=0,isensitivity=0,*konf=NULL,nbodyrhs,
    *iwork=NULL,nelt,lrgw,*igwk=NULL,itol,itmax,iter,ierr,iunit,ligw,
    mei[4]={0,0,0,0},*itreated=NULL,mscalmethod=-1,inoelfree,
    isiz=0,num_cpus,sys_cpus,ne1d2d=0,kchdep,nkftot,
    ifreesurface=0,*iponoelf=NULL,*inoelf=NULL,*iponoel=NULL,
    *damage_iponoel_trial=NULL,*damage_orphan_seen=NULL,*damage_addok=NULL,
    mortartrafoflag=0,*nelold=NULL,*nelnew=NULL,*nkold=NULL,*nknew=NULL,
    *ipompcf=NULL,*nodempcf=NULL,*nodebounf=NULL,*ndirbounf=NULL,
    *nelemloadf=NULL,*ipobodyf=NULL,nkf,nkonf,memmpcf,nbounf,nloadf,nmpcf,
    *ikbounf=NULL,*ilbounf=NULL,*ikmpcf=NULL,*ilmpcf=NULL,*iambounf=NULL,
    *iamloadf=NULL,*inotrf=NULL,*jqw=NULL,*iroww=NULL,nzsw,*jqtherm=NULL,
    *kslav=NULL,*lslav=NULL,*ktot=NULL,*ltot=NULL,nmasts,neqtot,
    intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
    intpointvart,*jqbi=NULL,*irowbi=NULL,*jqib=NULL,*irowib=NULL,
    idispfrdonly,*inumcp=NULL,nmethodold=*nmethod,
    idamage=0,iitsav=0,idamagereeq=0,ilocalsubstep=0,*ipkondamageini=NULL,
    *damage_de13_trigger_ip=NULL,*damage_ract=NULL,damage_dth_n=0,damage_dth_i=0,
    damage_release_probe=0,damage_release_armed=0,damage_release_pass=0,
    damage_release_rebuild=0,
    damage_release_nterm=0,damage_release_nother=0,
    damage_release_nisl=0,damage_release_ncoh=0,damage_release_iforbou=0,
    damage_batch=0,damage_scan_count=0,damage_nip_local=0,
    damage_mode=0,damage_predict_count=0,damage_event_cut=0,
    damage_active_pass=0,damage_soft_reeq=0,damage_fast_retry=0,damage_fast_used=0,
    damage_fast_recover=0,damage_fast_failures=0,
    damage_de12_enabled=0,damage_de12_matcount=0,damage_dm20_matcount=0,
    damage_tangent_mode=0,
    damage_cut_on=0,damage_cut_bad=0,damage_cut_narrow=0,
    damage_cut_exact=0,
    damage_topology_deferred_mode=0,damage_topology_rebuild=1,*damage_damcat=NULL,
    damage_snap_elem=0,damage_snap_bad=0,damage_batch_list=0,
    damage_float_new=0,damage_float_reach=0,damage_float_total=0,
    damage_float_coh=0,damage_float_isl=0,
    damage_conn=1,damage_conn_reach=0,damage_fracture_complete=0,
    /* bounded exactly like the DE1.3 terminal batch: once the hydride is
       gone every facet on its surface loses its plus side at the same
       instant, and removing all of them in one topology transaction is a
       far larger redistribution than the same-load solve can absorb */
    damage_float_batch=DAMAGE_DE13_BATCH_MAX,
    damage_dangle_new=0,damage_dangle_weak=0,damage_dangle_total=0,
    damage_dangle_max=0,damage_stiff_probe=0,
    damage_stab_maxdof=0,damage_stab_maxdead=0,*damage_stab_node=NULL,
    damage_fracture_link=0,damage_deadfacet=0,damage_facetdel=0,
    damage_facetdel_new=0,damage_facetdel_total=0,damage_spc_neg=0,damage_census_ok=1,
    damage_bare=0,damage_bare_rep=-1,damage_free_probe=0,
    damage_free_cnt=0,damage_free_rep=-1,*damage_free_nb=NULL,
    damage_free_worst=-1,damage_free_raw=0,damage_free_rawrep=-1,
    *damage_stiff_haz=NULL,
    damage_topology_orphans=0,damage_indexe=0,damage_lsr=0,
    damage_de13_transaction=0;

  double damage_cut_frac=0.,damage_cut_ref=-1.,damage_cut_now=0.,
    damage_cut_gmin=1.e-4;
  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cg=NULL,
    *cdn=NULL,*pslavsurfold=NULL,*fextload=NULL,
    *f=NULL,*fn=NULL,qa[4]={0.,0.,-1.,0.},qam[2]={0.,0.},dtheta,theta,
    err,ram[6]={0.,0.,0.,0.,0.,0.},*areaslav=NULL,
    *springarea=NULL,ram1[6]={0.,0.,0.,0.,0.,0.},
    ram2[6]={0.,0.,0.,0.,0.,0.},deltmx,ptime,smaxls,sminls,
    uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,*straight=NULL,
    *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
    *xforcact=NULL,*xloadact=NULL,*fext=NULL,*clearini=NULL,
    reltime,time,bet=0.,gam=0.,*aux2=NULL,dtime,*fini=NULL,
    *fextini=NULL,*veini=NULL,*accini=NULL,*xstateini=NULL,
    *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,pressureratio,
    *xbounini=NULL,dev,*xstiff=NULL,*stx=NULL,*stiini=NULL,
    *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
    *coefmpcref=NULL,*enerini=NULL,*emn=NULL,alpham,betam,
    *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
    *adview=NULL,*auview=NULL,*qfx=NULL,*cvini=NULL,*cv=NULL,
    *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
    *xbodyact=NULL,*cgr=NULL,dthetaref,dthetadamage,thetadamage,
    dthetarefdamage,theta_goal=0.,theta_local_start=0.,
    dtheta_restore=0.,dtheta_remaining, *vr=NULL,*vi=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
    *f_cm=NULL, *f_cs=NULL,*adc=NULL,*auc=NULL,*res=NULL,
    *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
    *xmastnor=NULL,*emeini=NULL,*tinc,*tper,*tmin,*tmax,*tincf,
    *doubleglob=NULL,*xnoels=NULL,*au=NULL,*resold=NULL,
    *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
    *x=NULL,*y=NULL,*z=NULL,*xo=NULL,sum1,sum2,flinesearch,
    *yo=NULL,*zo=NULL,*cdnr=NULL,*cdni=NULL,*fnext=NULL,*fnextini=NULL,
    allwk=0.,allwkini,energyini[4]={0.,0.,0.,0.},*cof=NULL,
    energyref,dtcont,dtvol,wavespeed[*nmat],emax,r_abs,
    enetoll,dampwk=0.,dampwkini=0.,temax,*tmp=NULL,energystartstep[4],
    sizemaxinc,*adblump=NULL,*adcpy=NULL,*aucpy=NULL,*rwork=NULL,
    *sol=NULL,*rgwk=NULL,tol,*sb=NULL,*sx=NULL,delcon,alea,
    *smscale=NULL,dtset,energym=0.,energymold=0.,*voldf=NULL,
    *coefmpcf=NULL,*xbounf=NULL,*xloadf=NULL,*xbounoldf=NULL,
    *xbounactf=NULL,*xloadoldf=NULL,*xloadactf=NULL,*auw=NULL,*volddof=NULL,
    *qb=NULL,*aloc=NULL,dtmin,*fric=NULL,*aubi=NULL,*auib=NULL,
    *fullgmatrix=NULL,*fullr=NULL,*alglob=NULL,*damn=NULL,*errn=NULL,
    *damdamageini=NULL,*damde1prev=NULL,*veolddamageini=NULL,
    *damage_de13_trigger_value=NULL,*damagebase=NULL,
    *damage_damjac=NULL,damage_snap_ratio=0.,
    *damage_damvisc=NULL,*damage_damviscini=NULL,damage_visc_eta=0.,
    *damage_frel=NULL,
    damage_dtheta_healthy=0.,
    damage_dth_ring[20]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
    *daba_res=NULL,*daba_v=NULL,*daba_stx=NULL,*daba_fn=NULL,
    *daba_f=NULL,*daba_dam=NULL,*daba_visc=NULL,*daba_xs=NULL,
    *daba_eme=NULL,*daba_stiff=NULL,*daba_qa=NULL,*daba_cam=NULL,
    damage_release_qa=0.,damage_release_qam=0.,
    damage_release_dt=0.,
    *damage_addiag=NULL,*damage_addiag0=NULL,
    damage_stiff_min=0.,
    *damage_free_g=NULL,damage_free_gm=0.,damage_free_dv=0.,
    damage_dmax=0.,
    damage_alphaevent=2.,damage_event_dtheta=0.,
    damage_event_raw=0.,damage_event_floor=0.;
  ITG damage_spc_force=0;

  /* [LOADCTL] who drives the load parameter: fifty-three locals that
     had no owner, and one mutual-exclusion rule.  See loadctl.c. */
  loadctl lc;

  /* [OPCHECK] and [WALLDIAG]/[DAMAGE RAY]/[DAMAGE ABA]: the driver
     state of the operator check and of the probes.  See opcheck.c
     and damdiag.c. */
  opcheckdrv opd;
  probedrv prb;

  /* [RESCUE] what happens when an increment will not converge:
     sixty-seven locals that had no owner.  See rescue.c. */
  rescue rsc;

  /* the Newton iteration budget, and the damage census it reads */
  slownewton slow;
  damstats de1;

  /* [PATHFOLLOW] the driver's state: sixty-three locals that had no
     owner.  The method itself is pathfollow.c's. */
  pathdrv pf;

  /* [DAMAGE CT] the continuation's own state: eighty-eight locals
     that had no owner.  See damcont.c. */
  damcont ct;

  /* [DAMAGE TR] the trust region's own state: forty-nine locals
     that had no owner.  See dogleg.c. */
  dogleg dog;

  /* [TRIAL] the residual evaluator's view of this frame; see trial.c */
  trialctx nlgt;

  /* [TOPOLOGY] the erosion transaction.  Nine locals with no owner became
     one object with one lifetime.  Six copies of "discard the marked set", in three
     different variants, are now one call. */
  topo_txn dtxn;

  /* [GLOBALIZE] the census of which globalization mechanism ever changes
     anything.  Counts and reports, decides
     nothing, changes no arithmetic. */
  glob_census damage_glob;
  double damage_nl_ell=0.;
  ITG damage_nl_mode=0;
  double damage_qam_floor=0.;
  converge damage_cvg;
  double damage_stab_alpha=0.,damage_spc_g=0.;
  /* [EROSION] who leaves the assembly, and what the run has taken so
     far: two objects instead of twelve locals with no owner. */
  erosion_policy damage_epol={DAMAGE_DE13_DELETE_D,1,NULL,0.,0.,
                              DAMAGE_DE13_BATCH_MAX};
  erosion_batch damage_ebatch;
  char *damage_stab_env=NULL,*damage_deadsole_env=NULL,
    *damage_deadall_env=NULL;

  /* ---- CCX_PATHFOLLOW: consistent dissipation path following ----------
     Everything below is inert unless CCX_PATHFOLLOW is set to a positive
     dissipation increment.  The state proper lives in pathfollow.c; these
     are only the handles the Newton loop needs. */

  /* ---- topology / deletion diagnostics (topodiag.c) -----------------
     Measurement only.  td_trace prints a content hash of the COMMITTED
     state at the top of every attempt and the sorted composition of every
     deletion batch, which is what turns "the same batch repeats from the
     same state" from an inference into a measurement.  td_from runs the
     connectivity/rank/residual report from that increment on. */

  ITG td_trace=0,td_from=0,*td_comp=NULL,*td_sort=NULL,td_armed=0;
  ITG td_nmode=0,td_qneq=0;
  double *td_qkeep=NULL;
  unsigned long long td_state=0ULL,td_batch=0ULL;
  topodiag_report td_rep;

  lsladder damage_lsl;

  /* ---- mixed-mode crack control (pf.codmode==2) ---------------------
     The control functional is rebuilt from the COMMITTED state at the top
     of every attempt and the constraint is written incrementally, so a
     redefinition between increments carries nothing over.  See
     crackcontrol.c. */

	 
  FILE *f1,*fdamage=NULL;

#ifdef SGI
  ITG token;
#endif
  
  /* declarations for mortar contact */

  ITG *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc=NULL,
    *islavactdof=NULL,*islavactini=NULL,*islavtie=NULL,
    *irowt=NULL,*jqt=NULL,*irowtinv=NULL,*jqtinv=NULL,
    *irowb=NULL,*jqb=NULL,*irowd=NULL,*jqd=NULL,*irowdtil=NULL,*jqdtil=NULL,
    *irowbtil=NULL,*jqbtil=NULL,*irowbhelp=NULL,*jqbhelp=NULL,
    *islavnodeinv=NULL,*islavquadel=NULL,*irowc2=NULL,*jqc2=NULL,nzsc2,
    *icolc2=NULL,
    *jqbd=NULL,*irowbd=NULL,*jqbdtil=NULL,*irowbdtil=NULL,*jqbdtil2=NULL,
    *irowbdtil2=NULL,
    *jqdd=NULL,*irowdd=NULL,*jqddtil=NULL,*irowddtil=NULL,*jqddtil2=NULL,
    *irowddtil2=NULL,
    *jqddinv=NULL,*irowddinv=NULL,*jqtemp=NULL,*irowtemp=NULL,*icoltemp=NULL,
    nzstemp[3];
  
  double *bp=NULL,*gap=NULL,*slavnor=NULL,
    *slavtan=NULL,*cdisp=NULL,*cstress=NULL,*cfs=NULL,*cfm=NULL,*cfsini=NULL,
    *cfstil=NULL,*bpini=NULL,
    *cstressini=NULL,*pslavdual=NULL,*aut=NULL,
    *autinv=NULL,*Bd=NULL,*Bdhelp=NULL,
    *Dd=NULL,*Ddtil=NULL,*Bdtil=NULL,*auc2=NULL,*adc2=NULL,*aubd=NULL,
    *audd=NULL,*auddtil=NULL,*auddtil2=NULL,*auddinv=NULL,*bhat=NULL,
    *aubdtil=NULL,*aubdtil2=NULL;

  /* end of declarations for mortar contact */

  icol=*icolp;irow=*irowp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;nelemload=*nelemloadp;
  iamload=*iamloadp;sideload=*sideloadp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  /* [TRIAL] bind the residual evaluator to this frame.  Every field is
     the ADDRESS of a local, so this is valid from here to the end of
     the function no matter how often the arrays are reallocated. */
  TRIAL_BIND(nlgt);
  dogleg_init(&dog);
  damcont_init(&ct);
  pathdrv_init(&pf);
  rescue_init(&rsc);
  opcheckdrv_init(&opd);
  probedrv_init(&prb);
  loadctl_init(&lc);
  slownewton_init(&slow);
  damstats_init(&de1);

  /* determining whether a node belongs to at least one element
     (needed in resultsforc.c) */
  
  NNEW(iponoel,ITG,*nk);
  FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));

  if(filab[4]!=' ') ne1d2d=1;

  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
    sys_cpus=atoi(envsys);
    if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
    sys_cpus=getSystemCPUs();
    if(sys_cpus<1) sys_cpus=1;
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
    if(env)
      num_cpus=atoi(env);
    if(num_cpus<1) {
      num_cpus=1;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }
  
  // MPADD: initialize rmin to the tolerance
  enetoll=0.02;
  r_abs=0.0;
  emax=0.0;
  // MPADD end

  delcon=ctrl[53];alea=ctrl[54];

  tinc=&timepar[0];
  tper=&timepar[1];
  tmin=&timepar[2];
  tmax=&timepar[3];
  tincf=&timepar[4];

  if(*ithermal==4){
    uncoupled=1;
    *ithermal=3;
  }else{
    uncoupled=0;
  }

  /* for massless explicit dynamics any other "nonlingeo" step in the same
     calculation (e.g. a static step or an implicit dynamics step) 
     is performed with node-to-face contact */
  
  if(*mortar!=1){
    if(*nintpoint!=0){
      maxprevcontel=0;
      *nintpoint=0;
      SFREE(pslavsurf);SFREE(clearini);
    }else{
      maxprevcontel=*nslavs;
    }
  }else if(*mortar==1){
    maxprevcontel=*nintpoint;
    if(*nstate_!=0){
      if(maxprevcontel!=0){
	MNEW(islavsurfold,ITG,2**ifacecount+2);
	MNEW(pslavsurfold,double,3**nintpoint);
	isiz=2**ifacecount+2;cpyparitg(islavsurfold,islavsurf,&isiz,&num_cpus);
	isiz=3**nintpoint;cpypardou(pslavsurfold,pslavsurf,&isiz,&num_cpus);
      }
    }
    nslavs_prev_step=*nslavs;
  }

  /* turbulence model 
     iturbulent==0: laminar
     iturbulent==1: k-epsilon
     iturbulent==2: q-omega
     iturbulent==3: BSL
     iturbulent==4: SST */
  
  iturbulent=(ITG)physcon[8];
  
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];

  sminls=ctrl[28];smaxls=ctrl[29];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];

  alpham=xmodal[0];
  betam=xmodal[1];

  /* check whether, for a dynamic calculation, damping is involved */
  
  if(*nmethod==4){
    if(*iexpl<=1){
	  
      /* implicit dynamics */
	  
      if((fabs(alpham)>1.e-30)||(fabs(betam)>1.e-30)){
	idamping=1;idampingwithoutcontact=1;
      }else{
	for(i=0;i<*ne;i++){
	  if(ipkon[i]<0) continue;
	  if(strcmp1(&lakon[i*8],"ED")==0){
	    idamping=1;idampingwithoutcontact=1;break;
	  }
	}
      }
    }else{
	  
      /* explicit dynamics */
	  
      if((fabs(alpham)>1.e-30)||((fabs(betam)>1.e-30)&&(*mortar==-1))){
	idamping=1;idampingwithoutcontact=1;
      }
      if((fabs(betam)>1.e-30)&&(*mortar!=-1)){
	printf
	  (" *ERROR in nonlingeo: in explicit dynamic calculations\n");
	printf
	  ("         without massless contact the damping is only\n");
	printf
	  ("         allowed to be mass proportional: the coefficient beta\n");
	printf("         of the stiffness proportional term must be zero\n");
	FORTRAN(stop,());
      }
    }
  }
  
  /* check whether a sensitivity step may follow (whether design variables
     were defined */

  for(i=0;i<*ntie;i++){
    if(strcmp1(&tieset[i*243+80],"D")==0){
      isensitivity=1;
      NNEW(adcpy,double,neq[1]);
      /* no asymmetric matrices allowed for sensitivity */
      NNEW(aucpy,double,nzs[1]);
      break;
    }
  }

  if((icascade==2)&&(*iexpl>1)){
    printf
      (" *ERROR in nonlingeo: linear and nonlinear MPC's depend on each other\n");
    printf("        This is not allowed in a explicit dynamic calculation\n");
    FORTRAN(stop,());
  }
      
  /* determining the global values to be used as boundary conditions
     for a submodel */

  ITG irefine=0;
  getglobalresults(&jobnamec[396],&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
		   ithermal,nk,t1,iamt1,&sigma,&irefine);
  
  if(iglob<0){
    printf(" *ERROR in nonlingeo: a submodel calculation for which\n");
    printf("        the global model results from a *FREQUENCY\n");
    printf("        calculation must be geometrically linear\n");
    FORTRAN(stop,());
  }

  /* reading temperatures from frd-file */
  
  if((itempuser[0]==2)&&(itempuser[1]!=itempuser[2])) {
    utempread(t1,&itempuser[2],jobnamec);
  }      
  
  /* invert nactdof */
  
  /*  NNEW(nactdofinv,ITG,mt**nk);
  MNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
			 SFREE(nodorig);*/
  
  /* allocating a field for the stiffness matrix */
  
  NNEW(xstiff,double,(long long)27*mi[0]**ne);
  
  /* allocating force fields */
  
  NNEW(f,double,neq[1]);
  NNEW(fext,double,neq[1]);
  
  NNEW(b,double,neq[1]);
  NNEW(vini,double,mt**nk);
  
  NNEW(aux,double,7*maxlenmpc);
  NNEW(iaux,ITG,2*maxlenmpc);
  
  /* allocating fields for the actual external loading */
  
  NNEW(xbounact,double,*nboun);
  NNEW(xbounini,double,*nboun);
  for(k=0;k<*nboun;++k){
    xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){
    xbodyact[k]=xbody[k];}
  
  /* assigning the body forces to the elements */ 
  
  if(*nbody>0){

    /* check whether there the previous step was in the relative
       system and a change to the absolute system was requested */
      
    if((*nmethod==4)&&(alpha[1]>1.)){
      NNEW(itreated,ITG,*nk);
      FORTRAN(velinireltoabs,(ibody,xbody,cbody,nbody,set,
			      istartset,iendset,ialset,nset,veold,mi,
			      ipkon,kon,lakon,co,itreated));
      SFREE(itreated);
    }
    if(*inewton==1){NNEW(cgr,double,4**ne);}
  }
  
  /* for mechanical calculations: updating boundary conditions
     calculated in a previous thermal step */
  
  if(*ithermal<2) FORTRAN(gasmechbc,(vold,nload,sideload,
				     nelemload,xload,mi));
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
    NNEW(itg,ITG,*nload+3**nflow);
    NNEW(ieg,ITG,*nflow);
    /* max 6 triangles per face, 4 entries per triangle */
    NNEW(kontri,ITG,24**nload);
    NNEW(nloadtr,ITG,*nload);
    NNEW(nacteq,ITG,4**nk);
    NNEW(nactdog,ITG,4**nk);
    NNEW(v,double,mt**nk);
    FORTRAN(envtemp,(itg,ieg,&ntg,&ntr,sideload,nelemload,
		     ipkon,kon,lakon,ielmat,ne,nload,
		     kontri,&ntri,nloadtr,nflow,ndirboun,nactdog,
		     nodeboun,nacteq,nboun,ielprop,prop,&nteq,
		     v,network,physcon,shcon,ntmat_,co,
		     vold,set,nshcon,rhcon,nrhcon,mi,nmpc,nodempc,
		     ipompc,labmpc,ikboun,&nasym,ttime,&time,
		     iaxial));
    SFREE(v);
      
    if((*mcs>0)&&(ntr>0)){
      NNEW(inocs,ITG,*nk);
      radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
	     iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
    }
    else{ntrit=ntri;}
      
    nzsrad=100*ntr;
    NNEW(mast1rad,ITG,nzsrad);
    NNEW(irowrad,ITG,nzsrad);
    NNEW(icolrad,ITG,ntr);
    NNEW(jqrad,ITG,ntr+1);
    NNEW(ipointerrad,ITG,ntr);
      
    if(ntr>0){
      mastructrad(&ntr,nloadtr,sideload,ipointerrad,
		  &mast1rad,&irowrad,&nzsrad,
		  jqrad,icolrad);
    }
      
    /* determine the network elements belonging to a given node (for usage
       in user subroutine film */

    if((*network>0)||(ntg>0)){
      NNEW(iponoeln,ITG,*nk);
      NNEW(inoeln,ITG,2**nkon);
      if(*network>0){
	FORTRAN(networkelementpernode,(iponoeln,inoeln,lakon,ipkon,kon,
				       &inoelnsize,nflow,ieg,ne,network));
	FORTRAN(checkforhomnet,(ieg,nflow,lakon,ipkon,kon,itg,&ntg,
				iponoeln,inoeln));
      }
      RENEW(inoeln,ITG,2*inoelnsize);
    }

    SFREE(ipointerrad);SFREE(mast1rad);
    RENEW(irowrad,ITG,nzsrad);
      
    RENEW(itg,ITG,ntg);
    NNEW(ineighe,ITG,ntg);
    RENEW(kontri,ITG,4*ntrit);
    RENEW(nloadtr,ITG,ntr);
      
    NNEW(adview,double,ntr);
    NNEW(auview,double,2*nzsrad);
    NNEW(tarea,double,ntr);
    NNEW(tenv,double,ntr);
    NNEW(fenv,double,ntr);
    NNEW(erad,double,ntr);
      
    NNEW(ac,double,nteq*nteq);
    NNEW(bc,double,nteq);
    NNEW(ipiv,ITG,nteq);
    NNEW(adrad,double,ntr);
    NNEW(aurad,double,2*nzsrad);
    NNEW(bcr,double,ntr);
    NNEW(ipivr,ITG,ntr);
  }
  
  /* check for fluid elements
     check for strain-less elements */
  
  NNEW(nactdoh,ITG,*ne);
  NNEW(nactdohinv,ITG,*ne);
  *nef=0;
  for(i=0;i<*ne;++i){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")==0){
      icfd=1;nactdohinv[*nef]=i+1;++*nef;nactdoh[i]=*nef;}
    if(istrainfree==0){
      if(ielmat[i]<0){istrainfree=1;}
    }
  }

  if(icfd==1){
    if(iturbulent>=20){

      /* CBS method for shallow water equations */
      
      iturbulent=iturbulent-20;
      ifreesurface=1;
      icfd=2;
    }else if(iturbulent>=10){

      /* CBS method for all other applications */
      
      iturbulent=iturbulent-10;
      icfd=2;
    }
  }
  
  if(icfd==1){
  }else if(icfd==2){
      SFREE(nactdoh);SFREE(nactdohinv);

      /* rearranging the fluid nodes and elements such that
	 no gaps occur */

      NNEW(ipkonf,ITG,*nef);
      NNEW(lakonf,char,8**nef);
      NNEW(ielmatf,ITG,mi[2]**nef);
      if(*norien>0) NNEW(ielorienf,ITG,mi[2]**nef);
      NNEW(nelold,ITG,*nef);
      NNEW(nelnew,ITG,*ne);
      NNEW(cof,double,3**nk);
      NNEW(voldf,double,mt**nk);
      NNEW(nkold,ITG,*nk);
      NNEW(nknew,ITG,*nk);
      NNEW(inotrf,ITG,2**nk);
      NNEW(konf,ITG,*nkon);
      NNEW(ipompcf,ITG,*nmpc);
      NNEW(ikmpcf,ITG,*nmpc);
      NNEW(ilmpcf,ITG,*nmpc);
      NNEW(nodempcf,ITG,3*memmpc_);
      NNEW(coefmpcf,double,memmpc_);
      NNEW(nodebounf,ITG,*nboun);
      NNEW(ndirbounf,ITG,*nboun);
      NNEW(ikbounf,ITG,*nboun);
      NNEW(ilbounf,ITG,*nboun);
      if(*nam>0) NNEW(iambounf,ITG,*nboun);
      NNEW(xbounf,double,*nboun);
      NNEW(xbounoldf,double,*nboun);
      NNEW(xbounactf,double,*nboun);
      NNEW(nelemloadf,ITG,2**nload);
      if(*nam>0) NNEW(iamloadf,ITG,2**nload);
      NNEW(xloadf,double,2**nload);
      NNEW(xloadoldf,double,2**nload);
      NNEW(xloadactf,double,2**nload);
      NNEW(sideloadf,char,20**nload);
      if(*nbody>0) NNEW(ipobodyf,ITG,2*(*ifreebody-1));

      FORTRAN(rearrangecfd,(ne,ipkon,lakon,ielmat,ielorien,norien,nef,ipkonf,
			    lakonf,ielmatf,ielorienf,mi,nelold,nelnew,nkold,
			    nknew,nk,&nkf,konf,&nkonf,nmpc,ipompc,nodempc,
			    coefmpc,&memmpc_,&nmpcf,ipompcf,nodempcf,coefmpcf,
			    &memmpcf,nboun,nodeboun,ndirboun,xboun,&nbounf,
			    nodebounf,ndirbounf,xbounf,nload,nelemload,
			    sideload,xload,&nloadf,nelemloadf,sideloadf,xloadf,
			    ipobody,ipobodyf,kon,&nkftot,co,cof,vold,voldf,
			    ikbounf,ilbounf,ikmpcf,ilmpcf,iambounf,iamloadf,
			    iamboun,iamload,xbounold,xbounoldf,xbounact,
			    xbounactf,xloadold,xloadoldf,xloadact,xloadactf,
			    inotr,inotrf,nam,ntrans,nbody));

      /* call rearrangecfd */
      
      RENEW(nkold,ITG,nkftot);
      RENEW(cof,double,3*nkftot);
      RENEW(voldf,double,mt*nkftot);
      NNEW(inotrf,ITG,2*nkftot);
      RENEW(konf,ITG,nkonf);
      RENEW(ipompcf,ITG,nmpcf);
      RENEW(ikmpcf,ITG,nmpcf);
      RENEW(ilmpcf,ITG,nmpcf);
      RENEW(nodempcf,ITG,3*memmpcf);
      RENEW(coefmpcf,double,memmpcf);
      RENEW(nodebounf,ITG,nbounf);
      RENEW(ndirbounf,ITG,nbounf);
      RENEW(ikbounf,ITG,nbounf);
      RENEW(ilbounf,ITG,nbounf);
      if(*nam>0) RENEW(iambounf,ITG,nbounf);
      RENEW(xbounf,double,nbounf);
      RENEW(xbounoldf,double,nbounf);
      RENEW(xbounactf,double,nbounf);
      RENEW(nelemloadf,ITG,2*nloadf);
      if(*nam>0) NNEW(iamloadf,ITG,2*nloadf);
      RENEW(xloadf,double,2*nloadf);
      RENEW(xloadoldf,double,2*nloadf);
      RENEW(xloadactf,double,2*nloadf);
      RENEW(sideloadf,char,20*nloadf);
      
      /* calculating topological properties for CFD */
      
      NNEW(sideface,char,6**nef);
      NNEW(nelemface,ITG,6**nef);
      NNEW(ipface,ITG,*nef);
      NNEW(ipoface,ITG,nkf);
      NNEW(nodface,ITG,5*6**nef);
      NNEW(ifreestream,ITG,nkf);
      NNEW(isolidsurf,ITG,nkf);
      NNEW(neighsolidsurf,ITG,nkf);
      NNEW(iponoelf,ITG,nkf);
      NNEW(inoelf,ITG,2*8**nef);
      NNEW(inomat,ITG,nkftot);
      FORTRAN(topocfdfem,(nelemface,sideface,&nface,ipoface,nodface,nef,ipkonf,
			  konf,lakonf,&nkf,isolidsurf,&nsolidsurf,ifreestream,
			  &nfreestream,neighsolidsurf,iponoelf,inoelf,
			  &inoelfree,cof,set,istartset,iendset,ialset,nset,
			  &iturbulent,inomat,ielmatf,ipface,nknew));
      RENEW(sideface,char,nface);
      RENEW(nelemface,ITG,nface);
      SFREE(ipoface);SFREE(nodface);
      RENEW(ifreestream,ITG,nfreestream);
      RENEW(isolidsurf,ITG,nsolidsurf);
      RENEW(neighsolidsurf,ITG,nsolidsurf);
      RENEW(inoelf,ITG,2*inoelfree);
      if(*ithermal==1){
	NNEW(qfx,double,3*mi[0]**ne);}
  }else{
    SFREE(nactdoh);SFREE(nactdohinv);
  }
  
  if(*ithermal>1){
    NNEW(qfx,double,3*mi[0]**ne);}
  
  /* contact conditions */
  
  inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	  lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
          &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
          mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,ifacecount,
	  iperturb,ikboun,nboun,co,istep,&xnoels);
  
  if(ncont!=0){
      
    NNEW(cg,double,3*ncont);
    NNEW(straight,double,16*ncont);
	  
    /* 11 instead of 10: last position is reserved for the
       local contact spring element number; needed as
       pointer into springarea */
      
    if(*mortar<=0){
      RENEW(kon,ITG,*nkon+11**nslavs);
      NNEW(springarea,double,2**nslavs);
      if((*nener==1)&&((maxprevcontel==0)&&(*nslavs!=0))){
	RENEW(ener,double,2*mi[0]*(*ne+*nslavs));
	DOUMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nslavs),0.);

	/* setting the entries for the friction contact energy to zero */

	/*	for(k=mi[0]*(2**ne+*nslavs);k<mi[0]*(*ne+*nslavs)*2;k++){
		ener[k]=0.;}*/
      }
      RENEW(ipkon,ITG,*ne+*nslavs);
      RENEW(lakon,char,8*(*ne+*nslavs));
	  
      if(*norien>0){
	RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
	for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++){
	  ielorien[k]=0;}
      }

      RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
      for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++){
	ielmat[k]=1;}

      if((maxprevcontel==0)&&(*nslavs!=0)){
	RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
	  xstate[k]=0.;
	}
      }
      maxprevcontel=*nslavs;

      NNEW(areaslav,double,*ifacecount);
    }else if(*mortar==1){
      NNEW(islavact,ITG,nslavnode[*ntie]);
      if((*istep==1)||(nslavs_prev_step==0))
	NNEW(clearini,double,3*9**ifacecount);

      /* check whether at least one contact definition involves true contact
	 and not just tied contact */

      FORTRAN(checktruecontact,(ntie,tieset,tietol,elcon,&itruecontact,
				ncmat_,ntmat_));
    }else if(*mortar>1){
      ismallsliding=1;
      NNEW(slavnor,double,3**nslavs);
      NNEW(slavtan,double,6**nslavs);
      inimortar(&ener,mi,ne,nslavs,nk,nener,&ipkon,&lakon,&kon,nkon,
		&maxprevcontel,&xstate,nstate_,&islavtie,&bp,&islavact,
		&gap,&cdisp,&cstress,&cfs,
		&bpini,&islavactini,&cstressini,ntie,
		tieset,nslavnode,islavnode,&islavnodeinv,&islavquadel,
		&pslavdual,&aut,&irowt,&jqt,&autinv,
		&irowtinv,&jqtinv,&Bd,&irowb,&jqb,&Bdhelp,&irowbhelp,
		&jqbhelp,&Dd,&irowd,&jqd,&Ddtil,&irowdtil,&jqdtil,&Bdtil,
		&irowbtil,&jqbtil,itiefac,islavsurf,nboun,
		nmpc,&nslavspc,&islavspc,&nslavmpc,&islavmpc,
		&nmastspc,&imastspc,&nmastmpc,&imastmpc,
		imastnode,nmastnode,&nasym,mortar,&ielmat,&ielorien,norien,
		ipompc,nodempc,ikboun,ilboun,ikmpc,ilmpc,jobnamef,set,
		co,vold,nset,&nslavquadel);
    }
    NNEW(xmastnor,double,3*nmastnode[*ntie]);
  }
  
  if(icascade==2){
    memmpcref_=memmpc_;mpcfreeref=mpcfree;maxlenmpcref=maxlenmpc;
    NNEW(nodempcref,ITG,3*memmpc_);
    for(k=0;k<3*memmpc_;k++){
      nodempcref[k]=nodempc[k];}
    NNEW(coefmpcref,double,memmpc_);
    for(k=0;k<memmpc_;k++){
      coefmpcref[k]=coefmpc[k];}
  }
  
  if((*ithermal==1)||(*ithermal>=3)){
    NNEW(t1ini,double,*nk);
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){
      t1act[k]=t1old[k];}
  }
  
  /* allocating a field for the instantaneous amplitude */
  
  NNEW(ampli,double,*nam);
  
  /* fini is also needed in static calculations if iforbou=1
     to get correct values of f after a divergent increment */

  NNEW(fini,double,neq[1]);
  
  /* allocating fields for nonlinear dynamics */
  
  if(*nmethod==4){
    mass[0]=1;
    mass[1]=1;
    NNEW(aux2,double,neq[1]);
    NNEW(fextini,double,neq[1]);
    NNEW(fnext,double,mt**nk);
    NNEW(fnextini,double,mt**nk);
    NNEW(veini,double,mt**nk);
    NNEW(accini,double,mt**nk);
    NNEW(adb,double,neq[1]);
    NNEW(aub,double,nzs[1]);
    NNEW(cvini,double,neq[1]);
    NNEW(cv,double,neq[1]);
  }

  //    if((*nstate_!=0)&&((*mortar!=1)||(ncont==0))){
  if((*nstate_!=0)&&(*mortar!=1)){
    NNEW(xstateini,double,*nstate_*mi[0]*(*ne+*nslavs));
    isiz=*nstate_*mi[0]*(*ne+*nslavs);cpypardou(xstateini,xstate,&isiz,&num_cpus);
    //    FORTRAN(stop,());
  }

  /* next lines: change on 8th of July 2023: initial state values
     for dynamic plastic calculations */
  
  if((*nstate_!=0)&&(*mortar==1)){
    NNEW(xstateini,double,*nstate_*mi[0]**ne);
    isiz=*nstate_*mi[0]**ne;cpypardou(xstateini,xstate,&isiz,&num_cpus);
  }
  
  NNEW(eei,double,6*mi[0]**ne);
  NNEW(stiini,double,6*mi[0]**ne);
  NNEW(emeini,double,6*mi[0]**ne);
  
  if(*nener==1){
    if((*mortar!=1)||(ncont==0)){
      NNEW(enerini,double,2*mi[0]*(*ne+*nslavs));
    }else{
      NNEW(enerini,double,2*mi[0]*(*ne+*nintpoint));
    }
    isiz=2*mi[0]**ne;cpypardou(enerini,ener,&isiz,&num_cpus);
  }
  
  qa[0]=qaold[0];
  qa[1]=qaold[1];
  
  /* normalizing the time */
  
  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,inext,&itp,istep,tper));
  dtheta=(*tinc)/(*tper);

  /* taking care of a small increment at the end of the step
     for face-to-face penalty contact */

  dthetaref=dtheta;
  if((dtheta<=1.e-6)&&(*iexpl<=1)){
    printf("\n *ERROR in nonlingeo\n");
    printf(" increment size smaller than one millionth of step size\n");
    printf(" increase increment size\n\n");
  }
  /* [SWITCHES] State the configuration of this run before anything acts on
     it, and name any CCX_* variable that is set but not read.  A flag that
     was never read is the difference between an A/B that is wrong and one
     that is uninformative, and the second kind costs a whole run to notice.
     Reports only; reads nothing, decides nothing. */
  ccxopt_report();

  opcheckdrv_configure_fd(&opd);

  /* [DAMAGE TMIN] statics.f:234-247 silently raises the deck's minimum
     increment to min(tinc,1e-6*tper) under automatic incrementation.  With
     tinc=1e-3 and tper=1 every deck in this project has therefore been run
     with tmin=1e-6, not the 1e-9 it asks for, and "increment size smaller
     than minimum" has been a CalculiX floor rather than a physical limit.
     Opt-in, physical units, applied before the normalisation below. */

  if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_TMIN"))!=NULL){
    double tmnew=atof(damage_de13_env);
    if(tmnew>0.){
      printf("[DAMAGE TMIN] minimum increment overridden: %.12e -> %.12e "
             "(physical units).  The statics.f 1e-6*tper floor is bypassed; "
             "the stock cutback machinery is otherwise untouched.%s",
             *tmin,tmnew,"\n");
      fflush(stdout);
      *tmin=tmnew;
    }
  }

  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;
  
  /* calculating an initial flux norm */
  
  if(*ithermal!=2){
    if(qau>1.e-10){qam[0]=qau;}
    else if(qa0>1.e-10){qam[0]=qa0;}
    else if(qa[0]>1.e-10){qam[0]=qa[0];}
    else {qam[0]=1.e-2;}
  }
  if(*ithermal>1){
    if(qau>1.e-10){
      qam[1]=qau;}
    else if(qa0>1.e-10){
      qam[1]=qa0;}
    else if(qa[1]>1.e-10){
      qam[1]=qa[1];}
    else {qam[1]=1.e-2;}
  }
  
  /* storing the element and topology information before introducing 
     contact elements */
  
  ne0=*ne;nkon0=*nkon;neold=*ne;

  /* DE1.2/DM2.0 progressive damage is enabled for either the historical
     four-constant Rice-Tracey evolution record or the new type-3 tabulated
     ductile fracture-locus record.  Both use the same Newton-integrated
     evolution and DE1.3.1 terminal topology backend. */
  if((*ndmat_>0)&&(*iexpl<=1)){
    for(i=0;i<*nmat;i++){
      if(damage_progressive_material(i+1,ndmcon,dmcon,*ndmat_,*ntmat_)){
        damage_de12_matcount++;
        if((ITG)dmcon[1+(*ndmat_+1)*(*ntmat_)*i]==3)
          damage_dm20_matcount++;
      }
    }
    if(damage_de12_matcount>0) damage_de12_enabled=1;
    if(damage_de12_enabled){
      /* FD_SYM - mode 1, the symmetric part of the rank-1 term folded
         into the 21-entry tangent - was DELETED.  It had
         two convention defects that made it destroy runs; corrected, it
         cost 13% more Newton iterations than applying no correction at
         all on both fast decks and moved the fracture by one element.
         The measurement carries the table.  The spelling
         is still recognised so that a deck or a script carrying it is
         TOLD, rather than silently running the stock tangent. */

      damage_tangent_env=ccxopt_getenv("CCX_DAMAGE_TANGENT");
      if((damage_tangent_env!=NULL)&&
         ((strcmp(damage_tangent_env,"FD_SYM")==0)||
          (strcmp(damage_tangent_env,"fd_sym")==0)||
          (strcmp(damage_tangent_env,"1")==0))){
        printf("*ERROR: CCX_DAMAGE_TANGENT=FD_SYM was removed.  The "
               "symmetric part of a nonsymmetric rank-1 term is not a "
               "Newton tangent; measured, it cost 13%% more iterations "
               "than no correction at all.  Use UNSYM for the consistent "
               "tangent, or leave the switch unset for g(D)*Cep.\n");
        FORTRAN(stop,());
      }else if((damage_tangent_env!=NULL)&&
               ((strcmp(damage_tangent_env,"UNSYM")==0)||
                (strcmp(damage_tangent_env,"unsym")==0)||
                (strcmp(damage_tangent_env,"2")==0))){

        /* Stage 1 of the consistent-tangent work: assemble and solve the
           bulk problem through the asymmetric path with the constitutive
           tangent left untouched.  resultsmech.f only ever tests
           de12tangent.eq.1, so mode 2 changes no material code. */

        damage_tangent_mode=2;
      }
      rsc.reeq_scale_env=ccxopt_getenv("CCX_DAMAGE_REEQ_SCALE");
      if((rsc.reeq_scale_env!=NULL)&&
         ((strcmp(rsc.reeq_scale_env,"PHYSICAL")==0)||
          (strcmp(rsc.reeq_scale_env,"physical")==0)||
          (strcmp(rsc.reeq_scale_env,"1")==0))){
        rsc.reeq_scale_mode=1;
      }
      /* J-09.  The same-load re-equilibration after a terminal deletion
         is tested with the STOCK RELATIVE correction criterion cam/uam, and
         REEQ_SCALE=PHYSICAL floors uam at the converged physical increment's
         displacement norm.  That floor SHRINKS WITH THE STEP, while the
         correction it is tested against does not: the perturbation is an
         element vanishing, not a load increment.  Measured on
         run_m12s3_rad - three consecutive failed attempts at dt = 2.64e-5,
         6.59e-6, 1.65e-6 (16x apart) give residuals agreeing to three
         digits, 0.548 / 0.550 / 0.551, while actual_uam holds at 1.43e-3
         and physical_ref collapses 2.37e-3 -> 1.48e-4.  Cutting the step is
         therefore the wrong instrument, and the loop can only end at tmin.

         This holds the reference at a fraction of the largest value it
         reached, exactly as CCX_DAMAGE_QAM_FLOOR does for the force side.
         The force side already had this; the displacement side did not.

         THIS CHANGES THE CONVERGENCE CRITERION AND THEREFORE THE ANSWER.
         A run that goes further with it is NOT thereby a success - that has
         to be shown on the physics (hydride consumption, failed facets,
         peakedness), and adopting it needs the full verify + ladder gate.
         Default 0 = off = bit-identical to the unpatched binary. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_REEQ_FLOOR"))!=NULL){
        rsc.reeq_uam_floor=atof(damage_de13_env);
        if(rsc.reeq_uam_floor<0.) rsc.reeq_uam_floor=0.;
        if(rsc.reeq_uam_floor>1.) rsc.reeq_uam_floor=1.;
      }

      probedrv_configure_aba(&prb);

      rescue_configure_backtrack(&rsc);

      rescue_configure_levels(&rsc,&lc,&prb);

      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_RESIDUAL_RAY"))!=NULL){
        prb.ray_probe=1;
        prb.ray_max=atoi(damage_de13_env);
        if(prb.ray_max<1) prb.ray_max=8;
        if(prb.ray_max>1000) prb.ray_max=1000;
        /* CCX_DAMAGE_RAY_INC="231" - walk only at these increments. */
        prb.ray_ninc=0;
        if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_RAY_INC"))!=NULL){
          char *rcp=damage_de13_env;
          while((*rcp!=0)&&(prb.ray_ninc<4)){
            while((*rcp==' ')||(*rcp==',')) rcp++;
            if(*rcp==0) break;
            prb.ray_inc[prb.ray_ninc++]=atoi(rcp);
            while((*rcp!=0)&&(*rcp!=',')) rcp++;
          }
        }
        printf("[DAMAGE RAY] DIAGNOSTIC: the residual is scanned along the "
               "Newton direction during same-load re-equilibration, at most "
               "%" ITGFORMAT " times.  Reads only; b is restored exactly and "
               "the scan ends on the full step, so no bit of the answer "
               "changes.%s",prb.ray_max,"\n");
        if(prb.ray_ninc>0){
          printf("[DAMAGE RAY] restricted to increments:");
          for(i=0;i<prb.ray_ninc;i++)
            printf(" %" ITGFORMAT,prb.ray_inc[i]);
          printf("%s","\n");
        }
        fflush(stdout);
      }

      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_RELEASE_PROBE"))!=NULL){
        damage_release_probe=atoi(damage_de13_env);
        if(damage_release_probe<0) damage_release_probe=0;
        if(damage_release_probe>2) damage_release_probe=2;
        if(damage_release_probe>0){
          printf("[DAMAGE RELEASE] DIAGNOSTIC level %" ITGFORMAT
                 ": the internal force released by each topology event is "
                 "measured at frozen displacement and reported split into "
                 "surviving and removed degrees of freedom.  Reads only; "
                 "changes no bit of the answer.%s", damage_release_probe,
                 "\n");
          fflush(stdout);
        }
      }
      rsc.linesearch_env=ccxopt_getenv("CCX_DAMAGE_LINESEARCH");
      if((rsc.linesearch_env!=NULL)&&
         ((strcmp(rsc.linesearch_env,"ADAPTIVE")==0)||
          (strcmp(rsc.linesearch_env,"adaptive")==0)||
          (strcmp(rsc.linesearch_env,"1")==0))){
        rsc.linesearch_mode=1;
      }
      /* DE1.3 keeps an element until it is degraded to 0.1% of its
         stiffness.  Between D=0.95 and that point the element carries
         almost no load but still owns equations, and in a phase that is
         failing over a whole region there can be thousands of them at
         once.  Making the threshold settable lets that be measured
         instead of argued about. */

      if((damage_de13_env=ccxopt_getenv("CCX_FRACTURE_LINK"))!=NULL){
        if((strcmp(damage_de13_env,"FACE")==0)||
           (strcmp(damage_de13_env,"face")==0)){
          damage_fracture_link=1;
        }else if((strcmp(damage_de13_env,"NODE")!=0)&&
                 (strcmp(damage_de13_env,"node")!=0)){
          printf("[FRACTURE TERMINATION] *WARNING: CCX_FRACTURE_LINK=%s is not NODE or FACE; keeping NODE\n",damage_de13_env);
        }
      }
      dogleg_configure(&dog,&rsc,&lc);

      damcont_configure(&ct,&dog,&rsc,&lc);

      /* A requested rescue flag that silently fails to arm has already
         cost two whole runs: CCX_DAMAGE_REEQ_RESCUE3 and then
         CCX_DAMAGE_RESCUE_CORRIDOR were each parsed inside a block whose
         outer condition did not list them, so the run executed as the
         plain control and read like a negative result.  Shout instead of
         pretending. */
      if((ccxopt_getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)&&
         (rsc.corr_mode==0)){
        printf("*ERROR: CCX_DAMAGE_RESCUE_CORRIDOR is set but the "
               "corridor did NOT arm; the run would silently be a plain "
               "control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)&&
         (lc.reg_nlam==0)){
        printf("*ERROR: CCX_DAMAGE_REEQ_RESCUE3 is set but the "
               "regularized level did NOT arm.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((ccxopt_getenv("CCX_DAMAGE_CONTINUATION")!=NULL)&&
         ((ct.mode==0)||(rsc.rescue_maxlevel!=4))){
        printf("*ERROR: CCX_DAMAGE_CONTINUATION is set but the continuation "
               "level did NOT arm; the run would silently be a plain "
               "Rescue2+dogleg control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((ccxopt_getenv("CCX_DAMAGE_TR_DOGLEG")!=NULL)&&
         ((dog.mode==0)||(rsc.rescue_maxlevel<3))){
        printf("*ERROR: CCX_DAMAGE_TR_DOGLEG is set but the dogleg level "
               "did NOT arm; the run would silently be a plain Rescue2 "
               "control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if(((ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE")!=NULL)||
          (ccxopt_getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL))&&
         (rsc.rescue_mode==0)){
        printf("*ERROR: a rescue flag is set but rescue did NOT arm.  "
               "Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }

      /* CCX_FRACTURE_CUT: judge the load path by AREA rather than by
         existence.  damconnect.f answers whether a chain of surviving
         elements still links the grips and says in its own comment why it
         can never answer more: "a nonzero stiffness is a load path, and
         calling it broken would be a modelling decision rather than a
         topological fact".  The modelling decision is exactly what is
         missing, and this is where the deck gets to make it - as a
         FRACTION of the load path the specimen started with, which is
         dimensionless and needs no knowledge of the section area.

         Measured on s3rad: the run ends with
         every topological rule saying CONNECTED and the two grips joined
         by ONE triangular face - 0.05% of a section - carrying 2.1% of
         peak load, after spending 18.4% of its wall clock on a specimen
         that was already finished. */

      damage_cut_frac=0.;
      damage_cut_ref=-1.;
      damage_cut_on=0;
      damage_cut_narrow=0;
      damage_cut_exact=0;
      {
        const char *cutenv=ccxopt_getenv("CCX_FRACTURE_CUT");
        if(cutenv!=NULL){
          damage_cut_frac=atof(cutenv);
          if((damage_cut_frac<=0.)||(damage_cut_frac>=1.)){
            printf("[LOADCUT] *ERROR: CCX_FRACTURE_CUT must be a fraction "
                   "strictly between 0 and 1; got %s.  Not armed.\n",cutenv);
            damage_cut_frac=0.;
          }else{
            damage_cut_on=1;

            /* The residual-stiffness floor is a DECK number, and it was
               hard coded here as 1e-4 while resultsmech.f and
               mafilldamas.f read CCX_DAMAGE_GMIN.  Three copies of one
               constant, two of them settable and one not, is how a
               measurement of the floor's own effect would have come out
               wrong.  Parsed here exactly as they parse it, bounds
               included. */

            {
              const char *ge=ccxopt_getenv("CCX_DAMAGE_GMIN");
              damage_cut_gmin=1.e-4;
              if(ge!=NULL){
                double gv=atof(ge);
                if((gv>=1.e-6)&&(gv<=0.2)) damage_cut_gmin=gv;
              }
              printf("[LOADCUT] residual-stiffness floor in the weighting: "
                     "%.6e\n",damage_cut_gmin);
            }
            if(ccxopt_getenv("CCX_FRACTURE_CUT_EXACT")!=NULL){
              damage_cut_exact=1;
              printf("[LOADCUT] the early exit is DISABLED: every batch "
                     "reports the true cut rather than a lower bound.  "
                     "Costs a full max-flow per committed batch and is "
                     "how the trajectory is measured before a threshold "
                     "is chosen.\n");
            }
            printf("[LOADCUT] armed: the run stops when the minimum cut "
                   "between the termination sets falls below %.4g of its "
                   "value at the first committed batch.\n",damage_cut_frac);
            printf("[LOADCUT] the cut is the smallest total (face area x "
                   "residual stiffness) that would have to break to "
                   "separate them - a width, not a yes/no.\n");
            fflush(stdout);
          }
        }
      }

      damage_fracture_env=ccxopt_getenv("CCX_FRACTURE_TERMINATION");
      if(damage_fracture_env!=NULL){
        damage_fracture_seta=strdup(damage_fracture_env);
        damage_fracture_setb=strchr(damage_fracture_seta,':');
        if(damage_fracture_setb!=NULL){
          *damage_fracture_setb=0;
          damage_fracture_setb++;
        }else{
          free(damage_fracture_seta);
          damage_fracture_seta=NULL;
          printf("[FRACTURE TERMINATION] *WARNING: expected "
                 "CCX_FRACTURE_TERMINATION=SETA:SETB; ignored\n");
        }
        if(damage_fracture_seta!=NULL){

          /* blank padded to 81 for the Fortran side: the FORTRAN macro
             passes no hidden string length, so an assumed-length dummy
             argument there would read garbage */

          memset(damage_fracture_a,' ',81);
          memset(damage_fracture_b,' ',81);
          memcpy(damage_fracture_a,damage_fracture_seta,
                 (strlen(damage_fracture_seta)<80)?
                 strlen(damage_fracture_seta):80);
          memcpy(damage_fracture_b,damage_fracture_setb,
                 (strlen(damage_fracture_setb)<80)?
                 strlen(damage_fracture_setb):80);
          printf("[FRACTURE TERMINATION] the run stops as soon as no "
                 "surviving element links %s to %s, elements conducting "
                 "through a shared %s\n",
                 damage_fracture_seta,damage_fracture_setb,
                 damage_fracture_link?"FACE":"node");
        }
      }
      /* Exclude a fully debonded cohesive facet from the TERMINATION
         connectivity.  A UC6 facet is never deleted - cohesive_uc6.f pins
         g at gmin - so once it has failed it still reads as a load path for
         ever, and on an interface-dominated fracture [FRACTURE COMPLETE] can
         then never fire.  Measured on DHC1 with facet viscosity: the specimen
         carries 0.0% of peak with 206 of 416 facets fully failed, and an
         offline replay puts severance at t=0.6100 while the run was driven on
         to 0.7124 (E-90).

         Failure is read from xstate slot 4, which resultsmech_uc6.f already
         writes as an explicit flag and which nothing else reads.  Like
         CCX_FRACTURE_LINK this deletes nothing and changes no equation - only
         the moment the run may stop - so it carries none of the risk that sank
         damfloatface (E-57, E-72).  Default OFF. */
      if((damage_de13_env=ccxopt_getenv("CCX_FRACTURE_DEADFACET"))!=NULL){
        damage_deadfacet=(strcmp(damage_de13_env,"0")==0)?0:1;
        if(damage_deadfacet){
          printf("[FRACTURE TERMINATION] a cohesive facet whose every "
                 "integration point has failed is excluded from the load "
                 "path\n");
        }
      }

      /* CCX_DAMAGE_FACET_DELETE - the same failure flag as
         CCX_FRACTURE_DEADFACET, but acted on in the EQUATIONS instead of only
         in the termination test.

         The bulk phases create a displacement discontinuity when they fail;
         the interface never does, because g = max(gmin, 1-dvisc) floors the
         facet stiffness and a fully debonded facet keeps conducting gmin*Kn
         for ever.  Measured directly: 45.3% of the radial twin's backed
         facets sit at Dmin >= 0.999 with ZERO removed (E-113).  E-111 item 1
         named this as the remaining gap once the internal length (E-114) and
         the penalty stiffness (E-113) had both been tested and excluded, and
         E-119/E-120 localised the unresolved route to DEBONDING - the one
         route whose crack this prevents.

         Unlike CCX_FRACTURE_DEADFACET this DOES change the equations, so it
         carries the risk that sank damfloatface (E-57, E-72) and is gated the
         same way: default OFF, batch-capped, and travelling the unchanged
         transactional path so rollback restores it (E-64).  It must not be
         adopted without verify 65/65 and the full ladder. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_FACET_DELETE"))!=NULL){
        damage_facetdel=(strcmp(damage_de13_env,"0")==0)?0:1;
        if(damage_facetdel){
          printf("[DAMAGE FACET DELETE] a cohesive facet whose every "
                 "integration point has failed is REMOVED from the equations "
                 "(its compressive penalty is removed with it)\n");
        }
      }

      if(ccxopt_getenv("CCX_DISSIPATION_REPORT")!=NULL) lc.diss_report=1;
      lc.diss_env=ccxopt_getenv("CCX_DISSIPATION_TARGET");
      if(lc.diss_env!=NULL){
        lc.diss_target=atof(lc.diss_env);
        if(lc.diss_target>0.) lc.diss_report=1;
      }
      if(ccxopt_getenv("CCX_DISSIPATION_PROBE")!=NULL) lc.diss_probe=1;
      if(ccxopt_getenv("CCX_DAMAGE_BATCH_LIST")!=NULL) damage_batch_list=1;
      /* The stress is scaled by 1-Dvis whenever the viscosity is
         on, so the deletion trigger has to read the same variable.
         Measured on DHC1: elements were being removed while still
         carrying up to 41% of their effective stress (batch_Dvis
         down to 0.592), and correcting it moved the run from
         lambda=0.3793 to 0.4651.  With the viscosity off damvisc
         is never allocated and the trigger falls back to D, so
         this is a no-op there - verified on SP1. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_DELETE_VISC"))!=NULL){
        damage_epol.delete_visc=(strcmp(damage_de13_env,"0")==0)?0:1;
      }
      if(ccxopt_getenv("CCX_DAMAGE_FREE_PROBE")!=NULL) damage_free_probe=1;
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NODE_DUMP"))!=NULL)
        prb.dump_node=atoi(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NODE_INC"))!=NULL)
        prb.dump_inc=atoi(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NULLVEC"))!=NULL)
        prb.null_inc=atoi(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NULLVEC_IT"))!=NULL)
        prb.null_nit=atoi(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_PATH_DROP"))!=NULL)
        lc.path_drop=atof(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_PATH_NSTEP"))!=NULL)
        lc.path_nstep=atoi(damage_de13_env);
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_PATH"))!=NULL){
        lc.path_on=atoi(damage_de13_env);
        if(lc.path_on<0) lc.path_on=0;
        if(lc.path_on>2) lc.path_on=2;
      }
      if(ccxopt_getenv("CCX_DAMAGE_STIFF_PROBE")!=NULL)
        damage_stiff_probe=1;
      if(ccxopt_getenv("CCX_DAMAGE_TANGENT_CENSUS")!=NULL)
        opd.unsym_census=1;
      /* resultsmech.f reads this name too and has no safe place to print
         from - it runs on several threads.  Reporting it here is what makes
         the difference between "the flag ran and changed nothing" and "the
         flag was never read", which is exactly the ambiguity that made the
         first CCX_DAMAGE_TANGENT_FULL A/B uninformative (J-10). */
      if(ccxopt_getenv("CCX_DAMAGE_TANGENT_FULL")!=NULL){
        opd.unsym_tanfull=1;
        printf("[DAMAGE TANGENT FULL] the consistent-tangent cut-off is "
               "taken on the SAME damage variable the stress uses instead "
               "of the stock D<0.999.  This CHANGES THE OPERATOR and "
               "therefore the answer\n");
        fflush(stdout);
      }
      /* Nonlocal damage: an internal length in the FORMULATION.
         Crack-band scales the dissipation and is correct, but it gives
         the localisation no width - three uniform bar meshes agree to
         0.94% before the peak and spread to 73-89% after it, and
         viscosity changes solvability without moving those curves at
         all (E-83).  Default off, so every stored baseline stands. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NONLOCAL"))!=NULL){
        damage_nl_ell=atof(damage_de13_env);
        if(damage_nl_ell<0.) damage_nl_ell=0.;
        if(damage_nl_ell>0.){
          printf("[DAMAGE NONLOCAL] the damage driving variable is averaged over a neighbourhood of radius 2*ell with ell=%.4e; the stress update stays local\n",damage_nl_ell);
        }
      }
      FORTRAN(damnonlocalset,(&damage_nl_ell));
      /* Which regularisation backend computes that average.

         INTEGRAL is the default and stays the default: it is what E-84
         and E-85 measured, so `CCX_DAMAGE_NONLOCAL=<ell>` alone must keep
         meaning exactly what those records say it means.

         GRADIENT solves  ebar - div(ell^2 grad ebar) = e  instead.  Same
         internal length, no neighbour list: the integral form stores
         O(N*(2*ell/h)^3) neighbours, which is 6.5 GB at ell=0.15 and
         58 GB at ell=0.45 on demo_realistic_clusters3 against 31.6 GB of
         machine, because one hydride element's neighbourhood holds
         thousands of tiny elements.  The PDE form is O(nnz) on any mesh
         and is the only one that can be used on a graded mesh at all. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_NONLOCAL_MODE"))!=NULL){
        if((strcmp(damage_de13_env,"GRADIENT")==0)||
           (strcmp(damage_de13_env,"gradient")==0)){
          damage_nl_mode=1;
        }else if((strcmp(damage_de13_env,"INTEGRAL")==0)||
                 (strcmp(damage_de13_env,"integral")==0)){
          damage_nl_mode=0;
        }else if((strcmp(damage_de13_env,"FROZEN")==0)||
                 (strcmp(damage_de13_env,"frozen")==0)){
          /* Not a backend: the CONTROL that separates the internal length
             from the integration scheme.  With ell>0 the driving variable
             is refreshed only at the predictor and the commit, so damage
             stops being integrated implicitly inside Newton - a change of
             scheme, not a lag.  Every local-vs-nonlocal comparison in this
             tree therefore moves two factors at once.  FROZEN supplies the
             same staggered update with NO averaging and NO length. */
          damage_nl_mode=2;
        }else{
          printf("*ERROR: CCX_DAMAGE_NONLOCAL_MODE must be INTEGRAL, GRADIENT or FROZEN, not %s\n",damage_de13_env);
          FORTRAN(stop,());
        }
      }
      if((damage_nl_ell>0.)&&(damage_nl_mode==1)){
        printf("[DAMAGE NONLOCAL] backend GRADIENT: implicit gradient, ebar - div(ell^2 grad ebar) = e, lumped mass, Jacobi-CG; no neighbour list\n");
      }
      if((damage_nl_ell>0.)&&(damage_nl_mode==2)){
        printf("[DAMAGE NONLOCAL] backend FROZEN: CONTROL ONLY - the staggered damage update WITHOUT spatial averaging.  ell is read but NOT applied; this arm regularises nothing and exists to separate the internal length from the integration scheme\n");
      }
      FORTRAN(damnonlocalmode,(&damage_nl_mode));
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_QAM_FLOOR"))!=NULL){
        damage_qam_floor=atof(damage_de13_env);
        if(damage_qam_floor<0.) damage_qam_floor=0.;
        if(damage_qam_floor>1.) damage_qam_floor=1.;
        if(damage_qam_floor>0.){
          printf("[DAMAGE QAM FLOOR] DIAGNOSTIC: the reference force qam is held at or above %.3e of its running maximum, so the RELATIVE force criterion cannot collapse as the specimen unloads.  This CHANGES THE CONVERGENCE CRITERION and therefore the answer; going further with it is not by itself a success\n",damage_qam_floor);
        }
      }
      if(ccxopt_getenv("CCX_DAMAGE_AUTOSPC_NEG")!=NULL) damage_spc_neg=1;
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_AUTOSPC"))!=NULL){
        damage_spc_g=atof(damage_de13_env);
        if(damage_spc_g<0.) damage_spc_g=0.;
        if(damage_spc_g>1.e-1) damage_spc_g=1.e-1;
        if(damage_spc_g>0.){
          damage_stiff_probe=1;
          printf("[DAMAGE AUTOSPC] a node whose assembled diagonal has fallen below %.1e of its own intact value is excluded from the DISPLACEMENT convergence norm; the force residual is untouched and nothing is deleted\n",damage_spc_g);
        }
      }
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_STIFF_MIN"))!=NULL){
        damage_stiff_min=atof(damage_de13_env);
        if(damage_stiff_min<0.) damage_stiff_min=0.;
        if(damage_stiff_min>0.1) damage_stiff_min=0.1;
        if(damage_stiff_min>0.) damage_stiff_probe=1;
      }
      /* [CENSUS] the same discipline, one level down: the census is a
         MEASUREMENT and not a judgement, but a measurement that is wrong is
         worse than none, because it is the number every reading in
         the census is taken against.  Prove it before
         printing it, on every run that arms it.

         Note what is suppressed on failure: the REPORT, and only the report.
         damage_stiff_probe also gates the assembly of damage_addiag, which
         the load-path judgement reads - clearing it would silently change
         the trajectory, which is exactly the kind of behaviour change a
         self test must not smuggle in.  Measured: doing it the wrong way
         round moved m.cvg on three gate cases. */
      /* Default 0 since E-22.  damdangle cannot distinguish a dangling
         sliver from the legitimate last element at a node, so it deletes
         healthy load-bearing material and on a small mesh removes every
         element, after which remastruct runs on an empty model and the
         process dies.  Measured cost: R2 and R3S PASS->FAIL,
         nc2_two_tet_wp and bk4_orphan_tet_wp access violation.  Measured
         benefit on DHC1: +0.0009, down from the +0.021 of E-04, which the
         E-05 deletion-on-Dtilde fix superseded.  Set to 1 only to
         reproduce the old behaviour. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_DANGLE"))!=NULL){
        damage_dangle_max=atoi(damage_de13_env);
        if(damage_dangle_max<0) damage_dangle_max=0;
        if(damage_dangle_max>4) damage_dangle_max=4;
      }
      damage_deadall_env=ccxopt_getenv("CCX_DAMAGE_DEADALL");
      if(damage_deadall_env!=NULL){
        damage_epol.deadall_g=atof(damage_deadall_env);
        if(damage_epol.deadall_g<0.) damage_epol.deadall_g=0.;
        if(damage_epol.deadall_g>0.5) damage_epol.deadall_g=0.5;
        if(damage_epol.deadall_g>0.){
          printf("[DAMAGE DEADALL] a node whose ENTIRE live support is dead "
                 "(g < %.1e) and which no cohesive facet holds has that "
                 "support deleted; only dead elements can be taken\n",
                 damage_epol.deadall_g);
        }
      }
      damage_deadsole_env=ccxopt_getenv("CCX_DAMAGE_DEADSOLE");
      if(damage_deadsole_env!=NULL){
        damage_epol.deadsole_g=atof(damage_deadsole_env);
        if(damage_epol.deadsole_g<0.) damage_epol.deadsole_g=0.;
        if(damage_epol.deadsole_g>0.5) damage_epol.deadsole_g=0.5;
        if(damage_epol.deadsole_g>0.){
          printf("[DAMAGE DEADSOLE] a dead element (g < %.1e) that is the sole "
                 "support of a node is deleted\n",damage_epol.deadsole_g);
        }
      }
      damage_stab_env=ccxopt_getenv("CCX_DAMAGE_STABILISE");
      if(damage_stab_env!=NULL){
        damage_stab_alpha=atof(damage_stab_env);
        if(damage_stab_alpha<0.) damage_stab_alpha=0.;
        if(damage_stab_alpha>1.e-1) damage_stab_alpha=1.e-1;
        if(damage_stab_alpha>0.){
          printf("[DAMAGE STABILISE] detached pieces held with alpha=%.3e "
                 "of the mean diagonal; nothing is deleted\n",
                 damage_stab_alpha);
        }
      }
      lc.diss_env=ccxopt_getenv("CCX_DISSIPATION_CONTROL");
      if((lc.diss_env!=NULL)&&(lc.diss_target>0.)){
        lc.diss_ctrl=(strcmp(lc.diss_env,"2")==0)?2:1;
        if(lc.diss_ctrl==2){
          NNEW(lc.diss_fhat,double,neq[1]);
          NNEW(lc.diss_uf,double,neq[1]);
        }
        lc.diss_report=1;
        printf("[DISSIPATION CONTROL] load factor solved from "
               "dG=%.6e per increment; theta is advanced only on "
               "acceptance\n",lc.diss_target);
      }

      /* Path following: give lambda an identity separate from theta.

         Everything in this file has so far identified the load factor WITH
         the step time: `lc.diss_lprev=theta`, and both dissipation
         schemes end by clamping `lamcur` back to `theta`
         (nonlingeo.c 5153 and 5635).  theta only ever advances, because
         dtime=dtheta*tper feeds the viscous update and must stay positive,
         so lambda could be pulled back INSIDE an increment and never end
         below the last converged value.  That is why the file's own comment
         says a genuine snap-back stays out of reach.

         The separation is the whole change:

             theta   monotone pseudo-time, still drives dtime
             lambda  the load factor, free to decrease

         and xbounact is built from lambda, which is legitimate because
         xbounold holds the value at the START of the step, so
         xbounold+(xboun-xbounold)*lambda is an absolute step fraction -
         the equivalence CCX_DAMAGE_PATH mode 1 was written to verify and
         did, at max|ramp-tempload| = 0.

         lambda still needs an EQUATION, and that is the bordered
         dissipation solve, so this requires CCX_DISSIPATION_CONTROL=2.
         Without it there is nothing to determine lambda and the flag is
         refused rather than silently doing something else.

         Measured motivation: soft50 with a released interface stops at 8.8%
         of peak on a force residual that fails at 81% NORMAL nodes, with a
         stable average force and a clean topology - a limit point, not an
         artefact (E-95).  m12_spc3 stops the same way at a healthy grip
         node (E-74). */
      /* The target does double duty: it sizes the step AND, through
         DAMAGE_DISS_ENGAGE, decides when the constraint takes over.  Those
         two pull in opposite directions - a target matched to the real
         dissipation rate engages while the response is still elastic, which
         is precisely the state the coupled solve cannot represent (dG is
         identically zero there).  Measured on the released-interface soft50:
         target 0.30 against a real post-peak rate of 0.75-2.2 per increment
         engaged at increment 9, lambda=0.0155, and died at increment 10.

         So give engagement its own control: stay in displacement control
         until the step time passes this value, whatever the dissipation is
         doing.  Set it past the load peak. */
      /* Separate the two halves of dissipation control.

         CCX_DISSIPATION_TARGET drives THREE things at once: it sizes the
         next increment, it gates engagement through DAMAGE_DISS_ENGAGE, and
         it is the right-hand side of the coupled constraint.  Measured on
         the released-interface soft50, a target matched to the real
         dissipation rate (0.75-2.2 per increment) makes the step controller
         grow the step by its 1.25 cap every increment and thrash into tmin
         at the load peak: step range 581x against the stock controller's
         64x, floor 2.83e-6 against 3.13e-5, and the run died with ZERO
         deletions before the constraint ever engaged (E-96).

         CCX_DISSIPATION_STEP=0 keeps the stock step controller and leaves
         the target to the constraint alone. */
      if((damage_de13_env=ccxopt_getenv("CCX_DISSIPATION_STEP"))!=NULL){
        lc.diss_step=(strcmp(damage_de13_env,"0")==0)?0:1;
        if(lc.diss_step==0){
          printf("[DISSIPATION] step sizing is OFF; the target drives the "
                 "constraint only, the stock controller sizes the step\n");
        }
      }
      if((damage_de13_env=ccxopt_getenv("CCX_DISSIPATION_ENGAGE_T"))!=NULL){
        lc.diss_engage_t=atof(damage_de13_env);
      }
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_ARCLENGTH"))!=NULL){
        lc.arc=(strcmp(damage_de13_env,"0")==0)?0:1;
      }
      if(lc.arc==1){
        if(lc.diss_ctrl!=2){
          printf("[PATH FOLLOWING] *WARNING: CCX_DAMAGE_ARCLENGTH needs "
                 "CCX_DISSIPATION_CONTROL=2 to have an equation for lambda; "
                 "it is ignored\n");
          lc.arc=0;
        }else{
          printf("[PATH FOLLOWING] lambda is decoupled from the step time "
                 "and may DECREASE; theta stays monotone for dtime\n");
        }
      }

      damage_epol.filter=ccxopt_getenv("CCX_DAMAGE_DELETE_MAT");

      damage_de13_env=ccxopt_getenv("CCX_DAMAGE_DELETE_D");
      if(damage_de13_env!=NULL){
        damage_epol.delete_d=atof(damage_de13_env);
        if((damage_epol.delete_d<0.5)||(damage_epol.delete_d>0.9999)){
          printf("[DAMAGE DE1.3.1] *WARNING: CCX_DAMAGE_DELETE_D=%s is "
                 "outside [0.5,0.9999]; keeping %.4f\n",
                 damage_de13_env,(double)DAMAGE_DE13_DELETE_D);
          damage_epol.delete_d=DAMAGE_DE13_DELETE_D;
        }
      }

      damage_topology_env=ccxopt_getenv("CCX_DAMAGE_TOPOLOGY");
      if((damage_topology_env!=NULL)&&
         ((strcmp(damage_topology_env,"DEFERRED")==0)||
          (strcmp(damage_topology_env,"deferred")==0)||
          (strcmp(damage_topology_env,"1")==0))){
        damage_topology_deferred_mode=1;
      }
    }
  }

  /* Path following is armed HERE, outside the *DAMAGE INITIATION block
     above.  It was inside it, which silently disabled it for any model
     without a damage material - including the cohesive benchmark, where
     the only nonlinearity is a UC6 interface and ndmat_ is zero.  The
     method has nothing to do with whether a bulk damage material exists. */

  /* ---- CCX_PATHFOLLOW ------------------------------------------
     A self-contained dissipation path following, independent of the
     CCX_DISSIPATION_CONTROL machinery above and of PARDISO.  It is
     armed only inside a domain where the bordered system is exactly
     the one pathfollow.c verifies, and it refuses to arm otherwise
     rather than degrade silently. */

  if(ccxopt_getenv("CCX_DAMAGE_BATCH_TRACE")!=NULL) td_trace=1;
  if(ccxopt_getenv("CCX_TOPODIAG")!=NULL) td_from=atoi(ccxopt_getenv("CCX_TOPODIAG"));
  if(ccxopt_getenv("CCX_DAMAGE_LS_PROBE")!=NULL) rsc.ls_probe=1;
  if(ccxopt_getenv("CCX_DAMAGE_WALL_NULLVEC")!=NULL) prb.wall_null=1;

  /* [WALLDIAG] Arm the increment-numbered probes by LOAD FACTOR instead.

     The wall is reproducible in theta and NOT in the increment index: the
     same binary on the same deck lands on it at increment 348 with two
     threads and 353 with four, because MKL_CBWR fixes the instruction set
     and not the reduction order.  Naming an increment therefore arms the
     diagnostics at a state that the next run does not have.  theta is the
     physical coordinate the wall is recorded in, so the gate takes it. */

  if(ccxopt_getenv("CCX_DAMAGE_WALL_THETA")!=NULL){
    prb.wall_theta=atof(ccxopt_getenv("CCX_DAMAGE_WALL_THETA"));
    printf("[WALLDIAG] armed at theta >= %.9e: from the first attempt at or "
           "beyond that load factor the linearisation check, the topology "
           "diagnostic, the deflated null-vector probe and the line-search "
           "ladder probe run on the CURRENT increment, whatever its "
           "number.%s",prb.wall_theta,"\n");
    fflush(stdout);
  }
  /* [WALLDIAG] MASKSTEP.  The wall measurement says 99.96% of |p_N|^2 sits
     on the three dofs of ONE node whose assembled diagonal has collapsed,
     while the force imbalance is two nodes away.  This probe walks the same
     eps ladder along p_N with every component on an AUTOSPC-masked node
     zeroed, so the question "is the step unusable BECAUSE it is spent on
     collapsed-diagonal nodes" becomes a number instead of a story.  It reads
     the mask AUTOSPC already builds, changes no solution path, and is off
     unless asked for. */
  /* [DAMAGE AUTOSPC-FORCE] AUTOSPC already decides which nodes have lost
     their load path - assembled diagonal below CCX_DAMAGE_AUTOSPC times the
     node's OWN intact value - and removes them from the DISPLACEMENT
     convergence norm only, deliberately leaving the force residual alone so
     that "a node that still carries load still blocks convergence".

     At the second wall that reasoning inverts.  The peak force residual sits
     on a node held by ONE live bulk element whose six cohesive facets have
     all failed (g -> 0), it needs a displacement of order one to be
     equilibrated, and it moves 0.1% per Newton iteration.  It does not carry
     load; it blocks convergence anyway, and it blocks the very increment in
     which its own last element would damage and delete.  That is a deadlock:
     the solve cannot advance because of the fragment, and the fragment
     cannot resolve because the solve cannot advance.

     This lifts the veto for exactly the AUTOSPC set and for nothing else.
     The dof is still assembled, still solved and still moved - only its
     veto over ram[0] is removed - and every excluded residual is printed, so
     it can never hide a growing imbalance.  OFF unless asked for. */
  if(ccxopt_getenv("CCX_DAMAGE_AUTOSPC_FORCE")!=NULL){
    damage_spc_force=1;
    printf("[DAMAGE AUTOSPC-FORCE] armed: a node already masked by AUTOSPC "
           "is excluded from the FORCE residual ram[0] as well as from "
           "cam[0].  Its dof is still solved and still moved; the largest "
           "excluded residual is reported on every increment.%s","\n");
    fflush(stdout);
  }
  if(ccxopt_getenv("CCX_DAMAGE_WALL_MASKSTEP")!=NULL){
    prb.wall_maskstep=1;
    printf("[WALLDIAG] MASKSTEP armed: the linearisation check gains a third "
           "pass along p_N with the AUTOSPC-masked nodes' components zeroed. "
           " Diagnostic only - nothing on a solution path reads it.%s","\n");
    fflush(stdout);
  }
  /* [CONVERGE] the reduction is not optional - every verdict in the run is
     made from the numbers it produces - so its self test runs on every run
     and a failure stops the job rather than degrading it quietly.  There is
     no fallback to degrade to: a wrong ram[0] is a wrong answer that looks
     like a right one. */
  /* [TOPOLOGY] the erosion transaction arms here, next to Convergence and
     under the same rule: the self test runs on every run and a failure
     stops the job.  A transaction whose discard leaks or double-frees
     corrupts a run that still prints plausible numbers. */
  glob_census_init(&damage_glob);
  glob_census_arm(&damage_glob);

  if(trial_check(&nlgt)!=0){
    printf("*ERROR: the trial-evaluation context has unbound fields.\n"
           "        Every field is the address of a local, so a NULL\n"
           "        means TRIAL_BIND and the struct have drifted\n"
           "        apart and results() would be called with one.\n");
    fflush(stdout);
    FORTRAN(stop,());
  }

  erosion_batch_init(&damage_ebatch);

  topo_txn_init(&dtxn);

  converge_init(&damage_cvg,damage_qam_floor,0,NULL,0);

  if((td_trace!=0)||(td_from>0)){
    {
      td_armed=1;
      NNEW(td_comp,ITG,*nk);
      NNEW(td_sort,ITG,(*ne>0)?*ne:1);
      printf("[TOPODIAG] armed: batch trace %s, topology report from "
             "increment %" ITGFORMAT ".  Nothing here changes the "
             "calculation.\n",td_trace?"ON":"off",td_from);
    }
  }

  pathdrv_configure(&pf,&lc,&nlgt,isolver,ncont);

  /* backup fields needed to roll back element deletion after
     an unsuccessful damage re-equilibration */

  if((*ndmat_>0)&&(*iexpl<=1)){
    NNEW(ipkondamageini,ITG,ne0);
    NNEW(damdamageini,double,mi[0]*ne0);
    NNEW(damde1prev,double,mi[0]*ne0);
    if(damage_de12_enabled){
      /* UNSYM stage 2 buffer: 6 effective-stress plus 6 dD/d(eps)
         components per integration point, filled by resultsmech.f and
         consumed by mafilldamas.f */

      if(damage_tangent_mode==2){
        NNEW(damage_damjac,double,12*mi[0]*ne0);
        /* why each element did or did not get the rank-1 term, kept so
           that the structural operator probe can attribute a discrepancy
           on a NAMED element to a named population */
        NNEW(damage_damcat,ITG,ne0);
      }

      /* R4 viscous regularisation buffers.  damage_damvisc is the trial
         relaxed degradation rebuilt from damage_damviscini on every
         results() call; damage_damviscini is the committed value, taken
         over at the start of the next physical increment exactly like
         damdamageini.  A rejected increment simply discards the trial. */

      damage_visc_env=ccxopt_getenv("CCX_DAMAGE_VISCOSITY");
      if(damage_visc_env!=NULL) damage_visc_eta=atof(damage_visc_env);
      if(damage_visc_eta<0.) damage_visc_eta=0.;
      if(damage_visc_eta>0.){
        NNEW(damage_damvisc,double,mi[0]*ne0);
        NNEW(damage_damviscini,double,mi[0]*ne0);
      }
      NNEW(damage_de13_trigger_value,double,ne0);
      NNEW(damage_de13_trigger_ip,ITG,ne0);
      for(i=0;i<ne0;i++){
        damage_de13_trigger_value[i]=-1.;
        damage_de13_trigger_ip[i]=0;
      }
    }
    if(*nmethod!=4) NNEW(veolddamageini,double,mt**nk);

    /* committed hard-deletion history; one file per CalculiX job */

    strcpy2(damagefilename,jobnamec,132);
    strcat(damagefilename,".damage");

    if(*istep==1){
      fdamage=fopen(damagefilename,"w");
    }else{
      fdamage=fopen(damagefilename,"a");
    }

    if(fdamage==NULL){
      printf(" *ERROR in nonlingeo: cannot open %s for writing\n",
             damagefilename);
      FORTRAN(stop,());
    }

    printf("[DAMAGE PATCH A3] adaptive-fast transactional damage enabled "
           "(batch_max=%d floor=%.2f direct_alpha=%.2f)\n",
           DAMAGE_FAST_BATCH_MAX,DAMAGE_FAST_MIN_FRACTION,
           DAMAGE_FAST_DIRECT_ALPHA);
    if(damage_de12_enabled){
      printf("[DAMAGE DE1.2] Newton-integrated trial damage enabled "
             "(%" ITGFORMAT " material(s)); BK1 one-pass local update / "
             "one global Newton solve; exact cell VTK enabled\n",
             damage_de12_matcount);
      if(damage_dm20_matcount>0){
        printf("[DAMAGE DM2.0] tabulated ductile initiation enabled "
               "(%" ITGFORMAT " material(s)); eps_f(eta) linear interpolation\n",
               damage_dm20_matcount);
      }
      printf("[DAMAGE DE1.3.1] terminal failure enabled "
             "(Ddelete=%.4f batch_max=%d); transactional A3 topology "
             "backend\n",damage_epol.delete_d,DAMAGE_DE13_BATCH_MAX);
      printf("[DAMAGE SOLVER NC1] adaptive slow-Newton controller enabled "
             "(stock_ic=%" ITGFORMAT " cap=%d max_extra=%d); "
             "final stock convergence criteria unchanged\n",
              (ITG)icref,DAMAGE_SLOW_NEWTON_MAX_ITERS,
              DAMAGE_SLOW_NEWTON_MAX_EXTRA);
      /* crack-band check: L*sigma_y/(u_f*E) must stay below 1 or the
         element softens faster than it can elastically unload */

      FORTRAN(damsnapcheck,(co,kon,ipkon,lakon,&ne0,ielmat,mi,ndmcon,
                            dmcon,ndmat_,ntmat_,elcon,ncmat_,plicon,
                            nplicon,npmat_,&damage_snap_ratio,
                            &damage_snap_elem,&damage_snap_bad));
      printf("[DAMAGE CRACK-BAND] max L*sigma_y/(u_f*E) = %.4f in "
             "element %" ITGFORMAT "\n",damage_snap_ratio,
             damage_snap_elem);
      if(damage_snap_bad>0){
        printf("[DAMAGE CRACK-BAND] *WARNING: %" ITGFORMAT " element(s) "
               "have r>1, i.e. the element-level softening branch is "
               "steeper than its own elastic unloading.  This is an "
               "indicator, not a proof: the structure may still find "
               "equilibrium by redistributing into neighbours.  It "
               "does mean convergence with a consistent tangent can "
               "be poor; refining the mesh or raising u_f removes "
               "it.\n",damage_snap_bad);
      }

      if(damage_tangent_mode==2){
        printf("[DAMAGE TANGENT UNSYM] stage 1: asymmetric assembly and "
               "solve with the unchanged symmetric g(D)*Cep tangent; "
               "results must match the symmetric path exactly\n");

        /* The rank-1 projection is the one piece of this mode the
           structural FD probe cannot exonerate on its own: a discrepancy
           it measures could equally be damjac, the projection, or the
           probe's reading of the CSR.  damrank1test settles the middle one
           offline, against e_c3d's own quadruple sum and against a finite
           difference of the element internal force.  It now runs in
           selftest_gate() with the other eighteen, before anything is
           armed, and a failure stops the run rather than disarming this
           mode - every number the rest of the run would produce comes from
           the same build. */
      }else{
        printf("[DAMAGE TANGENT BK2] secant g(D)*Cep baseline enabled - "
               "this is the DEFAULT and the measured best on both fast "
               "decks; set CCX_DAMAGE_TANGENT=UNSYM for the consistent "
               "tangent, which buys 0.2%% of iterations for 37%% of "
               "runtime\n");
      }
      if(rsc.reeq_scale_mode==1){
        printf("[DAMAGE SOLVER NC2] terminal same-load correction scale "
               "uses the converged physical-increment displacement norm; "
               "stock residual/correction tolerances retained\n");
      }else{
        printf("[DAMAGE SOLVER NC2] stock zero-load re-equilibration scale; "
               "set CCX_DAMAGE_REEQ_SCALE=PHYSICAL for NC2 trial\n");
      }
      if(rsc.reeq_uam_floor>0.){
        printf("[DAMAGE REEQ FLOOR] DIAGNOSTIC: the re-equilibration "
               "displacement reference is held at or above %.3e of its "
               "running maximum, so the RELATIVE correction criterion "
               "cannot collapse as the step is cut.  The perturbation of a "
               "same-load topology solve does not scale with dt, so cutting "
               "the step cannot converge it (J-09).  This CHANGES THE "
               "CONVERGENCE CRITERION and therefore the answer; going "
               "further with it is not by itself a success\n",
               rsc.reeq_uam_floor);
      }
      if(rsc.linesearch_mode==1){
        printf("[DAMAGE SOLVER BK3] adaptive line search enabled "
               "(growth=%.2f lambda=[%.2f,%.2f] trials=%d); full Newton retained "
               "while residual contracts\n",DAMAGE_LINESEARCH_GROWTH,
               DAMAGE_LINESEARCH_MIN,DAMAGE_LINESEARCH_MAX,
               DAMAGE_LINESEARCH_MAX_TRIALS);
      /* The damage line search is clamped to a compiled-in floor and trial
         cap.  On m12_field_soft50 both bind SIMULTANEOUSLY on every iteration
         of the failing increment - lambda pinned at exactly 0.100 with
         trials=3, seven iterations running, residual oscillating in
         1.70e-3..2.30e-3 with no contraction and a stable active set (E-81).
         A search that asks for a shorter step on every iteration and is
         refused on every iteration is a search whose floor is the binding
         constraint, not a search that has converged.

         The structural FD probe says the tangent there is NOT the problem:
         2.52e-3 median at the wall against a 1.33e-4 elastic noise floor on
         the same deck and SP1's healthy 1.92e-3 (E-91).  So make the two
         clamps measurable instead of assumed.  Defaults are the compiled-in
         values, so an unset environment is bit-identical. */
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_LS_MIN"))!=NULL){
        rsc.ls_min=atof(damage_de13_env);
        if(rsc.ls_min<1.e-6) rsc.ls_min=1.e-6;
        if(rsc.ls_min>DAMAGE_LINESEARCH_MAX) rsc.ls_min=DAMAGE_LINESEARCH_MAX;
      }
      if((damage_de13_env=ccxopt_getenv("CCX_DAMAGE_LS_TRIALS"))!=NULL){
        rsc.ls_trials=atoi(damage_de13_env);
        if(rsc.ls_trials<1) rsc.ls_trials=1;
        if(rsc.ls_trials>32) rsc.ls_trials=32;
      }

      /* ---- THE BACKTRACKING LADDER ---------------------------------
         Measured on the target with CCX_DAMAGE_LS_PROBE, over the 245
         line-search activations of a run to the recorded wall:

           in 114 of them (46.5%) the best step length lies BELOW the
           floor of 0.1, i.e. outside what the ladder {1.0, 0.5, 0.1}
           can reach;

           at the wall increment 353, iterations 12 to 16, the floor step
           RAISES |R|inf - 2.852e-03, 3.053e-03, 3.353e-03, 3.649e-03,
           3.822e-03 - while alpha=0.03 lowers it at every one of them.
           The iteration was converging until iteration 11 and the ladder
           turned it round.

         So the ladder is deepened (floor 1e-3, up to 8 halvings) and the
         search now falls back on the BEST alpha it measured rather than
         the last.  Neither changes any convergence criterion: what an
         increment must satisfy to be accepted is untouched, this only
         chooses how far along a direction to step.

         CCX_DAMAGE_LS_LEGACY=1 restores the old ladder exactly, which is
         what makes a controlled A/B possible from ONE binary. */

      /* The ladder is the thing that was wrong, so prove it is right
         before using it, on every run, the way pathfollow and
         crackcontrol do.  A failure here is not something to carry on
         through: the search would silently go back to handing back a
         step it had measured to be worse. */


      if(ccxopt_getenv("CCX_DAMAGE_LS_LEGACY")!=NULL){
        rsc.ls_legacy=1;
        rsc.ls_min=DAMAGE_LINESEARCH_MIN;
        rsc.ls_trials=DAMAGE_LINESEARCH_MAX_TRIALS;
        printf("[DAMAGE LINESEARCH] LEGACY ladder: floor %.2f, %d trials, "
               "and the last trial is taken whether or not it contracts.\n",
               DAMAGE_LINESEARCH_MIN,DAMAGE_LINESEARCH_MAX_TRIALS);
      }else{
        if(ccxopt_getenv("CCX_DAMAGE_LS_MIN")==NULL) rsc.ls_min=1.e-3;
        if(ccxopt_getenv("CCX_DAMAGE_LS_TRIALS")==NULL) rsc.ls_trials=8;
        printf("[DAMAGE LINESEARCH] ladder: floor %.4e, up to %"
               ITGFORMAT " trials, and the BEST alpha measured is taken "
               "when none contracts.\n",rsc.ls_min,rsc.ls_trials);
      }
      if((rsc.ls_min!=DAMAGE_LINESEARCH_MIN)||
         (rsc.ls_trials!=DAMAGE_LINESEARCH_MAX_TRIALS)){
        printf("[DAMAGE LINESEARCH] floor %.4e (default %.2f), "
               "max trials %" ITGFORMAT " (default %d)\n",
               rsc.ls_min,DAMAGE_LINESEARCH_MIN,
               rsc.ls_trials,DAMAGE_LINESEARCH_MAX_TRIALS);
      }

      }else{
        printf("[DAMAGE SOLVER BK3] disabled; set "
               "CCX_DAMAGE_LINESEARCH=ADAPTIVE for damage globalization\n");
      }
      if(damage_topology_deferred_mode==1){
        printf("[DAMAGE TOPOLOGY BK4] deferred sparse-structure compaction "
               "enabled; terminal elements leave the assembly immediately, "
               "remastruct is deferred until an active node becomes orphaned\n");
      }else{
        printf("[DAMAGE TOPOLOGY BK4] immediate remastruct baseline; set "
               "CCX_DAMAGE_TOPOLOGY=DEFERRED for deferred compaction\n");
      }
    }else{
      printf("[DAMAGE DE1.1] fixed-point controller available for legacy "
             "DE1 states; exact cell VTK enabled\n");
    }
    fflush(stdout);

    if(*istep==1){
      fprintf(fdamage,"# CalculiX committed damage deletion history v2\n");
      fprintf(fdamage,"# element step increment step_time total_time "
              "material damage critical_ip batch\n");
      fprintf(fdamage,"# DE1.3.1: damage is the captured terminal-trigger D "
              "at first topology marking\n");
      fflush(fdamage);
    }
  }
  
  /*********************************************************************/
  
  /* calculating of the acceleration due to force discontinuities
     (external - internal force) at the start of a step */
  
  /*********************************************************************/
  
  if((*nmethod==4)&&(*ithermal!=2)&&(icfd==0)){
    bet=(1.-alpha[0])*(1.-alpha[0])/4.;
    gam=0.5-alpha[0];
      
    /* calculating the stiffness and mass matrix */
      
    reltime=0.;
    time=0.;
    dtime=0.;
      
    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
		      xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,ampli,
		      &time,&reltime,ttime,&dtime,ithermal,nmethod,xbounold,
		      xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,nodeforc,
		      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,
		      ilboun,nelemload,sideload,mi,ntrans,trab,inotr,veold,
		      integerglob,doubleglob,tieset,istartset,iendset,ialset,
		      ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,ipobody,
		      iponoeln,inoeln,ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,
		      rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));
      
    time=0.;
    dtime=1.;
    
    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's)
	if contact arises the number of MPC's can also change */
      
    cam[0]=0.;cam[1]=0.;cam[2]=0.;
      
    if(icascade==2){
      memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
      RENEW(nodempc,ITG,3*memmpcref_);
      for(k=0;k<3*memmpcref_;k++){
	nodempc[k]=nodempcref[k];}
      RENEW(coefmpc,double,memmpcref_);
      for(k=0;k<memmpcref_;k++){
	coefmpc[k]=coefmpcref[k];}
    }

    newstep=0;
    FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		       nmpc,ikboun,ilboun,nboun,xbounold,aux,iaux,
		       &maxlenmpc,ikmpc,ilmpc,&icascade,
		       kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
		       &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));
    if(icascade==2){
      for(k=0;k<3*memmpc_;k++){
	nodempcref[k]=nodempc[k];}
      for(k=0;k<memmpc_;k++){
	coefmpcref[k]=coefmpc[k];}
    }

    /* recalculating the matrix structure */
    
    if(icascade>0){
      remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		 &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		 labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		 kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		 neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		 &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		 typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		 nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		 ielmat,matname);
    }

    /* invert nactdof */

    NNEW(nactdofinv,ITG,1);
      
    iout=-1;
    ielas=1;
      
    MNEW(fn,double,mt**nk);
    NNEW(stx,double,6*mi[0]**ne);
      
    if((*iexpl<=1)||(*mortar==-1)){intscheme=1;}
      
    if(ne1d2d==1)NNEW(inum,ITG,*nk);
    results(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xbounold,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,&bet,
	    &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,xstaten,
	    eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	    ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	    &ne0,thicke,shcon,nshcon,
	    sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
	    inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavquadel,aut,irowt,jqt,&mortartrafoflag,
	    &intscheme,physcon,dam,damn,iponoel);
      
    SFREE(fn);SFREE(stx);if(ne1d2d==1)SFREE(inum);
      
    if(*mortar<2){
      iout=0;
      ielas=0;
	  
      reltime=0.;
      time=0.;
      dtime=0.;
    }
      
    if(*iexpl>1){

      mscalmethod=0;
      nloadrhs=*nload;nbodyrhs=*nbody;
	  
      /* Explicit: Calculation of stable time increment according to
	 Courant's Law  Carlo Monjaraz Tec (CMT) and Selctive Mass Scaling CC*/

      /*Mass Scaling
	mscalmethod < 0: no explicit dynamics
	mscalmethod = 0: no mass scaling
	mscalmethod = 1: selective mass scaling for nonlinearity after 
	Olovsson et. al 2005

        mscalmethod=2 and mscalmethod=3 correspond to 0 and 1, 
        respectively with in addition contact scaling active; contact
        scaling is activated if the user time increment cannot be satisfied */

      dtset=*tmin*(*tper);
      NNEW(smscale,double,*ne);
	  
      FORTRAN(calcstabletimeincvol,(&ne0,elcon,nelcon,rhcon,nrhcon,alcon,
				    nalcon,orab,ntmat_,ithermal,alzero,plicon,
				    nplicon,plkcon,nplkcon,npmat_,mi,&dtime,
				    xstiff,ncmat_,vold,ielmat,t0,t1,matname,
				    lakon,wavespeed,nmat,ipkon,co,kon,&dtvol,
				    alpha,smscale,&dtset,&mscalmethod,mortar,
				    jobnamef,iperturb));

      printf(" Explicit time integration: Volumetric COURANT initial stable time increment:%e\n\n",dtvol);

      if(dtvol>(*tmax*(*tper))){
	*tinc=*tmax*(*tper);}
      else if(dtvol<dtset){
	*tinc=dtset;}
      else{
	*tinc=dtvol;
      }
	  
      dtheta=(*tinc)/(*tper);
      dthetaref=dtheta;
      printf(" SELECTED time increment (not considering penalty contact):%e\n\n",*tinc);
    }
      
    /* in mafillsm the stiffness and mass matrix are computed;
       The primary aim is to calculate the mass matrix (not 
       lumped for an implicit dynamic calculation, lumped for an
       explicit dynamic calculation). However:
       - for an implicit calculation the mass matrix is "doped" with
       a small amount of stiffness matrix, therefore the calculation
       of the stiffness matrix is needed.
       - for an explicit calculation the stiffness matrix is not 
       needed at all. Since the calculation of the mass matrix alone
       is not possible in mafillsm, the determination of the stiffness
       matrix is taken as unavoidable "ballast". */
      
    NNEW(ad,double,neq[1]);
    NNEW(au,double,nzs[1]);

    mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
		 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		 nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		 nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		 nmethod,ikmpc,ilmpc,ikboun,ilboun,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		 ielmat,ielorien,norien,orab,ntmat_,
		 t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		 xstiff,npmat_,&dtime,matname,mi,
		 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		 &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		 xstateini,xstate,thicke,integerglob,doubleglob,
		 tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		 pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		 iponoeln,inoeln,network,ntrans,inotr,trab,smscale,&mscalmethod,
		 set,nset,islavquadel,aut,irowt,jqt,&mortartrafoflag);
    
    if(*nmethod==0){
	  
      /* error occurred in mafill: storing the geometry in frd format */
	  
      ++*kode;
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	MNEW(ipneigh,ITG,*nk);
      }
	  
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti,damn,&errn);
	  
      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}      
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif	  
      if(nmethodold==0){FORTRAN(stopwithout201,());}else{FORTRAN(stop,());}
	  
    }

    /* massless contact: setting up the system matrices based on
       the stiffness and mass matrix; these matrices are not
       assumed to change during the step;
       factorization of the LHS matrix */
    
    if(*mortar==-1){
      if(ncont!=0){
	nmasts=nmastnode[*ntie];

	NNEW(kslav,ITG,3**nslavs);
	NNEW(lslav,ITG,3**nslavs);
	NNEW(ktot,ITG,3**nslavs+3*nmasts);
	NNEW(ltot,ITG,3**nslavs+3*nmasts);
	NNEW(fric,double,*nslavs);

	/*  Create set of slave and slave+master contact DOFS (sorted);
	    assign a friction coefficient to each slave node */
      
	FORTRAN(create_contactdofs,(kslav,lslav,ktot,ltot,nslavs,islavnode,
				    &nmasts,imastnode,nactdof,mi,&neqtot,
				    nslavnode,fric,tieset,tietol,ntie,elcon,
				    ncmat_,ntmat_));
      }else{
	neqtot=0;
      }
      
      /*   RENEW(kslav,ITG,3**nslavs);
      RENEW(lslav,ITG,3**nslavs);
      RENEW(ktot,ITG,neqtot);
      RENEW(ltot,ITG,neqtot);*/

      /* create RHS of system:  M/dt - (aM + bK)/2 and store in adc,auc */
      
      NNEW(adc,double,neq[0]);
      for(k=0;k<neq[0];k++){
	adc[k]=adb[k]/(*tinc)-(alpham*adb[k]+betam*ad[k])/2.0;
      }

      NNEW(auc,double,nzs[0]);
      for(k=0;k<nzs[0];k++){
	auc[k]=aub[k]/(*tinc)-(alpham*aub[k]+betam*au[k])/2.0;
      }

      /* create LHS of system:  M/dt + (aM + bK)/2  and store in adb,aub */

      for(k=0;k<neq[0];k++){
	adb[k]=adb[k]/(*tinc)+(alpham*adb[k]+betam*ad[k])/2.0;
      }

      for(k=0;k<nzs[0];k++){
	aub[k]=aub[k]/(*tinc)+(alpham*aub[k]+betam*au[k])/2.0;
      }

      /* reduce LHS and RHS by removing contact dofs (diagonal terms
         are set to 1, off-diagonal terms to 0 */

      if(ncont!=0){
	FORTRAN(reducematrix,(aub,adb,jq,irow,neq,&neqtot,ktot));
	FORTRAN(reducematrix,(auc,adc,jq,irow,neq,&neqtot,ktot));
      }

      /* factorize the LHS */

      if(*isolver==0){
#ifdef SPOOLES

	spooles_factor(adb,aub,adb,aub,&sigma,icol,irow,
		       &neq[0],&nzs[0],&symmetryflag,&inputformat,&nzs[0]);

#else
	printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],token);
#else
	printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(adb,&aub,adb,aub,&sigma,icol,&irow,&neq[0],&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
		       &symmetryflag,&inputformat,jq,&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			   &symmetryflag,&inputformat,jq,&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }

      // Storing contact force vector initial solution

      if(ncont!=0){
	NNEW(aloc,double,3**nslavs);
	NNEW(alglob,double,neqtot);
      }

      /* no nlgeom and no nonlinear material for massless explicit dynamics */
      
      if((iperturb[0]<3)&&(iperturb[1]==0)) masslesslinear=1;

      /* check whether the output consists of displacements only */

      FORTRAN(checkdispoutonly,(prlab,nprint,nlabel,filab,&idispfrdonly));

      if(idispfrdonly==1){
	NNEW(inumcp,ITG,*nk);
	strcpy1(&cflag[0],&filab[4],1);
	FORTRAN(createinum,(ipkon,inumcp,kon,lakon,nk,ne,&cflag[0],nelemload,
			    nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,
			    ielmat,ielprop,prop));
      }

      
    } //endif massless
      
    /* mass x acceleration = f(external)-f(internal) 
       only for the mechanical loading*/
      
    /* not needed for massless contact */
    
    if(*mortar!=-1){
      for(k=0;k<neq[0];++k){b[k]=fext[k]-f[k];}
    }
      
    if(*iexpl<=1){
	  
      /* a small amount of stiffness is added to the mass matrix
	 otherwise the system leads to huge accelerations in 
	 case of discontinuous load changes at the start of the step */
	  
      dtime=*tinc/10.;
      scal1=bet*dtime*dtime*(1.+alpha[0]);
      for(k=0;k<neq[0];++k){
	ad[k]=adb[k]+scal1*ad[k];
      }
      for(k=0;k<nzs[0];++k){
	au[k]=aub[k]+scal1*au[k];
      }
      if(*isolver==0){
#ifdef SPOOLES
	spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		&symmetryflag,&inputformat,&nzs[2]);
#else
	printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
#endif
      }
      else if((*isolver==2)||(*isolver==3)){
	preiter(ad,&au,b,&icol,&irow,&neq[0],&nzs[0],isolver,iperturb);
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
#else
	printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[0],&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		     &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
	printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
	printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
    }
      
    else{

      /* explicit dynamics; no selective mass scaling
         (at most spring scaling) */
      
      /* if massless contact: no acceleration needed */
      
      if((mscalmethod==0)||(mscalmethod==2)){
        if(*mortar!=-1){
	  for(k=0;k<neq[0];++k){b[k]=(fext[k]-f[k])/adb[k];}
	}
      }
	  
      else{

	/* explicit dynamics with selective mass scaling */

	inputformat=0;
	if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,&nzs[2]);
	  spooles_solve(b,&neq[0]);
#else
	  printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],token);
	  sgi_solve(b,token);
#else
	  printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(adb,&aub,adb,aub,&sigma,icol,&irow,&neq[0],&nzs[0]);
	  tau_solve(b,&neq[0]);
#else
	  printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[0]);

	  pardiso_solve(b,&neq[0],&symmetryflag,&inputformat,&nrhs);
#else
	  printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  pastix_factor_main(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			&symmetryflag,&inputformat,jq,&nzs[0]);

	  pastix_solve(b,&neq[0],&symmetryflag,&nrhs);
#else
	  printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
      }
    }
      
    /* for thermal loading the acceleration is set to zero */
      
    for(k=neq[0];k<neq[1];++k){
      b[k]=0.;
    }
      
    /* calculating the displacements, stresses and forces */
      
    if(*mortar!=-1){
      NNEW(v,double,mt**nk);
      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
      
      NNEW(stx,double,6*mi[0]**ne);
      MNEW(fn,double,mt**nk);
      
      /* setting a "special" time consisting of the first primes;
	 used to recognize the initial acceleration procedure
	 in file resultsini.f */

      if(ne1d2d==1)NNEW(inum,ITG,*nk);
      dtime=1.235711130e-20;
      trial_results(&nlgt);
      if(ne1d2d==1)SFREE(inum);
      dtime=0.;

      isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;	    
	cpypardou(sti,stx,&isiz,&num_cpus);
      }

      SFREE(v);SFREE(stx);SFREE(fn);
    }
    SFREE(ad);SFREE(au);
      
    /* the mass matrix is kept for subsequent calculations, therefore,
       no new mass calculation is necessary for the remaining iterations
       in the present step */
      
    mass[0]=0;intscheme=0;
    energyref=energy[0]+energy[1]+energy[2]+energy[3];

    if(*iexpl<=1){
	  
      NNEW(tmp,double,neq[1]);
      NNEW(adblump,double,neq[1]);
      for(k=0;k<neq[1];k++){
	tmp[k] = 1;
      }
      if(nasym==0){
	opmain(&neq[1],tmp,adblump,adb,aub,jq,irow); 
      }else{
	FORTRAN(opas,(&neq[1],tmp,adblump,adb,aub,jq,irow,nzs)); 
      }
      SFREE(tmp);
    }
  }

  /* warning: for C3D8R-elements the stiffness is needed in subroutine
              hgforce, therefore, in explicit dynamic steps
              with C3D8R-elements icmd should not be set to 3 */
  
  if(*iexpl>1) icmd=3;

  
  /**************************************************************/
  /* starting the loop over the increments                      */
  /**************************************************************/
  
  newstep=1;
	  
  //    MPADD start
  if((*nmethod==4)&&(*ithermal<2)&&(*iexpl<=1)){
    neini=*ne;
    for(k=0;k<4;k++){
      energystartstep[k]=energy[k];
    }
    emax=0.1*energyref;
    // Anti-stick at the beginning of simulation
  } 
  //    MPADD end

  /* saving the distributed loads (volume heating will be
     added because of friction heating) */

  if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
    nloadref=*nload;
    NNEW(nelemloadref,ITG,2**nload);
    if(*nam>0) NNEW(iamloadref,ITG,2**nload);
    NNEW(sideloadref,char,20**nload);
      
    isiz=2**nload;cpyparitg(nelemloadref,nelemload,&isiz,&num_cpus);
    if(*nam>0){
      isiz=2**nload;cpyparitg(iamloadref,iamload,&isiz,&num_cpus);
    }
    memcpy(&sideloadref[0],&sideload[0],sizeof(char)*20**nload);
  }
  
  while(((1.-theta>1.e-6)||(negpres==1))&&
        (damage_fracture_complete==0)){
      
    if((icutb==0)&&(idamagereeq==0)){

      /* ---- [DAMAGE CORR] one pass per converged increment ---------- */
      if(rsc.corr_mode==1){
        if((rsc.rec_used_in_inc==0)&&(lc.reg_on==0)){
          /* The crawl breaker needs a reference from HEALTHY operation.
             Taking it from the last clean increment before entry is
             wrong: that increment sits at the wall, where the stock
             controller has already cut dtheta hard.  Measured on
             run_c3_s3rad the reference came out 6.592e-06 while the run
             had been healthy at 5e-5..1.2e-4, so the mean inside the
             corridor (1.211e-05) was ABOVE it and the breaker could
             never fire.  Keep a ring of the last clean increments and
             use their MEDIAN instead. */
          damage_dtheta_healthy=dtheta;
          damage_dth_ring[damage_dth_i]=dtheta;
          damage_dth_i=(damage_dth_i+1)%20;
          if(damage_dth_n<20) damage_dth_n++;
        }
        if(rsc.corr_on==1){
          ITG cbrk=0;
          double cmean;
          rsc.corr_ninc++;
          rsc.corr_nsince++;
          rsc.corr_nwallstab=0;
          rsc.corr_dtsum+=dtheta;rsc.corr_dtn++;
          rsc.corr_nfact+=lc.reg_napply;
          /* A probe stays a PROBE until it has survived stableneed
             converged increments.  Clearing the flag after the FIRST
             success - the earlier draft did - makes a failure on the
             second or third increment look like a wall on the HELD
             lambda, so the wall handler would ESCALATE away from the
             proven value instead of falling back to it.  The flag and
             the promotion of lamstable must therefore happen together,
             and only once nsince has been reached. */
          if((rsc.corr_trial==1)&&
             (rsc.corr_nsince>=rsc.corr_stableneed)){
            rsc.corr_lamstable=rsc.corr_lam;
            rsc.corr_trial=0;
            printf("[DAMAGE CORR] lambda=%.3e has held %" ITGFORMAT
                   " increments and is now the PROVEN value%s",
                   rsc.corr_lam,rsc.corr_nsince,"\n");
          }
          if(rsc.corr_lam<=0.) rsc.corr_clean++;
          else rsc.corr_clean=0;
          cmean=(rsc.corr_dtn>0)?rsc.corr_dtsum/rsc.corr_dtn:0.;
          if(rsc.corr_clean>=rsc.corr_exit){
            printf("[DAMAGE CORR] EXIT: %" ITGFORMAT " consecutive "
                   "increments converged with NO help.  Corridor: %"
                   ITGFORMAT " increments, %" ITGFORMAT " intervention(s), "
                   "%" ITGFORMAT " regularized factorisation(s), step time "
                   "gained %.6e, mean dtime inside %.6e against reference "
                   "%.6e.  Stock Newton has the run back.%s",
                   rsc.corr_clean,rsc.corr_ninc,rsc.corr_nint,
                   rsc.corr_nfact,(theta-rsc.corr_theta0)**tper,
                   cmean**tper,rsc.corr_dtref**tper,"\n");
            fflush(stdout);
            rsc.corr_on=0;lc.reg_on=0;rsc.corr_lam=0.;
          }else{
            if(rsc.corr_ninc>rsc.corr_maxinc) cbrk=1;
            if((cbrk==0)&&(rsc.corr_ninc>rsc.corr_grace)&&
               (rsc.corr_dtref>0.)&&
               (cmean<rsc.corr_minfrac*rsc.corr_dtref)) cbrk=2;
            if(cbrk>0){
              printf("[DAMAGE CORR] CIRCUIT BREAKER (%s): %" ITGFORMAT
                     " increments, %" ITGFORMAT " intervention(s), %"
                     ITGFORMAT " regularized factorisation(s), step time "
                     "gained %.6e, mean dtime %.6e against reference %.6e."
                     "  Carried, not recovering: the corridor closes, the "
                     "mechanism DISARMS, and the next wall goes to the "
                     "original stock stop.%s",
                     (cbrk==1)?"length":"crawl",rsc.corr_ninc,
                     rsc.corr_nint,rsc.corr_nfact,
                     (theta-rsc.corr_theta0)**tper,cmean**tper,
                     rsc.corr_dtref**tper,"\n");
              fflush(stdout);
              rsc.corr_on=0;lc.reg_on=0;rsc.corr_lam=0.;
              rsc.rec_disarmed=1;
            }else{
              rsc.corr_try--;
              if((rsc.corr_try<=0)&&(rsc.corr_lam>0.)&&
                 (rsc.corr_trial==0)&&
                 (rsc.corr_nsince>=rsc.corr_stableneed)){
                rsc.corr_try=rsc.corr_tryevery;
                rsc.corr_lam=rsc.corr_lamstable*0.25;
                if(rsc.corr_lam<1.e-4) rsc.corr_lam=0.;
                rsc.corr_trial=1;
                rsc.corr_nsince=0;
                printf("[DAMAGE CORR] probing weaker help: lambda %.3e -> "
                       "%.3e (proven value kept for fallback)%s",
                       rsc.corr_lamstable,rsc.corr_lam,"\n");
                fflush(stdout);
              }
              lc.reg_lambda=rsc.corr_lam;
              lc.reg_on=(rsc.corr_lam>0.)?1:0;
            }
          }
          lc.reg_napply=0;
        }
      }

      /* [DAMAGE RESCUE] the previous increment converged.  If it was the
         rescue attempt, it is now committed together with any re-equilibration
         it needed: switch backtracking off at once and re-arm for the next
         INDEPENDENT wall.  theta has advanced, so this cannot loop. */

      /* [DAMAGE CT] acceptance needs BOTH the stock residual criteria and
         the constraint.  checkconvergence has already accepted on R; the
         constraint is verified here, before anything is committed. */
      if((ct.on==1)&&(fabs(ct.cprev)>ct.tolc)){
        printf("[DAMAGE CT] increment accepted on R but the CONSTRAINT is "
               "not satisfied: |c|=%.6e > tol_c=%.6e.  Requesting the clean "
               "PARTIAL exit rather than committing a state that does not "
               "lie on the constraint.%s",
               fabs(ct.cprev),ct.tolc,"\n");
        fflush(stdout);
        ct.partial=1;ct.on=0;
      }
      /* [DAMAGE CT] lambda >= 1 is NOT a completion in this MVP. */
      if((ct.on==1)&&(ct.lam>=1.)){
        printf("[DAMAGE CT] lambda has reached %.12e >= 1.  The MVP has no "
               "terminal landing, so this is a PARTIAL stop, not a "
               "completed step.%s",ct.lam,"\n");
        fflush(stdout);
        ct.partial=1;ct.on=0;
      }
      if(rsc.rec_used_in_inc==0){
        rsc.rec_healthy++;
        /* [DAMAGE TR] five consecutive clean increments, derived from
           INCREMENT NUMBERS rather than from a running flag.  The previous
           counter was incremented in one place and cleared in two, and its
           message carried no inc= token, so nothing reading the log could
           attribute it - both defects are fixed here.  dog.lasthelp is
           the increment of the last Rescue or trust-region firing of ANY
           level; iinc is the increment just committed, so the difference is
           the number of increments since then, and every one of them
           committed without help or lasthelp would have moved. */
        if((dog.narm>0)&&(dog.lasthelp>0)&&
           (iinc>dog.lasthelp)){
          dog.selfrec=iinc-dog.lasthelp;
          if((dog.selfrec>=5)&&
             (dog.recdone!=dog.lasthelp)){
            dog.recdone=dog.lasthelp;
            printf("[DAMAGE TR] SELF-RECOVERY inc=%" ITGFORMAT
                   ": increments %" ITGFORMAT "..%" ITGFORMAT
                   " - five consecutive - committed with PLAIN NEWTON, no "
                   "Rescue and no trust region, since the help at inc=%"
                   ITGFORMAT "; step time is now %.12e.  This, and not "
                   "t_end, is the criterion the method was built to meet.%s",
                   iinc,dog.lasthelp+1,iinc,dog.lasthelp,
                   theta**tper,"\n");
            fflush(stdout);
          }
        }
        if((rsc.rec_healthy==rsc.rec_window)&&
           (rsc.rec_unrec>0)){
          printf("[DAMAGE RESCUE] RECOVERED: %" ITGFORMAT " consecutive "
                 "increments converged with no intervention; the "
                 "un-recovered counter is cleared%s",
                 rsc.rec_window,"\n");
          fflush(stdout);
          rsc.rec_unrec=0;
        }
      }
      rsc.rec_used_in_inc=0;

      if(rsc.rescue_bt_on==1){
        rsc.rescue_bt_on=0;
        rsc.rescue_used=0;
        prb.evt_on=0;
        if(dog.on==1){
          printf("[DAMAGE TR] inc=%" ITGFORMAT " the TRUST-REGION attempt "
                 "CONVERGED at step time %.12e, judged by checkconvergence "
                 "on the UNMODIFIED residual.  The trust region switches OFF and the run "
                 "continues on plain Newton.  Steps kept in this attempt and "
                 "before it: %" ITGFORMAT " Newton, %" ITGFORMAT " Cauchy, %"
                 ITGFORMAT " dogleg; %" ITGFORMAT " iteration(s) accepted "
                 "nothing and kept the full Newton step; %" ITGFORMAT
                 " rejected trial(s); %" ITGFORMAT " residual evaluations; %"
                 ITGFORMAT " armed factorisation(s).  Whether the method "
                 "WORKED is decided by the next five increments, not by "
                 "this line.%s",
                 iinc,theta**tper,dog.nnewt,dog.ncau,
                 dog.ndog,
                 dog.nfail,dog.nrej,dog.neval,
                 dog.nfact,"\n");
          fflush(stdout);
        }
        dog.on=0;
        dog.delta=0.;
        if((lc.reg_on==1)&&(rsc.corr_mode==1)&&
           (rsc.corr_on==0)){
          rsc.corr_on=1;
          rsc.corr_lam=lc.reg_lambda;
          rsc.corr_lamstable=lc.reg_lambda;
          rsc.corr_ninc=0;rsc.corr_nint=1;
          rsc.corr_nfact=lc.reg_napply;
          rsc.corr_clean=0;rsc.corr_try=rsc.corr_tryevery;
          rsc.corr_dtsum=0.;rsc.corr_dtn=0;
          rsc.corr_trial=0;rsc.corr_nsince=0;
          rsc.corr_nwallstab=0;rsc.corr_nwalltot=0;
          rsc.corr_theta0=theta;
          {ITG mi1,mj1;double mtmp,msort[20];
           for(mi1=0;mi1<damage_dth_n;mi1++) msort[mi1]=damage_dth_ring[mi1];
           for(mi1=1;mi1<damage_dth_n;mi1++){
             mtmp=msort[mi1];
             for(mj1=mi1;(mj1>0)&&(msort[mj1-1]>mtmp);mj1--)
               msort[mj1]=msort[mj1-1];
             msort[mj1]=mtmp;
           }
           rsc.corr_dtref=(damage_dth_n>0)?msort[damage_dth_n/2]:
                                              damage_dtheta_healthy;
          }
          printf("[DAMAGE CORR] ENTERED at step time %.12e with "
                 "lambda=%.3e after %" ITGFORMAT " factorisation(s).  "
                 "Reference dtime: MEDIAN over the last %" ITGFORMAT
                 " clean increments = %.6e (the last clean increment "
                 "alone was %.6e and sits at the wall, so it is not a "
                 "healthy reference)%s",
                 theta**tper,rsc.corr_lam,lc.reg_napply,
                 damage_dth_n,rsc.corr_dtref**tper,
                 damage_dtheta_healthy**tper,"\n");
          fflush(stdout);
        }else if(lc.reg_on==1){
          printf("[DAMAGE REG] the regularized attempt CONVERGED after %"
                 ITGFORMAT " factorisation(s) at lambda=%.3e; the shift "
                 "is switched OFF and the run continues stock%s",
                 lc.reg_napply,lc.reg_lambda,"\n");
          fflush(stdout);
        }
        if(rsc.corr_on==0) lc.reg_on=0;
        rsc.rescue_nok++;
        printf("[DAMAGE RESCUE] ACCEPTED: the rescue increment converged at "
               "step time %.12e; backtracking switched OFF, rescue re-armed "
               "(%" ITGFORMAT " accepted of %" ITGFORMAT " fired)%s",
               theta**tper,rsc.rescue_nok,rsc.rescue_nfired,"\n");
        fflush(stdout);
      }
	  
      /* ---- [DAMAGE CT] ring, filled passively at the ONE commit point ---
         Six committed endpoints give five intervals.  Endpoint 0 is taken
         the first time this block runs; the next five commits complete the
         set.  Taken BEFORE vini<-vold and xstateini<-xstate, so vold/vini
         still bracket the increment that has just been committed and
         xstate still holds its converged values.  sti is used rather than
         stx for the compression flag: it is the committed stress copy and
         is alive here, whereas stx belongs to the iteration.
         Commit-only, so a cutback cannot corrupt it. */

      if(ct.mode==1){
        ITG cti,ctj,ctn;
        double ctd;
        if(ct.alloc==0){
          NNEW(ct.ring,double,18*mi[0]*ne0);
          NNEW(ct.fl,ITG,6*mi[0]*ne0);
          ct.alloc=1;ct.head=0;ct.nring=0;
          printf("[DAMAGE CT] ring allocated: 6 committed endpoints x %"
                 ITGFORMAT " integration points (%.1f MB).  It fills from "
                 "this commit onwards and changes nothing until arming.%s",
                 mi[0]*ne0,
                 (double)(18*mi[0]*ne0*8+6*mi[0]*ne0*4)/1048576.,"\n");
          fflush(stdout);
        }
        ctn=mi[0]*ne0;
        ct.head=(ct.nring==0)?0:((ct.head+1)%6);
        damcont_snap(co,kon,ipkon,lakon,vold,sti,xstate,ne0,mi[0],
                       *nstate_,mt,
                       &ct.ring[3*ctn*ct.head],
                       &ct.fl[ctn*ct.head]);
        if(ct.nring>0){
          ctd=0.;
          for(cti=0;cti<*nk;cti++)
            for(ctj=1;ctj<mt;ctj++)
              if(fabs(vold[mt*cti+ctj]-vini[mt*cti+ctj])>ctd)
                ctd=fabs(vold[mt*cti+ctj]-vini[mt*cti+ctj]);
          ct.duinf[ct.head]=ctd;
          ct.dt[ct.head]=dtime;
          ct.dlam[ct.head]=dtheta;
        }
        if(ct.nring<6) ct.nring++;

        /* ---- [DAMAGE CT] continuation commit lifecycle ----------------
           A committed increment while armed IS one accepted continuation
           step: advance lambda_c and delta_c, count the step, re-establish
           dtime = kappa*ds for the next one, and measure physical progress
           P1/P2/P3.  delta_c itself is re-anchored by the corrector at the
           first iteration of the next step (ct.newstep), so only the
           bookkeeping is done here. */
        if(ct.on==1){
          ITG cp1=0,cp2=0,cp3=0,cq;
          double cpd0;
          ct.step++;ct.ncommit++;
          ct.lamc=ct.lam;
          ct.newstep=1;
          for(cti=0;cti<ne0;cti++){
            if(ipkon[cti]<0){cp2++;continue;}
            if(lakon[8*cti]!='U') continue;
            if(ielprop[cti]<0) continue;
            cpd0=prop[ielprop[cti]+1];
            if(cpd0<=0.) continue;
            cpd0=cpd0/prop[ielprop[cti]];       /* d0 = Tn0/Kn */
            for(ctj=0;ctj<3;ctj++){
              if(ctj>=mi[0]) break;
              cq=mi[0]*cti+ctj;
              if(xstate[*nstate_*cq+3]>=0.5) cp1++;
              if(xstate[*nstate_*cq]>cpd0) cp3++;
            }
          }
          if((cp1>ct.p1r)||(cp2>ct.p2r)||
             (cp3>=ct.p3r+5)) ct.nprog=0;
          else ct.nprog++;
          printf("[DAMAGE CT] COMMIT step=%" ITGFORMAT " inc=%" ITGFORMAT
                 " lambda=%.12e theta=%.12e ds=%.6e dtime=%.6e ; "
                 "P1 failed facets %" ITGFORMAT " (was %" ITGFORMAT
                 "), P2 deleted elements %" ITGFORMAT " (was %" ITGFORMAT
                 "), P3 ip with dmax>d0 %" ITGFORMAT " (was %" ITGFORMAT
                 ") ; quiet windows %" ITGFORMAT "/10 ; evals %" ITGFORMAT
                 " fact %" ITGFORMAT "%s",
                 ct.step,iinc,ct.lam,theta,ct.ds,dtime,
                 cp1,ct.p1r,cp2,ct.p2r,cp3,ct.p3r,
                 ct.nprog,ct.neval,ct.nfact,"\n");
          fflush(stdout);
          ct.p1r=cp1;ct.p2r=cp2;ct.p3r=cp3;
          if((ct.nprog>=10)||
             (ct.step>=ct.maxstep)||
             (ct.neval>=ct.maxeval)||
             (ct.nfact>=ct.maxfact)){
            printf("[DAMAGE CT] segment ENDS: %s.  Requesting the clean "
                   "PARTIAL exit.%s",
                   (ct.nprog>=10)?
                   "no physical progress (P1=P2=P3 unchanged) over 10 "
                   "committed continuation steps":
                   ((ct.step>=ct.maxstep)?"step budget":
                    ((ct.neval>=ct.maxeval)?
                     "residual-evaluation budget":"factorisation budget")),
                   "\n");
            fflush(stdout);
            ct.partial=1;ct.on=0;
          }
        }
      }

      /* previous increment converged: update the initial values */
	  
      iinc++;
      jprint++;
      damage_active_pass=0;
      damage_soft_reeq=0;
      damage_fast_retry=0;
      damage_fast_used=0;
      damage_fast_recover=0;
      rsc.reeq_uam_ref[0]=0.;
      rsc.reeq_uam_ref[1]=0.;

      /* A new physical increment starts with no pending DE1.3 terminal
         triggers.  Entries are populated only at the exact topology-change
         event and are preserved across same-load active-set extensions. */
      if(damage_de13_trigger_value!=NULL){
        for(i=0;i<ne0;i++){
          damage_de13_trigger_value[i]=-1.;
          damage_de13_trigger_ip[i]=0;
        }
      }

      /* Establish the end point of the current normal CCX increment.
         If this increment later needs a cutback, theta_goal is kept
         fixed while reduced physical substeps traverse the same interval.
         Do not move the goal while local substepping is already active. */

      if((ilocalsubstep==0)&&(*ndmat_>0)&&(*iexpl<=1)&&
         (*nmethod!=4)&&(*idrct==0)){
        theta_local_start=theta;
        theta_goal=theta+dtheta;
        if(theta_goal>1.) theta_goal=1.;
        dtheta_restore=dtheta;
      }

      /* store number of elements (important for implicit dynamic
	 contact */

      neini=*ne;
	  
      /* ---- CCX_PATHFOLLOW: commit the increment just accepted --------
         This runs exactly once per accepted physical increment, and it
         must run BEFORE vini is overwritten with vold: the increment's
         displacement change is vold-vini, and after the copy that
         difference is identically zero.  A rejected attempt never reaches
         this block, which is what makes the constraint transactional. */

      if((pf.on==1)&&(pf.pending==1)&&(pathfollow_have()==1)){
        pf.pdu=pathfollow_project(pathfollow_fhat(),vold,vini,nactdof,*nk,mt);
        if(pf.codmode==1){
          pathfollow_cod_settarget(pathfollow_cod_target()+pf.dphi);
          pf.dgc=pf.cu;
        }else if(pf.codmode==2){

          /* CONSTRAINT RESIDUAL AT THE ACCEPTED STATE.  Everything else
             the driver prints is measured mid-iteration, before the last
             correction was applied, so it cannot answer "did the
             constraint converge".  This does: vold is the accepted
             state, vini is still the committed one, and the projection
             uses the same c the corrector used. */

          pf.phiacc=pathfollow_project(pathfollow_cod_c(),vold,vini,nactdof,*nk,mt);
          pf.gacc=pf.phiacc-pf.dphicur;
        }
        pathfollow_commit(pf.lam,pf.pdu,&pf.dgc);
        if(pf.engaged==1){
          pf.taucur*=1.4;
          if(pf.taucur>pf.tauv) pf.taucur=pf.tauv;
          pathfollow_settau(pf.taucur);
          if(pf.codmode==2){

            /* 1.4 was measured to be too greedy on the target: dphi grew
               past what the increment could take, attempt 1 failed, dphi
               was halved, attempt 2 converged - one wasted attempt per
               increment, every increment.  Grow gently instead. */

            pf.dphicur*=pf.ccgrow;
            if(pf.dphicur>pf.dphi) pf.dphicur=pf.dphi;
          }
        }
        pf.dlampred=pf.lam-pf.lamprev;
        pf.lamprev=pf.lam;
        if(!(fabs(pf.dlampred)>1.e-30)) pf.dlampred=dtheta;
        if((pf.engaged==0)&&(pf.dgc>0.2*pf.tauv)){
          pf.engaged=1;
          printf("[PATHFOLLOW] engaged at inc=%" ITGFORMAT " lambda=%.8f: "
                 "measured dG=%.6e has reached 0.2*tau=%.6e\n",
                 iinc,pf.lam,pf.dgc,0.2*pf.tauv);
        }
        if(pf.codmode==2){

          /* Census at the ACCEPTED state (vold, xstate), which is what
             "the front advanced" has to be measured on.  pf.cvec is
             refrozen from the committed state at the top of the next
             attempt, so using it as the scratch vector here costs
             nothing. */

          crackcontrol_build(pf.cvec,neq[1],pf.ccmode,co,kon,ipkon,lakon,
                             *ne,ielprop,prop,xstate,*nstate_,mi,vold,
                             nactdof,*nk,mt,&pf.cs);
          printf("[CRACKCTL] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "dphi=%.6e achieved=%.6e g=%.3e |R|=%.3e du=%.3e "
                 "zone=%" ITGFORMAT " load=%" ITGFORMAT " init=%"
                 ITGFORMAT " fail=%" ITGFORMAT " dead=%" ITGFORMAT
                 " deffmax=%.6e shear=%.4f w=%.6e fhat=%.6f/%.4f\n",
                 iinc,pf.lam,pf.dphicur,pf.phiacc,pf.gacc,ram[0],ram[1],
                 pf.cs.nzone,pf.cs.nload,pf.cs.ninit,pf.cs.nfail,
                 pf.cs.ndead,pf.cs.deffmax,pf.cs.shearfrac,pf.cs.weight,
                 pf.fhcos,pf.fhrat);
        }else if(pf.codmode==1){
          printf("[PATHFOLLOW] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "phi=%.6e target=%.6e\n",iinc,pf.lam,pf.cu,
                 pathfollow_cod_target());
        }else{
          printf("[PATHFOLLOW] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "dG=%.6e P=%.6e ff=%.6e engaged=%" ITGFORMAT " refusals=%"
                 ITGFORMAT "\n",iinc,pf.lam,pf.dgc,pathfollow_Pn(),
                 pathfollow_ff(),pf.engaged,pathfollow_refusals());
        }
        fflush(stdout);
      }
      pf.pending=pf.on;

      /* vold is copied into vini */
	  
      isiz=mt**nk;cpypardou(vini,vold,&isiz,&num_cpus);
	  
      isiz=*nboun;cpypardou(xbounini,xbounact,&isiz,&num_cpus);


      if((*ithermal==1)||(*ithermal>=3)){
	isiz=*nk;cpypardou(t1ini,t1act,&isiz,&num_cpus);
      }
      isiz=neq[1];cpypardou(fini,f,&isiz,&num_cpus);
      if(*nmethod==4){
	if(*iexpl<=1){
	  isiz=mt**nk;
	  cpypardou(veini,veold,&isiz,&num_cpus);
	  cpypardou(accini,accold,&isiz,&num_cpus);
	}
	isiz=mt**nk;cpypardou(fnextini,fnext,&isiz,&num_cpus);

	isiz=neq[1];
	cpypardou(fextini,fext,&isiz,&num_cpus);
	cpypardou(cvini,cv,&isiz,&num_cpus);
	      
	if(*ithermal<2){
	  allwkini=allwk;
	  // MPADD start
	  if(idamping==1)dampwkini = dampwk;
	  for(k=0;k<4;k++){
	    energyini[k]=energy[k];
	  }
	  // MPADD end
	}
      }
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;	    
	cpypardou(stiini,sti,&isiz,&num_cpus);
	cpypardou(emeini,eme,&isiz,&num_cpus);
      }

      /* the contact friction energy is stored at the slave nodes for
         mortar not equal to 1 */
      
      if(*nener==1){
	  isiz=2*mi[0]*ne0;
	cpypardou(enerini,ener,&isiz,&num_cpus);
      }
	      

      if(*mortar!=1){
	if(*nstate_!=0){
	  isiz=*nstate_*mi[0]*(ne0+*nslavs);
	  cpypardou(xstateini,xstate,&isiz,&num_cpus);
	}
      }

      /* store topology and damage at the start of the physical increment */

      if((*ndmat_>0)&&(*iexpl<=1)){
	isiz=ne0;
	cpyparitg(ipkondamageini,ipkon,&isiz,&num_cpus);
	isiz=mi[0]*ne0;
	cpypardou(damdamageini,dam,&isiz,&num_cpus);
	if(damage_damviscini!=NULL){
	  cpypardou(damage_damviscini,damage_damvisc,&isiz,&num_cpus);
	}
	if(*nmethod!=4){
	  isiz=mt**nk;
	  cpypardou(veolddamageini,veold,&isiz,&num_cpus);
	}
      }
	
      /* From this point on every results() call in the physical increment
         reconstructs DE1.2 trial damage from the immutable committed
         baseline.  The public results() ABI remains unchanged. */
      results_set_de12_context(damage_de12_enabled,damage_tangent_mode,
                               ndmat_,ndmcon,dmcon,damdamageini,
                               damage_damjac,damage_damvisc,
                               damage_damviscini,damage_visc_eta);

      if(*mortar>1){
	for (i=0;i<*ntie;i++){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    islavactini[j]=islavact[j];
	    bpini[j]=bp[j];
	    for(k=0;k<mt;k++){
	      cstressini[mt*j+k]=cstress[mt*j+k];
	    }		      
	  }    
	}
      }
    }
      
    /* check for max. # of increments */
      
    if(iinc>jmax[0]){
      printf(" *ERROR in nonlingeo: max. # of increments reached\n\n");
      FORTRAN(stop,());
    }

    if(*iexpl<=1){
      /* [GLOBALIZE] the previous attempt ends where the next is
         announced.  An attempt in which nothing was solved is not counted,
         so a cutback that never reaches a solve cannot dilute the
         denominator. */
      glob_attempt_end(&damage_glob);
      printf(" increment %" ITGFORMAT " attempt %" ITGFORMAT " \n",iinc,icutb+1);
      printf(" increment size= %e\n",dtheta**tper);
      printf(" sum of previous increments=%e\n",theta**tper);
      printf(" actual step time=%e\n",(theta+dtheta)**tper);
      printf(" actual total time=%e\n\n",*ttime+(theta+dtheta)**tper);
      
      printf(" iteration 1\n\n");
    }
      
    qamold[0]=qam[0];
    qamold[1]=qam[1];

    icntrl=0;

    /* restoring the distributed loading before adding the
       friction heating */

    if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
      *nload=nloadref;
      isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
      if(*nam>0){
	isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
      }
      memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
    }
      
    /* store the load level in case damage requires re-equilibration */

    if(idamagereeq==0){
      thetadamage=theta;
      dthetadamage=dtheta;
      dthetarefdamage=dthetaref;
    }

    /* determining the actual loads at the end of the new increment*/

    /* Once the constraint drives lambda, theta is only a monotone counter
       that keeps dtime positive and paces the output.  The stock
       controller still sizes it, but it is capped here: theta reaching 1
       ends the step, and an uncapped theta would end the run in the middle
       of the snap-back while lambda still has most of the branch to go. */

    /* CAP, not set.  The original code SET dtheta to a constant once the
       constraint took over, on the argument that theta is only a counter
       when lambda is decoupled from it, so shrinking it cannot make a
       retry easier.  That argument holds only for a rate-INDEPENDENT
       model.

       The s3rad target is not one: CCX_DAMAGE_VISCOSITY=1e-4 and the UC6
       law both relax with alpha = dt/(mu+dt), so dtime is a constitutive
       parameter and a theta cutback really does change the problem.
       Measured: the stock run needs dtime around 1e-6 at increment 250
       (alpha=0.0095) and recovers there by cutting it; pinning dtheta at
       1e-4 holds alpha at 0.5, i.e. removes almost all of the viscous
       regularisation, and the first engaged attempt failed through five
       cutbacks that could not touch dtime.

       Capping it and leaving the stock controller to size it below the
       cap was tried and MEASURED, and it does not work either, for a
       different reason: once engaged, an attempt that fails is fixed by
       halving dphi, not by halving dtheta, but the stock controller
       halves dtheta anyway and then refuses to grow it back because the
       accepted attempt took 29 iterations.  dtheta therefore ratchets
       down one notch per increment and the run dies on "increment size
       smaller than minimum" - measured, 14 accepted engaged increments
       and then that exact stop.

       So dtheta is SET again, and it is a CHOICE of time step: with
       lambda decoupled, theta is the clock the rate-dependent laws run
       on and nothing else, and holding it fixed makes the viscous
       relaxation per increment constant.  The continuation step size is
       dphi and that alone is what a cutback shrinks. */

    if((pf.on==1)&&(pf.engaged==1)){
      dtheta=pf.dtheta_eng;
      if(theta+dtheta>1.) dtheta=1.-theta;
    }

    reltime=theta+dtheta;

    /* Content hash of the COMMITTED state at the top of this attempt.
       vini/xstateini/damdamageini/ipkondamageini are exactly what a
       rollback restores, so two attempts that print the same hash really
       did start from the same state - which is the half of the
       "repeating event" claim that cannot be read off the existing log. */

    if((td_armed!=0)&&(td_trace!=0)){
      td_state=topodiag_hash_seed();
      td_state=topodiag_hash_d(vini,mt**nk,td_state);
      if(xstateini!=NULL)
        td_state=topodiag_hash_d(xstateini,*nstate_*mi[0]**ne,td_state);
      if(damdamageini!=NULL)
        td_state=topodiag_hash_d(damdamageini,mi[0]*ne0,td_state);
      if(ipkondamageini!=NULL)
        td_state=topodiag_hash_i(ipkondamageini,ne0,td_state);
      td_state=topodiag_hash_i(ipkon,*ne,td_state);
      printf("[BATCHTRACE] attempt inc=%" ITGFORMAT " icutb=%" ITGFORMAT
             " theta=%.12e dtheta=%.12e committed_state=%016llx\n",
             iinc,icutb,theta,dtheta,td_state);
      fflush(stdout);
    }

    /* [WALLDIAG] theta is the increment's starting load factor here, before
       dtheta is added, so the gate opens on the first attempt whose base
       state is at or beyond the recorded wall. */

    if((prb.wall_theta>=0.)&&(theta>=prb.wall_theta)){
      if(prb.wall_armed==0){
        prb.wall_armed=1;
        printf("[WALLDIAG] GATE OPEN at inc=%" ITGFORMAT " icutb=%" ITGFORMAT
               " theta=%.12e dtheta=%.12e%s",iinc,icutb,theta,dtheta,"\n");
        fflush(stdout);
      }
      dog.lincheck=iinc;
      rsc.ls_probe=1;
      if(td_from<=0) td_from=iinc;

      /* NOT the null-vector probe.  CCX_DAMAGE_NULLVEC ends with
         stopwithout201() by design - it is a terminal measurement, taken
         once, on a run that is not meant to continue.  Arming it from a
         load-factor gate killed the run at the gate instead of at the wall
         - measured, R0-early stopped at increment 92 with the ladder
         unmeasured.  It stays opt-in and its terminal nature is stated. */
      if((prb.wall_null!=0)&&(prb.null_inc<=0)){
        prb.null_inc=iinc;
        printf("[WALLDIAG] arming the deflated null-vector probe on this "
               "increment.  IT IS TERMINAL: the run stops after it "
               "prints.%s","\n");
        fflush(stdout);
      }
    }

    FORTRAN(uc6setinc,(&iinc));

    /* trial load factor for dissipation control; theta itself is
       left alone until the increment is accepted */

    if(lc.diss_ctrl>=1){
      lc.diss_lamcur=(lc.arc==1)?(lc.arc_lam+dtheta):(theta+dtheta);
      lc.diss_have=0;
    }
    time=reltime**tper;
    dtime=dtheta**tper;
      
    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
		      xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,ampli,
		      &time,&reltime,ttime,&dtime,ithermal,nmethod,xbounold,
		      xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,nodeforc,
		      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,
		      ilboun,nelemload,sideload,mi,ntrans,trab,inotr,veold,
		      integerglob,doubleglob,tieset,istartset,iendset,ialset,
		      ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,ipobody,
		      iponoeln,inoeln,ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,
		      rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));

    /* tempload has just built xbounact from reltime=theta+dtheta.  With
       lambda decoupled that is no longer the load factor, so rebuild the
       prescribed values from lambda itself.  xbounold is the value at the
       START of the step, so this is an absolute step fraction - exactly the
       ramp CCX_DAMAGE_PATH mode 1 verified against tempload at max|diff|=0. */

    if(lc.arc==1){
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*lc.diss_lamcur;
      }
    }

    /* ---- CCX_PATHFOLLOW: predictor for this attempt -----------------
       pathfollow_incstart discards whatever a rejected attempt left, so
       a cutback re-enters here from the committed state.  Before the
       constraint engages lambda simply follows the step time, which is
       the ordinary control the stock solver would have run. */

    if(pf.on==1){
      if(neq[1]!=pf.neqarm){
        if(pathfollow_resize(neq[1])==0){
          printf("[PATHFOLLOW] disarmed: could not follow the equation "
                 "count from %" ITGFORMAT " to %" ITGFORMAT "\n",
                 pf.neqarm,neq[1]);
          pf.on=0;
        }else{
          RENEW(pf.uf,double,neq[1]);
          printf("[PATHFOLLOW] equation count %" ITGFORMAT " -> %" ITGFORMAT
                 "; constraint origin restarted at the current state\n",
                 pf.neqarm,neq[1]);
          pf.neqarm=neq[1];

          /* The control functional lives in EQUATION space.  After a
             remastruct its entries refer to dofs that no longer exist,
             so a vector kept from before the change is not merely stale,
             it is read out of bounds.  The mixed-mode functional is
             rebuilt from the model a few lines below, every attempt; the
             absolute Mode-I functional cannot be rebuilt without redoing
             the geometry pass, so it is disarmed loudly rather than used
             wrong. */

          if(pf.codmode==2){
            RENEW(pf.cvec,double,neq[1]);
          }else if(pf.codmode==1){
            printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW_COD builds its "
                   "control functional once, from a fixed equation "
                   "numbering.  The equation count has changed, so the "
                   "vector no longer refers to the same dofs.  "
                   "Disarming rather than continuing with a stale "
                   "constraint; use CCX_CRACK_CONTROL, which refreezes "
                   "the functional every attempt.\n");
            pf.on=0;pf.engaged=0;pf.codmode=0;
          }
        }
      }
    }
    if(pf.on==1){

      /* The dissipation increment IS the step size here.  A stock cutback
         shrinks theta, but lambda no longer follows theta, so the retry
         would be the identical problem.  Shrink tau instead - that is the
         path-following cutback, and without it an engaged increment can
         only fail down to the minimum step time. */

      if((pf.engaged==1)&&(icutb>pf.icutbprev)){
        pf.taucur*=0.5;
        if(pf.taucur<1.e-6*pf.tauv) pf.taucur=1.e-6*pf.tauv;
        pf.ncut++;
        pathfollow_settau(pf.taucur);

        /* For crack control the CONTROL INCREMENT is the step size, for
           exactly the same reason: lambda no longer follows theta, so a
           theta cutback retries the identical problem. */

        if(pf.codmode==2){
          pf.dphicur*=0.5;
          if(pf.dphicur<1.e-8*pf.dphi){
            pf.dphicur=1.e-8*pf.dphi;

            /* TERMINATION.  With dtheta held fixed the stock stop
               "increment size smaller than minimum" can never fire, so
               without this the run has no termination criterion at all
               and a wall becomes an infinite loop instead of a reported
               result.  Measured on the target: the deletion batch at
               increment 348 leaves an orphan node, its same-load
               re-equilibration does not converge, the batch is rolled
               back, and the identical event repeats forever.

               dphi at its floor and still failing means the obstruction
               is not the size of the continuation step. */

            pf.ncutfloor++;
            if(pf.ncutfloor>=3){
              printf("[CRACKCTL] *WALL: the control increment is at its "
                     "floor (%.6e, 1e-8 of the requested %.6e) and the "
                     "attempt still fails.  The obstruction is not the "
                     "size of the continuation step.  Stopping.\n",
                     pf.dphicur,pf.dphi);
              fflush(stdout);
              FORTRAN(stop,());
            }
          }else{
            pf.ncutfloor=0;
          }
          printf("[CRACKCTL] cutback: dphi -> %.6e\n",pf.dphicur);
        }

        /* the predictor is proportional to tau, so it rescales itself */

        printf("[PATHFOLLOW] cutback %" ITGFORMAT ": tau -> %.6e\n",
               pf.ncut,pf.taucur);
        fflush(stdout);
      }
      pf.icutbprev=icutb;
      pathfollow_incstart();

      /* ---- mixed-mode crack control: refreeze the functional --------
         Built from the COMMITTED state (vini, xstateini), which is what
         a rollback restores, so a retry rebuilds the identical vector.
         The constraint is incremental, so the target is simply the
         control increment currently in force. */

      if(pf.codmode==2){
        pf.ccnw=crackcontrol_build(pf.cvec,neq[1],pf.ccmode,co,kon,ipkon,
                                   lakon,*ne,ielprop,prop,xstateini,
                                   *nstate_,mi,vini,nactdof,*nk,mt,&pf.cs);

        /* ZONE and DISS are supported on the process zone.  Before the
           first initiation, and again if every initiated point has
           failed, that set is empty and the functional carries no
           weight.  Fall back to the mean over all live facets for this
           attempt rather than handing the bordered row a zero vector. */

        if((pf.ccnw==0)&&(pf.ccmode!=0)){
          pf.ccnw=crackcontrol_build(pf.cvec,neq[1],0,co,kon,ipkon,
                                     lakon,*ne,ielprop,prop,xstateini,
                                     *nstate_,mi,vini,nactdof,*nk,mt,
                                     &pf.cs);
        }
        pathfollow_cod_arm(pf.cvec,neq[1]);
        pathfollow_cod_settarget(pf.dphicur);

        /* With the extrapolation suppressed, b/dlamjump at the first
           iteration IS -dR/dlambda, so there is no reason to keep a
           vector captured many increments and one remastruct ago.  The
           freeze existed to stop the dissipation constraint rescaling
           its own tau; the crack-control constraint does not depend on
           the scale of f_hat, only on its being the true derivative. */

        if(pf.engaged==1) pathfollow_unfreeze();

        if((pf.engaged==0)&&(pf.ccnw>0)){
          if(((pf.ccengage>0)&&(iinc>=pf.ccengage))||
             ((pf.ccengage<=0)&&(pf.cs.nzone>0))){
            pf.engaged=1;

            /* The cap above was applied before this decision was taken,
               so the attempt that engages would otherwise run with the
               uncapped stock dtheta - which is exactly the attempt where
               the constitutive relaxation must not jump.  Apply it here
               too.  Only the time-like quantities are rebuilt: xbounact
               is overwritten from pf.lam a few lines below anyway. */

            dtheta=pf.dtheta_eng;
            if(theta+dtheta>1.) dtheta=1.-theta;
            reltime=theta+dtheta;
            time=reltime**tper;
            dtime=dtheta**tper;
            printf("[CRACKCTL] engaged at inc=%" ITGFORMAT " lambda=%.8f: "
                   "process zone %" ITGFORMAT " ip(s) of which %"
                   ITGFORMAT " loading, initiated %" ITGFORMAT
                   ", failed %" ITGFORMAT ", shear fraction %.4f, "
                   "weight %.6e\n",iinc,pathfollow_lamn(),pf.cs.nzone,
                   pf.cs.nload,pf.cs.ninit,pf.cs.nfail,pf.cs.shearfrac,
                   pf.cs.weight);
            fflush(stdout);
          }
        }
      }

      if(pf.engaged==0){
        pf.lam=theta+dtheta;
      }else{

        /* Tangent predictor, equation (7) in pathfollow.c.  Its SIGN is
           not chosen here: it comes out of lambda_n*ff-P_n, which changes
           sign at the limit point.  That is what turns the branch. */

        if(pf.codmode==0){
          if(pathfollow_predictor(&pf.dlampred)==0){
            printf("[PATHFOLLOW] predictor degenerate at inc=%" ITGFORMAT
                   "; keeping dlambda=%+.4e\n",iinc,pf.dlampred);
          }
        }
        /* For crack-opening control the previous accepted dlambda is the
           predictor: it carries the sign of the branch, and the constraint
           corrects the length in one step because phi is linear in u.

           It must not be zero on the first engaged increment: f_hat is
           captured by dividing the first residual by this jump, so a zero
           predictor leaves f_hat undefined, the corrector never runs and
           lambda sits still while the step controller burns its cutbacks.
           Measured as "too many cutbacks" with no [PATHFOLLOW] it= line at
           all. */
        /* PROBE JUMP, not predictor jump, once the crack constraint is
           driving.

           f_hat is obtained as b/dlamjump at the first iteration.  With
           the displacement extrapolation suppressed that quotient is a
           SECANT of -dR/dlambda over the jump, and it is the tangent only
           in the limit of a small jump.  Using the previous accepted
           dlambda as that jump makes it a secant over the whole physical
           step: at increment 60 of the target the step is 2e-03, over
           which the boundary-adjacent material yields and damages, and
           an honest finite difference of the residual put
           |f_hat|/|dR/dlambda| at 1.43.

           A scale error there is not harmless.  The corrector adds
           dlambda*K^-1 f_hat to u and dlambda to lambda; if f_hat = s*q
           the DISPLACEMENT half is right and the LOAD half is a factor
           1/s short, so the prescribed dofs end up where the free dofs
           did not assume.  Measured on the target: the constraint row
           converged to 1e-19 every iteration while the equilibrium
           residual stalled at 0.17 - and it stalled on nodes next to the
           loaded face, which is exactly where -dR/dlambda lives.  With
           dlambda clamped to 1e-9 the same increment converged to 8e-04
           and was accepted, so the machinery is sound and the reference
           vector is what is wrong.

           So take the jump as a small PROBE of size pf.eps.  It costs
           nothing - no extra residual evaluation, no state to restore -
           because the first iteration of the increment has to evaluate a
           residual anyway.  The constraint then supplies the real
           predictor: its own dlambda closes the whole step in one
           corrector, since phi is linear in u. */

        if(pf.codmode==2) pf.dlampred=pf.eps;
        if(!(fabs(pf.dlampred)>1.e-12)) pf.dlampred=dtheta;
        pf.lam=pathfollow_lamn()+pf.dlampred;
      }
      pf.dlamjump=pf.lam-pathfollow_lamn();
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*pf.lam;
      }

      /* Bisection probe: is the model state at the START of an attempt
         really the committed one?  vold-vini must be identically zero
         here; anything else means the rollback did not restore what the
         constraint assumes it restored. */

      if(ccxopt_getenv("CCX_PATHFOLLOW_ACCUMCHECK")!=NULL){
        double pfd=0.,pfm=0.;
        ITG pfi,pfj,pfk;
        for(pfi=0;pfi<*nk;pfi++){
          for(pfj=1;pfj<mt;pfj++){
            pfk=nactdof[mt*pfi+pfj];
            if(pfk>0){
              pfd+=(vold[mt*pfi+pfj]-vini[mt*pfi+pfj])*
                   (vold[mt*pfi+pfj]-vini[mt*pfi+pfj]);
              if(fabs(vold[mt*pfi+pfj]-vini[mt*pfi+pfj])>pfm)
                pfm=fabs(vold[mt*pfi+pfj]-vini[mt*pfi+pfj]);
            }
          }
        }
        printf("[PF-START] inc=%" ITGFORMAT " attempt icutb=%" ITGFORMAT
               " |vold-vini|=%.6e max=%.6e lambda=%.8f\n",
               iinc,icutb,sqrt(pfd),pfm,pf.lam);
        fflush(stdout);
      }
    }


    /* ================= [DAMAGE CT] ARMING =========================
       Runs once, at the top of the attempt that follows the level-4
       firing.  The standard cutback rollback has already restored the
       increment-start state, so every quantity below is read from the last
       COMMITTED state.  ANY refusal leaves that state untouched, sets
       ct.refused and hands the wall to the ORIGINAL stock stop. */

    if(ct.arm==1){
      ITG ai,aj,ak,an,anp,anm,abest=-1,aip=-1,acnt=0,aok,ansel=0;
      ITG afl0,afl1,asupp[18],anod[18],adir[18];
      double aw[18],adl[3],armat[9],ash[3],ad0[3],ad1[3];
      double am[3],adeff,abeta,atn0,ats0,agc,atau,adf;
      double aq[5],akp[5],atmp,abestv=-1.,aclam,ag[3];
      ITG ah,ahp,as0,as1,ne_uc6=0;

      ct.arm=0;
      ct.refused=1;               /* pessimistic until every guard passes */

      /* ---- admissibility domain (SPEC FREEZE v1) ---- */
      aok=1;
      if(*nam!=0) aok=0;  if(*nforc!=0) aok=0;  if(*nload!=0) aok=0;
      if(*nbody!=0) aok=0; if(*nmpc!=0) aok=0;  if(ncont!=0) aok=0;
      if(*mortar>1) aok=0; if(*ithermal>=2) aok=0; if(*nmethod!=1) aok=0;
      if(*idrct!=0) aok=0; if(*iprestr!=0) aok=0; if(*nstate_<4) aok=0;
      if(*isolver!=7) aok=0;
      {
        char *ae=ccxopt_getenv("CCX_PARDISO_REUSE_SYMBOLIC");
        if((ae==NULL)||((strcmp(ae,"1")!=0)&&(strcmp(ae,"ON")!=0)&&
                        (strcmp(ae,"on")!=0)&&(strcmp(ae,"YES")!=0)&&
                        (strcmp(ae,"yes")!=0))) aok=0;
      }
      for(ai=0;ai<ne0;ai++)
        if((ipkon[ai]>=0)&&(lakon[8*ai]=='U')) ne_uc6++;
      if(ne_uc6==0) aok=0;
      if(ct.nring<6) aok=0;
      if(aok==0){
        printf("[DAMAGE CT] REFUSING TO ARM: admissibility domain or ring not "
               "satisfied (nam=%" ITGFORMAT " nforc=%" ITGFORMAT " nload=%"
               ITGFORMAT " nbody=%" ITGFORMAT " nmpc=%" ITGFORMAT " ncont=%"
               ITGFORMAT " mortar=%" ITGFORMAT " ithermal=%" ITGFORMAT
               " nmethod=%" ITGFORMAT " idrct=%" ITGFORMAT " iprestr=%"
               ITGFORMAT " nstate_=%" ITGFORMAT " isolver=%" ITGFORMAT
               " uc6_elements=%" ITGFORMAT " ring=%" ITGFORMAT "/6, "
               "PARDISO symbolic reuse required).  The last committed state "
               "is untouched; this wall goes to the ORIGINAL stock stop.%s",
               *nam,*nforc,*nload,*nbody,*nmpc,ncont,*mortar,*ithermal,
               *nmethod,*idrct,*iprestr,*nstate_,*isolver,ne_uc6,
               ct.nring,"\n");
        fflush(stdout);
      }else{

      /* ---- candidate scan over the five committed intervals ---- */
      ah=ct.head;
      an=mi[0]*ne0;
      printf("[DAMAGE CT] arming at inc=%" ITGFORMAT " theta=%.12e: scanning "
             "%" ITGFORMAT " UC6 elements over 6 committed endpoints.%s",
             iinc,theta,ne_uc6,"\n");
      for(ai=0;ai<ne0;ai++){
        if(ipkon[ai]<0) continue;
        if(lakon[8*ai]!='U') continue;
        if(ielprop[ai]<0) continue;
        atn0=prop[ielprop[ai]+1];
        ats0=prop[ielprop[ai]+2];
        agc =prop[ielprop[ai]+3];
        if((atn0<=0.)||(ats0<=0.)||(agc<=0.)) continue;
        abeta=(ats0/atn0)*(ats0/atn0);
        adf=2.*agc/atn0;
        atau=1.e-9*adf;
        for(aj=0;aj<3;aj++){
          ak=an*ah+ (mi[0]*ai+aj);
          afl0=ct.fl[an*ah+mi[0]*ai+aj];
          if(afl0<0) continue;                   /* not a live UC6 ip      */
          if(afl0&1) continue;                   /* already failed         */
          /* frozen mixed-mode direction m at the newest committed endpoint */
          ad1[0]=ct.ring[3*(an*ah+mi[0]*ai+aj)];
          ad1[1]=ct.ring[3*(an*ah+mi[0]*ai+aj)+1];
          ad1[2]=ct.ring[3*(an*ah+mi[0]*ai+aj)+2];
          adeff=(ad1[0]>0.?ad1[0]*ad1[0]:0.)
                +abeta*(ad1[1]*ad1[1]+ad1[2]*ad1[2]);
          adeff=(adeff>0.)?sqrt(adeff):0.;
          if(!(adeff>0.)) continue;
          if(fabs(ad1[0])<0.05*adeff) continue;  /* too close to the kink  */
          am[0]=(ad1[0]>0.?ad1[0]:0.)/adeff;
          am[1]=abeta*ad1[1]/adeff;
          am[2]=abeta*ad1[2]/adeff;
          /* five historical intervals, ALL must be admissible */
          aok=1;
          for(ak=0;ak<5;ak++){
            as1=(ah-ak+6)%6;  as0=(ah-ak-1+6)%6;
            afl1=ct.fl[an*as1+mi[0]*ai+aj];
            afl0=ct.fl[an*as0+mi[0]*ai+aj];
            if((afl1<0)||(afl0<0)){aok=0;break;}
            if((afl1&1)||(afl0&1)){aok=0;break;}        /* alive both ends */
            if((afl1&2)!=(afl0&2)){aok=0;break;}        /* ucomp unchanged */
            for(anp=0;anp<3;anp++){
              ad1[anp]=ct.ring[3*(an*as1+mi[0]*ai+aj)+anp];
              ad0[anp]=ct.ring[3*(an*as0+mi[0]*ai+aj)+anp];
            }
            if((ad1[0]>0.)!=(ad0[0]>0.)){aok=0;break;}  /* sign(dn) kept   */
            atmp=am[0]*(ad1[0]-ad0[0])+am[1]*(ad1[1]-ad0[1])
                +am[2]*(ad1[2]-ad0[2]);
            if(!(atmp>atau)){aok=0;break;}
            aq[ak]=atmp;
            if(!(ct.dt[as1]>0.)){aok=0;break;}
            akp[ak]=ct.dt[as1]/atmp;
            if(!(akp[ak]>0.)){aok=0;break;}
          }
          if(aok==0) continue;
          acnt++;
          if(aq[0]>abestv){
            abestv=aq[0];abest=ai;aip=aj;
            for(anp=0;anp<3;anp++) ct.m[anp]=am[anp];
            for(anp=0;anp<5;anp++){
              ct.dt[anp]=ct.dt[anp];      /* untouched      */
            }
            /* keep this candidate's five ratios and advances */
            for(anp=0;anp<5;anp++){ct.g[0]=0.;}
            for(anp=0;anp<5;anp++){aq[anp]=aq[anp];akp[anp]=akp[anp];}
            /* median of five: insertion sort on local copies */
            {
              double sq[5],sk[5];
              ITG bi,bj2;
              for(bi=0;bi<5;bi++){sq[bi]=aq[bi];sk[bi]=akp[bi];}
              for(bi=1;bi<5;bi++){
                atmp=sq[bi];for(bj2=bi;(bj2>0)&&(sq[bj2-1]>atmp);bj2--)
                  sq[bj2]=sq[bj2-1]; sq[bj2]=atmp;
                atmp=sk[bi];for(bj2=bi;(bj2>0)&&(sk[bj2-1]>atmp);bj2--)
                  sk[bj2]=sk[bj2-1]; sk[bj2]=atmp;
              }
              ct.ds0=sq[2];
              ct.kappa=sk[2];
              ct.tau=(sk[0]>0.)?sk[4]/sk[0]:1.e30;
            }
          }
        }
      }
      printf("[DAMAGE CT] candidates with five admissible intervals: %"
             ITGFORMAT "%s",acnt,"\n");
      if((acnt==0)||(abest<0)){
        printf("[DAMAGE CT] REFUSING TO ARM: no UC6 integration point has "
               "five consecutive admissible committed intervals at this "
               "state.  The last committed state is untouched; this wall "
               "goes to the ORIGINAL stock stop.%s","\n");
        fflush(stdout);
      }else if(!(ct.tau<=ct.kaptol)){
        printf("[DAMAGE CT] REFUSING TO ARM: kappa is not stable over the "
               "five intervals of the best candidate (max/min=%.4f > 1.5, "
               "median kappa=%.6e).  Refusing rather than freezing an "
               "unmeasured control-to-time coupling.%s",
               ct.tau,ct.kappa,"\n");
        fflush(stdout);
      }else{
        /* ---- freeze the control point, m, the stencil and c_lambda ---- */
        ct.elem=abest;ct.ip=aip;
        damcont_kin(co,kon,ipkon[abest],vold,mt,aip,adl,armat,ash);
        for(ak=0;ak<3;ak++){
          ag[ak]=ct.m[0]*armat[ak]+ct.m[1]*armat[3+ak]
                +ct.m[2]*armat[6+ak];
          ct.g[ak]=ag[ak];
          /* delta_c is taken from the SAME committed snapshot m came from
             - the newest ring endpoint - not recomputed from vold.  Mixing
             the two sources is what made the first constraint value differ
             from -ds by 83x on the short deck: the ring endpoint and vold
             are not the same instant when the wall is reached in a later
             increment than the last commit. */
          ct.dc[ak]=ct.ring[3*(mi[0]*ne0*ah+mi[0]*abest+aip)+ak];
        }
        printf("[DAMAGE CT] delta_c source check: ring=(%.9e,%.9e,%.9e) "
               "vold=(%.9e,%.9e,%.9e) m.(vold-ring)=%.6e%s",
               ct.dc[0],ct.dc[1],ct.dc[2],
               adl[0],adl[1],adl[2],
               ct.m[0]*(adl[0]-ct.dc[0])
              +ct.m[1]*(adl[1]-ct.dc[1])
              +ct.m[2]*(adl[2]-ct.dc[2]),"\n");
        ct.nsupp=0;aclam=0.;aok=1;
        for(ai=0;ai<3;ai++){
          anm=kon[ipkon[abest]+ai]-1;
          anp=kon[ipkon[abest]+ai+3]-1;
          for(ak=0;ak<3;ak++){
            anod[ct.nsupp]=anp;adir[ct.nsupp]=ak+1;
            aw[ct.nsupp]= ash[ai]*ag[ak];ct.nsupp++;
            anod[ct.nsupp]=anm;adir[ct.nsupp]=ak+1;
            aw[ct.nsupp]=-ash[ai]*ag[ak];ct.nsupp++;
          }
        }
        if(ct.w==NULL) NNEW(ct.w,double,18);
        for(ai=0;ai<ct.nsupp;ai++){
          ak=nactdof[mt*anod[ai]+adir[ai]];
          ct.w[ai]=aw[ai];
          ct.node[ai]=anod[ai];ct.dir[ai]=adir[ai];
          if(ak>0){
            asupp[ai]=ak-1;
          }else{
            asupp[ai]=-1;
            if(ak<0){
              if(ak!=2*(ak/2)) aok=0;          /* odd negative = MPC       */
            }
            /* SPC: find its boundary entry and add w * d(u_c)/d(lambda) */
            for(aj=0;aj<*nboun;aj++){
              if((nodeboun[aj]-1==anod[ai])&&(ndirboun[aj]==adir[ai])){
                aclam+=aw[ai]*(xboun[aj]-xbounold[aj]);
                break;
              }
            }
          }
          ct.supp[ai]=asupp[ai];
        }
        ct.clam=aclam;
        if(aok==0){
          printf("[DAMAGE CT] REFUSING TO ARM: an MPC-dependent degree of "
                 "freedom is in the control stencil; the MPC expansion is "
                 "not implemented.  ORIGINAL stock stop.%s","\n");
          fflush(stdout);
        }else{
          ct.lamc=theta;
          ct.lam=theta;
          ct.ds=ct.ds0;
          ct.dsmin=1.e-3*ct.ds0;
          ct.dsmax=20.*ct.ds0;
          ct.tolc=1.e-3*ct.ds0;
          if(ct.tolc<1.e-9) ct.tolc=1.e-9;
          /* healthy references for the absolute predictor bounds */
          {
            double sl[5],su[5];ITG bi,bj2;
            for(bi=0;bi<5;bi++){
              as1=(ah-bi+6)%6;
              sl[bi]=ct.dlam[as1];su[bi]=ct.duinf[as1];
            }
            for(bi=1;bi<5;bi++){
              atmp=sl[bi];for(bj2=bi;(bj2>0)&&(sl[bj2-1]>atmp);bj2--)
                sl[bj2]=sl[bj2-1]; sl[bj2]=atmp;
              atmp=su[bi];for(bj2=bi;(bj2>0)&&(su[bj2-1]>atmp);bj2--)
                su[bj2]=su[bj2-1]; su[bj2]=atmp;
            }
            ct.lamref=sl[2];ct.duref=su[2];
          }
          ct.on=1;ct.refused=0;ct.step=0;
          ct.it=0;ct.ncommit=0;
          dtheta=ct.kappa*ct.ds/(*tper);
          if(dtheta<=0.) dtheta=ct.dt[ah];
          dthetaref=dtheta;
          printf("[DAMAGE CT] ARMED.  control point: element %" ITGFORMAT
                 " ip %" ITGFORMAT "; m=(%.6f,%.6f,%.6f); |m.delta_c|=%.6e; "
                 "c_lambda=%.6e; support dofs %" ITGFORMAT " of 18 free; "
                 "kappa=%.6e (max/min=%.4f over 5 intervals); ds0=%.6e; "
                 "ds_min=%.6e ds_max=%.6e; tol_c=%.6e; lambda_c=%.12e; "
                 "dtheta=kappa*ds/tper=%.6e; healthy refs dlambda=%.6e "
                 "|du|inf=%.6e.  lambda now OWNS the boundary for the rest "
                 "of the step; every ending is PARTIAL.%s",
                 abest+1,aip+1,ct.m[0],ct.m[1],ct.m[2],
                 ct.m[0]*adl[0]+ct.m[1]*adl[1]
                 +ct.m[2]*adl[2],
                 ct.clam,ct.nsupp,ct.kappa,
                 ct.tau,ct.ds0,ct.dsmin,
                 ct.dsmax,ct.tolc,ct.lamc,dtheta,
                 ct.lamref,ct.duref,"\n");
          fflush(stdout);
        }
      }
      }
      if(ct.on==0) ct.refused=1;
    }

    /* [DAMAGE CT] once armed, lambda owns the boundary for the rest of the
       step.  xbounold is the value at the START of the step, so this is an
       absolute step fraction, and it is the same ramp CCX_DAMAGE_PATH mode 1
       verified against tempload at max|diff|=0.  theta stays monotone
       pseudo-time and carries dtime only. */
    if(ct.on==1){
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*ct.lam;
      }
      /* dtime = kappa*ds is re-established for every armed attempt, so a
         reduced ds reduces the viscous time by the same factor and the
         total viscous time across a given control advance is invariant. */
      if(ct.kappa>0.){
        dtheta=ct.kappa*ct.ds/(*tper);
        if(dtheta>0.) dthetaref=dtheta;
      }
    }

    /* Path control: give the load factor an identity of its own.

       tempload builds xbounact from reltime=theta+dtheta, and theta only
       ever advances, so lambda cannot be reduced.  That is the one
       capability missing when a deletion releases energy and the
       neighbouring equilibrium sits at a LOWER load.  xbounold holds the
       value at the start of the step (it is refreshed only after the
       increment loop), so xbounold+(xboun-xbounold)*lambda is an
       absolute step fraction and lambda is free to move either way.

       Mode 1 only measures: it rebuilds the ramp and reports the largest
       deviation from what tempload produced, which has to be zero before
       anything is allowed to rely on the equivalence.  Amplitudes would
       break it, and this reports that rather than assuming it. */

    if(lc.path_on>0){
      lc.path_lam=theta+dtheta;

      /* Diagnostic descent.  After a deletion transaction has been
         rolled back, the question is whether an equilibrium exists
         at all for the post-deletion topology, or only at a LOWER
         load.  theta cannot answer it: it only advances.  Here the
         retry is placed below the committed load factor, keeping
         theta as a monotone counter so dtime stays positive for the
         viscous update.  If Newton converges there, the stall is a
         limit point and path control is the answer; if it converges
         nowhere, it is not. */

      /* Arm on ANY repeated attempt, not only a rolled-back deletion.
         The first version armed on the deletion rollback alone and
         fired exactly once in a 533-increment run, never reaching the
         stall, which fails as slow convergence with no deletion in
         that increment at all. */
      /* Controlled descent, one accepted increment at a time.

         The first version dropped lambda by up to 80% inside a single
         retry.  That does not test whether an equilibrium exists lower
         down; it asks Newton to jump from a converged field at one load
         to a distant field at another, and it destroyed the state -
         displacement increment 229 in a specimen of size 4, residual
         3.4e8, damage saturated in 9195 elements.  A descent has to be
         walked, with each step accepted, exactly like a loading path. */

      if((lc.path_on>=2)&&(lc.path_desc>0)){
        lc.path_lam=lc.path_lamcom*(1.-lc.path_drop);
        if(lc.path_lam<1.e-6) lc.path_lam=1.e-6;

        /* The descent is armed when dtheta has already collapsed to
           the floor, so without this it fails the minimum-size test
           immediately and never gets to walk - three attempts and
           out.  With lambda decoupled, dtheta only carries time for
           the viscous update; the load step is the lambda drop.
           Give the descent a workable step and a fresh set of
           attempts, bounded so a descent that never converges cannot
           spin. */
        lc.path_att++;
        if(lc.path_att>4*lc.path_nstep){
          lc.path_desc=0;
          printf("[DAMAGE PATH] descent abandoned after %" ITGFORMAT
                 " attempts without an accepted step\n",
                 lc.path_att);
          fflush(stdout);
        }else{
          if(dtheta<0.5*dthetaref) dtheta=0.5*dthetaref;
          icutb=0;
        }
        printf("[DAMAGE PATH] inc=%" ITGFORMAT " descent step %"
               ITGFORMAT " left: lambda %.6f -> %.6f\n",
               iinc,lc.path_desc,lc.path_lamcom,lc.path_lam);
        fflush(stdout);
      }
      lc.path_dev=0.;
      for(k=0;k<*nboun;k++){
        lc.path_ref=xbounold[k]+(xboun[k]-xbounold[k])*lc.path_lam;
        if(fabs(lc.path_ref-xbounact[k])>lc.path_dev)
          lc.path_dev=fabs(lc.path_ref-xbounact[k]);
        if(lc.path_on>=2) xbounact[k]=lc.path_ref;
      }
      if(lc.path_dev>lc.path_devmax){
        lc.path_devmax=lc.path_dev;
        printf("[DAMAGE PATH] inc=%" ITGFORMAT " lambda=%.6f "
               "max|ramp-tempload|=%.6e (new maximum)\n",
               iinc,lc.path_lam,lc.path_dev);
        fflush(stdout);
      }
    }

    for(i=0;i<3;i++){
      cam[i]=0.;}
    for(i=3;i<5;i++){
      cam[i]=0.5;}
    if(*ithermal>1){
      radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
		  ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,
		  vold,
		  shcon,nshcon,ipkon,kon,co,
		  kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
		  nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
		  cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
		  ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
		  nactdog,nacteq,nodeboun,ndirboun,network,
		  rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,
		  ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,
		  ialset,nset,
		  ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,
		  iamload,
		  jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,ncocon,
		  iponoeln,
		  inoeln,nprop,amname,namta,amta,iexpl);
             
      /* check whether network iterations converged */

      if(qa[2]>0){
	checkdivergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
			kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
			xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
			ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
			vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
			&theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
			nactdof,b,tmin,ctrl,amta,namta,itpamp,inext,&dthetaref,
			&itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
			nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
			set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
			mortar,nmat,ielprop,prop,&ialeatoric,&kscale,
			energy, &allwk,&energyref,&emax,&r_abs,&enetoll,
			energyini,
			&allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			&dampwkini,energystartstep);

	/* the divergence is flagged by icntrl!=0
	   icutb is reset to zero in order to generate
	   regular contact elements etc.. */

	icutb--;
      }
    }
      
    if(icfd==2){
      compfluidfem(&cof,&nkf,&ipkonf,&konf,&lakonf,nef,&sideface,
		   ifreestream,&nfreestream,isolidsurf,neighsolidsurf,
		   &nsolidsurf,iponoelf,inoelf,nshcon,shcon,nrhcon,rhcon,
		   &voldf,ntmat_,nodebounf,ndirbounf,&nbounf,&ipompcf,
		   &nodempcf,&nmpcf,&ikmpcf,&ilmpcf,ithermal,ikbounf,ilbounf,
		   &iturbulent,isolver,iexpl,ttime,&time,&dtime,nodeforc,
		   ndirforc,xforc,nforc,nelemloadf,sideloadf,
		   xloadf,&nloadf,xbody,ipobodyf,nbody,&ielmatf,matname,mi,
		   ncmat_,physcon,istep,&iinc,ibody,xloadold,xbounf,&coefmpcf,
		   nmethod,xforcold,xforcact,iamforc,iamloadf,xbodyold,xbodyact,
		   t1old,t1,t1act,iamt1,amta,namta,nam,ampli,xbounold,xbounact,
		   iambounf,itg,&ntg,amname,t0,&nelemface,&nface,cocon,ncocon,
		   xloadact,tper,jmax,jout,set,nset,istartset,iendset,ialset,
		   prset,prlab,nprint,trab,inotr,ntrans,filab,&labmpc,sti,
		   norien,orab,jobnamef,tieset,ntie,mcs,ics,cs,nkon,&mpcfree,
		   &memmpc_,&fmpc,nef,&inomat,qfx,kode,ipface,ielprop,prop,
		   orname,tincf,&ifreesurface,&nkftot,ielorienf,nelold,nkold,
		   nknew,nelnew);

      for(i=0;i<nkftot;i++){
	for(j=0;j<mt;j++){
	  vold[mt*(nkold[i]-1)+j]=voldf[mt*i+j];
	}
      }
    }
      
    if(icascade==2){
      memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
      RENEW(nodempc,ITG,3*memmpcref_);
      isiz=3*memmpcref_;cpyparitg(nodempc,nodempcref,&isiz,&num_cpus);
      RENEW(coefmpc,double,memmpcref_);
      isiz=memmpcref_;cpypardou(coefmpc,coefmpcref,&isiz,&num_cpus);
    }

    /* generating contact elements */
      
    if((ncont!=0)&&(*mortar<=1)&&

       /*       for purely thermal calculations: determine contact integration
		points only at the start of a step */

       ((*ithermal!=2)||(iit==-1))){

      *ne=ne0;*nkon=nkon0;

      /* at start of new increment: 
	 - copy state variables (node-to-face)
	 - determine slave integration points (face-to-face)
	 - interpolate state variables (face-to-face) */

      if(icutb==0){
	if(*mortar==1){

	  if(*nstate_!=0){
	    if(maxprevcontel!=0){
	      if(iit!=-1){
		NNEW(islavsurfold,ITG,2**ifacecount+2);
		NNEW(pslavsurfold,double,3**nintpoint);
		isiz=2**ifacecount+2;
		cpyparitg(islavsurfold,islavsurf,&isiz,&num_cpus);
		isiz=3**nintpoint;
		cpypardou(pslavsurfold,pslavsurf,&isiz,&num_cpus);
	      }
	    }
	  }

	  *nintpoint=0;

	  /* determine the location of the slave integration
	     points */

	  precontact(&ncont,ntie,tieset,nset,set,istartset,
                     iendset,ialset,itietri,lakon,ipkon,kon,koncont,ne,
                     cg,straight,co,vold,istep,&iinc,&iit,itiefac,
                     islavsurf,islavnode,imastnode,nslavnode,nmastnode,
                     imastop,mi,ipe,ime,tietol,
		     nintpoint,&pslavsurf,xmastnor,cs,mcs,ics,clearini,
                     nslavs);
		  
	  /* changing the dimension of element-related fields */
		  
	  RENEW(kon,ITG,*nkon+22**nintpoint);
	  RENEW(springarea,double,2**nintpoint);
	  RENEW(pmastsurf,double,6**nintpoint);
		  
	  if(*nener==1){
	    RENEW(ener,double,mi[0]*(*ne+*nintpoint)*2);

	    /* setting the entries for the contact energy to zero */

	    DOUMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nintpoint),0.);

	  }
	  RENEW(ipkon,ITG,*ne+*nintpoint);
	  RENEW(lakon,char,8*(*ne+*nintpoint));
		  
	  if(*norien>0){
	    RENEW(ielorien,ITG,mi[2]*(*ne+*nintpoint));
	    ITGMEMSET(ielorien,mi[2]**ne,mi[2]*(*ne+*nintpoint),0);
	  }
	  RENEW(ielmat,ITG,mi[2]*(*ne+*nintpoint));
	  isiz=mi[2]**nintpoint;
	  ITGMEMSET(ielmat,mi[2]**ne,mi[2]*(*ne+*nintpoint),1);

	  /* interpolating the state variables */

	  if(*nstate_!=0){
	    if(maxprevcontel!=0){
	      RENEW(xstateini,double,
		    *nstate_*mi[0]*(ne0+maxprevcontel));
	      isiz=*nstate_*mi[0]*maxprevcontel;
	      cpypardou(&xstateini[*nstate_*mi[0]*ne0],
			&xstate[*nstate_*mi[0]*ne0],&isiz,&num_cpus);
	    }
		      
	    RENEW(xstate,double,*nstate_*mi[0]*(ne0+*nintpoint));
	    isiz=*nstate_*mi[0]**nintpoint;
	    DOUMEMSET(xstate,*nstate_*mi[0]*ne0,
		      *nstate_*mi[0]*(ne0+*nintpoint),0.);
		      
	    if((*nintpoint>0)&&(maxprevcontel>0)){
			  
	      /* interpolation of xstate */
			  
	      interpolatestatemain(ne,ipkon,kon,lakon,
				   &ne0,mi,xstate,pslavsurf,nstate_,
				   xstateini,islavsurf,islavsurfold,
				   pslavsurfold,tieset,ntie,itiefac);
			  
	    }

	    if(maxprevcontel!=0){
	      SFREE(islavsurfold);SFREE(pslavsurfold);
	    }

	    maxprevcontel=*nintpoint;

	    RENEW(xstateini,double,*nstate_*mi[0]*(ne0+*nintpoint));
	    isiz=*nstate_*mi[0]*(ne0+*nintpoint);
	    cpypardou(xstateini,xstate,&isiz,&num_cpus);
	  }

	}

	/* set the contact spring energy to zero at the start of
           an increment. The friction energy is summed in energy[3]
           based on energyini[3] */
	
	if(*nener==1){
	  if((*mortar!=1)||(ncont==0)){
	    DOUMEMSET(enerini,2*mi[0]**ne,2*mi[0]*(*ne+*nslavs),0.);
	  }else{
	    RENEW(enerini,double,2*mi[0]*(*ne+*nintpoint));
	    DOUMEMSET(enerini,2*mi[0]**ne,2*mi[0]*(*ne+*nintpoint),0.);
	  }
	}
	
      }

      /* massless contact: calculate matrix Wb */
      
      if(*mortar==-1){

	if((masslesslinear==0)||(iinc==1)){
	  nzsw=5*9**nslavs;
	  
	  /* 5 = 1 slave + maximal 4 master, so 5 terms in the equation times
	     3 dofs = 15 terms; for each slave node 3 dofs, so 3 equations */ 
	  
	  NNEW(auw,double,nzsw);
	  NNEW(jqw,ITG,3**nslavs+1);
	  NNEW(iroww,ITG,nzsw);
	}
	if((masslesslinear>0)&&(iinc==1)){
	  NNEW(fullgmatrix,double,9**nslavs**nslavs);
	  NNEW(fullr,double,3**nslavs);
	}
      }

      if((*mortar!=-1)||(masslesslinear==0)||(iinc==1)){
	contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
		co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
		&ne0,nmethod,
		iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,
		islavsurf,
		itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
		imastnode,nmastnode,xmastnor,filab,mcs,ics,&nasym,
		xnoels,mortar,pslavsurf,pmastsurf,clearini,&theta,
		xstateini,xstate,nstate_,&icutb,&ialeatoric,jobnamef,
		&alea,auw,jqw,iroww,&nzsw);
      }
   
      /* check whether, for a dynamic calculation, contact damping 
	 is involved */

      if(*nmethod==4){
	if(*iexpl<=1){
	  if(idampingwithoutcontact==0){
	    for(i=0;i<*ne;i++){
	      if(ipkon[i]<0) continue;
	      if(*ncmat_>=5){
		if(strcmp1(&lakon[i*8],"ES")==0){
		  if(strcmp1(&lakon[i*8+6],"C")==0){
		    imat=ielmat[i*mi[2]];
		    if(elcon[(*ncmat_+1)**ntmat_*(imat-1)+4]>0.){
		      idamping=1;break;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
	  
      if(*iexpl<=1) printf(" Number of contact spring elements=%"
			   ITGFORMAT "\n\n",*ne-ne0);
            
      /* carlo start */

      /* dynamic time step estimation for explicit dynamics under penalty 
	 contact(CMT) start */

      if((*iexpl>1)&&(*mortar!=-1)){
	      
	if((*ne-ne0)<ncontacts){

	  /* number of contact elements has decreased */
	  
	  ncontacts=*ne-ne0;
	  inccontact=0;
	}  
	else if((*ne-ne0)>ncontacts)  {

	  /* number of contact elements has increased */
	  
	  RENEW(smscale,double,*ne);
		  
	  FORTRAN(calcstabletimeinccont,(ne,lakon,kon,ipkon,mi,ielmat,elcon,
					 mortar,adb,alpha,nactdof,springarea,
					 &ne0,ntmat_,ncmat_,&dtcont,smscale,
					 &dtset,&mscalmethod));
	  if(dtcont<dtvol){
	    dtmin=dtcont;
	  }else{
	    dtmin=dtvol;
	  }
	  
	  if(dtmin>(*tmax*(*tper))){
	    dtime=*tmax*(*tper);}
	  else if(dtmin<dtset){
	    dtime=dtset;}
	  else {
	    dtime=dtmin;
	  }
	  
	  dtheta=(dtime)/(*tper);
	  reltime=theta+dtheta;
	  time=reltime**tper;
	  dthetaref=dtheta;
	  printf(" SELECTED time increment (based on contact):%e\n\n",dtime);

	  ncontacts=*ne-ne0; 
	  inccontact=0;
	}else if((inccontact==500)&&(ncontacts==0)){

          /* no contact elements (either no contact or massless contact) */

	  if(dtvol>(*tmax*(*tper))){
	    dtime=*tmax*(*tper);
	  }else if(dtvol<dtset){
	    dtime=dtset;
	  }else{
	    dtime=dtvol;
	  }
	  dtheta=(dtime)/(*tper);
	  reltime=theta+dtheta;
	  time=reltime**tper;
	  dthetaref=dtheta;
	  printf(" SELECTED time increment (based on contact):%e\n\n",*tinc);

	  dtcont=1.e30;
	}
	inccontact++;
      }
	  
      /* CMT end */

    }
      
    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's) */
      
    FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		       nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
		       &maxlenmpc,ikmpc,ilmpc,&icascade,
		       kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
		       &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));
      
    if(icascade==2){
      isiz=3*memmpc_;cpyparitg(nodempcref,nodempc,&isiz,&num_cpus);
      isiz=memmpc_;cpypardou(coefmpcref,coefmpc,&isiz,&num_cpus);
    }

    /* recalculating the matrix structure; only needed if:
       1) MPC's are cascaded
       2) contact occurs in an implicit calculation 
          (in penalty explicit no matrices are needed, in massless
           explicit no contact elements are generated) */
      
    if((icascade>0)||((ncont!=0)&&(*iexpl<=1)))
      remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		 &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		 labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		 kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		 neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		 &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		 typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		 nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		 ielmat,matname);

    /* invert nactdof (not for dynamic explicit calculations) */

    if(*iexpl<=1){
      SFREE(nactdofinv);
      NNEW(nactdofinv,ITG,mt**nk);
      MNEW(nodorig,ITG,*nk);
      FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			     ipkon,lakon,kon,ne));
      SFREE(nodorig);
    }
      
    /* check whether the forced displacements changed; if so, and
       if the procedure is static, the first iteration has to be
       purely linear elastic, in order to get an equilibrium
       displacement field; otherwise huge (maybe nonelastic)
       stresses may occur, jeopardizing convergence */
      
    iforbou=0;
      
    /* only for iinc=1 a linearized calculation is performed, since
       for iinc>1 a reasonable displacement field is predicted by using the
       initial velocity field at the end of the last increment */
      
    if((iinc==1)&&(*ithermal<2)&&((*nmethod!=4)||(*mortar==-1))){
      dev=0.;
      for(k=0;k<*nboun;++k){
	err=fabs(xbounact[k]-xbounini[k]);
	if(err>dev){dev=err;}
      }
      if(dev>1.e-5) iforbou=1;
    }
    if((*mortar==-1)&&(iforbou==1)){
      printf(" *ERROR in nonlingeo: nonzero boundary conditions are not allowed\n");
      printf("        in combination with massless contact\n\n");
      FORTRAN(stop,());
    }
      
    /* prediction of the kinematic vectors  */
      
    NNEW(v,double,mt**nk);
    
    /* for massless contact there is no need for prediction,
       since scheme is on velocity level */
      
    /* ---- CCX_PATHFOLLOW: the stock predictor is the wrong predictor ---
       prediction() extrapolates v = vold + dtime*veold, i.e. it repeats
       the previous increment's displacement change.  That is right when
       the step time IS the load parameter.  It is wrong here for two
       independent reasons, both measured on the target at increment 60:

       1. It moves the CONTROL COORDINATE before the corrector runs.  The
          extrapolation advanced phi to 1.65e-04 against a target of
          1.0e-05, so every corrector opened by having to undo 94% of the
          predictor, and a cutback made that worse rather than better
          because it shrinks dphi faster than it shrinks the predictor.

       2. It destroys the capture of f_hat.  f_hat is obtained as
          b/dlamjump at the first iteration, which is -dR/dlambda only if
          the ONLY thing that moved since the committed state is the
          prescribed pattern.  With an extrapolated u, b also contains
          K*du_pred - which is designed to cancel the load term - so the
          quotient is a difference of two nearly equal quantities.
          Measured against an honest finite difference of the residual at
          fixed u, the frozen f_hat came out at cos=+0.983 but
          |f_hat|/|q| = 1.43.  A scale error s there leaves the
          DISPLACEMENT correction right and the LOAD FACTOR step a factor
          1/s short, so the prescribed dofs end up where the free dofs did
          not assume: the constraint row converged to 1e-19 while the
          equilibrium row stalled with lambda marching down 6.7e-05 per
          iteration.  That is exactly the trace that was measured.

       prediction() skips the extrapolation when idiscon != 0 and resets
       the flag itself, so suppressing it is one assignment.  The method
       keeps its OWN predictor: the previous accepted dlambda, applied to
       lambda alone, which is the standard arc-length predictor. */

    if((pf.on==1)&&(pf.engaged==1)) idiscon=1;

    if (*mortar==-1){
      memcpy(&v[0],&vold[0],sizeof(double)*mt **nk);
    }else{
      prediction(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,accold,v,
		 &iinc,&idiscon,vold,nactdof,mi,&num_cpus);
    }
      
    MNEW(fn,double,mt**nk);
    NNEW(stx,double,6*mi[0]**ne);
      
    /* determining the internal forces at the start of the increment
	 
       for a static calculation with increased forced displacements
       the linear strains are calculated corresponding to
	 
       the displacements at the end of the previous increment, extrapolated
       if appropriate (for nondispersive media) +
       the forced displacements at the end of the present increment +
       the temperatures at the end of the present increment (this sum is
       v) -
       the displacements at the end of the previous increment (this is vold)
	 
       these linear strains are converted in stresses by multiplication
       with the tangent element stiffness matrix and converted into nodal
       forces. 
	 
       this boils down to the fact that the effect of forced displacements
       should be handled in a purely linear way at the
       start of a new increment, in order to speed up the convergence and
       (for dissipative media) guarantee smooth loading within the increment.
	 
       for all other cases the nodal force calculation is based on
       the true stresses derived from the appropriate strain tensor taking
       into account the extrapolated displacements at the end of the 
       previous increment + the forced displacements and the temperatures
       at the end of the present increment */
      
    iout=-1;
    if(istrainfree==1) iout=-2;
    iperturb_sav[0]=iperturb[0];
    iperturb_sav[1]=iperturb[1];
      
    /* first iteration in first increment: elastic tangent */
      
    if((*nmethod!=4)&&(iforbou==1)){
	  
      ielas=1;
	  
      iperturb[0]=-1;
      iperturb[1]=0;
	  
      isiz=neq[1];cpypardou(b,f,&isiz,&num_cpus);
      if(ne1d2d==1)NNEW(inum,ITG,*nk);
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t1ini,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd, ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
	      inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavquadel,aut,irowt,jqt,&mortartrafoflag,
	      &intscheme,physcon,dam,damn,iponoel);
      iperturb[0]=0;if(ne1d2d==1)SFREE(inum);
	  
      /* check whether any displacements or temperatures are changed
	 in the new increment */
	  
      for(k=0;k<neq[1];++k){
	f[k]=f[k]+b[k];}
	  
    }
    else{

      if(*mortar!=-1){
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	trial_results(&nlgt);
	if(ne1d2d==1)SFREE(inum);
	  
	isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
	  
	if(*ithermal!=2){
	  isiz=6*mi[0]*ne0;	    
	  cpypardou(sti,stx,&isiz,&num_cpus);
	}
      }  
    }
      
    ielas=0;
    iout=0;

    /* ---- [DAMAGE RELEASE] report ------------------------------------
       f now holds f_int on the NEW topology, at a displacement state that is
       bit-identical to the converged u*: the topology site set idiscon=1
       before looping back, and prediction() with idiscon!=0 copies vold into
       v with no extrapolation (prediction.c:99).  The difference below is
       therefore the released internal force and nothing else.

       Self-checks printed with the numbers, because a probe that is trusted
       without them is worse than none:
         - action=reuse-sparse-graph means nactdof did not change, so
           removed MUST be exactly zero;
         - anom counts DOF active AFTER but not BEFORE, which a deletion
           cannot produce - nonzero means the mapping is wrong;
         - iforbou=1 would mean f had a boundary term added to it
           (f[k]+=b[k] on the other branch above) and the reading is void. */
    if((damage_release_probe)&&(damage_release_armed)&&
       (damage_frel!=NULL)&&(damage_ract!=NULL)){
      ITG ri,rj,rk,rns=0,rds=0,rnr=0,rdr=0,ranom=0;
      double dfv,dfa,smax=0.,sl1=0.,sl2=0.,rmax=0.,rl1=0.,rl2=0.;
      for(ri=0;ri<*nk;ri++){
        for(rj=0;rj<mt;rj++){
          rk=nactdof[mt*ri+rj];
          dfa=(rk>0)?f[rk-1]:0.;
          if(damage_ract[mt*ri+rj]){
            dfv=damage_frel[mt*ri+rj]-dfa;
            if(rk>0){
              sl1+=fabs(dfv); sl2+=dfv*dfv;
              if(fabs(dfv)>smax){smax=fabs(dfv);rns=ri+1;rds=rj;}
            }else{
              rl1+=fabs(dfv); rl2+=dfv*dfv;
              if(fabs(dfv)>rmax){rmax=fabs(dfv);rnr=ri+1;rdr=rj;}
            }
          }else if(rk>0){
            ranom++;
          }
        }
      }
      sl2=sqrt(sl2); rl2=sqrt(rl2);
      printf("[DAMAGE RELEASE] inc=%" ITGFORMAT " pass=%" ITGFORMAT
             " time=%.12e action=%s dt=%.6e%s"
             "   surv:    dF_max=%.6e node=%" ITGFORMAT " dof=%" ITGFORMAT
             " dF_l1=%.6e dF_l2=%.6e%s"
             "   removed: dF_max=%.6e node=%" ITGFORMAT " dof=%" ITGFORMAT
             " dF_l1=%.6e dF_l2=%.6e%s"
             "   qa=%.6e qam=%.6e  dF_surv_max/qam=%.6e  dF_max/qam=%.6e%s"
             "   deleted: terminal=%" ITGFORMAT " deadall+deadsole=%" ITGFORMAT
             " islands=%" ITGFORMAT " cohfacets=%" ITGFORMAT "%s"
             "   selfcheck: anom=%" ITGFORMAT " iforbou=%" ITGFORMAT
             " removed_must_be_zero=%s%s",
             iinc,damage_release_pass,theta**tper,
             damage_release_rebuild?"remastruct":"reuse-sparse-graph",
             damage_release_dt,"\n",
             smax,rns,rds,sl1,sl2,"\n",
             rmax,rnr,rdr,rl1,rl2,"\n",
             damage_release_qa,damage_release_qam,
             (damage_release_qam>0.)?smax/damage_release_qam:-1.,
             (damage_release_qam>0.)?
               ((smax>rmax?smax:rmax)/damage_release_qam):-1.,"\n",
             damage_release_nterm,damage_release_nother,
             damage_release_nisl,damage_release_ncoh,"\n",
             ranom,damage_release_iforbou,
             damage_release_rebuild?"no":"YES","\n");
      fflush(stdout);
      damage_release_armed=0;
    }

      
    SFREE(fn);SFREE(v);
    if((*ithermal!=3)||(ncont==0)||(*mortar!=1)||(*ncmat_<11)) SFREE(stx);
      
    /***************************************************************/
    /* iteration counter and start of the loop over the iterations */
    /***************************************************************/

    if(*mortar>1){	  
      NNEW(bhat,double,neq[1]);
      NNEW(islavactdof,ITG,neq[1]);
    } 
      
    iit=1;

    /* change due to previous checkdivergence routine */

    if(icntrl!=0) icutb++;

    ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
    slow.extended=0;
    slow.camprev1=1.e300;
    slow.camprev2=1.e300;
    slow.maxiters=DAMAGE_SLOW_NEWTON_MAX_ITERS;
    if((ITG)icref>slow.maxiters)
      slow.maxiters=(ITG)icref;
    if(*nmethod!=4)NNEW(resold,double,neq[1]);
    if(uncoupled){
      *ithermal=2;
      NNEW(iruc,ITG,nzs[1]-nzs[0]);
      for(k=0;k<nzs[1]-nzs[0];k++){
	iruc[k]=irow[k+nzs[0]]-neq[0];}
    }

    /* Second half of the rollback probe: PF-START reports |vold-vini| at
       the top of the attempt, this reports it at the last statement before
       the Newton loop.  Anything nonzero here is state that moved between
       the two, which the constraint would otherwise attribute to the
       increment. */

    if((pf.on==1)&&(ccxopt_getenv("CCX_PATHFOLLOW_ACCUMCHECK")!=NULL)){
      double pfd=0.,pft;
      ITG pfi,pfj,pfk;
      for(pfi=0;pfi<*nk;pfi++){
        for(pfj=1;pfj<mt;pfj++){
          pfk=nactdof[mt*pfi+pfj];
          if(pfk>0){
            pft=vold[mt*pfi+pfj]-vini[mt*pfi+pfj];
            pfd+=pft*pft;
          }
        }
      }
      printf("[PF-PRELOOP] inc=%" ITGFORMAT " icutb=%" ITGFORMAT
             " |vold-vini|=%.6e\n",iinc,icutb,sqrt(pfd));
      fflush(stdout);
    }

    while(icntrl==0){

#ifdef COMPANY
      FORTRAN(uiter,(&iit));
#endif	  

      /*  updating the nonlinear mpc's (also affects the boundary
	  conditions through the nonhomogeneous part of the mpc's) */

      if((iit!=1)||((uncoupled)&&(*ithermal==1))){

	printf(" iteration %" ITGFORMAT "\n\n",iit);

	/* restoring the distributed loading before adding the
	   friction heating */
	  
	if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
	  *nload=nloadref;
	  isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
	  if(*nam>0){
	    isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
	  }
	  memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
	}
	  
	FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
			  xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
			  xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,
			  ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
			  xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
			  ndirboun,nodeforc,ndirforc,istep,&iinc,co,vold,itg,
			  &ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
			  ntrans,trab,inotr,veold,integerglob,doubleglob,
			  tieset,istartset,iendset,ialset,ntie,nmpc,ipompc,
			  ikmpc,ilmpc,nodempc,coefmpc,ipobody,iponoeln,inoeln,
			  ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,rhcon,
			  nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));

	for(i=0;i<3;i++){
	  cam[i]=0.;}
	for(i=3;i<5;i++){
	  cam[i]=0.5;}
	if(*ithermal>1){
	  radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,ac,bc,nload,
		      sideload,nelemload,xloadact,lakon,ipiv,ntmat_,vold,shcon,
		      nshcon,ipkon,kon,co,kontri,&ntri,nloadtr,tarea,tenv,
		      physcon,erad,&adview,&auview,nflow,ikboun,xbounact,nboun,
		      ithermal,&iinc,&iit,cs,mcs,inocs,&ntrit,nk,fenv,istep,
		      &dtime,ttime,&time,ilboun,ikforc,ilforc,xforcact,nforc,
		      cam,ielmat,&nteq,prop,ielprop,nactdog,nacteq,nodeboun,
		      ndirboun,network,rhcon,nrhcon,ipobody,ibody,xbodyact,
		      nbody,iviewfile,jobnamef,ctrl,xloadold,&reltime,nmethod,
		      set,mi,istartset,iendset,ialset,nset,ineighe,nmpc,
		      nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,iamload,
		      jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,ncocon,
		      iponoeln,inoeln,nprop,amname,namta,amta,iexpl);
             
	  /* check whether network iterations converged */

	  if(qa[2]>0){
	    checkdivergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod,kode,filab,
			    een,t1act,&time,epn,ielmat,matname,enern,xstaten,
			    nstate_,istep,&iinc,iperturb,ener,mi,output,
			    ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			    ielorien,norien,description,sti,&icutb,&iit,&dtime,
			    qa,vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,
			    &icntrl,&theta,&dtheta,veold,vini,idrct,tper,
			    &istab,tmax,nactdof,b,tmin,ctrl,amta,namta,itpamp,
			    inext,&dthetaref,&itp,&jprint,jout,&uncoupled,t1,
			    &iitterm,nelemload,nload,nodeboun,nboun,itg,
			    ndirboun,&deltmx,&iflagact,set,nset,istartset,
			    iendset,ialset,emn,thicke,jobnamec,mortar,nmat,
			    ielprop,prop,&ialeatoric,&kscale,energy,&allwk,
			    &energyref,&emax,&r_abs,&enetoll,energyini,
			    &allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			    &dampwkini,energystartstep);
	    continue;
	  }
	}

	if(icascade==2){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
	  RENEW(nodempc,ITG,3*memmpcref_);
	  isiz=3*memmpcref_;cpyparitg(nodempc,nodempcref,&isiz,&num_cpus);
	  RENEW(coefmpc,double,memmpcref_);
	  isiz=memmpcref_;cpypardou(coefmpc,coefmpcref,&isiz,&num_cpus);
	}

	if((ncont!=0)&&(*mortar<=1)&&(ismallsliding==0)&&
	   /*           for node-to-face contact: freeze contact elements for
			iterations 8 and higher */
	   ((iit<=8)||(*mortar==1))&&
	   /*           for purely thermal calculations: freeze contact elements
			during complete step */
	   ((*ithermal!=2)||(iit==-1))){

	  neold=*ne;
	  *ne=ne0;*nkon=nkon0;
	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
		  straight,nkon,co,vold,ielmat,cs,elcon,istep,
		  &iinc,&iit,ncmat_,ntmat_,&ne0,
		  nmethod,iperturb,
		  ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
		  itiefac,areaslav,iponoels,inoels,springarea,tietol,
		  &reltime,imastnode,nmastnode,xmastnor,
		  filab,mcs,ics,&nasym,xnoels,mortar,pslavsurf,pmastsurf,
		  clearini,&theta,xstateini,xstate,nstate_,&icutb,
		  &ialeatoric,jobnamef,&alea,auw,jqw,iroww,&nzsw);

	  /* check whether, for a dynamic calculation, contact damping is 
	     involved */
	      
	  if(*nmethod==4){
	    if(*iexpl<=1){
	      if(idampingwithoutcontact==0){
		for(i=0;i<*ne;i++){
		  if(ipkon[i]<0) continue;
		  if(*ncmat_>=5){
		    if(strcmp1(&lakon[i*8],"ES")==0){
		      if(strcmp1(&lakon[i*8+6],"C")==0){
			imat=ielmat[i*mi[2]];
			if(elcon[(*ncmat_+1)**ntmat_*(imat-1)+4]>0.){
			  idamping=1;break;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	      
	  if(*mortar==0){
	    if(*ne!=neold){iflagact=1;}
	  }else if(*mortar==1){
	    if(((*ne-ne0)<(neold-ne0)*(1.-delcon))||
	       ((*ne-ne0)>(neold-ne0)*(1.+delcon))){iflagact=1;}
	  }

	  printf(" Number of contact spring elements=%" ITGFORMAT "\n\n",
		 *ne-ne0);

	}
	  
	if(*ithermal==3){
	  for(k=0;k<*nk;++k){
	    t1act[k]=vold[mt*k];}
	}

	FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			   nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			   &maxlenmpc,ikmpc,ilmpc,&icascade,
			   kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,&iit,
			   &idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));

	if(icascade==2){
	  isiz=3*memmpc_;cpyparitg(nodempcref,nodempc,&isiz,&num_cpus);
	  isiz=memmpc_;cpypardou(coefmpcref,coefmpc,&isiz,&num_cpus);
	}

	/* recalculating the matrix structure */

	/* for face-to-face contact (mortar=1) this is only done if
	   the dependent term in nonlinear MPC's changed */
	
	if((icascade>0)||(ncont!=0)){
	  if((*mortar!=1)||(kchdep==1)){
	    remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		       &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,
		       ilboun,labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		       kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		       neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		       &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		       typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		       nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		       ielmat,matname);
	  }

	  /* invert nactdof */
	      
	  SFREE(nactdofinv);
	  NNEW(nactdofinv,ITG,mt**nk);
	  MNEW(nodorig,ITG,*nk);
	  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
				 ipkon,lakon,kon,ne));
	  SFREE(nodorig);
	      
	  MNEW(v,double,mt**nk);
	  NNEW(stx,double,6*mi[0]**ne);
	  MNEW(fn,double,mt**nk);
      
	  isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
	  iout=-1;
	      
	  if(ne1d2d==1)NNEW(inum,ITG,*nk);
	  trial_results(&nlgt);
	  
	  isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
	      
	  if(*ithermal!=2){
	    isiz=6*mi[0]*ne0;	    
	    cpypardou(sti,stx,&isiz,&num_cpus);
	  }
	      
	  SFREE(v);SFREE(fn);if(ne1d2d==1)SFREE(inum);
	  if((*ithermal!=3)||(ncont==0)||(*mortar!=1)||(*ncmat_<11)) SFREE(stx);
	  iout=0;
	}    
      }
	  
      /* add friction heating  */
      
      if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
	nload_=*nload+2*(*ne-ne0);

	RENEW(nelemload,ITG,2*nload_);
	ITGMEMSET(nelemload,2**nload,2*nload_,0);
	if(*nam>0){
	  RENEW(iamload,ITG,2*nload_);
	  ITGMEMSET(iamload,2**nload,2*nload_,0);
	}
	RENEW(xloadact,double,2*nload_);
	DOUMEMSET(xloadact,2**nload,2*nload_,0.);
	RENEW(sideload,char,20*nload_);
	DMEMSET(sideload,20**nload,20*nload_,'\0');

	MNEW(idefload,ITG,nload_);
	ITGMEMSET(idefload,0,nload_,1);
	FORTRAN(frictionheating,(&ne0,ne,ipkon,lakon,ielmat,mi,elcon,ncmat_,
				 ntmat_,kon,islavsurf,pmastsurf,springarea,co,
				 vold,veold,pslavsurf,xloadact,nload,&nload_,
				 nelemload,iamload,idefload,sideload,stx,nam,
				 &time,ttime,matname,istep,&iinc));
	SFREE(idefload);SFREE(stx);
      }

      /* calculate the stiffness matrix for:
         - implicit calculations
         - linear massless explicit calculations in the first increment
         - nonlinear massless explicit calculations */
      
      if((*iexpl<=1)||((*mortar==-1)&&((masslesslinear==0)||(iinc==1)))){

	/* calculating the local stiffness matrix and external loading */

	NNEW(ad,double,neq[1]);
	NNEW(au,double,nzs[1]);

	if(*nmethod==4){
	  DOUMEMSET(fnext,0,mt**nk,0.);
	}

	/* UNSYM stage 1: route the bulk assembly through the existing
	   asymmetric storage while the constitutive tangent is still the
	   stock symmetric g(D)*C_ep.

	   mafillsmmain writes the symmetric part into the lower triangle
	   exactly as before; mafillsmasmain then mirrors it into the upper
	   half at offset nzs[2] and adds asymmetric element contributions.
	   With no contact elements present it only performs the mirror, so
	   the assembled operator is numerically identical to the symmetric
	   one and PARDISO mtype=1 must reproduce the symmetric answer to
	   the last digit.  That is the acceptance test for this stage: it
	   validates the plumbing before any damage term is added to it.

	   Gates: implicit mechanical static path, no contact (contact owns
	   nasym itself), progressive damage material present. */

	opd.unsym_active=0;
	if((damage_tangent_mode==2)&&(damage_de12_enabled)&&
	   (ncont==0)&&(*iexpl<=1)&&(*nmethod!=4)&&(*nmethod!=5)&&
	   (*ithermal<2)&&(*mortar!=-1)){
	  opd.unsym_active=1;
	  if(nasym==0){
	    nasym=1;
	    printf("[DAMAGE TANGENT UNSYM] bulk assembly switched to the "
		   "asymmetric path (nasym=1); PARDISO symbolic reuse is "
		   "available here (mtype=1, structurally "
		   "symmetric) - set CCX_PARDISO_REUSE_SYMBOLIC=1\n");
	    fflush(stdout);
	  }
	}

	mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
		     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		     nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		     nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		     nmethod,ikmpc,ilmpc,ikboun,ilboun,
		     elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		     ielmat,ielorien,norien,orab,ntmat_,
		     t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		     nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		     xstiff,npmat_,&dtime,matname,mi,
		     ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		     physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		     &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		     xstateini,xstate,thicke,integerglob,doubleglob,
		     tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		     pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		     iponoeln,inoeln,network,ntrans,inotr,trab,smscale,
		     &mscalmethod,set,nset,islavquadel,aut,irowt,jqt,
		     &mortartrafoflag);
	//		     &nslavquadel);

	if(nasym==1){
	  RENEW(au,double,2*nzs[1]);
	  if(*nmethod==4){
	    RENEW(aub,double,2*nzs[1]);}
	  symmetryflag=2;
	  inputformat=1;

	  mafillsmasmain(co,nk,kon,ipkon,lakon,ne,nodeboun,
			 ndirboun,xbounact,nboun,
			 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
			 nforc,nelemload,sideload,xloadact,nload,xbodyact,
			 ipobody,
			 nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
			 nmethod,ikmpc,ilmpc,ikboun,ilboun,
			 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
			 ielmat,ielorien,norien,orab,ntmat_,
			 t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
			 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			 xstiff,npmat_,&dtime,matname,mi,
			 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
			 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,
			 &iinc,
			 &coriolis,ibody,xloadold,&reltime,veold,springarea,
			 nstate_,
			 xstateini,xstate,thicke,
			 integerglob,doubleglob,tieset,istartset,iendset,
			 ialset,ntie,&nasym,pslavsurf,pmastsurf,mortar,clearini,
			 ielprop,prop,&ne0,&kscale,iponoeln,inoeln,network,set,
			 nset);

	  /* UNSYM stage 2: the rank-1 damage correction

	         -sigma_eff (x) dD/d(eps)

	     is added on top of the symmetric operator that mafillsmasmain
	     has just mirrored into both halves of au.  Only integration
	     points that are actively softening carry a nonzero damjac, so
	     the extra work scales with the size of the process zone, not
	     with the model. */

	  if(opd.unsym_active==1){
	    opd.unsym_elems=0;
	    FORTRAN(mafilldamas,(co,kon,ipkon,lakon,&ne0,nactdof,jq,irow,
				 neq,nzs,au,ad,vold,mi,damage_damjac,nmpc,
				 &opd.unsym_elems,dam,damdamageini,
				 &opd.unsym_skip,&opd.unsym_adv,
				 &opd.unsym_hole,&opd.unsym_floor,
				 &opd.unsym_live,&opd.unsym_degen,
				 damage_damcat));
	    /* J-10: what the operator looks like on THIS iteration, not as a
	       running maximum.  A count that moves between iterations of a
	       SAME-LOAD solve means the operator is changing shape mid-Newton,
	       which no step size and no tolerance can repair. */
	    if((opd.unsym_census==1)&&(idamagereeq==1)){
	      printf("[DAMAGE TANGENT CENSUS] inc=%" ITGFORMAT " iit=%" ITGFORMAT
		     " assembled=%" ITGFORMAT " adv=%" ITGFORMAT
		     " gap=%" ITGFORMAT " hi=%" ITGFORMAT
		     " skip=%" ITGFORMAT " live=%" ITGFORMAT
		     " degen=%" ITGFORMAT "\n",
		     iinc,iit,opd.unsym_elems,opd.unsym_adv,
		     opd.unsym_hole,opd.unsym_floor,opd.unsym_skip,
		     opd.unsym_live,opd.unsym_degen);
	      fflush(stdout);
	    }
	    if(opd.unsym_adv>opd.unsym_advrep){
	      opd.unsym_advrep=opd.unsym_adv;
	      printf("[DAMAGE TANGENT HOLE] %" ITGFORMAT " element(s) with "
		     "ADVANCING damage carry no rank-1 term (new maximum); "
		     "%" ITGFORMAT " past initiation without one, "
		     "%" ITGFORMAT " assembled\n",
		     opd.unsym_adv,opd.unsym_skip,opd.unsym_elems);
	      fflush(stdout);
	    }

	    /* The number above conflates three populations and only one is a
	       defect.  dambase is the damage at the START OF THE INCREMENT, so
	       "advancing" there also covers every point that advanced earlier
	       in the increment and unloads elastically now - zero is correct
	       for those, and for points already on the residual-stiffness
	       floor, where dg/dD=0.  What is NOT correct is the band
	       D in [0.999, 1-gmin): resultsmech.f refuses the consistent
	       tangent at D>=0.999 while the floor only engages at 1-D<gmin,
	       so those elements carry g*C_ep with g in (1e-4,1e-3] - which
	       mpfd measures at rho=3443 with the WRONG SIGN on the dominant
	       entry already at D=0.994.  Report that band on its own so the
	       hypothesis can be falsified rather than argued. */
	    if(opd.unsym_hole>opd.unsym_holerep){
	      opd.unsym_holerep=opd.unsym_hole;
	      printf("[DAMAGE TANGENT GAP] inc=%" ITGFORMAT " %" ITGFORMAT
	    	 " element(s) in D=[0.999,1-gmin) carry g*C_ep with no rank-1 "
	    	 "term (new maximum); %" ITGFORMAT " on the residual-stiffness "
	    	 "floor where zero is correct; %" ITGFORMAT " advancing "
	    	 "without a term in total\n",
	    	 iinc,opd.unsym_hole,opd.unsym_floor,opd.unsym_adv);
	      fflush(stdout);
	    }
	    if(opd.unsym_elems>opd.unsym_report){
	      opd.unsym_report=opd.unsym_elems;
	      printf("[DAMAGE TANGENT UNSYM] rank-1 correction assembled "
		     "for %" ITGFORMAT " softening element(s) (new maximum)\n",
		     opd.unsym_elems);
	      fflush(stdout);
	    }
	  }
	}

	/* Per-node stiffness taken from the assembled diagonal.

	   This is the measure the topological rule in damdangle.f cannot
	   reach.  Node 501 keeps 3 of its original 8 elements and is only
	   the 63rd softest of 2387 nodes by mesh geometry - 62 softer nodes
	   cause no trouble at all - yet it is the one the solver throws
	   0.69 of a specimen length.  What separates it is the damage state
	   of the elements it still has, and ad[] carries g(D), so it reports
	   the stiffness the node actually has rather than the one its
	   connectivity suggests.

	   The reference is the node's OWN diagonal while the model was
	   intact, so the ratio is dimensionless and needs no median over the
	   mesh.  A constrained direction is excluded: a node held by a
	   boundary condition is never stranded. */

	/* The solver frees ad and au right after the factorisation, so a
	   probe that runs after results() would read an empty matrix -
	   which is exactly what the elastic control caught, reporting a
	   relative error of 1.0 because every coefficient read back was
	   zero.  Keep a copy while the operator still exists. */

	if((opd.fd_inc>0)&&(iinc>=opd.fd_inc)&&(iit>=opd.fd_it)&&
	   ((opd.fd_step<=0)||(*istep==opd.fd_step))){
	  if(opd.fd_ad==NULL){
	    NNEW(opd.fd_ad,double,neq[1]);
	    NNEW(opd.fd_au,double,(nasym+1)*nzs[1]);
	  }
	  memcpy(opd.fd_ad,ad,sizeof(double)*neq[1]);
	  memcpy(opd.fd_au,au,sizeof(double)*(nasym+1)*nzs[1]);
	}

	/* Addressed dump of one node.

	   Global counters have twice disagreed with each other here, so this
	   asks the question directly for a named node instead: every element
	   that touches it, its type, whether it is still assembled, and for a
	   cohesive element the damage at each integration point.  The
	   assembled diagonal is printed alongside, so the claim "this node has
	   lost its support" is backed by the operator rather than inferred
	   from connectivity. */

	if((prb.dump_node>0)&&(prb.dump_node<=*nk)&&
	   (iinc>=prb.dump_inc)){
	  prb.dump_n=prb.dump_node-1;
	  printf("[NODE DUMP] inc=%" ITGFORMAT " node %" ITGFORMAT
		 " at (%.4f, %.4f, %.4f)\n",
		 iinc,prb.dump_node,co[3*prb.dump_n],
		 co[3*prb.dump_n+1],co[3*prb.dump_n+2]);
	  for(idir=1;idir<=3;idir++){
	    k=nactdof[mt*prb.dump_n+idir];
	    if(k<=0){
	      printf("[NODE DUMP]   dof %" ITGFORMAT ": constrained\n",idir);
	    }else{
	      printf("[NODE DUMP]   dof %" ITGFORMAT ": ad=%.6e",idir,ad[k-1]);
	      if((damage_addiag0!=NULL)&&(damage_addiag0[prb.dump_n]>0.))
		printf("   intact reference %.6e   ratio %.4e",
		       damage_addiag0[prb.dump_n],
		       ad[k-1]/damage_addiag0[prb.dump_n]);
	      printf("\n");
	    }
	  }
	  prb.dump_nb=0;prb.dump_nu=0;
	  for(i=0;i<*ne;i++){
	    if(ipkon[i]==-1) continue;
	    prb.dump_alive=(ipkon[i]>=0)?1:0;
	    prb.dump_idx=(ipkon[i]>=0)?ipkon[i]:(-ipkon[i]-2);
	    if(prb.dump_idx<0) continue;
	    prb.dump_np=0;
	    if(lakon[8*i]=='C'){
	      prb.dump_np=(lakon[8*i+3]=='4')?4:0;
	    }else if(lakon[8*i]=='U'){
	      prb.dump_np=(ITG)((unsigned char)lakon[8*i+7]);
	      if((prb.dump_np<1)||(prb.dump_np>20)) prb.dump_np=0;
	    }
	    if(prb.dump_np<=0) continue;
	    prb.dump_hit=0;
	    for(j=0;j<prb.dump_np;j++){
	      if(kon[prb.dump_idx+j]-1==prb.dump_n) prb.dump_hit=1;
	    }
	    if(prb.dump_hit==0) continue;
	    if(lakon[8*i]=='C'){
	      if(prb.dump_alive) prb.dump_nb++;
	      printf("[NODE DUMP]   elem %-7" ITGFORMAT " C3D4  %-8s dam=%.6f\n",
		     i+1,prb.dump_alive?"alive":"DELETED",
		     (i<ne0)?dam[mi[0]*i]:-1.);
	    }else{
	      if(prb.dump_alive) prb.dump_nu++;
	      printf("[NODE DUMP]   elem %-7" ITGFORMAT " UC6   %-8s",
		     i+1,prb.dump_alive?"alive":"DELETED");
	      for(j=0;j<3;j++){
		prb.dump_dv=xstate[*nstate_*(mi[0]*i+j)+1];
		printf("  ip%" ITGFORMAT ": dvisc=%.6f g=%.6e",
		       j+1,prb.dump_dv,1.-prb.dump_dv);
	      }
	      printf("\n");
	    }
	  }
	  printf("[NODE DUMP]   live bulk=%" ITGFORMAT
		 "   live cohesive=%" ITGFORMAT "\n",
		 prb.dump_nb,prb.dump_nu);
	  /* The node's own 3x3 block of the assembled operator.
	
	     damage_addiag is the smallest AXIS-ALIGNED diagonal, which is only a
	     proxy for how compliant the node is: a node whose soft direction is
	     skew to x, y and z can carry three healthy diagonals while this block
	     has a small eigenvalue.  Three cohesive facets meeting at a node
	     define a plane whose normal is generally not an axis, so that is
	     exactly where the proxy should be weakest (E-79).  Printing the block
	     lets the eigenvalue be compared with the diagonal offline instead of
	     assumed.
	
	     Storage, from pardiso.c inputformat==3: the off-diagonal terms are a
	     full CSC of BOTH triangles, column by column, jq giving the 1-based
	     start of each column and irow the row.  A symmetric assembly stores
	     only the lower triangle, so it is mirrored. */
	  if((au!=NULL)&&(jq!=NULL)&&(irow!=NULL)){
	    ITG dmp_i,dmp_j,dmp_k,dmp_dof[3];
	    double dmp_blk[9];
	    for(dmp_i=0;dmp_i<9;dmp_i++) dmp_blk[dmp_i]=0.;
	    for(dmp_i=0;dmp_i<3;dmp_i++)
	      dmp_dof[dmp_i]=nactdof[mt*prb.dump_n+dmp_i+1];
	    for(dmp_j=0;dmp_j<3;dmp_j++){
	      if(dmp_dof[dmp_j]<=0) continue;
	      dmp_blk[dmp_j*3+dmp_j]=ad[dmp_dof[dmp_j]-1];
	      for(dmp_k=jq[dmp_dof[dmp_j]-1];dmp_k<=jq[dmp_dof[dmp_j]]-1;dmp_k++){
	        for(dmp_i=0;dmp_i<3;dmp_i++){
	          if(dmp_dof[dmp_i]<=0) continue;
	          if(irow[dmp_k-1]==dmp_dof[dmp_i]){
	            dmp_blk[dmp_i*3+dmp_j]=au[dmp_k-1];
	            /* au holds (nasym+1)*nzs[1] entries: the lower triangle in
	               [0,nzs) and, when nasym==1, its TRANSPOSE in [nzs,2nzs) -
	               the same layout the damping assembly uses a few hundred
	               lines below.  Reading only the first half produced a block
	               with a zero upper triangle, which is how this was caught. */
	            dmp_blk[dmp_j*3+dmp_i]=nasym?au[nzs[1]+dmp_k-1]:au[dmp_k-1];
	          }
	        }
	      }
	    }
	    printf("[NODE DUMP]   block3x3 nasym=%" ITGFORMAT,nasym);
	    for(dmp_i=0;dmp_i<9;dmp_i++) printf(" %.6e",dmp_blk[dmp_i]);
	    printf("\n");
	  }
	  fflush(stdout);
	}

	if(damage_stiff_probe>0){
	  /* [DAMSTATE] One owner for the judgement.  This block used to compute
	     the per-node diagonal, its intact reference and the AUTOSPC mask
	     inline; damstate.c now owns all three and the globals below are
	     views onto it, so every consumer - the displacement norm in
	     resultsini.c, the force norm here, the termination connectivity -
	     reads ONE decision instead of three separate ones.  Behaviour is
	     unchanged: damstate_update reproduces this loop exactly, including
	     that the minimum over the three dofs is a real minimum when a
	     negative diagonal is present, which its self test pins. */
	  if(damage_dstate.nk==0)
	    damstate_init(&damage_dstate,*nk,damage_spc_g,damage_spc_neg);
	  damstate_update(&damage_dstate,ad,nactdof,mt);
	  damage_addiag=damage_dstate.diag;
	  damage_addiag0=damage_dstate.diag0;
	  damage_addok=damage_dstate.ok;
	  if(damage_spc_g>0.){
	    damage_spc_mask=damage_dstate.dead;
	    damage_spc_nk=*nk;
	    damage_spc_count=damage_dstate.ndead;
	  }
	}

	iperturb[0]=iperturb_sav[0];
	iperturb[1]=iperturb_sav[1];

      }else{

	/* calculating the external loading 

	   This is only done once per increment. In reality, the
           external loading is a function of vold (specifically,
           the body forces and surface loading). This effect is
           neglected, since the increment size in dynamic explicit
           calculations is usually small */

	if((*mortar==-1)&&(masslesslinear==1)&&(iinc==2)){

	  /* check whether the distributed loading changes in this step 
	   (only for linear massless explicit dynamic calculations) */
	  
	  FORTRAN(checktempload,(iamload,nload,sideload,ibody,nbody,
				 &masslesslinear,&nloadrhs,&nbodyrhs,nam));

	  /* if no change: calculate the external force vector due to this
             loading only once at the start of the step */

	  if(masslesslinear==2){
	    NNEW(fextload,double,neq[1]);
	    nforcrhs=0;
	    rhsmain(co,nk,kon,ipkon,lakon,ne,
		    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		    &nforcrhs,nelemload,sideload,xloadact,nload,xbodyact,
		    ipobody,nbody,cgr,fextload,nactdof,&neq[1],
		    nmethod,ikmpc,ilmpc,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		    ielmat,ielorien,norien,orab,ntmat_,
		    t0,t1act,ithermal,iprestr,vold,iperturb,
		    iexpl,plicon,nplicon,plkcon,nplkcon,
		    npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		    xbodyold,&reltime,veold,matname,mi,ikactmech,
		    &nactmech,ielprop,prop,sti,xstateini,xstate,nstate_,
		    ntrans,inotr,trab,fnext);
	  }
	}

	/* call to rhsmain in every increment: if the distributed
           loading does not change only point forces are taken into account */
	
	rhsmain(co,nk,kon,ipkon,lakon,ne,
	  	ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	  	nforc,nelemload,sideload,xloadact,&nloadrhs,xbodyact,ipobody,
	  	&nbodyrhs,cgr,fext,nactdof,&neq[1],
	  	nmethod,ikmpc,ilmpc,
	  	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	  	ielmat,ielorien,norien,orab,ntmat_,
	  	t0,t1act,ithermal,iprestr,vold,iperturb,
	  	iexpl,plicon,nplicon,plkcon,nplkcon,
	  	npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	  	xbodyold,&reltime,veold,matname,mi,ikactmech,
	  	&nactmech,ielprop,prop,sti,xstateini,xstate,nstate_,
	        ntrans,inotr,trab,fnext);
	//   for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}

	/* adding fextload due to distributed loading */

	if(masslesslinear==2){
	  for(i=0;i<neq[1];i++){fext[i]+=fextload[i];}
	}

      }
      

      /* calculating the damping matrix for implicit dynamic
         calculations */

      if((idamping==1)&&(*iexpl<=1)){

	/* Rayleigh damping */

	MNEW(adc,double,neq[1]);DOUMEMSET(adc,neq[0],neq[1],0.);
	for(k=0;k<neq[0];k++){
	  adc[k]=alpham*adb[k]+betam*ad[k];}
	if(nasym==0){
	  MNEW(auc,double,nzs[1]);DOUMEMSET(auc,nzs[0],nzs[1],0.);
	  for(k=0;k<nzs[0];k++){
	    auc[k]=alpham*aub[k]+betam*au[k];}
	}else{
	  NNEW(auc,double,2*nzs[1]);DOUMEMSET(auc,2*nzs[0],2*nzs[1],0.);
	  for(k=0;k<2*nzs[0];k++){
	    auc[k]=alpham*aub[k]+betam*au[k];}
	}
 
	/* dashpots and contact damping */

	FORTRAN(mafilldm,(co,nk,kon,ipkon,lakon,ne,nodeboun,
			  ndirboun,xbounact,nboun,
			  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,
			  xforcact,
			  nforc,nelemload,sideload,xloadact,nload,xbodyact,
			  ipobody,nbody,cgr,
			  adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
			  ikmpc,ilmpc,ikboun,ilboun,
			  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
			  ielorien,norien,orab,ntmat_,
			  t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
			  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			  xstiff,npmat_,&dtime,matname,mi,ncmat_,
			  ttime,&time,istep,&iinc,ibody,clearini,mortar,
			  springarea,
			  pslavsurf,pmastsurf,&reltime,&nasym));
      }

      /* calculating the residual (RHS of equation system) */

      if(*mortar!=-1){
	trial_reduce(&nlgt,b);
      }else{
	NNEW(volddof,double,neq[0]);
	if(ncont!=0){NNEW(qb,double,neqtot);}
        massless(kslav,lslav,ktot,ltot,au,ad,auc,adc,jq,irow,neq,nzs,auw,jqw,
		 iroww,&nzsw,islavnode,nslavnode,nslavs,imastnode,nmastnode,
		 ntie,nactdof,mi,vold,volddof,veold,nk,fext,isolver,
		 &masslesslinear,co,springarea,&neqtot,qb,b,&dtime,aloc,fric,
		 iexpl,nener,ener,ne,&jqbi,&aubi,&irowbi,&jqib,&auib,&irowib,
		 &iclean,&iinc,fullgmatrix,fullr,alglob,&num_cpus,&ncont);
        if(masslesslinear==0){SFREE(ad);SFREE(au);} 
      }
      
      /*    for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
	    for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
	    for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	    for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
	    for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}*/

      /* mortar contact */

      if(*mortar>1){
	
	/* trafo u -> util */
	
	premortar(nzs,&nzsc2,&auc2,&adc2,
		  &irowc2,&icolc2,&jqc2,&aubd,&irowbd,&jqbd,&aubdtil,
		  &irowbdtil,&jqbdtil,&aubdtil2,&irowbdtil2,&jqbdtil2,
		  &audd,&irowdd,&jqdd,&auddtil,&irowddtil,&jqddtil,
		  &auddtil2,&irowddtil2,&jqddtil2,&auddinv,&irowddinv,
		  &jqddinv,&jqtemp,&irowtemp,&icoltemp,nzstemp,&iit,
		  icol,irow,jq,ikboun,ilboun,ikmpc,ilmpc,
		  imastnode,nmastnode,co,nk,kon,ipkon,lakon,ne,stn,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
		  iprestr,filab,eme,emn,een,iperturb,nactdof,&iout,qa,
		  vold,b,nodeboun,ndirboun,xbounact,xboun,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,neq,veold,accold,
		  &dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,stiini,vini,ener,enern,emeini,xstaten,
		  eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  nelemload,nload,istep,&iinc,springarea,&reltime,&ne0,
		  xforc,nforc,thicke,shcon,nshcon,sideload,xload,xloadold,
		  &icfd,inomat,islavquadel,islavsurf,iponoels,inoels,
		  mortar,nslavnode,
		  islavnode,nslavs,ntie,aut,irowt,jqt,autinv,
		  irowtinv,jqtinv,tieset,
		  itiefac,&rhsi,au,ad,&f_cm,&f_cs,t1act,cam,&bet,&gam,epn,
		  xloadact,nodeforc,ndirforc,xforcact,xbodyact,ipobody,
		  nbody,cgr,nzl,sti,iexpl,mass,&buckling,&stiffness,
		  &intscheme,physcon,&coriolis,ibody,integerglob,
		  doubleglob,&nasym,&alpham,&betam,pslavsurf,
		  pmastsurf,clearini,ielprop,prop,islavact,cdn,&memmpc_,
		  &idamping,&iforbou,iperturb_sav,
		  itietri,cg,straight,koncont,energyini,energy,&kscale,
		  iponoeln,inoeln,nener,orname,network,typeboun,&num_cpus,
		  t0g,t1g,smscale,&mscalmethod,&nslavquadel,iponoel);
	
	/* calculating coupling matrices and embedding weak 
	   contact conditions */ 
      
	contactmortar(&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,
		      itietri,lakon,ipkon,kon,koncont,ne,cg,straight,co,vold,
		      ielmat,elcon,istep,&iinc,&iit,ncmat_,ntmat_,&ne0,vini,
		      nmethod,neq,nzs,nactdof,itiefac,islavsurf,islavnode,
		      imastnode,nslavnode,nmastnode,ad,&au,b,&irow,icol,jq,
		      imastop,iponoels,inoels,&nzsc2,&auc2,adc2,&irowc2,jqc2,
		      islavact,gap,slavnor,slavtan,bhat,&irowbd,jqbd,&aubd,
		      &irowbdtil,jqbdtil,&aubdtil,&irowbdtil2,jqbdtil2,
		      &aubdtil2,&irowdd,jqdd,&audd,&irowddtil,jqddtil,&auddtil,
		      &irowddtil2,jqddtil2,&auddtil2,&irowddinv,jqddinv,
		      &auddinv,irowt,jqt,aut,irowtinv,jqtinv,
		      autinv,mi,ipe,ime,tietol,cstress,cstressini,
		      bp,nk,nboun,ndirboun,nodeboun,xbounact,nmpc,
		      ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		      nslavspc,islavspc,nslavmpc,islavmpc,
		      nmastmpc,imastmpc,
		      pslavdual,islavactdof,islavtie,
		      plicon,nplicon,npmat_,nelcon,&dtime,islavnodeinv,&Bd,
		      &irowb,jqb,&Bdhelp,&irowbhelp,jqbhelp,&Dd,&irowd,jqd,
		      &Ddtil,&irowdtil,jqdtil,&Bdtil,&irowbtil,jqbtil,
		      &bet,cfsini,
		      &reltime,ithermal,plkcon,nplkcon);
	  
	nzs[0]=nzs[1];
	nzs[2]=nzs[1];
	symmetryflag=2;
	inputformat=3; 
      }

      /* storing the residuum in resold (for line search) */


      if((((*mortar==1)&&(iit!=1)&&(*ne-ne0>0)&&(*nmethod!=4))||
          ((rsc.linesearch_mode==1)&&(damage_de12_enabled)&&
           (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
           (*ithermal<2)&&(*idrct==0)))){
	isiz=neq[1];cpypardou(resold,b,&isiz,&num_cpus);
      }
	  
      newstep=0;
      
      if(*nmethod==0){
	  
	/* error occurred in mafill: storing the geometry in frd format */
	  
	*nmethod=0;
	++*kode;
	NNEW(inum,ITG,*nk);ITGMEMSET(inum,0,*nk,1);
	if(strcmp1(&filab[1044],"ZZS")==0){
	  NNEW(neigh,ITG,40**ne);
	  MNEW(ipneigh,ITG,*nk);
	}
	  
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti,damn,&errn);

	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);} 
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
	SFREE(inum);
	if(nmethodold==0){FORTRAN(stopwithout201,());}else{FORTRAN(stop,());}
	  
      }
      
      /* implicit step (static or dynamic) */
      
      if(*iexpl<=1){
	if((*nmethod==4)&&(*mortar<2)){
	      
	  /* mechanical part */
	      
	  if(*ithermal!=2){
	    scal1=bet*dtime*dtime*(1.+alpha[0]);
	    for(k=0;k<neq[0];++k){
	      ad[k]=adb[k]+scal1*ad[k];
	    }
	    for(k=0;k<nzs[0];++k){
	      au[k]=aub[k]+scal1*au[k];
	    }
		  
	    /* upper triangle of asymmetric matrix */
		  
	    if(nasym>0){
	      for(k=nzs[2];k<nzs[2]+nzs[0];++k){
		au[k]=aub[k]+scal1*au[k];
	      }
	    }

	    /* damping */
		  
	    if(idamping==1){
	      scal1=gam*dtime*(1.+alpha[0]);
	      for(k=0;k<neq[0];++k){
		ad[k]+=scal1*adc[k];
	      }
	      for(k=0;k<nzs[0];++k){
		au[k]+=scal1*auc[k];
	      }
		      
	      /* upper triangle of asymmetric matrix */
		      
	      if(nasym>0){
		for(k=nzs[2];k<nzs[2]+nzs[0];++k){
		  au[k]+=scal1*auc[k];
		}
	      }
	    }

	  }
	      
	  /* thermal part */
	      
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      ad[k]=adb[k]/dtime+ad[k];
	    }
	    for(k=nzs[0];k<nzs[1];++k){
	      au[k]=aub[k]/dtime+au[k];
	    }
		  
	    /* upper triangle of asymmetric matrix */
		  
	    if(nasym>0){
	      for(k=nzs[2]+nzs[0];k<nzs[2]+nzs[1];++k){
		au[k]=aub[k]/dtime+au[k];
	      }
	    }
	  }
	}
      
	/*	for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
	for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
	for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
	for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}
	for(k=0;k<nzs[1];++k){printf("irow=%" ITGFORMAT ",%d\n",k,irow[k]);}
	for(k=0;k<neq[1]+1;++k){printf("jq=%" ITGFORMAT ",%d\n",k,jq[k]);}
	for(k=0;k<neq[1];++k){printf("icol=%" ITGFORMAT ",%d %d\n",k,icol[k],jq[k+1]-jq[k]);}*/
      

	/* Coupled dissipation control, part 1: capture the load direction.
	
	   f_hat = dR/dlambda needs no separate results() call.  At the first
	   iteration of an increment the prescribed values have just jumped
	   from lambda_n to lambda_n+dtheta and no correction has been applied
	   yet, so the residual standing in b is exactly the out-of-balance
	   that jump produced,
	
	       R_1 = -K_fp*u_hat*dtheta = f_hat*dtheta
	
	   and likewise P_1-P_n = k_pp*dtheta because du is still zero.  Both
	   fall out of the first iteration for free.
	
	   The staggered version failed here for a reason worth recording: it
	   could only move lambda once the residual was contracting, but at a
	   turning point the residual does not contract until lambda moves.
	   Solving both together removes that deadlock, because the constraint
	   becomes part of the linear system instead of a layer on top. */
	
	if((lc.diss_ctrl==2)&&(lc.diss_fhat!=NULL)){
	  if((iit==1)&&(dtheta>1.e-30)){
	    for(k=0;k<neq[1];k++) lc.diss_fhat[k]=b[k]/dtheta;
	    lc.diss_kpp=(lc.diss_p-lc.diss_pprev)/dtheta;
	    lc.diss_have=1;
	  }
	  if(lc.diss_have==1){
	    for(k=0;k<neq[1];k++) lc.diss_uf[k]=lc.diss_fhat[k];
	  }
	  if(lc.diss_probe==1){
	    printf("[DISS-GATE] it=%" ITGFORMAT " have=%" ITGFORMAT
	           " eng=%" ITGFORMAT " isolver=%" ITGFORMAT
	           " ithermal=%" ITGFORMAT " dtheta=%.4e\n",
	           iit,lc.diss_have,lc.diss_engaged,*isolver,
	           *ithermal,dtheta);
	    fflush(stdout);
	  }
	}



	/* ---- CCX_PATHFOLLOW_SOLVECHECK: keep the right-hand side -------
	   The method issues a SECOND solve with the same operator.  That is
	   only legitimate if the solver leaves ad/au intact and returns a
	   true solution both times.  Saving the rhs here lets both solves be
	   verified afterwards with one sparse mat-vec each. */

	if((pf.on==1)&&(ccxopt_getenv("CCX_PATHFOLLOW_SOLVECHECK")!=NULL)){
	  if(pf.rhs0==NULL){NNEW(pf.rhs0,double,neq[1]);NNEW(pf.y,double,neq[1]);}
	  isiz=neq[1];cpypardou(pf.rhs0,b,&isiz,&num_cpus);
	}

	/* ---- CCX_PATHFOLLOW: capture f_hat -----------------------------
	   b still holds fext-f = -R.  At the first iteration of an attempt
	   the ONLY thing that moved since the committed state is the
	   prescribed pattern, by pf.dlamjump, so -R = f_hat*pf.dlamjump and
	   the reference load vector conjugate to lambda comes out of the
	   iteration for free.  It is then frozen for the rest of the
	   increment, which is what makes the constraint exactly
	   differentiable. */

	/* ---- TOPOLOGY / RANK / RESIDUAL REPORT ------------------------
	   Runs at the first iteration of every attempt from CCX_TOPODIAG on.
	   b holds fext-f = -R and ad/au are assembled, so this is the real
	   operator and the real residual of the state the solver is about
	   to work on.  pardiso.c copies ad/au into its own array before
	   factorising, so the extra solve below cannot disturb the
	   factorisation the run then performs, and the scratch vector keeps
	   b untouched. */

	if((td_armed!=0)&&(td_from>0)&&(iinc>=td_from)&&(iit==1)&&
	   (*ithermal<2)){
	  double *td_x=NULL,*td_r0=NULL,*td_q=NULL,td_nx=0.,td_nb,td_pr;
	  ITG td_it,td_k,td_m,td_nm;
#define TD_NMODE 3

	  topodiag_run(&td_rep,td_comp,kon,ipkon,lakon,*ne,*nk,nactdof,mt,
	               nodeboun,ndirboun,*nboun,ipompc,nodempc,*nmpc,
	               ad,au,jq,irow,neq[1],nzs[0],b);
	  printf("[TOPODIAG] --- inc=%" ITGFORMAT " icutb=%" ITGFORMAT
	         " iter=%" ITGFORMAT " committed_state=%016llx ---\n",
	         iinc,icutb,iit,td_state);
	  topodiag_print(&td_rep,"assembled operator",iinc);

	  /* Softest mode by inverse iteration on the SAME operator, then
	     the residual's projection onto it.  If the residual lies along
	     a soft mode the obstruction is the topology; if it does not,
	     the topology is sound and the corrector is the problem.  That
	     is the whole point of the measurement. */

	  /* THREE soft directions, deflated against each other, not one.
	     A single inverse iteration converges to one vector, and a
	     residual can be orthogonal to that one while lying squarely in
	     a two- or three-dimensional soft SUBSPACE - which is exactly
	     what six rigid-body modes of a detached piece would look like.
	     The span projection below is the quantity that cannot be
	     fooled that way, and the random control says what "small"
	     means in 29501 dimensions. */

	  NNEW(td_x,double,neq[1]);
	  NNEW(td_r0,double,neq[1]);
	  NNEW(td_q,double,TD_NMODE*neq[1]);
	  isiz=neq[1];cpypardou(td_r0,b,&isiz,&num_cpus);
	  td_nm=0;
	  for(td_m=0;td_m<TD_NMODE;td_m++){
	    for(td_k=0;td_k<neq[1];td_k++)
	      td_x[td_k]=(double)(((td_k+7919*td_m)*2654435761U)%20011)
	                 /10005.-1.;
	    if(td_nm>0){
	      if(topodiag_deflate(td_x,td_q,td_nm,neq[1])==0) break;
	    }
	    for(td_it=0;td_it<3;td_it++){
	      td_nb=0.;
	      for(td_k=0;td_k<neq[1];td_k++) td_nb+=td_x[td_k]*td_x[td_k];
	      td_nb=sqrt(td_nb);
	      if(td_nb>0.) for(td_k=0;td_k<neq[1];td_k++) td_x[td_k]/=td_nb;
	      if(*isolver==0){
#ifdef SPOOLES
	        spooles(ad,au,adb,aub,&sigma,td_x,icol,irow,&neq[0],&nzs[0],
	                &symmetryflag,&inputformat,&nzs[2]);
#endif
	      }else if(*isolver==7){
#ifdef PARDISO
	        pardiso_main(ad,au,adb,aub,&sigma,td_x,icol,irow,&neq[0],
	                     &nzs[0],&symmetryflag,&inputformat,jq,&nzs[2],
	                     &nrhs);
#endif
	      }else{td_m=TD_NMODE;break;}
	      td_nx=0.;
	      for(td_k=0;td_k<neq[1];td_k++) td_nx+=td_x[td_k]*td_x[td_k];
	      td_nx=sqrt(td_nx);
	      if(td_nm>0) topodiag_deflate(td_x,td_q,td_nm,neq[1]);
	    }
	    if(td_m>=TD_NMODE) break;
	    if(topodiag_deflate(td_x,td_q,td_nm,neq[1])==0) break;
	    for(td_k=0;td_k<neq[1];td_k++) td_q[td_nm*neq[1]+td_k]=td_x[td_k];
	    td_nm++;
	    {
	      ITG td_bn=-1,td_bd=0,td_bc=-2,td_e,tb,tf;
	      double td_bv=0.;
	      for(td_k=0;td_k<*nk;td_k++){
	        for(td_it=1;td_it<mt;td_it++){
	          td_e=nactdof[mt*td_k+td_it];
	          if(td_e>0){
	            if(fabs(td_x[td_e-1])>td_bv){
	              td_bv=fabs(td_x[td_e-1]);td_bn=td_k+1;td_bd=td_it;
	              td_bc=td_comp[td_k];
	            }
	          }
	        }
	      }
	      topodiag_support(td_bn,kon,ipkon,lakon,*ne,&tb,&tf);
	      printf("[TOPODIAG]   soft mode %" ITGFORMAT
	             ": 1/sigma_min >= %.6e, peaks at node %" ITGFORMAT
	             " dir %" ITGFORMAT " (component %" ITGFORMAT
	             ", %" ITGFORMAT " live bulk element(s), %" ITGFORMAT
	             " live facet(s)), |cos(R,mode)|=%.3e\n",
	             td_nm,td_nx,td_bn,td_bd,td_bc,tb,tf,
	             topodiag_project(td_r0,td_x,neq[1]));
	    }
	  }
	  /* Keep the basis for the rest of the attempt: the decisive
	     question is not where the RESIDUAL points but where the
	     CORRECTION does.  With 1/sigma_min at 6e13 a residual component
	     along a soft mode of 1e-7 - which the cosine prints as zero -
	     becomes a correction of order one, and that is what a runaway
	     Newton step looks like. */

	  if((td_qkeep==NULL)||(td_qneq!=neq[1])){
	    if(td_qkeep!=NULL) SFREE(td_qkeep);
	    NNEW(td_qkeep,double,TD_NMODE*neq[1]);
	    td_qneq=neq[1];
	  }
	  isiz=TD_NMODE*neq[1];cpypardou(td_qkeep,td_q,&isiz,&num_cpus);
	  td_nmode=td_nm;

	  td_pr=topodiag_project_span(td_r0,td_q,td_nm,neq[1]);
	  for(td_k=0;td_k<neq[1];td_k++)
	    td_x[td_k]=(double)(((td_k+104729)*2246822519U)%20011)/10005.-1.;
	  printf("[TOPODIAG]   share of |R| in the span of the %" ITGFORMAT
	         " softest modes: %.6e   (random control %.6e)\n",
	         td_nm,td_pr,topodiag_project(td_r0,td_x,neq[1]));
	  {
	    ITG tb,tf;
	    topodiag_support(td_rep.resmaxnode,kon,ipkon,lakon,*ne,&tb,&tf);
	    printf("[TOPODIAG]   residual peak node %" ITGFORMAT
	           " has %" ITGFORMAT " live bulk element(s), %" ITGFORMAT
	           " live facet(s)\n",td_rep.resmaxnode,tb,tf);
	  }
	  SFREE(td_q);SFREE(td_r0);SFREE(td_x);
	  fflush(stdout);
	}

	if((pf.on==1)&&(iit==1)){

	  /* HOW FAR HAS THE FROZEN REFERENCE VECTOR DRIFTED?
	     f_hat is frozen at the first capture, and the corrector adds
	     dlambda*K^-1 f_hat to the displacement while adding dlambda to
	     lambda.  Those two are consistent only while f_hat is -dR/dlambda;
	     a scale error s makes the applied load-factor step 1/s of the one
	     the displacement correction assumed.  b/dlamjump IS the current
	     secant of -dR/dlambda, so the comparison costs one pass over the
	     vector and needs no extra residual evaluation.  Measured, not
	     assumed: on the mixed benchmark the ratio is still 0.997 at
	     increment 1500. */

	  if((pathfollow_have()==1)&&(fabs(pf.dlamjump)>1.e-8)){
	    const double *pffh=pathfollow_fhat();
	    double pfa=0.,pfb=0.,pfd=0.,pfs;
	    for(k=0;k<neq[1];k++){
	      pfs=b[k]/pf.dlamjump;
	      pfa+=pffh[k]*pffh[k];pfb+=pfs*pfs;pfd+=pffh[k]*pfs;
	    }
	    if((pfa>0.)&&(pfb>0.)){
	      pf.fhcos=pfd/sqrt(pfa*pfb);
	      pf.fhrat=sqrt(pfa/pfb);
	    }
	  }
	  pathfollow_capture(b,pf.dlamjump);
	  if(pathfollow_have()==1){
	    pathfollow_setPn(pathfollow_project(pathfollow_fhat(),vini,pf.uref,
	                                nactdof,*nk,mt));
	  }
	}

	/* Stabilisation of detached pieces.

	   A piece that has come loose by FACE connectivity still touches the
	   structure at a node or an edge, so it keeps free rotational modes
	   and the operator is singular in those directions.  At the
	   m12_field wall three independent pointers - the softest mode
	   (1/sigma_min = 230 against 3.7 elastic), the dof whose Newton
	   correction equalled the whole increment, and the largest residual -
	   all sat on nodes carrying a one-element floating piece (E-56 and
	   its follow-up).

	   Deleting those pieces was implemented and measured: it turns R0,
	   R1, R2 and the UC6 PASS disk from PASS to FAIL (E-57).  So hold
	   them instead.  Only nodes whose every surviving element is
	   unreached are touched, which leaves every load-carrying dof alone,
	   and the added stiffness is a fraction of the mean diagonal so the
	   piece is held without being welded back on.

	   DEFAULT OFF.  CCX_DAMAGE_STABILISE=<alpha> switches it on. */

	if((damage_stab_alpha>0.)&&(damage_de12_enabled)&&(*ithermal<2)){
	  ITG nstabnode=0;
	  NNEW(damage_stab_node,ITG,*nk);
	  FORTRAN(damfloatstab,(ipkon,kon,lakon,&ne0,nk,nodeboun,nboun,
				ipompc,nodempc,nmpc,damage_stab_node,
				&nstabnode));
	  /* Second criterion, and the one that matters.

	     `damfloatstab` flags nodes unreachable by FACE connectivity.  It
	     does NOT flag a node whose single supporting element is still
	     face-attached to the body but is DEAD - and that is the measured
	     case.  On m14_fine, node 776 had exactly one live element, that
	     element stood at D = 1.0000, the node travelled 9.63 mm while the
	     other three nodes of the same element moved 0.44, and the
	     element's longest edge went from 0.0974 to 9.809 - a stretch of
	     101x on a 4 mm specimen.

	     A node whose entire live support sits at g = 1-D near gmin is very
	     nearly free, whatever the connectivity says.  Flag it on the
	     degradation, which is what E-24 recorded as missing from every
	     topology predicate in this branch.  The threshold is deliberately
	     severe: the whole support must be below 1% of its stiffness. */
	  {
	    double *gsup=NULL;
	    const double *dsrc=(damage_visc_eta>0.&&damage_damvisc!=NULL)?
	      damage_damvisc:dam;
	    ITG nn,jj,ndead=0;
	    NNEW(gsup,double,*nk);
	    for(nn=0;nn<*nk;nn++) gsup[nn]=-1.;
	    for(i=0;i<ne0;i++){
	      if(ipkon[i]<0) continue;
	      if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
	      {
		double dmx=0.,g;
		ITG usedam=((damage_visc_eta>0.&&damage_damvisc!=NULL)?0:1);
		for(jj=0;jj<mi[0];jj++){
		  /* dam holds 1+D, visc holds D - see damage_de13_mark_deadsole */
		  double dd=usedam?(dsrc[mi[0]*i+jj]-1.):dsrc[mi[0]*i+jj];
		  if(dd<0.) dd=0.;
		  if(dd>1.) dd=1.;
		  if(dd>dmx) dmx=dd;
		}
		g=1.-dmx; if(g<0.) g=0.;
		for(jj=0;jj<4;jj++){
		  nn=kon[ipkon[i]+jj]-1;
		  if((nn<0)||(nn>=*nk)) continue;
		  if(g>gsup[nn]) gsup[nn]=g;
		}
	      }
	    }
	    for(nn=0;nn<*nk;nn++){
	      if(gsup[nn]<0.) continue;              /* no live support at all */
	      if(gsup[nn]>=DAMAGE_STAB_GDEAD) continue;
	      if(damage_stab_node[nn]==0){
		damage_stab_node[nn]=1;
		nstabnode++;
		ndead++;
	      }
	    }
	    SFREE(gsup);
	    if(ndead>damage_stab_maxdead){
	      damage_stab_maxdead=ndead;
	      printf("[DAMAGE STABILISE] inc=%" ITGFORMAT " nodes held on DEAD "
		     "support only: %" ITGFORMAT " (g < %.1e) (new maximum)\n",
		     iinc,ndead,DAMAGE_STAB_GDEAD);
	      fflush(stdout);
	    }
	  }

	  if(nstabnode>0){
	    double dref=0.;
	    ITG ndref=0;
	    for(k=0;k<neq[1];k++){
	      if(ad[k]>0.){dref+=ad[k];ndref++;}
	    }
	    if(ndref>0){
	      dref/=(double)ndref;
	      ITG nstabdof=0;
	      for(i=0;i<*nk;i++){
		if(damage_stab_node[i]==0) continue;
		for(idir=1;idir<=3;idir++){
		  k=nactdof[mt*i+idir];
		  if(k>0){
		    ad[k-1]+=damage_stab_alpha*dref;
		    nstabdof++;
		  }
		}
	      }
	      if(nstabdof>damage_stab_maxdof){
		damage_stab_maxdof=nstabdof;
		printf("[DAMAGE STABILISE] inc=%" ITGFORMAT " nodes=%"
		       ITGFORMAT " dof=%" ITGFORMAT " alpha=%.3e "
		       "mean_diag=%.6e (new maximum)\n",
		       iinc,nstabnode,nstabdof,damage_stab_alpha,dref);
		fflush(stdout);
	      }
	    }
	  }
	  SFREE(damage_stab_node);
	}

        /* ---- [DAMAGE REG] positive diagonal regularization, level-3
           attempts only.  Same diagonal the stabiliser above edits, same
           moment: after the whole assembly, before any solver sees it,
           and the solver frees ad right after the factorisation so
           nothing persists.  Only the DIRECTION changes - b holds the
           residual, and checkconvergence decides on ram/cam/qa/uam from
           results()/calcresidual, none of which sees the shift. */

        if((lc.reg_on==1)&&(*ithermal<2)){
          ITG rneg=0,rzero=0,nsum=0,nfl=0,nneg2=0;
          double rmin=1.e300,rmax=-1.e300,rabs,dsum=0.,dmean,dfloor;
          double ssum=0.,smax=0.,sh;
          for(k=0;k<neq[1];k++){
            if(ad[k]<0.) rneg++;
            if(ad[k]==0.) rzero++;
            if(ad[k]<rmin) rmin=ad[k];
            if(ad[k]>rmax) rmax=ad[k];
            rabs=(ad[k]<0.)?-ad[k]:ad[k];
            if(rabs>0.){dsum+=rabs;nsum++;}
          }
          dmean=(nsum>0)?dsum/nsum:0.;
          /* Cheap guard.  NaN fails every comparison, so !(dmean>0.)
             catches NaN, zero and negative alike; the second test
             catches infinity.  On any of them the shift is not applied
             and the attempt proceeds exactly as the stock one. */
          if((!(dmean>0.))||(dmean>1.e300)){
            printf("[DAMAGE REG] mean|ad| is not a usable scale "
                   "(zero, NaN or infinite); regularization NOT applied, "
                   "the attempt falls back to stock%s","\n");
            fflush(stdout);
            lc.reg_on=0;
          }else{
            dfloor=1.e-6*dmean;
            /* D must be STRICTLY positive.  Plain |ad| is not: it
               vanishes where the diagonal is zero, and for ad<0 the
               shift gives |ad|*(lambda-1), so lambda=1 lands exactly ON
               zero.  Hence the floor and a ladder that runs past 1. */
            for(k=0;k<neq[1];k++){
              rabs=(ad[k]<0.)?-ad[k]:ad[k];
              if(rabs<dfloor){rabs=dfloor;nfl++;}
              sh=lc.reg_lambda*rabs;
              ssum+=sh; if(sh>smax) smax=sh;
              ad[k]+=sh;
            }
            for(k=0;k<neq[1];k++) if(ad[k]<0.) nneg2++;
            lc.reg_napply++;
            if(lc.reg_napply==1){
              printf("[DAMAGE REG] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                     " diagonal sign census BEFORE the shift - a "
                     "DIAGNOSTIC of the diagonal, NOT evidence about the "
                     "definiteness of K: n=%" ITGFORMAT " negative=%"
                     ITGFORMAT " zero=%" ITGFORMAT " min=%.6e max=%.6e "
                     "mean|ad|=%.6e%s",
                     iinc,iit,neq[1],rneg,rzero,rmin,rmax,dmean,"\n");
              printf("[DAMAGE REG] applied ad[k] += lambda*D[k] with "
                     "D[k]=max(|ad[k]|,%.6e), lambda=%.3e (ladder step %"
                     ITGFORMAT " of %" ITGFORMAT ").  ACTUAL SHIFT: "
                     "mean=%.6e max=%.6e floored_dof=%" ITGFORMAT
                     " ; diagonal negatives %" ITGFORMAT " -> %" ITGFORMAT
                     ".  Acceptance remains on the UNMODIFIED residual%s",
                     dfloor,lc.reg_lambda,lc.reg_level+1,
                     lc.reg_nlam,ssum/neq[1],smax,nfl,rneg,nneg2,"\n");
              fflush(stdout);
            }
          }
        }

        /* ---- [DAMAGE TR] the residual, before the solver overwrites b ----
           b holds fext-f here (calcresidual.c:48-52), i.e. exactly -R with
           R=f_int-f_ext.  phi=1/2|R|^2=1/2|b|^2 is sign-blind; the GRADIENT
           is not, so the sign is carried explicitly where d=J^T b is formed.
           Nothing between the calcresidual that filled b and this point
           writes b on the mortar<=1 path. */
        dog.lc_due=0;
        if((dog.lincheck>0)&&(iinc==dog.lincheck)&&
           (iit>=dog.lc_it)&&(iit<dog.lc_it+dog.lc_nit)&&
           (dog.on==0)&&(damage_de12_enabled)&&(idamagereeq==0)&&
           (ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&
           (*idrct==0)&&(*mortar<=1)) dog.lc_due=1;
        if((dog.on==1)||(dog.lc_due==1)){
          if(dog.r0==NULL) NNEW(dog.r0,double,neq[1]);
          isiz=neq[1];cpypardou(dog.r0,b,&isiz,&num_cpus);
        }
        /* [DAMAGE CT] b holds fext-f = -R here; the solve overwrites it.
           The FULL read-before-write set is snapshotted HERE, at the current
           iterate, BEFORE the solve and before the stock results() applies
           z.  Taking it after that call mixed u_current with the history of
           u_current+z, so every finite difference was built on two
           different base states.  cam is RESTORED from this snapshot, never
           re-initialised: at this point it already holds the clean values
           the iteration top installed. */
        if(ct.on==1){
          ITG cq;
          if(ct.r0==NULL){
            NNEW(ct.r0,double,neq[1]);
            NNEW(ct.beps,double,neq[1]);
            NNEW(ct.y,double,neq[1]);
            NNEW(ct.z,double,neq[1]);
          }
          if(ct.dam==NULL){
            NNEW(ct.dam,double,mi[0]**ne);
            NNEW(ct.visc,double,mi[0]**ne);
            if(*nstate_>0) NNEW(ct.xs,double,*nstate_*mi[0]**ne);
            NNEW(ct.jac,double,12*mi[0]**ne);
            NNEW(ct.sgn,ITG,mi[0]*ne0);
          }
          isiz=neq[1];cpypardou(ct.r0,b,&isiz,&num_cpus);
          isiz=mi[0]**ne;cpypardou(ct.dam,dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(ct.visc,damage_damvisc,&isiz,&num_cpus);
          }
          if(*nstate_>0){
            isiz=*nstate_*mi[0]**ne;
            cpypardou(ct.xs,xstate,&isiz,&num_cpus);
          }
          if(damage_damjac!=NULL){
            isiz=12*mi[0]**ne;
            cpypardou(ct.jac,damage_damjac,&isiz,&num_cpus);
          }
          for(cq=0;cq<4;cq++) ct.qa[cq]=qa[cq];
          for(cq=0;cq<5;cq++) ct.cam[cq]=cam[cq];
          for(cq=0;cq<2;cq++) ct.uam[cq]=uam[cq];
          ct.lamsnap=ct.lam;
          ct.nfact++;
        }

	if(*isolver==0){
#ifdef SPOOLES
	  if(*ithermal<2){
	    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,&nzs[2]);
	    
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    spooles(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
		    &sigma,&b[neq[0]],&icol[neq[0]],iruc,
		    &n1,&n2,&symmetryflag,&inputformat,&nzs[2]);
	  }else{
	    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		    &symmetryflag,&inputformat,&nzs[2]);
	  }
#else
	  printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if((*isolver==2)||(*isolver==3)){
	  if(symmetryflag==2){
	    if(*isolver==3){
	      printf(" *WARNING in nonlingeo: the iterative Cholesky solver");
	      printf(" cannot be used for asymmetric matrices.\nThe");
	      printf(" iterative scaling solver will be used instead\n\n");
	    }
	    NNEW(rwork,double,neq[1]);
	    NNEW(sol,double,neq[1]);
	    RENEW(au,double,2*nzs[1]+neq[1]);
	    isiz=neq[1];cpypardou(&au[2*nzs[1]],ad,&isiz,&num_cpus);
	    nelt=2*nzs[1]+neq[1];
	    lrgw=131+16*neq[1];
	    isym=0;
	    NNEW(rgwk,double,lrgw);
	    NNEW(igwk,ITG,20);
	    for(i=0;i<neq[1];i++){
	      rwork[i]=1./ad[i];}
	    FORTRAN(predgmres_struct,(&neq[1],b,sol,&nelt,irow,jq,au,
				      &isym,&itol,&tol,&itmax,&iter,
				      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
				      &ligw,rwork,iwork));
	    isiz=neq[1];cpypardou(b,sol,&isiz,&num_cpus);
	    SFREE(rgwk);SFREE(igwk);SFREE(rwork);SFREE(sol);
	  }else{
	    preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
	  }
	}
	else if(*isolver==4){
#ifdef SGI
	  if(symmetryflag==2){
	    printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
	    FORTRAN(stop,());
	  }
	  token=1;
	  if(*ithermal<2){
	    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    sgi_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
		     &sigma,&b[neq[0]],&icol[neq[0]],iruc,
		     &n1,&n2,token);
	  }else{
	    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
	  }
#else
	  printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  if(symmetryflag==2){
	    printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	    FORTRAN(stop,());
	  }
	  tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
	  printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==6){
#ifdef MATRIXSTORAGE
	  matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
			ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
			lakon,kon,ne,mei,nboun,nmpc,cs,mcs,ithermal,nmethod);
	  strcpy2(fneig,jobnamec,132);
	  strcat(fneig,".frd");
	  if((f1=fopen(fneig,"ab"))==NULL){
	    printf(" *ERROR in nonlingeo: cannot open frd file for writing...");
	    exit(0);
	  }
	  fprintf(f1," 9999\n");
	  fclose(f1);
	  FORTRAN(stopwithout201,());
#else
	  printf(" *ERROR in nonlingeo: the MATRIXSTORAGE library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  if(*ithermal<2){
	    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    NNEW(jqtherm,ITG,n1+1);
	    for(i=0;i<n1+1;i++){
	      jqtherm[i]=jq[neq[0]+i]-nzs[0];}
	    pardiso_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			 &sigma,&b[neq[0]],&icol[neq[0]],iruc,
			 &n1,&n2,&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	    SFREE(jqtherm);
	  }else{
	    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }
#else
	  printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  if(*ithermal<2){
	    pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    NNEW(jqtherm,ITG,n1+1);
	    for(i=0;i<n1+1;i++){
	      jqtherm[i]=jq[neq[0]+i]-nzs[0];}
	    pastix_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			&sigma,&b[neq[0]],&icol[neq[0]],iruc,
			&n1,&n2,&symmetryflag,&inputformat,jqtherm,&nzs[2],&nrhs);
	    SFREE(jqtherm);
	  }else{
	    pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }
#else
	  printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	  

	/* Coupled dissipation control, part 2: second solve and the load
	   factor increment.
	
	       K du_R = R      (already done, du_R is in b)
	       K du_F = f_hat  (same matrix, second back substitution)
	       du     = du_R + dlambda*du_F
	
	   The constraint row for displacement control, with
	   dG = (P_n*lambda - P*lambda_n)/2 and g = dG - dtau, gives
	
	       dlambda = -[g + l_n/2*f_hat.du_R]
	                 /[l_n/2*f_hat.du_F + P_n/2 - l_n/2*k_pp]
	
	   PARDISO is called a second time on the same ad/au, so the
	   factorisation is repeated; correctness first, and the symbolic
	   analysis is already cached when CCX_PARDISO_REUSE_SYMBOLIC is set. */
	
	/* The constraint is unsatisfiable while the response is still
	   elastic: dG is identically zero there, so g = -dtau every
	   iteration and lambda is driven at the clamp forever.  That is
	   what the first coupled run did - residual 0.088 against an
	   average force of 0.023, diverging at increment 2.  Stay in
	   displacement control until the measured dissipation of an
	   accepted increment reaches a fraction of the target, then
	   engage and stay engaged. */

	if((lc.diss_ctrl==2)&&(lc.diss_have==1)&&
	   (lc.diss_engaged==1)&&
	   (*isolver==7)&&(*ithermal<2)){
#ifdef PARDISO
	  pardiso_main(ad,au,adb,aub,&sigma,lc.diss_uf,icol,irow,
	               &neq[0],&nzs[0],&symmetryflag,&inputformat,jq,
	               &nzs[2],&nrhs);
#endif
	  lc.diss_fr=0.;
	  lc.diss_ff=0.;
	  for(k=0;k<neq[1];k++){
	    lc.diss_fr+=lc.diss_fhat[k]*b[k];
	    lc.diss_ff+=lc.diss_fhat[k]*lc.diss_uf[k];
	  }
	  lc.diss_dgcur=0.5*(lc.diss_pprev*lc.diss_lamcur-
	                         lc.diss_p*lc.diss_lprev);
	  lc.diss_g=lc.diss_dgcur-lc.diss_target;
	  /* P is a function of u alone, so dg/dlambda is exactly 0.5*P_n and
	     every remaining lambda dependence already travels through du_F.
	     The earlier -0.5*lambda_n*k_pp counted dP/dlambda a second time,
	     which is what left the correction cancelling itself: the residual
	     sat at 0.043639 while the displacement correction was 5.7e-7. */
	  lc.diss_den=0.5*lc.diss_lprev*lc.diss_ff
	                 +0.5*lc.diss_pprev;
	  if(lc.diss_probe==1){
	    printf("[DISS-PROBE] it=%" ITGFORMAT " lam=%.6f dG=%.4e g=%.4e"
	           " fr=%.4e ff=%.4e den=%.4e\n",
	           iit,lc.diss_lamcur,lc.diss_dgcur,lc.diss_g,
	           lc.diss_fr,lc.diss_ff,lc.diss_den);
	    fflush(stdout);
	  }
	  if(fabs(lc.diss_den)>1.e-30){
	    lc.diss_dlam=-(lc.diss_g
	                       +0.5*lc.diss_lprev*lc.diss_fr)
	                     /lc.diss_den;
	    if(lc.diss_dlam>DAMAGE_DISS_DLAM*dthetaref)
	      lc.diss_dlam=DAMAGE_DISS_DLAM*dthetaref;
	    if(lc.diss_dlam<-DAMAGE_DISS_DLAM*dthetaref)
	      lc.diss_dlam=-DAMAGE_DISS_DLAM*dthetaref;
	    for(k=0;k<neq[1];k++)
	      b[k]+=lc.diss_dlam*lc.diss_uf[k];
	    lc.diss_lamcur+=lc.diss_dlam;
	    if(lc.diss_probe==1){
	      printf("[DISS-PROBE]      dlam=%.6e -> lam=%.6f\n",
	             lc.diss_dlam,lc.diss_lamcur);
	      fflush(stdout);
	    }
	    if(lc.arc==1){
	      if(lc.diss_lamcur<0.) lc.diss_lamcur=0.;
	    }else if(lc.diss_lamcur<=theta){
	      lc.diss_lamcur=theta+(*tmin);
	    }
	    if(lc.diss_lamcur>1.) lc.diss_lamcur=1.;
	    for(k=0;k<*nboun;k++){
	      xbounact[k]=xbounold[k]+
	        (xboun[k]-xbounold[k])*lc.diss_lamcur;
	    }
	  }
	}

	/* ---- CCX_PATHFOLLOW: close the bordered system -----------------
	   b now holds du_R.  Solve once more with the SAME operator for
	   du_F = K^-1 f_hat and close the second row.  The extra
	   factorisation is the honest price of staying solver agnostic: the
	   CalculiX solver wrappers factor and solve in one call, so a second
	   right-hand side cannot reuse the first factorisation without
	   changing them.  Correctness first; the cost is a constant factor
	   and is reported at the end of the step. */

	/* ---- where does the CORRECTION point? --------------------------
	   b now holds the ordinary Newton correction du_R = K^-1(-R), before
	   any constraint, line search or trust region touches it.  If its
	   norm is dominated by the soft subspace measured at iteration 1,
	   the runaway steps are the ill-conditioning being amplified, and no
	   amount of step-size control fixes that.  If it is not, the steps
	   are large for a constitutive reason. */

	if((td_armed!=0)&&(td_from>0)&&(iinc>=td_from)&&(td_qkeep!=NULL)&&
	   (td_qneq==neq[1])&&(td_nmode>0)&&(*ithermal<2)){
	  double tdn=0.,tds;
	  ITG tdk;
	  for(tdk=0;tdk<neq[1];tdk++) tdn+=b[tdk]*b[tdk];
	  tdn=sqrt(tdn);
	  tds=topodiag_project_span(b,td_qkeep,td_nmode,neq[1]);
	  printf("[TOPODIAG]   correction inc=%" ITGFORMAT " iter=%" ITGFORMAT
	         " |du|2=%.6e  share in the %" ITGFORMAT
	         " softest modes: %.6e\n",iinc,iit,tdn,td_nmode,tds);
	  fflush(stdout);
	}

	/* Reset here, not inside the block below: the line search has to be
	   able to tell "this iteration applied a constraint step" from "the
	   value left over by an earlier one". */

	pf.applied=0;pf.dlamit=0.;pf.lam0it=pf.lam;

	if((pf.on==1)&&(pathfollow_have()==1)&&(*ithermal<2)){
	  pathdrv_predictor(&pf,&damage_glob,&nlgt,xboun,xbounold,
	                    ad,au,icol,isolver,sigma,inputformat,nrhs,
	                    symmetryflag,iit);
	}



	/* Locate the softest mode of the assembled operator.

	   At the stall an unloading of 0.5% moves the structure by 2e-2 in a
	   corner far from the crack, and |du|/|r| has grown three orders from
	   its mid-run value.  That is a soft mode, and because it survives a
	   change of load it belongs to the topology, not to lambda - which is
	   why walking lambda downwards does not rescue the run.

	   Inverse iteration finds it: x <- K^-1 x, renormalised, converges to
	   the smallest-singular-value direction, and ||x||/||b|| after the
	   first solve estimates 1/sigma_min.  Reporting where that vector
	   lives says which part of the model is the mechanism, instead of
	   guessing at candidates one at a time.

	   The factorisation is already in PARDISO's hands, so the extra
	   solves are back-substitutions. */

	/* [WALLDIAG] Each deflated inverse iteration is a full factorisation
	   and solve.  Armed by hand on ONE named increment that is affordable;
	   armed by the load-factor gate it would run on every iteration of
	   every increment from the wall on, so under that gate it is taken
	   once per attempt.  CCX_DAMAGE_NULLVEC on its own is unchanged. */
	if((prb.null_inc>0)&&(iinc>=prb.null_inc)&&
	   ((prb.wall_armed==0)||(iit==1))&&
	   (*isolver==7)&&(*ithermal<2)&&(*mortar<=1)){
	  NNEW(prb.null_x,double,neq[1]);
	  for(k=0;k<neq[1];k++){
	    prb.null_seed=(1103515245*prb.null_seed+12345)&0x7fffffff;
	    prb.null_x[k]=(double)(prb.null_seed%20001-10000)/10000.;
	  }
	  for(prb.null_it=0;prb.null_it<prb.null_nit;prb.null_it++){
	    prb.null_nb=0.;
	    for(k=0;k<neq[1];k++)
	      prb.null_nb+=prb.null_x[k]*prb.null_x[k];
	    prb.null_nb=sqrt(prb.null_nb);
	    if(prb.null_nb>0.)
	      for(k=0;k<neq[1];k++) prb.null_x[k]/=prb.null_nb;
#ifdef PARDISO
	    pardiso_main(ad,au,adb,aub,&sigma,prb.null_x,icol,irow,
			 &neq[0],&nzs[0],&symmetryflag,&inputformat,jq,
			 &nzs[2],&nrhs);
#endif
	    prb.null_nx=0.;
	    for(k=0;k<neq[1];k++)
	      prb.null_nx+=prb.null_x[k]*prb.null_x[k];
	    prb.null_nx=sqrt(prb.null_nx);
	    if((prb.null_it<3)||(prb.null_it==prb.null_nit-1))
	    printf("[DAMAGE NULLVEC] inc=%" ITGFORMAT " iteration %" ITGFORMAT
		   " 1/sigma_min >= %.6e\n",iinc,prb.null_it+1,
		   prb.null_nx);
	  }
	  /* where does the mode live */
	  prb.null_amax=0.;
	  for(i=0;i<*nk;i++){
	    for(idir=1;idir<=3;idir++){
	      k=nactdof[mt*i+idir];
	      if(k<=0) continue;
	      if(fabs(prb.null_x[k-1])>prb.null_amax)
		prb.null_amax=fabs(prb.null_x[k-1]);
	    }
	  }
	  printf("[DAMAGE NULLVEC] nodes carrying the mode "
		 "(|component| > 0.2 of the maximum):\n");
	  prb.null_cnt=0;
	  for(i=0;i<*nk;i++){
	    prb.null_nn=0.;
	    for(idir=1;idir<=3;idir++){
	      k=nactdof[mt*i+idir];
	      if(k<=0) continue;
	      if(fabs(prb.null_x[k-1])>prb.null_nn)
		prb.null_nn=fabs(prb.null_x[k-1]);
	    }
	    if((prb.null_amax>0.)&&(prb.null_nn>0.2*prb.null_amax)){
	      prb.null_cnt++;
	      if(prb.null_cnt<=20){
		printf("[DAMAGE NULLVEC]   node %" ITGFORMAT
		       " at (%.3f, %.3f, %.3f) amplitude %.4f\n",
		       i+1,co[3*i],co[3*i+1],co[3*i+2],
		       prb.null_nn/prb.null_amax);
	      }
	    }
	  }
	  printf("[DAMAGE NULLVEC] %" ITGFORMAT " node(s) above the "
		 "threshold out of %" ITGFORMAT "\n",prb.null_cnt,*nk);
	  fflush(stdout);
	  SFREE(prb.null_x);
	  FORTRAN(stopwithout201,());
	}

	if(*mortar<=1){
	  if(isensitivity){
	    SFREE(adcpy);MNEW(adcpy,double,neq[1]);
	    SFREE(aucpy);MNEW(aucpy,double,(nasym+1)*nzs[1]);
	    isiz=neq[1];cpypardou(adcpy,ad,&isiz,&num_cpus);
	    isiz=(nasym+1)*nzs[1];cpypardou(aucpy,au,&isiz,&num_cpus);
	  }
	  /* ---- [DAMAGE TR] everything the dogleg needs from J, taken while
	     J is still allocated.  Two lines below, ad and au are freed, and
	     every trial residual evaluation happens after that - so a dogleg
	     that wanted a matrix-vector product per trial could not have one.
	     It does not need one:

	         J p_N = r0        (exactly: p_N is what the solver returned)
	         J d   = w         (computed here, once)

	     and every step this trust region can propose is p = pa*d + pb*p_N,
	     so J p = pa*w + pb*r0 is a two-scalar combination of vectors that
	     already exist.  No matrix, no second factorisation and no second
	     solve for the whole trial loop.

	     Storage (add_sm_st_as.f:29-56, mastruct.c:795-820, cross-checked
	     against opas.f:34-46): jq and irow are 1-based; column c (0-based)
	     owns au[jq[c]-1 .. jq[c+1]-2]; for slot k in that column with
	     r=irow[k]-1 > c,  au[k] = J(r,c)  and  au[nzs[2]+k] = J(c,r);
	     ad[i] = J(i,i).  sigma is 0 everywhere in nonlingeo, so the
	     factorised operator is (ad,au) with no shift.

	     THE TRANSPOSE IS PROVED, NOT ASSERTED.  d = J^T r0 and J p_N = r0
	     give dot(d,p_N) = r0^T J J^-1 r0 = |r0|^2 identically.  Swap the
	     two halves of au, or use J where J^T is meant, and the identity
	     fails at once.  It is checked on every armed iteration and the
	     mechanism REFUSES TO ARM when it is off by more than 1e-3. */
	  if(((dog.on==1)||(dog.lc_due==1))&&
	     (dog.r0!=NULL)&&(*mortar<=1)){
	    if((nasym!=1)||(symmetryflag!=2)||(*ithermal>=2)||
	       (nzs[2]!=nzs[1])||(neq[0]!=neq[1])){
	      if(dog.have>=0){
	        printf("[DAMAGE TR] REFUSING TO ARM: the assembled operator is "
	               "not the one this construction was proved on (nasym=%"
	               ITGFORMAT " symmetryflag=%" ITGFORMAT " ithermal=%"
	               ITGFORMAT " neq0=%" ITGFORMAT " neq1=%" ITGFORMAT
	               " nzs1=%" ITGFORMAT " nzs2=%" ITGFORMAT ").  The wall "
	               "goes to the original stock stop.%s",
	               nasym,symmetryflag,*ithermal,neq[0],neq[1],
	               nzs[1],nzs[2],"\n");
	        fflush(stdout);
	      }
	      dog.have=(dog.on==1)?-1:dog.have;
	      dog.on=0;dog.lc_due=0;
	    }else{
	      ITG dk,dc,dr;
	      double daden,dad;
	      if(dog.d==NULL){
	        NNEW(dog.d,double,neq[1]);
	        NNEW(dog.w,double,neq[1]);
	        NNEW(dog.pn,double,neq[1]);
	      }
	      /* d = J^T r0 */
	      for(dk=0;dk<neq[1];dk++)
	        dog.d[dk]=ad[dk]*dog.r0[dk];
	      for(dc=0;dc<neq[1];dc++){
	        for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	          dr=irow[dk]-1;
	          dog.d[dc]+=au[dk]*dog.r0[dr];
	          dog.d[dr]+=au[nzs[2]+dk]*dog.r0[dc];
	        }
	      }
	      /* w = J d   (the same loop with the two halves exchanged) */
	      for(dk=0;dk<neq[1];dk++)
	        dog.w[dk]=ad[dk]*dog.d[dk];
	      for(dc=0;dc<neq[1];dc++){
	        for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	          dr=irow[dk]-1;
	          dog.w[dr]+=au[dk]*dog.d[dc];
	          dog.w[dc]+=au[nzs[2]+dk]*dog.d[dr];
	        }
	      }
	      isiz=neq[1];cpypardou(dog.pn,b,&isiz,&num_cpus);
	      /* [WALLDIAG] MASKSTEP.  pm = p_N with every component that
	         belongs to an AUTOSPC-masked node zeroed, and wm = J pm by the
	         same proved loop.  Built HERE because ad and au are freed
	         before the linearisation check runs; pm and wm are the only
	         things that survive to it.  Nothing on a solution path reads
	         either. */
	      if((prb.wall_maskstep!=0)&&(damage_spc_mask!=NULL)&&
	         (damage_spc_nk>=*nk)){
	        ITG mi_,mj_,mk_;
	        if(dog.pm==NULL){
	          NNEW(dog.pm,double,neq[1]);
	          NNEW(dog.wm,double,neq[1]);
	        }
	        for(dk=0;dk<neq[1];dk++) dog.pm[dk]=dog.pn[dk];
	        for(mi_=0;mi_<*nk;mi_++){
	          if(damage_spc_mask[mi_]==0) continue;
	          for(mj_=1;mj_<mt;mj_++){
	            mk_=nactdof[mt*mi_+mj_];
	            if(mk_>0) dog.pm[mk_-1]=0.;
	          }
	        }
	        for(dk=0;dk<neq[1];dk++)
	          dog.wm[dk]=ad[dk]*dog.pm[dk];
	        for(dc=0;dc<neq[1];dc++){
	          for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	            dr=irow[dk]-1;
	            dog.wm[dr]+=au[dk]*dog.pm[dc];
	            dog.wm[dc]+=au[nzs[2]+dk]*dog.pm[dr];
	          }
	        }
	        dog.npm2=0.;
	        for(dk=0;dk<neq[1];dk++)
	          dog.npm2+=dog.pm[dk]*dog.pm[dk];
	      }
	      dog.nb2=0.;dog.nd2=0.;dog.nw2=0.;
	      dog.npn2=0.;dog.dtpn=0.;
	      for(dk=0;dk<neq[1];dk++){
	        dog.nb2+=dog.r0[dk]*dog.r0[dk];
	        dog.nd2+=dog.d[dk]*dog.d[dk];
	        dog.nw2+=dog.w[dk]*dog.w[dk];
	        dog.npn2+=dog.pn[dk]*dog.pn[dk];
	        dog.dtpn+=dog.d[dk]*dog.pn[dk];
	      }
	      /* how far from symmetric the operator actually IS - measured,
	         because the dogleg is only worth building if J^T != J */
	      dog.asym=0.;daden=0.;
	      for(dk=0;dk<nzs[1];dk++){
	        dad=fabs(au[dk]-au[nzs[2]+dk]);
	        if(dad>dog.asym) dog.asym=dad;
	        if(fabs(au[dk])>daden) daden=fabs(au[dk]);
	      }
	      if(daden>0.) dog.asym/=daden;
	      dog.ident=(dog.nb2>0.)?
	                       dog.dtpn/dog.nb2:0.;
	      if((dog.nb2<=0.)||(dog.nw2<=0.)||
	         (dog.nd2<=0.)||(fabs(dog.ident-1.)>1.e-3)){
	        printf("[DAMAGE TR] TRANSPOSE-CHECK FAILED at inc=%" ITGFORMAT
	               " iter=%" ITGFORMAT ": dot(J^T R,p_N)/|R|^2 = %.12e, "
	               "and it must be 1.  Either the transpose is not a "
	               "transpose or the solve is not a solve.  The dogleg "
	               "DISARMS and the wall goes to the original stock stop.%s",
	               iinc,iit,dog.ident,"\n");
	        fflush(stdout);
	        dog.have=-1;dog.on=0;dog.lc_due=0;
	      }else{
	        dog.tc=dog.nd2/dog.nw2;
	        dog.have=1;
	        if(dog.on==1) dog.nfact++;
	        if(dog.banner==0){
	          dog.banner=1;
	          printf("[DAMAGE TR] TRANSPOSE-CHECK inc=%" ITGFORMAT " iter=%"
	                 ITGFORMAT " dot(J^T R,p_N)/|R|^2 = %.12e (exact value "
	                 "1; a J used in place of J^T does not satisfy it).  "
	                 "OPERATOR ASYMMETRY max|au_L-au_U|/max|au_L| = %.6e, "
	                 "so J^T R is genuinely not J R here.%s",
	                 iinc,iit,dog.ident,dog.asym,"\n");
	          fflush(stdout);
	        }
	      }
	    }
	  }
	  SFREE(ad);SFREE(au);
	} 
      }
      
      /* explicit dynamic step */
      
      else{
	if(((mscalmethod==0)||(mscalmethod==2))&&(*mortar!=-1)){
	  //    for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	  if(*ithermal!=2){
	    isiz=neq[0];divparll(b,adb,&isiz,&num_cpus);
	  }
	  //      for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      b[k]=b[k]*dtime/adb[k];
	    }
	  }
	}
	else{
	  if(*ithermal!=2){
	    if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(b,&neq[0]);
#endif
	    }
	    else if(*isolver==4){
#ifdef SGI
	      sgi_solve(b,token);
#endif
	    }
	    else if(*isolver==5){
#ifdef TAUCS
	      tau_solve(b,&neq[0]);
#endif
	    }
	    else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(b,&neq[0],&symmetryflag,&inputformat,&nrhs);
#endif
	    }
	    else if(*isolver==8){
#ifdef PASTIX
	      pastix_solve(b,&neq[0],&symmetryflag,&nrhs);
#endif
	    }
	    if(*mortar==-1){
	      if(ncont!=0){
		if(iinc==1){
		  for(i=0;i<neqtot;i++){
		    k=floor(ltot[i]/10);
		    l=ltot[i]-10*k;
		    b[ktot[i]-1]=veold[mt*(k-1)+l];
		  }
		} else{

		  /* determine the velocity in the contact nodes */
	      
		  for(i=0;i<neqtot;++i){
		    b[ktot[i]-1]=(qb[i]-volddof[ktot[i]-1])/(dtime);
		  }
		}
		SFREE(qb);
	      }
	      SFREE(volddof);
	    }
	  }
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      b[k]=b[k]*dtime/adb[k];
	    }
	  }
	}
      }
      
      /* mortar */

      if(*mortar>1){	    
  
	/* restoring the structure of the original stiffness
	   matrix */

	for(i=0;i<3;i++){
	  nzs[i]=nzstemp[i];}
	for (i=0;i<neq[1];i++){jq[i]=jqtemp[i];icol[i]=icoltemp[i];}
	jq[neq[1]]=jqtemp[neq[1]];
	for (i=0;i<nzs[1];i++){irow[i]=irowtemp[i];}
	SFREE(jqtemp);SFREE(irowtemp);SFREE(icoltemp);

	/* trafo util->u , calculate cstress and update active set  */

	stressmortar(bhat,adc2,auc2,jqc2,irowc2,neq,gap,b,islavact,irowddinv,
		     jqddinv,auddinv,irowt,jqt,aut,irowtinv,
		     jqtinv,autinv,ntie,nslavnode,islavnode,nmastnode,
		     imastnode,slavnor,slavtan,nactdof,&iflagact,cstress,
		     cstressini,mi,cdisp,f_cs,f_cm,&iit,&iinc,vold,vini,bp,
		     nk,nboun,ndirboun,nodeboun,xboun,
		     nmpc,
		     ipompc,nodempc,coefmpc,
		     nslavmpc,islavmpc,
		     tieset,elcon,tietol,ncmat_,ntmat_,plicon,nplicon,npmat_,
		     nelcon,&dtime,cfs,cfm,islavnodeinv,Bd,irowb,jqb,Dd,
		     irowd,jqd,Ddtil,irowdtil,jqdtil,Bdtil,irowbtil,jqbtil,
		     nmethod,&bet,ithermal,
		     iperturb,labmpc,cam,veold,accold,&gam,
		     cfsini,cfstil,plkcon,nplkcon,filab,f,fn,qa,nprint,prlab,
		     xforc,nforc,iponoel);
	  
	SFREE(auc2);SFREE(adc2);SFREE(irowc2);SFREE(icolc2);SFREE(jqc2);
	SFREE(au);SFREE(ad);	  
      }

      /* calculating the displacements, stresses and forces */
      
      MNEW(v,double,mt**nk);
      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
      
      NNEW(stx,double,6*mi[0]**ne);
      MNEW(fn,double,mt**nk);

      /* for massless explicit dynamics without energy
         calculation only the displacements have to be calculated

         this does not work for explicit dynamics without massless
         contact since in that case f and fini have to be calculated  */
      
      if((*mortar==-1)&&(masslesslinear>0)&&(*nener==0)){
	resultsini(nk,v,ithermal,filab,iperturb,f,fn,
		   nactdof,&iout,qa,vold,b,nodeboun,ndirboun,
		   xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,
		   neq,veold,accold,&bet,&gam,&dtime,mi,vini,nprint,prlab,
		   &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
		   &ikin,&intpointvart,typeboun,&num_cpus,mortar,nener,iponoeln,
		   network);
      }else{
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	trial_results(&nlgt);
	if(ne1d2d==1)SFREE(inum);
      }

      /* implicit dynamics (Matteo Pacher) */

      if((*ne!=ne0)&&(*nmethod==4)&&(*ithermal<2)&&(*iexpl<=1)){
	FORTRAN(storecontactprop,(ne,&ne0,lakon,kon,ipkon,mi,ielmat,elcon,
				  mortar,adblump,nactdof,springarea,ncmat_,
				  ntmat_,stx,&temax));
      }

      /* Reaction conjugate to the prescribed pattern u_p = lambda*xboun.
         Taken on the main Newton path, where fn is the array just filled
         by results() and resultsini.c has set calcul_fn=1 for NLGEOM.
         Two earlier placements were wrong and both showed it plainly:
         after SFREE(fn) the sum was exactly zero, and at the output calls
         it was stale enough to give a negative dissipation increment. */

      if(lc.diss_report==1){
        lc.diss_p=0.;
        for(i=0;i<*nboun;i++){
          if((ndirboun[i]<1)||(ndirboun[i]>mi[1])) continue;
          lc.diss_p+=fn[mt*(nodeboun[i]-1)+ndirboun[i]]*xboun[i];
        }
      }

      opcheck_probe(&opd,&nlgt,ndmat_,damage_damcat,iit);


      /* Dissipation-controlled load factor.

         Shrinking the step cannot pass a limit point: with the target
         tightened threefold the stall moved by 0.0006 in step time, so at
         that lambda no equilibrium exists on the branch being tracked and
         the load factor itself has to become an unknown.

         This is the staggered form of the constraint rather than a fully
         coupled bordered solve.  Each iteration still solves K du = R at
         the current lambda, then moves lambda to drive

             g = dG - dtau = 0,    dG = (P_n*lambda - P*lambda_n)/2

         with the slope taken by secant from the previous iteration.  It
         converges more slowly than the coupled form, but it needs no
         second right-hand side and no f_hat vector, and it has the one
         property that matters here: lambda is free to decrease.

         xbounact is rebuilt directly, because tempload composes it as
         xbounold + (xboun-xbounold)*reltime and resultsini.c writes that
         value straight into v.  theta stays untouched until the increment
         is accepted; only then is dtheta set so checkconvergence lands on
         the lambda that was actually reached. */

      /* The constraint may only be driven by an equilibrium reaction.  In
         the coupled form that is automatic, because lambda and u are
         solved together and both hold at convergence.  Staggered, the
         early iterations carry a reaction that is not in equilibrium yet,
         and driving lambda with it produced a negative dG - impossible for
         a dissipation - and threw the load factor about near the turning
         point.  So lambda is held until the residual has contracted for
         two iterations, and a non-positive dG is treated as "not yet
         meaningful" rather than as a constraint violation. */

      lc.diss_ok=0;
      if((lc.diss_ctrl==1)&&(iit>=3)&&(lc.diss_init==1)&&
         (ram[0]<ram1[0])&&(ram1[0]<ram2[0])){
        lc.diss_ok=1;
      }

      if(lc.diss_ok==1){
        lc.diss_dgcur=0.5*(lc.diss_pprev*lc.diss_lamcur-
                               lc.diss_p*lc.diss_lprev);
        if(lc.diss_dgcur<=0.) lc.diss_ok=0;
      }

      if(lc.diss_ok==1){
        lc.diss_g=lc.diss_dgcur-lc.diss_target;

        if(lc.diss_have==1){
          lc.diss_slope=(lc.diss_dgcur-lc.diss_dgold)/
            (lc.diss_lamcur-lc.diss_lamold);
        }else{
          /* first pass: the only slope estimate available is the secant
             through the converged state */
          lc.diss_slope=0.5*lc.diss_pprev;
        }
        if(fabs(lc.diss_slope)<1.e-30) lc.diss_slope=
          (lc.diss_slope<0.)?-1.e-30:1.e-30;

        lc.diss_dgold=lc.diss_dgcur;
        lc.diss_lamold=lc.diss_lamcur;
        lc.diss_have=1;

        lc.diss_dlam=-lc.diss_g/lc.diss_slope;

        /* bound the move so a bad secant cannot throw lambda across the
           whole step */
        if(lc.diss_dlam>DAMAGE_DISS_DLAM*dthetaref)
          lc.diss_dlam=DAMAGE_DISS_DLAM*dthetaref;
        if(lc.diss_dlam<-DAMAGE_DISS_DLAM*dthetaref)
          lc.diss_dlam=-DAMAGE_DISS_DLAM*dthetaref;

        lc.diss_lamcur+=lc.diss_dlam;
        if(lc.arc==1){
          if(lc.diss_lamcur<0.) lc.diss_lamcur=0.;
        }else if(lc.diss_lamcur<theta){
          lc.diss_lamcur=theta;
        }
        if(lc.diss_lamcur>1.) lc.diss_lamcur=1.;

        for(i=0;i<*nboun;i++){
          xbounact[i]=xbounold[i]+
            (xboun[i]-xbounold[i])*lc.diss_lamcur;
        }
      }

      /* updating the external work (only for dynamic calculations) */

      if((*nmethod==4)&&(*ithermal<2)&&(*nener==1)){
	allwk=allwkini;
	worparll(&allwk,fnext,&mt,fnextini,v,vini,nk,&num_cpus);

        /* Work due to damping forces (cv and cvini) --> MPADD */

	if(idamping==1){
	  dampwk=dampwkini;
	  dam1parll(&mt,nactdof,aux2,v,vini,nk,&num_cpus);
	  dam2parll(&dampwk,cv,cvini,aux2,&neq[0],&num_cpus);
	}
        /* Damping forces --> MPADD */
      }

      /* line search (only for static surface-to-surface penalty contact)
         and not in the first iteration */

      if((*mortar==1)&&(iit!=1)&&(*ne-ne0>0)&&(*nmethod!=4)){

	SFREE(v);SFREE(stx);SFREE(fn);
      
	/* calculating the residual */
      
	NNEW(res,double,neq[1]);
	trial_reduce(&nlgt,res);

	/* calculating the line search factor */

	sum1=0.;sum2=0.;
	for(i=0;i<neq[1];i++){
	  sum1+=b[i]*resold[i];
	  sum2+=b[i]*res[i];
	}
	SFREE(res);

	if(fabs(sum1-sum2)<1.e-30){
	  flinesearch=1.;
	}else{
	  flinesearch=sum1/(sum1-sum2);
	  if(flinesearch>smaxls){
	    flinesearch=smaxls;
	  }else if(flinesearch<sminls){
	    flinesearch=sminls;
	  }
	}
	printf("line search factor=%f\n\n",flinesearch);

	/* update the solution */

	for(i=0;i<neq[1];i++){
	  b[i]*=flinesearch;}
      
	MNEW(v,double,mt**nk);
	isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
	  
	NNEW(stx,double,6*mi[0]**ne);
	MNEW(fn,double,mt**nk);
	  
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	trial_results(&nlgt);
	if(ne1d2d==1)SFREE(inum);

	/* Structural finite-difference check of the ASSEMBLED Jacobian.

	   mpfd validates the 6x6 material tangent and says it is right to
	   1e-3 even on a shear path.  Nothing validates what the element
	   loop and add_sm_st_as then make of it, and that is the only step
	   left between a correct constitutive law and the operator PARDISO
	   actually factorises.

	   K(i,j) = d fn_i / d u_j is measured by central differences on the
	   internal force, one column per perturbed degree of freedom, and
	   compared against the coefficient read back out of the sparse
	   structure.  Columns are taken at the nodes of the most damaged
	   element, because a tangent that is right in the elastic bulk and
	   wrong in the process zone is exactly the failure that would
	   survive every test written so far.

	   Storage convention, from add_sm_st_as.f: ad(i) is the diagonal;
	   K(i,j) with i>j sits in column j of au; with i<j it sits in
	   column i offset by nzs(3).

	   The probe perturbs v and calls results repeatedly, so the run is
	   stopped as soon as it reports.  This is a diagnostic build, not a
	   production one. */

      }

      /* BK3 adaptive damage line search.

         Evaluate the full Newton trial first.  A line search is admissible
         only in a contact-free, static physical increment, only while
         progressive damage is actually increasing from the committed
         baseline, and only when the mechanical infinity norm of the residual
         grows by more than DAMAGE_LINESEARCH_GROWTH.  Contracting Newton
         iterations therefore pay no extra constitutive evaluation and are
         never damped.  Terminal same-load topology solves stay on NC2. */



      /* ---- [DAMAGE TR LINCHECK] the sign convention, measured -----------

         Everything the trust region does rests on three conventions that are
         easy to get backwards and impossible to see in the output:

             b (before the solve) = fext - f = -R,   R = f_int - f_ext
             p_N = +b (after the solve),  so  J p_N = r0 = -R
             d   = J^T r0 = -g,  the steepest-DESCENT direction for phi

         This block measures them instead of asserting them.  On a healthy
         increment - the state must be smooth for the test to mean anything -
         it walks eps down a ladder and compares the residual actually
         returned by results()/calcresidual() against the linear model the
         dogleg uses:

             R(u + eps*p)  ->  R(u) + eps*J p        i.e.
             res(u+eps*p)  ->  r0 - eps*(J p)        (res = -R throughout)

         DEFECT = |res(u+eps p) - (r0 - eps Jp)| / (eps |Jp|)  must fall like
         O(eps).  It is run twice: once with p = p_N, where J p = r0 exactly
         and the ratio |res|/((1-eps)|r0|) must go to 1 - this pins the sign
         of b, of R and of the Newton correction together - and once with
         p = d scaled to |p_N|, where J p = w.  The second pass is the one
         that pins J^T: d came out of the transpose loop and w out of the
         forward loop, so a swapped pair cannot pass both.

         It ends by evaluating p_N, which is the state the unprobed code
         would have had, and puts cam/qa/uam back first. */

      if((dog.lc_due==1)&&(dog.have==1)){
        /* [WALLDIAG] The ladder has to reach BELOW the step the search
           actually takes, or it cannot tell a consistent tangent with a
           small radius of validity from an inconsistent one: the second
           wall's best rung is alpha~0.004, and the old ladder stopped at
           0.0156.  1/2^k down to k=14 puts six rungs below 0.004 while the
           defect is still far above the cancellation floor, which for
           |J p|=|r0| sits at about 1e-16/eps. */
        static const double lcE[15]={1.,0.5,0.25,0.125,0.0625,0.03125,
                                     0.015625,0.0078125,0.00390625,
                                     0.001953125,0.0009765625,
                                     0.00048828125,0.000244140625,
                                     0.0001220703125,0.00006103515625};
        ITG lnst,lii,ljj,lpass,lact,lspc,lmpcd;
        double le,lsc,lnum,lden,lnjp,lnres,lmv,ldv,linf,linf0;
        double lqas[4],luams[2];
        double *lp=NULL,*ljp=NULL;

        lnst=*nstate_;
        if(dog.res==NULL){
          NNEW(dog.res,double,neq[1]);
          NNEW(dog.dam,double,mi[0]**ne);
          NNEW(dog.visc,double,mi[0]**ne);
          if(lnst>0) NNEW(dog.xs,double,lnst*mi[0]**ne);
        }
        NNEW(lp,double,neq[1]);
        NNEW(ljp,double,neq[1]);

        lact=0;lspc=0;lmpcd=0;
        for(ljj=0;ljj<*nk;ljj++){
          for(lii=1;lii<mt;lii++){
            if(nactdof[mt*ljj+lii]>0) lact++;
            else if(nactdof[mt*ljj+lii]==0) lspc++;
            else lmpcd++;   /* negative: MPC-dependent, and in this tree
                               SPC-eliminated slots land here as well */
          }
        }
        printf("[DAMAGE TR LINCHECK] inc=%" ITGFORMAT " iter=%" ITGFORMAT
               " REDUCED SPACE: neq[0]=%" ITGFORMAT " neq[1]=%" ITGFORMAT
               "; mechanical slots %" ITGFORMAT " = active (nactdof>0) %"
               ITGFORMAT " + excluded, nactdof==0 %" ITGFORMAT
               " + excluded, nactdof<0 %" ITGFORMAT "; nboun=%" ITGFORMAT
               " nmpc=%" ITGFORMAT " nzs[1]=%" ITGFORMAT " nzs[2]=%" ITGFORMAT
               ".  active == neq[1] and excluded == nboun+MPC-dependent is "
               "the proof that b, f, fext, ad and au share ONE post-SPC/MPC "
               "numbering (nactdof>0, 1-based), so every vector the trust "
               "region forms already lives in the reduced space and no "
               "constrained degree of freedom is ever stepped.%s",
               iinc,iit,neq[0],neq[1],3*(*nk),lact,lspc,lmpcd,*nboun,*nmpc,
               nzs[1],nzs[2],"\n");
        printf("[DAMAGE TR LINCHECK] transpose identity dot(J^T R,p_N)/|R|^2 "
               "= %.12e (exact 1); operator asymmetry "
               "max|au_L-au_U|/max|au_L| = %.6e; |p_N|=%.6e |p_C|=%.6e "
               "|R|2=%.6e%s",
               dog.ident,dog.asym,sqrt(dog.npn2),
               dog.tc*sqrt(dog.nd2),sqrt(dog.nb2),"\n");
        fflush(stdout);

        isiz=mi[0]**ne;cpypardou(dog.dam,dam,&isiz,&num_cpus);
        if(damage_damvisc!=NULL){
          isiz=mi[0]**ne;
          cpypardou(dog.visc,damage_damvisc,&isiz,&num_cpus);
        }
        if((lnst>0)&&(dog.xs!=NULL)){
          isiz=lnst*mi[0]**ne;
          cpypardou(dog.xs,xstate,&isiz,&num_cpus);
        }
        for(ljj=0;ljj<4;ljj++) lqas[ljj]=qa[ljj];
        for(ljj=0;ljj<2;ljj++) luams[ljj]=uam[ljj];

        /* [WALLDIAG] The active set AT THE CURRENT ITERATE.  xstate, dam and
           stx here belong to u, because the iteration's own results() call
           built them and nothing has stepped yet.  Every rung below is
           compared against this one census, so a transition count is the
           number of integration points the step moved across a branch. */

        if(prb.wall_cat==NULL) NNEW(prb.wall_cat,ITG,mi[0]**ne);
        damage_ray_census(prb.wall_cat,xstate,xstateini,dam,damdamageini,
                          damage_damvisc,stx,ipkon,lakon,ne0,mi[0],*nstate_);
        {
          ITG wj,wnadv=0,wnlive=0;
          double wpinf=0.;
          for(wj=0;wj<mi[0]*ne0;wj++){
            if(prb.wall_cat[wj]&DAMCAT_USOFT) wnlive++;
            if(prb.wall_cat[wj]&DAMCAT_UADV) wnadv++;
          }
          for(wj=0;wj<neq[1];wj++)
            if(fabs(dog.pn[wj])>wpinf) wpinf=fabs(dog.pn[wj]);
          printf("[WALLDIAG] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                 " base state: UC6 points past initiation %" ITGFORMAT
                 ", of which ADVANCING (deff>dmax0) %" ITGFORMAT
                 "; |p_N|inf=%.6e |p_N|2=%.6e |R|2=%.6e%s",
                 iinc,iit,wnlive,wnadv,wpinf,sqrt(dog.npn2),
                 sqrt(dog.nb2),"\n");
          damage_wall_where("residual",dog.r0,neq[1],nactdof,mt,*nk,5,
                            ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,mi[0],
                            *nstate_);
          damage_wall_where("correction",dog.pn,neq[1],nactdof,mt,*nk,
                            5,ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,mi[0],
                            *nstate_);
          /* [WALLDIAG] STIFFNESS AT THE RESIDUAL PEAK.  The convergence test
             checkconvergence() applies is on max|R| over the mechanical
             block, not on |R|2, so the dof that decides the run is the peak
             one.  This reports, for the five largest residual dofs, the
             assembled diagonal AGAINST that node's own intact value and the
             displacement |R|/k that would be needed to null the residual
             locally.  A peak sitting on a node whose diagonal has collapsed,
             needing a displacement far larger than anything physical, is a
             different object from a peak on a healthy node, and only the
             second is a convergence problem in the ordinary sense. */
          {
            ITG *wsn=NULL,*wsd=NULL,wi,wj,wk,wt,wbest;
            double wa;
            NNEW(wsn,ITG,neq[1]);NNEW(wsd,ITG,neq[1]);
            for(wi=0;wi<neq[1];wi++){wsn[wi]=-1;wsd[wi]=0;}
            for(wi=0;wi<*nk;wi++)
              for(wj=1;wj<mt;wj++){
                wk=nactdof[mt*wi+wj];
                if((wk>0)&&(wk<=neq[1])){wsn[wk-1]=wi;wsd[wk-1]=wj;}
              }
            for(wt=0;wt<5;wt++){
              wbest=-1;wa=-1.;
              for(wi=0;wi<neq[1];wi++){
                if(wsn[wi]<0) continue;
                if(fabs(dog.r0[wi])>wa){wa=fabs(dog.r0[wi]);wbest=wi;}
              }
              if(wbest<0) break;
              wi=wsn[wbest];
              printf("[WALLDIAG]   Rpeak #%" ITGFORMAT ": node %" ITGFORMAT
                     " dir %" ITGFORMAT " R=%.6e  addiag=%.6e addiag0=%.6e "
                     "ratio=%.6e spc_masked=%" ITGFORMAT " need_du=|R|/k=%.6e"
                     "%s",wt+1,wi+1,wsd[wbest],dog.r0[wbest],
                     (damage_addiag!=NULL)?damage_addiag[wi]:0.,
                     (damage_addiag0!=NULL)?damage_addiag0[wi]:0.,
                     ((damage_addiag!=NULL)&&(damage_addiag0!=NULL)&&
                      (damage_addiag0[wi]>0.))?
                       damage_addiag[wi]/damage_addiag0[wi]:-1.,
                     ((damage_spc_mask!=NULL)&&(damage_spc_nk>wi))?
                       damage_spc_mask[wi]:-1,
                     ((damage_addiag!=NULL)&&(damage_addiag[wi]>0.))?
                       fabs(dog.r0[wbest])/damage_addiag[wi]:-1.,"\n");
              wsn[wbest]=-1;
            }
            SFREE(wsn);SFREE(wsd);
          }
          fflush(stdout);
        }

        {
        ITG lpassn=2;
        if((prb.wall_maskstep!=0)&&(dog.pm!=NULL)) lpassn=3;
        for(lpass=0;lpass<lpassn;lpass++){
          if(lpass==0){
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=dog.pn[ljj];
              ljp[ljj]=dog.r0[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 1: p = p_N, so J p = r0 "
                   "EXACTLY.  res(u+eps p)/((1-eps)|r0|) must go to 1 and "
                   "the defect to zero like O(eps).%s","\n");
          }else if(lpass==1){
            lsc=sqrt(dog.npn2/dog.nd2);
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=lsc*dog.d[ljj];
              ljp[ljj]=lsc*dog.w[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 2: p = %.6e * d (scaled to "
                   "|p_N|), so J p = %.6e * w.  d comes from the TRANSPOSE "
                   "loop and w from the forward loop, so a swapped pair "
                   "cannot pass this.%s",lsc,lsc,"\n");
          }else{
            /* [WALLDIAG] MASKSTEP.  p = p_N with the AUTOSPC-masked nodes'
               components zeroed.  UNSCALED on purpose: the question is not
               how this direction behaves at |p_N|, it is what the step
               actually does once the collapsed-diagonal dofs are taken out
               of it, so eps=1 here means "the rest of the Newton step". */
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=dog.pm[ljj];
              ljp[ljj]=dog.wm[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 3 (MASKSTEP): p = p_N with the "
                   "%" ITGFORMAT " AUTOSPC-masked nodes zeroed.  |p|=%.6e "
                   "against |p_N|=%.6e, so the mask removes %.4f%% of "
                   "|p_N|^2.  If the best rung here beats pass 1's best rung "
                   "(and the ladder's accepted alpha), the step is unusable "
                   "BECAUSE it is spent on collapsed-diagonal nodes; if it "
                   "does not, that hypothesis is dead.%s",
                   damage_spc_count,sqrt(dog.npm2),
                   sqrt(dog.npn2),
                   (dog.npn2>0.)?
                     100.*(1.-dog.npm2/dog.npn2):0.,"\n");
          }
          fflush(stdout);
          lnjp=0.;
          for(ljj=0;ljj<neq[1];ljj++) lnjp+=ljp[ljj]*ljp[ljj];
          lnjp=sqrt(lnjp);
          for(lii=0;lii<15;lii++){
            le=lcE[lii];
          /* [DAMAGE TR LINCHECK] cam starts clean for every probe, for the
             same reason as in the trust-region block. */
          for(ljj=0;ljj<3;ljj++) cam[ljj]=0.;
          for(ljj=3;ljj<5;ljj++) cam[ljj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,dog.dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,dog.visc,&isiz,&num_cpus);
          }
          if((lnst>0)&&(dog.xs!=NULL)){
            isiz=lnst*mi[0]**ne;
            cpypardou(xstate,dog.xs,&isiz,&num_cpus);
          }
          for(ljj=0;ljj<neq[1];ljj++) b[ljj]=le*lp[ljj];
          trial_residual(&nlgt,dog.res);

            dog.neval++;
            lnum=0.;lnres=0.;linf=0.;linf0=0.;
            if(prb.wall_def==NULL) NNEW(prb.wall_def,double,neq[1]);
            for(ljj=0;ljj<neq[1];ljj++){
              lmv=dog.r0[ljj]-le*ljp[ljj];
              ldv=dog.res[ljj]-lmv;
              prb.wall_def[ljj]=ldv;
              lnum+=ldv*ldv;
              lnres+=dog.res[ljj]*dog.res[ljj];
            }
            /* [WALLDIAG] The line search judges a trial by max|res| over the
               MECHANICAL block neq[0] (nonlingeo.c, damage_linesearch_*norm),
               while the direction it damps is only guaranteed to descend the
               2-norm merit - the transpose identity dot(J^T R,p_N)=|R|^2 is
               what makes that guarantee exact.  The two norms are therefore
               measured on the SAME rung here, or the disagreement between
               them stays an inference. */
            for(ljj=0;ljj<neq[0];ljj++){
              if(fabs(dog.res[ljj])>linf) linf=fabs(dog.res[ljj]);
              if(fabs(dog.r0[ljj])>linf0) linf0=fabs(dog.r0[ljj]);
            }
            lnum=sqrt(lnum);lnres=sqrt(lnres);
            lden=le*lnjp;
            printf("[DAMAGE TR LINCHECK] pass %" ITGFORMAT " eps=%.9f "
                   "|res(u+eps p)|=%.12e  linear-model defect=%.6e  "
                   "defect/eps=%.6e%s",lpass+1,le,lnres,
                   (lden>0.)?lnum/lden:0.,
                   (lden>0.)?lnum/(lden*le):0.,"\n");
            printf("[WALLDIAG] pass %" ITGFORMAT " eps=%.9f  TWO NORMS: "
                   "|res|2/|r0|2=%.12f  |res|inf/|r0|inf=%.12f  "
                   "(|r0|2=%.6e |r0|inf=%.6e)%s",lpass+1,le,
                   (dog.nb2>0.)?lnres/sqrt(dog.nb2):0.,
                   (linf0>0.)?linf/linf0:0.,sqrt(dog.nb2),linf0,"\n");
            if(lpass==0)
              printf("[DAMAGE TR LINCHECK]          ratio "
                     "|res|/((1-eps)|r0|) = %.12f  (must go to 1)%s",
                     ((1.-le)>0.)?lnres/((1.-le)*sqrt(dog.nb2)):0.,
                     "\n");
            {
              ITG wtot;
              wtot=damage_wall_setdiff(prb.wall_cat,xstate,xstateini,dam,
                                       damdamageini,damage_damvisc,stx,
                                       ipkon,lakon,ne0,mi[0],*nstate_,
                                       prb.wall_nb);
              printf("[WALLDIAG] pass %" ITGFORMAT " eps=%.9f  active-set "
                     "transitions %" ITGFORMAT " = UC6 loading/unloading %"
                     ITGFORMAT " + UC6 initiation %" ITGFORMAT
                     " + UC6 viscous %" ITGFORMAT " + UC6 failure %"
                     ITGFORMAT " + UC6 tension/compression %" ITGFORMAT
                     " + bulk plastic %" ITGFORMAT " + bulk initiation %"
                     ITGFORMAT " + bulk damage growth %" ITGFORMAT "%s",
                     lpass+1,le,wtot,prb.wall_nb[6],prb.wall_nb[3],
                     prb.wall_nb[4],prb.wall_nb[5],prb.wall_nb[7],
                     prb.wall_nb[0],prb.wall_nb[1],prb.wall_nb[2],
                     "\n");
            }
            /* [WALLDIAG] IS THE MISSING TERM THE DAMAGE RANK-1 TERM?

               J = K0 + E, where E is exactly what mafilldamas.f adds from
               damjac.  Assemble E ALONE into a zeroed pair and apply it to
               p_N.  If the true operator is A = K0 + (1+c)E - that is, if
               the rank-1 damage term is the right shape and the wrong size -
               then (J-A)p = -c*E*p, so the measured defect must be
               ANTI-PARALLEL to E*p and |defect|/(eps*|E p|) must equal |c|.
               A cosine of zero says the missing term is a different term,
               and no scaling of this one can supply it. */
            if(lii==14){
              ITG e1i,e1c,e1k,e1r;
              double *e1ad=NULL,*e1au=NULL,*e1y=NULL,e1n=0.,e1d=0.,e1w=0.;
              ITG e1nd,e1sk,e1ad2,e1ho,e1fl,e1lv,e1dg;
              ITG *e1cat=NULL;
              NNEW(e1ad,double,neq[1]);
              NNEW(e1au,double,(nasym+1)*nzs[1]);
              NNEW(e1y,double,neq[1]);
              e1nd=0;e1sk=0;e1ad2=0;e1ho=0;e1fl=0;e1lv=0;e1dg=0;
              NNEW(e1cat,ITG,ne0);
              FORTRAN(mafilldamas,(co,kon,ipkon,lakon,&ne0,nactdof,jq,irow,
                                   neq,nzs,e1au,e1ad,vold,mi,damage_damjac,
                                   nmpc,&e1nd,dam,damdamageini,
                                   &e1sk,&e1ad2,&e1ho,&e1fl,&e1lv,&e1dg,
                                   e1cat));
              SFREE(e1cat);
              for(e1k=0;e1k<neq[1];e1k++)
                e1y[e1k]=e1ad[e1k]*dog.pn[e1k];
              for(e1c=0;e1c<neq[1];e1c++){
                for(e1k=jq[e1c]-1;e1k<jq[e1c+1]-1;e1k++){
                  e1r=irow[e1k]-1;
                  e1y[e1r]+=e1au[e1k]*dog.pn[e1c];
                  e1y[e1c]+=e1au[nzs[2]+e1k]*dog.pn[e1r];
                }
              }
              for(e1i=0;e1i<neq[1];e1i++){
                e1n+=e1y[e1i]*e1y[e1i];
                e1d+=prb.wall_def[e1i]*prb.wall_def[e1i];
                e1w+=e1y[e1i]*prb.wall_def[e1i];
              }
              e1n=sqrt(e1n);e1d=sqrt(e1d);
              /* Is the rank-1 term small because dD/d(eps) is small, or
                 because the assembly loses it?  damjac slots 1..6 are the
                 effective stress and 7..12 are dD/d(eps); censusing both
                 separates "the derivative is tiny" from "the derivative is
                 there and the operator does not carry it". */
              {
                ITG qi,qn=0,qz=0;
                double qmax=0.,qsum=0.,tmax=0.;
                for(qi=0;qi<ne0;qi++){
                  ITG qk;double qa2=0.,qt2=0.;
                  if(ipkon[qi]<0) continue;
                  if(lakon[8*qi]!='C') continue;
                  for(qk=6;qk<12;qk++)
                    qa2+=damage_damjac[12*mi[0]*qi+qk]
                        *damage_damjac[12*mi[0]*qi+qk];
                  for(qk=0;qk<6;qk++)
                    qt2+=damage_damjac[12*mi[0]*qi+qk]
                        *damage_damjac[12*mi[0]*qi+qk];
                  if((qa2<=0.)&&(qt2<=0.)) continue;
                  qn++;
                  qa2=sqrt(qa2);qt2=sqrt(qt2);
                  if(qa2<=1.e-12) qz++;
                  if(qa2>qmax) qmax=qa2;
                  if(qt2>tmax) tmax=qt2;
                  qsum+=qa2;
                }
                printf("[WALLDIAG]   damjac census over %" ITGFORMAT
                       " elements with any entry: |dD/d(eps)| max=%.6e "
                       "mean=%.6e, of which %" ITGFORMAT
                       " have it BELOW 1e-12 (assembled anyway, because the "
                       "skip test sums the stress slots too); |sigma_eff| "
                       "max=%.6e%s",qn,qmax,(qn>0)?qsum/qn:0.,qz,tmax,"\n");
              }
              printf("[WALLDIAG]   damage rank-1 term applied to p_N: "
                     "elements %" ITGFORMAT ", |E p|=%.6e, |defect|/eps=%.6e,"
                     "  cos(defect,E p)=%+.6f  |defect|/(eps|E p|)=%.6f%s",
                     e1nd,e1n,e1d/le,
                     ((e1n>0.)&&(e1d>0.))?e1w/(e1n*e1d):0.,
                     (e1n>0.)?e1d/(le*e1n):0.,"\n");
              SFREE(e1ad);SFREE(e1au);SFREE(e1y);
            }
            if(lii==14){
              damage_wall_split("linear-model defect",prb.wall_def,neq[1],
                                nactdof,mt,*nk,ipkon,kon,lakon,dam,xstate,
                                *ne,ne0,mi[0],*nstate_);
              damage_wall_split("residual r0",dog.r0,neq[1],nactdof,mt,
                                *nk,ipkon,kon,lakon,dam,xstate,*ne,ne0,mi[0],
                                *nstate_);
              damage_wall_split("Newton step p_N",dog.pn,neq[1],nactdof,
                                mt,*nk,ipkon,kon,lakon,dam,xstate,*ne,ne0,
                                mi[0],*nstate_);
              damage_wall_where("defect",prb.wall_def,neq[1],nactdof,mt,
                                *nk,5,ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,
                                mi[0],*nstate_);
            }
            fflush(stdout);
          }
        }
        }

        /* leave the full Newton step, exactly as the unprobed code would */
        for(ljj=0;ljj<4;ljj++) qa[ljj]=lqas[ljj];
        for(ljj=0;ljj<2;ljj++) uam[ljj]=luams[ljj];
        le=1.;
        for(ljj=0;ljj<neq[1];ljj++) lp[ljj]=dog.pn[ljj];
          /* [DAMAGE TR LINCHECK] cam starts clean for every probe, for the
             same reason as in the trust-region block. */
          for(ljj=0;ljj<3;ljj++) cam[ljj]=0.;
          for(ljj=3;ljj<5;ljj++) cam[ljj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,dog.dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,dog.visc,&isiz,&num_cpus);
          }
          if((lnst>0)&&(dog.xs!=NULL)){
            isiz=lnst*mi[0]**ne;
            cpypardou(xstate,dog.xs,&isiz,&num_cpus);
          }
          for(ljj=0;ljj<neq[1];ljj++) b[ljj]=le*lp[ljj];
          trial_residual(&nlgt,dog.res);

        lnres=0.;
        for(ljj=0;ljj<neq[1];ljj++)
          lnres+=dog.res[ljj]*dog.res[ljj];
        printf("[DAMAGE TR LINCHECK] restored to the full Newton step; "
               "|res| there = %.12e.  Compare with the eps=1 line of pass 1: "
               "equal means every probe rolled the mutable state back "
               "completely.%s",sqrt(lnres),"\n");
        fflush(stdout);
        SFREE(lp);SFREE(ljp);
        dog.lc_due=0;
      }


      /* ============ [DAMAGE CT] bordered corrector ==================
         One factorisation per iteration, two solves against it:
             K z = -R = b        (b as the solver left it)
             K y = -q            (-q straight from two b vectors)
             den = c_lambda + c_u^T y
             dlambda = (-c - c_u^T z)/den ,   du = z + y*dlambda
         The factorisation has ALREADY been taken from the base state above,
         so the perturbation evaluations below cannot contaminate it: the
         cached LU is what pardiso_solve uses for y.
         c(u,lambda) = m.(delta - delta_c) - ds is exactly linear in u with
         a frozen m, so c is driven to zero every iteration by construction
         and its value is reported, never assumed. */

      ct.used=0;
      if((ct.on==1)&&(idamagereeq==0)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)&&
         (*mortar<=1)){
        damcont_corrector(&ct,&nlgt,xboun,xbounold,uam,damage_damjac,
                          damage_damvisc,inputformat,nrhs,symmetryflag);
      }

      /* ---- [DAMAGE TR] TRUST-REGION DOGLEG on the original residual -----

         Armed only as rescue level 3, i.e. only after the Rescue2 levels 1
         and 2 have both failed on the same wall, and only where the measured
         failure is a step-LENGTH failure.  At s3rad inc=569 BK3 reports, on
         all four identical attempts,

             res_old=1.942365e-03  res_full=2.520626e-02  lambda=0.100000
             trials=3  contracted=0  res_damped=2.001658e-03

         so the full Newton step overshoots 13x and the most damped step BK3
         is permitted to take still leaves the residual ABOVE where it
         started, because DAMAGE_LINESEARCH_MIN floors lambda at 0.1.  A
         trust region has no such floor, and it is not confined to the Newton
         direction.

             phi(u) = 1/2|R|^2,   R = f_int - f_ext = -r0
             p_N : J p_N = r0            (what PARDISO returned)
             d   = J^T r0 = -g           (steepest descent, formed as J^T)
             p_C = (|d|^2/|J d|^2) d
             p   = pa*d + pb*p_N         (Newton, Cauchy, or the blend)
             J p = pa*w + pb*r0          (so the model needs no matrix)
             rho = [phi(u)-phi(u+p)] / [phi(u) - 1/2|R+Jp|^2]

         Transactionality follows the accepted BT block: dam, damvisc and
         xstate are snapshotted once and restored before EVERY trial, and the
         loop always ends by re-evaluating the point it is going to keep, so
         what survives is that point own state and not a residue of the last
         trial.  cam, qa and uam are restored as well - cam[0] is a running
         MAXIMUM (resultsini.c:77-79) and uam[0]=max(uam[0],cam[0]) is never
         reset inside a step (nonlingeo.c) - so without this the rejected
         trials would permanently inflate the displacement criterion.  The
         kept point re-evaluation is compared against the value the trial
         pass measured for it: they must agree, and a mismatch is reported
         as an IMPURITY instead of being absorbed.

         Convergence of the INCREMENT is untouched: checkconvergence still
         decides on ram/cam/qa/qam from the unmodified residual. */

      dog.used=0;
      if((dog.on==1)&&(dog.have==1)&&(damage_de12_enabled)&&
         (ct.on==0)&&(ct.used==0)&&
         (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
         (*ithermal<2)&&(*idrct==0)&&(*mortar<=1)){
        /* what the region then does is dogleg.c's; the guard above -
           whether it may fire at all - is a decision about the increment
           and stays with the increment. */
        dogleg_rescue(&dog,&nlgt,&damage_glob,damage_damvisc,iit,icutb,uam);
      }

      rsc.linesearch_applied=0;
      if((rsc.linesearch_mode==1)&&(damage_de12_enabled)&&
         (dog.used==0)&&(dog.on==0)&&(ct.on==0)&&
         (ct.used==0)&&
         (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
         (*ithermal<2)&&(*idrct==0)){
        NNEW(res,double,neq[1]);
        trial_reduce(&nlgt,res);

        rsc.linesearch_oldnorm=0.;
        rsc.linesearch_fullnorm=0.;
        for(i=0;i<neq[0];i++){
          if(fabs(resold[i])>rsc.linesearch_oldnorm)
            rsc.linesearch_oldnorm=fabs(resold[i]);
          if(fabs(res[i])>rsc.linesearch_fullnorm)
            rsc.linesearch_fullnorm=fabs(res[i]);
        }

        rsc.linesearch_active=0;
        rsc.linesearch_nsoft=0;
        rsc.linesearch_maxdd=0.;
        if((isfinite(rsc.linesearch_oldnorm))&&
           (isfinite(rsc.linesearch_fullnorm))&&
           (rsc.linesearch_fullnorm>
              DAMAGE_LINESEARCH_GROWTH*rsc.linesearch_oldnorm)){
          rsc.linesearch_active=erosion_softening(
              dam,damdamageini,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,
              *ndmat_,*ntmat_,ne0,mi[0],&rsc.linesearch_nsoft,
              &rsc.linesearch_maxdd);
        }

        if(rsc.linesearch_active){
          sum1=0.;sum2=0.;
          for(i=0;i<neq[0];i++){
            sum1+=b[i]*resold[i];
            sum2+=b[i]*res[i];
          }

          if((!isfinite(sum1))||(!isfinite(sum2))||
             (fabs(sum1-sum2)<1.e-30)){
            flinesearch=DAMAGE_LINESEARCH_FALLBACK;
          }else{
            flinesearch=sum1/(sum1-sum2);
            if((!isfinite(flinesearch))||(flinesearch<=0.))
              flinesearch=DAMAGE_LINESEARCH_FALLBACK;
          }
          if(flinesearch>DAMAGE_LINESEARCH_MAX)
            flinesearch=DAMAGE_LINESEARCH_MAX;
          if(flinesearch<rsc.ls_min)
            flinesearch=rsc.ls_min;

          NNEW(rsc.linesearch_step,double,neq[1]);
          isiz=neq[1];cpypardou(rsc.linesearch_step,b,&isiz,&num_cpus);
          rsc.linesearch_contracted=0;
          glob_fired(&damage_glob,GLOB_LADDER);
          lsladder_start(&damage_lsl,flinesearch,
                         rsc.linesearch_oldnorm,rsc.ls_min,0.5,
                         rsc.ls_trials,rsc.ls_legacy);

          /* ---- BACKTRACKING LADDER PROBE (CCX_DAMAGE_LS_PROBE) --------
             The line search backtracks 1.0 -> 0.5 -> 0.1 and stops: three
             trials with a floor of 0.1.  At the increments that fail it
             reports contracted=0, i.e. even a tenth of the step makes the
             residual worse - by 14000x at increment 346.  Two things can
             cause that and they need different answers:

               the DIRECTION is usable and the floor is simply too coarse
               for a damage transition, in which case some smaller alpha
               contracts;

               the direction is wrong, in which case no alpha does.

             This walks a finer ladder and PRINTS the residual at each
             alpha.  It decides nothing: the normal trial loop below runs
             afterwards, unchanged, and its own final evaluation is what
             the solver keeps.  alpha=1 is evaluated TWICE so the log
             carries its own proof that a trial evaluation is
             reproducible - dam and damvisc are rebuilt from the committed
             baseline on every results() call, and this is what shows it
             rather than assuming it. */

          if(rsc.ls_probe!=0){
            static const double lsp[9]={1.,1.,0.5,0.2,0.1,0.03,0.01,
                                        0.003,0.001};
            ITG lpi;
            double lpr,lp0=-1.;
            for(lpi=0;lpi<9;lpi++){
              for(i=0;i<neq[1];i++)
                b[i]=lsp[lpi]*rsc.linesearch_step[i];
              trial_residual(&nlgt,res);
              lpr=0.;
              for(i=0;i<neq[0];i++) if(fabs(res[i])>lpr) lpr=fabs(res[i]);
              if(lpi==0) lp0=lpr;
              printf("[LS-PROBE] inc=%" ITGFORMAT " attempt=%" ITGFORMAT
                     " iter=%" ITGFORMAT " alpha=%.4f |R|inf=%.6e"
                     " (res_old=%.6e, ratio=%.3e)%s\n",
                     iinc,icutb+1,iit,lsp[lpi],lpr,
                     rsc.linesearch_oldnorm,
                     lpr/(rsc.linesearch_oldnorm+1.e-300),
                     ((lpi==1)&&(lp0>=0.))?
                     ((fabs(lpr-lp0)<=1.e-12*(lp0+1.e-300))?
                      "  [repeat of alpha=1: REPRODUCIBLE]":
                      "  [repeat of alpha=1: NOT REPRODUCIBLE]"):"");
            }
            fflush(stdout);
          }

          for(rsc.linesearch_trial=1;
              rsc.linesearch_trial<=rsc.ls_trials;
              rsc.linesearch_trial++){
            for(i=0;i<neq[1];i++)
              b[i]=flinesearch*rsc.linesearch_step[i];

            /* The bordered step is ONE Newton step in (u,lambda).
               Damping only its displacement half leaves the prescribed
               dofs at the undamped load factor, so the trial state is
               not on the ray the line search thinks it is contracting
               along - and the constraint it was solved to satisfy is not
               satisfied at any point of that ray.  Contract both. */

            if((pf.on==1)&&(pf.applied==1)){
              pf.lam=pf.lam0it+flinesearch*pf.dlamit;
              for(i=0;i<*nboun;i++)
                xbounact[i]=xbounold[i]+(xboun[i]-xbounold[i])*pf.lam;
            }

            trial_residual(&nlgt,res);
            rsc.linesearch_dampednorm=0.;
            for(i=0;i<neq[0];i++){
              if(fabs(res[i])>rsc.linesearch_dampednorm)
                rsc.linesearch_dampednorm=fabs(res[i]);
            }

            /* The ladder decides: accept, descend, or give up.  It is a
               separate unit (lsladder.c) with its own regression test,
               because both of the things it used to get wrong - the rungs
               it could reach and which rung it handed back - were the
               cause of the recorded wall. */

            damage_lsr=lsladder_step(&damage_lsl,
                                     rsc.linesearch_dampednorm);
            if(damage_lsr==1){
              rsc.linesearch_contracted=1;
              break;
            }
            if(damage_lsr<0) break;
            flinesearch=damage_lsl.alpha;
          }
          if(rsc.linesearch_trial>rsc.ls_trials)
            rsc.linesearch_trial=rsc.ls_trials;

          /* Nothing contracted.  Falling out with the LAST trial in b is
             the defect: the last trial is the floor, and the floor is not
             the best of what was measured.  Re-evaluate at the best alpha
             instead, so the step the solver takes is never worse than the
             best step this search actually saw. */

          if((rsc.linesearch_contracted==0)&&
             (fabs(lsladder_final(&damage_lsl)-flinesearch)>1.e-12)){
            flinesearch=lsladder_final(&damage_lsl);
            for(i=0;i<neq[1];i++)
              b[i]=flinesearch*rsc.linesearch_step[i];
            trial_residual(&nlgt,res);
            rsc.linesearch_dampednorm=damage_lsl.bestres;
            rsc.ls_bestused++;
          }
          SFREE(rsc.linesearch_step);
          rsc.linesearch_applied=1;
        }
        SFREE(res);
      }

      /* ---- [DAMAGE RAY] residual scan along the Newton direction --------
         Runs only where BK3 cannot: idamagereeq==1, the same-load solve after
         a topology change.  b still holds the Newton direction here; below,
         calcresidual overwrites b with the residual, so the direction is
         saved and restored explicitly.

         Two-point test first (alpha=0 and alpha=1) on every re-equilibration
         iteration - cheap, and alpha=0 is the reference this pass has no
         other way to obtain, resold being stale under this gate.  The full
         ray is walked only when the full step GREW the norm, which is the
         fatal signature, and at most prb.ray_max times.

         The walk ends on alpha=1, so v/stx/fn/dam/xstate are left in exactly
         the state the unprobed code would have produced. */
      prb.ray_incok=1;
      if(prb.ray_ninc>0){
        prb.ray_incok=0;
        for(i=0;i<prb.ray_ninc;i++)
          if(iinc==prb.ray_inc[i]) prb.ray_incok=1;
      }
      if((prb.ray_probe)&&(prb.ray_incok==1)&&
         (idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        /* 0 and 1 first; then the ladder from 1/1024 up, ending on the
           repeated 0.25 (purity pair, indices 10 and 12) and on alpha=1 so
           the state left behind is the unprobed one. */
        static const double rayA[]={0.,1.,
                                    0.0009765625,0.001953125,0.00390625,
                                    0.0078125,0.015625,0.03125,0.0625,0.125,
                                    0.25,0.5,0.25,1.};
        ITG rnpt,rii,rjj,rfull,rflip,rnpl[14],rnuc[14];
        double ra,rinf[14],rl2[14],rptr[14],rn2[14],rs;

        if(prb.ray_p==NULL){
          NNEW(prb.ray_p,double,neq[1]);
          NNEW(prb.ray_res,double,neq[1]);
        }
        isiz=neq[1];cpypardou(prb.ray_p,b,&isiz,&num_cpus);

        rs=(qam[0]>0.)?qam[0]:1.;
        rfull=0; rflip=0;
        rnpt=2;                       /* alpha=0 and alpha=1 only, at first */
        for(rii=0;rii<rnpt;rii++){
          ra=rayA[rii];
          for(rjj=0;rjj<neq[1];rjj++) b[rjj]=ra*prb.ray_p[rjj];
          trial_residual(&nlgt,prb.ray_res);
          rinf[rii]=0.; rl2[rii]=0.; rptr[rii]=0.;
          for(rjj=0;rjj<neq[0];rjj++){
            if(fabs(prb.ray_res[rjj])>rinf[rii])
              rinf[rii]=fabs(prb.ray_res[rjj]);
            rl2[rii]+=prb.ray_res[rjj]*prb.ray_res[rjj];
            rptr[rii]+=prb.ray_p[rjj]*prb.ray_res[rjj];
          }
          rn2[rii]=sqrt(rl2[rii]);
          rl2[rii]=rn2[rii]/rs;
          damage_ray_tally(xstate,xstateini,dam,damdamageini,damage_damvisc,
                           stx,ipkon,lakon,ne0,mi[0],*nstate_,
                           &rnpl[rii],&rnuc[rii]);

          /* alpha=0 is the reference: keep the WHOLE residual vector and the
             per-integration-point branch map.  Three scalars agreeing with the
             linear model prove nothing about the tangent VECTOR - that was an
             overreach in J-14 and it is corrected here by measuring the
             componentwise defect E(alpha)=R(alpha)-(1-alpha)R(0). */
          if(rii==0){
            if(prb.ray_r0==NULL) NNEW(prb.ray_r0,double,neq[1]);
            isiz=neq[1];cpypardou(prb.ray_r0,prb.ray_res,&isiz,&num_cpus);
            if(prb.ray_cat==NULL) NNEW(prb.ray_cat,ITG,mi[0]*ne0);
            damage_ray_census(prb.ray_cat,xstate,xstateini,dam,
                              damdamageini,damage_damvisc,stx,ipkon,lakon,
                              ne0,mi[0],*nstate_);
          }

          /* The fatal signature is decided on the RESIDUAL ALONE.  The
             budget must not enter this test: conflating "no overshoot" with
             "overshoot, but the walk budget is spent" makes the label lie,
             and it lied - the first build capped at 8 walks, all 8 were
             consumed by inc<=152, and inc=154 then printed flip=no for a
             step that had not been tested at all. */
          if((rii==1)&&(rinf[0]>0.)){
            if(rinf[1]>prb.ray_growth*rinf[0]) rflip=1;
            if((rflip)&&(prb.ray_shots<prb.ray_max)){
              rfull=1; rnpt=14; prb.ray_shots++;
            }
          }
        }

        printf("[DAMAGE RAY] inc=%" ITGFORMAT " iter=%" ITGFORMAT
               " time=%.12e dt=%.6e qam=%.6e"
               " R0=%.6e R1=%.6e ratio=%.3e flip=%s walked=%s%s",
               iinc,iit,theta**tper,dtime,qam[0],
               rinf[0],rinf[1],(rinf[0]>0.)?rinf[1]/rinf[0]:-1.,
               rflip?"YES":"no",rfull?"yes":"no",
               "\n");
        if(rfull){
          for(rii=0;rii<14;rii++){
            printf("   alpha=%.9f |R|inf=%.6e |R|2=%.6e phi=%.6e"
                   " |R|2/qam=%.6e pTR=%+.6e plastic_ip=%" ITGFORMAT
                   " ucomp_ip=%" ITGFORMAT "%s",
                   rayA[rii],rinf[rii],rn2[rii],0.5*rn2[rii]*rn2[rii],
                   rl2[rii],rptr[rii],rnpl[rii],rnuc[rii],
                   "\n");
          }
          /* Which measure would have accepted?  Armijo with c1=1e-4 on each
             of the three, against the alpha=0 reference. */
          for(rii=2;rii<12;rii++){
            printf("   ARMIJO alpha=%.9f  inf:%s  l2:%s  phi:%s%s",
                   rayA[rii],
                   (rinf[rii]<=(1.-1.e-4*rayA[rii])*rinf[0])?"PASS":"fail",
                   (rn2[rii] <=(1.-1.e-4*rayA[rii])*rn2[0]) ?"PASS":"fail",
                   (0.5*rn2[rii]*rn2[rii]<=
                    (1.-1.e-4*rayA[rii])*0.5*rn2[0]*rn2[0])?"PASS":"fail",
                   "\n");
          }
          /* componentwise defect against the linear model, and the branch
             census transitions.  This is the part that can actually tell a
             consistent tangent from a change of branch. */
          if(prb.ray_r0!=NULL){
            ITG dii,djj,dworst,dncat,dfirste,dfirstip,dfirsta,dfirstb;
            double de,dn2,di,r0n2,r0ni,dwv,dpred;
            r0n2=0.;r0ni=0.;
            for(djj=0;djj<neq[0];djj++){
              r0n2+=prb.ray_r0[djj]*prb.ray_r0[djj];
              if(fabs(prb.ray_r0[djj])>r0ni) r0ni=fabs(prb.ray_r0[djj]);
            }
            r0n2=sqrt(r0n2);
            if(r0n2<=0.) r0n2=1.e-300;
            if(r0ni<=0.) r0ni=1.e-300;
            for(dii=2;dii<14;dii++){
              for(djj=0;djj<neq[1];djj++)
                b[djj]=rayA[dii]*prb.ray_p[djj];
              trial_residual(&nlgt,prb.ray_res);
              dn2=0.;di=0.;dworst=-1;dwv=0.;dpred=0.;
              for(djj=0;djj<neq[0];djj++){
                de=prb.ray_res[djj]-(1.-rayA[dii])*prb.ray_r0[djj];
                dn2+=de*de;
                if(fabs(de)>di){
                  di=fabs(de);dworst=djj;dwv=prb.ray_res[djj];
                  dpred=(1.-rayA[dii])*prb.ray_r0[djj];
                }
              }
              dn2=sqrt(dn2);
              dncat=damage_ray_census_diff(prb.ray_cat,xstate,xstateini,
                        dam,damdamageini,damage_damvisc,stx,ipkon,lakon,ne0,
                        mi[0],*nstate_,&dfirste,&dfirstip,&dfirsta,&dfirstb);
              printf("   DEFECT a=%.6f |E|2/|R0|2=%.6e |E|inf/|R0|inf=%.6e"
                     " worstdof=%" ITGFORMAT " R=%.6e lin=%.6e E=%.6e"
                     "  switched_ip=%" ITGFORMAT,
                     rayA[dii],dn2/r0n2,di/r0ni,dworst,dwv,dpred,dwv-dpred,
                     dncat);
              if(dncat>0)
                printf("  first: el=%" ITGFORMAT " ip=%" ITGFORMAT
                       " cat %" ITGFORMAT "->%" ITGFORMAT,
                       dfirste,dfirstip,dfirsta,dfirstb);
              printf("%s","\n");
            }
            fflush(stdout);
          }
          /* purity: index 10 and index 12 are the same alpha, twice */
          printf("   PURITY alpha=%.6f evaluated twice: |R|inf %.17e vs %.17e"
                 "  %s%s",rayA[10],rinf[10],rinf[12],
                 (rinf[10]==rinf[12])?"BITWISE IDENTICAL":"*** DIFFERS ***",
                 "\n");
        }
        fflush(stdout);

        /* restore the Newton direction; v/stx/fn already hold the alpha=1
           state because the ray ends there */
        isiz=neq[1];cpypardou(b,prb.ray_p,&isiz,&num_cpus);
      }

      /* ---- [DAMAGE BT] transactional backtracking in re-equilibration ----
         THIS IS A SOLVER CHANGE, not a diagnostic.  It decides which point
         the Newton iteration continues from.

         BK3 cannot serve here: it is locked out of idamagereeq by two
         independent gates, its floor is 0.10, and on failure it silently
         keeps the floor-length step whose residual is known to be worse.
         This block restores exactly and hands the increment back to the
         standard cutback when no step is acceptable.

         Contract, one probe at a time:
           - topology is whatever the current remastruct produced; no deletion
             scan and no active-set change happens inside the search;
           - dam / damvisc / xstate are restored from the snapshot taken
             before the search, so every probe starts from the same state;
           - a rejected probe leaves nothing behind: the loop always ends by
             re-evaluating the point it is going to keep, so the surviving
             state is that point's own and not a residue of the last trial;
           - if nothing is acceptable the FULL step is restored and the normal
             divergence / cutback machinery takes over unchanged.

         Acceptance is Armijo on the mechanical infinity norm:
             |R(alpha)|inf <= (1 - c1*alpha) * |R(0)|inf ,  c1 = 1e-4.
         Default off; unset the flag and not one array is allocated. */
      if(((rsc.bt_mode==1)||(rsc.rescue_bt_on==1))&&
         (idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        /* what the ladder then does is rescue.c's; the guard above -
           whether it may run at all - is a decision about the increment. */
        rescue_backtrack(&rsc,&nlgt,&damage_glob,&prb,damage_damvisc,iit);
      }

      /* ---- [DAMAGE ABA] full-state purity test, fires once --------------- */
      prb.aba_hit=0;
      if((prb.aba_mode==1)&&(idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        if(prb.aba_ninc==0){
          if(prb.aba_done==0) prb.aba_hit=1;
        }else{
          for(i=0;i<prb.aba_ninc;i++)
            if((prb.aba_inc[i]>0)&&(iinc==prb.aba_inc[i])){
              prb.aba_hit=1; prb.aba_inc[i]=-prb.aba_inc[i];
            }
        }
      }
      if(prb.aba_hit==1){
        ITG apass,abad=0,anst,nstx,nxst,neme;
        double aal;
        anst=*nstate_;
        nstx=6*mi[0]**ne;
        nxst=27*mi[0]**ne;
        neme=6*mi[0]**ne;
        if(prb.ray_p==NULL){
          NNEW(prb.ray_p,double,neq[1]);
          NNEW(prb.ray_res,double,neq[1]);
        }
        isiz=neq[1];cpypardou(prb.ray_p,b,&isiz,&num_cpus);
        NNEW(daba_res,double,neq[1]);
        NNEW(daba_v,double,mt**nk);
        NNEW(daba_stx,double,nstx);
        NNEW(daba_fn,double,mt**nk);
        NNEW(daba_f,double,neq[1]);
        NNEW(daba_dam,double,mi[0]**ne);
        NNEW(daba_visc,double,mi[0]**ne);
        if(anst>0) NNEW(daba_xs,double,anst*mi[0]**ne);
        NNEW(daba_eme,double,neme);
        NNEW(daba_stiff,double,nxst);
        NNEW(daba_qa,double,4);
        NNEW(daba_cam,double,5);

        for(apass=0;apass<3;apass++){
          aal=(apass==1)?(1.-prb.aba_a):prb.aba_a;
          for(i=0;i<neq[1];i++) b[i]=aal*prb.ray_p[i];
          trial_residual(&nlgt,prb.ray_res);
          if(apass==0){
            isiz=neq[1];cpypardou(daba_res,prb.ray_res,&isiz,&num_cpus);
            isiz=mt**nk;cpypardou(daba_v,v,&isiz,&num_cpus);
            isiz=nstx;cpypardou(daba_stx,stx,&isiz,&num_cpus);
            isiz=mt**nk;cpypardou(daba_fn,fn,&isiz,&num_cpus);
            isiz=neq[1];cpypardou(daba_f,f,&isiz,&num_cpus);
            isiz=mi[0]**ne;cpypardou(daba_dam,dam,&isiz,&num_cpus);
            if(damage_damvisc!=NULL){
              isiz=mi[0]**ne;
              cpypardou(daba_visc,damage_damvisc,&isiz,&num_cpus);
            }
            if(anst>0){
              isiz=anst*mi[0]**ne;cpypardou(daba_xs,xstate,&isiz,&num_cpus);
            }
            isiz=neme;cpypardou(daba_eme,eme,&isiz,&num_cpus);
            isiz=nxst;cpypardou(daba_stiff,xstiff,&isiz,&num_cpus);
            for(i=0;i<4;i++) daba_qa[i]=qa[i];
            for(i=0;i<5;i++) daba_cam[i]=cam[i];
          }else if(apass==2){
            printf("[DAMAGE ABA] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " A=%.6f B=%.6f  full-state comparison of the two A"
                   " evaluations:%s",iinc,iit,prb.aba_a,
                   1.-prb.aba_a,"\n");
            damage_aba_cmp("residual",daba_res,prb.ray_res,neq[1],&abad);
            damage_aba_cmp("v",daba_v,v,mt**nk,&abad);
            damage_aba_cmp("stx",daba_stx,stx,nstx,&abad);
            damage_aba_cmp("fn",daba_fn,fn,mt**nk,&abad);
            damage_aba_cmp("f",daba_f,f,neq[1],&abad);
            damage_aba_cmp("dam",daba_dam,dam,mi[0]**ne,&abad);
            damage_aba_cmp("damvisc",daba_visc,damage_damvisc,
                           mi[0]**ne,&abad);
            damage_aba_cmp("xstate",daba_xs,xstate,anst*mi[0]**ne,&abad);
            damage_aba_cmp("eme",daba_eme,eme,neme,&abad);
            damage_aba_cmp("xstiff",daba_stiff,xstiff,nxst,&abad);
            damage_aba_cmp("qa",daba_qa,qa,4,&abad);
            damage_aba_cmp("cam",daba_cam,cam,5,&abad);
            printf("[DAMAGE ABA] VERDICT: %" ITGFORMAT
                   " array(s) differ -> %s%s",abad,
                   (abad==0)?"TRANSACTION IS CLEAN":
                             "STATE LEAK - a line search on this path compares"
                             " different physical states","\n");
            fflush(stdout);
          }
        }
        SFREE(daba_res);SFREE(daba_v);SFREE(daba_stx);SFREE(daba_fn);
        SFREE(daba_f);SFREE(daba_dam);SFREE(daba_visc);
        if(anst>0) SFREE(daba_xs);
        SFREE(daba_eme);SFREE(daba_stiff);SFREE(daba_qa);SFREE(daba_cam);
        prb.aba_done=1;

        /* leave the full step, exactly as the unprobed code would have */
        isiz=neq[1];cpypardou(b,prb.ray_p,&isiz,&num_cpus);
        trial_evaluate(&nlgt);
      }

      /* calculating the residual */

      // next line: change on 19072022
      if((*iexpl<=1)||(*nener==1)){
	trial_reduce(&nlgt,b);
      }

      if(rsc.linesearch_applied){
        rsc.linesearch_dampednorm=0.;
        for(i=0;i<neq[0];i++){
          if(fabs(b[i])>rsc.linesearch_dampednorm)
            rsc.linesearch_dampednorm=fabs(b[i]);
        }
        printf("[DAMAGE LINESEARCH BK3] inc=%" ITGFORMAT
               " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
               " soft_elements=%" ITGFORMAT " max_dD=%.6e "
               "res_old=%.6e res_full=%.6e lambda=%.6f "
               "trials=%" ITGFORMAT " contracted=%" ITGFORMAT
               " res_damped=%.6e\n",iinc,icutb+1,iit,
               rsc.linesearch_nsoft,rsc.linesearch_maxdd,
               rsc.linesearch_oldnorm,rsc.linesearch_fullnorm,
               flinesearch,rsc.linesearch_trial,
               rsc.linesearch_contracted,rsc.linesearch_dampednorm);
        fflush(stdout);
      }

      /* fix residuals for mortar contact, add contact forces */
      
      if(*mortar>1){
	for(k=0;k<neq[1];k++){
	  b[k]=b[k]-f_cs[k]-f_cm[k];}
      }	 

      isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
      if(*ithermal!=2){
	// next line: change on 19072022
	//if((*ithermal!=2)&&((*iexpl<=1)||(*nener==1))){
	for(k=0;k<6*mi[0]*ne0;++k){
	  sti[k]=stx[k];
	}
      
	/* calculating the ratio of the smallest to largest pressure
	   for face-to-face contact
	   only done at the end of a step */

	if((*mortar==1)&&(1.-theta-dtheta<=1.e-6)){
	  FORTRAN(negativepressure,(&ne0,ne,mi,stx,&pressureratio));
	}else{pressureratio=0.;}
      }

      SFREE(v);SFREE(stx);SFREE(fn);
      
      if((idamping==1)&&(*iexpl<=1)){SFREE(adc);SFREE(auc);}

      if(*iexpl<=1){
	  
	/* store the residual forces for the next iteration */

	/* [CONVERGE] the numbers the judgement is made from now have an
	   owner: converge.c.  The
	   AUTOSPC-FORCE exclusion used to be a `continue` in the middle of
	   the loop that computes ram[0]; it is now a named property of the
	   object that owns the reduction, and what it excluded is read back
	   below rather than left in three file-scope variables.  Pure code
	   movement: bit-identical, checked on the gate and on the deck. */

	damage_cvg.qam_floor=damage_qam_floor;
	damage_cvg.mask_force=damage_spc_force;
	damage_cvg.mask=damage_spc_mask;
	damage_cvg.mask_nk=damage_spc_nk;
	glob_iterate(&damage_glob,b,neq[0]);
	converge_norms(&damage_cvg,b,neq,nactdofinv,mt,*ithermal,*mortar,
	               *ne,ne0,neold,qa,qamold,jnz,qau,ea,
	               ram,ram1,ram2,cam,uam,qam);
	/* [CONVERGE] and presenting them is the second call.  nonlingeo()
	   now says "compute the norms" and "show them" where it had ninety
	   lines of arithmetic interleaved with eleven printf.  The three
	   file-scope variables that carried the excluded peak between the
	   two are gone: the object that excludes is the object that
	   reports. */

	converge_report(&damage_cvg,nactdofinv,mt,*ithermal,ctrl[18],
	                qa,qam,ram,cam,uam);

	FORTRAN(writecvg,(istep,&iinc,&icutb,&iit,ne,&ne0,ram,qam,cam,uam,
			  ithermal));

        /* NC1 damage-aware slow-convergence gate.

           Keep checkconvergence() itself stock.  ctrl[3] (ic) is raised only
           for this call when progressive damage is actually changing, or
           during the same-load equilibrium solve of a DE1.3 terminal
           deletion, and both the residual and correction have contracted
           twice.  If a
           previously extended solve loses that trend after stock ic, setting
           ic to the current iteration makes the existing "iit==ic" branch
           perform the normal dc cutback instead of risking an unbounded run.
           ctrl[3] is restored immediately after the call. */

        ctrl[3]=icref;
        slow.active=0;
        slow.allow=0;
        slow.nsoft=0;
        slow.maxdd=0.;
        slow.estres=DAMAGE_SLOW_NEWTON_INVALID_EST;
        slow.estcorr=DAMAGE_SLOW_NEWTON_INVALID_EST;
        slow.esttotal=DAMAGE_SLOW_NEWTON_INVALID_EST;
        slow.rratio=0.;
        slow.cratio=0.;

        if((damage_de12_enabled)&&(damdamageini!=NULL)&&
           (*iexpl<=1)&&(*nmethod!=4)&&(*nmethod!=5)&&
           (*ithermal<2)&&(*idrct==0)&&(ncont==0)&&
           (iit>=(ITG)irref)){
          /* During a terminal-deletion transaction the surviving elements
             need not acquire additional D: the topology change itself is the
             damage event being equilibrated.  Treat that same-load solve as
             damage-active, while retaining all monotonicity and cost gates
             below.  Legacy A3 transactions do not set
             damage_de13_transaction and therefore remain stock. */
          if((idamagereeq==1)&&(damage_de13_transaction==1)){
            slow.active=1;
            slow.nsoft=dtxn.count;
          }else{
            slow.active=erosion_softening(
                dam,damdamageini,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,
                *ndmat_,*ntmat_,ne0,mi[0],&slow.nsoft,
                &slow.maxdd);
          }
        }

        if(slow.active){
          slow.allow=slownewton_allow(
              iit,ram,ram1,ram2,cam,uam,slow.camprev1,
              slow.camprev2,qa,qam,ctrl,slow.maxiters,
              &slow.estres,&slow.estcorr,
              &slow.esttotal,&slow.rratio,
              &slow.cratio);
        }

        if((slow.allow)&&
           ((slow.esttotal>(ITG)icref)||slow.extended)){
          ctrl[3]=(double)slow.maxiters+0.5;
          if((slow.extended==0)||(iit==(ITG)icref)||
             ((iit%5)==0)){
            printf("[DAMAGE NEWTON EXTEND] inc=%" ITGFORMAT
                   " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
                   " soft_elements=%" ITGFORMAT " max_dD=%.6e "
                   "est_res=%" ITGFORMAT " est_corr=%" ITGFORMAT
                   " est_total=%" ITGFORMAT " cap=%" ITGFORMAT
                   " r_ratio=%.6f c_ratio=%.6f\n",
                   iinc,icutb+1,iit,slow.nsoft,slow.maxdd,
                   slow.estres,slow.estcorr,
                   slow.esttotal,slow.maxiters,
                   slow.rratio,slow.cratio);
            fflush(stdout);
          }
          slow.extended=1;
        }else if((slow.extended)&&(iit>(ITG)icref)){
          /* Avoid checkconvergence's fatal iit>ic guard: ic==iit reaches its
             ordinary too-slow cutback path on this iteration. */
          ctrl[3]=(double)iit+0.5;
          printf("[DAMAGE NEWTON STOP] inc=%" ITGFORMAT
                 " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
                 " reason=lost-contraction-or-cost-bound; stock cutback\n",
                 iinc,icutb+1,iit);
          fflush(stdout);
        }

        slow.camprev2=slow.camprev1;
        slow.camprev1=cam[0];

        /* NC2: a terminal topology solve has zero external load increment.
           Stock therefore normalizes its correction by the tiny displacement
           caused by deletion alone.  This can demand many factorizations even
           after the force residual satisfies the stock tolerance.  For this
           nested solve only, retain the displacement scale of the already
           converged physical increment.  No tolerance is relaxed: the same
           stock relative correction and force-residual tests are applied to
           the combined physical-increment/topology correction. */
        /* checkconvergence advances theta by dtheta on acceptance and
           sizes the next increment from dthetaref, which is untouched,
           so overriding dtheta here lands theta exactly on the load
           factor the constraint produced without disturbing the stock
           step controller. */

        /* While a descent is walking, lambda is no longer tied to
           theta, so the dissipation controller's handover would
           overwrite the step the descent just restored - which is
           why the descent only ever got three attempts. */
        if((lc.diss_ctrl>=1)&&(lc.path_desc==0)&&(lc.arc==0)){
          dtheta=lc.diss_lamcur-theta;
          if(dtheta<*tmin) dtheta=*tmin;
        }

        rsc.reeq_uam_actual[0]=uam[0];
        rsc.reeq_uam_actual[1]=uam[1];
        if((rsc.reeq_scale_mode==1)&&(idamagereeq==1)&&
           (damage_de13_transaction==1)){
          if(uam[0]<rsc.reeq_uam_ref[0])
            uam[0]=rsc.reeq_uam_ref[0];
          if(uam[1]<rsc.reeq_uam_ref[1])
            uam[1]=rsc.reeq_uam_ref[1];
          if(iit==1){
            printf("[DAMAGE REEQ NC2] inc=%" ITGFORMAT
                   " actual_uam=%.6e physical_ref=%.6e\n",
                   iinc,rsc.reeq_uam_actual[0],rsc.reeq_uam_ref[0]);
            fflush(stdout);
          }
        }

        /* [DAMAGE RESCUE] the step the failing attempt actually used is
           the last admissible one; checkconvergence is about to shrink it. */
        rsc.rescue_dtheta_last=dtheta;
        rsc.rescue_dthetaref_last=dthetaref;
        ccx_rescue_req=0;
        ccx_rescue_arm=0;
        if((rsc.rescue_mode==1)&&(rsc.rec_disarmed==0)&&
           (rsc.rescue_used<rsc.rescue_maxlevel)&&(ncont==0)&&
           (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
          ccx_rescue_arm=1;
          /* [DAMAGE CT] level 4 is granted only if the continuation could
             actually arm.  The candidate scan needs committed ring data
             only, so it is evaluated HERE, before checkconvergence defers
             the stop - a refusal then leaves the stock path byte-identical
             to the Rescue2+dogleg control. */
          if((ct.mode==1)&&(rsc.rescue_used>=3)&&
             (ct.on==0)){
            ITG cse,csi,csn,csr,csok=0;
            double csm[3],csk,csd,cst;
            if((ct.alloc==1)&&(ct.nring>=6)&&
               (ct.refused==0)){
              csok=damcont_select(ct.ring,ct.fl,ct.dt,
                                    ipkon,lakon,ielprop,prop,ne0,mi[0],
                                    ct.head,&cse,&csi,csm,&csk,&csd,
                                    &cst,&csn,&csr,ct.kaptol,ct.bl,ct.nbl);
            }else{csn=0;csr=1;}
            if(csok==0){
              printf("[DAMAGE CT] level 4 NOT granted at inc=%" ITGFORMAT
                     " (%s; candidates=%" ITGFORMAT ", ring=%" ITGFORMAT
                     "/6).  No continuation attempt is taken, the stock path "
                     "is untouched and this wall goes to the ORIGINAL stock "
                     "stop.%s",iinc,
                     (ct.alloc==0)?"ring not allocated":
                     ((ct.nring<6)?"ring incomplete":
                      ((ct.refused!=0)?"already refused once":
                       ((csr==2)?"kappa unstable over the five intervals":
                        "no candidate with five admissible intervals"))),
                     csn,ct.nring,"\n");
              fflush(stdout);
              ct.refused=1;
              ccx_rescue_arm=0;
            }
          }
        }

	lc.arc_theta0=theta;
	checkconvergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
			 kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
			 xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
			 ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			 ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
			 vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
			 &theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
			 nactdof,b,tmin,ctrl,amta,namta,itpamp,inext,&dthetaref,
			 &itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
			 nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
			 set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
			 mortar,nmat,ielprop,prop,&ialeatoric,&kscale,
			 energy,&allwk,&energyref,&emax,&r_abs,&enetoll,
			 energyini,
			 &allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			 &dampwkini,energystartstep);

        /* [DAMAGE RESCUE] checkconvergence deferred the stop.  It left
           exactly as an ordinary cutback leaves, so the standard rollback at
           `if(icutb!=0)` further down restores the start of this increment
           from the baselines saved when it began - no new snapshot.  Undo
           only its step reduction and arm BT for this one attempt. */
        rescue_attempt(&rsc,&dog,&ct,&lc,&prb,&damage_glob,&nlgt,
                       &dtheta,&dthetaref,theta,tper,iit,idamagereeq);

        /* checkconvergence advances theta itself, so an ACCEPTED increment is
           exactly one where theta moved.  Commit lambda there and nowhere
           else: the first version committed it in the dissipation-report
           block, gated on icutb==0 && idamagereeq==0, so every cutback and
           every re-equilibration advanced theta while leaving lambda behind.
           The drift accumulated - lambda 0.2733 against theta 0.2706 - and the
           run stalled at 98.6% of peak with 38 deletions. */
        if((lc.arc==1)&&(theta>lc.arc_theta0)){
          /* Until the constraint engages, lambda has no equation of its own
             and must track theta EXACTLY.  It cannot be taken from
             lc.diss_lamcur: that was built at the top of the increment
             from the TRIAL dtheta, while checkconvergence advances theta by
             whatever dtheta it settles on.  The two drift, the dG formula
             below then pairs a current theta with a stale lambda, dG comes
             out wrong and the engagement threshold never fires - measured,
             the arc arm never engaged while the identical run without it
             engaged and reached 70.9% of peak (E-96). */
          if(lc.diss_engaged==1){
            lc.arc_lam=lc.diss_lamcur;
          }else{
            lc.arc_lam+=theta-lc.arc_theta0;
          }
          lc.diss_lamcur=lc.arc_lam;
        }

        ctrl[3]=icref;
        uam[0]=rsc.reeq_uam_actual[0];
        uam[1]=rsc.reeq_uam_actual[1];

	if(*mortar>1){
	  SFREE(f_cs);SFREE(f_cm);
	} 
	  
      }else{

	/* explicit dynamics */

	icntrl=1;
	icutb=0;   

	theta=theta+dtheta;  
	if(dtheta>=1.-theta){
	  if(dtheta>1.-theta){
	    printf(" the increment size exceeds the remainder of the step and is decreased to %e\n\n",
		   dtheta**tper);
	  }
	  dtheta=1.-theta;
	  dthetaref=dtheta;
	}
	iflagact=0;
      }
    }

    if((*mortar==-1)&&(masslesslinear==0)&&(ncont!=0))
      {SFREE(auw);SFREE(jqw);SFREE(iroww);}

    if(*nmethod!=4)SFREE(resold);

    /*********************************************************/
    /*   end of the iteration loop                          */
    /*********************************************************/

    /* icutb=0 means that the iterations in the increment converged,
       icutb!=0 indicates that the increment has to be reiterated with
       another increment size (dtheta) */

    if(*mortar>1){
      SFREE(aubd);SFREE(jqbd);SFREE(irowbd);
      SFREE(aubdtil);SFREE(jqbdtil);SFREE(irowbdtil);
      SFREE(aubdtil2);SFREE(jqbdtil2);SFREE(irowbdtil2);
      SFREE(audd);SFREE(jqdd);SFREE(irowdd);
      SFREE(auddinv);SFREE(jqddinv);SFREE(irowddinv);
      SFREE(auddtil);SFREE(jqddtil);SFREE(irowddtil);
      SFREE(auddtil2);SFREE(jqddtil2);SFREE(irowddtil2);
      SFREE(bhat);	  
      SFREE(islavactdof);
    }

    if((icutb==0)&&(idamagereeq==1)){

      /* The same-load Newton solve on the current eroded topology has
         converged.  Recompute damage for every surviving element from the
         immutable physical-increment baseline.  This closes the coupling

              equilibrium -> damage -> topology -> equilibrium

         at one load level and, because dambase=damdamageini, avoids
         double-counting xstate-xstateini on repeated active-set passes. */

      damage_mode=1;
      idamage=0;
      damage_alphaevent=2.;
      damagebase=damdamageini;

      isiz=mi[0]*ne0;
      cpypardou(damde1prev,dam,&isiz,&num_cpus);

      FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
                          ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
                          ndmcon,dmcon,dam,damagebase,&dtime,sti,
                          ithermal,t1,xstate,xstateini,nstate_,vold,
                          &idamage,&damage_mode,&damage_alphaevent));

      damstats_element(dam,damde1prev,ipkon,lakon,ne0,mi[0],
                       &de1.nactive,&de1.gt01,
                       &de1.gt05,&de1.gt09,
                       &de1.nfull,&de1.nchanged,
                       &de1.dmax,&de1.maxdelta);

      /* DE1.3 terminal active-set extension.  The just-converged same-load
         equilibrium may have driven additional surviving DE1.2 elements
         to terminal degradation.  Mark only a bounded batch; the next
         same-load solve will expose any further terminal candidates. */
      damage_ebatch.marked=0;
      damage_ebatch.batch_dmax=0.;
      if(damage_de12_enabled){
        erosion_mark(&damage_epol,&damage_ebatch,
                     dam,damage_damvisc,ipkon,lakon,kon,ielmat,matname,
                     ndmcon,dmcon,*ndmat_,*ntmat_,*nk,*ne,ne0,mi[0],mi[2],
                     damage_de13_trigger_value,damage_de13_trigger_ip,
                     iinc,theta**tper);
        if(damage_ebatch.marked>0){
          damage_de13_transaction=1;
          idamage+=damage_ebatch.marked;
          printf("[DAMAGE DE1.3 EXTEND] inc=%" ITGFORMAT
                 " pass=%" ITGFORMAT " time=%.12e new_terminal=%"
                 ITGFORMAT " batch_Dmax=%.6e batch_Dvis=%.6e\n",
                 iinc,damage_active_pass+1,theta**tper,damage_ebatch.marked,
                 damage_ebatch.batch_dmax,damage_ebatch.batch_vmin);
          fflush(stdout);
        }
      }

      if(idamage>0){

        /* Redistribution on the current topology has caused additional
           elements to reach the damage limit at the SAME load level.
           These are not committed yet.  Extend the tentative topology,
           rebuild the equation structure and equilibrate once more. */

        damage_active_pass++;

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            stiffcensus damage_census;
            stiffcensus_take(&damage_census,*nk,damage_addiag,
                             damage_addiag0,damage_addok);
            if(damage_census_ok) monitor_stiffness(&damage_census,iinc,
                                                   theta**tper);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            /* H4 probe (E-24).  For every node that has lost all its bulk,
               assemble the local support operator from its LIVE cohesive
               facets

                 K = sum_ip (A/3) * shape_i(ip)^2 * g_ip * Kn *
                     ( n (x) n + beta*(t1 (x) t1 + t2 (x) t2) )

               and report lambda_min(K)/lambda_min(K_intact), where K_intact
               is the same sum with g=1.

               This is area- and direction-weighted and additive, which the
               previous version of this probe was not.  Max over a node's
               facets is the WRONG statistic (E-15 correction 2: one
               surviving integration point at 0.698 masked seventeen dead
               ones), and a plain sum cannot see rank collapse (H3): with one
               point bonded out of eighteen the operator has one stiff
               direction and two at gmin, which only an eigenvalue exposes.

               Reported as a distribution, because no threshold may be
               proposed until PASS and FAIL separate in one. */
            ITG h4i,h4j,h4k,h4m,h4n,h4np,h4a,h4b,h4p;
            double h4e1[3],h4e2[3],h4cr[3],h4nr[3],h4t1[3],h4t2[3],h4T[9];
            double h4nc,h4n1,h4area,h4kn,h4tn,h4ts,h4gm,h4bta,h4g,h4w,h4sh;
            double h4A[9],h4V[9],h4ev[3],h4r,h4off,h4th,h4c,h4s,h4t,h4tau;
            double *h4loc=NULL,*h4int=NULL,*h4rat=NULL;
            ITG h4cnt,h4tot;

            NNEW(damage_free_nb,ITG,*nk);
            NNEW(h4loc,double,9**nk);
            NNEW(h4int,double,9**nk);
            NNEW(h4rat,double,*nk);
            for(h4i=0;h4i<*nk;h4i++){damage_free_nb[h4i]=0;h4rat[h4i]=-1.;}
            for(h4i=0;h4i<9**nk;h4i++){h4loc[h4i]=0.;h4int[h4i]=0.;}

            /* live bulk per node */
            for(h4i=0;h4i<*ne;h4i++){
              if(ipkon[h4i]<0) continue;
              if(lakon[8*h4i]!='C') continue;
              h4np=4;
              if(lakon[8*h4i+3]=='6') h4np=6;
              else if(lakon[8*h4i+3]=='8') h4np=8;
              else if(lakon[8*h4i+3]=='1') h4np=10;
              for(h4j=0;h4j<h4np;h4j++){
                h4k=kon[ipkon[h4i]+h4j]-1;
                if((h4k>=0)&&(h4k<*nk)) damage_free_nb[h4k]++;
              }
            }

            /* local support operator from live cohesive facets */
            for(h4i=0;h4i<*ne;h4i++){
              if(ipkon[h4i]<0) continue;
              if(lakon[8*h4i]!='U') continue;
              h4a=kon[ipkon[h4i]]-1;
              h4b=kon[ipkon[h4i]+1]-1;
              h4n=kon[ipkon[h4i]+2]-1;
              if((h4a<0)||(h4b<0)||(h4n<0)) continue;
              for(h4j=0;h4j<3;h4j++){
                h4e1[h4j]=co[3*h4b+h4j]-co[3*h4a+h4j];
                h4e2[h4j]=co[3*h4n+h4j]-co[3*h4a+h4j];
              }
              h4cr[0]=h4e1[1]*h4e2[2]-h4e1[2]*h4e2[1];
              h4cr[1]=h4e1[2]*h4e2[0]-h4e1[0]*h4e2[2];
              h4cr[2]=h4e1[0]*h4e2[1]-h4e1[1]*h4e2[0];
              h4nc=sqrt(h4cr[0]*h4cr[0]+h4cr[1]*h4cr[1]+h4cr[2]*h4cr[2]);
              h4n1=sqrt(h4e1[0]*h4e1[0]+h4e1[1]*h4e1[1]+h4e1[2]*h4e1[2]);
              if((h4nc<1.e-30)||(h4n1<1.e-30)) continue;
              h4area=0.5*h4nc;
              for(h4j=0;h4j<3;h4j++){
                h4nr[h4j]=h4cr[h4j]/h4nc;
                h4t1[h4j]=h4e1[h4j]/h4n1;
              }
              h4t2[0]=h4nr[1]*h4t1[2]-h4nr[2]*h4t1[1];
              h4t2[1]=h4nr[2]*h4t1[0]-h4nr[0]*h4t1[2];
              h4t2[2]=h4nr[0]*h4t1[1]-h4nr[1]*h4t1[0];
              h4p=ielprop[h4i];
              h4kn=prop[h4p];h4tn=prop[h4p+1];h4ts=prop[h4p+2];h4gm=prop[h4p+4];
              if((h4kn<=0.)||(h4tn<=0.)) continue;
              h4bta=(h4ts/h4tn)*(h4ts/h4tn);
              for(h4j=0;h4j<3;h4j++)
                for(h4k=0;h4k<3;h4k++)
                  h4T[3*h4j+h4k]=h4nr[h4j]*h4nr[h4k]
                    +h4bta*(h4t1[h4j]*h4t1[h4k]+h4t2[h4j]*h4t2[h4k]);
              for(h4m=0;h4m<3;h4m++){
                h4g=1.-xstate[*nstate_*(mi[0]*h4i+h4m)+1];
                if(h4g<h4gm) h4g=h4gm;
                for(h4n=0;h4n<6;h4n++){
                  h4sh=(h4n%3==h4m)?(2./3.):(1./6.);
                  h4k=kon[ipkon[h4i]+h4n]-1;
                  if((h4k<0)||(h4k>=*nk)) continue;
                  h4w=h4area/3.*h4sh*h4sh*h4kn;
                  for(h4j=0;h4j<9;h4j++){
                    h4loc[9*h4k+h4j]+=h4w*h4g*h4T[h4j];
                    h4int[9*h4k+h4j]+=h4w*h4T[h4j];
                  }
                }
              }
            }

            /* lambda_min of each 3x3 by cyclic Jacobi, for bare nodes only */
            h4tot=0;
            for(h4i=0;h4i<*nk;h4i++){
              if(damage_free_nb[h4i]>0) continue;
              if(h4int[9*h4i]+h4int[9*h4i+4]+h4int[9*h4i+8]<=0.) continue;
              h4tot++;
              for(h4b=0;h4b<2;h4b++){
                for(h4j=0;h4j<9;h4j++)
                  h4A[h4j]=(h4b==0)?h4loc[9*h4i+h4j]:h4int[9*h4i+h4j];
                for(h4n=0;h4n<20;h4n++){
                  h4off=0.;
                  for(h4j=0;h4j<3;h4j++)
                    for(h4k=0;h4k<3;h4k++)
                      if(h4j!=h4k) h4off+=h4A[3*h4j+h4k]*h4A[3*h4j+h4k];
                  if(h4off<1.e-30) break;
                  for(h4j=0;h4j<2;h4j++){
                    for(h4k=h4j+1;h4k<3;h4k++){
                      if(fabs(h4A[3*h4j+h4k])<1.e-300) continue;
                      h4th=(h4A[3*h4k+h4k]-h4A[3*h4j+h4j])
                           /(2.*h4A[3*h4j+h4k]);
                      h4t=(h4th>=0.?1.:-1.)/(fabs(h4th)+sqrt(h4th*h4th+1.));
                      h4c=1./sqrt(h4t*h4t+1.);h4s=h4t*h4c;
                      for(h4m=0;h4m<3;h4m++){
                        h4tau=h4A[3*h4j+h4m];
                        h4A[3*h4j+h4m]=h4c*h4tau-h4s*h4A[3*h4k+h4m];
                        h4A[3*h4k+h4m]=h4s*h4tau+h4c*h4A[3*h4k+h4m];
                      }
                      for(h4m=0;h4m<3;h4m++){
                        h4tau=h4A[3*h4m+h4j];
                        h4A[3*h4m+h4j]=h4c*h4tau-h4s*h4A[3*h4m+h4k];
                        h4A[3*h4m+h4k]=h4s*h4tau+h4c*h4A[3*h4m+h4k];
                      }
                    }
                  }
                }
                h4ev[h4b]=h4A[0];
                if(h4A[4]<h4ev[h4b]) h4ev[h4b]=h4A[4];
                if(h4A[8]<h4ev[h4b]) h4ev[h4b]=h4A[8];
              }
              if(h4ev[1]>0.) h4rat[h4i]=h4ev[0]/h4ev[1];
            }

            /* distribution, not a threshold */
            if(h4tot>0){
              ITG h4bin[7];double h4edge[6];
              h4edge[0]=1.e-6;h4edge[1]=1.e-5;h4edge[2]=1.e-4;
              h4edge[3]=1.e-3;h4edge[4]=1.e-2;h4edge[5]=1.e-1;
              for(h4j=0;h4j<7;h4j++) h4bin[h4j]=0;
              h4r=2.;h4cnt=-1;
              for(h4i=0;h4i<*nk;h4i++){
                if(h4rat[h4i]<0.) continue;
                for(h4j=0;h4j<6;h4j++) if(h4rat[h4i]<h4edge[h4j]) break;
                h4bin[h4j]++;
                if(h4rat[h4i]<h4r){h4r=h4rat[h4i];h4cnt=h4i;}
              }
              printf("[H4 SUPPORT] inc=%" ITGFORMAT " time=%.6e bare_nodes=%"
                     ITGFORMAT "\n",iinc,theta**tper,h4tot);
              printf("[H4 SUPPORT]   lambda_min(K)/lambda_min(K_intact) bins:"
                     " <1e-6:%" ITGFORMAT " <1e-5:%" ITGFORMAT
                     " <1e-4:%" ITGFORMAT " <1e-3:%" ITGFORMAT
                     " <1e-2:%" ITGFORMAT " <1e-1:%" ITGFORMAT
                     " >=1e-1:%" ITGFORMAT "\n",
                     h4bin[0],h4bin[1],h4bin[2],h4bin[3],h4bin[4],h4bin[5],
                     h4bin[6]);
              if(h4cnt>=0)
                printf("[H4 SUPPORT]   worst node %" ITGFORMAT
                       " at (%.3f, %.3f, %.3f) ratio=%.4e\n",
                       h4cnt+1,co[3*h4cnt],co[3*h4cnt+1],co[3*h4cnt+2],h4r);
              fflush(stdout);
            }
            SFREE(h4rat);SFREE(h4int);SFREE(h4loc);SFREE(damage_free_nb);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
            damage_scan_count++;
          }
        }

        /* DE1 changes the constitutive stiffness but deliberately keeps
           the element topology and sparse pattern fixed.  If calcdamage
           requests another equilibrium and no topology changed, simply
           repeat Newton at the SAME load level. */
        if(damage_scan_count==0){
          damage_soft_reeq=1;

          printf("[DAMAGE DE1] pass=%" ITGFORMAT
                 " inc=%" ITGFORMAT " time=%.12e changed=%" ITGFORMAT
                 " action=same-load-reequilibrate\n",
                 damage_active_pass,iinc,theta**tper,idamage);
          printf("[DAMAGE DE1 STATUS] pass=%" ITGFORMAT
                 " active=%" ITGFORMAT " D>0.1=%" ITGFORMAT
                 " D>0.5=%" ITGFORMAT " D>0.9=%" ITGFORMAT
                 " Dfull=%" ITGFORMAT " Dmax=%.6e max_dD=%.6e "
                 "changed_tol=%" ITGFORMAT "\n",
                 damage_active_pass,de1.nactive,de1.gt01,
                 de1.gt05,de1.gt09,de1.nfull,
                 de1.dmax,de1.maxdelta,
                 de1.nchanged);
          fflush(stdout);

          /* DE1.1 closes the staggered fixed point on the largest scalar
             damage correction, not on the number of elements whose last
             bits changed.  The latter can stay O(1e5) for many passes even
             when every individual correction is negligible. */
          if(de1.maxdelta<=DAMAGE_DE1_FP_TOL){
            printf("[DAMAGE DE1 FP CONVERGED] inc=%" ITGFORMAT
                   " pass=%" ITGFORMAT " max_dD=%.6e tol=%.6e\n",
                   iinc,damage_active_pass,de1.maxdelta,
                   DAMAGE_DE1_FP_TOL);
            fflush(stdout);
            goto damage_active_set_closed;
          }

          if(damage_active_pass>=DAMAGE_DE1_MAX_PASSES){

            if(de1.maxdelta<=DAMAGE_DE1_FP_RELAX_TOL){
              printf("[DAMAGE DE1 FP ACCEPT] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT
                     " max_dD=%.6e relaxed_tol=%.6e\n",
                     iinc,damage_active_pass,de1.maxdelta,
                     DAMAGE_DE1_FP_RELAX_TOL);
              fflush(stdout);
              goto damage_active_set_closed;
            }

            /* Do not spend an unbounded number of PARDISO factorizations
               at one load level.  Reject this physical increment and let
               the existing transactional rollback retry a smaller step. */
            theta=thetadamage;
            dtheta=DAMAGE_DE1_CUTBACK_FACTOR*dthetadamage;
            if(dtheta<*tmin) dtheta=*tmin;
            dthetaref=dtheta;
            istab=0;
            icutb=1;

            printf("[DAMAGE DE1 CUTBACK] inc=%" ITGFORMAT
                   " pass=%" ITGFORMAT " max_dD=%.6e > %.6e "
                   "old_dt=%.6e retry_dt=%.6e\n",
                   iinc,damage_active_pass,de1.maxdelta,
                   DAMAGE_DE1_FP_RELAX_TOL,dthetadamage**tper,
                   dtheta**tper);
            fflush(stdout);
            goto damage_controller_done;
          }

          theta=thetadamage;
          dtheta=dthetadamage;
          dthetaref=dthetarefdamage;
          idiscon=1;

          /* idamagereeq stays set: after convergence calcdamagebase will
             recompute DE1 from the immutable physical-increment baseline
             and close the constitutive fixed point. */
          continue;
        }

        damage_soft_reeq=0;

        /* Rebuild the tentative history from the complete difference
           between the physical-increment baseline topology and the current
           active-set topology.  Earlier deleted elements retain the damage
           value stored when they were removed because calcdamage skips
           ipkon<0 elements. */

        topo_txn_discard(&dtxn);

        /* [TOPOLOGY] the block that was written out three times. */
        topo_txn_collect(&dtxn,*istep,iinc,theta**tper,*ttime+theta**tper,
                         ne0,ipkondamageini,ipkon,ielmat,mi,lakon,dam,
                         ndmcon,dmcon,*ndmat_,*ntmat_,
                         damage_de13_transaction,damage_de13_trigger_value,
                         damage_de13_trigger_ip);

        printf("[DAMAGE ACTIVESET] pass=%" ITGFORMAT
               " inc=%" ITGFORMAT " time=%.12e new_deleted=%" ITGFORMAT
               " total_tentative=%" ITGFORMAT "\n",
               damage_active_pass,iinc,theta**tper,idamage,
               dtxn.count);
        fflush(stdout);

        /* Which elements are actually in the batch.  Four hypotheses for
           the stall at lambda=0.3579 died for want of exactly this: the
           batch is what decides whether the surviving material still
           carries load, and the batch is the one thing the log never
           said. */
        if(damage_batch_list==1){
          printf("[DAMAGE BATCH] pass=%" ITGFORMAT " inc=%" ITGFORMAT
                 " elements:",damage_active_pass,iinc);
          for(k=0;k<dtxn.count;k++){
            printf(" %" ITGFORMAT,dtxn.elem[k]);
          }
          printf("\n");
          fflush(stdout);
        }

        /* The batch as a SET: sorted, hashed, and printed in full.  Two
           batches with the same hash are the same set of elements
           whatever order the scan produced them in. */

        if((td_armed!=0)&&(td_trace!=0)){
          td_batch=topodiag_hash_batch(dtxn.elem,dtxn.count,
                                       td_sort,topodiag_hash_seed());
          printf("[BATCHTRACE] batch inc=%" ITGFORMAT " icutb=%" ITGFORMAT
                 " pass=%" ITGFORMAT " n=%" ITGFORMAT
                 " committed_state=%016llx batch=%016llx sorted:",
                 iinc,icutb,damage_active_pass,dtxn.count,
                 td_state,td_batch);
          for(k=0;k<dtxn.count;k++)
            printf(" %" ITGFORMAT,td_sort[k]);
          printf("\n");
          fflush(stdout);
        }

        /* BK4 deferred sparse-structure compaction.  A terminal element is
           already absent from assembly because ipkon<0.  Keeping the old
           sparse graph is algebraically harmless while every previously
           active node still belongs to at least one active element: the
           obsolete entries are simply assembled as zero.  Rebuild as soon
           as an orphan node appears.  Detached but internally connected
           components remain the responsibility of the planned connectivity
           manager; a failed same-load solve still follows the full rollback
           path. */

        damage_topology_rebuild=1;
        damage_topology_orphans=0;
        if((damage_topology_deferred_mode==1)&&
           (damage_de13_transaction==1)){
          NNEW(damage_iponoel_trial,ITG,*nk);
          ITGMEMSET(damage_iponoel_trial,0,*nk,0);
          FORTRAN(nodebelongstoel,(damage_iponoel_trial,lakon,ipkon,kon,ne));
          NNEW(damage_orphan_seen,ITG,*nk);
          for(k=0;k<dtxn.count;k++){
            i=dtxn.elem[k]-1;
            if((i<0)||(i>=ne0)||(ipkon[i]>=0)) continue;
            if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
            damage_indexe=-ipkon[i]-2;
            for(j=0;j<4;j++){
              inode=kon[damage_indexe+j]-1;
              if((inode<0)||(inode>=*nk)) continue;
              if((damage_iponoel_trial[inode]==0)&&
                 (damage_orphan_seen[inode]==0)){
                for(idir=1;idir<=3;idir++){
                  if(nactdof[mt*inode+idir]>0){
                    damage_topology_orphans++;
                    damage_orphan_seen[inode]=1;
                    break;
                  }
                }
              }
            }
          }
          if(damage_topology_orphans==0){
            isiz=*nk;
            cpyparitg(iponoel,damage_iponoel_trial,&isiz,&num_cpus);
            damage_topology_rebuild=0;
          }
          SFREE(damage_orphan_seen);
          SFREE(damage_iponoel_trial);
        }

        printf("[DAMAGE TOPOLOGY BK4] inc=%" ITGFORMAT
               " pass=%" ITGFORMAT " tentative=%" ITGFORMAT
               " orphan_nodes=%" ITGFORMAT " action=%s\n",
               iinc,damage_active_pass,dtxn.count,
               damage_topology_orphans,
               damage_topology_rebuild?"remastruct":"reuse-sparse-graph");
        fflush(stdout);

        /* [DAMAGE RELEASE] arm.  f still holds f_int(u*) of the converged
           state - the elements were marked above but no results() has run
           since - and nactdof is still the pre-remastruct numbering.  Both
           are mapped into node space here so they survive remastruct. */
        if(damage_release_probe){
          ITG ai,aj,ak;
          if(damage_frel==NULL){
            NNEW(damage_frel,double,mt**nk);
            NNEW(damage_ract,ITG,mt**nk);
          }
          for(ai=0;ai<*nk;ai++){
            for(aj=0;aj<mt;aj++){
              ak=nactdof[mt*ai+aj];
              damage_frel[mt*ai+aj]=(ak>0)?f[ak-1]:0.;
              damage_ract[mt*ai+aj]=(ak>0)?1:0;
            }
          }
          damage_release_armed=1;
          damage_release_pass++;
          damage_release_qa=qa[0];
          damage_release_qam=qam[0];
          damage_release_dt=dtime;
          damage_release_rebuild=damage_topology_rebuild;
          damage_release_iforbou=iforbou;
          damage_release_nterm=damage_ebatch.terminal;
          damage_release_nother=damage_ebatch.marked-damage_ebatch.terminal;
          damage_release_nisl=damage_float_isl;
          damage_release_ncoh=damage_float_coh;

          /* level 2: one line per element in this batch.  The degradation it
             carried when it left is printed as a REFERENCE CHARACTERISTIC -
             it is a dimensionless multiplier, not a force, and it is NOT
             summed into any estimate of the released force.  The measured
             force is the surv/removed split above.

             It does classify the source for free, though: the terminal
             trigger cannot fire below its own threshold, so an element that
             left at Dvis < delete_d did NOT come from it - it came from
             damfloat, which reads no damage variable at all and therefore
             removes at whatever the element was carrying. */
          if(damage_release_probe>=2){
            ITG pe,pj,pnip,pel;
            double pd,pdv,pdmax,pdvmax,ptrig;
            for(pe=0;pe<dtxn.count;pe++){
              pel=dtxn.elem[pe]-1;
              if((pel<0)||(pel>=ne0)) continue;
              pnip=topo_element_nip(&lakon[8*pel],mi[0]);
              if(pnip<1) pnip=1;
              if(pnip>mi[0]) pnip=mi[0];
              pdmax=0.; pdvmax=0.;
              for(pj=0;pj<pnip;pj++){
                pd=dam[mi[0]*pel+pj]-1.;
                if(pd<0.) pd=0.;
                if(pd>1.) pd=1.;
                if(pd>pdmax) pdmax=pd;
                if(damage_damvisc!=NULL){
                  pdv=damage_damvisc[mi[0]*pel+pj];
                  if(pdv<0.) pdv=0.;
                  if(pdv>1.) pdv=1.;
                  if(pdv>pdvmax) pdvmax=pdv;
                }
              }
              ptrig=((damage_epol.delete_visc==1)&&(damage_damvisc!=NULL))
                    ?pdvmax:pdmax;
              /* The bulk damage variable exists ONLY for bulk elements.
                 calcdamage.f:135 skips every lakon(1:1) != 'C', so dam is
                 identically zero for a UC6 facet - printing it would read
                 as "left fully intact" when in truth the facet was
                 conducting gmin*Kn and its state lives in xstate, not here.
                 Both decks in play carry UC6, so this is not hypothetical:
                 damfloatcoh removals land in the same tentative batch. */
              if(lakon[8*pel]!='C'){
                printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                       " pass=%" ITGFORMAT " el=%" ITGFORMAT
                       " mat=%" ITGFORMAT " type=%.8s"
                       " D=n/a Dvis=n/a g_ref=n/a"
                       " note=non-bulk-element-damage-lives-in-xstate%s",
                       iinc,damage_release_pass,dtxn.elem[pe],
                       dtxn.mat[pe],&lakon[8*pel],
                       "\n");
                continue;
              }
              printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT " el=%" ITGFORMAT
                     " mat=%" ITGFORMAT " type=%.8s"
                     " D=%.6f Dvis=%.6f g_ref=%.4e"
                     " below_delete_d=%s%s",
                     iinc,damage_release_pass,dtxn.elem[pe],
                     dtxn.mat[pe],&lakon[8*pel],pdmax,pdvmax,1.-ptrig,
                     (ptrig<damage_epol.delete_d)?"YES-not-terminal":"no",
                     "\n");
            }
            fflush(stdout);
          }
        }

        if(damage_topology_rebuild){
          iitsav=iit;
          iit=-2;

          remastruct(ipompc,&coefmpc,&nodempc,nmpc,
                     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
                     ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
                     &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
                     &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
                     &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
                     mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
                     itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
                     tieset,ntie,&num_cpus,ielmat,matname);

          iit=iitsav;

          SFREE(nactdofinv);
          NNEW(nactdofinv,ITG,mt**nk);
          MNEW(nodorig,ITG,*nk);

          FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
                                 ipkon,lakon,kon,ne));

          SFREE(nodorig);

          ITGMEMSET(iponoel,0,*nk,0);
          FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

        /* checkconvergence() advanced theta when the just-finished
           re-equilibration converged.  Rewind to the saved beginning so
           the next Newton solve reaches exactly the same physical load. */

        theta=thetadamage;
        dtheta=dthetadamage;
        dthetaref=dthetarefdamage;
        idiscon=1;

        /* idamagereeq deliberately remains 1: the active set is still
           open and the next converged Newton solve must be checked again. */
        continue;
      }

damage_active_set_closed:

      /* No additional element failed after recomputing damage from the
         final post-deletion equilibrium.  The damage/topology active set
         is closed and the whole tentative deletion set may now be
         committed atomically to jobname.damage. */

      if(dtxn.count>0){
        damage_batch++;

        /* [TOPOLOGY] step C.  This block appeared once, so the reason to
           move it is not duplication: the nine-field record it writes is
           parsed by tools/ccxdiff.py, so the format had two
           implementations in this repository and an owner for neither.
           A writer that gains a field the reader does not expect makes
           ccxdiff compare the wrong columns - a silent failure in the
           instrument rather than in the thing measured. */

        topo_txn_commit(&dtxn,fdamage,damage_batch,
                        damage_de13_transaction,damage_active_pass);

        /* Configurable fracture termination.

           CCX_FRACTURE_TERMINATION="SETA:SETB" names two node sets whose
           load path is the thing being destroyed:

             tension   Face_X0_nset:Face_XL_nset
             Lame      InnerPressure_nset:OuterRadius_nset

           Once no chain of surviving elements links them, the specimen
           has separated and every further increment is a zero-load walk
           to the end of the step.  Stopping here is a demonstrated loss
           of the load-bearing path, not an increment-size failure. */

        if((damage_fracture_seta!=NULL)&&(damage_fracture_complete==0)){
          ITG *damage_ifacdead=NULL,ifd_i,ifd_n=0;
          NNEW(damage_ifacdead,ITG,*ne);
          if(damage_deadfacet&&(*nstate_>=4)){
            /* [DAMSTATE] the same judgement, asked of the one owner.  The
               three-point count and the reason it must not be mi[0] live in
               damstate_facet_dead, with a test on both branches. */
            for(ifd_i=0;ifd_i<*ne;ifd_i++){
              if(ipkon[ifd_i]<0) continue;
              if(lakon[8*ifd_i]!='U') continue;
              if(damstate_facet_dead(xstate,*nstate_,mi[0],ifd_i,3)){
                damage_ifacdead[ifd_i]=1;
                ifd_n++;
              }
            }
          }
          damage_conn=1;
          damage_conn_reach=0;
          FORTRAN(damconnectsets,(ipkon,kon,lakon,ne,nk,set,nset,
                                  istartset,iendset,ialset,
                                  damage_fracture_a,damage_fracture_b,
                                  &damage_conn,&damage_conn_reach,
                                  &damage_fracture_link,damage_ifacdead));

          /* the same question asked as a width.  It runs beside the
             boolean rather than instead of it, so the boolean's answer -
             and every baseline that depends on it - is untouched. */

          if((damage_cut_on==1)&&(damage_cut_bad==0)){
            ITG *lc_a=NULL,*lc_b=NULL,lc_na=0,lc_nb=0,lc_nf=0,lc_nel=0,
                lc_bel=0,lc_ex=1,lc_w=0;
            double lc_target=-1.;
            if((loadcut_sets(set,*nset,istartset,iendset,ialset,
                             damage_fracture_seta,&lc_a,&lc_na)==1)&&
               (loadcut_sets(set,*nset,istartset,iendset,ialset,
                             damage_fracture_setb,&lc_b,&lc_nb)==1)){
              if((damage_cut_ref>0.)&&(damage_cut_exact==0))
                lc_target=damage_cut_frac*damage_cut_ref;
              damage_cut_now=loadcut_width(co,ipkon,kon,lakon,*ne,*nk,
                                           lc_a,lc_na,lc_b,lc_nb,dam,mi,
                                           damage_cut_gmin,xstate,*nstate_,
                                           damage_ifacdead,lc_target,
                                           &lc_nf,&lc_nel,&lc_bel,&lc_ex,
                                           &lc_w);
              if(damage_cut_now<0.){
                damage_cut_bad=1;
              }else{
                if(damage_cut_ref<0.){
                  damage_cut_ref=damage_cut_now;
                  printf("[LOADCUT] reference cut at the first committed "
                         "batch: %.6e over %" ITGFORMAT " element(s) and %"
                         ITGFORMAT " shared face(s)\n",
                         damage_cut_ref,lc_nel,lc_nf);
                }
                /* an early exit returns a LOWER BOUND, not the cut, and
                   the report has to say which */
                printf("[LOADCUT] inc=%" ITGFORMAT " cut%s%.6e ratio%s%.6e "
                       "conn=%" ITGFORMAT "\n",iinc,
                       lc_ex?"=":">=",
                       lc_ex?damage_cut_now:lc_target,
                       lc_ex?"=":">=",
                       (damage_cut_ref>0.)?
                         (lc_ex?damage_cut_now/damage_cut_ref:damage_cut_frac)
                         :1.,
                       damage_conn);
                if((damage_cut_ref>0.)&&(damage_cut_narrow==0)&&
                   ((lc_bel==1)||((damage_cut_exact==1)&&
                     (damage_cut_now<damage_cut_frac*damage_cut_ref)))){
                  damage_cut_narrow=1;
                  damage_fracture_complete=1;
                  printf("\n[FRACTURE COMPLETE] inc=%" ITGFORMAT
                         " step_time=%.12e\n"
                         "                    the load path between %s and "
                         "%s has narrowed to %.4g of its\n"
                         "                    original width, below the "
                         "%.4g the deck asked for.  The topological\n"
                         "                    sweep still calls them %s.\n\n",
                         iinc,theta**tper,damage_fracture_seta,
                         damage_fracture_setb,
                         damage_cut_now/damage_cut_ref,damage_cut_frac,
                         damage_conn?"CONNECTED":"disconnected");
                }
                fflush(stdout);
              }
            }else{
              damage_cut_bad=1;
              printf("[LOADCUT] *WARNING: could not resolve the "
                     "termination node sets; the width is not reported.\n");
            }
            free(lc_a); free(lc_b);   /* plain malloc in loadcut.c */
          }

          if(damage_conn==0){
            damage_fracture_complete=1;
            printf("\n[FRACTURE COMPLETE] inc=%" ITGFORMAT
                   " step_time=%.12e\n"
                   "                    no surviving load path between "
                   "%s and %s\n"
                   "                    nodes reachable from the first "
                   "set: %" ITGFORMAT "\n\n",
                   iinc,theta**tper,damage_fracture_seta,
                   damage_fracture_setb,damage_conn_reach);
            if(ifd_n>0){
              printf("                    excluded %" ITGFORMAT
                     " fully failed cohesive facet(s)\n",ifd_n);
            }
            fflush(stdout);
          }
          SFREE(damage_ifacdead);
        }

        topo_txn_discard(&dtxn);
      }

      if(damage_soft_reeq==1){
        printf("[DAMAGE DE1 CONVERGED] inc=%" ITGFORMAT
               " time=%.12e passes=%" ITGFORMAT
               " active=%" ITGFORMAT " Dmax=%.6e max_dD=%.6e\n",
               iinc,theta**tper,damage_active_pass,
               de1.nactive,de1.dmax,
               de1.maxdelta);

        damstats_append(jobnamec,*istep,iinc,theta**tper,
                                *ttime+theta**tper,damage_active_pass,
                                de1.nactive,de1.gt01,
                                de1.gt05,de1.gt09,
                                de1.nfull,de1.dmax,
                                de1.maxdelta);

        damstats_write_vtk(jobnamec,co,vold,*nk,mt,kon,ipkon,lakon,
                             ielmat,mi[2],dam,mi[0],ne0,*istep,iinc,
                             theta**tper);

        printf("[DAMAGE DE1 OUTPUT] exact cell snapshot: %s.de1.vtk; "
               "history: %s.de1stats\n",jobnamec,jobnamec);
        fflush(stdout);
      }

      damage_soft_reeq=0;
      damage_de13_transaction=0;

      /* If a bounded fast event had to cut the physical increment,
         do not regrow dtheta through many 1.5x increments.  The local
         controller below will TRY the whole remaining distance to
         theta_goal in one Newton solve.  Failure is harmless: ordinary
         CalculiX cutback remains in charge. */
      if((damage_fast_used==1)&&(ilocalsubstep==1))
        damage_fast_recover=1;

      idamagereeq=0;

    }else if((icutb==0)&&(*ndmat_>0)){
      damage_event_cut=0;
      damage_mode=0;
      damage_predict_count=0;
      damage_alphaevent=2.;

      /* Use the immutable damage state at the beginning of the physical
         increment as the baseline.  In explicit dynamics no rollback
         baseline is allocated, so the current committed dam is used. */

      if((*iexpl<=1)&&(damdamageini!=NULL)){
        damagebase=damdamageini;
      }else{
        damagebase=dam;
      }

      /* First pass: predict whether the accumulated damage reaches the
         deletion limit inside this already-converged physical increment.
         This pass does not modify dam or ipkon. */

      FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
			  ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
			  ndmcon,dmcon,dam,damagebase,&dtime,sti,ithermal,t1,xstate,
			  xstateini,nstate_,vold,&damage_predict_count,
			  &damage_mode,&damage_alphaevent));

      /* Adaptive-fast event controller.

         A2 localized essentially every first crossing, which is robust but
         can collapse dtheta by factors of 50-100 when many C3D4 elements
         approach the limit one after another.

         A3 treats event placement as a TRIAL only:
           1) if the already-converged coarse increment predicts a modest
              batch and the first event is not too early, accept the coarse
              endpoint and let the transactional active-set test the batch;
           2) otherwise localize, but never below a fixed fraction of the
              current physical increment on the first fast attempt;
           3) if the post-deletion Newton solve fails, the existing rollback
              restores topology/state and damage_fast_retry forces the next
              attempt back to the exact A2 locator (no floor, no batching).

         Thus aggressive batching can cost at most one failed trial.  It can
         never be committed unless same-load equilibrium and active-set
         closure both converge. */

      /* damage_fast_used is intentionally NOT cleared here: when a
         bounded-localization trial is rolled back and re-solved, the flag
         must survive until the resulting deletion transaction commits so
         the fast recovery-to-theta_goal logic can be armed. */

      if((damage_predict_count>0)&&(*iexpl<=1)&&
         (damage_alphaevent>1.e-8)&&(damage_alphaevent<0.95)&&
         (dthetadamage>1.01**tmin)){

        damage_event_raw=1.02*damage_alphaevent*dthetadamage;
        damage_event_dtheta=damage_event_raw;

        if(damage_fast_retry==0){

          /* Fast path A: use the already-converged coarse endpoint for a
             bounded predicted batch.  Actual deletion is still based on the
             re-solved/current damage state, not on predictor membership. */

          if((damage_predict_count<=DAMAGE_FAST_BATCH_MAX)&&
             (damage_alphaevent>=DAMAGE_FAST_DIRECT_ALPHA)){

            damage_fast_used=1;

            printf("[DAMAGE FAST] inc=%" ITGFORMAT
                   " action=coarse-batch alpha=%.6e dt=%.6e "
                   "predicted=%" ITGFORMAT "\n",
                   iinc,damage_alphaevent,dthetadamage**tper,
                   damage_predict_count);
            fflush(stdout);

          }else{

            /* Fast path B: bounded localization.  This prevents a single
               early crossing from forcing an extremely small physical step.
               Any excessive topology jump is protected by transactional
               rollback and the exact-locator retry below. */

            damage_event_floor=DAMAGE_FAST_MIN_FRACTION*dthetadamage;
            if(damage_event_dtheta<damage_event_floor)
              damage_event_dtheta=damage_event_floor;

            if(damage_event_dtheta<*tmin) damage_event_dtheta=*tmin;
            if(damage_event_dtheta>0.95*dthetadamage)
              damage_event_dtheta=0.95*dthetadamage;

            if(damage_event_dtheta<0.999*dthetadamage){
              damage_fast_used=1;

              printf("[DAMAGE FAST] inc=%" ITGFORMAT
                     " action=bounded-localize alpha=%.6e old_dt=%.6e "
                     "raw_dt=%.6e trial_dt=%.6e predicted=%" ITGFORMAT
                     "\n",
                     iinc,damage_alphaevent,dthetadamage**tper,
                     damage_event_raw**tper,damage_event_dtheta**tper,
                     damage_predict_count);
              fflush(stdout);

              theta=thetadamage;
              dtheta=damage_event_dtheta;
              dthetaref=dtheta;
              istab=0;
              icutb=1;
              damage_event_cut=1;
            }
          }

        }else{

          /* Safety path after a failed fast post-deletion equilibrium:
             restore the exact A2 locator.  No artificial lower bound is
             applied. */

          if(damage_event_dtheta<*tmin) damage_event_dtheta=*tmin;
          if(damage_event_dtheta>0.95*dthetadamage)
            damage_event_dtheta=0.95*dthetadamage;

          if(damage_event_dtheta<0.999*dthetadamage){
            printf("[DAMAGE SAFE RETRY] inc=%" ITGFORMAT
                   " alpha=%.6e old_dt=%.6e exact_dt=%.6e "
                   "predicted=%" ITGFORMAT " fast_failures=%" ITGFORMAT
                   "\n",
                   iinc,damage_alphaevent,dthetadamage**tper,
                   damage_event_dtheta**tper,damage_predict_count,
                   damage_fast_failures);
            fflush(stdout);

            theta=thetadamage;
            dtheta=damage_event_dtheta;
            dthetaref=dtheta;
            istab=0;
            icutb=1;
            damage_event_cut=1;
          }
        }
      }

      if(damage_event_cut==0){

        /* Second pass: apply the damage increment and retain the original
           hard-deletion law. Because strongly interior crossings were
           localized above, the simultaneous deletion batch should now be
           much smaller and closer to the first physical damage event. */

        damage_mode=1;
        idamage=0;
        damage_alphaevent=2.;

        isiz=mi[0]*ne0;
        cpypardou(damde1prev,dam,&isiz,&num_cpus);

        FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
			    ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
			    ndmcon,dmcon,dam,damagebase,&dtime,sti,ithermal,t1,xstate,
			    xstateini,nstate_,vold,&idamage,
			    &damage_mode,&damage_alphaevent));

        damstats_element(dam,damde1prev,ipkon,lakon,ne0,mi[0],
                         &de1.nactive,&de1.gt01,
                         &de1.gt05,&de1.gt09,
                         &de1.nfull,&de1.nchanged,
                         &de1.dmax,&de1.maxdelta);

      if((idamage>0)&&(*iexpl<=1)){

        /* Build a tentative deletion batch by comparing the topology
           snapshot at the beginning of the physical increment with the
           topology immediately after calcdamage().  Nothing is written
           to disk here: the batch is committed only after successful
           same-load re-equilibration. */

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            stiffcensus damage_census;
            stiffcensus_take(&damage_census,*nk,damage_addiag,
                             damage_addiag0,damage_addok);
            if(damage_census_ok) monitor_stiffness(&damage_census,iinc,
                                                   theta**tper);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            NNEW(damage_free_nb,ITG,*nk);
            NNEW(damage_free_g,double,*nk);
            for(i=0;i<*nk;i++){damage_free_nb[i]=0;damage_free_g[i]=-1.;}
            for(i=0;i<*ne;i++){
              if(ipkon[i]<0) continue;
              if(lakon[8*i]=='C'){
                for(j=0;j<4;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k>=0)&&(k<*nk)) damage_free_nb[k]++;
                }
              }else if(lakon[8*i]=='U'){
                damage_free_gm=0.;
                for(j=0;j<3;j++){
                  damage_free_dv=xstate[*nstate_*(mi[0]*i+j)+1];
                  if(1.-damage_free_dv>damage_free_gm)
                    damage_free_gm=1.-damage_free_dv;
                }
                for(j=0;j<6;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k<0)||(k>=*nk)) continue;
                  if(damage_free_gm>damage_free_g[k])
                    damage_free_g[k]=damage_free_gm;
                }
              }
            }
            damage_free_cnt=0;
            for(i=0;i<*nk;i++){
              if(damage_free_nb[i]>0) continue;
              if(damage_free_g[i]<0.) continue;
              if(damage_free_g[i]<1.e-2) damage_free_cnt++;
            }
            if(damage_free_cnt>damage_free_rep){
              damage_free_rep=damage_free_cnt;
              printf("[DAMAGE FREE NODE] inc=%" ITGFORMAT " time=%.12e "
                     "%" ITGFORMAT " node(s) with no bulk AND every "
                     "remaining facet debonded (new maximum)\n",
                     iinc,theta**tper,damage_free_cnt);
              fflush(stdout);
            }
            SFREE(damage_free_nb);SFREE(damage_free_g);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
            damage_scan_count++;
          }
        }

        /* DE1: constitutive damage changed, but no element was physically
           deleted.  Keep the matrix graph intact and re-equilibrate the
           degraded material at the same physical load level. */
        if(damage_scan_count==0){
          damage_active_pass=1;
          damage_soft_reeq=1;
          dtxn.increment=iinc;
          dtxn.step_time=theta**tper;
          dtxn.total_time=*ttime+dtxn.step_time;

          printf("[DAMAGE DE1] pass=1 inc=%" ITGFORMAT
                 " time=%.12e changed=%" ITGFORMAT
                 " action=same-load-reequilibrate\n",
                 iinc,theta**tper,idamage);
          printf("[DAMAGE DE1 STATUS] pass=1 active=%" ITGFORMAT
                 " D>0.1=%" ITGFORMAT " D>0.5=%" ITGFORMAT
                 " D>0.9=%" ITGFORMAT " Dfull=%" ITGFORMAT
                 " Dmax=%.6e max_dD=%.6e changed_tol=%" ITGFORMAT "\n",
                 de1.nactive,de1.gt01,de1.gt05,
                 de1.gt09,de1.nfull,de1.dmax,
                 de1.maxdelta,de1.nchanged);
          fflush(stdout);

          theta=thetadamage;
          dtheta=dthetadamage;
          dthetaref=dthetarefdamage;
          idiscon=1;
          idamagereeq=1;
          continue;
        }

        damage_soft_reeq=0;

        /* A non-empty old buffer would indicate an internal logic error.
           Discard it rather than allowing trial data to leak into a later
           commit. */
        topo_txn_discard(&dtxn);

        /* [TOPOLOGY] the block that was written out three times. */
        topo_txn_collect(&dtxn,*istep,iinc,theta**tper,*ttime+theta**tper,
                         ne0,ipkondamageini,ipkon,ielmat,mi,lakon,dam,
                         ndmcon,dmcon,*ndmat_,*ntmat_,
                         damage_de13_transaction,damage_de13_trigger_value,
                         damage_de13_trigger_ip);

        /* First topology change at this physical load level. */
        damage_active_pass=1;

	iitsav=iit;
	iit=-2;

	remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
		   ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
		   &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
		   &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
		   &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
		   mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
		   itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
		   tieset,ntie,&num_cpus,ielmat,matname);

	iit=iitsav;

	SFREE(nactdofinv);
	NNEW(nactdofinv,ITG,mt**nk);
	MNEW(nodorig,ITG,*nk);

	FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			       ipkon,lakon,kon,ne));

	SFREE(nodorig);

	ITGMEMSET(iponoel,0,*nk,0);

	FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));

	theta=thetadamage;
	dtheta=dthetadamage;
	dthetaref=dthetarefdamage;
	idiscon=1;
	idamagereeq=1;

	continue;
      }
      }
    }

damage_controller_done:

    /* DE1.3 terminal failure manager.

       DE1.2 has already converged the continuous damage state inside Newton.
       Only now, at a converged physical load level, are near-fully degraded
       elements allowed to change topology.  The batch is tentative: after
       remastruct the model must converge again at exactly the same load.
       Failure enters the existing idamagereeq rollback path and restores the
       topology/damage state from the start of the physical increment. */
    /* Dissipation measure for path following.

       The loading here is displacement controlled, u_p = lambda*u_hat, so
       the Gutierrez dissipation collapses to a scalar built from the load
       factor and the reaction conjugate to the prescribed pattern:

           work      = integral P dlambda ~ (P_n+P_n+1)(l_n+1-l_n)/2
           stored    = (P_n+1 l_n+1 - P_n l_n)/2
           dissipated= (P_n l_n+1 - P_n+1 l_n)/2

       P = sum over the prescribed dofs of fn*xboun, because
       xbounact = lambda*xboun and fn is the reaction there.  resultsini.c
       sets calcul_fn=1 for the NLGEOM Newton path, so fn is current on
       every iteration and this costs one pass over nboun.

       This is the control variable a dissipation-controlled step would
       prescribe.  Measuring it first, before making it the control, keeps
       the two questions separate: is the quantity well behaved, and does
       driving it fix the limit point. */

    if((icutb==0)&&(idamagereeq==0)&&(*nmethod!=4)&&(lc.diss_report==1)){
      if(lc.diss_init==1){
        lc.diss_lamnow=(lc.arc==1)?lc.arc_lam:theta;
        lc.diss_dg=0.5*(lc.diss_pprev*lc.diss_lamnow-
                            lc.diss_p*lc.diss_lprev);
        lc.diss_total+=lc.diss_dg;

        if((lc.diss_ctrl==2)&&(lc.diss_engaged==0)&&
           (lc.diss_dg>DAMAGE_DISS_ENGAGE*lc.diss_target)&&
           ((lc.diss_engage_t<0.)||(theta>=lc.diss_engage_t))){
          lc.diss_engaged=1;
          printf("[DISSIPATION CONTROL] engaged at inc=%" ITGFORMAT
                 " lambda=%.6f dG=%.6e\n",iinc,theta,lc.diss_dg);
          fflush(stdout);
        }

        /* Dissipation-based step control.

           The measured history says the controller walks at full dtmax
           right up to the failure and never pre-emptively cuts: it has no
           way of knowing the dissipation is climbing, because the
           iteration count stays low until the step that cannot converge
           at all.  Sizing the next increment so its predicted dissipation
           stays near a target makes the step shrink as the crack starts
           to run, which is exactly the information the stock controller
           lacks.

           This is NOT path following: lambda still only increases, so a
           genuine snap-back remains out of reach.  It is the cheap half
           of the idea, worth measuring before rebuilding the Newton loop
           around a bordered system. */

        if((lc.diss_target>0.)&&(lc.diss_dg>1.e-30)&&
           (dtheta>0.)&&(*idrct==0)&&(lc.diss_step==1)){
          lc.diss_scale=lc.diss_target/lc.diss_dg;
          if(lc.diss_scale>DAMAGE_DISS_GROW)
            lc.diss_scale=DAMAGE_DISS_GROW;
          if(lc.diss_scale<DAMAGE_DISS_SHRINK)
            lc.diss_scale=DAMAGE_DISS_SHRINK;
          lc.diss_dtheta=dtheta*lc.diss_scale;
          if(lc.diss_dtheta>dthetaref) lc.diss_dtheta=dthetaref;
          if(lc.diss_dtheta<(*tmin)) lc.diss_dtheta=*tmin;
          if(lc.diss_dtheta<0.98*dtheta){
            printf("[DISSIPATION STEP] inc=%" ITGFORMAT
                   " dG=%.6e target=%.6e dtheta %.6e -> %.6e\n",
                   iinc,lc.diss_dg,lc.diss_target,dtheta,
                   lc.diss_dtheta);
            fflush(stdout);
          }
          dtheta=lc.diss_dtheta;
        }

        printf("[DISSIPATION] inc=%" ITGFORMAT " lambda=%.6f P=%.6e "
               "dG=%.6e G=%.6e\n",iinc,theta,lc.diss_p,
               lc.diss_dg,lc.diss_total);
        fflush(stdout);
      }
      lc.diss_init=1;
      lc.diss_lprev=lc.diss_lamnow;
      lc.diss_pprev=lc.diss_p;
    }

    if((icutb==0)&&(idamagereeq==0)){
      lc.path_retry=0;
      lc.path_lamcom=lc.path_lam;
      if(lc.path_desc>0) lc.path_desc--;
    /* Arm only when the step controller has actually run out of
       room.  Arming on icutb>=3 fired at increment 385 on an
       ordinary cutback the stock logic recovers from, derailed the
       path, and never reached the real stall at all. */
    }else if((lc.path_on>=2)&&(dtheta<10.*(*tmin))&&
             (lc.path_desc==0)){
      lc.path_desc=lc.path_nstep;
      lc.path_used=1;
      printf("[DAMAGE PATH] inc=%" ITGFORMAT " arming a %"
             ITGFORMAT "-increment descent at lambda=%.6f\n",
             iinc,lc.path_nstep,lc.path_lamcom);
      fflush(stdout);
    }
    if((damage_de12_enabled)&&(icutb==0)&&(idamagereeq==0)&&
       (damdamageini!=NULL)){

      erosion_mark(&damage_epol,&damage_ebatch,
                   dam,damage_damvisc,ipkon,lakon,kon,ielmat,matname,
                   ndmcon,dmcon,*ndmat_,*ntmat_,*nk,*ne,ne0,mi[0],mi[2],
                   damage_de13_trigger_value,damage_de13_trigger_ip,
                   iinc,theta**tper);

      if(damage_ebatch.marked>0){
        damage_de13_transaction=1;
        damage_soft_reeq=0;
        damage_active_pass=1;
        rsc.reeq_uam_ref[0]=uam[0];
        rsc.reeq_uam_ref[1]=uam[1];
        /* J-09: keep the reference from collapsing with the step.  The peak
           is a running maximum over the whole step, never reset per
           increment - the same shape as damage_qam_peak. */
        if(rsc.reeq_uam_floor>0.){
          if(rsc.reeq_uam_ref[0]>rsc.reeq_uam_peak[0])
            rsc.reeq_uam_peak[0]=rsc.reeq_uam_ref[0];
          if(rsc.reeq_uam_ref[1]>rsc.reeq_uam_peak[1])
            rsc.reeq_uam_peak[1]=rsc.reeq_uam_ref[1];
          if(rsc.reeq_uam_ref[0]<rsc.reeq_uam_floor*rsc.reeq_uam_peak[0])
            rsc.reeq_uam_ref[0]=rsc.reeq_uam_floor*rsc.reeq_uam_peak[0];
          if(rsc.reeq_uam_ref[1]<rsc.reeq_uam_floor*rsc.reeq_uam_peak[1])
            rsc.reeq_uam_ref[1]=rsc.reeq_uam_floor*rsc.reeq_uam_peak[1];
        }

        topo_txn_discard(&dtxn);

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            stiffcensus damage_census;
            stiffcensus_take(&damage_census,*nk,damage_addiag,
                             damage_addiag0,damage_addok);
            if(damage_census_ok) monitor_stiffness(&damage_census,iinc,
                                                   theta**tper);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            NNEW(damage_free_nb,ITG,*nk);
            NNEW(damage_free_g,double,*nk);
            for(i=0;i<*nk;i++){damage_free_nb[i]=0;damage_free_g[i]=-1.;}
            for(i=0;i<*ne;i++){
              if(ipkon[i]<0) continue;
              if(lakon[8*i]=='C'){
                for(j=0;j<4;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k>=0)&&(k<*nk)) damage_free_nb[k]++;
                }
              }else if(lakon[8*i]=='U'){
                damage_free_gm=0.;
                for(j=0;j<3;j++){
                  damage_free_dv=xstate[*nstate_*(mi[0]*i+j)+1];
                  if(1.-damage_free_dv>damage_free_gm)
                    damage_free_gm=1.-damage_free_dv;
                }
                for(j=0;j<6;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k<0)||(k>=*nk)) continue;
                  if(damage_free_gm>damage_free_g[k])
                    damage_free_g[k]=damage_free_gm;
                }
              }
            }
            damage_free_cnt=0;
            for(i=0;i<*nk;i++){
              if(damage_free_nb[i]>0) continue;
              if(damage_free_g[i]<0.) continue;
              if(damage_free_g[i]<1.e-2) damage_free_cnt++;
            }
            if(damage_free_cnt>damage_free_rep){
              damage_free_rep=damage_free_cnt;
              printf("[DAMAGE FREE NODE] inc=%" ITGFORMAT " time=%.12e "
                     "%" ITGFORMAT " node(s) with no bulk AND every "
                     "remaining facet debonded (new maximum)\n",
                     iinc,theta**tper,damage_free_cnt);
              fflush(stdout);
            }
            SFREE(damage_free_nb);SFREE(damage_free_g);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)) damage_scan_count++;
        }

        /* [TOPOLOGY] the block that was written out three times. */
        topo_txn_collect(&dtxn,*istep,iinc,theta**tper,*ttime+theta**tper,
                         ne0,ipkondamageini,ipkon,ielmat,mi,lakon,dam,
                         ndmcon,dmcon,*ndmat_,*ntmat_,
                         damage_de13_transaction,damage_de13_trigger_value,
                         damage_de13_trigger_ip);

        /* The batch as a SET: sorted, hashed and printed in full, next to
           the committed-state hash of the attempt that produced it.  Two
           batches with the same hash are the same set of elements
           whatever order the scan produced them in, and the PAIR
           (committed state, batch) is what would have to repeat for the
           batch/rollback loop to be real. */

        if((td_armed!=0)&&(td_trace!=0)&&(dtxn.count>0)){
          td_batch=topodiag_hash_batch(dtxn.elem,dtxn.count,
                                       td_sort,topodiag_hash_seed());
          printf("[BATCHTRACE] batch inc=%" ITGFORMAT " icutb=%" ITGFORMAT
                 " pass=%" ITGFORMAT " n=%" ITGFORMAT
                 " committed_state=%016llx batch=%016llx sorted:",
                 iinc,icutb,damage_active_pass,dtxn.count,
                 td_state,td_batch);
          for(k=0;k<dtxn.count;k++)
            printf(" %" ITGFORMAT,td_sort[k]);
          printf("\n");
          fflush(stdout);
        }

        damage_topology_rebuild=1;
        damage_topology_orphans=0;
        if(damage_topology_deferred_mode==1){
          NNEW(damage_iponoel_trial,ITG,*nk);
          ITGMEMSET(damage_iponoel_trial,0,*nk,0);
          FORTRAN(nodebelongstoel,(damage_iponoel_trial,lakon,ipkon,kon,ne));
          NNEW(damage_orphan_seen,ITG,*nk);
          for(k=0;k<dtxn.count;k++){
            i=dtxn.elem[k]-1;
            if((i<0)||(i>=ne0)||(ipkon[i]>=0)) continue;
            if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
            damage_indexe=-ipkon[i]-2;
            for(j=0;j<4;j++){
              inode=kon[damage_indexe+j]-1;
              if((inode<0)||(inode>=*nk)) continue;
              if((damage_iponoel_trial[inode]==0)&&
                 (damage_orphan_seen[inode]==0)){
                for(idir=1;idir<=3;idir++){
                  if(nactdof[mt*inode+idir]>0){
                    damage_topology_orphans++;
                    damage_orphan_seen[inode]=1;
                    break;
                  }
                }
              }
            }
          }
          if(damage_topology_orphans==0){
            isiz=*nk;
            cpyparitg(iponoel,damage_iponoel_trial,&isiz,&num_cpus);
            damage_topology_rebuild=0;
          }
          SFREE(damage_orphan_seen);
          SFREE(damage_iponoel_trial);
        }

        printf("[DAMAGE DE1.3 TERMINAL] inc=%" ITGFORMAT
               " time=%.12e new_terminal=%" ITGFORMAT
               " tentative=%" ITGFORMAT " batch_Dmax=%.6e "
               "batch_Dvis=%.6e "
               "orphan_nodes=%" ITGFORMAT " action=%s+same-load-Newton\n",
               iinc,theta**tper,damage_ebatch.marked,dtxn.count,
               damage_ebatch.batch_dmax,damage_ebatch.batch_vmin,
               damage_topology_orphans,
               damage_topology_rebuild?"remastruct":"reuse-sparse-graph");
        fflush(stdout);

        /* [DAMAGE RELEASE] arm.  f still holds f_int(u*) of the converged
           state - the elements were marked above but no results() has run
           since - and nactdof is still the pre-remastruct numbering.  Both
           are mapped into node space here so they survive remastruct. */
        if(damage_release_probe){
          ITG ai,aj,ak;
          if(damage_frel==NULL){
            NNEW(damage_frel,double,mt**nk);
            NNEW(damage_ract,ITG,mt**nk);
          }
          for(ai=0;ai<*nk;ai++){
            for(aj=0;aj<mt;aj++){
              ak=nactdof[mt*ai+aj];
              damage_frel[mt*ai+aj]=(ak>0)?f[ak-1]:0.;
              damage_ract[mt*ai+aj]=(ak>0)?1:0;
            }
          }
          damage_release_armed=1;
          damage_release_pass++;
          damage_release_qa=qa[0];
          damage_release_qam=qam[0];
          damage_release_dt=dtime;
          damage_release_rebuild=damage_topology_rebuild;
          damage_release_iforbou=iforbou;
          damage_release_nterm=damage_ebatch.terminal;
          damage_release_nother=damage_ebatch.marked-damage_ebatch.terminal;
          damage_release_nisl=damage_float_isl;
          damage_release_ncoh=damage_float_coh;

          /* level 2: one line per element in this batch.  The degradation it
             carried when it left is printed as a REFERENCE CHARACTERISTIC -
             it is a dimensionless multiplier, not a force, and it is NOT
             summed into any estimate of the released force.  The measured
             force is the surv/removed split above.

             It does classify the source for free, though: the terminal
             trigger cannot fire below its own threshold, so an element that
             left at Dvis < delete_d did NOT come from it - it came from
             damfloat, which reads no damage variable at all and therefore
             removes at whatever the element was carrying. */
          if(damage_release_probe>=2){
            ITG pe,pj,pnip,pel;
            double pd,pdv,pdmax,pdvmax,ptrig;
            for(pe=0;pe<dtxn.count;pe++){
              pel=dtxn.elem[pe]-1;
              if((pel<0)||(pel>=ne0)) continue;
              pnip=topo_element_nip(&lakon[8*pel],mi[0]);
              if(pnip<1) pnip=1;
              if(pnip>mi[0]) pnip=mi[0];
              pdmax=0.; pdvmax=0.;
              for(pj=0;pj<pnip;pj++){
                pd=dam[mi[0]*pel+pj]-1.;
                if(pd<0.) pd=0.;
                if(pd>1.) pd=1.;
                if(pd>pdmax) pdmax=pd;
                if(damage_damvisc!=NULL){
                  pdv=damage_damvisc[mi[0]*pel+pj];
                  if(pdv<0.) pdv=0.;
                  if(pdv>1.) pdv=1.;
                  if(pdv>pdvmax) pdvmax=pdv;
                }
              }
              ptrig=((damage_epol.delete_visc==1)&&(damage_damvisc!=NULL))
                    ?pdvmax:pdmax;
              /* The bulk damage variable exists ONLY for bulk elements.
                 calcdamage.f:135 skips every lakon(1:1) != 'C', so dam is
                 identically zero for a UC6 facet - printing it would read
                 as "left fully intact" when in truth the facet was
                 conducting gmin*Kn and its state lives in xstate, not here.
                 Both decks in play carry UC6, so this is not hypothetical:
                 damfloatcoh removals land in the same tentative batch. */
              if(lakon[8*pel]!='C'){
                printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                       " pass=%" ITGFORMAT " el=%" ITGFORMAT
                       " mat=%" ITGFORMAT " type=%.8s"
                       " D=n/a Dvis=n/a g_ref=n/a"
                       " note=non-bulk-element-damage-lives-in-xstate%s",
                       iinc,damage_release_pass,dtxn.elem[pe],
                       dtxn.mat[pe],&lakon[8*pel],
                       "\n");
                continue;
              }
              printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT " el=%" ITGFORMAT
                     " mat=%" ITGFORMAT " type=%.8s"
                     " D=%.6f Dvis=%.6f g_ref=%.4e"
                     " below_delete_d=%s%s",
                     iinc,damage_release_pass,dtxn.elem[pe],
                     dtxn.mat[pe],&lakon[8*pel],pdmax,pdvmax,1.-ptrig,
                     (ptrig<damage_epol.delete_d)?"YES-not-terminal":"no",
                     "\n");
            }
            fflush(stdout);
          }
        }

        if(damage_topology_rebuild){
          iitsav=iit;
          iit=-2;

          remastruct(ipompc,&coefmpc,&nodempc,nmpc,
                     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
                     ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
                     &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
                     &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
                     &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
                     mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
                     itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
                     tieset,ntie,&num_cpus,ielmat,matname);

          iit=iitsav;

          SFREE(nactdofinv);
          NNEW(nactdofinv,ITG,mt**nk);
          MNEW(nodorig,ITG,*nk);

          FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
                                 ipkon,lakon,kon,ne));

          SFREE(nodorig);

          ITGMEMSET(iponoel,0,*nk,0);
          FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

        theta=thetadamage;
        dtheta=dthetadamage;
        dthetaref=dthetarefdamage;
        idiscon=1;
        idamagereeq=1;
        continue;
      }
    }

    /* DE1.2 has already produced the converged trial damage inside the last
       Newton iteration.  There is no outer same-load damage fixed point to
       close.  Once the physical increment is accepted, report/serialize the
       exact integration-point state directly. */
    if((damage_de12_enabled)&&(icutb==0)&&(idamagereeq==0)&&
       (damdamageini!=NULL)){
      damstats_element(dam,damdamageini,ipkon,lakon,ne0,mi[0],
                       &de1.nactive,&de1.gt01,
                       &de1.gt05,&de1.gt09,
                       &de1.nfull,&de1.nchanged,
                       &de1.dmax,&de1.maxdelta);

      if(de1.nactive>0){
        printf("[DAMAGE DE1.2 COMMIT] inc=%" ITGFORMAT
               " time=%.12e active=%" ITGFORMAT
               " D>0.1=%" ITGFORMAT " D>0.5=%" ITGFORMAT
               " D>0.9=%" ITGFORMAT " Dfull=%" ITGFORMAT
               " Dmax=%.6e max_dD_inc=%.6e\n",
               iinc,theta**tper,de1.nactive,de1.gt01,
               de1.gt05,de1.gt09,de1.nfull,
               de1.dmax,de1.maxdelta);

        damstats_append(jobnamec,*istep,iinc,theta**tper,
                                *ttime+theta**tper,1,de1.nactive,
                                de1.gt01,de1.gt05,
                                de1.gt09,de1.nfull,
                                de1.dmax,de1.maxdelta);

        damstats_write_vtk(jobnamec,co,vold,*nk,mt,kon,ipkon,lakon,
                             ielmat,mi[2],dam,mi[0],ne0,*istep,iinc,
                             theta**tper);
        fflush(stdout);
      }
    }

    /* Local damage substepping.  checkconvergence() is left in charge
       of the actual cutback.  Here we only keep the original end point
       of the failed coarse increment and prevent the following successful
       reduced increments from stepping beyond it.  A physical substep is
       considered accepted only after any damage re-equilibration at its
       end has also converged. */

    if((icutb==0)&&(idamagereeq==0)&&(ilocalsubstep==1)){
      dtheta_remaining=theta_goal-theta;

      if((damage_fast_recover==1)&&(dtheta_remaining>1.e-12)&&
         (*itpamp==0)){

        /* Jump directly toward the saved coarse-increment goal instead of
           spending many accepted increments regrowing a tiny event step.
           This is only a proposed next step: Newton/checkconvergence can
           cut it back normally if the jump is too aggressive. */

        dtheta=dtheta_remaining;
        if(dtheta>dtheta_restore) dtheta=dtheta_restore;
        if(dtheta>*tmax) dtheta=*tmax;
        if(dtheta>1.-theta) dtheta=1.-theta;
        dthetaref=dtheta;
        istab=0;
        damage_fast_recover=0;

        printf("[DAMAGE FAST RECOVER] time=%e goal=%e next_dt=%e\n\n",
               theta**tper,theta_goal**tper,dtheta**tper);

      }else if(dtheta_remaining<=1.e-12){

        /* Eliminate round-off drift at the saved target. */

        theta=theta_goal;
        ilocalsubstep=0;

        /* Do not let the two-success CCX growth history accumulated by
           the internal substeps influence the first normal increment. */

        istab=0;

        printf(" Local damage substepping complete at step time %e.\n",
               theta**tper);

        /* TIME POINTS are handled inside checkconvergence().  Replacing
           its selected dtheta here could skip a requested time point, so
           only restore the pre-cutback coarse increment when no TIME
           POINTS sequence is active. */

        if((*itpamp==0)&&((1.-theta)>1.e-12)){
          dtheta=dtheta_restore;
          if(dtheta>*tmax) dtheta=*tmax;
          if(dtheta>1.-theta) dtheta=1.-theta;
          dthetaref=dtheta;
          printf(" Restoring normal increment size to %e.\n\n",
                 dtheta**tper);
        }else{
          printf(" Returning increment-size control to CalculiX.\n\n");
        }

      }else if(dtheta>dtheta_remaining){

        dtheta=dtheta_remaining;
        dthetaref=dtheta;
        printf(" Local damage substepping: next increment is clipped to "
               "%e to reach saved goal %e.\n\n",
               dtheta**tper,theta_goal**tper);
      }
    }

    /* printing the energies (only for dynamic calculations) */

    if((icutb==0)&&(*nmethod==4)&&(*ithermal<2)&&(jout[0]==jprint)&&
       (*nener==1)){

      printenergy(iexpl,ttime,&theta,tper,energy,ne,nslavs,ener,&energyref,
		  &allwk,&dampwk,&ea,&energym,&energymold,&jnz,&mscalmethod,
		  mortar,mi);

    }

    if(uncoupled){
      SFREE(iruc);
    }

    if(((qa[0]>ea*qam[0])||(qa[1]>ea*qam[1]))&&(icutb==0)){jnz++;}
    iit=0;

    if(icutb!=0){

      if(rsc.rescue_bt_on==1){
        printf("[DAMAGE RESCUE] STANDARD cutback rollback executed: vold, "
               "xbounact, f, sti, eme, ener, xstate, dam, damvisc and ipkon "
               "restored from the increment-start baselines; no additional "
               "snapshot was taken%s","\n");
        fflush(stdout);
      }

      /* A cutback may come either from the ordinary physical Newton
         solve or from the same-load damage re-equilibration.  In both
         cases checkconvergence() has already reduced dtheta.  Arm local
         substepping once and keep the target established at the start of
         the failed coarse increment.  The existing state/topology rollback
         below remains unchanged. */

      if((ilocalsubstep==0)&&(*ndmat_>0)&&(*iexpl<=1)&&
         (*nmethod!=4)&&(*idrct==0)&&
         (theta_goal>theta+1.e-12)){
        ilocalsubstep=1;
        if(damage_event_cut==1){
          printf(" Local damage substepping activated by damage-event "
                 "localization.\n");
        }else if(idamagereeq==1){
          printf(" Local damage substepping activated after damage "
                 "re-equilibration divergence.\n");
        }else{
          printf(" Local damage substepping activated after physical "
                 "Newton divergence.\n");
        }
        printf(" Saved interval: %e -> %e; reduced retry increment: %e.\n\n",
               theta_local_start**tper,theta_goal**tper,dtheta**tper);
      }

      isiz=mt**nk;cpypardou(vold,vini,&isiz,&num_cpus);

      isiz=*nboun;cpypardou(xbounact,xbounini,&isiz,&num_cpus);
      if((*ithermal==1)||(*ithermal>=3)){
	isiz=*nk;cpypardou(t1act,t1ini,&isiz,&num_cpus);
      }
      isiz=neq[1];cpypardou(f,fini,&isiz,&num_cpus);
      if(*nmethod==4){
	isiz=mt**nk;
	cpypardou(veold,veini,&isiz,&num_cpus);
	cpypardou(accold,accini,&isiz,&num_cpus);
	isiz=neq[1];
	cpypardou(fext,fextini,&isiz,&num_cpus);
	cpypardou(cv,cvini,&isiz,&num_cpus);
	if(*ithermal<2){
	  allwk=allwkini;
	  if(idamping==1)dampwk=dampwkini;
	  for(k=0;k<4;k++){
	    energy[k]=energyini[k];
	  }
	}
      }
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;
	cpypardou(sti,stiini,&isiz,&num_cpus);
	cpypardou(eme,emeini,&isiz,&num_cpus);
      }
      if(*nener==1){
	isiz=2*mi[0]*ne0;
	cpypardou(ener,enerini,&isiz,&num_cpus);
      }

      isiz=*nstate_*mi[0]*(ne0+maxprevcontel);cpypardou(xstate,xstateini,
							&isiz,&num_cpus);

      qam[0]=qamold[0];
      qam[1]=qamold[1];


      /* the rollback restores qam from before the failed attempt, which can
         be under the floor again; re-apply it so the criterion the retry
         faces is the same one the attempt faced (CCX_DAMAGE_QAM_FLOOR) */
      if((damage_qam_floor>0.)&&
         (qam[0]<damage_qam_floor*damage_cvg.qam_peak)){
        qam[0]=damage_qam_floor*damage_cvg.qam_peak;
      }

      if(*mortar>1){
	for (i=0;i<*ntie;i++){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    islavact[j]=islavactini[j];
	    bp[j]=bpini[j];
	    for(k=0;k<mt;k++){
	      cstress[mt*j+k]=cstressini[mt*j+k];
	    }
	  }    
	} 
      }
      /* if the failed attempt was a damage re-equilibration,
         restore the topology and damage state at the beginning
         of the physical increment before retrying with the smaller step */

      if(idamagereeq==1){

        if(damage_soft_reeq==1){
          printf("[DAMAGE DE1 ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e -> restore constitutive baseline\n",
                 iinc,theta**tper);
        }else if(damage_de13_transaction){
          printf("[DAMAGE DE1.3 ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e tentative_terminal=%" ITGFORMAT
                 " -> restore topology and constitutive baseline\n",
                 dtxn.increment,dtxn.step_time,
                 dtxn.count);
        }else{
          printf("[DAMAGE ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e tentative=%" ITGFORMAT "\n",
                 dtxn.increment,dtxn.step_time,
                 dtxn.count);
        }
        fflush(stdout);

        /* Trial deletions must never enter jobname.damage. */
        topo_txn_discard(&dtxn);
        if(damage_de13_transaction) lc.path_retry++;
        damage_de13_transaction=0;
        if(damage_de13_trigger_value!=NULL){
          for(i=0;i<ne0;i++){
            damage_de13_trigger_value[i]=-1.;
            damage_de13_trigger_ip[i]=0;
          }
        }

        /* A3 safety latch: if this failed transaction came from an
           aggressive fast trial, the retry uses the original exact A2
           event locator.  The latch survives the current increment retry
           and is reset only after that physical increment is accepted. */
        if(damage_fast_used==1){
          damage_fast_retry=1;
          damage_fast_used=0;
          damage_fast_recover=0;
          damage_fast_failures++;
          printf("[DAMAGE FAST BACKOFF] inc=%" ITGFORMAT
                 " failures=%" ITGFORMAT
                 " -> exact event localization on retry\n",
                 iinc,damage_fast_failures);
          fflush(stdout);
        }

	isiz=ne0;
	cpyparitg(ipkon,ipkondamageini,&isiz,&num_cpus);
	isiz=mi[0]*ne0;
	cpypardou(dam,damdamageini,&isiz,&num_cpus);

	if(*nmethod!=4){
	  isiz=mt**nk;
	  cpypardou(veold,veolddamageini,&isiz,&num_cpus);
	}

        if(damage_soft_reeq==0){
	iitsav=iit;
	iit=-2;

	remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
		   ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
		   &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
		   &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
		   &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
		   mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
		   itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
		   tieset,ntie,&num_cpus,ielmat,matname);

	iit=iitsav;

	SFREE(nactdofinv);
	NNEW(nactdofinv,ITG,mt**nk);
	MNEW(nodorig,ITG,*nk);

	FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			       ipkon,lakon,kon,ne));

	SFREE(nodorig);

	ITGMEMSET(iponoel,0,*nk,0);

	FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

	theta=thetadamage;
	idiscon=0;
	idamagereeq=0;
        damage_soft_reeq=0;

      }

    }
    /* ---- [DAMAGE CT] CLEAN PARTIAL EXIT ---------------------------
       Placed at the COMMON post-rollback join: the point reached once per
       outer pass, immediately after the complete if(icutb!=0) rollback
       block closes.  icutb!=0 here means checkconvergence() has just
       decided a cutback (it sets *icutb=0 on convergence and (*icutb)++
       otherwise), so the rollback above has just restored vold, xbounact,
       f, sti, eme, ener, xstate, dam, damvisc, ipkon and veold from the
       increment-start baselines: the reported state IS the last committed
       one, on EVERY rollback path.
       The earlier placement inside if(idamagereeq==1) was wrong twice
       over: that branch is only one of the rollback paths, so an ordinary
       cutback after a refusal never reached it, and the run ended at the
       stock stop with no PARTIAL report at all.  Placing it in the commit
       path was wrong for the same reason: a refusal never converges.
       It reports every quantity the record needs, states plainly that this
       is NOT a completed step and NOT a restart point, and leaves without
       the normal end-of-step bookkeeping.  Chosen deliberately (logged):
       the exit uses the stock 201 code, which already means "did not
       complete" and therefore cannot be mistaken for COMPLETED; a
       dedicated code would need a change outside nonlingeo.c, which this
       MVP is not allowed to make. */
    if((ct.partial==1)&&(icutb!=0)){
        ITG cq1=0,cq2=0,cq3=0,cqi,cqj,cqk;
        double cqd0,cqdl[3],cqrm[9],cqsh[3],cqop=0.;
        for(cqi=0;cqi<ne0;cqi++){
          if(ipkon[cqi]<0){cq2++;continue;}
          if(lakon[8*cqi]!='U') continue;
          if(ielprop[cqi]<0) continue;
          cqd0=prop[ielprop[cqi]+1];
          if(cqd0<=0.) continue;
          cqd0=cqd0/prop[ielprop[cqi]];
          for(cqj=0;cqj<3;cqj++){
            if(cqj>=mi[0]) break;
            cqk=mi[0]*cqi+cqj;
            if(xstate[*nstate_*cqk+3]>=0.5) cq1++;
            if(xstate[*nstate_*cqk]>cqd0) cq3++;
          }
        }
        if((ct.elem>=0)&&(ct.elem<ne0)&&
           (ipkon[ct.elem]>=0)&&(ct.w!=NULL)){
          damcont_kin(co,kon,ipkon[ct.elem],vold,mt,ct.ip,
                        cqdl,cqrm,cqsh);
          cqop=ct.m[0]*(cqdl[0]-ct.dc[0])
              +ct.m[1]*(cqdl[1]-ct.dc[1])
              +ct.m[2]*(cqdl[2]-ct.dc[2]);
        }
        printf("%s","\n");
        printf("[DAMAGE CT] ===== PARTIAL - bounded experimental "
               "continuation segment =====%s","\n");
        printf("[DAMAGE CT]   last COMMITTED state is preserved; this is NOT "
               "a completed CalculiX step and NOT a valid restart point.%s",
               "\n");
        printf("[DAMAGE CT]   theta (pseudo-time) = %.12e   time = %.12e   "
               "ttime = %.12e%s",theta,theta**tper,*ttime,"\n");
        printf("[DAMAGE CT]   lambda = %.12e   lambda_c = %.12e   "
               "ds = %.6e   kappa = %.6e%s",
               ct.lam,ct.lamc,ct.ds,ct.kappa,
               "\n");
        printf("[DAMAGE CT]   control point: element %" ITGFORMAT " ip %"
               ITGFORMAT "   m = (%.6f,%.6f,%.6f)   m.(delta-delta_c) = "
               "%.6e   |c| last = %.6e   tol_c = %.6e%s",
               ct.elem+1,ct.ip+1,ct.m[0],
               ct.m[1],ct.m[2],cqop,fabs(ct.cprev),
               ct.tolc,"\n");
        printf("[DAMAGE CT]   observables: P1 failed facets %" ITGFORMAT
               ", P2 deleted elements %" ITGFORMAT ", P3 ip with dmax>d0 %"
               ITGFORMAT "%s",cq1,cq2,cq3,"\n");
        printf("[DAMAGE CT]   cost: %" ITGFORMAT " accepted continuation "
               "commits, %" ITGFORMAT " corrector iterations, %" ITGFORMAT
               " factorisations, %" ITGFORMAT " residual evaluations, %"
               ITGFORMAT " blacklisted control point(s)%s",
               ct.ncommit,ct.ncorr,ct.nfact,
               ct.neval,ct.nbl,"\n");
        printf("[DAMAGE CT] ===== end PARTIAL report =====%s","\n");
        fflush(stdout);
        FORTRAN(stop,());
      }

    
    /* damage_event_cut is only a label for the rollback just consumed.
       Clear it before the next attempt so a later ordinary Newton
       cutback cannot be misclassified. */
    damage_event_cut=0;

    /* face-to-face penalty */

    if((*mortar==1)&&(icutb==0)&&(ncont!=0)){
	
      ntrimax=0;
      for(i=0;i<*ntie;i++){	    
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	  ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
      }
      MNEW(xo,double,ntrimax);	    
      MNEW(yo,double,ntrimax);	    
      MNEW(zo,double,ntrimax);	    
      MNEW(x,double,ntrimax);	    
      MNEW(y,double,ntrimax);	    
      MNEW(z,double,ntrimax);	   
      MNEW(nx,ITG,ntrimax);	   
      MNEW(ny,ITG,ntrimax);	    
      MNEW(nz,ITG,ntrimax);
      
      /*  Determination of active nodes (islavact) */
      
      FORTRAN(islavactive,(tieset,ntie,itietri,cg,straight,
			   co,vold,xo,yo,zo,x,y,z,nx,ny,nz,mi,
			   imastop,nslavnode,islavnode,islavact));

      SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
      SFREE(ny);SFREE(nz);

      if(*ithermal!=2){
	if(negpres==0){
	  if((*mortar==1)&&(1.-theta-dtheta<=1.e-6)&&(itruecontact==1)){
	    printf(" pressure ratio (smallest/largest pressure over all contact areas) =%e\n\n",pressureratio);
	    	    if(pressureratio<-0.05){
	    //	    if((pressureratio<-0.05)||((*nmethod==1)&&(iperturb[1]==1))){
	      printf(" zero-size increment is appended\n\n");
	      negpres=1;theta=1.-1.e-6;dtheta=1.e-6;
	    }
	  }
	}else{negpres=0;}
      }

    }

    /* output */

    if((jout[0]==jprint)&&(icutb==0)){

      jprint=0;

      /* calculating the displacements and the stresses and storing */
      /* the results in frd format  */
	
      MNEW(v,double,mt**nk);
      MNEW(fn,double,mt**nk);
      NNEW(stn,double,6**nk);
      if(*ithermal>1) NNEW(qfn,double,3**nk);
      NNEW(inum,ITG,*nk);
      NNEW(stx,double,6*mi[0]**ne);
      
      if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
      if(strcmp1(&filab[435],"PEEQ")==0) NNEW(epn,double,*nk);
      if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
      if(strcmp1(&filab[609],"SDV ")==0) NNEW(xstaten,double,*nstate_**nk);
      if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);
      if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
      if(strcmp1(&filab[4785],"DUCT")==0) NNEW(damn,double,*nk);

      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);

      if((*mortar==-1)&&(idispfrdonly==1)){

	/* nothing to do if massless explicit dynamics and all output
           consists of displacements */
	
	cpyparitg(inum,inumcp,nk,&num_cpus);

      }else{
      
	iout=2;
	icmd=3;
      
#ifdef COMPANY
	FORTRAN(uinit,());
#endif
	trial_results(&nlgt);
      
	isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);

      }
      
      iout=0;
      if(*iexpl<=1) icmd=0;
      
      ++*kode;
      if(*mcs!=0){
	ptime=*ttime+time;

	if(*mortar>1){
	  mortar_prefrd(ne,nslavs,mi,nk,nkon,&stx,cdisp,fn,cfs,cfm);       
	}
	
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
               nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
               ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
               cdn,mortar,nmat,qfx,ielprop,prop,damn,&errn);

	if(*mortar>1){
	  mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
	}
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
      }
      else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	  NNEW(neigh,ITG,40**ne);
	  MNEW(ipneigh,ITG,*nk);
	}

	ptime=*ttime+time;

	if(*mortar>1){
	  mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
	}
	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	    prop,sti,damn,&errn);
	if(*mortar>1){
	  mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
	}

	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
      }

      /* mesh refinement */
  
      if(strcmp1(&filab[4089],"RM")==0){
	refinemesh(nk,ne,co,ipkon,kon,v,veold,stn,een,emn,epn,enern,
		   qfn,errn,filab,mi,lakon,jobnamec,istartset,iendset,
		   ialset,set,nset,matname,ithermal,output,nmat,
		   nelemload,nload,sideload,nodeforc,
		   nforc,nodeboun,nboun,nodempc,ipompc,nmpc);

	/* free errn */
	
	if(((*nmethod!=5)||(mode==-1))&&
	   ((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2))) SFREE(errn);
      }
      
      SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);
      if(*ithermal>1){SFREE(qfn);}
      
      if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
      if(strcmp1(&filab[435],"PEEQ")==0) SFREE(epn);
      if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
      if(strcmp1(&filab[609],"SDV ")==0) SFREE(xstaten);
      if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);
      if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
      if(strcmp1(&filab[4785],"DUCT")==0) SFREE(damn);
    }
    
  }

  /*********************************************************/
  /*   end of the increment loop                          */
  /*********************************************************/

  if(jprint!=0){
    
    /* printing the energies (only for dynamic calculations) */

    if((*nmethod==4)&&(*ithermal<2)&&(*nener==1)){

      printenergy(iexpl,ttime,&theta,tper,energy,ne,nslavs,ener,&energyref,
		  &allwk,&dampwk,&ea,&energym,&energymold,&jnz,&mscalmethod,
		  mortar,mi);
    }

    /* calculating the displacements and the stresses and storing  
       the results in frd format */
  
    MNEW(v,double,mt**nk);
    MNEW(fn,double,mt**nk);
    NNEW(stn,double,6**nk);
    if(*ithermal>1) NNEW(qfn,double,3**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
    if(strcmp1(&filab[435],"PEEQ")==0) NNEW(epn,double,*nk);
    if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
    if(strcmp1(&filab[609],"SDV ")==0) NNEW(xstaten,double,*nstate_**nk);
    if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
    if(strcmp1(&filab[4785],"DUCT")==0) NNEW(damn,double,*nk);
    
    isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
    iout=2;
    icmd=3;

#ifdef COMPANY
    FORTRAN(uinit,());
#endif
    trial_results(&nlgt);
    
    isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);

    iout=0;
    if(*iexpl<=1) icmd=0;
    
    ++*kode;
    if(*mcs>0){
      ptime=*ttime+time;
      if(*mortar>1){
	mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
      }
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	     t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
             nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
             ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
             cdn,mortar,nmat,qfx,ielprop,prop,damn,&errn);
      if(*mortar>1){
	mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
      }
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif

    }else{
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	MNEW(ipneigh,ITG,*nk);
      }

      ptime=*ttime+time;
      if(*mortar>1){
	mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
      }
      frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	  prop,sti,damn,&errn);
      if(*mortar>1){
	mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
      }

      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
    }

    /* mesh refinement */
  
    if(strcmp1(&filab[4089],"RM")==0){
      refinemesh(nk,ne,co,ipkon,kon,v,veold,stn,een,emn,epn,enern,
		 qfn,errn,filab,mi,lakon,jobnamec,istartset,iendset,
		 ialset,set,nset,matname,ithermal,output,nmat,
		 nelemload,nload,sideload,nodeforc,
		 nforc,nodeboun,nboun,nodempc,ipompc,nmpc);

      /* free errn */
	
      if(((*nmethod!=5)||(mode==-1))&&
	 ((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2))) SFREE(errn);
    }

    SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);
    if(*ithermal>1){SFREE(qfn);}
    
    if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
    if(strcmp1(&filab[435],"PEEQ")==0) SFREE(epn);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    if(strcmp1(&filab[609],"SDV ")==0) SFREE(xstaten);
    if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
    if(strcmp1(&filab[4785],"DUCT")==0) SFREE(damn);

  }
    
  /* writing out the latest stiffness matrix for a subsequent
     sensitivity analysis */

  if(isensitivity){
      
    strcpy2(stiffmatrix,jobnamec,132);
    strcat(stiffmatrix,".stm");
      
    if((f1=fopen(stiffmatrix,"wb"))==NULL){
      printf(" *ERROR in nonlingeo: cannot open stiffness matrix file for writing...");
      exit(0);
    }
      
    /* storing the stiffness matrix */

    /* nzs,irow,jq and icol have to be stored too, since the static analysis
       can involve contact, whereas in the sensitivity analysis contact is not
       taken into account while determining the structure of the stiffness
       matrix (in mastruct.c)
    */
      
    if(fwrite(&nasym,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in nonlingeo saving the symmetry flag to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(nzs,sizeof(ITG),3,f1)!=3){
      printf(" *ERROR in nonlingeo saving the number of subdiagonal nonzeros to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(irow,sizeof(ITG),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR in nonlingeo saving irow to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(jq,sizeof(ITG),neq[1]+1,f1)!=neq[1]+1){
      printf(" *ERROR in nonlingeo saving jq to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(icol,sizeof(ITG),neq[1],f1)!=neq[1]){
      printf(" *ERROR in nonlingeo saving icol to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(adcpy,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in nonlingeo saving the diagonal of the stiffness matrix to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(aucpy,sizeof(double),(nasym+1)*nzs[2],f1)!=(nasym+1)*nzs[2]){
      printf(" *ERROR in nonlingeo saving the off-diagonal terms of the stiffness matrix to the stiffness matrix file...");
      exit(0);
    }
    fclose(f1);
    SFREE(adcpy);SFREE(aucpy);
  }
  
  /* restoring the distributed loading  */

  if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
    *nload=nloadref;
    RENEW(nelemload,ITG,2**nload);
    isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
    if(*nam>0){
      RENEW(iamload,ITG,2**nload);
      isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
    }
    RENEW(sideload,char,20**nload);memcpy(&sideload[0],&sideloadref[0],
					  sizeof(char)*20**nload);
      
    /* freeing the temporary fields */
      
    SFREE(nelemloadref);if(*nam>0){SFREE(iamloadref);};
    SFREE(sideloadref);
  }

  /* setting the velocity to zero at the end of a quasistatic or stationary
     step */

  if(abs(*nmethod)==1){
    for(k=0;k<mt**nk;++k){
      veold[k]=0.;}
  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){

    /* thermal boundary conditions are updated only if the
       step was thermal or thermomechanical */

    if(ndirboun[k]==0){
      if(*ithermal<2) continue;

      /* mechanical boundary conditions are updated only
	 if the step was not thermal or the node is a
	 network node */

    }else if((ndirboun[k]>0)&&(ndirboun[k]<4)){
      node=nodeboun[k];
      FORTRAN(nident,(itg,&node,&ntg,&id));
      networknode=0;
      if(id>0){
	if(itg[id-1]==node) networknode=1;
      }
      if((*ithermal==2)&&(networknode==0)) continue;
    }
    xbounold[k]=xbounact[k];
  }
  isiz=*nforc;cpypardou(xforcold,xforcact,&isiz,&num_cpus);
  isiz=2**nload;cpypardou(xloadold,xloadact,&isiz,&num_cpus);
  isiz=7**nbody;cpypardou(xbodyold,xbodyact,&isiz,&num_cpus);
  if(*ithermal==1){
    cpypardou(t1old,t1act,nk,&num_cpus);
    for(k=0;k<*nk;++k){
      vold[mt*k]=t1act[k];}
  }
  else if(*ithermal>1){
    for(k=0;k<*nk;++k){
      t1[k]=vold[mt*k];}
    if(*ithermal>=3){
      cpypardou(t1old,t1act,nk,&num_cpus);
    }
  }

  qaold[0]=qa[0];
  qaold[1]=qa[1];
  
  if(*iexpl>1){
    SFREE(smscale);

    if((mscalmethod==1)||(mscalmethod==3)||(*mortar==-1)){
      if(*isolver==0){
#ifdef SPOOLES
	spooles_cleanup();
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	sgi_cleanup(token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_cleanup(&neq[0],&symmetryflag,&inputformat);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
#endif
      }
    }
  }
  
  SFREE(f);SFREE(b);
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(xbodyact);

  if(*inewton==1){SFREE(cgr);}
  SFREE(fext);SFREE(ampli);SFREE(xbounini);SFREE(xstiff);
  if((*ithermal==1)||(*ithermal>=3)){SFREE(t1act);SFREE(t1ini);}

  if(*ithermal>1){
    SFREE(itg);SFREE(ieg);SFREE(kontri);SFREE(nloadtr);
    SFREE(nactdog);SFREE(nacteq);SFREE(ineighe);
    SFREE(tarea);SFREE(tenv);SFREE(fenv);SFREE(qfx);
    SFREE(erad);SFREE(ac);SFREE(bc);SFREE(ipiv);
    SFREE(bcr);SFREE(ipivr);SFREE(adview);SFREE(auview);SFREE(adrad);
    SFREE(aurad);SFREE(irowrad);SFREE(jqrad);SFREE(icolrad);
    if((*mcs>0)&&(ntr>0)){SFREE(inocs);}
    if((*network>0)||(ntg>0)){SFREE(iponoeln);SFREE(inoeln);}
    if(ntr>0){
    }
  }

  if(icfd==1){
  }else if(icfd==2){
    SFREE(sideface);SFREE(nelemface);SFREE(ifreestream);
    SFREE(isolidsurf);SFREE(neighsolidsurf);SFREE(iponoelf);SFREE(inoelf);
    SFREE(inomat);SFREE(ipface);
    if(*ithermal==1) SFREE(qfx);
	 
    SFREE(ipkonf);SFREE(lakonf);SFREE(ielmatf);SFREE(nelold);SFREE(nelnew);
    SFREE(cof);SFREE(voldf);SFREE(nkold);SFREE(nknew);SFREE(konf);
    SFREE(ipompcf);SFREE(nodempcf);SFREE(coefmpcf);SFREE(nodebounf);
    SFREE(ndirbounf);SFREE(xbounf);SFREE(nelemloadf);SFREE(xloadf);
    SFREE(sideloadf);SFREE(ikbounf);SFREE(ilbounf);
    SFREE(ikmpcf);SFREE(ilmpcf);SFREE(xbounoldf);SFREE(xbounactf);
    SFREE(xloadoldf);SFREE(xloadactf);SFREE(inotrf);
    if(*norien>0) SFREE(ielorienf);
    if(*nbody>0) SFREE(ipobodyf);
    if(*nam>0){SFREE(iambounf);SFREE(iamloadf);}
  }

  SFREE(fini);
  if(*nmethod==4){
    SFREE(aux2);SFREE(fextini);SFREE(veini);SFREE(accini);
    SFREE(adb);SFREE(aub);SFREE(cvini);SFREE(cv);SFREE(fnext);
    SFREE(fnextini);
  }
  SFREE(eei);SFREE(stiini);SFREE(emeini);
  if(*nener==1)SFREE(enerini);
  if(*nstate_!=0){SFREE(xstateini);}

  SFREE(aux);SFREE(iaux);SFREE(vini);

  if((*ndmat_>0)&&(*iexpl<=1)){
    /* The results context owns no storage.  Clear it before releasing the
       nonlingeo damage baselines so later steps/output calls cannot retain
       dangling pointers. */
    results_set_de12_context(0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0.);
    /* [GLOBALIZE] close the last attempt and report.  Unconditional:
       a census that has to be switched on is a census nobody reads,
       and this one exists to be read before a mechanism is deleted. */
    glob_attempt_end(&damage_glob);
    glob_census_report(&damage_glob);

    topo_txn_discard(&dtxn);

    if(fdamage!=NULL){
      fflush(fdamage);
      fclose(fdamage);
      fdamage=NULL;
    }

    SFREE(ipkondamageini);
    SFREE(damdamageini);
    SFREE(damde1prev);
    if(damage_damjac!=NULL) SFREE(damage_damjac);
    if(damage_fracture_seta!=NULL) free(damage_fracture_seta);
    if(lc.diss_fhat!=NULL) SFREE(lc.diss_fhat);
    if(lc.diss_uf!=NULL) SFREE(lc.diss_uf);
    if(damage_damvisc!=NULL) SFREE(damage_damvisc);
    if(damage_damviscini!=NULL) SFREE(damage_damviscini);
    if(rsc.bt_dam!=NULL) SFREE(rsc.bt_dam);
    if(rsc.bt_visc!=NULL) SFREE(rsc.bt_visc);
    if(rsc.bt_xs!=NULL) SFREE(rsc.bt_xs);
    if(prb.evt_sgn!=NULL) SFREE(prb.evt_sgn);
    if(prb.ray_r0!=NULL) SFREE(prb.ray_r0);
    if(prb.ray_cat!=NULL) SFREE(prb.ray_cat);
    if(prb.ray_p!=NULL) SFREE(prb.ray_p);
    if(prb.ray_res!=NULL) SFREE(prb.ray_res);
    if(ct.sgn!=NULL) SFREE(ct.sgn);
    if(ct.ring!=NULL) SFREE(ct.ring);
    if(ct.fl!=NULL) SFREE(ct.fl);
    if(ct.bl!=NULL) SFREE(ct.bl);
    if(ct.w!=NULL) SFREE(ct.w);
    if(ct.r0!=NULL) SFREE(ct.r0);
    if(ct.beps!=NULL) SFREE(ct.beps);
    if(ct.y!=NULL) SFREE(ct.y);
    if(ct.z!=NULL) SFREE(ct.z);
    if(ct.qh!=NULL) SFREE(ct.qh);
    if(ct.dam!=NULL) SFREE(ct.dam);
    if(ct.visc!=NULL) SFREE(ct.visc);
    if(ct.xs!=NULL) SFREE(ct.xs);
    if(ct.jac!=NULL) SFREE(ct.jac);
    if(dog.r0!=NULL) SFREE(dog.r0);
    if(dog.d!=NULL) SFREE(dog.d);
    if(dog.w!=NULL) SFREE(dog.w);
    if(dog.pn!=NULL) SFREE(dog.pn);
    if(dog.res!=NULL) SFREE(dog.res);
    if(dog.dam!=NULL) SFREE(dog.dam);
    if(dog.visc!=NULL) SFREE(dog.visc);
    if(dog.xs!=NULL) SFREE(dog.xs);
    if(dog.mode==1){
      printf("[DAMAGE TR] SUMMARY: armed %" ITGFORMAT " time(s); accepted "
             "steps %" ITGFORMAT " (Newton %" ITGFORMAT ", Cauchy %" ITGFORMAT
             ", dogleg %" ITGFORMAT "); iterations that accepted nothing %"
             ITGFORMAT "; rejected trials %" ITGFORMAT "; residual "
             "evaluations %" ITGFORMAT "; armed factorisations %" ITGFORMAT
             "; last transpose check %.12e; operator asymmetry %.6e\n",
             dog.narm,dog.nacc,dog.nnewt,dog.ncau,
             dog.ndog,dog.nfail,dog.nrej,dog.neval,
             dog.nfact,dog.ident,dog.asym);
      fflush(stdout);
    }
    if(damage_frel!=NULL) SFREE(damage_frel);
    if(damage_ract!=NULL) SFREE(damage_ract);
    if(damage_de13_trigger_value!=NULL) SFREE(damage_de13_trigger_value);
    if(damage_de13_trigger_ip!=NULL) SFREE(damage_de13_trigger_ip);
    if(*nmethod!=4) SFREE(veolddamageini);
  }

  if(icascade==2){
    memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
    RENEW(nodempc,ITG,3*memmpcref_);
    for(k=0;k<3*memmpcref_;k++){
      nodempc[k]=nodempcref[k];}
    RENEW(coefmpc,double,memmpcref_);
    for(k=0;k<memmpcref_;k++){
      coefmpc[k]=coefmpcref[k];}
    SFREE(nodempcref);SFREE(coefmpcref);
  }

  if(ncont!=0){
    *ne=ne0;*nkon=nkon0;
    if((*nener==1)&&(*mortar==1)){
      RENEW(ener,double,mi[0]**ne*2);
    }
    RENEW(ipkon,ITG,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(kon,ITG,*nkon);
    if(*norien>0){
      RENEW(ielorien,ITG,mi[2]**ne);
    }
    RENEW(ielmat,ITG,mi[2]**ne);

    if(*mortar>1){
      
      /// needed for next step coloumb friction
      
      for (i=0;i<*ntie;i++){
	if(tieset[i*(81*3)+80]=='C'){
	  if(*nstate_*mi[0]>0){
	    for(j=nslavnode[i];j<nslavnode[i+1];j++){	  	     
	      for(k=0;k<3;k++){	    	       
		xstate[*nstate_*mi[0]*(*ne+j)+k]=cstress[mt*j+k]; 
	      } 	    	       
	      xstate[*nstate_*mi[0]*(*ne+j)+3]=islavact[j]+0.5;
	    }
		      
	  }
	}
      }
    }
      
    SFREE(cg);SFREE(straight);
    SFREE(imastop);SFREE(itiefac);SFREE(islavnode);
    SFREE(nslavnode);SFREE(iponoels);SFREE(inoels);SFREE(imastnode);
    SFREE(nmastnode);SFREE(itietri);SFREE(koncont);SFREE(xnoels);
    SFREE(springarea);SFREE(xmastnor);

    if(*mortar==-1){
      if(ncont!=0){SFREE(kslav);SFREE(lslav);SFREE(ktot);SFREE(ltot);
	SFREE(aloc);SFREE(alglob);SFREE(areaslav);SFREE(fric);}
      SFREE(adc);SFREE(auc);
      if(idispfrdonly==1){SFREE(inumcp);}
      if(masslesslinear>0){
	SFREE(ad);SFREE(au);
	if(ncont!=0){SFREE(auw);SFREE(jqw);SFREE(iroww);
	  SFREE(fullgmatrix);SFREE(fullr);}
	iclean=1;
        massless(kslav,lslav,ktot,ltot,au,ad,auc,adc,jq,irow,neq,nzs,auw,jqw,
		 iroww,&nzsw,islavnode,nslavnode,nslavs,imastnode,nmastnode,
		 ntie,nactdof,mi,vold,volddof,veold,nk,fext,isolver,
		 &masslesslinear,co,springarea,&neqtot,qb,b,&dtime,aloc,fric,
		 iexpl,nener,ener,ne,&jqbi,&aubi,&irowbi,&jqib,&auib,&irowib,
		 &iclean,&iinc,fullgmatrix,fullr,alglob,&num_cpus,&ncont);
      }
      if(masslesslinear==2){SFREE(fextload);}

    }else if(*mortar==0){
      SFREE(areaslav);
    }else if(*mortar==1){
      SFREE(pmastsurf);SFREE(ipe);SFREE(ime);
      SFREE(islavact);
    }else if(*mortar>1){
      SFREE(islavact);SFREE(gap);SFREE(slavnor);SFREE(slavtan);
      SFREE(cstress);SFREE(ipe);SFREE(ime);SFREE(cfs);SFREE(cfm);
      SFREE(cdisp);SFREE(bp);SFREE(islavtie);
      SFREE(nslavspc);SFREE(islavspc);SFREE(nslavmpc);SFREE(islavmpc);
      SFREE(nmastspc);SFREE(imastspc);SFREE(nmastmpc);SFREE(imastmpc);
      SFREE(pslavdual);
      SFREE(cstressini);SFREE(bpini);SFREE(islavactini);
      SFREE(aut);SFREE(irowt);SFREE(jqt);
      SFREE(autinv);SFREE(irowtinv);SFREE(jqtinv);
      SFREE(Bd);SFREE(irowb);SFREE(jqb);
      SFREE(Bdhelp);SFREE(irowbhelp);SFREE(jqbhelp);
      SFREE(Dd);SFREE(irowd);SFREE(jqd);
      SFREE(Ddtil);SFREE(irowdtil);SFREE(jqdtil);
      SFREE(Bdtil);SFREE(irowbtil);SFREE(jqbtil);
      SFREE(islavnodeinv);SFREE(islavquadel);
    }
  }

  /* reset icascade */

  if(icascade==1){icascade=0;}

  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;

  if(iglob==1){SFREE(integerglob);SFREE(doubleglob);}

  *icolp=icol;*irowp=irow;*cop=co;*voldp=vold;

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*nelemloadp=nelemload;
  *iamloadp=iamload;*sideloadp=sideload;

  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;

  *islavsurfp=islavsurf;*pslavsurfp=pslavsurf;*clearinip=clearini;

  (*tmin)*=(*tper);
  (*tmax)*=(*tper);

  SFREE(nactdofinv);
  // MPADD start
  if((*nmethod==4)&&(*ithermal!=2)&&(*iexpl<=1)&&(icfd==0)){ SFREE(adblump);}
  // MPADD end
  
  (*ttime)+=(*tper);

  SFREE(iponoel);
  
  return;
}
