/*     CalculiX - A 3-dimensional finite element program                   */
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

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"
#include "pardiso.h"

/* next line is for the simultaneous use of PARDISO and PaStiX */

#include <mkl_service.h>

ITG *icolpardiso=NULL,*pointers=NULL,iparm[64];

/* The CSR REPACK cache.
   ------------------------------------------------------------------
   The mtype=1 branch below - structurally symmetric, numerically
   asymmetric, which is what CCX_DAMAGE_TANGENT=UNSYM produces and
   therefore what every run of this branch takes - rebuilds its whole CSR
   on EVERY factorisation: four allocations, a full sort of the lower
   triangle by row (isortiid over nzs entries), a sort of every row by
   column (isortid), and an interleave into neq+2*nzs slots.

   Measured on s3rad: 73 ms a call, 5344 calls, 9.4% of a 69-minute run -
   6.5 minutes spent rewriting the same numbers into the same places
   (research/01-PROFILING.md).

   But the permutation that sort produces depends ONLY on (icol,irow,jq),
   the sparsity pattern - the values are pure payload.  The mechanism that
   decides whether the pattern changed already exists and already runs:
   pardiso_structure_hash, which says it changes on 165 calls out of 5344.
   So on 97% of calls the sort reproduces itself exactly.

   pardiso_src[k] records where slot k came from:
       >= 0  ->  au[ pardiso_src[k] ]
       <  0  ->  ad[ -pardiso_src[k]-1 ]
   and refilling is then one indexed scatter.  It is built by running the
   ordinary construction with the SOURCE INDEX as payload instead of the
   value - the sorts carry a double and an index below 2^53 is exact in
   one - so there is exactly one copy of the ordering logic and the fast
   path cannot drift from it. */

static ITG *pardiso_src=NULL;
static ITG pardiso_src_n=0;
static ITG pardiso_repack_hits=0,pardiso_repack_builds=0;

void pardiso_repack_report(void);

static ITG pardiso_repack_on(void)
{
  static ITG init=0,on=0;
  if(!init){
    /* ON by default since 2026-09-12.  It only ever engages where
       symbolic reuse is already eligible, so it rides an existing opt-in
       and changes nothing where that is off.  The evidence for the flip:
       bit-identical on both fast decks and on the target deck - the
       deletion record byte-for-byte, all 3741, same ending increment -
       for 11.7x less repacking and 14% less wall clock at scale
       (research/01-PROFILING.md).  Set it to 0 to turn it off. */

    const char *e=ccxopt_getenv("CCX_PARDISO_REPACK");
    init=1;
    atexit(pardiso_repack_report);
    on=1;
    if((e!=NULL)&&((e[0]=='0')||(e[0]=='N')||(e[0]=='n')||
                   (e[0]=='F')||(e[0]=='f')||(e[0]==0))) on=0;
  }
  return on;
}

static ITG pardiso_repack_verify(void)
{
  static ITG init=0,on=0;
  if(!init){ init=1; on=(ccxopt_getenv("CCX_PARDISO_REPACK_VERIFY")!=NULL)?1:0; }
  return on;
}

static void pardiso_repack_free(void)
{
  if(pardiso_src!=NULL){ SFREE(pardiso_src); pardiso_src=NULL; }
  pardiso_src_n=0;
}

void pardiso_repack_report(void)
{
  if((pardiso_repack_builds+pardiso_repack_hits)<=0) return;
  printf("[PARDISO REPACK] rebuilt %" ITGFORMAT " time(s), refilled from the "
         "cached permutation %" ITGFORMAT " time(s)\n",
         pardiso_repack_builds,pardiso_repack_hits);
  fflush(stdout);
}
long long pt[64];
double *aupardiso=NULL;
/* double dparm[64];  not used */
ITG mthread_mkl=0;
char envMKL[32];

/* Optional symbolic-analysis cache for the symmetric nonlinear path.
   CalculiX normally calls PARDISO phase 12 and phase -1 for every Newton
   iteration even when the equation graph is unchanged.  With the opt-in
   environment setting below, retain PARDISO's symbolic state and execute
   phase 22 for subsequent numerical factorizations.  A hash of the actual
   sparse graph invalidates the cache after remastruct or any other structure
   change.  Unsymmetric storage remains on the untouched stock path. */
static ITG pardiso_reuse_mode=-1,pardiso_cache_valid=0,
  pardiso_cache_neq=0,pardiso_cache_nzs=0,
  pardiso_cache_symmetry=0,pardiso_cache_inputformat=0;
static unsigned long long pardiso_cache_hash=0;

static ITG pardiso_symbolic_reuse_requested(void)
{
  char *env;

  if(pardiso_reuse_mode>=0) return pardiso_reuse_mode;
  pardiso_reuse_mode=0;
  env=ccxopt_getenv("CCX_PARDISO_REUSE_SYMBOLIC");
  if((env!=NULL)&&
     ((strcmp(env,"1")==0)||(strcmp(env,"ON")==0)||
      (strcmp(env,"on")==0)||(strcmp(env,"YES")==0)||
      (strcmp(env,"yes")==0))) pardiso_reuse_mode=1;
  return pardiso_reuse_mode;
}

/* CCX_PARDISO_CGS=<L>: solve with preconditioned CGS using the LU that is
   already in the cache, instead of refactorising.  MKL does this inside a
   combined phase=23 call and falls back to a full numerical factorisation by
   itself when the Krylov iteration fails, so the downside is bounded by one
   wasted attempt.  L is the stopping tolerance exponent, 10^-L.

   E-101 removed the gate that kept the SYMBOLIC factorisation from being
   reused; on 84 000 elements that was worth only 2.3%, because the analysis
   is a nearly fixed cost while the NUMERICAL factorisation grows as ~N^1.5.
   This attacks the numerical phase instead, which is the part that dominates
   at the sizes the true-scale models need.  Measured justification: Newton
   here contracts the residual by a median 0.52-0.60 per iteration with only
   8-16% of steps quadratic, so the exact tangent is not buying quadratic
   convergence and a slightly stale operator should cost little. */

static ITG pardiso_cgs_mode=-1;
static double *pardiso_cgs_rhs=NULL;
static ITG pardiso_cgs_nrhs=1,pardiso_cgs_done=0;

/* Iterative refinement, counted.

   iparm(7) - iparm[6] here - is the number of refinement steps PARDISO
   ACTUALLY performed for a solve, as opposed to iparm(8) which only bounds
   them.  Nothing in this tree read it, and the complete s3rad profile
   (research/01-PROFILING.md) turned that into a question with no answer: the
   triangular solve steps from 28 ms a call to 87 ms across two 120-second
   windows while the numeric factorisation that produced those factors stays
   flat at 377 ms, the system size moves 0.3%, solves per factorisation stay
   at exactly 1.00, and the CGS path is not armed.  Refinement is the only
   candidate left that makes a back substitution three times more expensive
   without touching the factorisation.

   The falsifiable form: if this is refinement, the fast decks - which have
   no notch process zone holding thousands of elements at the residual
   stiffness floor - report a mean near zero, and s3rad reports a mean near
   two with the rise landing at increments 195-222.  If both report zero, the
   hypothesis is dead and the cost is somewhere else. */

static double pardiso_ref_solves=0.,pardiso_ref_steps=0.;
static ITG pardiso_ref_max=0,pardiso_ref_atexit=0;
static ITG pardiso_ref_every=200;

static void pardiso_refine_report(void){
  if(pardiso_ref_solves<=0.) return;
  printf("[PARDISO REFINE] solves=%.0f steps=%.0f mean=%.3f max=%" ITGFORMAT
         "\n",pardiso_ref_solves,pardiso_ref_steps,
         pardiso_ref_steps/pardiso_ref_solves,pardiso_ref_max);
  fflush(stdout);
}

static void pardiso_refine_note(ITG steps){
  if(!pardiso_ref_atexit){
    const char *e=ccxopt_getenv("CCX_PARDISO_REFINE_EVERY");
    if(e!=NULL){
      ITG v=(ITG)atoi(e);
      if(v>0) pardiso_ref_every=v;
    }
    pardiso_ref_atexit=1;
    atexit(pardiso_refine_report);
  }
  pardiso_ref_solves+=1.;
  pardiso_ref_steps+=(double)steps;
  if(steps>pardiso_ref_max) pardiso_ref_max=steps;
  if((pardiso_ref_every>0)&&
     (((ITG)pardiso_ref_solves)%pardiso_ref_every==0)) pardiso_refine_report();
}
static ITG pardiso_cgs_ok=0,pardiso_cgs_fail=0,pardiso_cgs_iter=0;

static ITG pardiso_cgs_level(void)
{
  char *env;
  ITG v;

  if(pardiso_cgs_mode>=0) return pardiso_cgs_mode;
  pardiso_cgs_mode=0;
  env=ccxopt_getenv("CCX_PARDISO_CGS");
  if(env!=NULL){
    v=atoi(env);
    if((v>0)&&(v<10)) pardiso_cgs_mode=v;
  }
  return pardiso_cgs_mode;
}

static void pardiso_cgs_report(ITG force)
{
  ITG n=pardiso_cgs_ok+pardiso_cgs_fail;

  if(n<=0) return;
  if((!force)&&(n%100!=0)) return;
  printf("[PARDISO CGS] %" ITGFORMAT " of %" ITGFORMAT " solves reused the LU "
         "(%.1f%%), mean %.1f CGS iterations; %" ITGFORMAT " fell back to a "
         "full factorisation\n",
         pardiso_cgs_ok,n,100.0*pardiso_cgs_ok/n,
         pardiso_cgs_ok>0?(double)pardiso_cgs_iter/pardiso_cgs_ok:0.0,
         pardiso_cgs_fail);
  fflush(stdout);
}

static ITG pardiso_reuse_eligible(ITG symmetryflag,ITG inputformat)
{
  /* Which matrix types may reuse the cached symbolic factorisation.

     The cache is keyed on neq, nzs and a hash of (icol,irow) - the
     lower-triangular pattern.  That describes the structure completely for
     the symmetric type (mtype=-2) and for the structurally symmetric,
     numerically asymmetric type (mtype=1).  mtype=1 is what the damage
     rank-1 tangent produces: mafilldamas writes into the pattern mafillsm
     has already built and the upper half mirrors it, so losing symmetry
     does not change the sparsity.  Until 2026-08-28 the gate here was
     (symmetryflag==0), so every run of the fracture branch - which is
     always CCX_DAMAGE_TANGENT=UNSYM, symmetryflag=2 - recomputed the
     ordering for every factorisation.

     inputformat==3 (mtype=11, structurally asymmetric - the contact path)
     builds its pattern from jq and nzs3, which the hash does NOT cover, so
     it stays excluded. */

  if(!pardiso_symbolic_reuse_requested()) return 0;
  if(symmetryflag==0) return 1;
  return (inputformat!=3);
}

static unsigned long long pardiso_structure_hash(const ITG *icol,
                                                 const ITG *irow,
                                                 ITG neq,ITG nzs)
{
  ITG i;
  unsigned long long h=1469598103934665603ULL;

  h^=(unsigned long long)neq;h*=1099511628211ULL;
  h^=(unsigned long long)nzs;h*=1099511628211ULL;
  for(i=0;i<neq;i++){
    h^=(unsigned long long)(unsigned int)icol[i];
    h*=1099511628211ULL;
  }
  for(i=0;i<nzs;i++){
    h^=(unsigned long long)(unsigned int)irow[i];
    h*=1099511628211ULL;
  }
  return h;
}

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
		    double *sigma,ITG *icol, ITG *irow, 
		    ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		    ITG *jq, ITG *nzs3){
  logview_begin_named("pardiso factor");

  char *env;
  /*  char env1[32]; */
  ITG i,j,k,l,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
    msglvl=0,error=0,*irowpardiso=NULL,kflag,kstart,n,ifortran,
    lfortran,index,id,k2,reuse_requested=0,reuse_current=0,cgs_active=0,
    repack_on=0,repack_fast=0,repack_verify=0;
  ITG ndim,nthread,nthread_v;
  double *b=NULL,*x=NULL,*repack_ref=NULL;
  unsigned long long structure_hash=0;

  repack_on=pardiso_repack_on();
  repack_verify=pardiso_repack_verify();
  reuse_requested=pardiso_reuse_eligible(*symmetryflag,*inputformat);
  if(reuse_requested){
    /* The price of the reuse mechanism, measured separately from what it
       buys.  This is an O(neq+nzs) pass taken on EVERY factorisation in
       order to decide whether to skip an analysis that is needed on about
       one call in seventy.  At s3rad scale the whole of pardiso_factor
       outside the solver call costs 82 ms against a 764 ms factorisation -
       6.3% of the run - and this event says how much of that is the hash. */
    logview_begin_named("pardiso structure hash");
    structure_hash=pardiso_structure_hash(icol,irow,*neq,*nzs);
    logview_end_named("pardiso structure hash");
    if((pardiso_cache_valid)&&
       (pardiso_cache_neq==*neq)&&(pardiso_cache_nzs==*nzs)&&
       (pardiso_cache_symmetry==*symmetryflag)&&
       (pardiso_cache_inputformat==*inputformat)&&
       (pardiso_cache_hash==structure_hash)) reuse_current=1;
  }

  if((pardiso_cache_valid)&&(!reuse_current)){
    ITG cache_neq=pardiso_cache_neq;
    ITG cache_symmetry=pardiso_cache_symmetry;
    ITG cache_inputformat=pardiso_cache_inputformat;
    pardiso_cleanup(&cache_neq,&cache_symmetry,&cache_inputformat);
  }

  if(reuse_current) phase=22;

  /* A cached LU plus a RHS handed in by pardiso_main is everything the
     combined factorise-and-solve phase needs. */
  if((pardiso_cgs_level()>0)&&(reuse_current)&&(pardiso_cgs_rhs!=NULL)){
    cgs_active=1;
    phase=23;
  }

  if(*symmetryflag==0){
    printf(" Factoring the system of equations using the symmetric pardiso solver\n");
  }else{
    printf(" Factoring the system of equations using the unsymmetric pardiso solver\n");
  }

  iparm[0]=0;
  iparm[1]=3;
  iparm[3]=0;
  if(cgs_active){
    /* iparm[0]=0 tells PARDISO to OVERWRITE iparm[1..63] with its defaults,
       so anything set here would be lost - which is also why the iparm[1]=3
       above has never taken effect.  After the first full call the array
       already holds those defaults, so switching to iparm[0]=1 preserves
       them and lets exactly one entry be changed. */
    iparm[0]=1;
    iparm[3]=10*pardiso_cgs_level()+1;
  }
  /* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
     must be done once  */
  if (mthread_mkl == 0) {
    nthread=1;
    env=getenv("MKL_NUM_THREADS");
    if(env) {
      nthread=atoi(env);}
    else {
      env=getenv("OMP_NUM_THREADS");
      if(env) {nthread=atoi(env);}
    }
    env=ccxopt_getenv("CCX_NPROC_EQUATION_SOLVER");
    if(env) {
      nthread_v=atoi(env);
      if (nthread_v <= nthread) {nthread=nthread_v;}
    }
    if (nthread < 1) {nthread=1;}
    sprintf(envMKL,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);  
    putenv(envMKL);
    mthread_mkl=nthread;
  }
    
  printf(" number of threads =% d\n\n",mthread_mkl);

  if(!reuse_current){
    for(i=0;i<64;i++){pt[i]=0;}
  }

  if(*symmetryflag==0){

    /* symmetric matrix; the subdiagonal entries are stored
       column by column in au, the diagonal entries in ad;
       pardiso needs the entries row per row */      

    mtype=-2;
      
    ndim=*neq+*nzs;

    if(!reuse_current){
      NNEW(pointers,ITG,*neq+1);
      NNEW(icolpardiso,ITG,ndim);
      NNEW(aupardiso,double,ndim);
    }
      
    k=ndim;
    l=*nzs;
      
    if(*sigma==0.){
      pointers[*neq]=ndim+1;
      for(i=*neq-1;i>=0;--i){
	for(j=0;j<icol[i];++j){
	  --k;--l;
	  if(!reuse_current) icolpardiso[k]=irow[l];
	  aupardiso[k]=au[l];
	}
	if(!reuse_current) pointers[i]=k;
	k--;
	if(!reuse_current) icolpardiso[k]=i+1;
	aupardiso[k]=ad[i];
      }
    }
    else{
      pointers[*neq]=ndim+1;
      for(i=*neq-1;i>=0;--i){
	for(j=0;j<icol[i];++j){
	  --k;--l;
	  if(!reuse_current) icolpardiso[k]=irow[l];
	  aupardiso[k]=au[l]-*sigma*aub[l];
	}
	if(!reuse_current) pointers[i]=k;
	k--;
	if(!reuse_current) icolpardiso[k]=i+1;
	aupardiso[k]=ad[i]-*sigma*adb[i];
      }
    }
  }else{

    if(*inputformat==3){

      /* off-diagonal terms  are stored column per
	 column from top to bottom in au;
	 diagonal terms are stored in ad  */

      /* structurally and numerically asymmetric */
	
      mtype=11;
	
      ndim=*neq+*nzs;
      NNEW(pointers,ITG,*neq+1);
      NNEW(irowpardiso,ITG,ndim);	  
      NNEW(icolpardiso,ITG,ndim);
      NNEW(aupardiso,double,ndim);
	  
      k=0;
      k2=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  if(au[k]>1.e-12||au[k]<-1.e-12){
	    icolpardiso[k2]=i+1;
	    irowpardiso[k2]=irow[k];
	    aupardiso[k2]=au[k];
	    k2++;		  
	  }
	  k++;	      
	}	  
      }  
      /* diagonal terms */  
      for(i=0;i<*neq;i++){
	icolpardiso[k2]=i+1;
	irowpardiso[k2]=i+1;
	aupardiso[k2]=ad[i];
	k2++;	  
      }
      ndim=k2;
	  
      /* pardiso needs the entries row per row; so sorting is
	 needed */ 
	  
      kflag=2;
      FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			&ndim,&kflag));
	  
      /* sorting each row */
	  
      k=0;
      pointers[0]=1;
      for(i=0;i<*neq;i++){
	j=i+1;
	kstart=k;
	do{
	  if(irowpardiso[k]!=j ){
	    n=k-kstart;		  
	    FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
			     &n,&kflag));
	    pointers[i+1]=k+1;
	    break;  
	  }else{
	    if(k+1==ndim){
	      n=k-kstart+1;	  
	      FORTRAN(isortid,(&icolpardiso[kstart],
			       &aupardiso[kstart],&n,&kflag));
	      break;	       
	    }else{
	      k++;	       
	    }  
	  }
	}while(1);
      }
      pointers[*neq]=ndim+1;
      SFREE(irowpardiso);

    }else if(*inputformat==1){
	  
      /* lower triangular matrix is stored column by column in
	 au, followed by the upper triangular matrix row by row;
	 the diagonal terms are stored in ad */

      /* structurally symmetric, numerically asymmetric */
	
      mtype=1;
	
      /* reordering lower triangular matrix */

      /* THIS BRANCH REBUILDS ITS CSC SCRATCH ON EVERY CALL, unlike the
         symmetric branch above, which guards the allocation with
         !reuse_current.  It cannot be guarded the same way: the buffers are
         allocated at nzs, filled, sorted, and only then RENEWed to
         neq+2*nzs, so the construction is multi-stage.

         Until E-101 that did not matter, because reuse was never eligible
         here and pardiso_cleanup - the only place that frees these three -
         ran after every solve.  Opening the cache to mtype=1 removed that
         cleanup and left the allocation, so EVERY factorisation leaked
         (neq+1)+nzs ITG plus nzs doubles.  Measured on a 190 000-equation
         model: ~120 MB per call, 4.75 GB after ninety seconds, 20.75 GB
         after forty minutes, against a flat 2.8 GB with the cache off.

         Freeing here costs nothing that matters - PARDISO's own symbolic
         state lives in pt and is untouched, which is where the reuse
         benefit actually is. */

      /* FAST PATH: the pattern has not changed, so pointers and
         icolpardiso are already the ones this matrix needs and only the
         values have to be put back in their slots. */

      if((repack_on)&&(reuse_current)&&(pardiso_src!=NULL)&&
         (pardiso_src_n==*neq+2**nzs)&&(aupardiso!=NULL)&&
         (icolpardiso!=NULL)&&(pointers!=NULL)){
        for(k=0;k<pardiso_src_n;k++){
          ITG sidx=pardiso_src[k];
          aupardiso[k]=(sidx>=0)?au[sidx]:ad[-sidx-1];
        }
        pardiso_repack_hits++;

        /* VERIFY: keep what the fast path produced, then rebuild the whole
           CSR from scratch anyway and compare slot by slot.  This is the
           only check that can actually fail - it compares the shortcut
           against the thing it is short-cutting, on real matrices, with no
           model of either.  CCX_PARDISO_REPACK_BREAK corrupts one entry of
           the permutation so that the check can be SEEN to go red. */

        if(repack_verify){
          NNEW(repack_ref,double,pardiso_src_n);
          for(k=0;k<pardiso_src_n;k++) repack_ref[k]=aupardiso[k];
        }else{
          repack_fast=1;
        }
      }

      if(!repack_fast){

      SFREE(pointers);
      SFREE(icolpardiso);
      SFREE(aupardiso);
      pointers=NULL;
      icolpardiso=NULL;
      aupardiso=NULL;
      pardiso_repack_free();

      ndim=*nzs;
      NNEW(pointers,ITG,*neq+1);
      NNEW(irowpardiso,ITG,ndim);
      NNEW(icolpardiso,ITG,ndim);
      NNEW(aupardiso,double,ndim);
	  
      k=0;
      for(i=0;i<*neq;i++){
	for(j=0;j<icol[i];j++){
	  icolpardiso[k]=i+1;
	  irowpardiso[k]=irow[k];
	  aupardiso[k]=repack_on?(double)k:au[k];
	  k++;
	}
      }
	  
      /* pardiso needs the entries row per row; so sorting is
	 needed */
	  
      kflag=2;
      FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			&ndim,&kflag));
	  
      /* sorting each row */
	  
      k=0;
      pointers[0]=1;
      if(ndim>0){
	for(i=0;i<*neq;i++){
	  j=i+1;
	  kstart=k;
	  do{

	    /* end of row reached */

	    if(irowpardiso[k]!=j){
	      n=k-kstart;
	      FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
			       &n,&kflag));
	      pointers[i+1]=k+1;
	      break;
	    }else{

	      /* end of last row reached */

	      if(k+1==ndim){
		n=k-kstart+1;
		FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				 &n,&kflag));
		break;
	      }else{

		/* end of row not yet reached */

		k++;
	      }
	    }
	  }while(1);
	}
      }
      pointers[*neq]=ndim+1;
      SFREE(irowpardiso);

      /* composing the matrix: lower triangle + diagonal + upper triangle */

      ndim=*neq+2**nzs;
      RENEW(icolpardiso,ITG,ndim);
      RENEW(aupardiso,double,ndim);
      k=ndim;
      for(i=*neq-1;i>=0;i--){
	l=k+1;
	for(j=jq[i+1]-1;j>=jq[i];j--){
	  icolpardiso[--k]=irow[j-1];
	  aupardiso[k]=repack_on?(double)(j+*nzs3-1):au[j+*nzs3-1];
	}
	icolpardiso[--k]=i+1;
	aupardiso[k]=repack_on?(double)(-(i+1)):ad[i];
	for(j=pointers[i+1]-1;j>=pointers[i];j--){
	  icolpardiso[--k]=icolpardiso[j-1];
	  aupardiso[k]=aupardiso[j-1];
	}
	pointers[i+1]=l;
      }
      pointers[0]=1;

      /* Materialise: the payload is currently the source index of every
         slot.  Record it, then put the real values through it - which is
         also the first exercise of the fast path's own arithmetic, so a
         wrong index shows up immediately rather than on the next call. */

      if(repack_on){
        pardiso_repack_free();
        NNEW(pardiso_src,ITG,ndim);
        pardiso_src_n=ndim;
        for(k=0;k<ndim;k++) pardiso_src[k]=(ITG)aupardiso[k];
        for(k=0;k<ndim;k++){
          ITG sidx=pardiso_src[k];
          aupardiso[k]=(sidx>=0)?au[sidx]:ad[-sidx-1];
        }
        pardiso_repack_builds++;

        /* The break is applied to the CACHE only, after this rebuild has
           already produced the correct matrix.  Corrupting it earlier
           would corrupt both sides of the comparison equally and the
           check would pass, which is precisely the failure mode a
           deliberate-break switch exists to rule out. */

        {
          const char *brk=ccxopt_getenv("CCX_PARDISO_REPACK_BREAK");
          if((brk!=NULL)&&(ndim>1)){
            ITG kb=(ITG)atoi(brk);
            if((kb<0)||(kb>=ndim)) kb=ndim/2;
            pardiso_src[kb]=pardiso_src[(kb+1)%ndim];
            printf("[PARDISO REPACK] DELIBERATELY BROKEN: cache slot %"
                   ITGFORMAT " now reads its neighbour's source.  The "
                   "verify check must go red on the next reuse.\n",kb);
            fflush(stdout);
          }
        }
      }

      if(repack_ref!=NULL){
        ITG nbad=0;
        double worst=0.;
        for(k=0;k<ndim;k++){
          double d=fabs(repack_ref[k]-aupardiso[k]);
          if(d>0.){ nbad++; if(d>worst) worst=d; }
        }
        if(nbad==0){
          printf("[PARDISO REPACK] verify: %" ITGFORMAT " slot(s) identical "
                 "to a full rebuild\n",ndim);
        }else{
          printf("[PARDISO REPACK] verify *ERROR: %" ITGFORMAT " of %"
                 ITGFORMAT " slot(s) differ from a full rebuild, worst "
                 "%.6e.  The cached permutation does not describe this "
                 "matrix.\n",nbad,ndim,worst);
        }
        fflush(stdout);
        SFREE(repack_ref); repack_ref=NULL;
      }

      }   /* end of !repack_fast */
    }
  }

/* next line is for the simulateous use of PARDISO and PaStiX */

  mkl_domain_set_num_threads(mthread_mkl,MKL_DOMAIN_PARDISO);

  if(cgs_active){
    b=pardiso_cgs_rhs;
    nrhs=pardiso_cgs_nrhs;
    NNEW(x,double,nrhs**neq);
  }

  /* Name the event by the PHASE actually executed, so one run measures what
     the symbolic analysis costs instead of needing an A/B.  phase 12 is
     analyse+factorise, 22 is factorise against a retained analysis, 23 is
     the combined factorise-and-solve the CGS path uses.  The difference
     between the means of 12 and 22 IS the analysis, and that is the number
     the fixed-sparsity candidate in 08-OBJECT-MODEL.md section 4 turns on. */
  {
    const char *ph=(phase==12)?"pardiso phase 12 (analyse+numeric)":
                   ((phase==22)?"pardiso phase 22 (numeric only)":
                                "pardiso phase 23 (numeric+solve)");
    logview_begin_named(ph);
    FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		     pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                     b,x,&error));
    logview_end_named(ph);
  }

  if(cgs_active){
    if(error==0){
      for(i=0;i<nrhs**neq;i++){pardiso_cgs_rhs[i]=x[i];}
      pardiso_cgs_done=1;
      /* iparm[19]>0 is the number of CGS iterations that converged;
         <0 means the Krylov iteration failed and PARDISO refactorised. */
      if(iparm[19]>0){
        pardiso_cgs_ok++;
        pardiso_cgs_iter+=iparm[19];
      }else{
        pardiso_cgs_fail++;
      }
      pardiso_cgs_report(0);
    }
    SFREE(x);
    x=NULL;
  }

  if(reuse_requested){
    if(error==0){
      if(!reuse_current){
        printf("[PARDISO CACHE] symbolic analysis retained "
               "(neq=%" ITGFORMAT " nzs=%" ITGFORMAT ")\n",
               *neq,*nzs);
      }
      pardiso_cache_valid=1;
      pardiso_cache_neq=*neq;
      pardiso_cache_nzs=*nzs;
      pardiso_cache_symmetry=*symmetryflag;
      pardiso_cache_inputformat=*inputformat;
      pardiso_cache_hash=structure_hash;
    }else{
      pardiso_cache_valid=0;
    }
  }

  logview_end_named("pardiso factor");
  return;
}

void pardiso_solve(double *b, ITG *neq,ITG *symmetryflag,ITG *inputformat,
		   ITG *nrhs){
  logview_begin_named("pardiso solve");

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
    msglvl=0,i,error=0;
  double *x=NULL;

  /*  if(*symmetryflag==0){
    printf(" Solving the system of equations using the symmetric pardiso solver\n");
  }else{
    printf(" Solving the system of equations using the unsymmetric pardiso solver\n");
    }*/

  if(*symmetryflag==0){
    mtype=-2;
  }else{
    if(*inputformat==3){
      mtype=11;
    }else{
      mtype=1;
    }
  }
  iparm[1]=3;
  
  /* pardiso_factor has been called before, MKL_NUM_THREADS=mthread_mkl is set*/

  //  printf(" number of threads =% d\n\n",mthread_mkl);

  NNEW(x,double,*nrhs**neq);

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,nrhs,iparm,&msglvl,
                   b,x,&error));

  pardiso_refine_note(iparm[6]);

  for(i=0;i<*nrhs**neq;i++){b[i]=x[i];}
  SFREE(x);

  logview_end_named("pardiso solve");
  return;
}

void pardiso_cleanup(ITG *neq,ITG *symmetryflag,ITG *inputformat){
  logview_begin_named("pardiso cleanup");

  ITG maxfct=1,mnum=1,phase=-1,*perm=NULL,nrhs=1,mtype,
    msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0){
    mtype=-2;
  }else{
    if(*inputformat==3){
      mtype=11;
    }else{
      mtype=1;
    }
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  SFREE(icolpardiso);
  SFREE(aupardiso);
  SFREE(pointers);
  icolpardiso=NULL;
  aupardiso=NULL;
  pointers=NULL;
  pardiso_cache_valid=0;
  pardiso_cache_hash=0;
  pardiso_repack_free();   /* the permutation describes a CSR that is gone */
  pardiso_cgs_report(1);
  pardiso_cgs_ok=0;pardiso_cgs_fail=0;pardiso_cgs_iter=0;

  logview_end_named("pardiso cleanup");
  return;
}

void pardiso_main(double *ad, double *au, double *adb, double *aub, 
		  double *sigma,double *b, ITG *icol, ITG *irow, 
		  ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		  ITG *jq, ITG *nzs3,ITG *nrhs){

  if(*neq==0) return;

  /* Offer the RHS to pardiso_factor: if a cached LU is usable it will run the
     combined phase 23 and solve there, and pardiso_cgs_done says so. */
  pardiso_cgs_rhs=b;
  pardiso_cgs_nrhs=*nrhs;
  pardiso_cgs_done=0;

  pardiso_factor(ad,au,adb,aub,sigma,icol,irow,
		 neq,nzs,symmetryflag,inputformat,jq,nzs3);

  pardiso_cgs_rhs=NULL;

  if(!pardiso_cgs_done){
    pardiso_solve(b,neq,symmetryflag,inputformat,nrhs);
  }

  if(!pardiso_reuse_eligible(*symmetryflag,*inputformat)){
    pardiso_cleanup(neq,symmetryflag,inputformat);
  }

  return;
}

#endif
