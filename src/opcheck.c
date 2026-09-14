/*     Operator: does the assembled tangent differ the residual?
 *
 *     This is the named hole in the diagnostics, and the Operator row
 *     carries the note
 *     "nothing verifies it".  Both were wrong about the tooling and right
 *     about the state of knowledge: a column-by-column finite-difference
 *     check has been sitting in nonlingeo.c the whole time, behind three
 *     environment names that no test set and no document mentioned.  It has
 *     an owner now, and a self test, and one thing it did not have before.
 *
 *     WHY A CENTRAL DIFFERENCE IS NOT ENOUGH, which is the whole point of
 *     this file.
 *
 *     The first run of that probe found the assembled tangent exact in the
 *     elastic bulk (5.4e-08 relative) and wrong by up to 9.6e-02 in the
 *     process zone, with the answer FLAT over h from 1e-11 to 1e-8.  A
 *     converged plateau looks like proof that the difference is real and not
 *     truncation.  It is not proof of anything of the kind:
 *
 *       at a point of non-smoothness the central difference converges to the
 *       MEAN of the one-sided derivatives.
 *
 *     That mean is a stable, converged number, different from both branches,
 *     and its plateau in h is indistinguishable from the plateau of a
 *     genuinely wrong tangent.  Three kinks are documented in this very
 *     process zone - damage initiation (deff>d0), loading against unloading
 *     (deff against dmax), and crack-face closure (deltal(1)=0, a measured
 *     factor of 1e+06 in the normal slope) - so the ambiguity is not
 *     hypothetical here, it is the default expectation.
 *
 *     THE DISCRIMINATING TEST.  Take the two one-sided differences
 *     separately:
 *
 *       K_fwd = ( f(u + h e) - f(u) ) / h        K_bwd = ( f(u) - f(u - h e) ) / h
 *
 *       a KINK:           K_fwd and K_bwd disagree with EACH OTHER, and one
 *                         of them agrees with K_asm - the tangent is right,
 *                         it is right on one branch, and Newton is standing
 *                         on the switching surface.
 *       a WRONG TANGENT:  K_fwd and K_bwd agree with each other, and both
 *                         disagree with K_asm.  There is nothing non-smooth
 *                         here and the operator is simply not the
 *                         differential.
 *       both:             one-sided differences disagree AND neither matches
 *                         K_asm.  Non-smooth, and the tangent is not on
 *                         either branch.
 *
 *     Everything is measured relative to the SCALE of the column - the
 *     largest |K| in it - because an absolute error on a coefficient that is
 *     itself a millionth of the column means nothing.
 *
 *     This file classifies.  It reads no environment, decides no step, and
 *     is on no solution path.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "CalculiX.h"

/* One coefficient of the assembled operator, out of CalculiX's own storage.

   Returns 0 for a structurally absent coefficient, which is the correct
   value - the sparsity pattern is a superset of the assembled entries and a
   missing slot means the two degrees of freedom share no element.  Moved
   here verbatim from nonlingeo.c, where it had no other caller. */
static double opcheck_coeff(const double *ad,const double *au,
                            const ITG *jq,const ITG *irow,const ITG *nzs,
                            ITG nasym,ITG i,ITG j)
{
  ITG col,want,lo,hi,mid,off;

  if(i==j) return ad[i-1];

  /* nasym says whether au carries an upper triangle at all.  When it does
     not - the stock symmetric path, which is what a deck without
     CCX_DAMAGE_TANGENT=UNSYM runs - the coefficient above the diagonal is
     the one below it, and reading at offset nzs[2] would run off the end of
     the array and return whatever is there.  That is not a hypothetical:
     the first run of this check on test/pathfollow/close.inp, which is
     symmetric, reported the cohesive tangent wrong by up to 1.0 of the
     column scale, and every one of those numbers came from past the end of
     au. */
  if(i>j){col=j;want=i;off=0;}
  /* i<j and no upper triangle: the value wanted is its mirror K(j,i), which
     lives in the LOWER triangle of column i - the smaller index is always
     the column there. */
  else if(nasym==0){col=i;want=j;off=0;}
  else   {col=i;want=j;off=nzs[2];}

  lo=jq[col-1];hi=jq[col]-1;
  while(lo<=hi){
    mid=(lo+hi)/2;
    if(irow[mid-1]==want) return au[mid-1+off];
    if(irow[mid-1]<want) lo=mid+1; else hi=mid-1;
  }
  return 0.;
}

/* The verdict for one coefficient.  Everything else in this file is
   bookkeeping around these four lines, so they are the ones the self test
   exercises directly. */
ITG opcheck_classify(double kfwd,double kbwd,double kasm,double scale,
                     double tol){
  double t=tol*scale;
  double dside=fabs(kfwd-kbwd);
  double dctr =fabs(0.5*(kfwd+kbwd)-kasm);
  double dfwd =fabs(kfwd-kasm);
  double dbwd =fabs(kbwd-kasm);
  double dbest=(dfwd<dbwd)?dfwd:dbwd;

  if(dside<=t) return (dctr<=t)?OPCHECK_OK:OPCHECK_WRONG;
  return (dbest<=t)?OPCHECK_KINK:OPCHECK_BOTH;
}

void opcheck_begin(opcheck *o,double tol){
  o->n=0;o->nok=0;o->nkink=0;o->nwrong=0;o->nboth=0;
  o->scale=0.;o->relctr=0.;o->relside=0.;o->relbest=0.;
  o->worst=-1;o->worstdir=0;o->wfwd=0.;o->wbwd=0.;o->wasm=0.;
  o->tol=tol;
}

/*  One column of the operator, compared against the two one-sided
 *  differences of the internal force.
 *
 *  Two passes on purpose.  The scale of a column is the largest |K| IN it,
 *  and it cannot be known until the column has been read; classifying
 *  against a running maximum - which is what the original inline code did
 *  for its "bad" count - makes the verdict depend on the order the rows
 *  happen to come in.  Two passes over a few thousand doubles costs nothing
 *  next to the three residual evaluations that produced them.
 */
void opcheck_column(opcheck *o,ITG nk,ITG mt,const ITG *nactdof,
                    const double *f0,const double *fp,const double *fm,
                    double h,ITG col,
                    const double *ad,const double *au,
                    const ITG *jq,const ITG *irow,const ITG *nzs,ITG nasym){
  ITG k,d,row,kind;
  double kfwd,kbwd,kasm,ctr,r;

  o->scale=0.;
  for(k=0;k<nk;k++){
    for(d=1;d<=3;d++){
      if(nactdof[mt*k+d]<=0) continue;
      ctr=0.5*((fp[mt*k+d]-f0[mt*k+d])+(f0[mt*k+d]-fm[mt*k+d]))/h;
      if(fabs(ctr)>o->scale) o->scale=fabs(ctr);
    }
  }
  if(o->scale<=0.) return;

  for(k=0;k<nk;k++){
    for(d=1;d<=3;d++){
      row=nactdof[mt*k+d];
      if(row<=0) continue;
      kfwd=(fp[mt*k+d]-f0[mt*k+d])/h;
      kbwd=(f0[mt*k+d]-fm[mt*k+d])/h;
      kasm=opcheck_coeff(ad,au,jq,irow,nzs,nasym,row,col);
      o->n++;
      kind=opcheck_classify(kfwd,kbwd,kasm,o->scale,o->tol);
      if(kind==OPCHECK_OK)         o->nok++;
      else if(kind==OPCHECK_KINK)  o->nkink++;
      else if(kind==OPCHECK_WRONG) o->nwrong++;
      else                         o->nboth++;

      r=fabs(kfwd-kbwd)/o->scale;
      if(r>o->relside) o->relside=r;
      r=fabs(0.5*(kfwd+kbwd)-kasm)/o->scale;
      if(r>o->relctr){
        o->relctr=r;o->worst=k+1;o->worstdir=d;
        o->wfwd=kfwd;o->wbwd=kbwd;o->wasm=kasm;
      }
      r=fabs(kfwd-kasm);
      if(fabs(kbwd-kasm)<r) r=fabs(kbwd-kasm);
      r/=o->scale;
      if(r>o->relbest) o->relbest=r;
    }
  }
}

/*  Self test.  The classifier is four comparisons and every one of them is a
 *  claim about physics, so each is exercised in both directions.  A verdict
 *  that cannot be shown to change is not a verdict. */
/* The coefficient reader, against a hand-built matrix in CalculiX's own
   storage.  This is the part that was wrong, and it was wrong in the way
   that is hardest to notice: it returned a plausible number read from past
   the end of an array.  Four equations; column 1 holds rows 2 and 3, column
   2 holds row 3, columns 3 and 4 hold nothing. */
static ITG opcheck_coeff_selftest(void){
  ITG bad=0;
  ITG jq[5]={1,3,4,4,4};
  ITG irow[3]={2,3,3};
  ITG nzs[3]={3,3,3};
  double ad[4]={11.,22.,33.,44.};
  double au[6]={21.,31.,32., 12.,13.,23.};   /* lower, then upper at nzs[2] */
  struct{ITG i,j,nasym;double want;const char *what;}c[]={
    {1,1,1,11.,"the diagonal"},
    {2,1,1,21.,"below the diagonal"},
    {3,1,1,31.,"below the diagonal, second row of the column"},
    {3,2,1,32.,"below the diagonal, second column"},
    {1,2,1,12.,"ABOVE the diagonal, unsymmetric storage"},
    {1,3,1,13.,"above the diagonal, second row"},
    {2,3,1,23.,"above the diagonal, second column"},
    {1,2,0,21.,"above the diagonal, SYMMETRIC storage: the mirror"},
    {2,3,0,32.,"above the diagonal, symmetric, second column"},
    {4,3,1,0. ,"a structurally absent coefficient"},
    {4,3,0,0. ,"a structurally absent coefficient, symmetric"},
    {0,0,0,0.,NULL}};
  ITG n;
  for(n=0;c[n].what!=NULL;n++){
    double got=opcheck_coeff(ad,au,jq,irow,nzs,c[n].nasym,c[n].i,c[n].j);
    if(got!=c[n].want){
      printf("[OPCHECK] *ERROR: %s: K(%" ITGFORMAT ",%" ITGFORMAT
             ") nasym=%" ITGFORMAT " read %g, want %g\n",
             c[n].what,c[n].i,c[n].j,c[n].nasym,got,c[n].want);bad++;}
  }
  return bad;
}

ITG opcheck_selftest(void){
  ITG bad=0,k;
  const double s=1000.,tol=1.e-4;   /* threshold is 0.1 in absolute terms */

  bad+=opcheck_coeff_selftest();

  /* smooth and correct: both one-sided differences and the tangent agree */
  k=opcheck_classify(500.,500.,500.,s,tol);
  if(k!=OPCHECK_OK){
    printf("[OPCHECK] *ERROR: an exact column was not called OK (%" ITGFORMAT
           ")\n",k);bad++;}

  /* smooth and wrong: the one-sided differences agree with EACH OTHER and
     both differ from the tangent.  This is the case the central difference
     alone cannot separate from the next one. */
  k=opcheck_classify(500.,500.,300.,s,tol);
  if(k!=OPCHECK_WRONG){
    printf("[OPCHECK] *ERROR: a smooth mismatch was not called WRONG (%"
           ITGFORMAT ")\n",k);bad++;}

  /* a kink with the tangent sitting on the FORWARD branch.  Note that the
     central difference here is 400, which differs from the tangent by
     exactly as much as the WRONG case above - the two are indistinguishable
     to it, and that is the whole reason this file exists. */
  k=opcheck_classify(300.,500.,300.,s,tol);
  if(k!=OPCHECK_KINK){
    printf("[OPCHECK] *ERROR: a kink on the forward branch was not called "
           "KINK (%" ITGFORMAT ")\n",k);bad++;}
  /* and on the backward branch, because a test that only covers one side
     would pass with the min() written as either argument */
  k=opcheck_classify(300.,500.,500.,s,tol);
  if(k!=OPCHECK_KINK){
    printf("[OPCHECK] *ERROR: a kink on the backward branch was not called "
           "KINK (%" ITGFORMAT ")\n",k);bad++;}

  /* non-smooth AND on neither branch */
  k=opcheck_classify(300.,500.,900.,s,tol);
  if(k!=OPCHECK_BOTH){
    printf("[OPCHECK] *ERROR: a kink with the tangent on neither branch was "
           "not called BOTH (%" ITGFORMAT ")\n",k);bad++;}

  /* the scale matters: the SAME absolute discrepancy is a mismatch in a
     small column and noise in a large one.  Without this the classifier
     would report the whole process zone as broken on rounding alone. */
  k=opcheck_classify(500.,500.,500.05,s,tol);
  if(k!=OPCHECK_OK){
    printf("[OPCHECK] *ERROR: a discrepancy below the column scale was not "
           "tolerated\n");bad++;}
  k=opcheck_classify(0.5,0.5,0.55,1.,tol);
  if(k!=OPCHECK_WRONG){
    printf("[OPCHECK] *ERROR: the same discrepancy in a small column was "
           "tolerated\n");bad++;}

  /* Exactly at the threshold is not over it.  Written in powers of two on
     purpose: with s=1000 and tol=1e-4 the product is 0.1, which is not
     representable, and 500.1-500. comes out as 0.10000000000002274 - so the
     obvious version of this check fails on the arithmetic rather than on
     the logic.  It did, on the first run, which is the self test earning
     its keep at the cost of one build. */
  k=opcheck_classify(512.,512.,513.,1024.,1./1024.);
  if(k!=OPCHECK_OK){
    printf("[OPCHECK] *ERROR: a value exactly at the threshold was called a "
           "mismatch\n");bad++;}
  k=opcheck_classify(512.,512.,514.,1024.,1./1024.);
  if(k!=OPCHECK_WRONG){
    printf("[OPCHECK] *ERROR: a value at twice the threshold was tolerated\n");
    bad++;}
  return bad;
}

/* The values the thirty-nine driver locals carried at their declarations
   in nonlingeo().  The probe is opt-in and every one of them is inert
   until a CCX_STRUCT_FD_* or CCX_DAMAGE_TANGENT_* switch arms it. */
void opcheckdrv_init(opcheckdrv *p)
{
  memset(p,0,sizeof(*p));
  p->fd_ad=NULL;
  p->fd_au=NULL;
  p->fd_f0=NULL;
  p->fd_fm=NULL;
  p->fd_fp=NULL;
  p->fd_vsav=NULL;
  p->fd_vtrue=NULL;
  p->fd_el=-1;
  p->fd_it=1;
  p->fd_tel=-1;
  p->fd_uel=-1;
  p->fd_h=1.e-7;
  p->fd_udmax=-1.;
  p->unsym_advrep=-1;
}

/* ---- arming ------------------------------------------------------------

   Read unconditionally, beside the object it configures.  The comment the
   block carries is the reason it is unconditional, and it is worth more
   than the six lines it explains. */

void opcheckdrv_configure_fd(opcheckdrv *p)
{
  const char *e;

  /* The operator check reads its configuration HERE, unconditionally.

     It used to be parsed inside the block gated by damage_de12_enabled -
     that is, only on a deck carrying a progressive BULK damage material -
     for no reason except that it was written next to the code that needed
     that gate.  The consequence was measured rather than argued: on
     test/pathfollow/close.inp, the one deck in this tree that isolates the
     crack-face closure kink, the probe could not be armed at all.  Every
     one of its switches was reported by [SWITCHES LEFT] as set and never
     read, which is exactly the failure that report exists to catch, on its
     first real use.

     "A responsibility with no home ends up nested inside whatever code
     happened to be nearby" - describing a
     different instance of the same thing. */
  if((e=ccxopt_getenv("CCX_STRUCT_FD_INC"))!=NULL)
    p->fd_inc=atoi(e);
  if((e=ccxopt_getenv("CCX_STRUCT_FD_ITER"))!=NULL)
    p->fd_it=atoi(e);
  if((e=ccxopt_getenv("CCX_STRUCT_FD_H"))!=NULL)
    p->fd_h=atof(e);
  if((e=ccxopt_getenv("CCX_STRUCT_FD_STEP"))!=NULL)
    p->fd_step=atoi(e);
  if((e=ccxopt_getenv("CCX_STRUCT_FD_ELEM"))!=NULL)
    p->fd_pick=atoi(e);
  if((e=ccxopt_getenv("CCX_STRUCT_FD_BASE"))!=NULL)
    p->fd_base=((strcmp(e,"VOLD")==0)||
                    (strcmp(e,"vold")==0))?1:0;
}

/* ---- the probe ---------------------------------------------------------

   The measurement this file exists for, moved beside the classification it
   feeds.  238 lines that compare the ASSEMBLED tangent, column by column,
   against a central difference of the internal force, and split the answer
   by population: all columns, the one chosen column, the bulk, the
   cohesive facets.

   It is a DIAGNOSTIC and it ends the run.  It perturbs the displacement
   repeatedly and restores it, so nothing downstream may depend on what it
   leaves behind; the guard that decides whether it fires at all stays with
   the caller, like every other mechanism here.

   damcat is the rank-1 census's per-element population, filled elsewhere
   and read here, which is why it is a parameter and not state: this file
   does not decide which population an element is in, it reports the
   operator error per population.                                      */

void opcheck_probe(opcheckdrv *o,const trialctx *t,const ITG *ndmat_,
                   const ITG *damcat,ITG iit)
{
  double *v=*(t->v),*vold=*(t->vold),*fn=*(t->fn),*dam=*(t->dam);
  double *xstate=*(t->xstate);
  ITG *nk=*(t->nk),*mi=*(t->mi),*ipkon=*(t->ipkon),*kon=*(t->kon);
  ITG *nactdof=*(t->nactdof),*nstate_=*(t->nstate_);
  ITG *irow=*(t->irow),*jq=*(t->jq),*nzs=*(t->nzs),*istep=*(t->istep);
  char *lakon=*(t->lakon);
  ITG iinc=*(t->iinc),ne0=*(t->ne0),nasym=*(t->nasym),mt=mi[1]+1;
  ITG i,k,idir;

  /* Aliasing v and fn is safe HERE and would not be everywhere: this probe
     only ever calls trial_results(), which reads and writes through the
     scratch arrays that are there.  trial_evaluate() and trial_residual()
     free and reallocate them, and an alias taken before one of those is
     dangling after it.  CalculiX.h says which is which. */

  if((o->fd_inc>0)&&(iinc>=o->fd_inc)&&(iit>=o->fd_it)&&
     ((o->fd_step<=0)||(*istep==o->fd_step))){

    opcheck oc_all,oc_col_o,oc_bulk,oc_coh;

    NNEW(o->fd_vsav,double,mt**nk);
    NNEW(o->fd_vtrue,double,mt**nk);
    NNEW(o->fd_fp,double,mt**nk);
    NNEW(o->fd_fm,double,mt**nk);
    NNEW(o->fd_f0,double,mt**nk);
    /* v must come back exactly as it was, whatever the probe uses as its
       base, or the solve continues from a state the probe invented. */
    memcpy(o->fd_vtrue,v,sizeof(double)*mt**nk);

    /* WHICH STATE TO DIFFERENTIATE AROUND.  This decides whether the
       comparison means anything at all.

       The matrix was assembled from xstiff at the top of this Newton
       iteration, where the iterate is vold.  The solve then produced b,
       and the results() call just above evaluated the residual at
       v = vold + b.  Probing around v therefore compares K(vold) against
       dR/du(vold+b) - an offset of one Newton step, which would look
       exactly like a wrong tangent: smooth, h-independent, and growing
       with how fast the state moves.

       vold is still untouched here; nonlingeo advances it 2600 lines
       further down.  So CCX_STRUCT_FD_BASE=VOLD differentiates around the
       state the matrix actually came from, and the difference between the
       two settings IS the size of that objection. */
    memcpy(o->fd_vsav,
           (o->fd_base==1)?vold:v,sizeof(double)*mt**nk);

    /* columns are taken at the nodes of the most damaged element: a
       tangent that is right in the elastic bulk and wrong in the
       process zone is precisely the error every earlier test would
       have missed */

    /* Which population is each element in?  damcat is filled by
       mafilldamas and is the only thing that can attribute a measured
       discrepancy on a NAMED element to a named reason.  It is written
       only under CCX_DAMAGE_TANGENT=UNSYM; in the other modes there is
       no rank-1 term anywhere and the census would be a row of zeros
       pretending to mean something, so say that instead. */
    if(damcat!=NULL){
      ITG oc_cc[8],oc_c;
      for(oc_c=0;oc_c<8;oc_c++) oc_cc[oc_c]=0;
      for(i=0;i<ne0;i++){
        oc_c=damcat[i];
        if((oc_c>=0)&&(oc_c<8)) oc_cc[oc_c]++;
      }
      printf("[OPCHECK] rank-1 term by population: assembled=%" ITGFORMAT
             " pre-initiation=%" ITGFORMAT " not-advancing=%" ITGFORMAT
             " on-floor=%" ITGFORMAT " terminal-hole=%" ITGFORMAT
             " advancing-live=%" ITGFORMAT " degenerate=%" ITGFORMAT "\n",
             oc_cc[1],oc_cc[2],oc_cc[3],
             oc_cc[4],oc_cc[5],oc_cc[6],
             oc_cc[7]);
      for(oc_c=1;oc_c<8;oc_c++){
        ITG shown=0;
        if((oc_c==2)||(oc_c==3)) continue;
        if(oc_cc[oc_c]==0) continue;
        printf("[OPCHECK]   category %" ITGFORMAT " elements:",oc_c);
        for(i=0;(i<ne0)&&(shown<10);i++){
          if(damcat[i]!=oc_c) continue;
          printf(" %" ITGFORMAT,i+1);shown++;
        }
        printf("%s\n",(oc_cc[oc_c]>shown)?" ...":"");
      }
      fflush(stdout);
    }else{
      printf("[OPCHECK] no rank-1 census: mafilldamas runs only under "
             "CCX_DAMAGE_TANGENT=UNSYM, so in this mode NO element "
             "carries the -sigma_eff (x) dD/d(eps) term at all\n");
    }

    o->fd_el=-1;o->fd_dmax=-1.;
    /* A NAMED element needs no damage array: the decisive comparison is
       the same element probed on a deck with the damage module and on one
       without it, and the second has no dam[] at all. */
    for(i=0;(i<ne0)&&((o->fd_pick>0)||((*ndmat_>0)&&(dam!=NULL)));i++){
      if(ipkon[i]<0) continue;
      if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
      if(o->fd_pick>0){
        if(i+1!=o->fd_pick) continue;
        o->fd_dmax=(dam!=NULL)?dam[mi[0]*i]:-1.;
        o->fd_el=i;break;
      }
      if(dam[mi[0]*i]>o->fd_dmax){
        o->fd_dmax=dam[mi[0]*i];o->fd_el=i;
      }
    }
    /* And the most damaged COHESIVE facet, because that is where the
       three documented kinks live - damage initiation (deff>d0),
       loading against unloading (deff against dmax) and crack-face
       closure (deltal(1)=0, a measured factor of 1e+06 in the normal
       slope).  Probing only the bulk cannot see any of them, and a
       verdict of "no kink" from bulk columns alone would be a statement
       about where the probe looked rather than about the model.
       UC6 damage is xstate[nstate_*(mi[0]*elem+ip)+1] over 3 points. */
    o->fd_uel=-1;o->fd_udmax=-1.;
    if(*nstate_>1){
      for(i=0;i<ne0;i++){
        if(ipkon[i]<0) continue;
        if((lakon[8*i]!='U')||(lakon[8*i+1]!='C')||(lakon[8*i+2]!='6'))
          continue;
        for(k=0;k<3;k++){
          if(xstate[*nstate_*(mi[0]*i+k)+1]>o->fd_udmax){
            o->fd_udmax=xstate[*nstate_*(mi[0]*i+k)+1];
            o->fd_uel=i;
          }
        }
      }
    }
    if((o->fd_el<0)&&(o->fd_uel<0)){
      printf("[OPCHECK] no active element to probe\n");
    }else{
      printf("[OPCHECK] inc=%" ITGFORMAT " iter=%" ITGFORMAT
  	   " step=%" ITGFORMAT " h=%.3e base=%s\n",iinc,iit,*istep,
         o->fd_h,(o->fd_base==1)?
         "VOLD (the state the matrix was assembled at)":
         "V (the post-solve iterate)");
      if(o->fd_el>=0)
        printf("[OPCHECK]   bulk     element %" ITGFORMAT " dam=%.6f "
               "population=%" ITGFORMAT "%s\n",
               o->fd_el+1,o->fd_dmax,
               (damcat!=NULL)?damcat[o->fd_el]:-1,
               (damcat!=NULL)?
               ((damcat[o->fd_el]==1)?" (ASSEMBLED)":
                (damcat[o->fd_el]==4)?" (on the floor - zero is correct)":
                (damcat[o->fd_el]==5)?" (terminal hole)":
                (damcat[o->fd_el]==6)?" (advancing, no term)":
                (damcat[o->fd_el]==7)?" (degenerate, dropped)":
                " (no term needed)"):" (no census in this mode)");
      if(o->fd_uel>=0)
        printf("[OPCHECK]   cohesive element %" ITGFORMAT " dv=%.6f "
               "(g=%.3e)\n",o->fd_uel+1,o->fd_udmax,
               1.-o->fd_udmax);
      printf("[OPCHECK] one-sided differences taken SEPARATELY: a "
             "central difference converges to their mean at a kink, so "
             "it cannot tell a kink from a wrong tangent.\n");
      fflush(stdout);
    }
    opcheck_begin(&oc_all,1.e-4);
    opcheck_begin(&oc_bulk,1.e-4);
    opcheck_begin(&oc_coh,1.e-4);

    o->fd_ncol=0;
    /* target 0 is the bulk element (4 nodes), target 1 the facet (6) */
    for(o->fd_t=0;o->fd_t<2;o->fd_t++){
    o->fd_tel=(o->fd_t==0)?o->fd_el:o->fd_uel;
    o->fd_tnn=(o->fd_t==0)?4:6;
    for(o->fd_j=0;(o->fd_j<o->fd_tnn)&&(o->fd_tel>=0);
        o->fd_j++){
      o->fd_node=kon[ipkon[o->fd_tel]+o->fd_j]-1;
      if((o->fd_node<0)||(o->fd_node>=*nk)) continue;

      for(idir=1;idir<=3;idir++){
        o->fd_col=nactdof[mt*o->fd_node+idir];
        if(o->fd_col<=0) continue;

        /* central difference on the internal force */

        /* three evaluations, not two: the base state is needed to
           split the central difference into its two one-sided halves,
           which is the only reading that separates a kink from a
           tangent that is not the differential. */
        for(o->fd_s=0;o->fd_s<3;o->fd_s++){
  	memcpy(v,o->fd_vsav,sizeof(double)*mt**nk);
  	if(o->fd_s==0)      v[mt*o->fd_node+idir]+=o->fd_h;
  	else if(o->fd_s==1) v[mt*o->fd_node+idir]-=o->fd_h;
  	trial_results(t);
  	if(o->fd_s==0)
  	  memcpy(o->fd_fp,fn,sizeof(double)*mt**nk);
  	else if(o->fd_s==1)
  	  memcpy(o->fd_fm,fn,sizeof(double)*mt**nk);
  	else
  	  memcpy(o->fd_f0,fn,sizeof(double)*mt**nk);
        }

        /* Classify the column: opcheck.c owns the comparison and the
           verdict, this loop owns only the three evaluations that feed
           it. */

        opcheck_begin(&oc_col_o,1.e-4);
        opcheck_column(&oc_col_o,*nk,mt,nactdof,
                       o->fd_f0,o->fd_fp,o->fd_fm,
                       o->fd_h,o->fd_col,
                       o->fd_ad,o->fd_au,jq,irow,nzs,nasym);
        monitor_opcheck(&oc_col_o,iinc,iit,o->fd_node+1,idir,
                        o->fd_h);
        oc_all.n+=oc_col_o.n;
        oc_all.nok+=oc_col_o.nok;
        oc_all.nkink+=oc_col_o.nkink;
        oc_all.nwrong+=oc_col_o.nwrong;
        oc_all.nboth+=oc_col_o.nboth;
        if(o->fd_t==0){
          oc_bulk.n+=oc_col_o.n;
          oc_bulk.nok+=oc_col_o.nok;
          oc_bulk.nkink+=oc_col_o.nkink;
          oc_bulk.nwrong+=oc_col_o.nwrong;
          oc_bulk.nboth+=oc_col_o.nboth;
        }else{
          oc_coh.n+=oc_col_o.n;
          oc_coh.nok+=oc_col_o.nok;
          oc_coh.nkink+=oc_col_o.nkink;
          oc_coh.nwrong+=oc_col_o.nwrong;
          oc_coh.nboth+=oc_col_o.nboth;
        }
        fflush(stdout);
        o->fd_ncol++;
      }
    }
    }

    memcpy(v,o->fd_vtrue,sizeof(double)*mt**nk);
    SFREE(o->fd_vsav);SFREE(o->fd_vtrue);
    SFREE(o->fd_fp);SFREE(o->fd_fm);SFREE(o->fd_f0);
    if(oc_bulk.n>0){
      printf("[OPCHECK] --- bulk columns ---\n");
      monitor_opcheck_total(&oc_bulk,iinc,iit,o->fd_el+1,
                            o->fd_dmax,o->fd_h,o->fd_ncol);
    }
    if(oc_coh.n>0){
      printf("[OPCHECK] --- cohesive columns, where the three documented "
             "kinks live ---\n");
      monitor_opcheck_total(&oc_coh,iinc,iit,o->fd_uel+1,
                            o->fd_udmax,o->fd_h,o->fd_ncol);
    }
    printf("[OPCHECK] stopping: the probe perturbed v repeatedly and "
           "this run is diagnostic only\n");
    fflush(stdout);
    FORTRAN(stopwithout201,());
  }
}
