/*     Operator: does the assembled tangent differ the residual?
 *
 *     07-RESEARCH-AGENDA.md rank 3 calls this "the named hole in the
 *     diagnostics", and 08-OBJECT-MODEL.md gives the Operator row the note
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
