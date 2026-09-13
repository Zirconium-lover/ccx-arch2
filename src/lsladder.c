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
  THE BACKTRACKING LADDER OF THE ADAPTIVE DAMAGE LINE SEARCH
  ==========================================================

  This is the decision the line search makes, on its own, with no model
  attached: given the residual measured at a trial step length, what step
  length to try next, and which one to finally take.  It lives here rather
  than inline in the Newton loop because it was WRONG and the way it was
  wrong is worth a regression test.

  WHAT WAS MEASURED
  -----------------
  On the s3rad target, over the 245 line-search activations of a run to
  the recorded wall (CCX_DAMAGE_LS_PROBE, which walks a finer ladder and
  prints the residual at each step without changing what the solver does):

    the ladder was {1.0, 0.5, 0.1} - three trials with a floor of 0.1;

    in 114 of the 245 activations (46.5%) the BEST step length lay below
    0.1, i.e. outside anything the ladder could reach;

    when nothing contracted, the loop fell out with the LAST trial - the
    floor - in hand, whatever it had done to the residual.  At the wall
    increment, iterations 12 to 16, that floor step raised the residual
    every time (2.852e-03, 3.053e-03, 3.353e-03, 3.649e-03, 3.822e-03)
    while alpha = 0.03 lowered it every time.  The iteration had been
    converging up to iteration 11; the ladder turned it round, and the
    step-time controller then stopped the run for "too slow convergence".

  So there are two defects, and they are independent:

    1. the ladder cannot reach the step lengths that work;
    2. it takes the last step it tried rather than the best it measured.

  Neither is a convergence criterion.  What an increment must satisfy to
  be accepted is untouched; this only chooses how far to step along a
  direction the solver has already computed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CalculiX.h"

/* alpha0 is the FIRST rung, and it is not 1: the caller computes a secant
   step sum1/(sum1-sum2) and clamps it into [floor, 0.80] before the ladder
   ever runs.  Overwriting it with 1 changes the search - measured, the
   trajectory diverges at increment 100 - so it is an argument here rather
   than an assumption. */

void lsladder_start(lsladder *l,double alpha0,double oldres,double flr,
                    double ratio,ITG ntrial,ITG legacy){

  l->alpha=((alpha0>0.)&&(alpha0<=1.))?alpha0:1.;
  l->best=-1.;
  l->bestres=0.;
  l->oldres=oldres;
  l->floor=(flr>0.)?flr:1.e-3;
  l->ratio=((ratio>0.)&&(ratio<1.))?ratio:0.5;
  l->ntrial=(ntrial>=1)?ntrial:1;
  l->itrial=0;
  l->contracted=0;
  l->legacy=(legacy!=0)?1:0;
}

/* Record the residual measured at the current alpha.

   Returns  1  accept: this trial contracted, stop here
            0  continue: alpha has been lowered to the next rung
           -1  exhausted: no rung contracted, take lsladder_final()      */

ITG lsladder_step(lsladder *l,double res){

  l->itrial++;

  if((l->best<0.)||(res<l->bestres)){
    l->best=l->alpha;
    l->bestres=res;
  }

  if(res<=l->oldres){
    l->contracted=1;
    return 1;
  }

  if(l->itrial>=l->ntrial) return -1;

  /* Already standing on the floor: descending again would re-test the
     same step length and burn the remaining rungs on an answer the
     search already has. */

  if(l->alpha<=l->floor*(1.+1.e-12)) return -1;

  /* The legacy ladder jumped straight to the floor on its last rung,
     which is why {1.0, 0.5, 0.1} has three rungs and not four.  Keep
     that shape when reproducing it. */

  if(l->itrial==l->ntrial-1){
    l->alpha=l->floor;
  }else{
    l->alpha*=l->ratio;
    if(l->alpha<l->floor) l->alpha=l->floor;
  }
  return 0;
}

/* The step length to actually take.  On a contraction that is the trial
   that contracted.  Otherwise it is the BEST rung measured - never the
   last one, unless the two coincide or the legacy shape is asked for. */

double lsladder_final(const lsladder *l){

  if(l->contracted!=0) return l->alpha;
  if(l->legacy!=0) return l->alpha;
  return (l->best>0.)?l->best:l->alpha;
}

/* ------------------------------------------------------------------ */
/* Regression test for the defect this file exists to fix.             */
/* ------------------------------------------------------------------ */

ITG lsladder_selftest(void){

  ITG nbad=0,i,r;
  lsladder l;
  double a[16];

  printf("[LSLADDER] self test\n");

  /* ---- A: the legacy ladder is exactly {1.0, 0.5, 0.1} ------------- */

  lsladder_start(&l,1.,1.,0.10,0.5,3,1);
  for(i=0;i<3;i++){a[i]=l.alpha;r=lsladder_step(&l,1.e30);}
  printf("   %-30s %.4f %.4f %.4f  %s\n","A legacy rungs",a[0],a[1],a[2],
         ((fabs(a[0]-1.)<1.e-12)&&(fabs(a[1]-0.5)<1.e-12)&&
          (fabs(a[2]-0.1)<1.e-12)&&(r==-1))?"ok":"FAIL");
  if(!((fabs(a[0]-1.)<1.e-12)&&(fabs(a[1]-0.5)<1.e-12)&&
       (fabs(a[2]-0.1)<1.e-12)&&(r==-1))) nbad++;

  /* ---- B: legacy takes the LAST rung, which is the defect ---------- */
  /* residuals 9, 3, 7 against oldres 1: nothing contracts, the best rung
     is 0.5 with 3, and the legacy rule hands back 0.1 anyway. */

  lsladder_start(&l,1.,1.,0.10,0.5,3,1);
  lsladder_step(&l,9.);lsladder_step(&l,3.);lsladder_step(&l,7.);
  printf("   %-30s best=%.4f taken=%.4f  %s\n","B legacy takes the last",
         l.best,lsladder_final(&l),
         ((fabs(l.best-0.5)<1.e-12)&&(fabs(lsladder_final(&l)-0.1)<1.e-12))?
         "ok (reproduces the defect)":"FAIL");
  if(!((fabs(l.best-0.5)<1.e-12)&&(fabs(lsladder_final(&l)-0.1)<1.e-12)))
    nbad++;

  /* ---- C: the fixed rule takes the BEST rung ----------------------- */

  lsladder_start(&l,1.,1.,1.e-3,0.5,8,0);
  lsladder_step(&l,9.);lsladder_step(&l,3.);lsladder_step(&l,7.);
  printf("   %-30s best=%.4f taken=%.4f  %s\n","C fixed takes the best",
         l.best,lsladder_final(&l),
         (fabs(lsladder_final(&l)-0.5)<1.e-12)?"ok":"FAIL");
  if(!(fabs(lsladder_final(&l)-0.5)<1.e-12)) nbad++;

  /* ---- D: the fixed ladder reaches below 0.1 ----------------------- */
  /* 46.5% of the measured activations had their optimum below the old
     floor, so reaching there is the whole point. */

  lsladder_start(&l,0.8,1.,1.e-3,0.5,8,0);
  for(i=0;i<8;i++){a[i]=l.alpha;if(lsladder_step(&l,1.e30)!=0) break;}
  printf("   %-30s rungs",("D fixed ladder depth"));
  for(r=0;r<8;r++) printf(" %.4g",a[r]);
  printf("  %s\n",(a[7]<0.1)?"ok":"FAIL");
  if(!(a[7]<0.1)) nbad++;

  /* ---- E: a contracting rung is accepted at once ------------------- */

  lsladder_start(&l,1.,1.,1.e-3,0.5,8,0);
  r=lsladder_step(&l,0.5);
  printf("   %-30s return=%" ITGFORMAT " alpha=%.4f  %s\n",
         "E contraction accepted",r,lsladder_final(&l),
         ((r==1)&&(fabs(lsladder_final(&l)-1.)<1.e-12))?"ok":"FAIL");
  if(!((r==1)&&(fabs(lsladder_final(&l)-1.)<1.e-12))) nbad++;

  /* ---- F: alpha is monotone and never below the floor -------------- */

  lsladder_start(&l,1.,1.,0.02,0.5,12,0);
  {
    double prev=2.;
    ITG ok=1;
    for(i=0;i<12;i++){
      if(!(l.alpha<prev)&&(i>0)) ok=0;
      if(l.alpha<0.02-1.e-15) ok=0;
      prev=l.alpha;
      if(lsladder_step(&l,1.e30)!=0) break;
    }
    printf("   %-30s %s\n","F monotone, respects the floor",
           ok?"ok":"FAIL");
    if(!ok) nbad++;
  }

  /* ---- G: the first rung is the caller's alpha0, not 1 -------------- */
  /* The caller clamps a secant step into [floor, 0.80]; forcing the ladder
     to start at 1 instead was measured to change the trajectory. */

  lsladder_start(&l,0.548186,1.,1.e-3,0.5,8,0);
  a[0]=l.alpha;
  lsladder_step(&l,1.e30);
  a[1]=l.alpha;
  printf("   %-30s first two rungs %.6f %.6f  %s\n","G alpha0 is honoured",
         a[0],a[1],
         ((fabs(a[0]-0.548186)<1.e-12)&&(fabs(a[1]-0.274093)<1.e-9))?
         "ok":"FAIL");
  if(!((fabs(a[0]-0.548186)<1.e-12)&&(fabs(a[1]-0.274093)<1.e-9))) nbad++;

  printf("[LSLADDER] self test %s (%" ITGFORMAT " failure(s))\n",
         (nbad==0)?"PASSED":"FAILED",nbad);
  return nbad;
}
