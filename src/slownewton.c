/*     CalculiX - damage/fracture extension                              */
/*     slownewton.c: how many Newton iterations this increment may have.  */

/* Why this module exists
   ----------------------
   This is a DECISION, and it was two file-statics and thirteen locals with
   no owner.  The stock controller starts estimating convergence at ir and
   cuts the increment at ic - normally sixteen - even when the residual and
   the displacement correction are both contracting monotonically.  A
   damaged tangent is approximate, so linear convergence is what one
   expects while D evolves, and the stock rule throws away increments that
   were going to converge.

   The extension is BOUNDED and it is conditional on softening actually
   being under way, which is why erosion_softening() exists: "there is a
   damage material in the model" is not a reason to spend forty iterations.

   Contract
   --------
     - slownewton_estimate() answers one question: at the contraction
       measured over the last step, how many iterations in total would this
       quantity need to reach its target?  It is a pure function of three
       doubles.  Invalid or non-contracting data returns a large sentinel,
       deliberately, so that bad data POSTPONES nothing;
     - slownewton_allow() is the decision.  It says yes only while the
       iteration is genuinely contracting, only within MAX_ITERS, and only
       within MAX_EXTRA beyond the estimate;
     - the stock checkconvergence() criteria and every one of its
       divergence paths are untouched.  This extends a budget; it does not
       relax a tolerance, and nothing here can turn a diverging iteration
       into an accepted one.                                            */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void slownewton_init(slownewton *s)
{
  memset(s,0,sizeof(*s));
  s->camprev1=1.e300;
  s->camprev2=1.e300;
  s->maxiters=DAMAGE_SLOW_NEWTON_MAX_ITERS;
}

/* Geometric estimate of the total iteration number required for value to
   reach target at the contraction measured over the latest Newton step.
   Invalid/non-contracting data deliberately returns a large sentinel so the
   stock cutback is not postponed. */
ITG slownewton_estimate(ITG iit,double value,
                                       double previous,double target)
{
  double ratio,nextra,iest;

  if((!isfinite(value))||(!isfinite(previous))||(!isfinite(target))||
     (value<0.)||(previous<=0.)||(target<=0.))
    return DAMAGE_SLOW_NEWTON_INVALID_EST;
  if(value<=target) return iit;
  if(value>=previous) return DAMAGE_SLOW_NEWTON_INVALID_EST;

  ratio=value/previous;
  if((ratio<=0.)||(ratio>=1.)) return DAMAGE_SLOW_NEWTON_INVALID_EST;
  nextra=ceil(log(target/value)/log(ratio));
  if((!isfinite(nextra))||(nextra<0.)||
     (nextra>(double)DAMAGE_SLOW_NEWTON_INVALID_EST-iit))
    return DAMAGE_SLOW_NEWTON_INVALID_EST;

  iest=(double)iit+nextra;
  if(iest>(double)DAMAGE_SLOW_NEWTON_INVALID_EST)
    return DAMAGE_SLOW_NEWTON_INVALID_EST;
  return (ITG)iest;
}

/* Decide whether it is worth paying for more PARDISO factorizations.  The
   direct stock force and correction targets are used.  Requiring two
   consecutive contractions rejects residual oscillation and one-step noise;
   the hard cap and remaining-iteration cap bound the cost. */
ITG slownewton_allow(ITG iit,const double *ram,
                                    const double *ram1,const double *ram2,
                                    const double *cam,const double *uam,
                                    double camprev1,double camprev2,
                                    const double *qa,const double *qam,
                                    const double *ctrl,ITG maxiters,
                                    ITG *iestres,ITG *iestcorr,
                                    ITG *iesttotal,double *rratio,
                                    double *cratio)
{
  ITG ip;
  double ea,c1,c2,targetres,targetcorr;

  *iestres=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *iestcorr=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *iesttotal=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *rratio=0.;
  *cratio=0.;

  if((iit<3)||(iit>=maxiters)) return 0;
  if((ram[0]<=0.)||(ram1[0]<=0.)||(ram2[0]<=0.)||
     (cam[0]<=0.)||(camprev1<=0.)||(camprev2<=0.)||(uam[0]<=0.))
    return 0;

  /* Strict two-step monotonicity.  ram2 still contains the genuine
     two-iterations-old residual here; checkconvergence() has not yet folded
     it into its running minimum. */
  if(!((ram[0]<ram1[0])&&(ram1[0]<ram2[0])&&
       (cam[0]<camprev1)&&(camprev1<camprev2))) return 0;

  *rratio=ram[0]/ram1[0];
  *cratio=cam[0]/camprev1;

  ea=ctrl[23];
  ip=(ITG)ctrl[2];
  if(qa[0]>ea*qam[0]){
    c1=(iit<=ip)?ctrl[18]:ctrl[22];
    c2=ctrl[19];
  }else{
    c1=ea;
    c2=ctrl[24];
  }

  targetres=c1*qam[0];
  targetcorr=c2*uam[0];
  *iestres=slownewton_estimate(iit,ram[0],ram1[0],targetres);
  *iestcorr=slownewton_estimate(iit,cam[0],camprev1,targetcorr);
  *iesttotal=(*iestres>*iestcorr)?*iestres:*iestcorr;

  if(*iesttotal>maxiters) return 0;
  if(*iesttotal-iit>DAMAGE_SLOW_NEWTON_MAX_EXTRA) return 0;
  return 1;
}
