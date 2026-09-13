/*     CalculiX - damage/fracture extension                              */
/*     converge.c: the numbers the convergence judgement is made from.   */

/* Why this module exists
   ----------------------
   Ask this tree "what counts as converged here" and you arrive at five
   places (handover/10-CONVERGENCE.md section 1).  This file is the first of
   them to get an owner: the REDUCTION - what ram, ram1, ram2, uam and qam
   are, and which degrees of freedom are allowed to contribute to them.

   That is a real responsibility and it had no home.  The clearest symptom:
   CCX_DAMAGE_AUTOSPC_FORCE excludes dofs from the force residual by a
   `continue` in the middle of the loop that computes it.  There is no
   object whose answer that is, so it lives where the loop happened to be -
   and the thing it modifies has no name, which is why it could be measured
   for two days before anyone established that on the target deck it changes
   nothing at all (research/14-THE-TRAP.md).

   Contract
   --------
     - one call per Newton iteration, after the solve and before the
       verdict;
     - it reads the solution vector and the load measures and writes the
       norms; it decides nothing and prints nothing;
     - what it EXCLUDED is part of its state, not a local.  The declaration
       of CCX_DAMAGE_AUTOSPC_FORCE promises the excluded residual is
       reported next to the criterion it was excluded from, and a promise
       like that belongs to the thing that does the excluding;
     - the module refuses to arm if its self test fails, the discipline
       lsladder.c set.

   What this is NOT
   ----------------
   Not the verdict (checkconvergence.c still owns the boolean), not the
   increment control, not the printing.  Those are steps B and C and the
   Increment and Monitor objects; 10-CONVERGENCE.md sections 3 and 4 say
   why they are held back rather than done here.

   Provenance
   ----------
   Transcribed from nonlingeo.c lines 12509-12601 and the two `ram < 1e-6`
   cut-offs that followed them, unchanged.  Bit-identity is the standard for
   this step and it is checked, not assumed: gate and target deck.        */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void converge_init(converge *c,double qam_floor,ITG mask_force,
                   const ITG *mask,ITG mask_nk)
{
  c->qam_floor=qam_floor;
  c->mask_force=mask_force;
  c->mask=mask;
  c->mask_nk=mask_nk;
  c->qam_peak=0.;
  c->excl_max=0.;
  c->excl_node=0;
  c->excl_count=0;
}

/* One Newton iteration's worth of norms.
 *
 * ram[0],ram[1]  largest |residual| over the mechanical / thermal dofs
 * ram[2],ram[3]  the equation it was found at, carried as k+0.5
 * ram[4],ram[5]  the face-to-face contact divergence pair
 * ram1,ram2      the two previous iterations, for the divergence test
 * uam            high-water mark of the correction
 * qam            running average of the load measure, optionally floored
 *
 * The order is load-bearing and is the order it had inline: the mortar pair
 * is formed from the UNCUT ram[0], and only then is ram cut off at 1e-6. */
void converge_norms(converge *c,const double *b,const ITG *neq,
                    const ITG *nactdofinv,ITG mt,ITG ithermal,ITG mortar,
                    ITG ne,ITG ne0,ITG neold,
                    const double *qa,const double *qamold,ITG jnz,
                    double qau,double ea,
                    double *ram,double *ram1,double *ram2,
                    const double *cam,double *uam,double *qam)
{
  ITG k;
  double err;

  /* the correction high-water mark and the load reference */

  if(ithermal!=2){
    if(cam[0]>uam[0]){
      uam[0]=cam[0];}
    if(qau<1.e-10){
      if(qa[0]>ea*qam[0]){
        qam[0]=(qamold[0]*jnz+qa[0])/(jnz+1);}
      else {
        qam[0]=qamold[0];}
    }
    /* CCX_DAMAGE_QAM_FLOOR - DIAGNOSTIC, default off.

       The force criterion is RELATIVE: the verdict tests ram[0] against
       c1[0]*qam[0], and qam[0] is a running average of the internal force
       over increments.  On a specimen that is unloading as it breaks,
       qam[0] follows the load down and the absolute tolerance collapses
       with it.  Measured on run_s0_fine_grad14 at its wall: qa=0.001552,
       qam=0.001604, residual 0.005474 - a miss by 3.4x - while the
       specimen was still carrying 2.33 N, 22.7% of its peak, and the
       residual had just fallen 14x in one iteration.  The run died on a
       vanishing reference, not on a growing residual.

       Flooring qam at a fraction of the largest value it ever reached keeps
       the tolerance tied to the load the specimen ONCE carried.

       THIS CHANGES THE CONVERGENCE CRITERION AND THEREFORE THE ANSWER.  It
       is a diagnostic for whether a given wall is criterial or physical.  A
       run that goes further with it is NOT thereby a success - that has to
       be shown on the physics (E-108) - and adopting it needs the full
       verify + ladder gate. */
    if(c->qam_floor>0.){
      if(qam[0]>c->qam_peak){c->qam_peak=qam[0];}
      if(qam[0]<c->qam_floor*c->qam_peak){
        qam[0]=c->qam_floor*c->qam_peak;}
    }
  }
  if(ithermal>1){
    if(cam[1]>uam[1]){
      uam[1]=cam[1];}
    if(qau<1.e-10){
      if(qa[1]>ea*qam[1]){
        qam[1]=(qamold[1]*jnz+qa[1])/(jnz+1);}
      else {
        qam[1]=qamold[1];}
    }
  }

  /* the two previous residuals, then the new one */

  for(k=0;k<2;++k){
    ram2[k]=ram1[k];
    ram1[k]=ram[k];
    ram[k]=0.;
  }

  if(ithermal!=2){
    c->excl_max=0.;c->excl_node=0;c->excl_count=0;
    for(k=0;k<neq[0];++k){
      err=fabs(b[k]);
      /* CCX_DAMAGE_AUTOSPC_FORCE: a node damstate.c has judged to have lost
         its load path does not get to veto the increment through the force
         norm.  Excluded, not hidden - see excl_* below. */
      if((c->mask_force!=0)&&(c->mask!=NULL)&&(nactdofinv!=NULL)){
        ITG spcnd=nactdofinv[k]/mt;
        if((spcnd>=0)&&(spcnd<c->mask_nk)&&(c->mask[spcnd]!=0)){
          c->excl_count++;
          if(err>c->excl_max){
            c->excl_max=err;c->excl_node=spcnd+1;}
          continue;
        }
      }
      if(err>ram[0]){
        ram[0]=err;
        ram[2]=k+0.5;}
    }
  }
  if(ithermal>1){
    for(k=neq[0];k<neq[1];++k){
      err=fabs(b[k]);
      if(err>ram[1]){
        ram[1]=err;
        ram[3]=k+0.5;}
    }
  }

  /* Divergence criteria for face-to-face penalty is different */

  if(mortar==1){
    for(k=4;k<6;++k){
      ram2[k]=ram1[k];
      ram1[k]=ram[k];
    }
    ram[4]=ram[0]+ram1[0];
    ram[5]=(ne-ne0)-(neold-ne0)+0.5;
  }

  /* next line is inserted to cope with stress-less temperature
     calculations.  AFTER the mortar pair, which is formed from the uncut
     value - that order was in the original and is not incidental. */

  if(ithermal!=2){
    if(ram[0]<1.e-6){
      ram[0]=0.;}
  }
  if(ithermal>1){
    if(ram[1]<1.e-6){
      ram[1]=0.;}
  }
}

/* Present the norms.  Separate from computing them, because "produce a
 * diagnostic" and "present it" are two jobs and this file only does the
 * second one here under protest: the printing belongs to Monitor, and this
 * is a way-station so that nonlingeo() has two named calls where it had
 * ninety lines of arithmetic interleaved with output.
 *
 * Byte-for-byte the text that was inline, including the two blank lines
 * after the correction and the position of fflush.  The AUTOSPC-FORCE line
 * is printed from the object's own excl_* state, which is the whole point:
 * the declaration of that switch promises the excluded residual appears
 * next to the criterion it was excluded from, and now one object both
 * excludes and reports.
 *
 * ran is ctrl[18], the coefficient the force criterion uses, so the printed
 * tolerance is the one the verdict will actually apply. */
void converge_report(const converge *c,const ITG *nactdofinv,ITG mt,
                     ITG ithermal,double ran,
                     const double *qa,const double *qam,const double *ram,
                     const double *cam,const double *uam)
{
  ITG inode,idir;

  if(ithermal!=2){
    printf(" average force= %f\n",qa[0]);
    printf(" time avg. forc= %f\n",qam[0]);
    if((c->mask_force!=0)&&(c->excl_count>0)){
      /* Auditable by construction: the excluded peak is printed next to
         the criterion it was excluded from, so a masked residual that
         starts to grow is visible in the same place ram[0] is read. */
      printf("[DAMAGE AUTOSPC-FORCE] excluded %" ITGFORMAT " dof(s) on "
             "AUTOSPC-masked nodes from ram[0]; largest excluded "
             "residual %.6e at node %" ITGFORMAT " (ram[0]=%.6e, "
             "tolerance %.6e = %.4f x qam)%s",
             c->excl_count,c->excl_max,c->excl_node,
             ram[0],ran*qam[0],
             (qam[0]>0.)?c->excl_max/qam[0]:0.,"\n");
    }
    if((ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1==0){
      printf(" largest residual force= %f\n",
             ram[0]);
    }else{
      inode=(ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1;
      idir=nactdofinv[(ITG)ram[2]]-mt*(inode-1);
      printf(" largest residual force= %f in node %" ITGFORMAT
             " and dof %" ITGFORMAT "\n",
             ram[0],inode,idir);
    }
    printf(" largest increment of disp= %e\n",uam[0]);
    if((ITG)cam[3]==0){
      printf(" largest correction to disp= %e\n\n",
             cam[0]);
    }else{
      inode=(ITG)((double)nactdofinv[(ITG)cam[3]]/mt)+1;
      idir=nactdofinv[(ITG)cam[3]]-mt*(inode-1);
      printf(" largest correction to disp= %e in node %" ITGFORMAT
             " and dof %" ITGFORMAT "\n\n",cam[0],inode,idir);
    }
  }
  if(ithermal>1){
    printf(" average flux= %f\n",qa[1]);
    printf(" time avg. flux= %f\n",qam[1]);
    if((ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1==0){
      printf(" largest residual flux= %f\n",
             ram[1]);
    }else{
      inode=(ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1;
      idir=nactdofinv[(ITG)ram[3]]-mt*(inode-1);
      printf(" largest residual flux= %f in node %" ITGFORMAT
             " and dof %" ITGFORMAT "\n",ram[1],inode,idir);
    }
    printf(" largest increment of temp= %e\n",uam[1]);
    if((ITG)cam[4]==0){
      printf(" largest correction to temp= %e\n\n",
             cam[1]);
    }else{
      inode=(ITG)((double)nactdofinv[(ITG)cam[4]]/mt)+1;
      idir=nactdofinv[(ITG)cam[4]]-mt*(inode-1);
      printf(" largest correction to temp= %e in node %" ITGFORMAT
             " and dof %" ITGFORMAT "\n\n",cam[1],inode,idir);
    }
  }
  fflush(stdout);
}

/* ------------------------------------------------- what counts as converged

   Step B of handover/10-CONVERGENCE.md.  Transcribed from
   checkconvergence.c lines 128-205, clause for clause, in the same order.

   ONE DELIBERATE DEVIATION from the sketch in the design note, and it is
   the reason this is worth doing: the sketch short-circuited like the
   original, so a leaf past the first failure was never evaluated and could
   not be reported.  Every leaf here is a comparison of two doubles with no
   side effect, so evaluating all of them cannot change the combination -
   and then the table is complete even for the iteration that failed on the
   first clause.  "Which clause is holding this increment back" is the
   question this object exists to answer; a table with holes in it does not
   answer it.

   The boolean is unchanged.  That is the standard for this step and it is
   checked on the gate, on stdout and on the target deck, not asserted. */

static ITG leaf(cvg_verdict *v,const char *name,ITG ok,
                double value,double thresh)
{
  if(v->nclause<CVG_MAXCLAUSE){
    v->clause[v->nclause].name=name;
    v->clause[v->nclause].value=value;
    v->clause[v->nclause].thresh=thresh;
    v->clause[v->nclause].status=ok?CVG_PASS:CVG_FAIL;
    v->nclause++;
  }
  return ok;
}

/* record the first failing leaf of a top-level AND chain */
static void blocked_by(cvg_verdict *v,const char *name,double value,
                       double thresh)
{
  if(v->blocker==NULL){
    v->blocker=name;
    v->blocker_value=value;
    v->blocker_thresh=thresh;
  }
}

void cvg_tol_from_ctrl(cvg_tol *t,const double *ctrl)
{
  t->ip   =(ITG)ctrl[2];
  t->ran  =ctrl[18];
  t->can  =ctrl[19];
  t->rap  =ctrl[22];
  t->ea   =ctrl[23];
  t->cae  =ctrl[24];
  t->ral  =ctrl[25];
  t->cetol=ctrl[39];
}

ITG converge_verdict(cvg_verdict *v,const cvg_tol *t,ITG ithermal,ITG iit,
                     ITG nmethod,ITG iflagact,ITG ntg,double deltmx,
                     const double *ram,const double *ram1,double *ram2,
                     const double *cam,const double *uam,
                     const double *qa,const double *qam,
                     double *c1,double *c2)
{
  ITG a1,a2,a3,a4,b1,b2,b3,b4,b5,b6,mech,ther;

  v->converged=0; v->nclause=0; v->blocker=NULL;
  v->blocker_value=0.; v->blocker_thresh=0.;

  /* which tolerance applies.  This is part of the verdict and not a
     preamble to it: whether the increment is judged at the loose working
     tolerance or the tight end-of-load one is decided here, from whether
     the load increment is still large. */

  if(ithermal!=2){
    if(qa[0]>t->ea*qam[0]){
      if(iit<=t->ip){c1[0]=t->ran;}
      else{c1[0]=t->rap;}
      c2[0]=t->can;
    }
    else{
      c1[0]=t->ea;
      c2[0]=t->cae;
    }
    if(ram1[0]<ram2[0]){ram2[0]=ram1[0];}
  }
  if(ithermal>1){
    if(qa[1]>t->ea*qam[1]){
      if(iit<=t->ip){c1[1]=t->ran;}
      else{c1[1]=t->rap;}
      c2[1]=t->can;
    }
    else{
      c1[1]=t->ea;
      c2[1]=t->cae;
    }
    if(ram1[1]<ram2[1]){ram2[1]=ram1[1];}
  }

  /* mechanical */

  if(ithermal<2){
    a1=leaf(v,"IterationsAtLeast2",iit>1,(double)iit,1.);
    a2=leaf(v,"ForceResidual",ram[0]<=c1[0]*qam[0],ram[0],c1[0]*qam[0]);
    a3=leaf(v,"ContactSetStable",iflagact==0,(double)iflagact,0.);
    a4=leaf(v,"CreepTolerance",(nmethod!=-1)||(qa[3]<=t->cetol),
            qa[3],t->cetol);
    b1=leaf(v,"SolutionChange",cam[0]<=c2[0]*uam[0],cam[0],c2[0]*uam[0]);
    b2=leaf(v,"Projected",ram[0]*cam[0]<c2[0]*uam[0]*ram2[0],
            ram[0]*cam[0],c2[0]*uam[0]*ram2[0]);
    b3=leaf(v,"ResidualWellBelow",ram[0]<=t->ral*qam[0],ram[0],t->ral*qam[0]);
    b4=leaf(v,"LoadIncrementSmall",qa[0]<=t->ea*qam[0],qa[0],t->ea*qam[0]);
    b5=leaf(v,"NoGasNetwork",ntg==0,(double)ntg,0.);
    b6=leaf(v,"CorrectionFloor",cam[0]<1.e-8,cam[0],1.e-8);

    if(a1&&a2&&a3&&a4&&(b1||(((b2||b3||b4))&&b5)||b6)) v->converged=1;

    if(!a1)      blocked_by(v,"IterationsAtLeast2",(double)iit,1.);
    else if(!a2) blocked_by(v,"ForceResidual",ram[0],c1[0]*qam[0]);
    else if(!a3) blocked_by(v,"ContactSetStable",(double)iflagact,0.);
    else if(!a4) blocked_by(v,"CreepTolerance",qa[3],t->cetol);
    else if(!v->converged)
      blocked_by(v,"SolutionChange",cam[0],c2[0]*uam[0]);
  }

  /* thermal */

  if(ithermal==2){
    a2=leaf(v,"FluxResidual",ram[1]<=c1[1]*qam[1],ram[1],c1[1]*qam[1]);
    a3=leaf(v,"TempChangeLimit",cam[2]<deltmx,cam[2],deltmx);
    b1=leaf(v,"TempChange",cam[1]<=c2[1]*uam[1],cam[1],c2[1]*uam[1]);
    b2=leaf(v,"TempProjected",ram[1]*cam[1]<c2[1]*uam[1]*ram2[1],
            ram[1]*cam[1],c2[1]*uam[1]*ram2[1]);
    /* change 25.11.2017 in the stock source: this leaf alone carries an
       extra iit>1, which the mechanical and thermomechanical forms of the
       same leaf do not.  Transcribed as found; see handover/05-DEBT.md. */
    b3=leaf(v,"TempResidualWellBelow",(ram[1]<=t->ral*qam[1])&&(iit>1),
            ram[1],t->ral*qam[1]);
    b4=leaf(v,"TempLoadIncrementSmall",qa[1]<=t->ea*qam[1],qa[1],
            t->ea*qam[1]);
    b5=leaf(v,"NoGasNetwork",ntg==0,(double)ntg,0.);
    b6=leaf(v,"TempCorrectionFloor",cam[1]<1.e-8,cam[1],1.e-8);

    if(a2&&a3&&(b1||(((b2||b3||b4))&&b5)||b6)) v->converged=1;

    if(!a2)      blocked_by(v,"FluxResidual",ram[1],c1[1]*qam[1]);
    else if(!a3) blocked_by(v,"TempChangeLimit",cam[2],deltmx);
    else if(!v->converged)
      blocked_by(v,"TempChange",cam[1],c2[1]*uam[1]);
  }

  /* thermomechanical */

  if(ithermal==3){
    a1=leaf(v,"IterationsAtLeast2",iit>1,(double)iit,1.);
    a2=leaf(v,"ForceResidual",ram[0]<=c1[0]*qam[0],ram[0],c1[0]*qam[0]);
    /* no ContactSetStable leaf here: the stock thermomechanical form omits
       iflagact, which the mechanical form requires.  See 05-DEBT.md. */
    a4=leaf(v,"CreepTolerance",(nmethod!=-1)||(qa[3]<=t->cetol),
            qa[3],t->cetol);
    b1=leaf(v,"SolutionChange",cam[0]<=c2[0]*uam[0],cam[0],c2[0]*uam[0]);
    b2=leaf(v,"Projected",ram[0]*cam[0]<c2[0]*uam[0]*ram2[0],
            ram[0]*cam[0],c2[0]*uam[0]*ram2[0]);
    b3=leaf(v,"ResidualWellBelow",ram[0]<=t->ral*qam[0],ram[0],t->ral*qam[0]);
    b4=leaf(v,"LoadIncrementSmall",qa[0]<=t->ea*qam[0],qa[0],t->ea*qam[0]);
    b5=leaf(v,"NoGasNetwork",ntg==0,(double)ntg,0.);
    b6=leaf(v,"CorrectionFloor",cam[0]<1.e-8,cam[0],1.e-8);
    mech=(a1&&a2&&a4&&(b1||(((b2||b3||b4))&&b5)||b6));

    a3=leaf(v,"FluxResidual",ram[1]<=c1[1]*qam[1],ram[1],c1[1]*qam[1]);
    ther=leaf(v,"TempChangeLimit",cam[2]<deltmx,cam[2],deltmx);
    {
      ITG d1,d2,d3,d4,d6;
      d1=leaf(v,"TempChange",cam[1]<=c2[1]*uam[1],cam[1],c2[1]*uam[1]);
      d2=leaf(v,"TempProjected",ram[1]*cam[1]<c2[1]*uam[1]*ram2[1],
              ram[1]*cam[1],c2[1]*uam[1]*ram2[1]);
      /* and here the same leaf has NO iit>1 - the asymmetry with the
         pure-thermal form above, visible only once both have names. */
      d3=leaf(v,"TempResidualWellBelow",ram[1]<=t->ral*qam[1],
              ram[1],t->ral*qam[1]);
      d4=leaf(v,"TempLoadIncrementSmall",qa[1]<=t->ea*qam[1],qa[1],
              t->ea*qam[1]);
      d6=leaf(v,"TempCorrectionFloor",cam[1]<1.e-8,cam[1],1.e-8);
      ther=(a3&&ther&&(d1||(((d2||d3||d4))&&b5)||d6));
    }

    if(mech&&ther) v->converged=1;

    if(!a1)       blocked_by(v,"IterationsAtLeast2",(double)iit,1.);
    else if(!a2)  blocked_by(v,"ForceResidual",ram[0],c1[0]*qam[0]);
    else if(!a4)  blocked_by(v,"CreepTolerance",qa[3],t->cetol);
    else if(!mech)blocked_by(v,"SolutionChange",cam[0],c2[0]*uam[0]);
    else if(!a3)  blocked_by(v,"FluxResidual",ram[1],c1[1]*qam[1]);
    else if(!ther)blocked_by(v,"TempChangeLimit",cam[2],deltmx);
  }

  return v->converged;
}

void converge_verdict_print(const cvg_verdict *v)
{
  ITG k;
  printf("[CONVERGE VERDICT] %s\n",v->converged?"CONVERGED":"not yet");
  for(k=0;k<v->nclause;k++){
    printf("   %-24s %-13s value=%-14.6e thresh=%-14.6e\n",
           v->clause[k].name,
           v->clause[k].status==CVG_PASS?"pass":
           (v->clause[k].status==CVG_FAIL?"FAIL":"-"),
           v->clause[k].value,v->clause[k].thresh);
  }
  if((v->converged==0)&&(v->blocker!=NULL)){
    printf("   held back by: %s (value=%.6e vs %.6e)\n",
           v->blocker,v->blocker_value,v->blocker_thresh);
  }
  fflush(stdout);
}

/* --------------------------------------------------- why the run stopped

   Step C.  Nothing here changes a decision: the four stop sites keep their
   stock message and their stock exit.  What is added is the one line that
   says WHICH of them fired and what the last verdict was blocked on, so
   that "increment size smaller than minimum" stops being a symptom shared
   by three unrelated situations.

   The verdict is passed in rather than remembered in a file-scope
   variable: three of the four sites are in the same function as the
   verdict, and the fourth (checkdivergence.c) honestly has none - it
   passes NULL and the report says so instead of inventing one.        */

const char *converge_reason_name(cvg_reason r)
{
  switch(r){
  case CVG_DIVERGED_MINSTEP_EXTERNAL:      return "DIVERGED_MINSTEP_EXTERNAL";
  case CVG_DIVERGED_MINSTEP_TOO_SLOW:      return "DIVERGED_MINSTEP_TOO_SLOW";
  case CVG_DIVERGED_MINSTEP_ON_DIVERGENCE: return "DIVERGED_MINSTEP_ON_DIVERGENCE";
  case CVG_DIVERGED_MINSTEP_AFTER_CONV:    return "DIVERGED_MINSTEP_AFTER_CONVERGENCE";
  case CVG_ITERATING:                      return "ITERATING";
  case CVG_CONVERGED_CRITERIA:             return "CONVERGED_CRITERIA";
  }
  return "UNKNOWN";
}

const char *converge_reason_explain(cvg_reason r)
{
  switch(r){
  case CVG_DIVERGED_MINSTEP_EXTERNAL:
    return "divergence was detected outside checkconvergence and the "
           "cut-back increment fell below tmin";
  case CVG_DIVERGED_MINSTEP_TOO_SLOW:
    return "the increment was converging too slowly - the estimated "
           "iteration count exceeded the limit - and the cut-back "
           "increment fell below tmin";
  case CVG_DIVERGED_MINSTEP_ON_DIVERGENCE:
    return "the residual diverged and the cut-back increment fell below "
           "tmin";
  case CVG_DIVERGED_MINSTEP_AFTER_CONV:
    return "the increment CONVERGED, but it took enough iterations that "
           "the next increment was decreased below tmin - the run stops "
           "on a step size, not on a failure to converge";
  case CVG_ITERATING:        return "still iterating";
  case CVG_CONVERGED_CRITERIA:
    return "every clause of the convergence criterion is satisfied";
  }
  return "unknown reason";
}

void converge_stop_report(cvg_reason r,const cvg_verdict *v)
{
  printf("\n[CONVERGE STOP] %s\n",converge_reason_name(r));
  printf("   %s\n",converge_reason_explain(r));
  if(v==NULL){
    printf("   last verdict: not available at this site\n");
  }else if(v->converged){
    printf("   last verdict: converged\n");
  }else if(v->blocker!=NULL){
    printf("   last verdict: held back by %s (%.6e vs %.6e)\n",
           v->blocker,v->blocker_value,v->blocker_thresh);
  }else{
    printf("   last verdict: not yet, no clause recorded\n");
  }
  fflush(stdout);
}

/* ---------------------------------------------------------------- tests */

static ITG cvg_chk(const char *name,double got,double want,double tol,
                   ITG *nbad)
{
  ITG ok=(fabs(got-want)<=tol);
  /* %.17g, not %g: the ordering check below differs from its expectation
     by 1e-7 in a number of order 5, and a report that prints "got=5 want=5
     FAIL" reads as a broken test rather than a caught defect. */
  printf("   %-36s got=%-22.17g want=%-22.17g %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

static ITG cvg_chki(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-36s got=%-22" ITGFORMAT " want=%-22" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

static ITG cvg_chks(const char *name,const char *got,const char *want,
                    ITG *nbad)
{
  ITG ok=((got!=NULL)&&(strcmp(got,want)==0));
  printf("   %-36s got=%-22s want=%-22s %s%s",
         name,got?got:"(none)",want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

/* A verdict case: everything the decision reads, so a test reads as the
   situation it describes rather than as nine array writes. */
static ITG cvg_case(cvg_verdict *v,const cvg_tol *t,ITG ithermal,ITG iit,
                    ITG iflagact,ITG ntg,
                    double ram0,double cam0,double uam0,double qa0,
                    double qam0,double ram2_0)
{
  double ram[6],ram1[6],ram2[6],cam[5],uam[2],qa[4],qam[2],c1[2],c2[2];
  ITG k;
  for(k=0;k<6;k++){ram[k]=0.;ram1[k]=1.e30;ram2[k]=1.e30;}
  for(k=0;k<5;k++) cam[k]=0.;
  for(k=0;k<4;k++) qa[k]=0.;
  uam[0]=uam[1]=0.; qam[0]=qam[1]=1.; c1[0]=c1[1]=0.; c2[0]=c2[1]=0.;
  ram[0]=ram0; ram2[0]=ram2_0; cam[0]=cam0; uam[0]=uam0;
  qa[0]=qa0; qam[0]=qam0;
  return converge_verdict(v,t,ithermal,iit,1,iflagact,ntg,1.e30,
                          ram,ram1,ram2,cam,uam,qa,qam,c1,c2);
}

ITG converge_selftest(void)
{
  ITG nbad=0,neq[2],nactdofinv[24],mt=4,k;
  double b[8],qa[2],qamold[2],ram[6],ram1[6],ram2[6],cam[2],uam[2],qam[2];
  ITG mask[8];
  converge c;

  printf("[CONVERGE] self test%s","\n");

  /* eight mechanical equations on nodes 0..5, mt=4 so node = dof/mt */
  neq[0]=8; neq[1]=8;
  for(k=0;k<8;k++) nactdofinv[k]=mt*(k/2)+1;   /* two dofs per node */
  for(k=0;k<8;k++) b[k]=0.;
  b[1]=-3.0;      /* node 0 */
  b[4]= 7.0;      /* node 2 - the largest                              */
  b[6]= 2.0;      /* node 3 */
  qa[0]=10.; qa[1]=0.;
  qamold[0]=4.; qamold[1]=0.;
  for(k=0;k<6;k++){ram[k]=0.;ram1[k]=0.;ram2[k]=0.;}
  cam[0]=0.5; cam[1]=0.;
  uam[0]=0.2; uam[1]=0.;
  qam[0]=4.;  qam[1]=0.;

  /* A: no mask.  The reduction finds |b| max and where it was. */
  converge_init(&c,0.,0,NULL,0);
  converge_norms(&c,b,neq,nactdofinv,mt,1,0,0,0,0,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  cvg_chk("A largest residual",ram[0],7.0,0.,&nbad);
  cvg_chk("A found at equation 4",ram[2],4.5,0.,&nbad);
  cvg_chk("A uam takes the correction",uam[0],0.5,0.,&nbad);
  cvg_chk("A qam running average",qam[0],(4.*1+10.)/(1+1),1.e-14,&nbad);
  cvg_chk("A nothing excluded",(double)c.excl_count,0.,0.,&nbad);

  /* B: mask node 2 - the one carrying the largest residual. */
  for(k=0;k<8;k++) mask[k]=0;
  mask[2]=1;
  qam[0]=4.;uam[0]=0.2;
  converge_init(&c,0.,1,mask,8);
  converge_norms(&c,b,neq,nactdofinv,mt,1,0,0,0,0,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  cvg_chk("B masked node does not veto",ram[0],3.0,0.,&nbad);
  cvg_chk("B excluded peak reported",c.excl_max,7.0,0.,&nbad);
  cvg_chk("B excluded node is 1-based",(double)c.excl_node,3.,0.,&nbad);
  cvg_chk("B excluded dof count",(double)c.excl_count,2.,0.,&nbad);

  /* C: the mask is inert unless AUTOSPC_FORCE armed it */
  qam[0]=4.;
  converge_init(&c,0.,0,mask,8);
  converge_norms(&c,b,neq,nactdofinv,mt,1,0,0,0,0,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  cvg_chk("C mask off: full residual",ram[0],7.0,0.,&nbad);

  /* D: the history shifts by exactly one iteration */
  cvg_chk("D ram1 is the previous ram",ram1[0],3.0,0.,&nbad);
  cvg_chk("D ram2 is the one before",ram2[0],7.0,0.,&nbad);

  /* E: the 1e-6 cut-off, and that it is applied AFTER the mortar pair */
  for(k=0;k<8;k++) b[k]=0.;
  b[3]=1.e-7;
  for(k=0;k<6;k++){ram[k]=0.;ram1[k]=0.;ram2[k]=0.;}
  ram[0]=5.;                 /* shifts into ram1[0] inside the call       */
  qam[0]=4.;
  converge_init(&c,0.,0,NULL,0);
  converge_norms(&c,b,neq,nactdofinv,mt,1,1,12,10,10,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  cvg_chk("E tiny residual is cut to zero",ram[0],0.,0.,&nbad);
  /* ram[4]=ram[0]+ram1[0] is formed BEFORE the cut-off, so it must carry
     the 1e-7 and not a bare 5.0.  This is the one ordering in the block
     that is easy to get wrong on a move and impossible to see afterwards. */
  cvg_chk("E mortar pair used the UNCUT value",ram[4],5.+1.e-7,1.e-12,&nbad);
  cvg_chk("E mortar element count",ram[5],2.5,0.,&nbad);

  /* F: the qam floor holds the reference at a fraction of its peak */
  qam[0]=100.;qamold[0]=100.;qa[0]=1.e-6;
  converge_init(&c,0.25,0,NULL,0);
  converge_norms(&c,b,neq,nactdofinv,mt,1,0,0,0,0,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  qamold[0]=1.;
  converge_norms(&c,b,neq,nactdofinv,mt,1,0,0,0,0,qa,qamold,1,0.,1.e-2,
                 ram,ram1,ram2,cam,uam,qam);
  cvg_chk("F qam floored at 0.25 of its peak",qam[0],25.,1.e-12,&nbad);

  /* ----------------------------------------------------------- reason */
  {
    cvg_reason rs[6]={CVG_DIVERGED_MINSTEP_EXTERNAL,
                      CVG_DIVERGED_MINSTEP_TOO_SLOW,
                      CVG_DIVERGED_MINSTEP_ON_DIVERGENCE,
                      CVG_DIVERGED_MINSTEP_AFTER_CONV,
                      CVG_ITERATING,CVG_CONVERGED_CRITERIA};
    ITG a,bq,distinct=1;

    /* every reason has a name and a sentence, and no two share either -
       the whole complaint against rc=201 is that one name covered four
       sites, so a duplicate here would reproduce the defect */
    for(a=0;a<6;a++){
      for(bq=a+1;bq<6;bq++){
        if(strcmp(converge_reason_name(rs[a]),
                  converge_reason_name(rs[bq]))==0) distinct=0;
        if(strcmp(converge_reason_explain(rs[a]),
                  converge_reason_explain(rs[bq]))==0) distinct=0;
      }
      if(strcmp(converge_reason_name(rs[a]),"UNKNOWN")==0) distinct=0;
    }
    cvg_chki("all six reasons named and distinct",distinct,1,&nbad);

    /* PETSc's sign convention: negative diverged, zero iterating,
       positive converged */
    cvg_chki("diverged reasons are negative",
             (CVG_DIVERGED_MINSTEP_EXTERNAL<0)&&
             (CVG_DIVERGED_MINSTEP_TOO_SLOW<0)&&
             (CVG_DIVERGED_MINSTEP_ON_DIVERGENCE<0)&&
             (CVG_DIVERGED_MINSTEP_AFTER_CONV<0),1,&nbad);
    cvg_chki("iterating is zero",CVG_ITERATING==0,1,&nbad);
    cvg_chki("converged is positive",CVG_CONVERGED_CRITERIA>0,1,&nbad);

    /* the one that is easiest to get wrong when reading the log: this
       stop follows a CONVERGED increment and is about step size */
    cvg_chks("after-convergence stop names convergence",
             strstr(converge_reason_explain(CVG_DIVERGED_MINSTEP_AFTER_CONV),
                    "CONVERGED")?"yes":"no","yes",&nbad);
  }

  /* ---------------------------------------------------------- verdict */
  {
    cvg_verdict v; cvg_tol t;
    double ctrl[60],ram[6],ram1[6],ram2[6],cam[5],uam[2],qa[4],qam[2];
    double c1[2],c2[2];
    ITG j,got;

    for(j=0;j<60;j++) ctrl[j]=0.;
    ctrl[2]=4.;        /* ip    */
    ctrl[18]=0.005;    /* ran   - the working force tolerance   */
    ctrl[19]=0.01;     /* can   */
    ctrl[22]=0.02;     /* rap   */
    ctrl[23]=1.e-5;    /* ea    */
    ctrl[24]=1.e-3;    /* cae   */
    ctrl[25]=1.e-4;    /* ral   */
    ctrl[39]=1.e30;    /* cetol */
    cvg_tol_from_ctrl(&t,ctrl);
    cvg_chk("tol.ran from ctrl[18]",t.ran,0.005,0.,&nbad);
    cvg_chk("tol.ral from ctrl[25]",t.ral,1.e-4,0.,&nbad);

    /* the ordinary converged iteration: residual under c1*qam and the
       solution change under c2*uam */
    got=cvg_case(&v,&t,0,3,0,0, 1.e-3, 1.e-4, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech converged",got,1,&nbad);

    /* first iteration never converges, whatever the numbers say */
    got=cvg_case(&v,&t,0,1,0,0, 1.e-3, 1.e-4, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech iit=1 blocked",got,0,&nbad);
    cvg_chks("  blocker",v.blocker,"IterationsAtLeast2",&nbad);

    /* residual too large - and the report names the clause, which is the
       whole reason this object exists */
    got=cvg_case(&v,&t,0,3,0,0, 1.0, 1.e-4, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech residual blocked",got,0,&nbad);
    cvg_chks("  blocker",v.blocker,"ForceResidual",&nbad);

    /* a changing contact set blocks it even with both norms converged */
    got=cvg_case(&v,&t,0,3,1,0, 1.e-3, 1.e-4, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech contact-change blocked",got,0,&nbad);
    cvg_chks("  blocker",v.blocker,"ContactSetStable",&nbad);

    /* solution change too big, but the residual is far below ral*qam:
       the OR branch carries it, and only because there is no gas network */
    got=cvg_case(&v,&t,0,3,0,0, 1.e-9, 1.0, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech OR-branch converged",got,1,&nbad);
    got=cvg_case(&v,&t,0,3,0,1, 1.e-9, 1.0, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("mech OR-branch + gas network",got,0,&nbad);
    cvg_chks("  blocker",v.blocker,"SolutionChange",&nbad);

    /* the 1e-8 floor accepts a correction nobody can call significant */
    got=cvg_case(&v,&t,0,3,0,1, 1.e-3, 1.e-9, 1.e30, 1.0, 1.0, 1.0);
    cvg_chki("mech correction floor",got,1,&nbad);

    /* every leaf is recorded even though the first one failed: a table
       with holes cannot answer "what else is not yet satisfied" */
    got=cvg_case(&v,&t,0,1,0,0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    cvg_chki("all leaves evaluated at iit=1",v.nclause,10,&nbad);

    /* which tolerance applies is part of the verdict.  Load increment
       still large and iit<=ip -> ran; past ip -> rap; increment small
       -> the tight end-of-load pair. */
    for(j=0;j<6;j++){ram[j]=0.;ram1[j]=1.e30;ram2[j]=1.e30;}
    for(j=0;j<5;j++) cam[j]=0.;
    for(j=0;j<4;j++) qa[j]=0.;
    uam[0]=uam[1]=1.; qam[0]=qam[1]=1.;
    qa[0]=1.0;                       /* qa[0] > ea*qam[0] */
    converge_verdict(&v,&t,0,2,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,qam,
                     c1,c2);
    cvg_chk("c1 at iit<=ip is ran",c1[0],0.005,0.,&nbad);
    cvg_chk("c2 at iit<=ip is can",c2[0],0.01,0.,&nbad);
    converge_verdict(&v,&t,0,9,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,qam,
                     c1,c2);
    cvg_chk("c1 past ip is rap",c1[0],0.02,0.,&nbad);
    qa[0]=1.e-9;                     /* qa[0] <= ea*qam[0] */
    converge_verdict(&v,&t,0,9,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,qam,
                     c1,c2);
    cvg_chk("c1 at small increment is ea",c1[0],1.e-5,0.,&nbad);
    cvg_chk("c2 at small increment is cae",c2[0],1.e-3,0.,&nbad);

    /* ram2 ratchets down to the best residual seen, and never up */
    ram1[0]=3.0; ram2[0]=5.0;
    converge_verdict(&v,&t,0,9,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,qam,
                     c1,c2);
    cvg_chk("ram2 ratchets down to ram1",ram2[0],3.0,0.,&nbad);
    ram1[0]=8.0;
    converge_verdict(&v,&t,0,9,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,qam,
                     c1,c2);
    cvg_chk("ram2 does not ratchet up",ram2[0],3.0,0.,&nbad);

    /* THE TWO ASYMMETRIES, pinned so that a later tidy-up has to be a
       decision rather than an accident.  Both are in the stock source and
       both are invisible until the leaves have names; 05-DEBT.md carries
       them.

       1. iflagact gates the mechanical form but NOT the thermomechanical
          one: the same changing contact set that blocks ithermal<2
          converges at ithermal==3. */
    for(j=0;j<6;j++){ram[j]=0.;ram1[j]=1.e30;ram2[j]=1.e30;}
    for(j=0;j<5;j++) cam[j]=0.;
    for(j=0;j<4;j++) qa[j]=0.;
    uam[0]=uam[1]=1.; qam[0]=qam[1]=1.;
    ram[0]=ram[1]=1.e-9; cam[0]=cam[1]=1.e-9; qa[0]=qa[1]=1.e-9;
    got=converge_verdict(&v,&t,3,3,1,1,0,1.e30,ram,ram1,ram2,cam,uam,qa,
                         qam,c1,c2);
    cvg_chki("thermomech ignores iflagact",got,1,&nbad);
    got=converge_verdict(&v,&t,0,3,1,1,0,1.e30,ram,ram1,ram2,cam,uam,qa,
                         qam,c1,c2);
    cvg_chki("  mech with same numbers blocked",got,0,&nbad);

    /*  2. the ral leaf carries an extra iit>1 in the pure-thermal form
           (a dated change in the stock source) and not in the thermo-
           mechanical one.  At iit==1 with only that leaf available the
           two forms disagree by construction. */
    for(j=0;j<6;j++){ram[j]=0.;ram1[j]=1.e30;ram2[j]=1.e30;}
    for(j=0;j<5;j++) cam[j]=0.;
    for(j=0;j<4;j++) qa[j]=0.;
    uam[0]=uam[1]=1.; qam[0]=qam[1]=1.;
    ram[1]=1.e-9;      /* under ral*qam  -> the OR branch's third leaf   */
    cam[1]=1.0;        /* over  c2*uam   -> the primary leaf fails       */
    qa[1]=1.0;         /* over  ea*qam   -> LoadIncrementSmall fails     */
    ram2[1]=0.;        /* projected leaf fails: 0 < 0 is false           */
    got=converge_verdict(&v,&t,2,1,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,
                         qam,c1,c2);
    cvg_chki("thermal ral leaf gated on iit>1",got,0,&nbad);
    for(j=0;j<6;j++){ram1[j]=1.e30;ram2[j]=1.e30;}
    ram2[1]=0.;
    got=converge_verdict(&v,&t,2,2,1,0,0,1.e30,ram,ram1,ram2,cam,uam,qa,
                         qam,c1,c2);
    cvg_chki("  same numbers at iit=2 converge",got,1,&nbad);
  }


  printf("[CONVERGE] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
