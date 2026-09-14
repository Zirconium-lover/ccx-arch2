/*     CalculiX - damage/fracture extension                              */
/*     release.c: what did deleting those elements release?              */

/* Why this module exists
   ----------------------
   When a batch of elements leaves the assembly, the internal force they
   were carrying has to go somewhere.  Where it goes is the question this
   file answers, and it is the only question it answers.

   Before this file, it was answered in five places inside nonlingeo() and
   owned by none of them.  Twelve locals - damage_release_armed, _pass,
   _probe, _qa, _qam, _dt, _rebuild, _iforbou, _nterm, _nother, _nisl,
   _ncoh - plus two arrays, spread across a configuration site, a report,
   and TWO ARM BLOCKS OF NINETY-ONE LINES THAT WERE CHARACTER-IDENTICAL.

   That last part is the whole argument.  docs/ARCHITECTURE.md says
   duplication is evidence rather than a verdict - "the forty-line deletion
   block appearing twice mattered because both copies were one decision, not
   because they were identical".  Here the copies are one decision: capture
   f_int at the converged state before remastruct renumbers anything.  They
   sit 1,120 lines apart, and nothing but reading both told you they were
   the same.

   How the measurement works, and why it is trustworthy
   ---------------------------------------------------
   The arm runs at the converged state, before the topology changes.  It
   maps f_int out of equation numbering into NODE space, because nactdof is
   about to be renumbered and an equation index would mean something else
   afterwards.  The report runs after the rebuild, at a displacement state
   that is bit-identical to the converged one - the topology site sets
   idiscon=1, and prediction() with idiscon!=0 copies vold into v with no
   extrapolation (prediction.c:99).  The difference is therefore the
   released internal force and nothing else.

   Three self-checks are printed with the numbers, because a probe trusted
   without them is worse than none: `action=reuse-sparse-graph' means
   nactdof did not change, so `removed' MUST be zero; `anom' counts degrees
   of freedom active after but not before, which a deletion cannot produce;
   `iforbou' would mean a boundary term had been added to f.

   What this is NOT
   ----------------
   It reads and it prints.  It changes no equation, takes no decision, and
   nothing downstream branches on its output - the banner says so in as many
   words, and it is why this file sits among the observers.        */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"
#include "ccxfork.h"
#include "ccxopt.h"

void release_init(release *r){ memset(r,0,sizeof(*r)); }

/* CCX_DAMAGE_RELEASE_PROBE=1 measures; =2 adds one line per deleted
   element.  Returns nonzero when armed, so the caller's banner stays the
   caller's. */
ITG release_configure(release *r)
{
  const char *e=ccxopt_getenv("CCX_DAMAGE_RELEASE_PROBE");
  if(e==NULL) return 0;
  r->probe=atoi(e);
  if(r->probe<0) r->probe=0;
  if(r->probe>2) r->probe=2;
  return r->probe;
}

void release_free(release *r)
{
  if(r->frel!=NULL) SFREE(r->frel);
  if(r->ract!=NULL) SFREE(r->ract);
  r->frel=NULL; r->ract=NULL;
}

/* Capture f_int(u*) in node space, plus the bookkeeping of what is leaving.
   Called at the converged state, after the batch is marked and before any
   results() has run - both of which the caller guarantees by where it puts
   the call. */
void release_arm(release *r,const trialctx *mdl,const nlstate *n,
                 const erosion_batch *b,const topo_txn *txn,
                 const erosion_policy *pol,const double *damvisc,
                 ITG rebuild,ITG iforbou,ITG nisl,ITG ncoh)
{
  const ITG *nk=*(mdl->nk),*mi=*(mdl->mi),*nactdof=*(mdl->nactdof);
  const ITG mt=mi[1]+1,ne0=*(mdl->ne0),iinc=*(mdl->iinc);
  const double *f=*(mdl->f),*qa=mdl->qa,*dam=*(mdl->dam);
  const char *lakon=*(mdl->lakon);
  ITG ai,aj,ak;

  if(!r->probe) return;
  if(r->frel==NULL){
    NNEW(r->frel,double,mt**nk);
    NNEW(r->ract,ITG,mt**nk);
  }
  for(ai=0;ai<*nk;ai++){
    for(aj=0;aj<mt;aj++){
      ak=nactdof[mt*ai+aj];
      r->frel[mt*ai+aj]=(ak>0)?f[ak-1]:0.;
      r->ract[mt*ai+aj]=(ak>0)?1:0;
    }
  }
  r->armed=1;
  r->pass++;
  r->qa=qa[0];
  r->qam=n->qam[0];
  r->dt=*(mdl->dtime);
  r->rebuild=rebuild;
  r->iforbou=iforbou;
  r->nterm=b->terminal;
  r->nother=b->marked-b->terminal;
  r->nisl=nisl;
  r->ncoh=ncoh;

  /* level 2: one line per element in this batch.  The degradation it
     carried when it left is printed as a REFERENCE CHARACTERISTIC - it is a
     dimensionless multiplier, not a force, and it is NOT summed into any
     estimate of the released force.  The measured force is the surv/removed
     split in the report.

     It does classify the source for free, though: the terminal trigger
     cannot fire below its own threshold, so an element that left at
     Dvis < delete_d did NOT come from it - it came from damfloat, which
     reads no damage variable at all and therefore removes at whatever the
     element was carrying. */
  if(r->probe>=2){
    ITG pe,pj,pnip,pel;
    double pd,pdv,pdmax,pdvmax,ptrig;
    for(pe=0;pe<txn->count;pe++){
      pel=txn->elem[pe]-1;
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
        if(damvisc!=NULL){
          pdv=damvisc[mi[0]*pel+pj];
          if(pdv<0.) pdv=0.;
          if(pdv>1.) pdv=1.;
          if(pdv>pdvmax) pdvmax=pdv;
        }
      }
      ptrig=((pol->delete_visc==1)&&(damvisc!=NULL))?pdvmax:pdmax;
      /* The bulk damage variable exists ONLY for bulk elements.
         calcdamage.f:135 skips every lakon(1:1) != 'C', so dam is
         identically zero for a UC6 facet - printing it would read as "left
         fully intact" when in truth the facet was conducting gmin*Kn and
         its state lives in xstate, not here.  Both decks in play carry UC6,
         so this is not hypothetical: damfloatcoh removals land in the same
         tentative batch. */
      if(lakon[8*pel]!='C'){
        printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
               " pass=%" ITGFORMAT " el=%" ITGFORMAT
               " mat=%" ITGFORMAT " type=%.8s"
               " D=n/a Dvis=n/a g_ref=n/a"
               " note=non-bulk-element-damage-lives-in-xstate%s",
               iinc,r->pass,txn->elem[pe],
               txn->mat[pe],&lakon[8*pel],
               "\n");
        continue;
      }
      printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
             " pass=%" ITGFORMAT " el=%" ITGFORMAT
             " mat=%" ITGFORMAT " type=%.8s"
             " D=%.6f Dvis=%.6f g_ref=%.4e"
             " below_delete_d=%s%s",
             iinc,r->pass,txn->elem[pe],
             txn->mat[pe],&lakon[8*pel],pdmax,pdvmax,1.-ptrig,
             (ptrig<pol->delete_d)?"YES-not-terminal":"no",
             "\n");
    }
    fflush(stdout);
  }
}

/* Read the captured force against the rebuilt one and say what moved.
   Disarms itself: one arm, one report. */
void release_report(release *r,const trialctx *mdl,double steptime)
{
  const ITG *nk=*(mdl->nk),*mi=*(mdl->mi),*nactdof=*(mdl->nactdof);
  const ITG mt=mi[1]+1,iinc=*(mdl->iinc);
  const double *f=*(mdl->f);
  ITG ri,rj,rk,rns=0,rds=0,rnr=0,rdr=0,ranom=0;
  double dfv,dfa,smax=0.,sl1=0.,sl2=0.,rmax=0.,rl1=0.,rl2=0.;

  if(!((r->probe)&&(r->armed)&&(r->frel!=NULL)&&(r->ract!=NULL))) return;

  for(ri=0;ri<*nk;ri++){
    for(rj=0;rj<mt;rj++){
      rk=nactdof[mt*ri+rj];
      dfa=(rk>0)?f[rk-1]:0.;
      if(r->ract[mt*ri+rj]){
        dfv=r->frel[mt*ri+rj]-dfa;
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
         iinc,r->pass,steptime,
         r->rebuild?"remastruct":"reuse-sparse-graph",
         r->dt,"\n",
         smax,rns,rds,sl1,sl2,"\n",
         rmax,rnr,rdr,rl1,rl2,"\n",
         r->qa,r->qam,
         (r->qam>0.)?smax/r->qam:-1.,
         (r->qam>0.)?((smax>rmax?smax:rmax)/r->qam):-1.,"\n",
         r->nterm,r->nother,r->nisl,r->ncoh,"\n",
         ranom,r->iforbou,
         r->rebuild?"no":"YES","\n");
  fflush(stdout);
  r->armed=0;
}
