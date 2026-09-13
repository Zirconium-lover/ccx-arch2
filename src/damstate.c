/*     CalculiX - damage/fracture extension                              */
/*     damstate.c: THE judgement about what no longer carries load.      */

/* Why this module exists
   ----------------------
   The code decides, per node, that a degree of freedom has lost its load
   path.  Before this module that judgement was computed inline and then told
   to consumers PIECEMEAL: the displacement norm learned it through
   CCX_DAMAGE_AUTOSPC, the force norm through CCX_DAMAGE_AUTOSPC_FORCE, the
   termination connectivity through CCX_FRACTURE_DEADFACET, and the solve
   itself was never told at all.

   The second s3rad wall was exactly the gap between the first two: AUTOSPC
   knew node 1246 had lost its load path and told the displacement norm,
   while the force norm - which is what actually vetoed the increment - was
   never told.  Closing that gap passed the wall and immediately exposed the
   next consumer.  Every wall of that family is the same physical state
   surfacing through a consumer nobody wired up, so the judgement belongs in
   one place with one owner.

   This module owns the computation.  It changes no equation and takes no
   decision of its own; consumers ask it and act.

   Contract
   --------
     - one update per assembled operator, before anything reads it;
     - a node is judged only when all three of its translational dofs are
       active, because a constrained dof cannot appear in the norms anyway;
     - a node counts against its OWN intact diagonal, so the ratio is
       dimensionless and no mesh-wide median is needed;
     - the module refuses to arm if its self test fails, the discipline
       lsladder.c set.                                                    */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

/* A non-positive assembled diagonal is not a small number, it is a different
   object: the node sits on a descending branch.  Treating it as dead is a
   bigger step than treating a collapsed one as dead, so it is separate and
   off by default. */

void damstate_init(damstate *s,ITG nk,double g,ITG allowneg)
{
  s->nk=nk; s->g=g; s->allowneg=allowneg; s->ndead=0;
  NNEW(s->ok,ITG,nk);
  NNEW(s->dead,ITG,nk);
  NNEW(s->diag,double,nk);
  NNEW(s->diag0,double,nk);
}

void damstate_free(damstate *s)
{
  if(s->ok!=NULL) SFREE(s->ok);
  if(s->dead!=NULL) SFREE(s->dead);
  if(s->diag!=NULL) SFREE(s->diag);
  if(s->diag0!=NULL) SFREE(s->diag0);
  s->ok=NULL; s->dead=NULL; s->diag=NULL; s->diag0=NULL; s->nk=0; s->ndead=0;
}

/* Recompute the per-node diagonal and the judgement from the operator that
   is about to be solved.  `diag0` keeps the FIRST positive value ever seen
   for a node, which is its intact reference; it never drifts downward. */
void damstate_update(damstate *s,const double *ad,const ITG *nactdof,ITG mt)
{
  ITG i,idir,k,first;
  double dmin;

  s->ndead=0;
  for(i=0;i<s->nk;i++){
    s->ok[i]=1; dmin=0.; first=1;
    for(idir=1;idir<=3;idir++){
      k=nactdof[mt*i+idir];
      if(k<=0){s->ok[i]=0;break;}
      /* A real minimum.  The historical defect here was
         `if((dmin<0.)||(ad<dmin))`, which stops being a minimum as soon as
         one negative value enters - the first clause is then true on every
         later dof, so the stored number is the LAST negative diagonal, not
         the smallest.  Self test case F pins that. */
      if(first||(ad[k-1]<dmin)){dmin=ad[k-1];first=0;}
    }
    if(s->ok[i]==0){s->diag[i]=0.;s->dead[i]=0;continue;}
    s->diag[i]=dmin;
    if((s->diag0[i]<=0.)&&(dmin>0.)) s->diag0[i]=dmin;

    s->dead[i]=0;
    if(s->g<=0.) continue;
    if(s->diag0[i]<=0.) continue;
    if(dmin<=0.){
      if(s->allowneg){s->dead[i]=1;s->ndead++;}
      continue;
    }
    if(dmin < s->g*s->diag0[i]){s->dead[i]=1;s->ndead++;}
  }
}

ITG damstate_dead(const damstate *s,ITG node)
{
  if((s==NULL)||(s->dead==NULL)) return 0;
  if((node<0)||(node>=s->nk)) return 0;
  return s->dead[node];
}

/* A cohesive facet is dead only when EVERY one of its integration points has
   failed.  The flag is xstate slot 4, which resultsmech_uc6.f writes
   explicitly and nothing else reads.  UC6 carries exactly three points, and
   looping to the allocated stride mi[0] instead would average in
   uninitialised slots - the defect E-84 found in the nonlocal average. */
ITG damstate_facet_dead(const double *xstate,ITG nstate,ITG mi0,
                        ITG elem,ITG nip)
{
  ITG j;
  if((xstate==NULL)||(nstate<4)||(nip<1)) return 0;
  if(nip>mi0) nip=mi0;
  for(j=0;j<nip;j++){
    if(xstate[3+nstate*(j+mi0*elem)]<0.5) return 0;
  }
  return 1;
}

/* How many cohesive facets hold each node.
 *
 * One owner for a count that used to be an inlined loop in nonlingeo.c.
 * CCX_DAMAGE_DEADALL skips any node a cohesive facet holds, on the argument
 * that a node tied to the other side of an interface is not free however
 * dead its bulk is.  That argument is right: 286 of 318 cohesive-only nodes
 * on m12_epsf50 are well supported by it.
 *
 * A NARROWED COUNT WAS TRIED AND RETIRED.  The count includes facets that
 * have SEPARATED - every integration point failed - and the argument for
 * skipping them is easy to make: a separated facet ties nothing.  It was
 * implemented as CCX_DAMAGE_DEADALL_FACET, with a self test that could be
 * shown failing, bit-identical when off.  Seven s3rad arms on one binary
 * (research/14-THE-TRAP.md) then said:
 *
 *   - it buys nothing.  The arm carrying it lands where the arm without it
 *     lands, to four figures on the terminal grip reaction;
 *   - in one configuration it is the ONLY difference between a run that
 *     reaches theta 0.5569 and one that stops at 0.2550, because the twelve
 *     extra elements it deletes are enough to put the trajectory into a
 *     trap that a one-element perturbation decides.
 *
 * By CLAUDE.md's test - name the failure it addresses and the gate case that
 * would go red without it - it had neither, so it is gone.  This comment is
 * here so the next reader does not re-derive it: the idea is plausible, it
 * was measured, and the measurement did not support it.
 *
 * Returns the number of facets counted; nfac[nk] receives the per-node count
 * and is zeroed first. */
ITG damstate_facet_support(const ITG *ipkon,const char *lakon,const ITG *kon,
                           ITG ne,ITG nk,ITG *nfac)
{
  ITG i,j,n,nope,ntot=0;

  if((nfac==NULL)||(nk<=0)) return 0;
  for(n=0;n<nk;n++) nfac[n]=0;
  if((ipkon==NULL)||(lakon==NULL)||(kon==NULL)) return 0;

  /* user elements carry their node count in byte 7 of lakon, the idiom the
     node dump uses. */
  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    nope=(ITG)((unsigned char)lakon[8*i+7]);
    if((nope<1)||(nope>20)) continue;
    ntot++;
    for(j=0;j<nope;j++){
      n=kon[ipkon[i]+j]-1;
      if((n>=0)&&(n<nk)) nfac[n]++;
    }
  }
  return ntot;
}

/* ---------------------------------------------------------------- tests */

static ITG dst_chk(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-34s got=%-4" ITGFORMAT " want=%-4" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG damstate_selftest(void)
{
  ITG nbad=0,mt=4,nk=5,i;
  ITG nactdof[20];
  double ad[16];
  damstate s;

  printf("[DAMSTATE] self test%s","\n");

  /* five nodes, three active dofs each except node 3 which has one
     constrained; equations are numbered 1..12 then 13..15 */
  for(i=0;i<nk*mt;i++) nactdof[i]=0;
  for(i=0;i<nk;i++){
    nactdof[mt*i+1]=3*i+1;
    nactdof[mt*i+2]=3*i+2;
    nactdof[mt*i+3]=3*i+3;
  }
  nactdof[mt*3+2]=0;                       /* node 3: one dof constrained */

  for(i=0;i<15;i++) ad[i]=1000.;

  damstate_init(&s,nk,1.e-3,0);

  /* first pass: everything intact, nothing dead, references captured */
  damstate_update(&s,ad,nactdof,mt);
  dst_chk("A intact: none dead",s.ndead,0,&nbad);
  dst_chk("A node 3 not judged (dof held)",s.ok[3],0,&nbad);
  dst_chk("A intact reference captured",(s.diag0[0]==1000.),1,&nbad);

  /* B: node 1 collapses below the threshold */
  ad[3]=0.5; ad[4]=0.6; ad[5]=0.7;         /* node 1, min 0.5 < 1e-3*1000 */
  damstate_update(&s,ad,nactdof,mt);
  dst_chk("B collapsed node is dead",damstate_dead(&s,1),1,&nbad);
  dst_chk("B only that one",s.ndead,1,&nbad);
  dst_chk("B healthy node still alive",damstate_dead(&s,0),0,&nbad);

  /* C: exactly at the threshold is NOT dead (strict inequality) */
  ad[3]=1.0; ad[4]=1.0; ad[5]=1.0;         /* 1.0 == 1e-3*1000 */
  damstate_update(&s,ad,nactdof,mt);
  dst_chk("C exactly at threshold is alive",damstate_dead(&s,1),0,&nbad);

  /* D: the intact reference does not drift down with the collapse */
  dst_chk("D reference unchanged",(s.diag0[1]==1000.),1,&nbad);

  /* E: a non-positive diagonal is not dead unless allowneg */
  ad[3]=-5.; ad[4]=100.; ad[5]=100.;
  damstate_update(&s,ad,nactdof,mt);
  dst_chk("E negative diagonal not dead by default",damstate_dead(&s,1),0,&nbad);
  {
    damstate t; ITG j;
    damstate_init(&t,nk,1.e-3,1);
    for(j=0;j<nk;j++) t.diag0[j]=1000.;
    damstate_update(&t,ad,nactdof,mt);
    dst_chk("E negative diagonal dead when asked",damstate_dead(&t,1),1,&nbad);
    damstate_free(&t);
  }

  /* F: the minimum is a REAL minimum with a negative present.  The old
     inline code returned the LAST negative value, not the smallest. */
  ad[3]=-5.; ad[4]=-50.; ad[5]=-1.;
  damstate_update(&s,ad,nactdof,mt);
  dst_chk("F min over dofs is the smallest",(s.diag[1]==-50.),1,&nbad);

  /* G: threshold zero disarms the judgement entirely */
  {
    damstate t; ITG j;
    damstate_init(&t,nk,0.,0);
    for(j=0;j<nk;j++) t.diag0[j]=1000.;
    ad[3]=1.e-9; ad[4]=1.e-9; ad[5]=1.e-9;
    damstate_update(&t,ad,nactdof,mt);
    dst_chk("G threshold 0 judges nobody",t.ndead,0,&nbad);
    damstate_free(&t);
  }

  /* H: facet is dead only when EVERY point has failed */
  {
    ITG nstate=5,mi0=3,e=2,j;
    double xs[64];
    for(j=0;j<64;j++) xs[j]=0.;
    for(j=0;j<3;j++) xs[3+nstate*(j+mi0*e)]=1.;
    dst_chk("H all points failed -> dead",
            damstate_facet_dead(xs,nstate,mi0,e,3),1,&nbad);
    xs[3+nstate*(1+mi0*e)]=0.;
    dst_chk("H one point alive -> not dead",
            damstate_facet_dead(xs,nstate,mi0,e,3),0,&nbad);
  }

  /* I: per-node facet support.  Three elements: one bulk (must be ignored),
     two six-node facets sharing nodes 4,5,6. */
  {
    ITG ne3=3,nk10=10,j,ntot;
    ITG ipkon[3],kon[64],nfac[10];
    char lak[24];

    for(j=0;j<24;j++) lak[j]=' ';
    strcpy1(&lak[0],"C3D4",4); lak[7]=(char)4;
    lak[8]='U'; lak[8+1]='C'; lak[8+2]='6'; lak[8+7]=(char)6;
    lak[16]='U'; lak[16+1]='C'; lak[16+2]='6'; lak[16+7]=(char)6;
    ipkon[0]=0; ipkon[1]=4; ipkon[2]=10;
    for(j=0;j<4;j++)  kon[j]=j+1;            /* bulk on nodes 1..4       */
    for(j=0;j<6;j++)  kon[4+j]=j+1;          /* facet 1 on nodes 1..6    */
    for(j=0;j<6;j++)  kon[10+j]=j+4;         /* facet 2 on nodes 4..9    */

    ntot=damstate_facet_support(ipkon,lak,kon,ne3,nk10,nfac);
    dst_chk("I two facets counted",ntot,2,&nbad);
    dst_chk("I shared node holds twice",nfac[3],2,&nbad);
    dst_chk("I private node holds once",nfac[0],1,&nbad);
    dst_chk("I bulk element is not a facet",nfac[8],1,&nbad);

    ipkon[2]=-ipkon[2]-2;                            /* facet 2 deleted   */
    ntot=damstate_facet_support(ipkon,lak,kon,ne3,nk10,nfac);
    dst_chk("I deleted facet is not counted",ntot,1,&nbad);
    dst_chk("I deleted facet holds nothing",nfac[8],0,&nbad);
    dst_chk("I its shared node keeps the live one",nfac[3],1,&nbad);
  }

  damstate_free(&s);
  printf("[DAMSTATE] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
