/*     CalculiX - damage/fracture extension                              */
/*     erosion.c: which elements leave the assembly, and why.            */

/* Why this module exists
   ----------------------
   "Where is it decided that an element is deleted?" had no answer that fit
   in one place.  The decision was five file-static functions inside
   nonlingeo.c and, on top of them, a forty-line block - terminal scan,
   then DEADALL, then DEADSOLE, each against the same batch budget, each
   with its own banner - that appeared TWICE, at the two points where a
   converged state is scanned.  The two copies differ by indentation and by
   nothing else; the diff is in the commit that created this file.

   That is the same shape topology.c was written for, and the same cause:
   nobody chose to write the scan twice, there was simply nothing to call.

   Contract
   --------
     - erosion_mark() is the ONE place that decides who leaves.  It applies
       the three rules in their established order against one batch budget
       and reports what it took.  A rule that is off costs nothing;
     - the three rules are file-static here, so "who else calls DEADALL"
       has the answer `nobody, by construction' rather than a grep;
     - it MARKS - ipkon -> -ipkon-2 - and nothing else.  Rebuilding the
       equation structure, committing the batch to .damage and rolling it
       back are topology.c's and nonlingeo's, and this file contains no
       code that pretends otherwise;
     - the policy is a struct read once from the configuration, not eleven
       arguments threaded through two call sites that can drift apart.

   What this is NOT
   ----------------
   Not the damage model: what D is at an integration point is calcdamage's
   answer and this file only reads it.  Not a criterion change - the three
   rules, their order, their thresholds and their printed banners are the
   ones that were there, which is why the acceptance test is bit identity
   on the whole gate.                                                    */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "CalculiX.h"

/* What counts as softening in the present trial.  Same 1.e-12 as before,
   under the name of the question it answers rather than of the one caller
   that happened to be written first (it was EROSION_SOFTENING_D_TOL). */
#define EROSION_SOFTENING_D_TOL 1.e-12

/* Progressive damage material classifier shared by DE1 and DM2.0.
   Rice-Tracey + Evolution=Displacement keeps the historical four-constant
   signature.  DM2.0 is identified by model type 3 and a variable-length
   constant count 3+2*NPOINTS (NPOINTS>=2). */
ITG damage_progressive_material(ITG imat,const ITG *ndmcon,
                                const double *dmcon,ITG ndmat,
                                ITG ntmat)
{
  ITG nconst,type,off;

  if(imat<1) return 0;
  nconst=ndmcon[2*(imat-1)];
  if((dmcon==NULL)||(ndmat<1)||(ntmat<1)) return 0;

  off=1+(ndmat+1)*ntmat*(imat-1);
  type=(ITG)dmcon[off];
  if((type==1)&&(nconst==4)) return 1;
  if((type==3)&&(nconst>=7)&&(((nconst-3)%2)==0)) return 1;

  return 0;
}

/* Detect actual progressive softening in the present Newton trial.  Merely
   having a DE1/DM2.0 material in the model is not enough: at least one active
   integration point must have D_trial>D_committed.  Comparing degradation D
   rather than the overloaded raw dam value also handles initiation crossing
   (omega<1 -> dam=1+D) without a false large jump. */
ITG erosion_softening(const double *dam,
                                      const double *dambase,
                                      const ITG *ipkon,const char *lakon,
                                      const ITG *ielmat,ITG mi2,
                                      const ITG *ndmcon,const double *dmcon,
                                      ITG ndmat,ITG ntmat,ITG ne0,ITG mi0,
                                      ITG *nsoft,double *maxdd)
{
  ITG i,j,nip,imat,elementsoft;
  double dtrial,dbase,dd;

  *nsoft=0;
  *maxdd=0.;
  if((dam==NULL)||(dambase==NULL)) return 0;

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='C') continue;

    imat=ielmat[mi2*i];
    if(!damage_progressive_material(imat,ndmcon,dmcon,ndmat,ntmat))
      continue;

    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    elementsoft=0;

    for(j=0;j<nip;j++){
      dtrial=dam[mi0*i+j]-1.;
      dbase=dambase[mi0*i+j]-1.;
      if(dtrial<0.) dtrial=0.;
      if(dbase<0.) dbase=0.;
      dd=dtrial-dbase;
      if(dd>EROSION_SOFTENING_D_TOL){
        elementsoft=1;
        if(dd>*maxdd) *maxdd=dd;
      }
    }
    if(elementsoft) (*nsoft)++;
  }

  return (*nsoft>0)?1:0;
}

/* Mark a bounded batch of DE1.2 C3D4 elements for terminal deletion.
   DE1 uses dam = 1 + D after initiation, so D is recovered locally without
   changing the public history layout.  Both Rice-Tracey DE1 and DM2.0
   tabulated ductile progressive materials are eligible.  The element is only
   marked here; remastruct and transactional
   commit/rollback remain under nonlingeo's existing A3 machinery. */
/* Terminal deletion can be restricted to named materials.

   The TP1 ladder needs its rungs to differ ONLY in which phase is allowed
   to erode, with the constitutive routines untouched:

     CCX_DAMAGE_DELETE_MAT=NONE   continuous damage, no topology change
     CCX_DAMAGE_DELETE_MAT=ZrH    only the hydride erodes
     CCX_DAMAGE_DELETE_MAT=ALL    stock behaviour (default when unset)

   The list is comma separated and matched case-insensitively against the
   *MATERIAL names.  A rejected element is removed from the terminal scan
   only; its damage keeps evolving exactly as before, so the comparison
   isolates the topology change and nothing else. */

static ITG erosion_material_allowed(ITG imat,const char *matname,
                                 const char *filter)
{
  const char *p,*q;
  char nm[81];
  ITG i,n;

  if(filter==NULL) return 1;
  if((strcmp(filter,"ALL")==0)||(strcmp(filter,"all")==0)) return 1;
  if((strcmp(filter,"NONE")==0)||(strcmp(filter,"none")==0)) return 0;
  if((matname==NULL)||(imat<1)) return 1;

  n=0;
  for(i=0;i<80;i++){
    if(matname[80*(imat-1)+i]==' ') break;
    nm[i]=matname[80*(imat-1)+i];
    n++;
  }
  nm[n]=0;
  if(n==0) return 1;

  p=filter;
  while(*p!=0){
    q=strchr(p,',');
    if(q==NULL) q=p+strlen(p);
    if((ITG)(q-p)==n){
      for(i=0;i<n;i++){
        if(tolower((unsigned char)p[i])!=tolower((unsigned char)nm[i])) break;
      }
      if(i==n) return 1;
    }
    p=(*q==0)?q:q+1;
  }
  return 0;
}

/* Deletes a DEAD element that is the SOLE support of a node.

   Measured configuration (E-61): on m14_fine node 776 had exactly one live
   element, that element stood at D = 1.0000, and the node travelled 9.63 mm
   while the other three nodes of the same element moved 0.44.  The element's
   longest edge went from 0.0974 to 9.809 - a stretch of 101x on a 4 mm
   specimen.  An element below a per cent of its stiffness carries almost
   nothing, so a node whose whole support is one such element is very nearly
   free and Newton solves for it.

   Holding that node with a diagonal term was tried and does not cure it
   (E-61): the added stiffness would have to be non-perturbative to compete
   with the element itself.  Removing the dead element instead is physically
   near-free - it was carrying under 1% - and it hands the node to the
   existing BK4 / damfloat machinery, which is built for a node that has lost
   all its bulk.

   This is NOT damdangle.  damdangle counted support without ever reading
   degradation, so it deleted healthy load-bearing material and emptied small
   meshes (E-22).  Here the element must itself be dead, and it must be the
   only thing a node has.  Both conditions are necessary and both are
   measured, not assumed.

   The batch is bounded exactly like the terminal batch, and marking uses the
   same ipkon -> -ipkon-2 convention. */
static ITG erosion_mark_deadsole(const double *dam,const double *visc,
                                     ITG usevisc,ITG *ipkon,const char *lakon,
                                     const ITG *kon,ITG nk,ITG ne0,ITG mi0,
                                     double gdead,ITG batchmax)
{
  ITG i,j,n,nip,nnew=0,usedam,*nlive=NULL;
  const double *src;
  double dmx,g;

  if((nk<=0)||(batchmax<=0)) return 0;
  usedam=((usevisc&&(visc!=NULL))?0:1);
  src=usedam?dam:visc;

  NNEW(nlive,ITG,nk);
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    /* kon is 0-based here: the nodes of element i are kon[ipkon[i]+0] to
       kon[ipkon[i]+nope-1].  frd.c:1684 takes the last node as
       kon[ipkon[i]+nope-1], and the VTK writer below (which uses j=0..3)
       reproduces the deck connectivity exactly - 27359 of 27359 cells on
       m12_eta15.  An earlier version of this loop ran j=1..4, which skips
       node 1 and picks up node 1 of the NEXT element instead. */
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n>=0)&&(n<nk)) nlive[n]++;
    }
  }

  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    dmx=0.;
    for(j=0;j<nip;j++){
      /* dam holds 1+D (the integer part flags initiation); visc holds D
         directly.  Reading dam raw gives 1.0024 for D=0.0024 and deletes a
         healthy element - measured, and it is exactly how damdangle went
         wrong.  mark_terminal decodes it the same way. */
      double dd=usedam?(src[mi0*i+j]-1.):src[mi0*i+j];
      if(dd<0.) dd=0.;
      if(dd>1.) dd=1.;
      if(dd>dmx) dmx=dd;
    }
    g=1.-dmx; if(g<0.) g=0.;
    if(g>=gdead) continue;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      if(nlive[n]==1){
        printf("[DAMAGE DEADSOLE]   element %" ITGFORMAT " g=%.6e sole "
               "support of node %" ITGFORMAT "\n",
               i+1,g,n+1);
        ipkon[i]=-ipkon[i]-2;
        nnew++;
        /* the node counts are now stale for this element's nodes; drop them
           so a second element at the same node cannot also be taken in the
           same pass */
        for(j=0;j<4;j++){
          n=kon[(-ipkon[i]-2)+j]-1;
          if((n>=0)&&(n<nk)) nlive[n]=0;
        }
        break;
      }
    }
  }
  SFREE(nlive);
  return nnew;
}

/* DEADALL - the Class A artefact (E-69).
 *
 * DEADSOLE above asks whether a dead element is the SOLE support of a node.
 * The three runs that died by divergence died on a node with TWO, not one:
 *
 *   uf50   node 2471  2 live tets, 0 facets, both D=1.0000, edge stretch 22.2
 *   eta1e4 node 3354  2 live tets, 0 facets, both D=1.0000, edge stretch 49.8
 *   eta15  node 1581  2 live tets, 0 facets, both D=1.0000, edge stretch 35.5
 *
 * In each of those models exactly ONE node in the whole mesh had its entire
 * live support dead - out of 17, 145 and 65 low-support nodes respectively -
 * and it was the node the solver threw.  Three out of three, with a
 * selectivity of one in seventeen to one in a hundred and forty-five.
 *
 * The node is a free swinging point: everything holding it carries only gmin,
 * so Newton solves it almost unconstrained and draws the elements out into a
 * needle.  Note that a needle PRESERVES VOLUME - it stretches along one
 * direction and collapses across the other two - so V/V0 does not see it
 * (1.2 on eta15 while an edge went 0.120 -> 4.151 mm).  Do not look for this
 * with a volume ratio.
 *
 * SAFETY.  A node qualifies only when every live element at it is dead, so
 * every element this routine deletes is itself dead: element i is live and
 * touches a qualifying node n, and n qualifies only if all its live elements
 * are dead.  Nothing load-bearing can be removed.  That is the property
 * damdangle did not have (E-22), and it is what makes this different.
 *
 * The facet guard is the second half.  A node still tied to the other side of
 * an interface is not free, whatever its bulk looks like, and the
 * cohesive-only node is Class B - a CONDITIONING defect (E-67, E-68) that must
 * not be answered by deleting material.  286 of 318 cohesive-only nodes on
 * m12_epsf50 are perfectly well supported.
 *
 * A NARROWED GUARD WAS TRIED AND RETIRED.  The count includes facets that
 * have SEPARATED, and a separated facet ties nothing - so counting only
 * LIVE facets looks like a correction rather than a loosening.  It was
 * implemented as CCX_DAMAGE_DEADALL_FACET and measured on seven s3rad arms
 * built from one binary.  It bought nothing: the
 * arm carrying it lands where the arm without it lands, to four figures on
 * the terminal grip reaction.  Worse, in one configuration it was the only
 * difference between a run reaching theta 0.5569 and one stopping at
 * 0.2550 - the twelve extra elements it deleted put the trajectory into a
 * trap that a ONE-element perturbation decides.  It is gone.  The count
 * lives in damstate_facet_support() with its own self test; this comment
 * exists so the idea is not re-derived without the measurement.
 */
static ITG erosion_mark_deadall(const double *dam,const double *visc,
                                    ITG usevisc,ITG *ipkon,const char *lakon,
                                    const ITG *kon,ITG nk,ITG ne,ITG ne0,
                                    ITG mi0,double gdead,ITG batchmax,
                                    ITG *nnodes)
{
  ITG i,j,n,nip,nnew=0,usedam;
  ITG *nlive=NULL,*ndead=NULL,*nfac=NULL,*take=NULL;
  const double *src;
  double dmx,g;

  if(nnodes!=NULL) *nnodes=0;
  if((nk<=0)||(batchmax<=0)) return 0;
  usedam=((usevisc&&(visc!=NULL))?0:1);
  src=usedam?dam:visc;

  NNEW(nlive,ITG,nk);
  NNEW(ndead,ITG,nk);
  NNEW(nfac,ITG,nk);

  /* cohesive support per node - damstate.c owns the count and its self test */
  damstate_facet_support(ipkon,lakon,kon,ne,nk,nfac);

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    dmx=0.;
    for(j=0;j<nip;j++){
      /* dam holds 1+D, visc holds D - same decode as mark_deadsole. */
      double dd=usedam?(src[mi0*i+j]-1.):src[mi0*i+j];
      if(dd<0.) dd=0.;
      if(dd>1.) dd=1.;
      if(dd>dmx) dmx=dd;
    }
    g=1.-dmx; if(g<0.) g=0.;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      nlive[n]++;
      if(g<gdead) ndead[n]++;
    }
  }

  NNEW(take,ITG,nk);
  for(n=0;n<nk;n++){
    if(nfac[n]>0) continue;
    if(nlive[n]<1) continue;
    if(ndead[n]!=nlive[n]) continue;
    take[n]=1;
    if(nnodes!=NULL) (*nnodes)++;
  }

  /* The qualifying set is a SNAPSHOT.  Deleting an element also takes support
     away from its other three nodes, and resolving that here would need a
     fixed point; the next increment picks it up instead, which is how every
     other batch in this file behaves. */
  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      if(take[n]==0) continue;
      printf("[DAMAGE DEADALL]   element %" ITGFORMAT " deleted: node %"
             ITGFORMAT " has %" ITGFORMAT " live element(s), all dead, "
             "no cohesive facet\n",i+1,n+1,nlive[n]);
      ipkon[i]=-ipkon[i]-2;
      nnew++;
      break;
    }
  }

  SFREE(take);SFREE(nfac);SFREE(ndead);SFREE(nlive);
  return nnew;
}

static ITG erosion_mark_terminal(double *dam,ITG *ipkon,
                                     const char *lakon,const ITG *ielmat,
                                     ITG mi2,const ITG *ndmcon,
                                     const double *dmcon,ITG ndmat,ITG ntmat,
                                     ITG ne0,ITG mi0,double ddelete,
                                     ITG batchmax,double *batch_dmax,
                                     double *trigger_value,ITG *trigger_ip,
                                     const char *matname,
                                     const char *delfilter,
                                     const double *visc,ITG usevisc,
                                     double *batch_vmin)
{
  ITG i,j,nip,imat,nnew=0;
  double de,demax,dv,dvmax,dtrig;

  *batch_dmax=0.;
  if(batch_vmin!=NULL) *batch_vmin=1.;

  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;

    imat=ielmat[mi2*i];
    if(imat<1) continue;
    if(!damage_progressive_material(imat,ndmcon,dmcon,ndmat,ntmat)) continue;
    if(!erosion_material_allowed(imat,matname,delfilter)) continue;

    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;

    demax=0.;
    dvmax=0.;
    for(j=0;j<nip;j++){
      de=dam[mi0*i+j]-1.;
      if(de<0.) de=0.;
      if(de>1.) de=1.;
      if(de>demax) demax=de;
      if(visc!=NULL){
        dv=visc[mi0*i+j];
        if(dv<0.) dv=0.;
        if(dv>1.) dv=1.;
        if(dv>dvmax) dvmax=dv;
      }
    }

    /* The stress is scaled by 1-Dvis, not by 1-D: resultsmech.f uses
       damvisc whenever the viscosity is on.  Triggering deletion on D
       therefore removes an element that is still carrying 1-Dvis of its
       effective stress, and releases that force in one increment at
       constant load.  It also explains why cutting the step never helped
       and sometimes hurt: beta=dt/(eta+dt), so a smaller step makes Dvis
       lag further behind D and the deleted element carries more.
       usevisc makes the trigger read the same variable the stress does. */
    dtrig=demax;
    if((usevisc==1)&&(visc!=NULL)) dtrig=dvmax;

    if(dtrig>=ddelete){
      if(batch_vmin!=NULL){
        if(dvmax<*batch_vmin) *batch_vmin=dvmax;
      }
      /* DE1.3.1: capture the terminal trigger at the instant the element
         first changes topology.  Later same-load redistribution can alter
         dam for already deleted elements, so .damage must not reconstruct
         this value from the current constitutive state. */
      if((trigger_value!=NULL)&&(trigger_ip!=NULL)){
        trigger_value[i]=demax;
        trigger_ip[i]=1;
      }
      ipkon[i]=-ipkon[i]-2;
      nnew++;
      if(demax>*batch_dmax) *batch_dmax=demax;
    }
  }

  return nnew;
}


/* ---- the one place that decides who leaves -----------------------------

   The three rules in the order they have always run, against ONE batch
   budget, with the banners they have always printed.  Both sites that used
   to carry this block now call it.                                      */

ITG erosion_mark(const erosion_policy *p,erosion_batch *b,
                 double *dam,const double *damvisc,
                 ITG *ipkon,const char *lakon,const ITG *kon,
                 const ITG *ielmat,const char *matname,
                 const ITG *ndmcon,const double *dmcon,
                 ITG ndmat,ITG ntmat,ITG nk,ITG ne,ITG ne0,ITG mi0,ITG mi2,
                 double *trigger_value,ITG *trigger_ip,
                 ITG iinc,double steptime)
{
  ITG nda,nds;

  b->deadall=0; b->deadsole=0; b->deadall_nodes=0;
  b->batch_dmax=0.;

  b->marked=erosion_mark_terminal(
      dam,ipkon,lakon,ielmat,mi2,ndmcon,dmcon,ndmat,ntmat,
      ne0,mi0,p->delete_d,p->batchmax,
      &b->batch_dmax,trigger_value,trigger_ip,matname,p->filter,
      damvisc,p->delete_visc,&b->batch_vmin);
  b->terminal=b->marked;

  if((p->deadall_g>0.)&&(b->marked<p->batchmax)){
    nda=erosion_mark_deadall(
        dam,damvisc,p->delete_visc,ipkon,lakon,kon,nk,
        ne,ne0,mi0,p->deadall_g,
        p->batchmax-b->marked,&b->deadall_nodes);
    if(nda>0){
      b->marked+=nda;
      b->deadall=nda;
      b->total_deadall+=nda;
      printf("[DAMAGE DEADALL] inc=%" ITGFORMAT " time=%.12e "
             "nodes=%" ITGFORMAT " deleted=%" ITGFORMAT " total=%"
             ITGFORMAT "\n",iinc,steptime,b->deadall_nodes,
             nda,b->total_deadall);
      fflush(stdout);
    }
  }
  if((p->deadsole_g>0.)&&(b->marked<p->batchmax)){
    nds=erosion_mark_deadsole(
        dam,damvisc,p->delete_visc,ipkon,lakon,kon,nk,
        ne0,mi0,p->deadsole_g,
        p->batchmax-b->marked);
    if(nds>0){
      b->marked+=nds;
      b->deadsole=nds;
      b->total_deadsole+=nds;
      printf("[DAMAGE DEADSOLE] inc=%" ITGFORMAT " time=%.12e "
             "deleted=%" ITGFORMAT " total=%" ITGFORMAT "\n",
             iinc,steptime,nds,b->total_deadsole);
      fflush(stdout);
    }
  }
  return b->marked;
}

void erosion_batch_init(erosion_batch *b)
{
  b->marked=0; b->terminal=0; b->deadall=0; b->deadsole=0;
  b->deadall_nodes=0; b->batch_dmax=0.; b->batch_vmin=1.;
  b->total_deadall=0; b->total_deadsole=0;
}

/* ---- self test ---------------------------------------------------------

   Four decisions are checked, and each one is a decision that had no test
   at all while it lived inside nonlingeo():

     1. the material filter.  A prefix must NOT match - a deck with
        materials ZrH and Zr must be able to erode one and not the other,
        and the length comparison in erosion_material_allowed() is the only
        thing standing between those two;
     2. the trigger variable.  CCX_DAMAGE_DELETE_VISC decides whether the
        threshold is read against D or against Dvis, and the whole argument
        for its default is that the two give DIFFERENT answers.  So an
        element with D past the threshold and Dvis short of it must survive
        with the switch on and go with it off;
     3. the batch budget.  More terminal candidates than batchmax must
        leave exactly batchmax of them, not all;
     4. the marking convention.  ipkon -> -ipkon-2 has to be invertible,
        because the rollback reads ipkon back through it; and an already
        marked element must not be marked twice.

   A tiny mesh is built here rather than read from a deck, so the test runs
   in every job with no input file.                                      */

static void erosion_chki(const char *what,ITG got,ITG want,ITG *nbad)
{
  if(got!=want) (*nbad)++;
  printf("[EROSION]   %-52s %6" ITGFORMAT " (want %" ITGFORMAT ") %s\n",
         what,got,want,(got==want)?"PASS":"FAIL");
}

ITG erosion_selftest(void)
{
  ITG nbad=0,i;
  /* four C3D4 elements, one integration point each, materials 1 and 2 */
  char lakon[4*8+1]="C3D4    C3D4    C3D4    C3D4    ";
  char matname[2*80+1];
  ITG ipkon[4]={0,4,8,12};
  ITG kon[16]={1,2,3,4, 2,3,4,5, 3,4,5,6, 4,5,6,7};
  ITG ielmat[4]={1,1,2,2};
  /* Three materials, not two: the test asks about material 3 (the one
     with no damage record), and both tables are indexed by the material
     number - ndmcon at 2*(imat-1) and dmcon at 1+(ndmat+1)*ntmat*(imat-1),
     which is 11 for imat=3.  Sized for two, the question read off the end
     of both arrays; -Wall said so and it was right. */
  ITG ndmcon[6]={4,0,4,0,0,0};      /* nconst per material              */
  double dmcon[16];                 /* off=1+(ndmat+1)*ntmat*(imat-1)   */
  double dam[4],visc[4],dmax,vmin;
  double tv[4]; ITG tip[4];
  erosion_policy pol;
  erosion_batch bat;

  memset(matname,' ',sizeof(matname)-1); matname[sizeof(matname)-1]=0;
  memcpy(matname,"ZrH",3);            /* material 1 */
  memcpy(matname+80,"Zr",2);          /* material 2 - a prefix of the first */
  for(i=0;i<16;i++) dmcon[i]=0.;
  dmcon[1]=1.;                        /* material 1: Rice-Tracey, type 1 */
  dmcon[6]=1.;                        /* material 2: same              */

  erosion_chki("progressive material is recognised",
               damage_progressive_material(1,ndmcon,dmcon,4,1),1,&nbad);
  erosion_chki("a material with no damage record is not",
               damage_progressive_material(3,ndmcon,dmcon,4,1),0,&nbad);

  /* 1. the filter.  "Zr" must take material 2 and leave material 1. */
  erosion_chki("filter unset takes every material",
               erosion_material_allowed(1,matname,NULL),1,&nbad);
  erosion_chki("filter ALL takes every material",
               erosion_material_allowed(1,matname,"ALL"),1,&nbad);
  erosion_chki("filter NONE takes none",
               erosion_material_allowed(1,matname,"ALL")&&
               !erosion_material_allowed(1,matname,"NONE"),1,&nbad);
  erosion_chki("filter Zr does NOT match ZrH (prefix)",
               erosion_material_allowed(1,matname,"Zr"),0,&nbad);
  erosion_chki("filter Zr matches Zr",
               erosion_material_allowed(2,matname,"Zr"),1,&nbad);
  /* the other direction, and the one the length comparison in
     erosion_material_allowed() is the ONLY guard for: without it "ZrH"
     matches the first two characters of "Zr" and the wrong phase erodes.
     Removing that comparison was tried; this is the line that went red. */
  erosion_chki("filter ZrH does NOT match Zr (the length guard)",
               erosion_material_allowed(2,matname,"ZrH"),0,&nbad);
  erosion_chki("filter is case insensitive",
               erosion_material_allowed(1,matname,"zrh"),1,&nbad);
  erosion_chki("a list takes any member",
               erosion_material_allowed(2,matname,"Steel,Zr,Cu"),1,&nbad);
  erosion_chki("a list refuses a non-member",
               erosion_material_allowed(1,matname,"Steel,Zr,Cu"),0,&nbad);

  /* 2. D against Dvis.  dam carries 1+D; visc carries D. */
  for(i=0;i<4;i++){ ipkon[i]=4*i; dam[i]=1.; visc[i]=0.; tv[i]=-1.; tip[i]=0; }
  dam[0]=1.9995; visc[0]=0.990;     /* past on D, short on Dvis */
  erosion_chki("trigger on Dvis leaves an element short on Dvis",
               erosion_mark_terminal(dam,ipkon,lakon,ielmat,1,ndmcon,dmcon,4,1,
                                     4,1,0.999,64,&dmax,tv,tip,matname,NULL,
                                     visc,1,&vmin),0,&nbad);
  erosion_chki("...and it is still assembled",(ipkon[0]>=0),1,&nbad);
  erosion_chki("trigger on D takes the same element",
               erosion_mark_terminal(dam,ipkon,lakon,ielmat,1,ndmcon,dmcon,4,1,
                                     4,1,0.999,64,&dmax,tv,tip,matname,NULL,
                                     visc,0,&vmin),1,&nbad);

  /* 4. the marking convention, and no double marking. */
  erosion_chki("marking is ipkon -> -ipkon-2, invertible",
               (ipkon[0]<0)&&((-ipkon[0]-2)==0),1,&nbad);
  erosion_chki("the trigger value is captured at the marking",
               (tv[0]>0.9994)&&(tv[0]<0.9996)&&(tip[0]==1),1,&nbad);
  erosion_chki("an element already marked is not marked again",
               erosion_mark_terminal(dam,ipkon,lakon,ielmat,1,ndmcon,dmcon,4,1,
                                     4,1,0.999,64,&dmax,tv,tip,matname,NULL,
                                     visc,0,&vmin),0,&nbad);

  /* 3. the batch budget, through the entry point the solver uses. */
  for(i=0;i<4;i++){ ipkon[i]=4*i; dam[i]=2.0; visc[i]=1.0; tv[i]=-1.; tip[i]=0; }
  pol.delete_d=0.999; pol.delete_visc=1; pol.filter=NULL;
  pol.deadall_g=0.; pol.deadsole_g=0.; pol.batchmax=2;
  erosion_batch_init(&bat);
  erosion_chki("four candidates, a budget of two, two leave",
               erosion_mark(&pol,&bat,dam,visc,ipkon,lakon,kon,ielmat,matname,
                            ndmcon,dmcon,4,1,8,4,4,1,1,tv,tip,1,0.),2,&nbad);
  erosion_chki("...and they are the first two, in element order",
               (ipkon[0]<0)&&(ipkon[1]<0)&&(ipkon[2]>=0)&&(ipkon[3]>=0),
               1,&nbad);
  erosion_chki("the rest are reported as terminal, not as a load-path rule",
               (bat.terminal==2)&&(bat.deadall==0)&&(bat.deadsole==0),1,&nbad);

  /* the filter, applied through the entry point: material 2 only. */
  for(i=0;i<4;i++){ ipkon[i]=4*i; tv[i]=-1.; tip[i]=0; }
  pol.filter="Zr"; pol.batchmax=64;
  erosion_batch_init(&bat);
  erosion_chki("with filter Zr only the two Zr elements leave",
               erosion_mark(&pol,&bat,dam,visc,ipkon,lakon,kon,ielmat,matname,
                            ndmcon,dmcon,4,1,8,4,4,1,1,tv,tip,1,0.),2,&nbad);
  erosion_chki("...and they are elements 3 and 4",
               (ipkon[0]>=0)&&(ipkon[1]>=0)&&(ipkon[2]<0)&&(ipkon[3]<0),
               1,&nbad);

  printf("[EROSION] self test %s (%" ITGFORMAT " failure(s))\n",
         nbad?"FAILED":"PASSED",nbad);
  fflush(stdout);
  return nbad;
}
