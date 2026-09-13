/*     The census: how far the assembled nodal stiffness has fallen.
 *
 *     This computation appeared in nonlingeo.c THREE TIMES, character for
 *     character - 29 lines each, at lines 13076, 14074 and 14572 of a 16000
 *     line file.  It was written three times because it had no owner: each
 *     time a new place needed the answer, the loop was pasted rather than
 *     called, and nothing said where the answer was supposed to come from.
 *
 *     That is the whole argument for this file.  A quantity with an owner is
 *     computed once and asked for; a quantity without one is recomputed
 *     wherever somebody happens to need it, and the copies drift.  These
 *     three had not drifted yet.  The next one would have.
 *
 *     WHAT IT MEASURES.  For every node, the ratio of its assembled diagonal
 *     to its OWN first positive value - so it is dimensionless, needs no
 *     mesh-wide median, and means the same thing on any deck.  Counts below
 *     1e-3, 1e-2 and 1e-1, the number of non-positive diagonals, and the
 *     extreme.  02-DIAGNOSTICS.md section 3 reads it: 3.7e-07 is a node that
 *     has lost ALL its bulk and hangs on failed facets; 0.04 has lost its
 *     bulk but its facets are intact; 0.1 to 0.3 still has live elements.
 *     The number tells you the topology.
 *
 *     It decides nothing.  damstate.c owns the JUDGEMENT of what has lost
 *     its load path; this owns the MEASUREMENT that the judgement is read
 *     against, and the two are deliberately separate objects.
 */
#include <stdio.h>
#include <stdlib.h>
#include "CalculiX.h"

/* Take the census.  Pure: it reads three arrays and writes one struct.
   addok may be NULL, in which case every node is considered - which is what
   the three original copies did when damage_addok had not been allocated. */
void stiffcensus_take(stiffcensus *c,ITG nk,const double *addiag,
                      const double *addiag0,const ITG *addok){
  ITG i;
  double r;
  c->below1=0;c->below2=0;c->below3=0;c->nonpositive=0;
  c->worst=-1;c->worstratio=2.;c->nnode=0;
  if((addiag==NULL)||(addiag0==NULL)) return;
  for(i=0;i<nk;i++){
    if((addok!=NULL)&&(addok[i]==0)) continue;
    if(addiag0[i]<=0.) continue;
    c->nnode++;
    if(addiag[i]<=0.){c->nonpositive++;continue;}
    r=addiag[i]/addiag0[i];
    if(r<1.e-3) c->below1++;
    if(r<1.e-2) c->below2++;
    if(r<1.e-1) c->below3++;
    if(r<c->worstratio){c->worstratio=r;c->worst=i+1;}
  }
}

/* Self test.  What has to be established is that the three thresholds are
   STRICT and nested, that a non-positive diagonal is counted and NOT
   compared, that a node with no intact reference is skipped entirely, and
   that the extreme is the extreme.  Every one of those is a line the three
   copies had in common and that a fourth copy could have got wrong. */
ITG stiffcensus_selftest(void){
  ITG bad=0;
  stiffcensus c;
  /* six nodes: one healthy, one at each threshold exactly, one below all
     three, one non-positive, one with no intact reference */
  double d0[6]={1.,1.,1.,1.,1.,0.};
  double d [6]={1.,1.e-3,1.e-2,1.e-4,-5.,1.};
  ITG    ok[6]={1,1,1,1,1,1};

  stiffcensus_take(&c,6,d,d0,ok);
  /* 1.e-3 is NOT below 1.e-3: the comparison is strict, and a node exactly
     at the threshold has repeatedly been the one under discussion */
  if(c.below1!=1){
    printf("[CENSUS] *ERROR: below_1e-3 is %" ITGFORMAT ", want 1 "
           "(is the comparison strict?)\n",c.below1);bad++;}
  if(c.below2!=2){
    printf("[CENSUS] *ERROR: below_1e-2 is %" ITGFORMAT ", want 2\n",
           c.below2);bad++;}
  if(c.below3!=3){
    printf("[CENSUS] *ERROR: below_1e-1 is %" ITGFORMAT ", want 3\n",
           c.below3);bad++;}
  if(c.nonpositive!=1){
    printf("[CENSUS] *ERROR: nonpositive is %" ITGFORMAT ", want 1\n",
           c.nonpositive);bad++;}
  if(c.worst!=4){
    printf("[CENSUS] *ERROR: the worst node is %" ITGFORMAT ", want 4\n",
           c.worst);bad++;}
  if(c.worstratio!=1.e-4){
    printf("[CENSUS] *ERROR: the worst ratio is %.4e, want 1.0000e-04\n",
           c.worstratio);bad++;}
  /* the node with no intact reference is not a node, and a non-positive
     diagonal must not contribute a ratio */
  if(c.nnode!=5){
    printf("[CENSUS] *ERROR: %" ITGFORMAT " nodes counted, want 5 "
           "(is a node with no intact reference being counted?)\n",
           c.nnode);bad++;}

  /* a masked node is not counted at all, however bad it is: masking the
     worst node must remove it from every count AND hand the extreme to the
     next node down */
  ok[3]=0;
  stiffcensus_take(&c,6,d,d0,ok);
  if((c.below1!=0)||(c.below2!=1)||(c.below3!=2)||(c.worst!=2)||(c.nnode!=4)){
    printf("[CENSUS] *ERROR: a masked node still reached the census "
           "(below_1e-3=%" ITGFORMAT " below_1e-2=%" ITGFORMAT
           " below_1e-1=%" ITGFORMAT " worst=%" ITGFORMAT " nodes=%"
           ITGFORMAT ")\n",c.below1,c.below2,c.below3,c.worst,c.nnode);bad++;}

  /* and with nothing to measure, the census is empty rather than wrong */
  stiffcensus_take(&c,6,NULL,NULL,NULL);
  if((c.nnode!=0)||(c.worst!=-1)||(c.below1!=0)){
    printf("[CENSUS] *ERROR: an empty census is not empty\n");bad++;}
  return bad;
}
