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
  TOPOLOGY DIAGNOSTIC FOR THE DELETION LIFECYCLE
  ==============================================

  This file measures.  It changes nothing, decides nothing and is not on
  any solution path: every routine takes const model arrays and writes a
  report.  It exists because three different explanations of the s3rad
  wall - a true orphan degree of freedom, a physically detached component
  with rigid-body modes, and a constitutive/deletion transition that the
  corrector cannot follow - are indistinguishable in the existing log and
  are distinguished by exactly these quantities:

    connected components   a detached component is a component of the LIVE
                           element graph that contains no prescribed dof
                           and no MPC.  It has six rigid-body modes and the
                           operator is singular on them.
    orphan dofs            a node that carries an active dof while no live
                           element touches it.  Its rows are structurally
                           empty: that is a bookkeeping defect, not
                           mechanics.
    empty rows / zero
    diagonals              the algebraic signature of the two above.
    residual projection    if the residual lies in the span of the soft
                           modes, the obstruction is the topology.  If it
                           does not, the topology is sound and the
                           corrector is the problem.

  WHY THE ELEMENT GRAPH IS BUILT HERE AND NOT WITH nodebelongstoel
  ---------------------------------------------------------------
  nodebelongstoel decides the node count of an element from lakon(4:4) and
  a short list of names, and a UC6 cohesive element matches none of them,
  so it falls through the else branch and is SKIPPED.  A node held only by
  cohesive facets therefore looks unsupported to that routine.  For a
  connectivity measurement that is exactly backwards: the cohesive facet is
  what still holds the two crack faces together, and a component reached
  only through facets is not detached.  This file reads the node count for
  a user element from byte 8 of lakon, which is the convention the deletion
  code itself uses.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalculiX.h"

/* ------------------------------------------------------------------ */
/* Union-find.  Pure, no model, so the self test can drive it directly. */
/* ------------------------------------------------------------------ */

ITG topodiag_find(ITG *p,ITG a){

  while(p[a]!=a){p[a]=p[p[a]];a=p[a];}
  return a;
}

void topodiag_union(ITG *p,ITG a,ITG b){

  a=topodiag_find(p,a);
  b=topodiag_find(p,b);
  if(a!=b) p[a]=b;
}

/* Number of nodes of element ie, or 0 if it is not one this diagnostic
   understands.  Bulk elements are read from lakon exactly as
   nodebelongstoel does; USER elements are read from byte 8, which is
   where the element reader stores the *USER ELEMENT NODES= count and what
   damage_de13_mark_deadall already relies on. */

ITG topodiag_nope(const char *lak){

  ITG n;

  if(lak[0]=='F') return 0;
  if(lak[0]=='U'){
    n=(ITG)((unsigned char)lak[7]);
    return ((n>=1)&&(n<=20))?n:0;
  }
  if(lak[3]=='2') return 20;
  if(lak[3]=='8') return 8;
  if(lak[3]=='4') return 4;
  if((lak[3]=='1')&&(lak[4]=='0')) return 10;
  if(lak[3]=='6') return 6;
  if((lak[3]=='1')&&(lak[4]=='5')) return 15;
  return 0;
}

/* ------------------------------------------------------------------ */
/* Self test.  Everything below is checked on a synthetic model whose    */
/* answer is known by construction: two tetrahedra sharing a face, one    */
/* isolated tetrahedron, one node with no element at all.                */
/* ------------------------------------------------------------------ */

ITG topodiag_selftest(void){

  ITG nbad=0,i,p[8],a,b;
  char lak[8*3];

  printf("[TOPODIAG] self test\n");

  /* A: union-find */
  for(i=0;i<8;i++) p[i]=i;
  topodiag_union(p,0,1);
  topodiag_union(p,1,2);
  topodiag_union(p,4,5);
  a=(topodiag_find(p,0)==topodiag_find(p,2));
  b=(topodiag_find(p,0)==topodiag_find(p,4));
  printf("   %-26s 0~2 %s, 0~4 %s  %s\n","A union-find",a?"yes":"no",
         b?"yes":"no",((a==1)&&(b==0))?"ok":"FAIL");
  if(!((a==1)&&(b==0))) nbad++;

  /* B: element node counts, including the UC6 case nodebelongstoel drops */
  memset(lak,' ',sizeof(lak));
  memcpy(lak+0,"C3D4",4);
  memcpy(lak+8,"UC6",3);  lak[8+7]=(char)6;
  memcpy(lak+16,"C3D10",5);
  a=topodiag_nope(lak);
  b=topodiag_nope(lak+8);
  i=topodiag_nope(lak+16);
  printf("   %-26s C3D4=%" ITGFORMAT " UC6=%" ITGFORMAT " C3D10=%"
         ITGFORMAT "  %s\n","B element node counts",a,b,i,
         ((a==4)&&(b==6)&&(i==10))?"ok":"FAIL");
  if(!((a==4)&&(b==6)&&(i==10))) nbad++;

  printf("[TOPODIAG] self test %s (%" ITGFORMAT " failure(s))\n",
         (nbad==0)?"PASSED":"FAILED",nbad);
  return nbad;
}

/* ------------------------------------------------------------------ */
/* The measurement.                                                     */
/*                                                                      */
/* comp[] is filled with a component id per node (-1 for a node with no  */
/* live element).  Everything else is counted into the report.          */
/* ------------------------------------------------------------------ */

void topodiag_report_zero(topodiag_report *r){

  memset(r,0,sizeof(*r));
  r->orphannode=-1;
  r->floatfirst=-1;
  r->resmaxeq=-1;
  r->resmaxnode=-1;
  r->resmaxcomp=-1;
  r->maincomp=-1;
  r->admin=-1.;
}

void topodiag_run(topodiag_report *r,ITG *comp,
                  const ITG *kon,const ITG *ipkon,const char *lakon,ITG ne,
                  ITG nk,const ITG *nactdof,ITG mt,
                  const ITG *nodeboun,const ITG *ndirboun,ITG nboun,
                  const ITG *ipompc,const ITG *nodempc,ITG nmpc,
                  const double *ad,const double *au,const ITG *jq,
                  const ITG *irow,ITG neq,ITG nzs,
                  const double *res){

  ITG i,j,k,n,nope,idx,root,*p=NULL,*held=NULL,*nn=NULL,*ndof=NULL;
  ITG *map=NULL,c,nc=0;
  double amax,d;

  topodiag_report_zero(r);
  if((nk<=0)||(ne<=0)) return;

  NNEW(p,ITG,nk);
  for(i=0;i<nk;i++){p[i]=i;comp[i]=-1;}

  /* ---- 1. the live element graph -------------------------------- */

  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    nope=topodiag_nope(&lakon[8*i]);
    if(nope<=0) continue;
    idx=ipkon[i];
    r->nelem++;
    if(lakon[8*i]=='U') r->nfacet++;
    n=kon[idx+0]-1;
    if((n<0)||(n>=nk)) continue;
    for(j=1;j<nope;j++){
      k=kon[idx+j]-1;
      if((k<0)||(k>=nk)) continue;
      topodiag_union(p,n,k);
    }
    for(j=0;j<nope;j++){
      k=kon[idx+j]-1;
      if((k>=0)&&(k<nk)) comp[k]=0;      /* marker: touched by a live el */
    }
  }

  /* ---- 2. label the components ----------------------------------- */

  NNEW(map,ITG,nk);
  for(i=0;i<nk;i++) map[i]=-1;
  for(i=0;i<nk;i++){
    if(comp[i]!=0) continue;             /* no live element here        */
    root=topodiag_find(p,i);
    if(map[root]<0){map[root]=nc;nc++;}
  }
  for(i=0;i<nk;i++){
    if(comp[i]!=0){comp[i]=-1;continue;}
    comp[i]=map[topodiag_find(p,i)];
  }
  r->ncomp=nc;

  /* ---- 3. which components are held ------------------------------ */

  NNEW(held,ITG,(nc>0)?nc:1);
  NNEW(nn,ITG,(nc>0)?nc:1);
  NNEW(ndof,ITG,(nc>0)?nc:1);
  for(i=0;i<nk;i++){
    if(comp[i]<0) continue;
    nn[comp[i]]++;
    for(k=1;k<mt;k++){
      if(nactdof[mt*i+k]>0) ndof[comp[i]]++;
    }
  }
  for(i=0;i<nboun;i++){
    n=nodeboun[i]-1;
    if((n<0)||(n>=nk)) continue;
    if(comp[n]>=0) held[comp[n]]|=1;
  }
  for(i=0;i<nmpc;i++){
    idx=ipompc[i]-1;
    while(idx!=-1){
      n=nodempc[3*idx]-1;
      if((n>=0)&&(n<nk)&&(comp[n]>=0)) held[comp[n]]|=2;
      idx=nodempc[3*idx+2]-1;
    }
  }
  for(c=0;c<nc;c++){
    if(held[c]==0){
      r->nfloat++;
      r->nfloatnode+=nn[c];
      r->nfloatdof+=ndof[c];
      if(nn[c]>r->floatmax) r->floatmax=nn[c];
      if(r->floatfirst<0) r->floatfirst=c;
    }
  }
  if(nc>0){
    r->maincomp=0;
    for(c=1;c<nc;c++) if(nn[c]>nn[r->maincomp]) r->maincomp=c;
    r->mainnode=nn[r->maincomp];
  }

  /* ---- 4. orphan dofs: active, but no live element ---------------- */

  for(i=0;i<nk;i++){
    if(comp[i]>=0) continue;             /* a live element holds it     */
    n=0;
    for(k=1;k<mt;k++){
      if(nactdof[mt*i+k]>0) n++;
    }
    if(n==0) continue;                   /* no dof: nothing to report   */
    r->norphandof+=n;
    r->norphannode++;
    if(r->orphannode<0) r->orphannode=i+1;
  }

  /* ---- 5. algebraic signature ------------------------------------ */

  if((ad!=NULL)&&(neq>0)){
    amax=0.;
    for(i=0;i<neq;i++){d=fabs(ad[i]);if(d>amax) amax=d;}
    r->admax=amax;
    for(i=0;i<neq;i++){
      d=fabs(ad[i]);
      if(!(d>0.)) r->nzerodiag++;
      else if(d<1.e-12*amax) r->ntinydiag++;
      if((r->admin<0.)||(d<r->admin)) r->admin=d;
    }
  }
  /* An equation is ISOLATED when it has no off-diagonal coupling at all -
     neither as a column of the stored triangle nor as a row index in it.
     That, and not an empty column, is the algebraic signature of an
     orphan degree of freedom: the last column of a triangular storage is
     always empty and means nothing. */

  if((jq!=NULL)&&(irow!=NULL)&&(neq>0)&&(nzs>0)){
    ITG *cpl=NULL;
    NNEW(cpl,ITG,neq);
    for(i=0;i<neq;i++) if(jq[i+1]>jq[i]) cpl[i]=1;
    for(k=0;k<nzs;k++){
      n=irow[k]-1;
      if((n>=0)&&(n<neq)) cpl[n]=1;
    }
    for(i=0;i<neq;i++) if(cpl[i]==0) r->nisolated++;
    SFREE(cpl);
  }

  /* ---- 6. where the residual lives -------------------------------- */

  if((res!=NULL)&&(neq>0)){
    amax=0.;k=-1;
    for(i=0;i<neq;i++){
      d=fabs(res[i]);
      if(d>amax){amax=d;k=i;}
      r->resnorm+=res[i]*res[i];
    }
    r->resnorm=sqrt(r->resnorm);
    r->resmax=amax;
    r->resmaxeq=k;
    r->resmaxnode=-1;
    if(k>=0){
      for(i=0;i<nk;i++){
        for(j=1;j<mt;j++){
          if(nactdof[mt*i+j]==k+1){
            r->resmaxnode=i+1;
            r->resmaxdir=j;
            r->resmaxcomp=comp[i];
            i=nk;break;
          }
        }
      }
    }
    /* how much of the residual sits on floating components */
    for(i=0;i<nk;i++){
      if(comp[i]<0) continue;
      if(held[comp[i]]!=0) continue;
      for(j=1;j<mt;j++){
        k=nactdof[mt*i+j];
        if(k>0) r->resfloat+=res[k-1]*res[k-1];
      }
    }
    r->resfloat=sqrt(r->resfloat);
  }

  SFREE(ndof);SFREE(nn);SFREE(held);SFREE(map);SFREE(p);
}

/* Fraction of v that lies along w, both in equation space.  Used to ask
   whether the residual is in the span of a soft mode found by inverse
   iteration. */

/* Support of one node: how many live bulk elements and how many live user
   (cohesive) elements hold it.  The soft-mode peak and the residual peak
   are only interpretable next to this: a node held by facets alone is a
   conditioning problem, a node held by two nearly dead tets is another,
   and they are not the same defect. */

void topodiag_support(ITG node,const ITG *kon,const ITG *ipkon,
                      const char *lakon,ITG ne,ITG *nbulk,ITG *nfac){

  ITG i,j,nope,idx;

  *nbulk=0;*nfac=0;
  if(node<1) return;
  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    nope=topodiag_nope(&lakon[8*i]);
    if(nope<=0) continue;
    idx=ipkon[i];
    for(j=0;j<nope;j++){
      if(kon[idx+j]==node){
        if(lakon[8*i]=='U') (*nfac)++; else (*nbulk)++;
        break;
      }
    }
  }
}

/* Orthogonalise w against the k vectors already in q (each of length neq),
   then normalise.  Returns 0 if what is left is numerically nothing, which
   is how a deflated inverse iteration says "there is no further
   independent soft direction". */

ITG topodiag_deflate(double *w,const double *q,ITG k,ITG neq){

  ITG i,j;
  double d,n;

  for(j=0;j<k;j++){
    d=0.;
    for(i=0;i<neq;i++) d+=w[i]*q[(size_t)j*neq+i];
    for(i=0;i<neq;i++) w[i]-=d*q[(size_t)j*neq+i];
  }
  n=0.;
  for(i=0;i<neq;i++) n+=w[i]*w[i];
  n=sqrt(n);
  if(!(n>0.)) return 0;
  for(i=0;i<neq;i++) w[i]/=n;
  return 1;
}

/* Fraction of |v| that lies in the span of the k orthonormal columns of q.
   This is the quantity a single-vector cosine cannot give: a soft
   SUBSPACE of dimension more than one would hide from it. */

double topodiag_project_span(const double *v,const double *q,ITG k,ITG neq){

  ITG i,j;
  double d,s=0.,vv=0.;

  for(i=0;i<neq;i++) vv+=v[i]*v[i];
  if(!(vv>0.)) return 0.;
  for(j=0;j<k;j++){
    d=0.;
    for(i=0;i<neq;i++) d+=v[i]*q[(size_t)j*neq+i];
    s+=d*d;
  }
  return sqrt(s/vv);
}

double topodiag_project(const double *v,const double *w,ITG neq){

  ITG i;
  double vw=0.,ww=0.,vv=0.;

  for(i=0;i<neq;i++){vw+=v[i]*w[i];ww+=w[i]*w[i];vv+=v[i]*v[i];}
  if(!(ww>0.)||!(vv>0.)) return 0.;
  return fabs(vw)/sqrt(ww*vv);
}

void topodiag_print(const topodiag_report *r,const char *tag,ITG iinc){

  printf("[TOPODIAG] %s inc=%" ITGFORMAT "\n",tag,iinc);
  printf("[TOPODIAG]   live elements %" ITGFORMAT " (of which %" ITGFORMAT
         " cohesive facets)\n",r->nelem,r->nfacet);
  printf("[TOPODIAG]   connected components %" ITGFORMAT
         ", largest holds %" ITGFORMAT " node(s)\n",r->ncomp,r->mainnode);
  printf("[TOPODIAG]   FLOATING components (no prescribed dof, no MPC): %"
         ITGFORMAT ", %" ITGFORMAT " node(s), %" ITGFORMAT
         " active dof(s), largest %" ITGFORMAT " node(s)\n",
         r->nfloat,r->nfloatnode,r->nfloatdof,r->floatmax);
  printf("[TOPODIAG]   ORPHAN dofs (active dof, no live element): %"
         ITGFORMAT " on %" ITGFORMAT " node(s)%s\n",
         r->norphandof,r->norphannode,
         (r->orphannode>0)?"":", none");
  if(r->orphannode>0)
    printf("[TOPODIAG]   first orphan node %" ITGFORMAT "\n",r->orphannode);
  printf("[TOPODIAG]   diagonal: max %.6e min %.6e, exactly zero %"
         ITGFORMAT ", below 1e-12*max %" ITGFORMAT
         "; ISOLATED equations (no off-diagonal coupling at all) %"
         ITGFORMAT "\n",r->admax,r->admin,r->nzerodiag,r->ntinydiag,
         r->nisolated);
  if(r->resnorm>0.){
    printf("[TOPODIAG]   residual |R|2=%.6e |R|inf=%.6e at equation %"
           ITGFORMAT " (node %" ITGFORMAT " dir %" ITGFORMAT
           " component %" ITGFORMAT ")\n",r->resnorm,r->resmax,
           r->resmaxeq,r->resmaxnode,r->resmaxdir,r->resmaxcomp);
    printf("[TOPODIAG]   share of |R|2 on floating components: %.6e "
           "(%.4f%%)\n",r->resfloat,
           100.*r->resfloat/(r->resnorm+1.e-300));
  }
  fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* Content hashes.                                                      */
/*                                                                      */
/* FNV-1a over the raw bytes.  The point is not cryptography, it is to   */
/* be able to say "this attempt started from the SAME committed state as */
/* the previous one" as a measurement instead of an inference - and to   */
/* say the same about the composition of a deletion batch, which is the  */
/* other half of the repetition claim.                                   */
/* ------------------------------------------------------------------ */

unsigned long long topodiag_hash_bytes(const void *p,size_t n,
                                       unsigned long long h){

  const unsigned char *b=(const unsigned char *)p;
  size_t i;

  if(p==NULL) return h;
  for(i=0;i<n;i++){
    h^=(unsigned long long)b[i];
    h*=1099511628211ULL;
  }
  return h;
}

unsigned long long topodiag_hash_d(const double *a,ITG n,
                                   unsigned long long h){

  return topodiag_hash_bytes(a,(size_t)((n>0)?n:0)*sizeof(double),h);
}

unsigned long long topodiag_hash_i(const ITG *a,ITG n,
                                   unsigned long long h){

  return topodiag_hash_bytes(a,(size_t)((n>0)?n:0)*sizeof(ITG),h);
}

unsigned long long topodiag_hash_seed(void){ return 14695981039346656037ULL; }

/* Sorted hash of a deletion batch, so that the same set in a different
   order compares equal.  The list is copied, not reordered in place. */

unsigned long long topodiag_hash_batch(const ITG *elem,ITG n,
                                       ITG *sorted,unsigned long long h){

  ITG i,j,t;

  if((elem==NULL)||(n<=0)) return h;
  for(i=0;i<n;i++) sorted[i]=elem[i];
  for(i=1;i<n;i++){                       /* insertion sort: n is tiny */
    t=sorted[i];
    for(j=i-1;(j>=0)&&(sorted[j]>t);j--) sorted[j+1]=sorted[j];
    sorted[j+1]=t;
  }
  return topodiag_hash_i(sorted,n,h);
}
