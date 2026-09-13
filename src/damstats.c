/*     CalculiX - damage/fracture extension                              */
/*     damstats.c: what the damage state looks like, reported.           */

/* Why this module exists
   ----------------------
   Three functions that report the damage field - the per-element census,
   the appended per-increment record and the VTK series - plus the eight
   counters they fill, were inside nonlingeo().  None of them decides
   anything: they read dam, they count, and they write files that nothing
   in the solver reads back.

   Contract
   --------
     - damstats_element() is the census: for each live element the maximum
       degradation over its active integration points, and the counts that
       follow from it.  Pure;
     - damstats_append() and damstats_write_vtk() write.  They open, they
       write, they close, and the solver never reads what they wrote;
     - damstats is the eight counters as one object, so that "what did the
       last census say" has one answer instead of eight.

   The VTK writer's connectivity is worth the comment it carries in the
   body: an earlier version ran j=1..4 over a four-node element and
   silently produced a mesh that was wrong by one node everywhere.     */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void damstats_init(damstats *d){ memset(d,0,sizeof(*d)); }

/* Exact DE1 element statistics.
   For each element the maximum degradation over its active integration
   points is used.  In the present DE1 implementation only C3D4 is enabled,
   therefore this is exactly the single integration-point value. */
void damstats_element(const double *dam,const double *damold,
                             const ITG *ipkon,const char *lakon,
                             ITG ne0,ITG mi0,
                             ITG *nactive,ITG *ngt01,ITG *ngt05,
                             ITG *ngt09,ITG *nfull,ITG *nchanged,
                             double *dmax,double *maxdelta)
{
  ITG i,j,nip;
  double de,dold,delta,demax,delmax;

  *nactive=0;*ngt01=0;*ngt05=0;*ngt09=0;*nfull=0;*nchanged=0;
  *dmax=0.;*maxdelta=0.;

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='C') continue;

    nip=topo_element_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;

    demax=0.;delmax=0.;
    for(j=0;j<nip;j++){
      de=dam[mi0*i+j]-1.;
      if(de<0.) de=0.;
      if(de>1.) de=1.;
      if(de>demax) demax=de;

      if(damold!=NULL){
        dold=damold[mi0*i+j]-1.;
        if(dold<0.) dold=0.;
        if(dold>1.) dold=1.;
        delta=fabs(de-dold);
        if(delta>delmax) delmax=delta;
      }
    }

    if(demax>0.){
      (*nactive)++;
      if(demax>0.1) (*ngt01)++;
      if(demax>0.5) (*ngt05)++;
      if(demax>0.9) (*ngt09)++;
      if(demax>=0.999) (*nfull)++;
      if(demax>*dmax) *dmax=demax;
    }

    if(delmax>DAMAGE_DE1_FP_TOL) (*nchanged)++;
    if(delmax>*maxdelta) *maxdelta=delmax;
  }
}

/* Append one compact accepted-state record.  This file is intentionally
   independent of the legacy .damage hard-deletion history. */
void damstats_append(const char *jobnamec,ITG istep,ITG iinc,
                                    double steptime,double totaltime,
                                    ITG passes,ITG nactive,ITG ngt01,
                                    ITG ngt05,ITG ngt09,ITG nfull,
                                    double dmax,double maxdelta)
{
  char fname[200]="";
  FILE *f=NULL;
  static ITG de1stats_initialized=0;

  strcpy2(fname,jobnamec,132);
  strcat(fname,".de1stats");

  /* DE1.3.1: one solver process == one fresh diagnostic history.
     The first accepted damage state truncates stale data from an older run;
     subsequent accepted states append normally. */
  if(de1stats_initialized==0){
    f=fopen(fname,"w");
  }else{
    f=fopen(fname,"a");
  }
  if(f==NULL) return;

  if(de1stats_initialized==0){
    fprintf(f,"# CalculiX DE1.3.1 accepted damage states\n");
    fprintf(f,"# step increment step_time total_time passes active "
              "Dgt0.1 Dgt0.5 Dgt0.9 Dfull Dmax max_dD\n");
    de1stats_initialized=1;
  }

  fprintf(f,"%" ITGFORMAT " %" ITGFORMAT " %.15e %.15e "
            "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " "
            "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " "
            "%.15e %.15e\n",
          istep,iinc,steptime,totaltime,passes,nactive,ngt01,ngt05,
          ngt09,nfull,dmax,maxdelta);
  fclose(f);
}

/* Write the latest accepted DE1 state as an exact element-cell VTK snapshot.
   The file is overwritten, so long calculations do not accumulate large
   post-processing files.  All active C3D4 cells are written.  DE1_D and
   DUCT_IP are CELL_DATA taken directly from the solver integration-point
   history; no extrapolation or nodal averaging is involved.  Coordinates
   are written in the current deformed configuration. */
void damstats_write_vtk(const char *jobnamec,
                                 const double *co,const double *vold,
                                 ITG nk,ITG mt,const ITG *kon,
                                 const ITG *ipkon,const char *lakon,
                                 const ITG *ielmat,ITG mi2,
                                 const double *dam,ITG mi0,ITG ne0,
                                 ITG istep,ITG iinc,double steptime)
{
  char fname[200]="",seq[32]="";
  FILE *f=NULL;
  ITG i,j,indexe,node,ncell=0;
  double d,duct,x,y,z;
  static ITG vtkseries=-1,vtkcount=0;

  /* CCX_DAMAGE_VTK_SERIES keeps every accepted state as
     <job>.de1.NNNNN.vtk instead of overwriting a single frame.  Needed to
     show a sequence of damage rather than one final picture. */

  if(vtkseries<0){
    vtkseries=(ccxopt_getenv("CCX_DAMAGE_VTK_SERIES")!=NULL)?1:0;
  }

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)==0) ncell++;
  }

  strcpy2(fname,jobnamec,132);
  if(vtkseries==1){
    sprintf(seq,".de1.%5.5d.vtk",(int)vtkcount);
    strcat(fname,seq);
    vtkcount++;
  }else{
    strcat(fname,".de1.vtk");
  }
  f=fopen(fname,"w");
  if(f==NULL) return;

  fprintf(f,"# vtk DataFile Version 3.0\n");
  fprintf(f,"CalculiX DE1.3.1 exact cell damage step=%" ITGFORMAT
            " inc=%" ITGFORMAT " time=%.12e\n",istep,iinc,steptime);
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET UNSTRUCTURED_GRID\n");

  fprintf(f,"POINTS %" ITGFORMAT " double\n",nk);
  for(i=0;i<nk;i++){
    x=co[3*i];
    y=co[3*i+1];
    z=co[3*i+2];
    if(vold!=NULL){
      x+=vold[mt*i+1];
      y+=vold[mt*i+2];
      z+=vold[mt*i+3];
    }
    fprintf(f,"%.15e %.15e %.15e\n",x,y,z);
  }

  fprintf(f,"CELLS %" ITGFORMAT " %" ITGFORMAT "\n",ncell,5*ncell);
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    indexe=ipkon[i];
    fprintf(f,"4");
    for(j=0;j<4;j++){
      node=kon[indexe+j]-1;
      fprintf(f," %" ITGFORMAT,node);
    }
    fprintf(f,"\n");
  }

  fprintf(f,"CELL_TYPES %" ITGFORMAT "\n",ncell);
  for(i=0;i<ncell;i++) fprintf(f,"10\n");

  fprintf(f,"CELL_DATA %" ITGFORMAT "\n",ncell);

  fprintf(f,"SCALARS DE1_D double 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    d=dam[mi0*i]-1.;
    if(d<0.) d=0.;
    if(d>1.) d=1.;
    fprintf(f,"%.15e\n",d);
  }

  fprintf(f,"SCALARS DUCT_IP double 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    duct=dam[mi0*i];
    fprintf(f,"%.15e\n",duct);
  }

  fprintf(f,"SCALARS ELEMENT_ID int 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    fprintf(f,"%" ITGFORMAT "\n",i+1);
  }

  fprintf(f,"SCALARS MATERIAL_ID int 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    fprintf(f,"%" ITGFORMAT "\n",ielmat[mi2*i]);
  }

  fclose(f);
}
