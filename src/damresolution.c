/*     [RESOLUTION]: the resolution and scale levers of a damage deck, printed
 *     once at the start of a step.  Output only - nothing here is read back.
 *
 *     WHY.  The decks this module is used on are rescaled: specimen size,
 *     hydride size and morphology, material data.  Each of the mechanisms
 *     has a window it works in, stated in its own units - the internal
 *     length against the element size, the failure displacement against the
 *     element's elastic stretch, the cohesive zone against the facet size,
 *     the penalty against the bulk stiffness, the viscosities against the
 *     step - and a rescale moves the deck across those windows without any
 *     message.  This table puts every one of them on one screen as a
 *     dimensionless number with the bound it is judged against, so the
 *     lever that left its window is visible before the run, not after it.
 *
 *     The bounds are the literature's where there is one (Turon et al.
 *     2007 for the cohesive zone and the penalty; Bazant and Jirasek 2002
 *     for the resolved length) and the code's own where there is not (the
 *     soft resolution criterion in damnonlocal.f).  The groups follow the
 *     library note notes/scaling-levers.md.
 *
 *     ELEMENT SIZE h is (6V)^(1/3) for a linear tetrahedron, the legacy
 *     crack-band width, and the mean corner-edge length for every other
 *     family: a statistic for judging a window, not a width any law uses.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalculiX.h"

static int dres_cmp(const void *a,const void *b){
  double x=*(const double *)a,y=*(const double *)b;
  return (x<y)?-1:((x>y)?1:0);
}

/* element size of a volume element, see the header */
static double dres_h(const double *co,const ITG *kon,ITG indexe,
                     const char *lakonl){
  ITG nc,i,j,k,n[8],ne=0;
  double s=0.,d,a[3],b[3],c[3],det;
  if(lakonl[3]=='4'||(lakonl[3]=='1'&&lakonl[4]=='0')) nc=4;
  else if(lakonl[3]=='6'||(lakonl[3]=='1'&&lakonl[4]=='5')) nc=6;
  else if(lakonl[3]=='8'||lakonl[3]=='2') nc=8;
  else return 0.;
  for(i=0;i<nc;i++) n[i]=kon[indexe+i]-1;
  if(nc==4){
    for(k=0;k<3;k++){
      a[k]=co[3*n[1]+k]-co[3*n[0]+k];
      b[k]=co[3*n[2]+k]-co[3*n[0]+k];
      c[k]=co[3*n[3]+k]-co[3*n[0]+k];
    }
    det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
       +a[2]*(b[0]*c[1]-b[1]*c[0]);
    return pow(fabs(det),1./3.);
  }
  for(i=0;i<nc;i++){
    for(j=i+1;j<nc;j++){
      d=0.;
      for(k=0;k<3;k++) d+=pow(co[3*n[i]+k]-co[3*n[j]+k],2);
      s+=sqrt(d);ne++;
    }
  }
  /* the mean over ALL corner pairs overstates the edge by the diagonals;
     scale to the edge of a cube (12 edges, 12 face and 4 body diagonals:
     mean = (12+12*sqrt2+4*sqrt3)/28 of the edge) or a wedge (unit
     right prism: 15 pairs, close enough to the same factor) */
  return (ne>0)?(s/ne)/((12.+12.*sqrt(2.)+4.*sqrt(3.))/28.):0.;
}

void damage_resolution_report(const double *co,const ITG *kon,
     const ITG *ipkon,const char *lakon,ITG ne,const ITG *ielmat,
     const ITG *mi,const char *matname,ITG nmat,const double *elcon,
     const ITG *nelcon,ITG ncmat_,ITG ntmat_,const double *plicon,
     const ITG *nplicon,ITG npmat_,const ITG *ndmcon,const double *dmcon,
     ITG ndmat_,const ITG *ielprop,const double *prop,double dt0,
     double visceta){
  ITG i,m,n,cnt,kind,iok,mode,nres,nsoft,ncard;
  double *h=NULL,E,sy,uf,gf,ell,ellsave,ellmax,hmed,hmin,hmax,Ebulk=0.,
    hbulk=0.;
  ITG nbulkmax=0;
  char name[81];

  FORTRAN(damnonlocalget,(&ellsave));
  FORTRAN(damnlellmax,(&ellmax));
  FORTRAN(damnonlocalmodeget,(&mode));
  FORTRAN(damnlellcardcount,(&ncard));

  printf("\n[RESOLUTION] scale and resolution levers of this step "
         "(output only; bound in brackets, '!!' = outside it)\n");

  NNEW(h,double,(ne>0)?ne:1);
  for(m=1;m<=nmat;m++){
    if(ndmcon[2*(m-1)]<=0) continue;           /* no damage card */
    cnt=0;
    for(i=0;i<ne;i++){
      if(ipkon[i]<0) continue;
      if(lakon[8*i]!='C') continue;
      if(ielmat[mi[2]*i]!=m) continue;
      h[cnt]=dres_h(co,kon,ipkon[i],&lakon[8*i]);
      if(h[cnt]>0.) cnt++;
    }
    if(cnt==0) continue;
    qsort(h,cnt,sizeof(double),dres_cmp);
    hmin=h[0];hmax=h[cnt-1];hmed=h[cnt/2];

    strncpy(name,&matname[80*(m-1)],80);name[80]='\0';
    for(n=79;n>=0&&name[n]==' ';n--) name[n]='\0';

    /* isotropic elasticity: nelcon 2, or -50..-54 for the isotropic-
       elastic plasticity models, which keep E in the same slot */
    E=((nelcon[2*(m-1)]==2)||
       ((nelcon[2*(m-1)]<=-50)&&(nelcon[2*(m-1)]>=-54)))?
      elcon[(ncmat_+1)*ntmat_*(m-1)+1]:0.;
    sy=(nplicon[(ntmat_+1)*(m-1)]>0)?
      plicon[(2*npmat_+1)*ntmat_*(m-1)+1]:0.;
    if(cnt>nbulkmax){nbulkmax=cnt;Ebulk=E;hbulk=hmed;}

    printf("[RESOLUTION] material %s: %" ITGFORMAT " elements, h min/median/"
           "max = %.3e / %.3e / %.3e\n",name,cnt,hmin,hmed,hmax);

    /* the internal length this material is averaged with */
    ell=0.;
    if(ellsave>0.) ell=ellsave;
    else{
      FORTRAN(damnlellgetmat,(&m,&ell,&iok));
      if((iok==0)&&(ncard>0)) ell=ellmax;     /* inherits, see damnlellel */
    }
    if((ell>0.)&&(mode!=2)){
      nres=0;nsoft=0;
      for(i=0;i<cnt;i++){
        if(h[i]<=0.5*ell) nres++;              /* strict: h <= ell/2       */
        if(0.19*h[i]<=ell) nsoft++;            /* the code's soft criterion */
      }
      printf("[RESOLUTION]   ell = %.3e: ell/h(median) = %.3f; resolved "
             "h <= ell/2 for %.1f %% of elements [100]%s, soft ell >= 0.19h "
             "for %.1f %%\n",ell,ell/hmed,100.*nres/cnt,
             (nres<cnt)?" !!":"",100.*nsoft/cnt);
    }else if(ell>0.){
      printf("[RESOLUTION]   FROZEN control: no length is applied\n");
    }else{
      printf("[RESOLUTION]   local damage: the band width follows the mesh "
             "(crack band)\n");
    }

    /* failure displacement against the element's elastic stretch */
    FORTRAN(damcbevolget,(&m,&kind));
    if((int)dmcon[(ndmat_+1)*ntmat_*(m-1)+1]==1)
      uf=dmcon[(ndmat_+1)*ntmat_*(m-1)+4];
    else
      uf=dmcon[(ndmat_+1)*ntmat_*(m-1)+3];
    gf=0.;
    if(kind==2){gf=uf;uf=(sy>0.)?2.*gf/sy:0.;}
    if((uf>0.)&&(E>0.)&&(sy>0.)){
      double sb=hmed*sy/(uf*E),sbmax=hmax*sy/(uf*E);
      printf("[RESOLUTION]   u_f = %.3e%s: u_f/h(median) = %.3e; "
             "snap-back h*sigma_y/(u_f*E) median %.3e, max %.3e [<1]%s\n",
             uf,(kind==2)?" (2G_f/sigma_y)":"",uf/hmed,sb,sbmax,
             (sbmax>=1.)?" !!":"");
    }
  }
  SFREE(h);

  /* cohesive sections, one line per distinct property set */
  {
    ITG last=-1,ip,nf;
    for(;;){
      ITG next=-1;
      for(i=0;i<ne;i++){
        if(ipkon[i]<0) continue;
        if(strncmp(&lakon[8*i],"UC6",3)!=0) continue;
        if(ielprop[i]<0) continue;
        if((ielprop[i]>last)&&((next<0)||(ielprop[i]<next))) next=ielprop[i];
      }
      if(next<0) break;
      last=next;ip=next;
      {
        double kn=prop[ip],tn0=prop[ip+1],ts0=prop[ip+2],gc=prop[ip+3],
          gmin=prop[ip+4],mu=prop[ip+5],d0,df,lcz,smed,*s=NULL;
        nf=0;
        for(i=0;i<ne;i++)
          if((ipkon[i]>=0)&&(ielprop[i]==ip)&&
             (strncmp(&lakon[8*i],"UC6",3)==0)) nf++;
        NNEW(s,double,nf);
        n=0;
        for(i=0;i<ne;i++){
          ITG k,a=0,b=1,c=2,idx;
          double e1[3],e2[3],cr[3];
          if((ipkon[i]<0)||(ielprop[i]!=ip)||
             (strncmp(&lakon[8*i],"UC6",3)!=0)) continue;
          idx=ipkon[i];
          for(k=0;k<3;k++){
            e1[k]=co[3*(kon[idx+b]-1)+k]-co[3*(kon[idx+a]-1)+k];
            e2[k]=co[3*(kon[idx+c]-1)+k]-co[3*(kon[idx+a]-1)+k];
          }
          cr[0]=e1[1]*e2[2]-e1[2]*e2[1];
          cr[1]=e1[2]*e2[0]-e1[0]*e2[2];
          cr[2]=e1[0]*e2[1]-e1[1]*e2[0];
          s[n++]=sqrt(0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]));
        }
        qsort(s,n,sizeof(double),dres_cmp);
        smed=(n>0)?s[n/2]:0.;
        SFREE(s);
        if(!(kn>0.)||!(tn0>0.)) continue;
        d0=tn0/kn;df=2.*gc/tn0;
        printf("[RESOLUTION] UC6 section (prop %" ITGFORMAT "): %" ITGFORMAT
               " facets, size sqrt(A) median %.3e; delta_0 = %.3e, "
               "delta_f = %.3e (delta_f/delta_0 = %.1f)\n",
               ip+1,nf,smed,d0,df,df/d0);
        if((Ebulk>0.)&&(smed>0.)){
          lcz=Ebulk*gc/(tn0*tn0);
          printf("[RESOLUTION]   cohesive zone l_cz = E*Gc/Tn0^2 = %.3e: "
                 "l_cz/facet = %.2f [>=3, Turon 2007]%s\n",
                 lcz,lcz/smed,(lcz/smed<3.)?" !!":"");
        }
        if((Ebulk>0.)&&(hbulk>0.)){
          printf("[RESOLUTION]   penalty Kn*h/E = %.1f [>=50, Turon 2007]%s "
                 "(E, h of the largest damage material)\n",
                 kn*hbulk/Ebulk,(kn*hbulk/Ebulk<50.)?" !!":"");
        }
        printf("[RESOLUTION]   residual gmin*Kn*delta_f/Tn0 = %.3e [<<1]%s; "
               "shear/normal strength Ts0/Tn0 = %.3f\n",
               gmin*kn*df/tn0,(gmin*kn*df/tn0>0.1)?" !!":"",ts0/tn0);
        if((mu>0.)&&(dt0>0.)){
          printf("[RESOLUTION]   viscosity mu/dt0 = %.3e (initial step "
                 "%.3e)\n",mu/dt0,dt0);
        }
      }
    }
  }

  if((visceta>0.)&&(dt0>0.)){
    printf("[RESOLUTION] bulk damage viscosity eta/dt0 = %.3e (initial step "
           "%.3e)\n",visceta/dt0,dt0);
  }
  printf("\n");
  fflush(stdout);
}
