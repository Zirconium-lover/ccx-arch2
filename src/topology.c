/*     CalculiX - damage/fracture extension                              */
/*     topology.c: the transaction that commits an erosion.              */

/* Why this module exists
   ----------------------
   handover/11-TOPOLOGY.md has the argument; the short form is a count.
   nonlingeo.c held 139 references to nine `damage_tent_*` locals, and with
   no object to call, two blocks got written more than once:

     - "discard the transaction" - four SFREE, four =NULL and count=0 -
       appeared FOUR times;
     - "collect what was eroded" - forty-five lines - appeared TWICE,
       differing by one blank line and one closing brace.

   Nobody chose to write the collection loop twice.  It happened because
   there was nothing to call.

   Contract
   --------
     - the marked set is one object with one lifetime, on the DMLabel
       pattern (PETSc include/petscdmlabel.h): created once, passed around,
       discarded by one call.  Freeing three of the four arrays and leaving
       the fourth is no longer expressible;
     - discard is idempotent and safe on a transaction that was never
       opened, because four copies of a lifetime is exactly the shape that
       produces a double free;
     - it records; it does not decide.  What counts as eroded is the damage
       model's answer and this object observes it.

   What this is NOT
   ----------------
   Not a rollback.  Discarding frees the marked set; it does not restore
   ipkon, and there is no code here that pretends it does.  Not the
   sparsity decision either - a committed batch is when the symbolic
   factorization stops being valid, which is a behaviour change with its
   own hypothesis and not part of an extraction.                        */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

void topo_txn_init(topo_txn *t)
{
  t->elem=NULL; t->mat=NULL; t->ip=NULL; t->value=NULL;
  t->count=0; t->step=0; t->increment=0;
  t->step_time=0.; t->total_time=0.;
}

/* Release the marked set.  This is the one that was written out four
   times.  Safe on a transaction that holds nothing, and safe twice: the
   pointers are nulled as they go, which is the property the four copies
   each had to remember separately. */
void topo_txn_discard(topo_txn *t)
{
  if(t->elem!=NULL){
    SFREE(t->elem);  t->elem=NULL;
    SFREE(t->mat);   t->mat=NULL;
    SFREE(t->ip);    t->ip=NULL;
    SFREE(t->value); t->value=NULL;
  }
  t->count=0;
}

/* How many integration points does this element type carry?

   Moved here from a static in nonlingeo.c, where sixteen call sites shared
   it and nothing could test it.  It is a property of the mesh, which is
   what this object owns, and it is a pure function of the element label -
   so it is the one part of this module that can be tested exhaustively
   rather than by example. */
ITG topo_element_nip(const char *lakonel,ITG mi0)
{
  if((lakonel[6]=='L')&&(lakonel[7]=='C')) return mi0;
  if(strncmp(lakonel,"C3D20RB",7)==0) return mi0;
  if(strncmp(lakonel,"C3D8R",5)==0) return 1;
  if(strncmp(lakonel,"C3D8I",5)==0) return 8;
  if(strncmp(lakonel,"C3D20R",6)==0) return 8;
  if(strncmp(lakonel,"C3D20",5)==0) return 27;
  if(strncmp(lakonel,"C3D10",5)==0) return 4;
  if(strncmp(lakonel,"C3D4",4)==0) return 1;
  if(strncmp(lakonel,"C3D15",5)==0) return 9;
  if(strncmp(lakonel,"C3D6",4)==0) return 2;
  if(strncmp(lakonel,"C3D8",4)==0) return 8;
  return mi0;
}

/* Collect the elements this increment eroded.

   This is the block that was written out THREE times - 13472, 14382 and
   14850 in the pre-extraction file, identical once whitespace and comments
   are stripped.  Two passes over the same predicate, `was alive at the
   baseline and is dead now`: one to size the marked set, one to fill it.

   The DE1.3 override at the end is deal.II's prepare_coarsening_and_-
   refinement: the marked set is adjusted before it is committed, because
   the DE1.3 path already knows the value and the integration point that
   triggered the deletion and they are better than what the scan found.
   This tree had that phase and no name for it.

   Discards first, so that "collect" means "the marked set is now exactly
   this".  Two of the three sites discarded immediately before; the third
   did not, and was safe only because nothing allocates between its path's
   discard and it - an argument the reader had to make from a hundred and
   sixty lines away.  Now nobody has to make it.

   damage_progressive_material is called directly and stays defined in
   nonlingeo.c: whether a material's damage is progressive is a Material
   question and that object does not exist yet.  Building a per-element
   predicate array to avoid the dependency would add an allocation and a
   loop the original never did - paying real work to hide a honest
   coupling.  The coupling is named here instead.                       */
void topo_txn_collect(topo_txn *t,ITG step,ITG increment,
                      double step_time,double total_time,
                      ITG ne0,const ITG *ipkondamageini,const ITG *ipkon,
                      const ITG *ielmat,const ITG *mi,const char *lakon,
                      const double *dam,
                      const ITG *ndmcon,const double *dmcon,
                      ITG ndmat_,ITG ntmat_,
                      ITG de13_transaction,const double *de13_trigger_value,
                      const ITG *de13_trigger_ip)
{
  ITG i,j,nip,cap=0;
  double dmax;

  topo_txn_discard(t);

  for(i=0;i<ne0;i++){
    if((ipkondamageini[i]>=0)&&(ipkon[i]<0)) cap++;
  }
  if(cap<=0) return;

  NNEW(t->elem,ITG,cap);
  NNEW(t->mat,ITG,cap);
  NNEW(t->ip,ITG,cap);
  NNEW(t->value,double,cap);

  t->step=step;
  t->increment=increment;
  t->step_time=step_time;
  t->total_time=total_time;
  t->count=0;

  for(i=0;i<ne0;i++){
    if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
      t->elem[t->count]=i+1;
      t->mat[t->count]=ielmat[mi[2]*i];

      nip=topo_element_nip(&lakon[8*i],mi[0]);
      if(nip<1) nip=1;
      if(nip>mi[0]) nip=mi[0];

      dmax=dam[mi[0]*i];
      t->ip[t->count]=1;
      for(j=1;j<nip;j++){
        if(dam[mi[0]*i+j]>dmax){
          dmax=dam[mi[0]*i+j];
          t->ip[t->count]=j+1;
        }
      }
      if((t->mat[t->count]>0)&&
         damage_progressive_material(t->mat[t->count],ndmcon,dmcon,
                                     ndmat_,ntmat_)&&
         (dmax>1.)) dmax-=1.;

      if((de13_transaction)&&(de13_trigger_value!=NULL)&&
         (de13_trigger_value[i]>=0.)){
        t->value[t->count]=de13_trigger_value[i];
        t->ip[t->count]=de13_trigger_ip[i];
      }else{
        t->value[t->count]=dmax;
      }
      t->count++;
    }
  }
}

/* ------------------------------------------------------------ the commit

   Step C.  Unlike discard and collect this appears ONCE, so the argument
   for moving it is not deduplication - it is that the record format had
   two implementations in this repository and no owner for either.

   The nine fields written here are parsed by tools/ccxdiff.py
   (read_damage), which is the comparison tool every bit-identity claim in
   this project rests on.  If the writer gains a field and the reader does
   not, ccxdiff keeps comparing and starts comparing the wrong columns -
   a silent failure in the instrument rather than in the thing measured,
   which is the expensive kind.  The order is therefore stated once, here,
   and the self test writes a batch and reads it back to prove the writer
   still produces it.

   The printing is Monitor's job and this is a way-station, the same
   admission converge_report makes: the batch line is here so that the
   commit is one call, not so that this file owns output.              */

/* element step increment step_time total_time material damage ip batch */
#define TOPO_HISTORY_FIELDS 9

/* Writing the record and announcing it are two things, and the self test
   is what made that obvious: a commit that always prints puts a
   [DAMAGE COMMIT] line into every run's log at start-up, from a test. */
void topo_txn_write_history(const topo_txn *t,FILE *fdamage,ITG batch)
{
  ITG i;

  if(t->count<=0) return;

  if(fdamage!=NULL){
    for(i=0;i<t->count;i++){
      fprintf(fdamage,
              "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT
              " %.15e %.15e %" ITGFORMAT " %.15e %" ITGFORMAT
              " %" ITGFORMAT "\n",
              t->elem[i],t->step,
              t->increment,t->step_time,
              t->total_time,t->mat[i],
              t->value[i],t->ip[i],batch);
    }
    fflush(fdamage);
  }
}

void topo_txn_commit(const topo_txn *t,FILE *fdamage,ITG batch,
                     ITG de13_transaction,ITG active_pass)
{
  if(t->count<=0) return;

  topo_txn_write_history(t,fdamage,batch);

  printf("[DAMAGE COMMIT] batch=%" ITGFORMAT
         " inc=%" ITGFORMAT " time=%.12e deleted=%" ITGFORMAT
         " active_passes=%" ITGFORMAT "\n",
         batch,t->increment,
         t->step_time,t->count,
         active_pass);
  if(de13_transaction){
    printf("[DAMAGE DE1.3 COMMIT] inc=%" ITGFORMAT
           " time=%.12e terminal_deleted=%" ITGFORMAT
           " active_passes=%" ITGFORMAT "\n",
           t->increment,t->step_time,
           t->count,active_pass);
  }
  fflush(stdout);
}

/* ---------------------------------------------------------------- tests */

static ITG topo_chki(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-36s got=%-10" ITGFORMAT " want=%-10" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG topo_selftest(void)
{
  ITG nbad=0;
  topo_txn t;

  printf("[TOPOLOGY] self test%s","\n");

  /* the element map, exhaustively: eleven branches and a fallback, and a
     pure function of the label is the one thing here that can be tested
     completely rather than by example */
  {
    struct{const char *lab;ITG mi0,want;}m[]={
      {"C3D8R   ",27,1},{"C3D8I   ",27,8},{"C3D20R  ",27,8},
      {"C3D20   ",27,27},{"C3D10   ",27,4},{"C3D4    ",27,1},
      {"C3D15   ",27,9},{"C3D6    ",27,2},{"C3D8    ",27,8},
      {"C3D20RB ",27,27},          /* the RB form defers to mi0 */
      {"C3D8  LC",27,27},          /* composite layer: also mi0.  lakon is
                                      EXACTLY eight characters and the L/C
                                      marker sits at 6 and 7, so a nine-
                                      character label is not a composite -
                                      which is what this expectation said
                                      on the first attempt, and the test
                                      caught its own author. */
      {"UNKNOWN ",13,13}};         /* fallback is mi0, not a guess */
    ITG k,allok=1;
    for(k=0;k<12;k++){
      if(topo_element_nip(m[k].lab,m[k].mi0)!=m[k].want){
        printf("   nip(%s) got=%" ITGFORMAT " want=%" ITGFORMAT "%s",
               m[k].lab,topo_element_nip(m[k].lab,m[k].mi0),m[k].want,"\n");
        allok=0;
      }
    }
    topo_chki("element nip map, all 12 forms",allok,1,&nbad);
  }

  topo_txn_init(&t);
  topo_chki("a fresh transaction holds nothing",
            (t.elem==NULL)&&(t.count==0),1,&nbad);

  /* discarding a transaction that was never opened is a no-op, not a
     crash - the case four copies of the block each had to get right */
  topo_txn_discard(&t);
  topo_chki("discard of an empty transaction is safe",
            (t.elem==NULL)&&(t.count==0),1,&nbad);

  NNEW(t.elem,ITG,4); NNEW(t.mat,ITG,4); NNEW(t.ip,ITG,4);
  NNEW(t.value,double,4);
  t.count=4; t.step=2; t.increment=17;
  t.elem[0]=11; t.value[3]=0.5;
  topo_chki("an opened transaction holds its marks",
            (t.elem!=NULL)&&(t.count==4)&&(t.elem[0]==11),1,&nbad);

  topo_txn_discard(&t);
  topo_chki("discard releases every array",
            (t.elem==NULL)&&(t.mat==NULL)&&(t.ip==NULL)&&(t.value==NULL),
            1,&nbad);
  topo_chki("discard zeroes the count",t.count,0,&nbad);

  /* the stamp survives discard: what was committed and when is not
     invalidated by releasing the marked set */
  topo_chki("discard leaves the stamp alone",
            (t.step==2)&&(t.increment==17),1,&nbad);

  /* twice, because four copies of a lifetime is the shape that produces a
     double free and this is the check that would have caught it */
  topo_txn_discard(&t);
  topo_chki("discard is idempotent",
            (t.elem==NULL)&&(t.count==0),1,&nbad);

  /* collect, on a four-element mesh where two elements died this
     increment.  mat=0 everywhere so the progressive-material branch is
     not reached and no material tables are needed - the point here is the
     scan, the worst-integration-point search and the DE1.3 override. */
  {
    ITG ipkondamageini[4]={0,0,0,-1};   /* element 3 was already dead */
    ITG ipkon[4]={-1,0,-1,-1};          /* 0 and 2 died this increment */
    ITG ielmat[4]={0,0,0,0};
    ITG mi[3]={2,0,1};                  /* two integration points       */
    char lakon[33]="C3D8    C3D8    C3D8    C3D8    ";
    double dam[8]={0.3,0.7,  0.0,0.0,  0.9,0.1,  0.0,0.0};
    double trig[4]={-1.,-1.,-1.,-1.};
    ITG trigip[4]={0,0,0,0};
    topo_txn c;

    topo_txn_init(&c);
    topo_txn_collect(&c,3,42,0.25,7.25,4,ipkondamageini,ipkon,ielmat,mi,
                     lakon,dam,NULL,NULL,0,0,0,NULL,NULL);
    topo_chki("collect finds both newly dead elements",c.count,2,&nbad);
    topo_chki("  and not the one dead at the baseline",
              (c.count==2)&&(c.elem[0]==1)&&(c.elem[1]==3),1,&nbad);
    topo_chki("  worst integration point of element 1 is 2",
              (c.count==2)&&(c.ip[0]==2),1,&nbad);
    topo_chki("  worst integration point of element 3 is 1",
              (c.count==2)&&(c.ip[1]==1),1,&nbad);
    topo_chki("  and the stamp is the one it was given",
              (c.step==3)&&(c.increment==42),1,&nbad);

    /* the DE1.3 override - deal.II's prepare phase: the marked set is
       adjusted before commit, because the DE1.3 path knows the value and
       the point that triggered the deletion and the scan does not */
    trig[2]=0.55; trigip[2]=2;
    topo_txn_collect(&c,3,42,0.25,7.25,4,ipkondamageini,ipkon,ielmat,mi,
                     lakon,dam,NULL,NULL,0,0,1,trig,trigip);
    topo_chki("DE1.3 override replaces value and point",
              (c.count==2)&&(c.ip[1]==2),1,&nbad);
    topo_chki("  and leaves the element it did not flag alone",
              (c.count==2)&&(c.ip[0]==2),1,&nbad);

    /* collect replaces the marked set rather than appending to it: the
       third call must not find four elements */
    topo_txn_collect(&c,3,43,0.25,7.25,4,ipkondamageini,ipkon,ielmat,mi,
                     lakon,dam,NULL,NULL,0,0,0,NULL,NULL);
    topo_chki("collect replaces, it does not append",c.count,2,&nbad);
    topo_txn_discard(&c);
  }

  /* the record format, round-tripped.  tools/ccxdiff.py parses these nine
     fields positionally, so the contract this test defends is the ORDER
     and the COUNT: a writer that gains a field the reader does not expect
     does not make ccxdiff fail, it makes ccxdiff compare the wrong
     columns and keep going. */
  {
    topo_txn c;
    FILE *f;
    char line[512];
    double f3,f4,f6;
    long f0,f1,f2,f5,f7,f8;
    ITG nf=0,nrow=0,fieldsok=0;
    const char *path="topo_selftest_history.tmp";

    topo_txn_init(&c);
    NNEW(c.elem,ITG,2); NNEW(c.mat,ITG,2); NNEW(c.ip,ITG,2);
    NNEW(c.value,double,2);
    c.count=2; c.step=3; c.increment=42;
    c.step_time=0.25; c.total_time=7.25;
    c.elem[0]=11; c.mat[0]=2; c.ip[0]=5; c.value[0]=0.75;
    c.elem[1]=13; c.mat[1]=2; c.ip[1]=1; c.value[1]=0.5;

    f=fopen(path,"w");
    if(f!=NULL){
      topo_txn_write_history(&c,f,9);
      fclose(f);
      f=fopen(path,"r");
      if(f!=NULL){
        while(fgets(line,sizeof(line),f)!=NULL){
          nf=sscanf(line,"%ld %ld %ld %lf %lf %ld %lf %ld %ld",
                    &f0,&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8);
          nrow++;
          if(nrow==1){
            fieldsok=((nf==TOPO_HISTORY_FIELDS)&&
                      (f0==11)&&(f1==3)&&(f2==42)&&
                      (f3==0.25)&&(f4==7.25)&&(f5==2)&&
                      (f6==0.75)&&(f7==5)&&(f8==9));
          }
        }
        fclose(f);
      }
      remove(path);
    }
    topo_chki("history: one row per marked element",nrow,2,&nbad);
    topo_chki("history: nine fields, in the order ccxdiff reads",
              fieldsok,1,&nbad);
    topo_txn_discard(&c);
  }

  printf("[TOPOLOGY] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
