/*     Options: the one place that knows what this binary can be told to do.
 *
 *     The state this replaces, measured rather than asserted (docs/SWITCHES.md
 *     is generated): 140 CCX_* names read from 177 getenv calls in C and a
 *     handful more in Fortran; 62 with nothing written about them anywhere
 *     but the line that reads them; 123 that no test in this tree sets.
 *     Nothing declared a type, a default, a range or a legal spelling, so
 *     "what may this be set to" had no answer short of reading the parse.
 *
 *     PETSc's options database is the model - PetscOptionsBegin/End, typed
 *     getters, -help producing the complete documented list, and
 *     -options_left warning about anything that was SET and never QUERIED.
 *     src/damswitch.c, which this replaces, had reinvented about a tenth of
 *     that: it could say which names are set and which are unknown, and
 *     nothing about what any of them means or whether the run read it.
 *
 *     THE MIGRATION, and why it is shaped this way.  Every call site keeps
 *     its own parsing.  This is deliberate and it is not laziness: the tree
 *     has at least four incompatible notions of "true" -
 *
 *         CCX_DAMAGE_TR_DOGLEG      set to ANYTHING, including "0"
 *         CCX_FRACTURE_DEADFACET    anything except the string "0"
 *         CCX_DAMAGE_LINESEARCH     "ADAPTIVE", "adaptive" or "1"
 *         CCX_PARDISO_REUSE_SYMBOLIC  "1", "ON", "on", "YES" or "yes"
 *
 *     - and a single ccxopt_bool() would silently change the behaviour of
 *     whichever sites it did not match.  So the first pass routes every read
 *     through ccxopt_getenv(), which returns exactly what getenv() returned
 *     and is therefore bit-identical by construction, while the registry
 *     gains the two things the tree never had: a DECLARATION per option
 *     (type, default, range, legal spellings, one line of prose) and a
 *     record of which options this run actually READ.
 *
 *     The declarations are in ccxopt_decl.h and cover the options that any
 *     test or the target-deck script sets.  Everything else is reported as
 *     UNDECLARED, and that list is the retirement queue, produced
 *     mechanically instead of from memory.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"
#include "ccxopt_list.h"
#include "ccxopt_decl.h"

extern char **environ;

/*     Names that start with CCX_ but are not options, so that the unknown
 *     report stays worth reading.  CCX_EXE is set by the run scripts to
 *     choose the binary; CCX_JOBNAME_GETJOBNAME is set by ccx_2.23.c itself
 *     with putenv so that getjobname can find the job name.  Neither is ever
 *     read as an option, and warning about them every run would train the
 *     reader to ignore the warning, which is the only way this report can
 *     fail. */
static const char *const ccxopt_notanoption[]={
  "CCX_EXE",
  "CCX_JOBNAME_GETJOBNAME",
  NULL
};

/* which known names this run actually queried, indexed like ccxopt_known_name */
static char ccxopt_queried[CCXOPT_KNOWN_COUNT];
static ITG  ccxopt_reported=0;
/* The self test proves validation by feeding it values it must REJECT.  The
   rejections are the passing case, so they must not print: a run log that
   says *ERROR when nothing is wrong trains the reader to ignore *ERROR, and
   the gate reads that word as a failure. */
static ITG  ccxopt_quiet=0;

static ITG ccxopt_ignored(const char *name,size_t n){
  ITG i;
  for(i=0;ccxopt_notanoption[i]!=NULL;i++){
    if((strlen(ccxopt_notanoption[i])==n)&&
       (strncmp(ccxopt_notanoption[i],name,n)==0)) return 1;
  }
  return 0;
}

/* the generated list is sorted, so this is a binary search; it runs on every
   option read and there are tens of thousands of those in a long run */
static ITG ccxopt_index(const char *name,size_t n){
  ITG lo=0,hi=CCXOPT_KNOWN_COUNT-1;
  while(lo<=hi){
    ITG mid=(lo+hi)/2,c;
    size_t m=strlen(ccxopt_known_name[mid]);
    c=strncmp(ccxopt_known_name[mid],name,n);
    if(c==0) c=(m<n)?-1:((m>n)?1:0);
    if(c==0) return mid;
    if(c<0) lo=mid+1; else hi=mid-1;
  }
  return -1;
}

static ITG ccxopt_known(const char *name,size_t n){
  return ccxopt_index(name,n)>=0;
}

static const ccxopt_decl *ccxopt_declared(const char *name){
  ITG i;
  for(i=0;i<CCXOPT_DECL_COUNT;i++){
    if(strcmp(ccxopt_decl_table[i].name,name)==0) return &ccxopt_decl_table[i];
  }
  return NULL;
}

/* ------------------------------------------------------------------ read */

/*  Exactly getenv(), plus the note that this run read this option.
 *
 *  Returning getenv()'s own pointer rather than a copy is what makes the
 *  migration provably free: every call site that used to hold that pointer
 *  still holds the same one, with the same lifetime.  */
const char *ccxopt_getenv(const char *name){
  ITG i=ccxopt_index(name,strlen(name));
  if(i>=0) ccxopt_queried[i]=1;
  return getenv(name);
}

const char *ccxopt_type_name(ITG t){
  switch(t){
  case CCXOPT_BOOL:   return "flag";
  case CCXOPT_INT:    return "integer";
  case CCXOPT_REAL:   return "real";
  case CCXOPT_ENUM:   return "one of";
  case CCXOPT_STRING: return "text";
  default:            return "?";
  }
}

/* ------------------------------------------------------------ validation */

/* Is this value admissible for this declaration?  Returns 0 if it is, and
   otherwise prints why.  Checked once, at the report, which happens before
   the first solve - so an inadmissible option is named before it can matter,
   and the gate's self-test reader turns red on the *ERROR. */
static ITG ccxopt_validate(const ccxopt_decl *d,const char *value){
  char *end=NULL;
  double x;
  if(d->deprecated!=NULL){
    if(!ccxopt_quiet) printf("[SWITCHES] *ERROR: %s is deprecated: %s\n",d->name,d->deprecated);
    return 1;
  }
  if(d->type==CCXOPT_ENUM){
    const char *p=d->choices;
    size_t n=strlen(value);
    while(p!=NULL&&*p!='\0'){
      const char *bar=strchr(p,'|');
      size_t len=(bar!=NULL)?(size_t)(bar-p):strlen(p);
      if((len==n)&&(strncmp(p,value,n)==0)) return 0;
      p=(bar!=NULL)?bar+1:NULL;
    }
    if(!ccxopt_quiet) printf("[SWITCHES] *ERROR: %s=%s is not one of %s\n",
           d->name,value,d->choices);
    return 1;
  }
  if((d->type==CCXOPT_INT)||(d->type==CCXOPT_REAL)){
    x=strtod(value,&end);
    if((end==value)||((end!=NULL)&&(*end!='\0'))){
      if(!ccxopt_quiet) printf("[SWITCHES] *ERROR: %s=%s is not a number\n",d->name,value);
      return 1;
    }
    if((d->lo<=d->hi)&&((x<d->lo)||(x>d->hi))){
      if(!ccxopt_quiet) printf("[SWITCHES] *ERROR: %s=%s is outside [%g,%g]\n",
             d->name,value,d->lo,d->hi);
      return 1;
    }
  }
  return 0;
}

/* ------------------------------------------------------------- self test */

ITG ccxopt_selftest(void){
  ITG i,bad=0;
  ccxopt_decl d;
  if(CCXOPT_KNOWN_COUNT<=0){
    printf("[SWITCHES] *ERROR: the generated list is empty\n");return 1;}
  for(i=0;i<CCXOPT_KNOWN_COUNT;i++){
    if((ccxopt_known_name[i]==NULL)||(ccxopt_known_name[i][0]=='\0')){
      printf("[SWITCHES] *ERROR: entry %" ITGFORMAT " is empty\n",i);bad++;continue;}
    if(strncmp(ccxopt_known_name[i],"CCX_",4)!=0){
      printf("[SWITCHES] *ERROR: %s does not start with CCX_\n",
             ccxopt_known_name[i]);bad++;}
    if((i>0)&&(strcmp(ccxopt_known_name[i-1],ccxopt_known_name[i])>=0)){
      printf("[SWITCHES] *ERROR: %s and %s are out of order or duplicated\n",
             ccxopt_known_name[i-1],ccxopt_known_name[i]);bad++;}
  }
  /* the lookup must recognise a name that is in the list and reject one that
     is not, which is the only behaviour anything downstream depends on.  It
     is a binary search over a generated table, so it is worth testing both
     ends and not only the middle. */
  if(!ccxopt_known(ccxopt_known_name[0],strlen(ccxopt_known_name[0]))){
    printf("[SWITCHES] *ERROR: the first listed name is not recognised\n");bad++;}
  if(!ccxopt_known(ccxopt_known_name[CCXOPT_KNOWN_COUNT-1],
                   strlen(ccxopt_known_name[CCXOPT_KNOWN_COUNT-1]))){
    printf("[SWITCHES] *ERROR: the last listed name is not recognised\n");bad++;}
  if(!ccxopt_known(ccxopt_known_name[CCXOPT_KNOWN_COUNT/2],
                   strlen(ccxopt_known_name[CCXOPT_KNOWN_COUNT/2]))){
    printf("[SWITCHES] *ERROR: a middle listed name is not recognised\n");bad++;}
  if(ccxopt_known("CCX_THIS_IS_NOT_A_SWITCH",24)){
    printf("[SWITCHES] *ERROR: an unlisted name is recognised\n");bad++;}
  /* a PREFIX of a real name must not match it: the length comparison in
     ccxopt_index is the only thing preventing that, and getting it wrong
     would make the unknown-name report silently useless */
  if(ccxopt_known(ccxopt_known_name[0],4)){
    printf("[SWITCHES] *ERROR: a prefix matched a full name\n");bad++;}
  if(!ccxopt_ignored("CCX_EXE",7)){
    printf("[SWITCHES] *ERROR: the ignore list does not cover CCX_EXE\n");bad++;}
  if(ccxopt_ignored("CCX_THIS_IS_NOT_A_SWITCH",24)){
    printf("[SWITCHES] *ERROR: the ignore list swallows an unknown name\n");
    bad++;}
  for(i=0;ccxopt_notanoption[i]!=NULL;i++){
    if(ccxopt_known(ccxopt_notanoption[i],strlen(ccxopt_notanoption[i]))){
      printf("[SWITCHES] *ERROR: %s is both read and ignored\n",
             ccxopt_notanoption[i]);bad++;}
  }
  /* every declaration must name an option this binary actually reads, or
     the generated documentation would describe something that cannot be
     set; and the ranges must be the right way round */
  for(i=0;i<CCXOPT_DECL_COUNT;i++){
    const ccxopt_decl *e=&ccxopt_decl_table[i];
    if(!ccxopt_known(e->name,strlen(e->name))){
      printf("[SWITCHES] *ERROR: %s is declared but this binary never reads it\n",
             e->name);bad++;}
    if((e->doc==NULL)||(e->doc[0]=='\0')){
      printf("[SWITCHES] *ERROR: %s is declared with no documentation\n",
             e->name);bad++;}
    if((e->type==CCXOPT_ENUM)&&((e->choices==NULL)||(e->choices[0]=='\0'))){
      printf("[SWITCHES] *ERROR: %s is an enum with no choices\n",e->name);bad++;}
  }
  /* validation itself has to be shown working, in both directions, or a
     green report means only that nothing was checked */
  ccxopt_quiet=1;
  d.name="CCX_SELFTEST";d.type=CCXOPT_REAL;d.dflt="unset";
  d.lo=0.;d.hi=1.;d.choices=NULL;d.doc="self test";d.deprecated=NULL;
  if(ccxopt_validate(&d,"0.5")!=0){
    printf("[SWITCHES] *ERROR: a value inside the range was rejected\n");bad++;}
  if(ccxopt_validate(&d,"2.0")==0){
    printf("[SWITCHES] *ERROR: a value outside the range was accepted\n");bad++;}
  if(ccxopt_validate(&d,"banana")==0){
    printf("[SWITCHES] *ERROR: a value that is not a number was accepted\n");
    bad++;}
  d.type=CCXOPT_ENUM;d.choices="NODE|FACE";
  if(ccxopt_validate(&d,"FACE")!=0){
    printf("[SWITCHES] *ERROR: a listed choice was rejected\n");bad++;}
  if(ccxopt_validate(&d,"FAC")==0){
    printf("[SWITCHES] *ERROR: a prefix of a choice was accepted\n");bad++;}
  if(ccxopt_validate(&d,"EDGE")==0){
    ccxopt_quiet=0;
    printf("[SWITCHES] *ERROR: an unlisted choice was accepted\n");bad++;}
  ccxopt_quiet=0;
  return bad;
}

/* ---------------------------------------------------------------- report */

static void ccxopt_left(void);

void ccxopt_report(void){
  char **e;
  ITG narmed=0,nunknown=0,nbadvalue=0;
  /* nonlingeo() is entered once per *STEP, and the environment cannot change
     between them, so a second report would be noise that trains the reader
     to skip the block - which is the one way this can fail. */
  if(ccxopt_reported) return;
  ccxopt_reported=1;
  if(ccxopt_selftest()!=0){
    printf("[SWITCHES] the registry self test failed; reporting nothing "
           "rather than reporting something that may be wrong\n");
    fflush(stdout);
    return;
  }
  printf("[SWITCHES] this binary reads %d CCX_* names, %d of them declared "
         "with a type, a default, a range and a line of documentation.  "
         "In force now:\n",CCXOPT_KNOWN_COUNT,CCXOPT_DECL_COUNT);
  for(e=environ;*e!=NULL;e++){
    char *eq=strchr(*e,'=');
    size_t n;
    if(eq==NULL) continue;
    n=(size_t)(eq-*e);
    if((n<4)||(strncmp(*e,"CCX_",4)!=0)) continue;
    if(ccxopt_known(*e,n)){
      char name[128];
      const ccxopt_decl *d;
      if(n>=sizeof(name)) n=sizeof(name)-1;
      memcpy(name,*e,n); name[n]='\0';
      d=ccxopt_declared(name);
      if(d!=NULL){
        printf("[SWITCHES]   %s = %s   [%s%s%s, default %s]\n",
               name,eq+1,ccxopt_type_name(d->type),
               (d->type==CCXOPT_ENUM)?" ":"",
               (d->type==CCXOPT_ENUM)?d->choices:"",d->dflt);
        nbadvalue+=ccxopt_validate(d,eq+1);
      }else{
        printf("[SWITCHES]   %s = %s   [UNDECLARED: no type, no range, "
               "no documented default]\n",name,eq+1);
      }
      narmed++;
    }else if(!ccxopt_ignored(*e,n)){
      nunknown++;
    }
  }
  if(narmed==0) printf("[SWITCHES]   (none - stock behaviour)\n");
  if(nunknown>0){
    printf("[SWITCHES] *WARNING: %" ITGFORMAT " CCX_* name(s) are set that "
           "this binary does NOT read.  A run configured with one of these "
           "is not testing what it looks like it is testing:\n",nunknown);
    for(e=environ;*e!=NULL;e++){
      char *eq=strchr(*e,'=');
      size_t n;
      if(eq==NULL) continue;
      n=(size_t)(eq-*e);
      if((n<4)||(strncmp(*e,"CCX_",4)!=0)) continue;
      if((!ccxopt_known(*e,n))&&(!ccxopt_ignored(*e,n)))
        printf("[SWITCHES]   UNKNOWN %.*s\n",(int)n,*e);
    }
  }
  fflush(stdout);
  atexit(ccxopt_left);
}

/*  PETSc's -options_left, and the failure it prevents is the one this
 *  project has actually paid for more than once: an A/B that is not wrong
 *  but UNINFORMATIVE, because the option under test was set and never read.
 *  It has to run at the END of the run - an option read only on a branch the
 *  run reaches at increment 400 has not been "left" at increment 1.
 *
 *  Options read from Fortran are excluded, and said to be excluded: the
 *  Fortran sites do not route through here yet, so claiming they were never
 *  read would be exactly the confident nonsense this report exists to stop.
 */
static void ccxopt_left(void){
  char **e;
  ITG nleft=0,nfortran=0;
  for(e=environ;*e!=NULL;e++){
    char *eq=strchr(*e,'=');
    size_t n; ITG i;
    if(eq==NULL) continue;
    n=(size_t)(eq-*e);
    if((n<4)||(strncmp(*e,"CCX_",4)!=0)) continue;
    i=ccxopt_index(*e,n);
    if(i<0) continue;
    if(ccxopt_queried[i]) continue;
    if(ccxopt_known_fortran[i]) nfortran++; else nleft++;
  }
  if((nleft==0)&&(nfortran==0)) return;
  if(nleft>0){
    printf("[SWITCHES LEFT] *WARNING: %" ITGFORMAT " option(s) were set and "
           "never read by this run.  Whatever they were meant to change, "
           "they did not:\n",nleft);
    for(e=environ;*e!=NULL;e++){
      char *eq=strchr(*e,'=');
      size_t n; ITG i;
      if(eq==NULL) continue;
      n=(size_t)(eq-*e);
      if((n<4)||(strncmp(*e,"CCX_",4)!=0)) continue;
      i=ccxopt_index(*e,n);
      if((i<0)||(ccxopt_queried[i])||(ccxopt_known_fortran[i])) continue;
      printf("[SWITCHES LEFT]   %.*s\n",(int)n,*e);
    }
  }
  if(nfortran>0){
    printf("[SWITCHES LEFT] %" ITGFORMAT " further option(s) are read from "
           "Fortran, where reads are not yet tracked, so nothing is claimed "
           "about them.\n",nfortran);
  }
  fflush(stdout);
}
