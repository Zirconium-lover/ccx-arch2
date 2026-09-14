/*     CalculiX - damage/fracture extension                              */
/*     selftest_main.c: run every self test without a finite element deck. */

/* Why this program exists
   -----------------------
   Every self test in this tree was compiled into the solver and reached only
   from inside nonlingeo().  The consequence was not that the tests were bad;
   it was that running one cost a finite element analysis.  A forty-line
   material classifier could not be proved without a mesh, a load history and
   a linear solver, and the price was paid twice over: in the minutes it took
   to ask a question, and in the fifty-six to eighty-two lines of self test
   output that every production job printed before it reached increment 1.

   A test that costs minutes is a test that gets run once a day, and once a
   day is not often enough to catch the mistake while the hand is still on
   it.  The evidence is not hypothetical.  dammat_selftest() was written the
   hour this file was planned, and being cheap to run it immediately found
   two defects IN ITSELF - an offset that made three of its assertions
   unfailable, and a missing case that let a `nconst==4' rule relax to
   `nconst>=4' unnoticed.  Neither would have been seen at deck speed,
   because neither would have been looked for.

   What this is
   ------------
   A main() and a table.  It links against the same archive the solver does,
   so the code under test is the code that ships - not a copy, not a subset
   compiled with different flags.  Adding a test is one line in the table.

   What this is NOT
   ----------------
   Not a replacement for the gate.  These are unit tests of pure decisions;
   whether the solver still produces the same answer on a real deck is what
   test/regress/run.py is for, and neither can substitute for the other.

   Not a place for tests that need a live solver state.  trial_check() stays
   in nonlingeo(), where it belongs: it verifies that a context has been
   bound to the caller's locals, which is a statement about a running
   solver and is meaningless here.

   Usage
   -----
     ccx_selftest            run everything
     ccx_selftest NAME...    run only the tests whose name contains NAME
     ccx_selftest --list     print the table and exit
     ccx_selftest --twice    run the table twice in one process

   --twice is not a convenience.  These tests measure with the real
   machinery and then put it back, so each one is a claim that it restored
   what it borrowed; running the table twice is the only thing that checks
   the claim, and a test that passes alone and fails after itself has left
   the solver's state dirty.  It found logview_selftest() the first time it
   was used: that test armed the instrument by assignment rather than
   through logview_enabled(), left the interim-report clock at zero, and its
   own logview_end() then re-entered it through logview_table().

   Exit status is the number of failing TESTS, capped at 125, so it is
   usable from a hook and from test/regress/run.py's preflight.          */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

typedef ITG (*selftest_fn)(void);

/* damrank1test is Fortran and returns its count through an argument, so it
   gets a C shim rather than a special case in the loop below. */
static ITG damrank1_selftest(void)
{
  ITG nbad=0;
  FORTRAN(damrank1test,(&nbad));
  return nbad;
}

static const struct { const char *name; selftest_fn fn; } TESTS[]={
  {"ccxopt",        ccxopt_selftest},
  {"logview",       logview_selftest},
  {"stiffcensus",   stiffcensus_selftest},
  {"dammat",        dammat_selftest},
  {"topology",      topo_selftest},
  {"damstate",      damstate_selftest},
  {"topodiag",      topodiag_selftest},
  {"opcheck",       opcheck_selftest},
  {"loadcut",       loadcut_selftest},
  {"globalize",     glob_selftest},
  {"converge",      converge_selftest},
  {"crackcontrol",  crackcontrol_selftest},
  {"damcont",       damcont_selftest},
  {"erosion",       erosion_selftest},
  {"lsladder",      lsladder_selftest},
  {"dogleg",        dogleg_selftest},
  {"pathfollow",    pathfollow_selftest},
  {"pathfollow_legacy", pathfollow_legacycheck},
  {"damrank1",      damrank1_selftest},
};
#define NTEST ((ITG)(sizeof(TESTS)/sizeof(TESTS[0])))

static ITG wanted(int argc,char **argv,const char *name)
{
  int i;
  int nsel=0;
  for(i=1;i<argc;i++) if(argv[i][0]!='\0') nsel++;
  if(nsel==0) return 1;          /* no filter, or only flags: run everything */
  for(i=1;i<argc;i++)
    if((argv[i][0]!='\0')&&(strstr(name,argv[i])!=NULL)) return 1;
  return 0;
}

int main(int argc,char **argv)
{
  ITG i,nbad,nrun=0,nfail=0,pass,npass=1;

  if((argc>1)&&(strcmp(argv[1],"--list")==0)){
    for(i=0;i<NTEST;i++) printf("%s\n",TESTS[i].name);
    return 0;
  }
  for(i=1;i<argc;i++) if(strcmp(argv[i],"--twice")==0){ npass=2; argv[i]=""; }

  for(pass=0;pass<npass;pass++){
    if(npass>1) printf("\n########## pass %" ITGFORMAT " ##########\n",pass+1);
    for(i=0;i<NTEST;i++){
      if(!wanted(argc,argv,TESTS[i].name)) continue;
      printf("\n===== %s =====\n",TESTS[i].name);
      fflush(stdout);
      nbad=TESTS[i].fn();
      nrun++;
      if(nbad!=0){
        nfail++;
        printf("[SELFTEST] %s FAILED with %" ITGFORMAT " failure(s)%s\n",
               TESTS[i].name,nbad,(pass>0)?" ON THE SECOND PASS - it did not"
               " restore what it borrowed":"");
      }
    }
  }

  printf("\n[SELFTEST] %" ITGFORMAT " test(s) run, %" ITGFORMAT " failed -- %s\n",
         nrun,nfail,nfail?"FAILED":"PASSED");
  fflush(stdout);
  /* A run that selected nothing is a mistyped filter, not a pass. */
  if(nrun==0){
    printf("[SELFTEST] no test matched; --list shows the table\n");
    return 126;
  }
  return (int)(nfail>125?125:nfail);
}
