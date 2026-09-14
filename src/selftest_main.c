/*     CalculiX - damage/fracture extension                              */
/*     selftest_main.c: run every self test without a finite element deck. */

/* Why this program exists
   -----------------------
   Every self test in this tree used to be reachable only from inside
   nonlingeo(), so running one cost a mesh, a load history and a linear
   solve.  A test that costs minutes is a test that runs once a day, and
   once a day is too late to catch the mistake while the hand is still on
   it.  The evidence arrived faster than the argument: this runner found
   two real defects in its first fifteen minutes of existing - a logview
   test that could not be run outside a solver at all, and an assertion in
   a test written an hour earlier that could not fail.

   The table itself lives in selftest.c, because the solver needs it too:
   nonlingeo() runs the same nineteen through selftest_gate() as an
   interlock.  One table, two callers, no chance of a test that only one
   of them knows about.

   This program links the same archive the solver does, so what it
   exercises is what ships - not a copy, not a subset built with different
   flags.

   What this is NOT
   ----------------
   Not a replacement for the gate.  These are unit tests of pure decisions;
   whether the solver still produces the same answer on a real deck is what
   test/regress/run.py is for, and neither substitutes for the other.

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
   the solver's state dirty for whatever runs next.  It found
   logview_selftest() the first time it was used.

   Exit status is the number of failing tests, capped at 125, so it is
   usable from a hook and from test/regress/run.py's preflight.          */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

static ITG wanted(int argc,char **argv,const char *name)
{
  int i,nsel=0;
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
    for(i=0;i<NSELFTEST;i++) printf("%s\n",SELFTESTS[i].name);
    return 0;
  }
  for(i=1;i<argc;i++) if(strcmp(argv[i],"--twice")==0){ npass=2; argv[i]=""; }

  for(pass=0;pass<npass;pass++){
    if(npass>1) printf("\n########## pass %" ITGFORMAT " ##########\n",pass+1);
    for(i=0;i<NSELFTEST;i++){
      if(!wanted(argc,argv,SELFTESTS[i].name)) continue;
      printf("\n===== %s =====\n",SELFTESTS[i].name);
      fflush(stdout);
      nbad=SELFTESTS[i].fn();
      nrun++;
      if(nbad!=0){
        nfail++;
        printf("[SELFTEST] %s FAILED with %" ITGFORMAT " failure(s)%s\n",
               SELFTESTS[i].name,nbad,(pass>0)?" ON THE SECOND PASS - it did"
               " not restore what it borrowed":"");
      }
    }
  }

  printf("\n[SELFTEST] %" ITGFORMAT " test(s) run, %" ITGFORMAT " failed -- %s\n",
         nrun,nfail,nfail?"FAILED":"PASSED");
  fflush(stdout);
  if(nrun==0){
    printf("[SELFTEST] no test matched; --list shows the table\n");
    return 126;
  }
  return (int)(nfail>125?125:nfail);
}
