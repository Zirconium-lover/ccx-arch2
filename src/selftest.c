/*     CalculiX - damage/fracture extension                              */
/*     selftest.c: the table of self tests, and one place that runs them. */

/* Why this module exists
   ----------------------
   Ten places inside nonlingeo() ran a self test.  Each was an INTERLOCK
   rather than a test invocation - it stopped the run, or disarmed the
   mechanism it guarded, when the test came back bad - and that is a
   property worth keeping: a solver that carries on with a mechanism known
   to be wrong produces numbers nobody can use.

   What was not worth keeping was the cost.  Four of the ten ran on every
   job; the other six ran whenever their switch armed; and between them a
   production run printed 56 to 82 lines of self test output before it
   reached increment 1.  The interlock does not need any of that output.  It
   needs a verdict.

   So the ten sites become one call, made once, before anything is armed.
   The table below is the thing: adding a test is a line in it.

   Silent unless it matters
   ------------------------
   The run is made with stdout redirected to /dev/null and restored
   afterwards, so a passing build prints NOTHING.  If anything failed, the
   failing tests are then run AGAIN with output, so the operator sees
   exactly which assertion went and why.  That is the right way round: the
   common case is free, and the rare case is fully explained rather than
   summarised.

   What changed, and what did not
   ------------------------------
   Every test that ran before still runs, and the four that stopped the run
   still stop it.  Two things are deliberately different:

     - the six tests that used to run only when their own switch armed now
       run always.  A test that fails is a broken build, and finding that
       out depends on which switches a particular job happened to set was
       never a property anybody wanted.  Verified they do not depend on
       arming: all nineteen already pass from ccx_selftest, which sets no
       switches at all;

     - a failure in one of those six used to DISARM the mechanism and carry
       on.  It now stops the run, like the other four.  Degrading quietly is
       the right answer when a MEASUREMENT is unavailable; it is the wrong
       answer when the code is wrong, because every number the rest of the
       run produces was computed by the same build.

   trial_check() is not here and does not belong here.  It verifies that a
   context has been bound to a live caller's locals - a statement about a
   running solver, not about a pure decision - and it stays at its site in
   nonlingeo().                                                           */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "CalculiX.h"

/* damrank1test is Fortran and returns its count through an argument, so it
   gets a shim rather than a special case in the loop. */
static ITG damrank1_selftest(void)
{
  ITG nbad=0;
  FORTRAN(damrank1test,(&nbad));
  return nbad;
}

const selftest_entry SELFTESTS[]={
  {"ccxopt",            ccxopt_selftest},
  {"logview",           logview_selftest},
  {"stiffcensus",       stiffcensus_selftest},
  {"dammat",            dammat_selftest},
  {"topology",          topo_selftest},
  {"damstate",          damstate_selftest},
  {"topodiag",          topodiag_selftest},
  {"opcheck",           opcheck_selftest},
  {"loadcut",           loadcut_selftest},
  {"globalize",         glob_selftest},
  {"converge",          converge_selftest},
  {"crackcontrol",      crackcontrol_selftest},
  {"damcont",           damcont_selftest},
  {"erosion",           erosion_selftest},
  {"lsladder",          lsladder_selftest},
  {"dogleg",            dogleg_selftest},
  {"pathfollow",        pathfollow_selftest},
  {"pathfollow_legacy", pathfollow_legacycheck},
  {"damrank1",          damrank1_selftest},
};
const ITG NSELFTEST=(ITG)(sizeof(SELFTESTS)/sizeof(SELFTESTS[0]));

/* Redirect stdout at the FILE DESCRIPTOR, not by reopening the stream: the
   tests print through printf and the Fortran one writes through its own
   unit, and only the fd is common to both.  Flush before and after or the
   buffered output crosses the swap and lands in the wrong place. */
static int quiet_begin(void)
{
  int saved,devnull;
  fflush(stdout);
  saved=dup(fileno(stdout));
  if(saved<0) return -1;
  devnull=open("/dev/null",O_WRONLY);
  if(devnull<0){ close(saved); return -1; }
  dup2(devnull,fileno(stdout));
  close(devnull);
  return saved;
}

static void quiet_end(int saved)
{
  if(saved<0) return;
  fflush(stdout);
  dup2(saved,fileno(stdout));
  close(saved);
}

/* Run the whole table.  verbose!=0 prints everything; verbose==0 prints
   nothing unless a test fails, and then prints that test in full.
   Returns the number of FAILING TESTS. */
ITG selftest_run_all(ITG verbose)
{
  ITG i,nfail=0,bad[64],nbad=0;
  int saved=-1;

  if(!verbose){
    saved=quiet_begin();
    /* If the redirect could not be set up, run loudly rather than not at
       all: noise is a nuisance, a skipped interlock is a wrong answer. */
    if(saved<0) verbose=1;
  }

  for(i=0;i<NSELFTEST;i++){
    if(verbose){ printf("\n===== %s =====\n",SELFTESTS[i].name); fflush(stdout); }
    if(SELFTESTS[i].fn()!=0){
      nfail++;
      if(nbad<(ITG)(sizeof(bad)/sizeof(bad[0]))) bad[nbad++]=i;
    }
  }

  if(!verbose){
    quiet_end(saved);
    /* Something failed while nobody was listening.  Say which, and run it
       again with its output so the reason is on the record. */
    for(i=0;i<nbad;i++){
      printf("[SELFTEST] %s FAILED; re-running it with output:\n",
             SELFTESTS[bad[i]].name);
      fflush(stdout);
      SELFTESTS[bad[i]].fn();
    }
  }
  return nfail;
}

/* The interlock nonlingeo() calls.  Silent on a sound build; on a bad one
   it says what broke and stops, because every number the rest of the run
   would produce was computed by the same build. */
void selftest_gate(void)
{
  if(selftest_run_all(0)==0) return;
  printf("\n*ERROR: the solver's own self tests failed.  The decisions they\n"
         "        check - which elements leave the assembly, whether an\n"
         "        increment converged, how long a step to take - are the\n"
         "        ones every number in this run would be made from, so the\n"
         "        run stops here rather than reporting a result computed\n"
         "        by a rule that is not the one that was tested.\n"
         "        Run ccx_selftest for the whole table.\n");
  fflush(stdout);
  FORTRAN(stop,());
}
