/*     CalculiX - damage/fracture extension                              */
/*     nlstate.c: where the Newton solve IS.                             */

/* Why this module exists
   ----------------------
   Twenty-seven of this fork's public functions take more than eight
   arguments; the widest take twenty-two.  A function with twenty-two
   arguments does not have an interface - it has a copy of its caller's
   locals, and every local the caller adds is an argument it grows.

   Reading those twenty-seven says the lists got long for exactly three
   reasons, and they have three different fixes.  This file is the second
   of the three.  The others are scalars DERIVED from the model and passed
   by hand - mi[0], mi[1]+1, thirty-three hand-offs - and model state the
   evaluator's context never learned about.

   This one is state that NO OBJECT OWNS.  Where the Newton solve is - the
   iteration number, the residual norms this iteration and the two before
   it, the correction norms, the reference force, the control table - is
   asked for by the convergence test, the slow-Newton gate, the continuation
   corrector and the line searches.  Every one of them had to take it apart
   as loose arguments, because there was nowhere to take it FROM.

   The shape is trialctx's, and for the same reason
   -----------------------------------------------
   Every field holds the ADDRESS of the caller's storage, never a copy.  The
   iteration counter changes every iteration and the norms are rewritten in
   place, so a context holding values would be stale before its first use.
   For the fixed-size arrays of nonlingeo() - ram[6], cam[5], qa[4] - the
   array IS the address, so the rule reads the same way for both and
   nlstate_check() can walk the whole struct as an array of pointers.

   Deliberately small
   ------------------
   It holds what the functions being converted actually ask for, and not one
   field more.  A context that tries to anticipate its callers becomes the
   second copy of nonlingeo()'s prologue, which is the thing being fixed.
   Adding a field is one line here and one in NLSTATE_BIND.

   What this is NOT
   ----------------
   Not the model.  What the mesh is, what the material does and where the
   damage stands are trialctx's, and this file knows nothing about them.
   Not a decision: it holds numbers and takes no view.  Everything that
   reads it decides for itself.                                          */

#include <stdio.h>
#include "CalculiX.h"
#include "ccxfork.h"

/* Every field is the address of one of the caller's locals, so every field
   is non-NULL.  A field the bind macro forgot is therefore a NULL, and this
   loop is the only thing between that and a wrong answer three hundred
   increments later.  Walking the struct as an array of pointers is legal
   here because every member is a pointer and C guarantees no padding
   between members of identical type - the same argument trial_check()
   makes, for the same layout. */
ITG nlstate_check(const nlstate *n)
{
  const void *const *p=(const void *const *)n;
  ITG k=(ITG)(sizeof(nlstate)/sizeof(void *)),i,nbad=0;

  for(i=0;i<k;i++) if(p[i]==NULL) nbad++;
  printf("[NLSTATE] bound: %" ITGFORMAT " field(s), %" ITGFORMAT
         " unbound -- %s\n",k,nbad,nbad?"FAILED":"PASSED");
  fflush(stdout);
  return nbad;
}
