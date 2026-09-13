/*     CalculiX - damage/fracture extension                              */
/*     rescue.c: the state of what happens when an increment fails.      */

/* Why this module exists
   ----------------------
   Six clusters of locals in nonlingeo() - sixty-seven of them - are one
   mechanism seen from six angles: the damage line search, its probe,
   transactional backtracking, the same-load re-equilibration that follows
   a deletion, the bounded recovery corridor, and the levels that order
   them.  Together they are the answer to "what does this solver do when an
   increment will not converge", and that answer had no place to live.

   The fields keep their cluster's name on purpose.  Four of the six had a
   `mode'; flattening the names would have merged four distinct variables
   into one field, compiled cleanly, and changed the answer.

   Contract
   --------
     - rescue is the state: what each level is configured to do, what it
       has spent, and the snapshots it restores from.  One object, one
       lifetime, initialised once;
     - the LADDER - which level fires, in what order, and what a failure at
       the top does - is still nonlingeo's, because it cuts increments,
       rolls back topology and ends steps.  Moving the state is what makes
       moving the ladder possible later; it is not the same thing and this
       file does not claim to be it.

   The measured reason any of this exists is in globalize.c: six
   mechanisms stacked in a fixed order with no interface, three
   consecutive attempts producing bit-identical residual sequences, and two
   rescue levels that had run and changed nothing.                      */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

/* The values the sixty-seven locals carried at their declarations. */
void rescue_init(rescue *r)
{
  memset(r,0,sizeof(*r));
  r->bt_dam=NULL; r->bt_visc=NULL; r->bt_xs=NULL;
  r->linesearch_env=NULL; r->linesearch_step=NULL; r->reeq_scale_env=NULL;
  r->bt_floor=0.015625;
  r->bt_growth=1.;
  r->bt_window=1;
  r->corr_exit=5;
  r->corr_grace=10;
  r->corr_maxesc=3;
  r->corr_maxinc=60;
  r->corr_maxwall=30;
  r->corr_stableneed=3;
  r->corr_tryevery=5;
  r->corr_minfrac=0.05;
  r->ls_trials=DAMAGE_LINESEARCH_MAX_TRIALS;
  r->ls_min=DAMAGE_LINESEARCH_MIN;
  r->rescue_maxlevel=1;
}
