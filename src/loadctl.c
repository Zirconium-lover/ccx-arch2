/*     CalculiX - damage/fracture extension                              */
/*     loadctl.c: who drives the load parameter.                         */

/* Why this module exists
   ----------------------
   Four mechanisms answer the same question - what sets lambda this
   increment - and all four refuse to run beside each other: dissipation
   control, path control, the arc-length boundary and the regularisation
   ladder.  They were fifty-three locals of nonlingeo() with no owner.

   The reason they are one object is not that they are one algorithm; they
   are not.  It is that the mutual exclusion is ONE RULE, and it was
   written out as four loose integers in every other mechanism's arming
   block.  The dogleg refuses to arm when the corridor or the
   regularisation is on; the continuation refuses when any of the four is.
   Those conditions now read one object and ask it one question.

   Contract
   --------
     - loadctl is the state of all four, with each field keeping its
       mechanism's prefix so that nothing is merged by accident;
     - loadctl_driving() is the question the rest of the solver asks: is
       anything here holding the load parameter?  One place to change when
       a fifth mechanism is added, instead of every arming block;
     - the mechanisms' own loops are still in nonlingeo(), because they set
       boundary conditions and end increments.  This is their state and
       their exclusion rule.

   What this is NOT
   ----------------
   Not a load-parameter abstraction.  Nothing here decides what lambda is;
   it records which mechanism has the right to.                        */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "CalculiX.h"

#include "ccxfork.h"
void loadctl_init(loadctl *c)
{
  memset(c,0,sizeof(*c));
  c->diss_env=NULL; c->diss_fhat=NULL; c->diss_uf=NULL;
  c->diss_engage_t=-1.;
  c->diss_scale=1.;
  c->diss_step=1;
  c->path_devmax=-1.;
  c->path_drop=0.25;
  c->path_nstep=20;
  c->reg_lam[0]=1.e-2; c->reg_lam[1]=1.e-1; c->reg_lam[2]=1.e0;
  c->reg_lam[3]=4.e0;  c->reg_lam[4]=1.6e1;
}

/* The one question the other mechanisms ask.  Written here so that adding
   a fifth driver is one edit rather than one per arming block. */
ITG loadctl_driving(const loadctl *c)
{
  if(c->arc==1) return 1;
  if(c->diss_ctrl>=1) return 1;
  if(c->path_on>0) return 1;
  if(c->reg_nlam>0) return 1;
  return 0;
}
