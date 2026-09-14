/*     CalculiX - damage/fracture extension                              */
/*     ccxopt.h: the switch registry: one place reads CCX_*. */

#ifndef CCX_CCXOPT_H
#define CCX_CCXOPT_H

#include "CalculiX.h"

/* ccxopt.c - Options: the one place that knows what this binary can be told
   to do.  Every CCX_* read goes through ccxopt_getenv, which is getenv plus
   the note that this run read it; ccxopt_decl.h declares type, default,
   range, legal spellings and prose for the options that have an owner;
   ccxopt_report states the configuration of the run, validates the declared
   ones, names anything set that is not read, and at exit names anything set
   that was never read.  Replaces damswitch.c, which did about a tenth of
   this.  PETSc's options database is the model. */
const char *ccxopt_getenv(const char *name);
void ccxopt_report(void);
ITG  ccxopt_selftest(void);

#endif /* CCX_CCXOPT_H */
