/*     CalculiX - damage/fracture extension                              */
/*     logview.h: where the wall clock went. */

#ifndef CCX_LOGVIEW_H
#define CCX_LOGVIEW_H

#include "CalculiX.h"

/* logview.c - where the run spends its time.  Named events with a call
   count, an inclusive time and a self time, in the shape of PETSc's
   -log_view.  Measurement only: off unless CCX_LOG_VIEW is set, on no
   solution path, and suppressed entirely if its own self test fails. */
ITG  logview_enabled(void);
ITG  logview_event(const char *name);
void logview_begin(ITG id);
void logview_end(ITG id);
void logview_begin_named(const char *name);
void logview_end_named(const char *name);
void logview_report(double totalseconds);
ITG  logview_selftest(void);

#endif /* CCX_LOGVIEW_H */
