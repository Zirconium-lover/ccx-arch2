/*     Monitor: presenting a quantity, which is not the same job as
 *     producing one.
 *
 *     The diagnostics in this tree are printf calls interleaved with the
 *     solve, in ad-hoc formats, gated by environment variables.
 *     02-DIAGNOSTICS.md is eleven readings a human performs by eye, and it
 *     cannot be anything else while the output is prose.
 *
 *     PETSc's separation is the model: a monitor PRODUCES nothing and a
 *     viewer DECIDES nothing.  Here that means census.c computes and this
 *     file prints, and it prints twice - the line the tree already had, byte
 *     for byte, and a machine-readable record beside it.
 *
 *     The legacy line is preserved exactly on purpose.  test/regress/run.py
 *     reads below_1e-3 and worst=node_N_at_R out of it to pin two of the
 *     nine gate cases; changing that format would have made this extraction
 *     a behaviour change wearing a refactor's clothes.
 */
#include <stdio.h>
#include "CalculiX.h"

/*  One JSON object per line, prefixed so it can be grepped out of a log that
 *  also carries the solver's own chatter:
 *
 *      grep '^\[MON\] ' run.log | jq -s '...'
 *
 *  This is what turns 02-DIAGNOSTICS.md from a document into a program.
 */
void monitor_stiffness(const stiffcensus *c,ITG iinc,double time){
  /* the line as it has always been - the gate parses it */
  printf("[DAMAGE STIFFNESS] inc=%" ITGFORMAT " time=%.12e "
         "below_1e-3=%" ITGFORMAT " below_1e-2=%" ITGFORMAT
         " below_1e-1=%" ITGFORMAT " nonpositive=%" ITGFORMAT
         " worst=node_%" ITGFORMAT "_at_%.4e\n",
         iinc,time,c->below1,c->below2,c->below3,c->nonpositive,
         c->worst,c->worstratio);
  /* and the same numbers in a form nothing has to be taught to read */
  printf("[MON] {\"kind\":\"stiffness\",\"inc\":%" ITGFORMAT
         ",\"time\":%.12e,\"below_1e-3\":%" ITGFORMAT
         ",\"below_1e-2\":%" ITGFORMAT ",\"below_1e-1\":%" ITGFORMAT
         ",\"nonpositive\":%" ITGFORMAT ",\"nodes\":%" ITGFORMAT
         ",\"worst_node\":%" ITGFORMAT ",\"worst_ratio\":%.6e}\n",
         iinc,time,c->below1,c->below2,c->below3,c->nonpositive,
         c->nnode,c->worst,c->worstratio);
  fflush(stdout);
}

/*  The operator check, per column and in total.
 *
 *  The old [STRUCT-FD] line reported one relative error, which - now that
 *  the one-sided differences are separated - is known to have been an
 *  ambiguous number.  It is kept, so a reader who has seen the old output
 *  can still find their bearings, and the discriminating columns are added
 *  beside it: how far the two one-sided differences are from EACH OTHER,
 *  how far the better of them is from the tangent, and the classification.
 */
void monitor_opcheck(const opcheck *o,ITG iinc,ITG iit,ITG node,ITG dir,
                     double h){
  printf("[OPCHECK] node %-6" ITGFORMAT " dir %" ITGFORMAT
         "  scale=%-11.4e |ctr-asm|=%-10.3e |fwd-bwd|=%-10.3e "
         "min|side-asm|=%-10.3e  ok=%-6" ITGFORMAT " kink=%-6" ITGFORMAT
         " wrong=%-6" ITGFORMAT " both=%-6" ITGFORMAT "\n",
         node,dir,o->scale,o->relctr,o->relside,o->relbest,
         o->nok,o->nkink,o->nwrong,o->nboth);
  printf("[MON] {\"kind\":\"opcheck_column\",\"inc\":%" ITGFORMAT
         ",\"iter\":%" ITGFORMAT ",\"node\":%" ITGFORMAT ",\"dir\":%"
         ITGFORMAT ",\"h\":%.6e,\"scale\":%.6e,\"rel_ctr\":%.6e,"
         "\"rel_side\":%.6e,\"rel_best\":%.6e,\"n\":%" ITGFORMAT
         ",\"ok\":%" ITGFORMAT ",\"kink\":%" ITGFORMAT ",\"wrong\":%"
         ITGFORMAT ",\"both\":%" ITGFORMAT ",\"worst_row\":%" ITGFORMAT
         ",\"worst_dir\":%" ITGFORMAT ",\"worst_fwd\":%.6e,"
         "\"worst_bwd\":%.6e,\"worst_asm\":%.6e}\n",
         iinc,iit,node,dir,h,o->scale,o->relctr,o->relside,o->relbest,
         o->n,o->nok,o->nkink,o->nwrong,o->nboth,o->worst,o->worstdir,
         o->wfwd,o->wbwd,o->wasm);
  fflush(stdout);
}

void monitor_opcheck_total(const opcheck *o,ITG iinc,ITG iit,ITG elem,
                           double dam,double h,ITG ncol){
  ITG bad=o->nkink+o->nwrong+o->nboth;
  printf("[OPCHECK] inc=%" ITGFORMAT " iter=%" ITGFORMAT " element=%"
         ITGFORMAT " dam=%.6f h=%.3e %" ITGFORMAT " column(s), %" ITGFORMAT
         " coefficient(s)\n",iinc,iit,elem,dam,h,ncol,o->n);
  printf("[OPCHECK]   %-8" ITGFORMAT " agree with the tangent\n",o->nok);
  printf("[OPCHECK]   %-8" ITGFORMAT " KINK   - the one-sided differences "
         "disagree and the tangent is on one branch\n",o->nkink);
  printf("[OPCHECK]   %-8" ITGFORMAT " WRONG  - the one-sided differences "
         "AGREE and the tangent is on neither\n",o->nwrong);
  printf("[OPCHECK]   %-8" ITGFORMAT " BOTH   - non-smooth, and the tangent "
         "is on neither branch\n",o->nboth);
  if(bad>0){
    printf("[OPCHECK]   verdict: %.1f%% of the disagreement is a kink the "
           "tangent is right about; %.1f%% is the tangent being wrong where "
           "the residual is smooth\n",
           100.*(double)o->nkink/(double)bad,
           100.*(double)o->nwrong/(double)bad);
  }else{
    printf("[OPCHECK]   verdict: the assembled tangent IS the differential "
           "of the residual here\n");
  }
  printf("[MON] {\"kind\":\"opcheck\",\"inc\":%" ITGFORMAT ",\"iter\":%"
         ITGFORMAT ",\"element\":%" ITGFORMAT ",\"dam\":%.6f,\"h\":%.6e,"
         "\"columns\":%" ITGFORMAT ",\"n\":%" ITGFORMAT ",\"ok\":%"
         ITGFORMAT ",\"kink\":%" ITGFORMAT ",\"wrong\":%" ITGFORMAT
         ",\"both\":%" ITGFORMAT "}\n",
         iinc,iit,elem,dam,h,ncol,o->n,o->nok,o->nkink,o->nwrong,o->nboth);
  fflush(stdout);
}
