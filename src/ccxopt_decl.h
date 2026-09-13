/*     The declarations: what an option is, not merely that it exists.
 *
 *     Hand written, and deliberately not generated - a declaration scraped
 *     from the line that reads it is the same information the tree already
 *     had.  Each entry states the type, the legal spellings, the range, the
 *     default and one line of prose, and docs/SWITCHES.md is generated FROM
 *     here rather than scraped from getenv.
 *
 *     Scope of this first batch: every option that any test in this tree or
 *     test/s3rad/run_s3rad.sh sets - that is, every option whose behaviour
 *     is actually pinned by something.  The other 116 names the binary reads
 *     are reported as UNDECLARED, and that list is the retirement queue.
 *
 *     RANGES record what the code does, including where it does it silently.
 *     CCX_DAMAGE_AUTOSPC is clamped at 1.e-1 inside nonlingeo.c with no
 *     message (05-DEBT.md item 7): a deck whose worst node sits at 1.04e-01
 *     masks nobody and gives no indication why.  Declaring the range does
 *     not fix that; it makes it visible in one place instead of none.
 *
 *     SPELLINGS are recorded as they are, not as they should be.  Four
 *     mutually incompatible notions of "true" are in use here.  Writing them
 *     down is the first step to having one.
 */
#ifndef CCXOPT_DECL_H
#define CCXOPT_DECL_H

typedef enum{
  CCXOPT_BOOL, CCXOPT_INT, CCXOPT_REAL, CCXOPT_ENUM, CCXOPT_STRING
}ccxopt_type;

typedef struct{
  const char *name;
  ITG type;
  const char *dflt;        /* the default, as text                        */
  double lo,hi;            /* inclusive range; lo>hi means unbounded      */
  const char *choices;     /* "A|B|C" for CCXOPT_ENUM, else NULL          */
  const char *doc;         /* one line, and it has to be true             */
  const char *deprecated;  /* replacement or reason, or NULL              */
}ccxopt_decl;

#define CCXOPT_UNBOUNDED 1.,0.

static const ccxopt_decl ccxopt_decl_table[]={

/* ---- the load-path judgement (damstate.c) ------------------------- */
{"CCX_DAMAGE_AUTOSPC",CCXOPT_REAL,"unset (no mask)",0.,1.e-1,NULL,
 "fraction of a node's OWN intact assembled diagonal below which it is "
 "judged to have lost its load path and is excluded from the displacement "
 "residual; silently clamped to 1.e-1 in nonlingeo.c, so a deck whose worst "
 "node is at 1.04e-01 can never exercise it",NULL},

{"CCX_DAMAGE_AUTOSPC_FORCE",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "extend the same judgement to the FORCE residual.  Set to anything.  Prints "
 "the excluded residual next to the criterion: if that number stops "
 "returning to zero the mechanism has become a fiction and is hiding a real "
 "imbalance (02-DIAGNOSTICS.md section 4)",NULL},

{"CCX_DAMAGE_DEADALL",CCXOPT_REAL,"unset (off)",0.,0.5,NULL,
 "a node whose entire live support is dead below this fraction is treated as "
 "having none; clamped to [0,0.5] in nonlingeo.c",NULL},

/* ---- erosion and topology ----------------------------------------- */
{"CCX_DAMAGE_DELETE_MAT",CCXOPT_STRING,"unset (no filter)",CCXOPT_UNBOUNDED,NULL,
 "restrict terminal deletion to elements of these materials; ALL means every "
 "material",NULL},

{"CCX_DAMAGE_DELETE_D",CCXOPT_REAL,"0.999",0.5,0.9999,NULL,
 "the terminal deletion threshold: a live C3D4 of a progressive material "
 "leaves the assembly when its worst integration point reaches this damage.  "
 "Outside [0.5,0.9999] the value is refused with a warning and 0.999 kept.  "
 "READ WITH CCX_DAMAGE_DELETE_VISC: by default the trigger is the VISCOUS "
 "damage, not the instantaneous one, so an element at D=1.0000 whose Dvis is "
 "0.998 stays.  MEASURED 2026-09-13 on the wrapped fast deck: the outcome is "
 "not monotone in this number.  0.9990-0.9950 wall at theta 0.158; "
 "0.9920-0.9850 run to theta 1.0; 0.9800-0.9600 wall again.  The band that "
 "runs to the end is a PHANTOM - the exact minimum cut reaches 2.9e-04 of "
 "its reference by increment 56 and stays there while the solver applies "
 "load for 480 more increments.  Lowering this threshold does not make the "
 "specimen break; it makes the solver stop noticing that it already has.  "
 "ON THE TARGET DECK IT IS WORSE STILL: 0.99 with CCX_FRACTURE_CUT armed "
 "dies at increment 168, theta 0.1914, with the grip reaction at 95.4 "
 "percent of peak and 278 elements deleted, against theta 0.2588 and 3749 "
 "deletions at the default.  The run is killed just past peak load.  The "
 "mechanism is the one written at the trigger itself: an element removed at "
 "Dvis=0.99 is still carrying one percent of its effective stress and "
 "releases it in one increment at constant load, and the residual there "
 "DIVERGES (7.2 -> 16.9) instead of creeping.  Raising this number is a "
 "change to how much force a deletion dumps, not a change to when the "
 "specimen breaks",NULL},

{"CCX_DAMAGE_DELETE_VISC",CCXOPT_BOOL,"1 (on)",CCXOPT_UNBOUNDED,NULL,
 "judge terminal deletion by the VISCOUS damage rather than the "
 "instantaneous one.  On by default and the string 0 is the only value that "
 "turns it off.  The reason it is on: resultsmech.f scales the stress by "
 "1-Dvis whenever the viscosity is armed, so a trigger reading D removes an "
 "element still carrying 1-Dvis of its effective stress and releases that "
 "force in one increment at constant load - measured on DHC1, elements left "
 "while carrying up to 41 percent (batch_Dvis down to 0.592), and correcting "
 "it moved the run from lambda=0.3793 to 0.4651.  The cost is a lag that "
 "GROWS as the step shrinks, beta=dt/(eta+dt): the s3rad trap strands "
 "element 19535 at D=1.0000 with Dvis short of the threshold "
 "(research/16-TRAP-ANATOMY.md).  Turning it off does not clear that wall "
 "either - measured, the wrapped deck still stops at theta 0.1575.  With the "
 "viscosity off damvisc is never allocated and the trigger falls back to D, "
 "so this is a no-op there",NULL},

{"CCX_DAMAGE_TOPOLOGY",CCXOPT_ENUM,"immediate",CCXOPT_UNBOUNDED,
 "DEFERRED|deferred|1",
 "DEFERRED batches topology changes into one transaction committed at the "
 "end of the increment instead of applying each deletion as it is found",NULL},

{"CCX_FRACTURE_TERMINATION",CCXOPT_STRING,"unset (never terminates)",
 CCXOPT_UNBOUNDED,NULL,
 "SETA:SETB - stop the run when no load path remains between these two node "
 "sets.  This is the only thing standing between a run and the phantom "
 "regime of 05-DEBT.md item 1",NULL},

{"CCX_FRACTURE_CUT",CCXOPT_REAL,"unset (off)",0.,1.,NULL,
 "stop when the load path between the CCX_FRACTURE_TERMINATION sets has "
 "narrowed to this FRACTION of the width it had at the first committed "
 "deletion batch.  The width is the minimum cut - the smallest total "
 "(shared face area x residual stiffness) that would have to break to "
 "separate them - so this is a quantity where the other four termination "
 "switches ask a yes/no.  Measured on s3rad: every topological rule says "
 "CONNECTED at the end of a run whose two grips are joined by ONE "
 "triangular face, 0.05 percent of a section, carrying 2.1 percent of "
 "peak load (research/09-SEVERANCE.md).  Off by default and bit-identical "
 "when off; the measure has a self test and refuses to arm if it fails",
 NULL},

{"CCX_FRACTURE_CUT_EXACT",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "with CCX_FRACTURE_CUT armed, disable the early exit so every committed "
 "batch reports the TRUE minimum cut instead of a lower bound.  Costs a "
 "full max-flow per batch - about 1.6 percent of an s3rad run - and is how "
 "the cut trajectory is measured before a stopping fraction is chosen.  "
 "MEASURED EXPENSIVE at s3rad scale: 274 increments in 66 minutes against "
 "599 in 69 minutes without it, so roughly half the throughput.  A "
 "diagnostic, not a production setting - the early exit is what makes the "
 "measure free.  "
 "Without it the reported number is printed as cut>= and ratio>=, because "
 "that is what it is",NULL},

{"CCX_FRACTURE_DEADFACET",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "exclude a cohesive facet whose every integration point has failed from the "
 "load path used by the termination test.  Any value except the string 0 "
 "means on",NULL},

{"CCX_FRACTURE_LINK",CCXOPT_ENUM,"NODE",CCXOPT_UNBOUNDED,
 "NODE|node|FACE|face",
 "what counts as a connection when the termination test walks the live bulk: "
 "sharing a node, or sharing a whole face",NULL},

/* ---- globalization ------------------------------------------------ */
{"CCX_DAMAGE_LINESEARCH",CCXOPT_ENUM,"off",CCXOPT_UNBOUNDED,
 "ADAPTIVE|adaptive|1",
 "arm the adaptive damage line-search ladder (lsladder.c).  Measured cost: "
 "the ladder and the rescues together consume 4.6% of the run "
 "(research/01-PROFILING.md)",NULL},

{"CCX_DAMAGE_REEQ_RESCUE2",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "arm Rescue level 2 with an event step.  Set to anything, INCLUDING 0",NULL},

{"CCX_DAMAGE_REEQ_SCALE",CCXOPT_ENUM,"off",CCXOPT_UNBOUNDED,
 "PHYSICAL|physical|1",
 "scale the re-equilibration step by a physical length rather than by the "
 "residual norm",NULL},

{"CCX_DAMAGE_TR_DOGLEG",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "arm the dogleg trust region as rescue level 3.  Set to anything, INCLUDING "
 "0.  Refuses to arm without CCX_DAMAGE_REEQ_RESCUE2: it is a level on top "
 "of Rescue 2, not a replacement",NULL},

/* ---- path following ----------------------------------------------- */
{"CCX_PATHFOLLOW",CCXOPT_REAL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "dissipation path following: tau per increment.  Must be strictly positive; "
 "the code refuses to arm otherwise",NULL},

{"CCX_PATHFOLLOW_DTHETA",CCXOPT_REAL,"1.e-3",CCXOPT_UNBOUNDED,NULL,
 "the load-factor increment at which path following engages; a non-positive "
 "value falls back to the default",NULL},

{"CCX_CRACK_CONTROL",CCXOPT_REAL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "crack control: the control increment per step, which must be strictly "
 "positive or the mechanism refuses to arm.  Needs UC6 elements present and "
 "its own kinematics self test to pass.  Mutually exclusive with "
 "CCX_PATHFOLLOW_COD; the code refuses to arm if both are set",NULL},

{"CCX_CRACK_CONTROL_ENGAGE",CCXOPT_INT,"0",CCXOPT_UNBOUNDED,NULL,
 "increment at which crack control engages",NULL},

{"CCX_DISSIPATION_CONTROL",CCXOPT_ENUM,"unset (off)",CCXOPT_UNBOUNDED,"1|2",
 "dissipation control mode; 2 assembles the f_hat = dR/dlambda vector.  Inert "
 "unless CCX_DISSIPATION_TARGET is positive",NULL},

{"CCX_DISSIPATION_TARGET",CCXOPT_REAL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "dissipation target per increment; a positive value also turns the "
 "dissipation report on",NULL},

/* ---- the operator ------------------------------------------------- */
{"CCX_DAMAGE_TANGENT",CCXOPT_ENUM,"unset - g(D)*Cep, the measured best",
 CCXOPT_UNBOUNDED,"UNSYM|unsym|2|FD_SYM|fd_sym|1",
 "which tangent to assemble.  Unset is the secant g(D)*Cep and is the "
 "default.  UNSYM (2) adds the consistent rank-1 correction "
 "-sigma_eff (x) dD/d(eps) through the asymmetric assembly; on fast-plain "
 "it produces a byte-identical fracture for 0.2 percent fewer Newton "
 "iterations and 37 percent MORE wall time, so it is off by default "
 "(research/06-TANGENT-VERDICT.md).  FD_SYM (1) was DELETED on 2026-09-11 "
 "and now stops the run with a message rather than being ignored.  The "
 "rank-1 projection itself is verified by damrank1test in src/damrank1.f",
 NULL},

{"CCX_DAMAGE_GMIN",CCXOPT_REAL,"1.e-4",1.e-6,0.2,NULL,
 "the residual-stiffness floor: a fully damaged element keeps this fraction "
 "g=1-D of its tangent until the terminal scan deletes it.  Out-of-range "
 "values are silently replaced by the default in resultsmech.f and "
 "mafilldamas.f, so the range here is the real one.  COST, measured: the "
 "conditioning it produces forces PARDISO into 8062 iterative-refinement "
 "steps over 5344 solves on s3rad, about 5.6 percent of the run "
 "(research/01-PROFILING.md).  BENEFIT, measured: on fast-wrapped raising "
 "it from 1e-04 to 1e-02 REMOVES the wall entirely - the run completes the "
 "load history with the identical 64-element deletion set and half the "
 "Newton iterations per increment - while 1e-06 through 1e-03 are "
 "indistinguishable and 1e-01 is worse than the default.  The useful band "
 "is narrow and the default sits four decades below it "
 "(research/12-GMIN.md, gate case fast-wrapped-gmin)",NULL},

{"CCX_DAMAGE_VISCOSITY",CCXOPT_REAL,"0 (off)",CCXOPT_UNBOUNDED,NULL,
 "viscous regularisation eta for the damage evolution; negative values are "
 "clamped to zero",NULL},

{"CCX_DAMAGE_TANGENT_DUMP",CCXOPT_INT,"0 (off)",CCXOPT_UNBOUNDED,NULL,
 "print, for this element, the two factors of the rank-1 correction - "
 "|dD/d(eps)| and |sigma_eff| - and their product, so the size of the "
 "correction can be compared against the operator error the structural "
 "probe measures at the same point",NULL},

{"CCX_DAMAGE_TANGENT_H",CCXOPT_REAL,"1.e-7",1.e-12,1.e-3,NULL,
 "perturbation of the forward difference that builds dD/d(eps) for the "
 "rank-1 term.  A forward difference carries O(h) truncation and O(eps/h) "
 "roundoff, so the operator error against h is a curve with a minimum; the "
 "shipped value had never been placed on it",NULL},

/* ---- the constitutive law ----------------------------------------- */
{"CCX_UC6_CONTACT_SMOOTH",CCXOPT_REAL,"unset (sharp law)",CCXOPT_UNBOUNDED,NULL,
 "penetration band over which the crack-face closure kink is blended.  The "
 "kink is a measured factor of 1/gmin = 1e+06 in the normal tangent; the "
 "blend is a byte-for-byte no-op wherever no DAMAGED facet closes",NULL},

/* ---- the linear solver -------------------------------------------- */
{"CCX_PARDISO_REUSE_SYMBOLIC",CCXOPT_ENUM,"off",CCXOPT_UNBOUNDED,
 "1|ON|on|YES|yes",
 "keep PARDISO's symbolic factorisation across numerical factorisations, "
 "keyed on a hash of the sparsity pattern.  MEASURED: it retains the "
 "analysis for 603 of 606 factorisations on fast-wrapped and buys nothing "
 "outside run-to-run noise (handover/04-REFUTED.md)",NULL},

/* ---- operator verification ----------------------------------------
   Found, not built: the directional-derivative check that
   07-RESEARCH-AGENDA.md rank 3 calls "the named hole in the diagnostics"
   was already here, under three undocumented names that no test set.  It is
   PETSc's -snes_test_jacobian in all but the name.  Declared rather than
   retired, because a tool nobody could find is not the same thing as a tool
   that does not exist.  What it says: research/03-OPERATOR.md. */
{"CCX_STRUCT_FD_INC",CCXOPT_INT,"0 (off)",CCXOPT_UNBOUNDED,NULL,
 "arm the operator check at this increment: compare the ASSEMBLED tangent, "
 "column by column, against a central difference of the internal force.  The "
 "run stops afterwards - the probe perturbs the displacement repeatedly, so "
 "the run is diagnostic only",NULL},

{"CCX_STRUCT_FD_ITER",CCXOPT_INT,"0",CCXOPT_UNBOUNDED,NULL,
 "the Newton iteration at which the operator check runs",NULL},

{"CCX_STRUCT_FD_STEP",CCXOPT_INT,"0 (any step)",CCXOPT_UNBOUNDED,NULL,
 "restrict the operator check to this *STEP.  iinc restarts at 1 in every "
 "step, so without this the probe can only ever fire in the first one - "
 "which is the wrong one whenever the interesting state is reached by "
 "unloading, as the crack-face closure benchmark is",NULL},

{"CCX_STRUCT_FD_ELEM",CCXOPT_INT,"0 (the most damaged C3D4)",
 CCXOPT_UNBOUNDED,NULL,
 "probe this bulk element instead of the most damaged one, so a measured "
 "discrepancy can be attributed to a chosen population of the rank-1 "
 "census rather than to whichever element happened to be worst",NULL},

{"CCX_STRUCT_FD_BASE",CCXOPT_ENUM,"V",CCXOPT_UNBOUNDED,"V|VOLD|vold",
 "which state the operator check differentiates around.  V is the iterate "
 "the residual was last evaluated at, which is one Newton step AFTER the "
 "state the matrix was assembled from; VOLD is that state itself.  The "
 "difference between the two settings is the size of the "
 "offset-of-one-iterate objection to any discrepancy the check reports",NULL},

{"CCX_STRUCT_FD_H",CCXOPT_REAL,"1.e-7",CCXOPT_UNBOUNDED,NULL,
 "the perturbation of the central difference.  Measured: the answer is flat "
 "over 1e-11 to 1e-8 and has moved by 1e-3, so the default sits inside the "
 "converged plateau",NULL},

/* ---- diagnostics that were in the retirement queue only because nobody
   had written them down.  Each answers a question 02-DIAGNOSTICS.md asks a
   human to answer by eye. */
{"CCX_DAMAGE_TANGENT_CENSUS",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "per-iteration census of the damage tangent: how many elements got the "
 "rank-1 term -sigma_eff (x) dD/d(eps), how many with ADVANCING damage did "
 "not, and how many were skipped.  This is the report that names WHY the "
 "assembled tangent is not the differential of the residual "
 "(research/03-OPERATOR.md); it was in the retirement queue",NULL},

{"CCX_DAMAGE_NODE_DUMP",CCXOPT_INT,"0 (off)",CCXOPT_UNBOUNDED,NULL,
 "dump one named node every increment: every element that touches it, its "
 "type, whether it is still assembled, the damage at each integration point "
 "of a cohesive element, and the assembled diagonal beside its intact "
 "reference.  Global counters have twice disagreed with each other here; "
 "this asks the question directly instead",NULL},

{"CCX_DAMAGE_NODE_INC",CCXOPT_INT,"1",CCXOPT_UNBOUNDED,NULL,
 "the increment from which CCX_DAMAGE_NODE_DUMP starts printing",NULL},

{"CCX_DAMAGE_STIFF_PROBE",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "print the stiffness census - how many nodes have lost what fraction of "
 "their OWN intact assembled diagonal (02-DIAGNOSTICS.md section 3).  Armed "
 "automatically whenever CCX_DAMAGE_AUTOSPC is",NULL},

{"CCX_DAMAGE_STIFF_MIN",CCXOPT_REAL,"0 (off)",0.,0.1,NULL,
 "report nodes whose assembled diagonal has fallen below this fraction of "
 "their intact value; silently clamped to [0,0.1] like CCX_DAMAGE_AUTOSPC",
 NULL},

{"CCX_DAMAGE_AUTOSPC_NEG",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "count a NON-POSITIVE assembled diagonal as dead as well.  The census "
 "reports nonpositive= separately, and on every deck measured so far it is "
 "zero",NULL},

/* ---- measurement -------------------------------------------------- */
{"CCX_PARDISO_REPACK",CCXOPT_BOOL,"ON by default; set 0 to disable",
 CCXOPT_UNBOUNDED,NULL,
 "cache the CSR permutation and refill only the values when the sparsity "
 "pattern has not changed.  The structurally-symmetric asymmetric branch "
 "(mtype=1), which is what CCX_DAMAGE_TANGENT=UNSYM makes every run of this "
 "branch take, otherwise rebuilds its whole CSR on EVERY factorisation - "
 "four allocations, a full sort of the lower triangle by row, a sort of "
 "every row by column, and an interleave - measured at 73 ms a call and "
 "9.4 percent of an s3rad run, 6.5 minutes of 69.  The permutation depends "
 "only on the pattern, and the hash that detects a pattern change already "
 "runs on every call: it fires on 165 of 5344.  Acceptance is BIT "
 "IDENTITY, because the same values land in the same slots - verified on "
 "both fast decks and on s3rad, where the deletion record is byte-for-byte "
 "identical, all 3741, for 11.7x less repacking and 14 percent less wall "
 "clock.  It only engages where symbolic reuse is already eligible",NULL},

{"CCX_PARDISO_REPACK_VERIFY",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "with CCX_PARDISO_REPACK armed, rebuild the whole CSR anyway on every "
 "reuse call and compare it slot by slot against what the cached "
 "permutation produced.  This is the only check that can actually fail - "
 "it compares the shortcut against the thing it short-cuts, on real "
 "matrices.  Expensive by construction: it does both",NULL},

{"CCX_PARDISO_REPACK_BREAK",CCXOPT_INT,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "corrupt one entry of the cached permutation so that the verify check can "
 "be SEEN to go red.  A check that has never been observed failing is not "
 "a check.  DIAGNOSTIC ONLY - it makes the matrix wrong",NULL},

{"CCX_PARDISO_REFINE_EVERY",CCXOPT_INT,"200 (solves between reports)",
 CCXOPT_UNBOUNDED,NULL,
 "how often [PARDISO REFINE] prints the running count of iterative "
 "refinement steps PARDISO actually performed - iparm(7), which nothing read "
 "until 2026-09-11. A cumulative line every N solves so that consecutive "
 "lines can be differenced the way the LOGVIEW interim tables are. A "
 "non-positive value reports only at exit",NULL},

{"CCX_LOG_VIEW_EVERY",CCXOPT_REAL,"600 (ten minutes)",CCXOPT_UNBOUNDED,NULL,
 "seconds between INTERIM profile reports; 0 prints only at exit.  Its "
 "defender: two 2.3-hour runs of the target deck were killed part way "
 "through and produced no profile at all, because the table was printed "
 "from atexit.  A profiler that reports only at the end is useless on "
 "exactly the runs it exists for",NULL},

{"CCX_CONVERGE_EXPLAIN",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "print the convergence verdict as a table of named clauses - each with "
 "the value it tested, the threshold it was tested against and pass/FAIL "
 "- and name the clause holding the increment back.  Its defender: the "
 "criterion is eight clauses of && and || written out three times, and "
 "until this existed nothing could say WHICH one an increment was failing. "
 "Diagnostic only; the verdict does not read it and the run is "
 "bit-identical with it on",NULL},

{"CCX_LOG_VIEW",CCXOPT_BOOL,"unset (off)",CCXOPT_UNBOUNDED,NULL,
 "print where the run spent its time: named events with a call count, an "
 "inclusive and a self time, and who called whom.  Measurement only; the run "
 "is bit-identical with it on",NULL},

};

#define CCXOPT_DECL_COUNT ((ITG)(sizeof(ccxopt_decl_table)/sizeof(ccxopt_decl_table[0])))

#endif
