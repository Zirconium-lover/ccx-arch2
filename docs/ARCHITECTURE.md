# The architecture of the damage/fracture solver

This describes `src/nonlingeo.c` and the files around it: what the objects
are, where the boundaries run, what is still unowned and why, and how a
change to any of it is accepted.

Numbers in this document are produced by `tools/arch.py`, not asserted.
Run it.

## The problem, as a measurement

`nonlingeo()` is the nonlinear solve: the Newton loop, the step control,
the convergence verdict, the removal of failed elements, the stabilisation
mechanisms and the diagnostics. This fork added a damage and fracture model
to it, and the model was added *inside* it.

At the start of this work:

| | |
|---|---|
| `src/nonlingeo.c` | 16,120 lines |
| `nonlingeo()` itself | 14,402 lines |
| locals that function declares | **1,068** |
| of them, stock CalculiX's | 479 |
| of them, this fork's | **589** |
| other functions in the same file | 32, all file-static |

The locals are the number that matters. Each one is reachable and writable
from every line of the function, so *every* piece of state is in scope for
*every* decision. That is what "there is no architecture" means here in a
form you can count, and it is also the mechanism: a block cannot become a
function while it names two hundred locals, so nothing could be moved out,
so everything new was written where everything else already was.

The second measurement is locality. For each named subsystem — the tags the
run prints, `[DAMAGE TR]`, `[WALLDIAG]`, `[PATHFOLLOW]` — `tools/arch.py`
counts how many separate **sites** in `nonlingeo.c` a reader has to find,
and how far apart the first and last are:

```
  DAMAGE TR   20 sites over 15,498 lines    no file owns it
  DAMAGE CT   18 sites over 15,144 lines    no file owns it
  WALLDIAG    12 sites over  9,923 lines    no file owns it
  LSLADDER     0 sites                      lsladder.c
```

The last row is the target the others are measured against. A subsystem
with a file has no sites, and "where is this decided" has an answer.

## The rule

Something is an object here when it can answer one question, and answer it
alone. Concretely, all four of:

1. **It names a question somebody asks.** "Which elements leave the
   assembly?" "Did this increment converge?" "How long a step?" A file
   named after a question is findable; one named after a layer is not.
2. **Its state has one lifetime.** The fields are created together, they
   are meaningless apart, and one call initialises them. If half the fields
   make sense without the other half, it is two objects.
3. **Its boundary is stateable in one sentence, including what it does
   *not* do.** `erosion.c` marks; it does not rebuild the equation
   structure. `damdiag.c` observes; it cannot change an answer.
   `trial.c` evaluates; it does not decide. The "not" half is the half that
   stops the object growing.
4. **It can be wrong in a way a test can catch.** Every object that makes a
   decision carries a self test that runs before the first increment, and
   the self test has been *observed failing* — see "How a change is
   accepted" below.

Two things that are **not** reasons to make an object: that some lines are
long, and that some lines look alike. Duplication is evidence, not a
verdict; the forty-line deletion block appearing twice mattered because
both copies were one decision, not because they were identical.

## The objects

Each file answers the question in its own first paragraph. Sizes are
`wc -l`.

### The atom everything is built from

- **`trial.c`** — *evaluate the residual at a trial state.* Put a step in
  `b`, evaluate the model, read the residual. Every globalisation mechanism
  in this solver is built from it, and it had been written out by hand
  seventeen times: 28 calls to `results()` of which 25 had a
  character-identical 130-argument list, 20 calls to `calcresidual()` of
  which 19 differed only in the destination array. `trialctx` is not an
  abstraction of the model — it is that argument list written down once,
  every field holding the *address* of a caller's local so the binding
  survives every `NNEW`, `SFREE` and `remastruct`. `trial_check()` walks
  it at arming and reports 180 fields, 0 unbound.

  This one comes first because it is what unblocked the rest. A block that
  needs a residual needs those 130 locals; while the only way to name them
  was to be inside `nonlingeo()`, every mechanism had to be written inside
  `nonlingeo()` too.

### Decisions

- **`erosion.c`** — *which elements leave the assembly, and why.* Three
  rules (terminal damage, DEADALL, DEADSOLE) against one batch budget.
  `erosion_mark()` is the single entry point; the three rules are
  file-static, so "what else can delete an element" is answered by
  construction. It marks and does no more.
- **`converge.c`** — *the numbers the convergence judgement is made from.*
  Owns the reduction into `ram`/`uam`/`qam` and which dofs may contribute.
  It decides nothing: `checkconvergence.c` still owns the verdict.
- **`slownewton.c`** — *how many Newton iterations this increment may
  have.* A bounded extension, conditional on damage actually advancing.
  It extends a budget; it relaxes no tolerance and closes no divergence
  path.
- **`dogleg.c`** — *the trust-region step.* `dogleg_pick()` is the choice
  alone: a radius and five scalars in, a branch and two coefficients out,
  which is what makes it testable against a closed form.
- **`lsladder.c`** — *what step length to try next, and which to take.*
- **`loadcut.c`**, **`crackcontrol.c`**, **`pathfollow.c`** — the load
  parameter: when to cut it, how to drive it by crack opening, how to
  follow the dissipation.
- **`topology.c`** — *the transaction that commits an erosion.* One marked
  set, one lifetime, one discard.

### State with an owner but a driver still in `nonlingeo()`

These files hold a cluster's state and its arithmetic; the loop that
orders them is still in `nonlingeo()`, because it cuts increments, rolls
back topology and ends steps.

- **`damcont.c`** — the bounded local continuation (88 locals → one
  object; 12 of the 88 were declared and never read, and are gone).
- **`rescue.c`** — what happens when an increment will not converge: the
  line search, its probe, transactional backtracking, same-load
  re-equilibration, the recovery corridor and the levels that order them
  (67 locals, six clusters, one object).
- **`pathfollow.c`** — `pathdrv`, the driver state of the path follower
  (63 locals), beside the method it drives.

### Observation

- **`damdiag.c`** — *the probes. They look; they never touch.* Ten pure
  functions of state handed to them explicitly, plus `probedrv`, the
  driver state (50 locals). A probe that cannot change an answer cannot
  change one by being moved, which is why this extraction was free.
- **`opcheck.c`** — is the assembled tangent the differential of the
  residual? Plus `opcheckdrv` (39 locals).
- **`globalize.c`** — which globalisation mechanism actually did anything.
  Six mechanisms were stacked in a fixed order with no interface; this
  turns "two rescue levels ran and changed nothing" from an anecdote into
  a counter.
- **`topodiag.c`**, **`damstate.c`**, **`damstats.c`**, **`logview.c`**,
  **`ccxopt.c`** — the topology/rank report, the load-path judgement, the
  damage census and its VTK series, the profile, the switch registry.

## What is still unowned, and why

`tools/arch.py` prints this; it is not a list to be kept by hand. At the
time of writing, `nonlingeo()` still declares about two hundred of this
fork's locals, in clusters of thirty and fewer: dissipation control, path
control, the release probe, the topology-diagnostic handles, the free/float
stabilisation knobs, the termination test's state.

More importantly, *loops* are still inline, and for the same reason in each
case: they do not only compute, they **commit**. The rescue ladder and the
continuation driver decide a cutback, roll topology back, or end a step.
Moving one of those out means giving `nonlingeo()`'s control flow —
`icutb`, `idamagereeq`, `dtheta`, the increment's rollback baselines — an
owner too. That object has not been built.

One loop has moved, as the proof that it can be done. The trust-region
dogleg was 239 lines naming twenty things outside itself; before `trial.c`
and the `dogleg` object, thirteen of those twenty were raw locals of
`nonlingeo()` and the other seven were the thirty-line
`results()`/`calcresidual()` pair written out by hand, so there was no
signature to give it. There is one now, of seven arguments, and the body
is `dogleg_rescue()` in `dogleg.c`. The **guard** stayed behind: whether
the region may fire at all — not thermal, not dynamic, no contact, no
continuation running — is a decision about the increment, and it belongs
where the increment is. That split is the pattern for the loops that
remain.

## How a change is accepted

Three checks, in this order. None of them is optional and none of them is a
matter of opinion.

1. **The gate.** `CCX_EXE=... test/regress/run.py` — sixteen cases in about
   six minutes on one core, covering the load-path judgement, the
   crack-face kink, bulk damage with deletion and cutbacks, and the
   analytical mixed-mode branch. It checks measured scalars against
   `cases.json`, byte identity between cases that must agree, that every
   switch a case asked for actually reached the binary, and that every
   required self test reported PASSED.

2. **Byte identity against the binary before the change.**
   `tools/abruns.py A_RUNDIR B_RUNDIR` compares two whole gate runs file by
   file: 16 cases, 144 files. Exactly two things may differ, both clocks
   rather than arithmetic — the `UTIME` record of a `.frd`, and
   `provenance.txt`, which records the binary's hash on purpose. Every
   extraction described here was accepted on **0 differ**.

   An extraction that changes an answer is not an extraction. If a change
   is *meant* to move a number, this check is the wrong one and the number
   needs a measurement of its own; that is what `cases.json` records.

3. **The ratchet.** `tools/arch.py --check` fails when any measured number
   — the function's length, its locals, or a subsystem's site count — is
   worse than `tools/arch_budget.json`. A number allowed to drift back up
   is not a measurement.

For a new object, two more:

- **The field set must match exactly.** When a cluster of locals becomes a
  struct, the set of fields declared and the set used must be equal in both
  directions. This is what catches the two mistakes this transformation
  actually makes, both of which compile cleanly: a function name swept up
  by the prefix rename, and two different clusters' `mode` merged into one
  field. It is also what found the twelve continuation locals and the one
  tangent-census local that were declared and never read.
- **The self test must have been observed failing.** A check nobody has
  seen go red is not a check. `erosion_selftest()`'s material-filter case
  was written, passed, and was then *proved worthless* by removing the
  length comparison it claimed to guard — the test survived the mutation.
  The case was rewritten to test the other direction, the mutation was run
  again, and it went red. That sequence is the standard, not the anecdote.
