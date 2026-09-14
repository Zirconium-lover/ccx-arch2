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

Where the shapes come from. None of this is invented here, and the files
name their sources the way `topology.c` already named PETSc's `DMLabel`
and `globalize.c` named `SNESLineSearchReason`:

- `trialctx` is PETSc's SNES application context (`SNESSetFunction`'s
  `ctx`) and deal.II's `ScratchData` — both exist so that an evaluation
  routine does not take many individual arguments, and both are documented
  as the way those libraries avoid global state.
- the policy / workspace / census split inside `dogleg` and `damcont` is
  the ordinary shape of a trust-region or continuation object: what it is
  configured to do, what it works in, what it has spent.
- keeping the **guard** with the caller and moving only the body is what
  lets a mechanism stay opt-in without the object having to know about
  increments at all.

Two things that are **not** reasons to make an object: that some lines are
long, and that some lines look alike. Duplication is evidence, not a
verdict; the forty-line deletion block appearing twice mattered because
both copies were one decision, not because they were identical.

## The other half of the rule

Everything above is bottom-up: it says what may be pulled out of a long
function. It was enough to get twenty-three files out of `nonlingeo()` and
it is not enough to keep them worth having. A tree can satisfy all four
criteria on every file and still be miserable to work in, because the four
say nothing about what it COSTS to change one. Four more, and each is a cost
somebody pays on every change rather than a matter of taste:

5. **A signature is an interface, or it is a copy of the caller's locals.**
   A function taking twenty-two arguments has no boundary: every local the
   caller adds is an argument it grows. The budget is six, and when a
   signature will not fit, the thing to look for is the missing object, not
   a shorter name. Measured here, the lists got long for three separable
   reasons, with three different fixes: scalars DERIVED from the model and
   passed by hand (`mi[0]`, `mi[1]+1` — thirty-three hand-offs across the
   worst twenty-seven), iteration state that no object owns (`ram`, `uam`,
   `qam`, `iit`), and model state the evaluator's context never learned
   about (`dambase`, `dmcon`). None of the three is fixed by editing a
   signature.

6. **Dependencies go one way, and the way is written down.** Six layers,
   listed in `tools/arch.py`; a module may call strictly down. The table was
   not designed — it was derived from the call graph that already existed,
   and it is checked rather than believed. This matters most for the files
   that are meant to be harmless: an observer that can only call downward
   cannot change an answer, and that is a stronger statement than
   `damdiag.c`'s promise not to.

   The rule earned its keep the hour it was written. `topology.c`, which
   says what the model IS, was calling `damage_progressive_material()` in
   `erosion.c`, which decides what LEAVES the model — a question about a
   material card, answered inside a deletion policy because deletion was the
   caller that got written first. `erosion.c` had already had the right
   instinct and applied it to the wrong thing: the constant beside that
   function carries a comment saying it was renamed "under the name of the
   question it answers rather than of the one caller that happened to be
   written first". The name moved; the file did not. It lives in
   `dammat.c` now, below both callers.

7. **A test that needs a deck is a test that runs once a day.** All
   nineteen self tests used to be compiled into the solver and reached only
   from inside `nonlingeo()`, so proving a forty-line classifier cost a
   finite element analysis — and every production job paid for it too, 56 to
   82 lines of self test output before increment 1. `ccx_selftest` runs all
   nineteen in 14 ms now, and `selftest_gate()` keeps the interlock from
   `main()`; production output is one line. Criterion 4 above says an
   object must be able to fail a test; this one says the test has to be
   cheap enough that it actually got run. `dammat_selftest()` found two
   bugs in its own first hour — an offset that made three assertions
   unfailable, and a missing case that let `nconst==4` relax to `nconst>=4`
   unnoticed — because it could be compiled and run in a second.

8. **The rebuild is part of the interface.** Changing any extension
   declaration used to recompile 208 of 208 `.c` files, because they all
   lived in `CalculiX.h`. That is the price of every experiment, paid by
   whoever tries anything, and it is not a question of style. They live in
   `ccxfork.h` now — 49 files — with `ccxopt.h` and `logview.h` split off
   because they are platform rather than mechanism.

### The limit of a context, found by measuring

Nine of the twenty-nine widest functions are called by a self test, and they
include the widest three — `erosion_mark` (22), `loadcut_width` (22),
`converge_norms` (21). That is not a coincidence. A self test builds a
synthetic four-element mesh and three materials; it cannot build a 180-field
`trialctx` bound to the locals of a running solver, and it should not have
to. **A function that takes plain arrays is a function you can test; a
function that takes the context is one you can only run.**

So the rule is not "pass the context everywhere". It is:

> Convert at the boundary the driver calls. Keep a plain-array core wherever
> a test needs one. Where both are wanted, the context-taking function is a
> thin wrapper over the testable core.

`tools/arch.py` knows this — it exempts any public function a self test
calls, and prints how many. Without that exemption the metric would have
rewarded deleting the only checks on the three widest functions in the tree,
which is the exact opposite of what it exists for. A metric that pushes
against criterion 4 is a worse metric than none.

### Where these numbers are now

| | at the start | now |
|---|---|---|
| widest public signature | 22 | 22 |
| signatures over the 6-argument budget | 36 | 35 |
| files recompiled by an interface change | 210 | 49 |
| self tests runnable without a deck | 0 of 19 | 19 of 19 |
| self test lines in a production run | 56–82 | 1 |
| layering violations | 1 | 0 |
| `nonlingeo()` lines / locals | 11,313 / 631 | 11,252 / 593 |

The first row is the honest one: the signatures are barely touched. One of
the twenty-seven is converted and it was the cheapest. What moves that
number is the three contexts — `trialctx` existed, `nlstate` now exists, and
the derived-scalar accessors do not.

Everything above that row is infrastructure, and it was taken first on
purpose: `ccx_selftest` pays for itself on every later step, and an
interface change that recompiles 49 files instead of 210 is what makes
trying one affordable.

Where these come from, in the same spirit as the shapes above: the layering
rule is the ordinary package-dependency discipline (Martin's stable-
dependencies principle, stated as a table rather than as advice); the
signature budget is why PETSc passes one `ctx` instead of an argument list
and why deal.II's `WorkStream` passes `ScratchData`; and the separate test
binary is the arrangement every one of those libraries already has and this
tree does not.

### Why these four and not others

Each has a direction that is not arguable. Nobody wants a wider signature, a
dependency cycle, a slower test or a bigger rebuild, so each can be a
ratchet — `tools/arch.py --check` fails when any of them drifts up. Numbers
that trade off against each other cannot be ratcheted, and a number that is
allowed to drift back up is not a measurement, it is a mood.

The layering entry is a floor rather than a ratchet: the table was adopted on
a tree with zero violations, so any violation at all is new.

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
- **`dogleg.c`** — *the trust region.* `dogleg_pick()` is the step choice
  alone: a radius and five scalars in, a branch and two coefficients out,
  which is what makes it testable against a closed form.
  `dogleg_rescue()` is the loop that uses it.
- **`lsladder.c`** — *what step length to try next, and which to take.*
- **`loadcut.c`**, **`crackcontrol.c`**, **`pathfollow.c`** — the load
  parameter: when to cut it, how to drive it by crack opening, how to
  follow the dissipation.
- **`topology.c`** — *the transaction that commits an erosion.* One marked
  set, one lifetime, one discard.
- **`loadctl.c`** — *who drives the load parameter.* Dissipation control,
  path control, the arc-length boundary and the regularisation ladder: 53
  locals, four mechanisms, and one rule — they refuse to run beside each
  other. That refusal was written out as four loose integers in every
  other mechanism's arming block; `loadctl_driving()` is the one question
  they ask now.

### State with an owner, and drivers partly still in `nonlingeo()`

These files hold a cluster's state and its arithmetic. Where a loop has
moved it is named; where it has not, it is because the loop cuts
increments, rolls topology back or ends steps.

- **`rescue.c`** — what happens when an increment will not converge: the
  line search, its probe, transactional backtracking, same-load
  re-equilibration, the recovery corridor and the levels that order them
  (73 locals, seven clusters, one object). `rescue_backtrack()` is the
  backtracking loop; `rescue_attempt()` is the decision about which level
  gets the one extra attempt a deferred stop buys. The three
  `ccx_rescue_*` flags — the ladder's handshake with
  `checkconvergence.c`'s verdict — are defined here too.

  `rescue_attempt()` takes **thirteen** parameters and the count is
  deliberate: six objects, and seven pieces of `nonlingeo()`'s own control
  flow. Bundling those seven would shorten the signature, leave the
  coupling exactly as it is, and put stock CalculiX locals into a
  fork-specific type — which costs something real on the next upstream
  merge. The length is the measurement; it is a number to reduce, not to
  hide.
- **`damcont.c`** — the bounded local continuation (88 locals → one
  object; 12 of the 88 were declared and never read, and are gone). The
  driver is still inline.
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

## Arming: each object reads its own switches

The configuration was 141 `ccxopt_getenv` calls in one function. Most of
them now sit in the file that owns what they configure, as a
`<module>_configure()` — `dogleg_configure()`, `damcont_configure()`,
`rescue_configure_backtrack()`, `rescue_configure_levels()`,
`probedrv_configure_aba()`, `opcheckdrv_configure_fd()`. This is PETSc's
`XXXSetFromOptions`: one per object, called in the order the objects are
created.

Two properties of these blocks decide how they had to move.

They **refuse rather than degrade**, and they refuse *on each other's
state*: the dogleg will not arm without Rescue2, the continuation will not
arm without the dogleg, and neither will arm while anything is driving the
load parameter. So the order in which they run is part of the behaviour,
and every `_configure()` call sits at exactly the line its block occupied.
Nothing about which switch is read first has changed.

They **write across objects**: the corridor sets the regularisation
ladder's length, two rescue levels reset the event census. The signatures
say so — `rescue_configure_levels(rescue*, loadctl*, probedrv*)`, with the
two it writes to non-const. A signature that hid that would be worse than
the inline block.

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

Two loops have moved, as the proof that it can be done.

The trust-region dogleg was 239 lines naming twenty things outside itself;
before `trial.c` and the `dogleg` object, thirteen of those twenty were raw
locals of `nonlingeo()` and the other seven were the thirty-line
`results()`/`calcresidual()` pair written out by hand, so there was no
signature to give it. There is one now, of seven arguments, and the body is
`dogleg_rescue()` in `dogleg.c`. Transactional backtracking followed: 148
lines, nineteen of its twenty-seven external names now arriving through
`trialctx`, and `rescue_backtrack()` takes six arguments.

In both cases the **guard** stayed behind. Whether the mechanism may fire
at all — not thermal, not dynamic, no contact, no continuation running — is
a decision about the increment, and it belongs where the increment is. What
moves is what the mechanism then does. That split is the pattern for the
loops that remain.

Moving a block also makes the compiler able to see it. `-Wall` on the small
functions found an out-of-bounds read in a self test's own data, two unused
locals, and a `may be used uninitialized` that is a false positive — and
that last one is worth the note it carries in `rescue.c`, because the
identical lines produced no warning while they were inside a 13,000-line
function. That was checked by compiling the pre-move file, not assumed.

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

   Two things this check has been wrong about, both of them the check's
   fault and not the code's, both now fixed in `tools/abruns.py`: it
   compared the `.frd` `UDATE` record, so a comparison run either side of
   midnight reported sixteen differing files, one per case, all identical
   in length; and the build it runs against only rejected `error:`, so an
   integer passed where a `const loadctl *` was expected compiled with a
   warning and killed eleven of sixteen cases with SIGSEGV. The build now
   carries `-Werror=int-conversion -Werror=incompatible-pointer-types` —
   two classes, both of them the ones C lets through silently, and the
   whole tree including stock CalculiX builds clean with them.

3. **The ratchet.** `tools/arch.py --check` fails when any measured number
   — the function's length, its locals, or a subsystem's site count — is
   worse than `tools/arch_budget.json`. A number allowed to drift back up
   is not a measurement.

### What the gate covers, and what it does not

Of the loops that have moved, three are executed by gate cases and one is
not, and the difference is worth stating case by case rather than letting
"144 files identical" carry weight it has not earned:

| moved loop | exercised by | evidence |
|---|---|---|
| `dogleg_rescue()` | `fast-wrapped` and 2 more | 79 `[DAMAGE TR]` lines, same before and after |
| `rescue_backtrack()` | `fast-wrapped` | 66 `[DAMAGE BT]` lines |
| `rescue_attempt()` | 3 cases | 8 `FIRED`, 2 `ACCEPTED` each |
| `pathdrv_predictor()` | both `mixed-analytic` cases | 5000 accepted increments |
| `opcheck_probe()` | nothing — **now `fast-plain-opcheck`** | 47 probe lines diffed by hand across the move, then given a case |
| `damcont_corrector()` | **nothing, and it cannot be** | see below |

`damcont_corrector()` is the honest bad case. Level 4 arms on every deck in
this tree and is then refused by `damcont_select()`; raising the curvature
tolerance only moves the refusal to another clause. The mechanism was built
for a target deck that is not here. What was checked is the arming and
selection path, which does run: the whole `[DAMAGE CT]` output of the
wrapped deck, 28 lines, identical across the move. The body is covered by
the compiler and by reading, and `damcont.c` says so at the function.

Moving it was still right: the code is untestable either way, and after the
move it is untestable *in the file that owns it* rather than untestable in
the middle of a thirteen-thousand-line function.

There is a limit to what check 2 can say, and it has to be said out loud:
**byte identity over the gate proves nothing about code the gate does not
execute.** 133 of the 153 switches are set by no case, so every mechanism
behind one of them is moved on the strength of the compiler and of reading.
When a mechanism like that is moved, run it directly, before and after, and
diff its own output — that is how the operator check was verified (47
`[OPCHECK]` lines identical across the move) — and then, if the mechanism
is worth keeping, give it a case so the next person does not have to. The
operator check has one now: `fast-plain-opcheck`, asserting 27 columns, 0
wrong and 0 kink, and observed going red when the expectation is moved.

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
