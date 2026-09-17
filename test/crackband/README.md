# Crack-band width: what is here and how to check it

Eleven files, no obvious entry point, so this is the map.  Every number below
is stamped with the commit it was measured on and the command that
reproduces it, because a number written into prose drifts silently when
somebody else's change moves it.  That failure has already happened once
in this tree, to a case description of mine, so the numbers here are quoted
*with* their command rather than asserted.

## The defect

`calcdamage.f` derived the characteristic length from the element volume:

```fortran
det6v=dabs(ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx))
if(det6v.gt.1.d-30) charlen=det6v**(1.d0/3.d0)
```

Two problems, and they are separate.

**It is the wrong length.**  Crack-band theory fixes the dissipation per
unit volume at `G_f/L`, where `L` is the width of the band measured
*across* it (Bazant and Oh 1983).  `(6V)^(1/3)` is a volume-equivalent
size, not a width across anything, and it has no idea which way the band
runs.  Jirasek and Bauer 2012 evaluate exactly this family in section
5.3.1 and conclude the volume-based estimate "cannot be recommended".

**It only existed for one element family.**  A volume-equivalent cube edge
is at least defined for any shape, but the code accepted C3D4 alone and
called `exit(201)` for everything else, so hexahedra could not run the
damage model at all.

## The replacement

The width is the projection of the element onto the major principal strain
direction, taken at the **element centre**, frozen at damage initiation
(Jirasek and Bauer 2012, section 7).  Three details of that sentence are
load-bearing:

* **major principal direction**, because the band runs normal to it, so
  this is the only direction in which "across the band" means anything;
* **at the element centre**, not per Gauss point - the same section reports
  that per-point projections give "excessively large or even infinite
  estimates";
* **frozen at initiation**, because the band's width is a property of the
  band, and the principal direction keeps rotating after it forms.

No orientation factor is applied; section 7 does not use one.

Selected with `CCX_DAMAGE_CHARLEN` (`0` = legacy, `1` or `PROJ` =
projection).  The legacy path is still the default and still refuses
non-tetrahedra, so nothing that ran before changes answer unless the
switch is set.

## Files

| file | what it is |
|---|---|
| `mkcross.py` | generates the objectivity decks - a fixed specimen whose **cross-section** is refined, so the band stays one slice thick while the mesh changes |
| `check_objectivity.py` | compares the runs of a sweep: spread of the external work, the fit to the `1/L` law, and the pointwise spread of `F(u)` split at the peak |
| `run_objectivity.sh` | drives the sweep: C3D4 on both arms, C3D8 on the projection arm only |
| `cbwidth_test.f` | eight hand-computable geometries put through `damcbwidth` directly |
| `run_cbwidth_test.sh` | builds and runs that unit test against `ccx_2.23.a` |
| `run_energy_equiv.sh` | checks that `EVOLUTION=ENERGY` and `EVOLUTION=DISPLACEMENT` agree when `u_f = 2 G_f / sigma_0` |
| `run_nlobjectivity.sh` | the other sweep - refines **along** the axis at fixed `ell`, so the band must choose its own width, which is the sweep an internal length has to be objective against |
| `run_nlwidth_scaling.sh` | holds the mesh and moves `ell`, which separates what the refinement sweep cannot: whether the width follows the internal length, and whether the dissipation follows with it |
| `bandwidth.py` | the band width as `sum(D*V)/A`, with no threshold - it replaced a count over `D>0.5` that gave three verdicts for the three thresholds `de1stats` writes |
| `run_nlwidth_gate.sh` | runs that sweep with `CCX_DAMAGE_NLWIDTH` off and on and requires the `ell` dependence to come down; a relative criterion, so there is no ceiling to pick |
| `run_xsection_floor.sh` | varies the bar's cross-section to find what stops the local band collapsing - it is the specimen, not the code |

## Why the sweep refines the cross-section

Refining *along* the band changes the band's width, so a sweep that does
it is measuring two things at once and cannot separate them.  Refining the
cross-section holds the width fixed - one slice, thickness 1.0, on every
mesh - and varies only how that slice is discretised.  Under the
projection `L` is then literally the same number on all three meshes,
which makes the sweep a clean test of everything *except* the width.

That property is also a trap, and I fell into it: I attributed the
projection arm's residual spread to the band width, which is impossible,
since a quantity that does not move across the sweep cannot produce a
spread across it.  `check_objectivity.py` now says so in place of the fit
whenever `L` is constant.

## Reproducing

```
CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_objectivity.sh /tmp/sweep
```

Measured on `c9f7ab4`:

| arm | `L` across the sweep | spread of `W` | worst departure from `W = A + B/L` |
|---|---|---|---|
| C3D4 legacy | 1.000 -> 0.630 -> 0.481 | 73.57 % | 5.31 % |
| C3D4 projection | 1.000 on every mesh | 5.97 % | law predicts zero |
| C3D8 projection | 1.000 on every mesh | 0.71 % (0.72 % re-measured) | law predicts zero |

Read those three rows together rather than separately.

The legacy row is the case *for* the change: a 73.57 % spread that follows
`1/L` to within 5.31 %, so the error really is the length and replacing
the length really is the fix.  The law accounts for the bulk of it and not
all of it, which is what a 5.31 % departure means and all it means.

The projection rows are **not** evidence that the remaining 5.97 % and
0.71 % are a better-estimated width.  They cannot be: `L` does not move on
those arms.  Whatever causes them is something else - on a deck where
every slice carries the damage law and the internal length is on, the
deleted volume comes out at exactly one slice on all three meshes and the
increment sequence is identical, so it is not the failure zone and not
time discretisation either.  The residual is how the same physical slice
redistributes load internally when it is six tetrahedra rather than
fifty-four, and crack-band scaling makes no claim about that.

The C3D8 arm has no legacy counterpart to compare against, because the
legacy path refuses the family outright.  There the projection is not a
more objective option - it is the only one.

That row was worth re-measuring rather than trusting, and for a specific
reason: its three arms stop at step times of 0.954, 0.994 and 0.962, and the
`work()` they were first measured with integrated each curve to wherever it
ended, so an arm was credited with less work for stopping earlier.  Re-measured
with the windowed `work(c, umax)` the spread is **0.72 %** against the published
0.71 %, and the pre-peak spread is exactly 0.00 % with 2.74 % after the peak.
The legacy and C3D4-projection rows need no such check: all six of those arms
run the step to completion.

## Does the internal length give mesh objectivity

Two sweeps, and the division between them is the point.  `run_objectivity.sh`
refines the cross-section, which holds the band one slice thick and so asks
an internal length nothing at all.  `run_nlobjectivity.sh` fixes the
specimen at length 6 and refines **along** the axis, so the band sizes
itself and the width is the model's rather than the deck's.

```
CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_nlobjectivity.sh /tmp/nl
CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_nlwidth_scaling.sh /tmp/sc
```

Three conditions have to hold before either is worth reading, and all three
are checked rather than assumed.

* **One specimen on every mesh.**  The narrowing is imposed on nodes, so a
  step notch is sampled differently by each mesh and sharpens as the mesh is
  refined.  The preflight interpolates the coarse deck's surface and compares
  it against every node of the finer ones; it fails the old configuration by
  `5.8e-2` of the cross-section and passes the tent by zero.
* **The same stage of failure on every arm.**  Not approximately: every arm
  runs to complete severance, `Dmax = 1`, elements deleted, `F = 0` well
  before the end of the step.  That is a stage nothing can be halfway into.
* **The stabiliser is not the answer.**  Viscosity is needed - with it off,
  four of six arms stop before rupture - so it is applied equally to every
  arm and its own effect is measured rather than assumed away.  Tripling it
  to `3.e-3` moves each arm's total work by 1.0 to 6.1 % and *reduces* both
  spreads, by 0.9 and 2.5 points: the disagreement below is not the
  stabiliser.  The width measure is less indifferent to it - the coarsest
  nonlocal arm moves from 1.000 to 1.333 - so that arm is quoted but not
  leaned on, and the verdict is the same at both viscosities.

Measured on `4d06196`, `ell = 0.5`, viscosity `1.e-3`, width by
`bandwidth.py`:

| quantity | local | `NONLOCAL=0.5` |
|---|---|---|
| band width, `h` = 1.0 / 0.5 / 0.25 | 1.049 -> 0.619 -> 0.592 | 1.407 -> 1.066 -> 1.261 |
| nonlocal / local | | 1.34 -> 1.72 -> 2.13 |
| peak force | 373.8 / 368.2 / 372.1 | 374.9 / 370.6 / 373.3 |
| work before the peak | spread 1.07 % | spread 1.03 % |
| work after the peak | spread 145 % | spread 124 % |

**A RETRACTED ROW BELONGS HERE.**  The first version of this table measured
the width as the count of integration points over `D>0.5` and reported
`1.000 -> 0.667 -> 0.417` against `1.000 -> 0.833 -> 1.083`, concluding that
the local band collapses while the nonlocal one holds.  `de1stats` writes
three thresholds, so testing that choice cost nothing and was not done.  On
the same runs `D>0.1` says the local band holds too and `D>0.9` says the
nonlocal band collapses too.  The verdict was a property of the number 0.5.

With the threshold-free measure the local band does **not** collapse: it falls
and nearly stops, and a sharper notch on finer meshes (span 1.0, `h` = 0.5 /
0.25 / 0.125) leaves the LOCAL arm at 0.634 / 0.545 / 0.642 - so the floor is
not the notch either.  Only the local arm: that sweep's nonlocal arm at
`h`=0.125 stopped on an unconverged attempt, so its width and the ratio built
on it are not comparable and are withheld.  The refutation rests on the local
row, which is where it always rested.  What holds it there is **the specimen's cross-section**, and
`run_xsection_floor.sh` measures that rather than arguing it (see below).  What survives is the
ratio, which needs no absolute scale and rises monotonically at both
viscosities (1.34 -> 1.72 -> 2.13 at `1.e-3`, 1.39 -> 1.83 -> 2.08 at
`3.e-3`).  The dissipation disagrees by more than a factor of two either way,
and the regularisation improves it by 20 points out of 145, which is nothing.

The refinement sweep cannot say why, because every candidate cause moves
with `h`.  `run_nlwidth_scaling.sh` holds `h` = 0.25 and moves `ell`:

| `ell` | width | width / `2*ell` | layers | `W_pre` | `W_post` |
|---|---|---|---|---|---|
| 0 (local) | 0.5918 | - | 2.37 | 6.3839 | 8.8725 |
| 0.25 | 0.8267 | 1.653 | 3.31 | 6.3851 | 11.2615 |
| 0.5 | 1.2608 | 1.261 | 5.04 | 6.3859 | 16.5533 |
| 1.0 | 2.1690 | 1.084 | 8.68 | 6.3860 | 29.9855 |

`W_pre` is identical to 0.03 %, so nothing outside the fracture process
moves.  The width follows `ell` and the dissipated work follows the width.

That is predictable with nothing fitted.  An element whose `charlen` is its
own size `h` reaches `D = 1` when `h*eps_p = u_f`, so it dissipates
`G_f = sigma_0 u_f / 2` per unit area - once, whatever its size, which is the
property Bazant and Oh 1983 build the crack band for.  A band `w` wide holds
`w/h` such layers, so

```
    W_post  ~  G_f * w / charlen        charlen = h  ->  G_f * w/h
                                        charlen = w  ->  G_f, always
```

Against `G_f = 4.08` read from the deck's own header, the ratio
`W_post / (G_f * w/h)` comes out at 0.92, 0.84, 0.80, 0.85 - nearly constant,
which is the statement that matters: **the dissipation is proportional to the
number of element layers in the band, not to the band's width.**

An earlier version of this section claimed the prediction held to `-8.0 /
-6.4 / +5.0 %` with nothing fitted; that accuracy belonged to the `D>0.5`
threshold and is withdrawn.

**The coefficient is 1, and the way to see that is to stop fitting.**  Fitting
`W_post = A + k G_f (w/h)` across the `ell` arms gave slopes of `0.83 G_f` at
`h`=0.25 and `1.39 G_f` at `h`=0.5, and that discrepancy was reported here as
the coefficient failing to transfer between meshes.  It is a two-parameter fit
to three or four scattered points: the free intercept absorbs whatever the
slope does not, so neither number measures `k`.

Measure `k` directly instead, one arm at a time, on the **local** arm - where
the substitution is inert and no internal length enters, so it is a statement
about the crack-band law alone - across the three meshes of the sharp-notch
sweep:

| `h` | `w` | `w/h` | `W_post` | `k = W_post/(G_f w/h)` |
|---|---|---|---|---|
| 0.500 | 0.6337 | 1.27 | 4.158 | 0.804 |
| 0.250 | 0.5454 | 2.18 | 10.699 | 1.202 |
| 0.125 | 0.6417 | 5.13 | 21.983 | 1.049 |

Mean 1.02, scatter 39 % of the mean, **no trend** across a factor of four in
element size.  So the crack-band law charges one `G_f` per element layer and
the layer count is `w/h`, with the coefficient 1 to within the scatter.  The
0.83-against-1.39 reading is superseded by this rather than explained by it.

**So the gradient backend is not failing objectivity - it is doing its job,
and the softening law is not.**  `calcdamage.f:1973` replaces the local
plastic-strain increment with the nonlocal average; `calcdamage.f:1983` still
multiplies it by `charlen`, the element's own width.  The averaging was given
a length and the softening law was not, so the fracture energy of the
combination is `G_f * w/charlen` - a material constant only when the band is
one element wide, which is exactly the case the crack band was derived for
and exactly the case an internal length abolishes.  Jirasek and Bauer 2012
section 5 state the requirement directly: the width entering the softening
law must be the width of the band that actually forms.

## What stops the local band collapsing

```
CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_xsection_floor.sh /tmp/xs
```

`--h` is the bar's transverse dimension, and it is one variable here: the
projected width of every Kuhn tetrahedron is the slice thickness along the
axis, so `charlen` does not move when the cross-section does, and the notch is
a fixed FRACTION of the section, so the relative stress profile along the bar
is identical too.

| `hcross` | `w` at `h`=0.5 | `/hcross` | `w` at `h`=0.25 | `/hcross` |
|---|---|---|---|---|
| 1.0 | 0.6186 | 0.619 | 0.5918 | 0.592 |
| 0.5 | 0.5672 | 1.134 | 0.3239 | 0.648 |
| 0.25 | 0.6315 | 2.526 | 0.2906 | 1.162 |

**Two floors, and they cross over.**  A band cannot be shorter than the element
carrying it, so at `h`=0.5 that floor hides everything else and `w` barely
moves with the cross-section - 0.619, 0.567, 0.632 across a factor of four.  At
`h`=0.25 the cross-section is free to act and halving it *halves* `w`: 0.5918
-> 0.3239, a factor of 0.547 for a factor of 0.5.  Quarter it and the element
floor takes over again, 0.2906 being 1.16 elements.

That decides between the two candidates rather than fitting either: no
dependence on the cross-section would have predicted 0.592 against 1.18, and
proportionality predicts 0.592 against 0.648.  The floor is the **specimen**.

The mechanism had two candidates that both predict a floor of order the
transverse dimension - Saint-Venant, which is geometric, and necking
triaxiality, which needs the material to contract laterally - and they differ
on Poisson's ratio.  **That test has now been run**, `NU=0.3` against the decks'
`nu=0`, everything else fixed, at `h`=0.25:

| `hcross` | `w` at `nu`=0 | `w` at `nu`=0.3 |
|---|---|---|
| 1.0 | 0.5918 | 0.6213 |
| 0.5 | 0.3239 | 0.3226 |
| `w/hcross` | 0.592 / 0.648 | 0.621 / 0.645 |

The floor and its proportionality survive unchanged, within 5 %, so the floor
is **not** the lateral contraction.  That eliminates triaxiality and leaves the
geometric reading standing.  It is still not a positive measurement of a
Saint-Venant decay length - nothing here measures one - but the alternative is
gone rather than merely unfavoured.

The consequence for everything above is the one worth carrying: on this
specimen the crack band is **never** one element wide, at any mesh, on any arm.
The premise the crack-band scaling is derived from does not hold here, which is
why the dissipation tracks the layer count in the local model too and not only
in combination with an internal length.

## The substitution, and what it did

`CCX_DAMAGE_NLWIDTH=1` makes the width in the softening law the width the
averaging forms, `2*ell`, wherever the mesh resolves it, and leaves the
element's own width elsewhere.  Off by default, so every earlier answer is
unchanged to the bit and the gate does not move.

Two guards, and both make it a no-op rather than a guess.  The length comes
from `damnlelleff(iel,ell,iok)` - the length that element was actually
averaged with, including the material factor and the localising `g(D)`, not a
global maximum - and it is used only at `iok=1`, by the criterion that lives
in `damnonlocal.f` and is deliberately not repeated in `calcdamage.f`, so the
two copies cannot drift apart.  And it never goes below the element's own
width, since a band cannot be narrower than the element carrying it.

Measured on the same sweep, `h` = 0.25, the same binary, the decks unchanged:

| | `NLWIDTH=0` | `NLWIDTH=1` |
|---|---|---|
| `W_post` at `ell` = 0.25 / 0.5 / 1.0 | 11.26 / 16.55 / 29.99 | 5.69 / 4.42 / 6.12 |
| `ell` dependence of `W_post` | **x2.66** | **x1.38** |
| `W_post / (G_f * w/h)` | 0.84 / 0.80 / 0.85 | 0.48 / 0.25 / 0.24 |
| local arm, unchanged by design | 8.87 | 8.87 |

Read the third row first.  Off, that ratio is constant - the dissipation is
one `G_f` per element layer, which is the defect.  On, it is no longer
constant, so the per-layer charge is gone, and the `ell` dependence falls
from 2.66 to 1.38.

### The gradient backend, asked the same question

`MODE=GRADIENT` runs the same sweep through the PDE form, which could not be
asked before `4d1c250` because a card-only length never reached its assembly.
Two things came out, and the second one is about this directory's own tooling.

| `ell/h` | `w / 2*ell`, integral | `w / 2*ell`, gradient |
|---|---|---|
| 1 | 1.65 | - |
| 2 | 1.26 | 1.20 |
| 4 | 1.08 | 0.83 |

So the band reaches `2*ell` at an `ell/h` of roughly 2 to 4 on both backends,
and the participation criterion that admits `ell >= 0.39 h` is more permissive
than that.  **Two comparable arms per backend is thin**, and the gradient rows
are consistent with the integral ones rather than an independent confirmation
of them.

The tooling point: the first gradient sweep came out non-monotone in `ell` -
0.882, 0.741, 1.582 - and the middle arm had stopped at `theta=6U` with nothing
deleted.  A width read off a run that never ruptured is the width of a
half-formed band, and nothing in the number says so.  `work()` has a common
window for exactly this, but **a width has no window**, so the script now
checks two independent signs of the same stage - the step completed AND
something was deleted - and marks any arm that fails them, excluding it from
the `ell` dependence.  Raising the viscosity does not fix the gradient arms
either: at `3.e-3` a *different* arm fails, which is worth knowing before
anyone treats viscosity as the way to make that backend comparable.

Measured again on a second mesh, `h` = 0.5, through `run_nlwidth_gate.sh`
end to end: **3.34 -> 1.35**.  The residual is therefore the same on both
meshes, 1.35 and 1.38, which says it is not a discretisation artefact - a
mesh-dependent leftover would not land twice on the same number.

**The residual is not diffuse - it is one condition, and it is a condition on
the specimen.**  With the switch on, the law's width is `2*ell`, so the
prediction to test is `G_f w/(2 ell)` and not `G_f w/h`.  At `h`=0.25, with
`ell`=0.375 excluded for not converging:

| `ell` | `w` | `w/(2 ell)` | predicted | measured | off |
|---|---|---|---|---|---|
| ~~0.25~~ | ~~0.7280~~ | ~~1.456~~ | ~~5.940~~ | ~~5.710~~ | ~~-4 %~~ WITHDRAWN |
| 0.5 | 1.0725 | 1.072 | 4.376 | 4.442 | +1 % |
| 0.75 | 1.3193 | 0.880 | 3.588 | 5.093 | +42 % |
| 1.0 | 1.5352 | 0.768 | 3.132 | 6.135 | +95 % |

**The first row is withdrawn**, and by the check in this directory rather than
by anything external: `ell`=0.25 stops at a step time of 0.147 with two elements
deleted, so it never ruptured.  It was published as the strongest point of the
table an hour before the stage check was rewritten to read the step time
instead of the attempt counter - the corrected check then caught a number that
the broken one had passed.  Its exclusion leaves **one** arm inside the domain,
not two.

The error is **monotone in `w/(2 ell)`** across the three arms that do reach
rupture, and crosses zero where the band stops being as wide as the length the
law charges.  At the one valid arm with `w >= 2 ell` the substitution is right
to 1 %; outside the domain the law asserts a width the specimen does not
sustain, and the accounting fails by the overshoot.  One point inside a domain
is a consistency, not a verification, and the honest way to strengthen it is an
arm at `w/(2 ell)` well above 1 that also ruptures - which on this specimen
means a smaller `ell`, and `ell`=0.375 already fails to converge.

So the substitution's validity condition is **not** `iok` - the mesh resolves
`ell` on all four of those arms - but `w >= 2 ell`, a statement about the
specimen.  On this bar the effective width inferred from the measured energy
saturates near 1.0, which is the cross-section, and the cross-section is what
sets the floor too; that reading is inside the coefficient's own 20 % scatter,
so it is named and not established.

What remains genuinely open is smaller than it was: `W_post` still varies with
`ell` inside the valid domain, because `w/(2 ell)` is 1.46 and 1.07 rather than
1 - the law is satisfied and the band the averaging forms simply is not exactly
`2*ell`.  And `ell`=0.375 failed to converge at a step time of 0.13 between two
values that both converge, which is not explained.  `run_nlwidth_gate.sh`
therefore judges by comparison - the same sweep with the switch off and on,
and the dependence must come down - because a fixed ceiling would have to be
picked above whatever today's value happens to be, and picking a threshold to
pass is the mistake the retracted width measure above already made once.

## The unit test

```
test/crackband/run_cbwidth_test.sh
```

Eight geometries whose answer can be worked out by hand: a unit cube pulled
along an axis and along a face diagonal, an elongated brick pulled along its
long side and then across its thin one, a corner tetrahedron, and three
degenerate cases where no band normal exists - hydrostatic loading, where
every eigenvalue is equal and no direction is the major one.  One case is
a C3D10 whose midside node has been pushed far outside the element, which
checks two things at once: that the width is unmoved, which it can only be
if the midside nodes are really excluded, and that the higher-order warning
fires.  A midside node left on its straight edge would hide that bug, since
it lies inside the corner hull and changes no extent.

Each case prints a verdict and a disagreement sets the exit code.

**The gate runs it**, as a preflight, and this paragraph used to say it was a
live proposal not yet landed - which was true when written and stopped being
true without the text noticing.  Check rather than trust either statement:

```
grep -n crackband test/regress/run.py
```

A hit means a wrong width turns the gate red by name.  No hit would mean the
opposite and is worth knowing: no regression case pins the geometry directly,
so a width error would only be caught once it moved a trajectory, arriving as
a mystery instead of a named failure.

Two things about how it got there are worth keeping.  The gate's own comment
records that each preflight self test was shown to go red **by the agent who
did not write it** - and that one of the two was proposed for the preflight
while it still could not fail at all, its author retracting the claim.  That
author was me: my self test reported disagreement and returned zero, so it
would have sat in the gate as decoration.  A test nobody has watched fail is
a claim, not a test.

## The energy equivalence

```
CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_energy_equiv.sh /tmp/eq
```

A `G_f` card and the equivalent `u_f` card must describe the same material,
since `G_f = sigma_0 u_f / 2` for the linear law.  The script runs both at
two strain levels a factor of ten apart, **bounds both**, and decides.  It
does not print a number for a reader to interpret, and it does not name a
cause.

Measured on `37d0d53`: 0.0822 % at working strain against a 0.5 % ceiling,
0.0082 % at one tenth of it against 0.05 %.  Remove the peak rule in
`damcbufset` and both go over - 2.9774 % and 0.1707 % - and the script
exits 1.

**The ratio is printed and does not decide, and that is the part worth
reading.**  It used to decide.  The argument was that an exact conversion
leaves only a residual that shrinks with the strain, so a big ratio meant
the conversion was sound.  That argument named a cause - a finite-strain
stress-measure effect - and the cause was wrong: the residual was
`sigma_0` sampled after the point had begun to unload, which scales with
strain in exactly the same way.  The ratio sat at 17.44 against a
threshold of 5 for as long as the defect existed.

Two things follow, and both cost something to learn here:

* **a check is blind to whatever its criterion is not a function of.**
  The ratio was a function of how the disagreement scales, so an error
  that scaled the same way stayed invisible however large it grew.
* **a retracted diagnosis leaves working code behind.**  When the cause
  was retracted, the ratio branch kept its power to exit green *and kept
  asserting the retracted cause while doing so* - so the one path
  reachable when something had gone wrong would have explained it away.
  Agent 1 found that in review; the fix was to take the verdict off it
  entirely rather than to reword it.
