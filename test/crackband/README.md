# Crack-band width: what is here and how to check it

Eight files, no obvious entry point, so this is the map.  Every number below
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
| C3D8 projection | 1.000 on every mesh | 0.71 % | law predicts zero |

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

Measured on `4d06196`, `ell = 0.5`, viscosity `1.e-3`:

| quantity | local | `NONLOCAL=0.5` |
|---|---|---|
| band width, `h` = 1.0 / 0.5 / 0.25 | 1.000 -> 0.667 -> 0.417 | 1.000 -> 0.833 -> 1.083 |
| the same in element layers | 1.00 -> 1.33 -> 1.67 | 1.00 -> 1.67 -> 4.33 |
| peak force | 373.8 / 368.2 / 372.1 | 374.9 / 370.6 / 373.3 |
| work before the peak | spread 1.07 % | spread 1.03 % |
| work after the peak | spread 145 % | spread 124 % |

**The width is held and the energy is not, and those are two different
verdicts about the same run.**  The local band is one to two elements wide on
every mesh - a width proportional to `h`, which is the pathology - while the
nonlocal band spans 1 to 4.3 elements at a length that stays near `2*ell`.
So the averaging does set the width.  The dissipation meanwhile disagrees by
more than a factor of two, and the regularisation improves it by 20 points
out of 145, which is nothing.

The refinement sweep cannot say why, because every candidate cause moves
with `h`.  `run_nlwidth_scaling.sh` holds `h` = 0.25 and moves `ell`:

| `ell` | width | width / `2*ell` | layers | `W_pre` | `W_post` |
|---|---|---|---|---|---|
| 0 (local) | 0.4167 | - | 1.67 | 6.3839 | 8.8725 |
| 0.25 | 0.7500 | 1.500 | 3.00 | 6.3851 | 11.2615 |
| 0.5 | 1.0833 | 1.083 | 4.33 | 6.3859 | 16.5533 |
| 1.0 | 1.7500 | 0.875 | 7.00 | 6.3860 | 29.9855 |

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

Against `G_f = 4.08` read from the deck's own header: `-8.0 %`, `-6.4 %`,
`+5.0 %` on the three nonlocal arms, and `+30.5 %` on the local one, where
the band is 1.67 elements and a count of points over a threshold is at its
coarsest.

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

The substitution that removes it is one line in the width cache, and it is
not taken here: it changes the answer of every existing nonlocal run, the
length has to arrive from the nonlocal side, and which length it should be -
`2*ell`, the measured width, or `u_f` scaled by `h/w` instead - is a
modelling decision and not a local one.  The numbers above are what such a
change has to be judged against, and `run_nlwidth_scaling.sh` is the shape
of the gate: after the fix `W_post` must stop depending on `ell`, which it
currently does by a factor of 2.66.

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

Whether the gate runs it is worth checking rather than assuming, since it is
a live proposal and not mine to land:

```
grep -n crackband test/regress/run.py
```

No hit means a wrong width would not turn the gate red on its own: no
regression case pins the geometry directly, so only a width error big enough
to move a trajectory gets caught, and then it arrives as a mystery rather
than as a named failure.

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
