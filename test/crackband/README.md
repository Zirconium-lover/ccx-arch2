# Crack-band width: what is here and how to check it

Six files, no obvious entry point, so this is the map.  Every number below
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

## The unit test

```
test/crackband/run_cbwidth_test.sh
```

Eight geometries whose answer can be worked out by hand: a unit cube pulled
along an axis and along a face diagonal, an elongated brick pulled along its
long side and then across its thin one, a corner tetrahedron, and three
degenerate cases where no band normal exists - hydrostatic loading, where
every eigenvalue is equal and no direction is the major one.  It returns
non-zero on disagreement, checked by changing an expected value.

Whether the gate runs it is worth checking rather than assuming, since it is
a live proposal and not mine to land:

```
grep -n crackband test/regress/run.py
```

No hit means a wrong width would not turn the gate red on its own: no
regression case pins the geometry directly, so only a width error big enough
to move a trajectory gets caught, and then it arrives as a mystery rather
than as a named failure.
