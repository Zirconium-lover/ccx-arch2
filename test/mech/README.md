# test/mech - mechanics tests with a reference

    CCX_EXE=$PWD/src/ccx_2.23_pardiso python3 test/mech/run.py [name ...] [--work DIR] [--list]

One directory per mechanism.  `<name>/test.py` generates its decks, runs the
solver and returns checks; `run.py` prints every check (measured, reference,
tolerance, why) and a pass/fail table, and exits with the number of failed
tests.  Every test runs in under 60 s on one thread.  The gate runs the whole set as the case `mech` in test/regress/cases.json, except a test that is red by design (output_invariance, until its fix).

Statuses: PASS / FAIL; XFAIL is a known limitation of the code that is still
failing as documented (counts as passing), XPASS means it stopped failing and
must be revisited (counts as failing); INFO is a measured number with no
criterion yet.

`MECH_RED_ENV="NAME=VALUE ..."` is applied over every solver call, so a test
can be shown red on a wrong variant without editing it.  Each test's
docstring records its red path: what was run and which checks went red.

| test | mechanism | reference |
|---|---|---|
| `nl_objectivity` | nonlocal GRADIENT/INTEGRAL: W and band width converge in h at fixed ell; ell->0 is the unaveraged scheme; h > ell (XFAIL) | convergence between meshes; local model as red path |
| `lag_order` | lagged nonlocal scheme: FROZEN -> implicit LOCAL as dt -> 0, order; INNER=1 | self-convergence against the implicit solution |
| `uc6_law` | UC6 facet driven through the solver: T(d) modes I/II/45 deg, work = Gc, secant unloading, compression Kn, shear under compression, gmin residual, Duvaut-Lions viscosity; rigid rotation (XFAIL) | closed-form bilinear law, Duvaut-Lions ODE |
| `bulk_law` | bulk DE1/DM2 law on one element: omega against eps_f(eta) increment by increment (Rice-Tracey and table; uniaxial, equibiaxial, pure shear, mixed), dD/d eps_p = L/u_f with CHARLEN 0/1 on elements turned 45 deg, W = (V/L) G_f, EVOLUTION=ENERGY; eta and sigma_0 from PK2 (XFAIL) | the law as coded, Rice & Tracey 1969, Jirasek & Bauer 2012, Bazant & Oh 1983 |
| `damper` | CCX_DAMAGE_VISCOUS_DAMPING: dashpot force law through the external work of an elastic bar (exact), the energy audit (two flaws, XFAIL), W -> static as c -> 0 | closed form for a linear field; self-convergence |
| `output_invariance` | .sta and .damage byte-identical at output FREQUENCY 1 and 1000; bulk LOCAL/FROZEN/GRADIENT, viscosity on/off, UC6. **Red by design until the output results() stops rewriting the committed damage** (developer's acceptance test) | none needed: a correct solver does not read the output frequency |
| `nl_benchmark` | nonlocal model against the exact 1-D bar B1 of the library: W/G_f = 1.98 (C3D8), 1.15 (C3D8I) against 1.008 - the plastic flow spreads past the weak element (lateral constraint, C3D8 locking); XFAIL with the mechanism, anatomy pinned | library notes/nonlocal-benchmarks.md sec. 1; Bazant & Jirasek 2002 eqs. 61-63 |
