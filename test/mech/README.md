# test/mech - mechanics tests with a reference

    CCX_EXE=$PWD/src/ccx_2.23_pardiso python3 test/mech/run.py [name ...] [--work DIR] [--list]

One directory per mechanism.  `<name>/test.py` generates its decks, runs the
solver and returns checks; `run.py` prints every check (measured, reference,
tolerance, why) and a pass/fail table, and exits with the number of failed
tests.  Every test runs in under 60 s on one thread, except those marked LONG (`--long`, or named on the command line), which are not in the gate.  The gate runs the whole set as the case `mech` in test/regress/cases.json, except a test that is red by design (output_invariance, until its fix).

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
| `nl_benchmark` (LONG) | nonlocal model against the exact 1-D bar B1 of the library: W/G_f = 1.98 (C3D8), 1.15 (C3D8I) against 1.008 - the plastic flow spreads past the weak element (lateral constraint, C3D8 locking); XFAIL with the mechanism, anatomy pinned | library notes/nonlocal-benchmarks.md sec. 1; Bazant & Jirasek 2002 eqs. 61-63 |
| `neck_dissipation` | C3D4 bar vs the 1-D B1 reference: the local crack band's W/G_f grows as h is refined (2.11 / 8.46 at h = 0.125 for CHARLEN 0 / 1), the nonlocal one converges to ~1.35 G_f; known 3-D property, pinned + XFAIL | library notes/nonlocal-benchmarks.md sec. 1 |
| `scale_invariance` | x1 vs x1/25: bulk law and UC6 facet identical when the numerical parameters follow their groups; element-level size effect u_f/L x k; UC6 overshoot against mu v/delta0 (disk 0.067: T_max +5 %, cell 1.25: +70 %) | library notes/scaling-levers.md; exact Duvaut-Lions (cohesive_ref.py) |
| `scale_bar` (LONG) | fracturing C3D4 bar x1 vs x1/25: byte-identical when u_f, ell and Delta scale with the size; with the physics kept, strength unchanged and the size effect recorded (local u_f/h 0.4 -> 10; INTEGRAL ductility x80) | library notes/scaling-levers.md sec. 1, 3, 5 |
| `size_effect` (LONG) | notched UC6 plate, D/l_ch 0.1-10: net-section strength at small size (0.978), LEFM approached from below at large size (0.924, c_f = 0.52 l_ch), Bazant fit | net section; LEFM M(T) (Tada, Feddersen); Bazant 1984 |
| `uc6_failed_open` | ACCEPTANCE (not in the gate until merged): CCX_UC6_FAILED_OPEN=1 - floor removed (g = max(0, 1 - D)): T = 0 in tension and shear past delta_f, Kn in contact, W = Gc, the size_effect plate breaks at g_min = 1e-3 | closed-form bilinear law without the floor; size_effect |
| `uc6_bk` | ACCEPTANCE (not in the gate until merged): B-K mixed mode, 8 constants - W = G_c(B) at B 0..1, peaks Tn0/Ts0, reduction to the old law at G_IIc = G_Ic | Benzeggagh & Kenane 1996; Turon 2006/2010; library cohesive_ref.py |
| `bugfix_accept` | ACCEPTANCE for pending bug fixes: neighbour count scale-free at 2 ell = 8 h (XFAIL), deleting the last element (XFAIL, SIGSEGV), dashpot audit finite (PASS) | reproducers |
