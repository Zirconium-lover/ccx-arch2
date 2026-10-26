# test/mech - mechanics tests with a reference

    CCX_EXE=$PWD/src/ccx_2.23_pardiso python3 test/mech/run.py [name ...] [--work DIR] [--list]

One directory per mechanism.  `<name>/test.py` generates its decks, runs the
solver and returns checks; `run.py` prints every check (measured, reference,
tolerance, why) and a pass/fail table, and exits with the number of failed
tests.  Every test runs in under 60 s on one thread.

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
