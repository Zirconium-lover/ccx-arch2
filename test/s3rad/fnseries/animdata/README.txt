Packed results of the Fn series for animations and plots
=========================================================

One .npz per run of the final table (results_base/SERIES_FINAL.txt), made by
tools/extract_anim.py from the raw m.inp / m.frd / m.dat / m.damage /
m.de1.vtk / m.de1stats (0.7 GB per run, not kept).  The key list is in the
docstring of extract_anim.py; anim/fndata.py loads them:

    import sys; sys.path.insert(0, 'test/s3rad/fnseries/anim')
    from fndata import Run, SERIES
    r = Run(0.4); s = r.at(0.05)      # state at the grip strain 5 %

    run          nominal Fn   radial plates k   k/12
    fn0          0            0                 0
    fn20         0.20         2                 0.167
    fn25v1e-3    0.25         3                 0.25
    fn33v1e-3    0.33         4                 0.333
    fn40         0.40         5                 0.417
    fn60         0.60         7                 0.583
    fn80v1e-3    0.80         10                0.833
    fn100        1            12                1

The decks place 12 equal plates, so the radial fraction of a run is k/12 (the
nominal Fn was rounded to whole plates when the decks were made).

Fields are at the output frames (every 10 increments, ~0.5 % of grip
strain), in float16; deletion times of elements are exact.  omega and D are
rebuilt from the PEEQ and stress history with the law of calcdamage.f
(Run.check_omega(): against the solver's own values at the last state the
median ratio is 0.96-1.01, correlation 0.73-0.80, omega >= 1 agrees for
81-88 % of the tets): good for pictures, not for numbers.

data_preview.png (anim/data_preview.py): curves, mechanisms at the end, 3D
crack with the plates, sections of omega / PEEQ / sigma_xx.
