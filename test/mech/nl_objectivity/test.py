"""Nonlocal regularisation: mesh objectivity at a fixed internal length.

The bar is refined ALONG its axis with ell fixed, h = ell/2, ell/4, ell/8,
on both backends.  What a regularised model owes is a dissipated energy and a
band width that stop depending on h; the local model owes nothing and is run
on the same meshes as the red path of the criterion - it must FAIL the same
spread bound, or the bound is not able to tell the two apart.

Measured on 71ed70d (ell=0.5; h = 0.25 / 0.125 / 0.0625):
  W_total  GRADIENT 13.99 / 14.11 / 14.23   INTEGRAL 11.52 / 11.74 / 11.81
           LOCAL    15.29 / 26.92 / 37.03   (spread 147 %)
  width    GRADIENT 1.230 / 1.271 / 1.205   INTEGRAL 1.041 / 1.096 / 1.077
  h = 2 ell (nslice 6): GRADIENT 19.59 (+38 %), INTEGRAL 13.22 (+12 %)
  W_post/(G_f A): GRADIENT 2.0, INTEGRAL 1.6, local one-element band 1.00-1.05

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_NLWIDTH=0 (the element width in the
softening law instead of l_d) turns six checks red: W_total spread 44 %
(GRADIENT) and 49 % (INTEGRAL), width spread 15 % / 23 %, and the h=ell/8
arms of both backends stop rc=201.

Tolerances.  5 % on W_total: the converged arms spread 1.7 % (GRADIENT) and
2.5 % (INTEGRAL) and the local arm 147 %, so 5 % sits a factor two above the
noise of a regularised model and thirty below the failure it must catch.  The
width carries its own noise (it is read from the damage field of the last
snapshot, and the separation end stops each arm at a slightly different
stage), measured at 5.3 %, so it gets 10 %.
"""
import os
from mechlib import Check, cross_bar, spread, BAR_GF, BAR_AREA

WHAT = "nonlocal (GRADIENT, INTEGRAL): W and band width converge in h at fixed ell; ell->0 is the unaveraged scheme"
SOURCE = ("convergence between meshes (Peerlings et al. 1996 sec. 5 bar); "
          "G_F = G_f via l_d: Bazant & Jirasek 2002 eqs. 61-63 (recorded, not asserted)")

ELL = 0.5
MESHES = (24, 48, 96)          # h = ell/2, ell/4, ell/8 on a bar of length 6


def run(work, exe):
    ck = []
    res = {}
    for mode in ("GRADIENT", "INTEGRAL", "LOCAL"):
        res[mode] = [cross_bar(os.path.join(work, "%s_%d" % (mode, n)), n,
                               None if mode == "LOCAL" else ELL, mode, exe=exe)
                     for n in MESHES]
    for mode in ("GRADIENT", "INTEGRAL"):
        r = res[mode]
        ck.append(Check("%s all arms finish (rc)" % mode,
                        max(abs(x["rc"]) for x in r), 0, tol=0))
        ck.append(Check("%s W_total spread over h=l/2..l/8" % mode,
                        spread([x["wtot"] for x in r]), 0.05, cmp="le"))
        ck.append(Check("%s band width spread" % mode,
                        spread([x["width"] for x in r]), 0.10, cmp="le"))
        ck.append(Check("%s W_post/(G_f*A) at h=l/8" % mode,
                        r[-1]["wpost"] / (BAR_GF * BAR_AREA), 1.0, cmp="info",
                        why="l_d theory predicts 1; the local one-element band gives 1.00-1.05"))
    ck.append(Check("RED PATH: LOCAL W_total spread must exceed 5 %",
                    spread([x["wtot"] for x in res["LOCAL"]]), 0.05, cmp="ge",
                    why="the same bound applied to the unregularised model"))
    # ell -> 0: the regularisation must switch itself off and leave exactly
    # the staggered scheme it runs in, which is FROZEN (the lag without the
    # average).  Its distance from LOCAL is the lag error, test lag_order.
    r0 = {}
    for mode in ("GRADIENT", "INTEGRAL", "FROZEN"):
        r0[mode] = cross_bar(os.path.join(work, "%s_ell0" % mode), 24, 0.01, mode, exe=exe)
    for mode in ("GRADIENT", "INTEGRAL"):
        ck.append(Check("%s ell=0.01: W_total equals FROZEN" % mode,
                        r0[mode]["wtot"], r0["FROZEN"]["wtot"], rtol=1e-9,
                        why="ell is 1/25 of h: every element is unresolved and gets its local value"))
    ck.append(Check("FROZEN vs LOCAL W_total at h=0.25 (lag error)",
                    r0["FROZEN"]["wtot"], res["LOCAL"][0]["wtot"], cmp="info",
                    why="measured +3.75 %; see test lag_order"))
    # h > ell: outside the resolved range.  Known limitation, documented.
    for mode in ("GRADIENT", "INTEGRAL"):
        rc = cross_bar(os.path.join(work, "%s_coarse" % mode), 6, ELL, mode, exe=exe)
        ck.append(Check("%s h=2*ell: W_total within 5 %% of converged" % mode,
                        rc["wtot"], res[mode][-1]["wtot"], rtol=0.05, expect_fail=True,
                        why="h > ell does not resolve the length; measured +38 % GRADIENT, +12 % INTEGRAL"))
    return ck
