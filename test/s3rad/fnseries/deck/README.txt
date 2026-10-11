The talk "Роль ориентации гидридов в разрушении циркониевых сплавов по данным
конечно-элементного моделирования" as a PPTX (16:9, dark, Cambria).

    python formulas.py assets                       # LaTeX formulas -> white PNG
    python ../anim/fracture_surface.py <in>/fs      # fracture surfaces, slide 11
    ANIM=../anim SCRATCH=<in> python build_deck.py out.pptx

SCRATCH holds what is not kept in this public repository: fn_measure.gif and
fn_block.gif renders, the article figures of slide 7 (fracto_rings.png,
brittle_fn.png, elong_fn.png, cropped from Saburov et al., YaFI 2025, and
Petrov et al., Letters on Materials 2025) and fs/fs_<run>.png.  The deck
itself is not committed for the same reason.

Slides 3-6 carry one GIF each with its text; no bullets anywhere (plain
paragraphs), the conclusions are numbered.  Slides 9 and 12 use native
charts, each headed above.  The damage law on slide 9 is written as in the
literature (ABAQUS ductile damage, Kim and Lee 2024): onset by omega, then
D = u/u_f with du = L d(eps_pl), which is what the solver does.

Size: the GIFs of slides 3-6 are re-encoded for the deck at 1344 px
(they fill 70 % of the slide width; texture and fn_measure also at 20 fps,
64 / 112 colours), GIFS=<dir> points the builder at them; the deck is then
~28 MB.
