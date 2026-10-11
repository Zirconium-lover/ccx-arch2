The talk "Роль ориентации гидридов в разрушении циркониевых сплавов по данным
конечно-элементного моделирования" as a PPTX (16:9, dark, Cambria).

    python formulas.py assets                       # LaTeX formulas -> white PNG
    python ../anim/fracture_surface.py <in>/fs      # fracture surfaces, slide 9
    ANIM=../anim SCRATCH=<in> python build_deck.py out.pptx

SCRATCH holds what is not kept in this public repository: fn_measure.gif and
fn_block.gif renders, the article figures of slide 5 (fracto_rings.png,
brittle_fn.png, elong_fn.png, cropped from Saburov et al., YaFI 2025, and
Petrov et al., Letters on Materials 2025) and fs/fs_<run>.png.  The deck
itself is not committed for the same reason.

Slides 3 and 4 swap one GIF (with its theses) for another on a click
(p:timing written by click_swap); slides 7 and 10 use native charts.
