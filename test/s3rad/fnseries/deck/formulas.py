"""LaTeX formulas as white-on-transparent PNG for the dark slides.

    python formulas.py out_dir
Each formula: pdflatex (standalone) -> pdftoppm at 600 dpi -> alpha from the
grey level, colour white (or the given colour)."""
import os, subprocess, sys, tempfile
import numpy as np
from PIL import Image

FORMULAS = {
    'omega': r'\omega=\int\limits_0^{\bar\varepsilon^{\,p}}\frac{d\bar\varepsilon^{\,p}}{\varepsilon_f(\eta)}',
    'eta': r'\eta=\frac{\sigma_m}{\sigma_{eq}}',
    'damage': r'D=\frac{L_e\,\Delta\bar\varepsilon^{\,p}}{u_f}',
    'sigma': r'\sigma=(1-D)\,\tilde\sigma',
    'gc': r'G_c=\int\limits_0^{\delta_f}T\,d\delta',
}
TEMPLATE = r'''\documentclass[border=2pt]{standalone}
\usepackage{amsmath}
\begin{document}
$\displaystyle %s$
\end{document}
'''


def render(name, tex, out, color=(255, 255, 255), dpi=600):
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, 'f.tex'), 'w').write(TEMPLATE % tex)
        subprocess.run(['pdflatex', '-interaction=nonstopmode', 'f.tex'], cwd=d, check=True,
                       stdout=subprocess.DEVNULL)
        subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', 'f.pdf', 'f'], cwd=d, check=True)
        g = np.asarray(Image.open(os.path.join(d, 'f.png')).convert('L'), float)
    a = (255 - g).clip(0, 255).astype(np.uint8)
    rgba = np.zeros(g.shape + (4,), np.uint8)
    rgba[..., :3] = color
    rgba[..., 3] = a
    Image.fromarray(rgba).save(os.path.join(out, 'f_%s.png' % name))
    return g.shape


if __name__ == '__main__':
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for k, v in FORMULAS.items():
        print(k, render(k, v, out))
