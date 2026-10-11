"""Figures of the experiment slide, cut from the two papers and made to sit on
the dark slides.  The pictures are third-party: they go to SCRATCH, never to
the repository.

    python article_figs.py saburov2025.pdf petrov2025.pdf out_dir

Saburov et al., YaFI 2025, 16 (4) 442-454:
    fig. 5 (page 6, embedded raster): fractographs, (a) ductile, (b) brittle.
        Cut square, inside the photo edges, above the scale-bar box (the bar
        is drawn on the slide: 10 um = SCALE_PX pixels of the raster).
    fig. 9b (page 8, rendered: its text is vector): elongation vs Fn.
Petrov et al., Letters on Materials 2025, 15 (1) 49-54:
    fig. 6 (page 4, embedded raster with its soft mask): brittle share vs Fn;
        the legend is translated.
The two charts are inverted: light ink on transparent.
"""
import os, subprocess, sys, tempfile
import numpy as np
from PIL import Image, ImageDraw, ImageFont

SCALE_PX = 76                       # 10 um in the fig. 5 raster (bar measured on a page render)
INK = (0xD8, 0xDD, 0xE3)
SERIF = '/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf'


def images(pdf, page, d):
    subprocess.run(['pdfimages', '-png', '-f', str(page), '-l', str(page), pdf, os.path.join(d, 'p')], check=True)
    return sorted(os.path.join(d, f) for f in os.listdir(d) if f.startswith('p-'))


def inverted(grey, out):
    """Black ink on white -> INK on transparent."""
    a = 255 - np.asarray(grey, float)
    rgba = np.zeros(a.shape + (4,), np.uint8)
    rgba[..., :3] = INK
    rgba[..., 3] = a.clip(0, 255).astype(np.uint8)
    Image.fromarray(rgba).save(out)


def main(sab, pet, out):
    os.makedirs(out, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        # fractographs: two photos side by side in one raster, white between
        im = Image.open(images(sab, 6, d)[0]).convert('L')
        a = np.asarray(im)
        on = (a < 235).mean(0) > 0.5
        edges = np.flatnonzero(np.diff(on.astype(int)))          # [l0, r0, l1, r1]
        top = np.flatnonzero((a < 235).mean(1) > 0.3)[0]
        side = 585                                               # square, clear of the scale-bar box
        for name, left in (('fracto_ductile.png', edges[0]), ('fracto_brittle.png', edges[2])):
            x0, y0 = left + 5, top + 4
            im.crop((x0, y0, x0 + side, y0 + side)).save(os.path.join(out, name))
    with tempfile.TemporaryDirectory() as d:
        # fig. 9b: page rendered at 600 dpi (axis text is vector), panel cut, its letter removed
        subprocess.run(['pdftoppm', '-r', '600', '-f', '8', '-l', '8', '-png', '-singlefile', sab,
                        os.path.join(d, 'pg')], check=True)
        g = Image.open(os.path.join(d, 'pg.png')).convert('L')
        panel = g.crop((2505, 4235, 4295, 5840))
        ImageDraw.Draw(panel).rectangle((3343 - 2505 - 6, 0, 3436 - 2505 + 6, 4319 - 4235 + 6), fill=255)
        inverted(panel, os.path.join(out, 'elong_fn.png'))
    with tempfile.TemporaryDirectory() as d:
        # fig. 6: raster + soft mask on white; "Manual annotation" -> "ручная разметка"
        p = images(pet, 4, d)
        rgb = np.asarray(Image.open(p[2]).convert('L'), float)
        al = np.asarray(Image.open(p[3]).convert('L'), float) / 255
        g = Image.fromarray((rgb * al + 255 * (1 - al)).clip(0, 255).astype(np.uint8))
        dr = ImageDraw.Draw(g)
        font = ImageFont.truetype(SERIF, 24)
        for base in (71, 143):                                   # baselines of the two entries
            dr.rectangle((214, base - 19, 405, base + 6), fill=255)
            dr.text((217, base), 'ручная разметка', font=font, fill=0, anchor='ls')
        inverted(g, os.path.join(out, 'brittle_fn.png'))
    for f in sorted(os.listdir(out)):
        print(f, Image.open(os.path.join(out, f)).size)


if __name__ == '__main__':
    main(*sys.argv[1:4])
