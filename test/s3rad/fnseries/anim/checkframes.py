#!/usr/bin/env python3
"""checkframes.py <video.mp4> [n=8]: sample n frames, check that nothing is
drawn in the title band or the side/bottom margins (CONTENT of style.py),
and write <video>_sheet.png - the frames with the grid drawn on them."""
import subprocess, sys, tempfile, os
import numpy as np
from PIL import Image, ImageDraw

v, n = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 8
dur = float(subprocess.check_output(['ffprobe', '-v', 'error', '-show_entries', 'format=duration',
                                     '-of', 'csv=p=0', v]))
tmp = tempfile.mkdtemp()
frames = []
for i in range(n):
    t = dur * (i + 0.5) / n
    f = os.path.join(tmp, '%d.png' % i)
    subprocess.run(['ffmpeg', '-v', 'error', '-y', '-ss', '%.3f' % t, '-i', v, '-frames:v', '1', f], check=True)
    frames.append((t, Image.open(f).convert('RGB')))
W, H = frames[0][1].size
u = H / 8.0
def px(x, y): return (W / 2 + x * u, H / 2 - y * u)
X0, Y1 = px(-6.71, 2.85); X1, Y0 = px(6.71, -3.65)
bad = 0
for t, im in frames:
    a = np.asarray(im).max(axis=2) > 40
    ys, xs = np.nonzero(a)
    if not len(xs):
        continue
    l, r, top, b = xs.min(), xs.max(), ys.min(), ys.max()
    msg = []
    if top < Y1 - 2: msg.append('title band (top at %d px, band ends %d)' % (top, Y1))
    if l < X0 - 2: msg.append('left margin (%d < %d)' % (l, X0))
    if r > X1 + 2: msg.append('right margin (%d > %d)' % (r, X1))
    if b > Y0 + 2: msg.append('bottom margin (%d > %d)' % (b, Y0))
    print('t=%5.2f  ink x %4d..%4d  y %4d..%4d  %s' % (t, l, r, top, b, '; '.join(msg) or 'ok'))
    bad += bool(msg)
cols = 2
sheet = Image.new('RGB', (W // 2 * cols, H // 2 * ((n + cols - 1) // cols)))
for i, (t, im) in enumerate(frames):
    g = im.copy(); d = ImageDraw.Draw(g)
    d.rectangle([X0, Y1, X1, Y0], outline=(80, 80, 80), width=2)
    sheet.paste(g.resize((W // 2, H // 2)), ((i % cols) * W // 2, (i // cols) * H // 2))
out = os.path.splitext(v)[0] + '_sheet.png'
sheet.save(out)
print('FRAMES %d with problems; sheet %s' % (bad, out))
