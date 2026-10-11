"""The talk as a PPTX: "Роль ориентации гидридов в разрушении циркониевых
сплавов по данным конечно-элементного моделирования".

Dark slides (the GIFs are black), Cambria, 16:9.  The animations go in as
GIFs; slides 3 and 4 swap one GIF (with its theses) for another on a click.
Charts of slides 7 and 10 are native PowerPoint charts (editable, in the
deck's font); formulas are LaTeX pictures (formulas.py).

    python build_deck.py out.pptx

Inputs: the GIFs (ANIM), the formula and figure pictures (ASSETS), the packed
runs (fndata).  Third-party figures (micrograph of fn_measure, article
figures of slide 5) are read from SCRATCH and never go to the public repo.
"""
import copy, os, sys
import numpy as np
from lxml import etree
from pptx import Presentation
from pptx.chart.data import XyChartData, CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_MARKER_STYLE, XL_TICK_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt, Emu

HERE = os.path.dirname(os.path.abspath(__file__))
ANIM = os.environ.get('ANIM', os.path.join(HERE, '..', 'anim'))
ASSETS = os.environ.get('ASSETS', os.path.join(HERE, 'assets'))
SCRATCH = os.environ.get('SCRATCH', ASSETS)
GIFS = os.environ.get('GIFS', ANIM)            # the GIFs (a lighter re-encode may stand in)
sys.path.insert(0, ANIM)

FONT = 'Cambria'
W, H = 13.333, 7.5
BG = RGBColor(0, 0, 0)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GREY = RGBColor(0xB8, 0xBE, 0xC6)
DIM = RGBColor(0x7C, 0x84, 0x8E)
ACC = RGBColor(0xFF, 0xD2, 0x4A)
RAD = RGBColor(0xE8, 0x60, 0x4C)
TAN = RGBColor(0x2C, 0xB1, 0xD9)
MECH = {'ductile': RGBColor(0xA9, 0xB4, 0xC2), 'cleavage': RGBColor(0xFF, 0xA0, 0x40),
        'debond': RGBColor(0x3C, 0xC9, 0xA0)}
CARD = RGBColor(0x14, 0x19, 0x20)


def hexrgb(c):
    return RGBColor.from_string(c.lstrip('#').upper())


def fn_rgb(f):
    """Colour of a run by Fn: blue (0) -> light -> red (1), as in the GIFs."""
    a, m, b = np.array([0x2C, 0xB1, 0xD9]), np.array([0xE6, 0xE6, 0xE6]), np.array([0xE8, 0x60, 0x4C])
    c = a + (m - a) * f / 0.5 if f <= 0.5 else m + (b - m) * (f - 0.5) / 0.5
    return RGBColor(*[int(round(x)) for x in c])


def dec(x, nd=2):
    s = ('%.*f' % (nd, x)).rstrip('0').rstrip('.') if nd else '%d' % round(x)
    return s.replace('.', ',')


# ---------------------------------------------------------------- text
def runs(par, text, size, color=WHITE, bold=False):
    """Mini markup: *x* italic, _{..} subscript, ^{..} superscript, **..** bold."""
    i, it, bd = 0, False, bold
    buf = ''

    def flush(t, base=None, it_=None, bd_=None):
        if not t:
            return
        r = par.add_run()
        r.text = t
        f = r.font
        f.name, f.size, f.color.rgb = FONT, Pt(size), color
        f.italic = it if it_ is None else it_
        f.bold = bd if bd_ is None else bd_
        rPr = r._r.get_or_add_rPr()
        for tag in ('a:latin', 'a:cs', 'a:ea'):
            e = rPr.find(qn(tag))
            if e is None:
                e = etree.SubElement(rPr, qn(tag))
            e.set('typeface', FONT)
        if base:
            rPr.set('baseline', base)
    while i < len(text):
        if text.startswith('**', i):
            flush(buf); buf = ''; bd = not bd; i += 2
        elif text[i] == '*':
            flush(buf); buf = ''; it = not it; i += 1
        elif text[i] in '_^' and i + 1 < len(text) and text[i + 1] == '{':
            flush(buf); buf = ''
            j = text.index('}', i)
            flush(text[i + 2:j], '-25000' if text[i] == '_' else '30000')
            i = j + 1
        else:
            buf += text[i]; i += 1
    flush(buf)


def textbox(slide, x, y, w, h, paras, size=18, color=WHITE, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP,
            spacing=1.0, space_after=6, bullets=False):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for side in ('margin_left', 'margin_right', 'margin_top', 'margin_bottom'):
        setattr(tf, side, Inches(0.02))
    for k, p in enumerate(paras):
        if isinstance(p, str):
            p = dict(text=p)
        par = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        par.alignment = p.get('align', align)
        par.line_spacing = spacing
        par.space_after = Pt(p.get('after', space_after))
        if bullets or p.get('bullet'):
            pPr = par._p.get_or_add_pPr()
            pPr.set('marL', str(Emu(Inches(0.28)))); pPr.set('indent', str(-Emu(Inches(0.28))))
            bu = etree.SubElement(pPr, qn('a:buFont')); bu.set('typeface', 'Arial')
            ch = etree.SubElement(pPr, qn('a:buChar')); ch.set('char', '•')
        runs(par, p['text'], p.get('size', size), p.get('color', color), p.get('bold', False))
    return tb


def title(slide, text, sub=None):
    textbox(slide, 0.45, 0.22, 12.4, 0.7, [dict(text=text, size=30)], anchor=MSO_ANCHOR.TOP)
    if sub:
        textbox(slide, 0.47, 0.82, 12.4, 0.4, [dict(text=sub, size=16, color=GREY)])


def page_no(slide, n):
    textbox(slide, 12.3, 7.05, 0.8, 0.3, [dict(text=str(n), size=11, color=DIM, align=PP_ALIGN.RIGHT)])


def picture(slide, path, x, y, w=None, h=None):
    kw = {}
    if w is not None:
        kw['width'] = Inches(w)
    if h is not None:
        kw['height'] = Inches(h)
    return slide.shapes.add_picture(path, Inches(x), Inches(y), **kw)


def rect(slide, x, y, w, h, fill=CARD, line=None, rounded=True):
    s = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE,
                               Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid(); s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line; s.line.width = Pt(1.25)
    if rounded:
        s.adjustments[0] = 0.06
    s.shadow.inherit = False
    return s


def card_picture(slide, path, x, y, w, h, pad=0.08):
    """A white figure from a paper on a light card (keeps its aspect, centred)."""
    rect(slide, x, y, w, h, fill=RGBColor(0xFF, 0xFF, 0xFF))
    from PIL import Image
    iw, ih = Image.open(path).size
    k = min((w - 2 * pad) / iw, (h - 2 * pad) / ih)
    pw, ph = iw * k, ih * k
    return picture(slide, path, x + (w - pw) / 2, y + (h - ph) / 2, w=pw)


def new_slide(prs):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = BG
    return s


# ---------------------------------------------------------------- click swap
def click_swap(slide, out_shapes, in_shapes):
    """One click: out_shapes fade out, in_shapes fade in (in_shapes hidden
    before).  Builds the slide's p:timing."""
    P = 'http://schemas.openxmlformats.org/presentationml/2006/main'
    ids = iter(range(3, 1000))

    def eff(spid, entr, first):
        cid = next(ids)
        node = 'clickEffect' if first else 'withEffect'
        cls = 'entr' if entr else 'exit'
        vis = ('<p:set><p:cBhvr><p:cTn id="%d" dur="1" fill="hold"><p:stCondLst><p:cond delay="%d"/>'
               '</p:stCondLst></p:cTn><p:tgtEl><p:spTgt spid="%d"/></p:tgtEl><p:attrNameLst>'
               '<p:attrName>style.visibility</p:attrName></p:attrNameLst></p:cBhvr><p:to><p:strVal val="%s"/>'
               '</p:to></p:set>') % (next(ids), 0 if entr else 499, spid, 'visible' if entr else 'hidden')
        fade = ('<p:animEffect transition="%s" filter="fade"><p:cBhvr><p:cTn id="%d" dur="500"/>'
                '<p:tgtEl><p:spTgt spid="%d"/></p:tgtEl></p:cBhvr></p:animEffect>') % (
                    'in' if entr else 'out', next(ids), spid)
        body = vis + fade if entr else fade + vis
        return ('<p:par><p:cTn id="%d" presetID="10" presetClass="%s" presetSubtype="0" fill="hold" '
                'grpId="0" nodeType="%s"><p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>%s'
                '</p:childTnLst></p:cTn></p:par>') % (cid, cls, node, body)
    effs = []
    first = True
    for sh in out_shapes:
        effs.append(eff(sh.shape_id, False, first)); first = False
    for sh in in_shapes:
        effs.append(eff(sh.shape_id, True, False))
    blds = ''.join('<p:bldP spid="%d" grpId="0"/>' % sh.shape_id
                   for sh in list(out_shapes) + list(in_shapes) if sh.has_text_frame)
    xml = ('<p:timing xmlns:p="%s"><p:tnLst><p:par><p:cTn id="1" dur="indefinite" restart="never" '
           'nodeType="tmRoot"><p:childTnLst><p:seq concurrent="1" nextAc="seek"><p:cTn id="2" '
           'dur="indefinite" nodeType="mainSeq"><p:childTnLst><p:par><p:cTn id="900" fill="hold">'
           '<p:stCondLst><p:cond delay="indefinite"/></p:stCondLst><p:childTnLst><p:par><p:cTn id="901" '
           'fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>%s</p:childTnLst>'
           '</p:cTn></p:par></p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn><p:prevCondLst>'
           '<p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:prevCondLst>'
           '<p:nextCondLst><p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond>'
           '</p:nextCondLst></p:seq></p:childTnLst></p:cTn></p:par></p:tnLst>%s</p:timing>') % (
               P, ''.join(effs), ('<p:bldLst>%s</p:bldLst>' % blds) if blds else '')
    slide._element.append(etree.fromstring(xml))


# ---------------------------------------------------------------- charts
def style_chart(chart, size=12):
    chart.font.name = FONT
    chart.font.size = Pt(size)
    chart.font.color.rgb = GREY
    cs = chart._chartSpace
    # no fill, no border on the chart and plot areas
    for parent in (cs, cs.find(qn('c:chart')).find(qn('c:plotArea'))):
        sp = parent.find(qn('c:spPr'))
        if sp is None:
            sp = etree.SubElement(parent, qn('c:spPr'))
        for ch in list(sp):
            sp.remove(ch)
        etree.SubElement(sp, qn('a:noFill'))
        ln = etree.SubElement(sp, qn('a:ln')); etree.SubElement(ln, qn('a:noFill'))
    # chartSpace spPr must come before txPr
    sp = cs.find(qn('c:spPr')); tx = cs.find(qn('c:txPr'))
    if tx is not None and sp is not None and list(cs).index(sp) > list(cs).index(tx):
        cs.remove(sp); tx.addprevious(sp)


def style_axis(ax, title=None, size=12, fmt=None, minimum=None, maximum=None, major=None):
    ax.has_major_gridlines = False
    ax.has_minor_gridlines = False
    ax.format.line.color.rgb = DIM
    ax.format.line.width = Pt(1.25)
    ax.tick_labels.font.size = Pt(size)
    ax.tick_labels.font.color.rgb = GREY
    ax.tick_labels.font.name = FONT
    if fmt:
        ax.tick_labels.number_format = fmt
        ax.tick_labels.number_format_is_linked = False
    if minimum is not None:
        ax.minimum_scale = minimum
    if maximum is not None:
        ax.maximum_scale = maximum
    if major is not None:
        ax.major_unit = major
    if title:
        ax.has_title = True
        tf = ax.axis_title.text_frame
        tf.text = ''
        runs(tf.paragraphs[0], title, size + 1, GREY)


def xy_chart(slide, x, y, w, h, series, xt, yt, xr=None, yr=None, xfmt='0.0', yfmt='0', markers=True,
             legend=False, size=12):
    data = XyChartData()
    for name, xs, ys, _ in series:
        s = data.add_series(name)
        for a, b in zip(xs, ys):
            s.add_data_point(a, b)
    gf = slide.shapes.add_chart(XL_CHART_TYPE.XY_SCATTER_LINES if markers else XL_CHART_TYPE.XY_SCATTER_LINES_NO_MARKERS,
                                Inches(x), Inches(y), Inches(w), Inches(h), data)
    ch = gf.chart
    style_chart(ch, size)
    ch.has_legend = legend
    if legend:
        ch.legend.position = XL_LEGEND_POSITION.TOP
        ch.legend.include_in_layout = False
        ch.legend.font.size = Pt(size); ch.legend.font.color.rgb = GREY; ch.legend.font.name = FONT
    for s, (_, _, _, col) in zip(ch.plots[0].series, series):
        s.smooth = False
        s.format.line.color.rgb = col
        s.format.line.width = Pt(2.75)
        if markers:
            s.marker.style = XL_MARKER_STYLE.CIRCLE
            s.marker.size = 8
            s.marker.format.fill.solid(); s.marker.format.fill.fore_color.rgb = col
            s.marker.format.line.color.rgb = col
    style_axis(ch.category_axis, xt, size, xfmt, *(xr or (None, None, None)))
    style_axis(ch.value_axis, yt, size, yfmt, *(yr or (None, None, None)))
    return ch


# ---------------------------------------------------------------- the slides
def build(out):
    from fndata import Run, SERIES
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    n = 0

    def nxt():
        nonlocal n
        n += 1
        s = new_slide(prs)
        if n > 1:
            page_no(s, n)
        return s

    # 1. title
    s = nxt()
    textbox(s, 0.8, 1.55, 11.7, 2.2, [dict(text='Роль ориентации гидридов в разрушении циркониевых сплавов '
                                             'по данным конечно-элементного моделирования', size=36)],
            spacing=1.05)
    textbox(s, 0.8, 4.1, 11.7, 0.5, [dict(text='**Р.А. Ильясов**, М.И. Петров, М.Г. Исаенкова, В.В. Михальчик',
                                          size=22)])
    textbox(s, 0.8, 4.75, 11.7, 0.9, [dict(text='Национальный исследовательский ядерный университет «МИФИ», Москва',
                                           size=18, color=GREY, after=2),
                                      dict(text='il.re0310@mail.ru', size=16, color=TAN)])
    textbox(s, 0.8, 6.55, 6, 0.4, [dict(text='Москва, 2026', size=16, color=GREY)])

    # 2. introduction
    s = nxt()
    title(s, 'Введение')
    textbox(s, 0.8, 1.45, 11.8, 5.6, [
        'Циркониевые сплавы широко применяются в качестве конструкционного материала активной зоны '
        'водо-водяных энергетических реакторов (ВВЭР).',
        'При эксплуатации элементы активной зоны контактируют с водным теплоносителем и насыщаются водородом. '
        'При превышении предела растворимости водород выделяется в виде хрупких δ-гидридов циркония, '
        'которые снижают пластичность материала и ускоряют деградацию конструктивных элементов.',
        'Существенное влияние на пластичность оказывает ориентация гидридов. Механические испытания '
        'в сочетании с численным моделированием позволяют установить, какие механизмы разрушения '
        'определяют снижение пластичности при разной ориентации гидридов.',
        dict(text='**Цель работы:** установить влияние ориентации гидридов на характер разрушения '
                  'и пластичность циркониевого сплава методом конечных элементов.', color=WHITE),
    ], size=20, space_after=14, spacing=1.08)

    # 3. delta hydrides: deck of cards -> (click) -> Fn measurement
    s = nxt()
    gw = 9.35
    gy = 0.55
    a1 = picture(s, os.path.join(GIFS, 'deck_of_cards.gif'), 0.1, gy, w=gw)
    t1 = textbox(s, 9.6, 1.55, 3.45, 5.0, [
        'Нанопластины δ-ZrH_{1,66} выделяются в базисной плоскости α-Zr',
        'Мезомасштабный гидрид собран из нанопластин ступенями, как колода карт, и лежит вдоль плоскости '
        '(10̅1̅7) под углом 14,7° к базисной',
        'На шлифах видны именно мезомасштабные гидриды'], size=17, bullets=True, space_after=12)
    b1 = picture(s, os.path.join(SCRATCH, 'fn_measure.gif'), 0.1, gy, w=gw)
    t2 = textbox(s, 9.6, 1.55, 3.45, 5.0, [
        'Ориентацию гидридов характеризует коэффициент *F_{n}*',
        '*F_{n}*: доля суммарной длины гидридов, отклонённых от радиального направления трубы не более '
        'чем на 45°',
        '*F_{n}* = 0: все гидриды тангенциальные, *F_{n}* = 1: все радиальные'], size=17, bullets=True,
        space_after=12)
    title(s, 'δ-гидриды циркония')
    click_swap(s, [a1, t1], [b1, t2])

    # 4. what sets the orientation: texture -> (click) -> stress reorientation
    s = nxt()
    a1 = picture(s, os.path.join(GIFS, 'texture.gif'), 0.1, gy, w=gw)
    t1 = textbox(s, 9.6, 1.55, 3.45, 5.0, [
        'Гидриды выделяются вдоль базисных плоскостей α-Zr, поэтому их ориентация следует текстуре',
        'В трубах базисные полюса направлены близко к радиальному направлению, поэтому гидриды выделяются преимущественно тангенциально'],
        size=17, bullets=True, space_after=12)
    b1 = picture(s, os.path.join(GIFS, 'reorient.gif'), 0.1, gy, w=gw)
    t2 = textbox(s, 9.6, 1.55, 3.45, 5.0, [
        'При нагреве гидриды растворяются, при охлаждении выделяются вновь',
        'Растягивающее окружное напряжение σ_{θ} при охлаждении разворачивает выделяющиеся гидриды радиально',
        'После снятия нагрузки радиальные гидриды сохраняются'], size=17, bullets=True, space_after=12)
    title(s, 'Что влияет на ориентацию гидридов?')
    click_swap(s, [a1, t1], [b1, t2])

    # 5. orientation and fracture in experiments (draft)
    s = nxt()
    title(s, 'Ориентация гидридов и характер разрушения', 'черновик: данные испытаний кольцевых образцов')
    fr = os.path.join(SCRATCH, 'fracto_rings.jpg')
    card_picture(s, fr, 0.45, 1.45, 6.6, 2.55)
    textbox(s, 0.45, 4.05, 3.3, 0.4, [dict(text='вязкий излом', size=16, color=GREY, align=PP_ALIGN.CENTER)])
    textbox(s, 3.75, 4.05, 3.3, 0.4, [dict(text='хрупкий излом', size=16, color=GREY, align=PP_ALIGN.CENTER)])
    card_picture(s, os.path.join(SCRATCH, 'brittle_fn.png'), 7.35, 1.45, 2.65, 2.65)
    card_picture(s, os.path.join(SCRATCH, 'elong_fn.png'), 10.2, 1.45, 2.7, 2.65)
    textbox(s, 7.35, 4.15, 2.65, 0.7, [dict(text='доля хрупкой составляющей излома', size=14, color=GREY,
                                            align=PP_ALIGN.CENTER)])
    textbox(s, 10.2, 4.15, 2.7, 0.7, [dict(text='относительное удлинение колец', size=14, color=GREY,
                                           align=PP_ALIGN.CENTER)])
    textbox(s, 0.8, 5.05, 11.8, 2.0, [
        'С ростом *F_{n}* излом переходит от вязкого (ямочного) к хрупкому',
        'Переход резкий: при *F_{n}* ≈ 0,25–0,3 доля хрупкой составляющей растёт от 0,1 до 0,9',
        'Удлинение колец 23–27 % при *F_{n}* < 0,2; при *F_{n}* > 0,25 оно снижается вплоть до нуля'],
        size=19, bullets=True, space_after=8)

    # 6. the numerical experiment
    s = nxt()
    picture(s, os.path.join(SCRATCH, 'fn_block.gif'), 0.55, 0.1, w=11.0)
    title(s, 'Постановка численного эксперимента')
    r0 = Run(0.0)
    textbox(s, 0.8, 6.32, 11.8, 1.0, [
        'Блок 4,0 × 3,6 × 1,6 мм (*T* × *R* × *L*), 12 одинаковых дисков гидрида диаметром 0,6 мм и толщиной '
        '0,09 мм; меняется только ориентация дисков',
        'Растяжение вдоль *T* (окружное напряжение σ_{θ}) заданным перемещением грани; '
        'тетраэдральная сетка ≈ %d тыс. элементов, когезионные элементы на границе фаз' % round(len(r0.T) / 1000)],
        size=15, space_after=3, color=GREY)

    # 7. materials
    s = nxt()
    title(s, 'Свойства материалов')
    cards = [('Матрица α-Zr', ['*E* = 108,8 ГПа', 'σ_{т} = 500 МПа'], RGBColor(0xC9, 0xD1, 0xDB)),
             ('Гидрид δ-ZrH_{1,66}', ['*E* = 130 ГПа', 'σ_{т} = 600 МПа'], MECH['cleavage']),
             ('Граница фаз', ['*T*_{0} = 300 МПа', '*G*_{c} = 2,4 кДж/м^{2}'], MECH['debond'])]
    for k, (head, lines, col) in enumerate(cards):
        y = 1.3 + k * 1.17
        rect(s, 0.45, y, 3.2, 1.07, fill=CARD, line=col)
        textbox(s, 0.6, y + 0.06, 3.0, 0.4, [dict(text=head, size=17, color=col)])
        textbox(s, 0.6, y + 0.47, 1.45, 0.55, [dict(text=lines[0], size=15)])
        textbox(s, 2.0, y + 0.47, 1.6, 0.55, [dict(text=lines[1], size=15)])
    eta = [0, 0.33, 0.5, 0.67, 1.0, 1.5]
    efm = [30, 12, 7, 4.25, 2.5, 1.5]
    efh = [12, 7, 4.5, 2.8, 1.5, 1.0]
    textbox(s, 4.0, 1.25, 4.6, 0.4, [dict(text='Деформация начала повреждения', size=16, color=GREY,
                                          align=PP_ALIGN.CENTER)])
    xy_chart(s, 3.95, 1.6, 4.7, 2.9, [('матрица', eta, efm, RGBColor(0xC9, 0xD1, 0xDB)),
                                       ('гидрид', eta, efh, MECH['cleavage'])],
             'трёхосность η', 'ε\u0304_{D}^{pl}, %', xr=(0, 1.5, 0.5), yr=(0, 30, 10), xfmt='0.0', yfmt='0',
             legend=True)
    # T(delta): a drawn schematic
    ox, oy, sx, sy = 9.35, 4.2, 3.3, 2.3
    p0, pk, pf = (ox, oy), (ox + 0.12 * sx, oy - 0.85 * sy), (ox + 0.85 * sx, oy)
    E = lambda q: (int(Inches(q[0])), int(Inches(q[1])))      # whole EMU: the builder rounds its units
    ff = s.shapes.build_freeform(*E(p0), scale=1.0)
    ff.add_line_segments([E(pk), E(pf), E(p0)])
    tri = ff.convert_to_shape()
    tri.fill.solid(); tri.fill.fore_color.rgb = MECH['debond']
    sf = tri.fill._xPr.find(qn('a:solidFill'))[0]
    al = etree.SubElement(sf, qn('a:alpha')); al.set('val', '35000')
    tri.line.color.rgb = MECH['debond']; tri.line.width = Pt(2.75)
    for (x0, y0, x1, y1) in ((ox, oy, ox + sx + 0.25, oy), (ox, oy, ox, oy - sy - 0.2)):
        c = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x0), Inches(y0), Inches(x1), Inches(y1))
        c.line.color.rgb = DIM; c.line.width = Pt(1.5)
        te = etree.SubElement(c.line._get_or_add_ln(), qn('a:tailEnd')); te.set('type', 'triangle')
    textbox(s, ox - 0.45, oy - sy - 0.35, 0.4, 0.4, [dict(text='*T*', size=17, color=GREY)])
    textbox(s, ox + sx + 0.25, oy - 0.2, 0.4, 0.4, [dict(text='δ', size=17, color=GREY)])
    textbox(s, ox - 0.55, pk[1] - 0.18, 0.5, 0.4, [dict(text='*T*_{0}', size=15, color=GREY)])
    textbox(s, pf[0] - 0.15, oy + 0.03, 0.5, 0.4, [dict(text='δ_{f}', size=15, color=GREY)])
    textbox(s, ox + 0.25 * sx, oy - 0.42 * sy, 0.8, 0.4, [dict(text='*G*_{c}', size=18, color=WHITE)])
    textbox(s, 8.85, 1.25, 4.3, 0.4, [dict(text='Когезионный закон границы фаз', size=16, color=GREY,
                                           align=PP_ALIGN.CENTER)])
    # the damage laws, as written in the literature; no evolution law for D
    fy = 4.75
    for x, img, hgt, cap in (
            (0.6, 'f_omega.png', 1.0, 'ω: исчерпание пластичности; при ω = 1 начинается повреждение'),
            (4.05, 'f_eta.png', 1.2, 'трёхосность напряжённого состояния; σ_{1}, σ_{2}, σ_{3}: главные напряжения'),
            (9.75, 'f_sigma.png', 0.33, '*D*: повреждённость, растёт после ω = 1; при *D* = 1 элемент удаляется')):
        picture(s, os.path.join(ASSETS, img), x, fy + min(0., 1.0 - hgt) + max(0., 1.0 - hgt) / 2, h=hgt)
        textbox(s, x, fy + 1.12, 5.2 if img == 'f_eta.png' else 3.3, 0.9, [dict(text=cap, size=14, color=GREY)])
    textbox(s, 0.6, 6.85, 12, 0.4, [dict(text='ε\u0304^{pl}: эквивалентная пластическая деформация; '
                                              'ε\u0304_{D}^{pl}(η): её значение в начале повреждения; '
                                              'σ\u0304: напряжение без повреждения', size=13, color=DIM)])

    # 8. the series
    s = nxt()
    picture(s, os.path.join(GIFS, 'res_sections_seq.gif'), 0, 0, w=W)
    title(s, 'Серия расчётов с различным *F_{n}*')

    # 9. fracture surfaces + mechanism GIFs as the colour key
    s = nxt()
    title(s, 'Поверхности разрушения')
    runs9 = [(0.0, 'fn0'), (0.25, 'fn25v1e-3'), (0.58, 'fn60'), (1.0, 'fn100')]
    for k, (f, name) in enumerate(runs9):
        x = 0.35 + k * 3.18
        picture(s, os.path.join(SCRATCH, 'fs', 'fs_%s.png' % name), x, 1.35, w=3.1)
        lab = '0' if f == 0 else '1' if f == 1 else dec(f)
        textbox(s, x, 1.05, 3.1, 0.4, [dict(text='*F_{n}* = ' + lab, size=18, color=fn_rgb(f),
                                            align=PP_ALIGN.CENTER)])
    for k, (gif, key, name) in enumerate((('mech1_ductile.gif', 'ductile', 'вязкий разрыв матрицы'),
                                          ('mech2_cleave.gif', 'cleavage', 'скол гидрида'),
                                          ('mech3_debond.gif', 'debond', 'отслоение по границе фаз'))):
        x = 0.9 + k * 4.15
        picture(s, os.path.join(GIFS, gif), x, 4.05, h=2.55)
        rect(s, x + 2.7, 4.95, 0.32, 0.32, fill=MECH[key], rounded=False)
        textbox(s, x + 2.7, 5.35, 1.35, 1.0, [dict(text=name, size=15, color=MECH[key])])

    # 10. numbers
    s = nxt()
    title(s, 'Количественные характеристики')
    runs10 = [Run(fn) for fn, _ in SERIES]
    F = [r.k_radial / 12 for r in runs10]
    ef = [100 * r.eps_f10 for r in runs10]
    sb = [r.peak for r in runs10]
    wk = []
    for r in runs10:
        c = r.curve
        m = c[:, 0] <= r.eps_f10 + 1e-9
        x_, y_ = np.r_[0, c[m, 0]], np.r_[0, c[m, 1]]
        wk.append(float(np.sum(0.5 * (y_[1:] + y_[:-1]) * np.diff(x_))))
    blue = RGBColor(0x8F, 0xC8, 0xEA)
    xy_chart(s, 0.35, 1.2, 6.2, 2.85, [('', F, ef, blue)], '*F_{n}*', 'ε_{f}, %',
             xr=(0, 1, 0.2), yr=(0, 14, 2), xfmt='0.0', yfmt='0')
    xy_chart(s, 6.8, 1.2, 6.2, 2.85, [('', F, wk, blue)], '*F_{n}*', '*W*, МДж/м^{3}',
             xr=(0, 1, 0.2), yr=(0, 80, 20), xfmt='0.0', yfmt='0')
    xy_chart(s, 0.35, 4.3, 6.2, 2.85, [('', F, sb, blue)], '*F_{n}*', 'σ_{в}, МПа',
             xr=(0, 1, 0.2), yr=(500, 600, 25), xfmt='0.0', yfmt='0')
    # brittle share of the fracture surface seen along T (projection on the R-L plane, as in
    # fractography): debonding + hydride cleavage; Fn = 0.33 from its completed rerun
    from fracture_surface import projected_brittle
    from pptx.enum.dml import MSO_LINE_DASH_STYLE
    names = [n_ if fn != 0.33 else 'fn33v4e-3' for fn, n_ in SERIES]
    br = [projected_brittle(Run(n_))[0] for n_ in names]
    fx = np.linspace(0, 0.5, 26)                        # rings, Petrov et al. 2025 (fig. 6): logistic fit
    fexp = 0.05 + 0.87 / (1 + np.exp(-(fx - 0.265) / 0.018))
    # the longer series first: LibreOffice drops a first series that has fewer points
    ch = xy_chart(s, 6.8, 4.3, 6.2, 2.85, [('эксперимент', list(fx), list(fexp), ACC), ('модель', F, br, blue)],
                  '*F_{n}*', 'доля хрупкой составляющей', xr=(0, 1, 0.2), yr=(0, 1, 0.2), xfmt='0.0',
                  yfmt='0.0', legend=True)
    exp = ch.plots[0].series[0]
    exp.marker.style = XL_MARKER_STYLE.NONE
    exp.format.line.dash_style = MSO_LINE_DASH_STYLE.DASH
    exp.format.line.width = Pt(2.25)
    textbox(s, 6.8, 7.12, 5.4, 0.3, [dict(text='излом в проекции на плоскость *R*–*L*; хрупкое: отслоение '
                                              'и скол гидрида', size=12, color=DIM, align=PP_ALIGN.CENTER)])

    # 11. conclusions (draft)
    s = nxt()
    title(s, 'Выводы', 'черновик')
    i0, i1 = 0, len(runs10) - 1
    textbox(s, 0.8, 1.45, 11.8, 5.6, [
        'Построена конечно-элементная модель разрушения циркониевого сплава с пластинчатыми гидридами: '
        'вязкое разрушение матрицы и скол гидрида описываются исчерпанием пластичности и повреждённостью, '
        'граница фаз описывается когезионным законом',
        'С ростом *F_{n}* от 0 до 1 деформация разрушения снижается с %s до %s %%, работа разрушения с %s до '
        '%s МДж/м^{3}, а прочность меняется лишь на %s %%' % (
            dec(ef[i0], 1), dec(ef[i1], 1), dec(wk[i0], 0), dec(wk[i1], 0),
            dec(100 * (sb[i0] - sb[i1]) / sb[i0], 0)),
        'Тангенциальные гидриды раскалываются, а трещина идёт вязко по матрице; радиальные гидриды '
        'отслаиваются, и трещина проходит от пластины к пластине по перемычкам',
        'Доля хрупкой составляющей излома (в проекции на плоскость *R*–*L*) растёт с %s при *F_{n}* = 0 до %s '
        'при *F_{n}* = 1; в испытаниях колец рост резче и достигает 0,9: в модели 12 отдельных пластин не '
        'образуют сквозной сетки гидридов' % (dec(br[i0], 2), dec(br[i1], 2)),
        'Модель качественно воспроизводит наблюдаемое в испытаниях снижение пластичности с ростом *F_{n}*'],
        size=19, bullets=True, space_after=12)

    # 12. thanks
    s = nxt()
    textbox(s, 0.8, 2.6, 11.7, 1.0, [dict(text='Спасибо за внимание!', size=40, align=PP_ALIGN.CENTER)])
    textbox(s, 0.8, 4.0, 11.7, 1.2, [
        dict(text='**Р.А. Ильясов**, М.И. Петров, М.Г. Исаенкова, В.В. Михальчик', size=18, align=PP_ALIGN.CENTER),
        dict(text='НИЯУ МИФИ, il.re0310@mail.ru', size=16, color=GREY, align=PP_ALIGN.CENTER)])

    prs.save(out)
    print('saved', out, n, 'slides')


if __name__ == '__main__':
    build(sys.argv[1] if len(sys.argv) > 1 else 'deck.pptx')
