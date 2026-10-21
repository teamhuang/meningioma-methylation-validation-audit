from pathlib import Path
from io import BytesIO
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(r'E:\脑膜瘤npj')
OUT = ROOT / 'npjDM_submission_candidate_v6'
FINAL_PDF = OUT / 'Supplementary_Information.pdf'

figs = [
 ('Supplementary Figure S1.', 'Training-reference fingerprint audit. Per-cohort distances to the nearest reference profile and identifier checks; all 866 external specimens are shown against the exact-match thresholds used for the author releases.'),
 ('Supplementary Figure S2.', 'Exploratory decision curves with corrected censoring weights, presented as relative comparisons without a clinical threshold.'),
 ('Supplementary Figure S3.', 'Robust-panel transfer. The corrected panel, shared-probe counts and rejection flags are shown alongside the documented probe-selection error and its correction.'),
 ('Supplementary Figure S4.', 'Longitudinal ancillary series. Per-cluster labels are compared with the constant-class comparator; paired probe-difference values and donor-fingerprint conflict flags are also shown. This computational grouping is not treated as patient-confirmed.'),
]
data = [
 ('Supplementary Data 1.', 'Corrected robust-panel profile-level results and probe-level metadata. The compressed CSV contains profile-level transfer results and probe-level fields used to document the corrected panel.'),
 ('Supplementary Data 2.', 'Longitudinal specimen-to-cluster mapping, donor-fingerprint conflict flags and the constant-class baseline for the 42-specimen ancillary series.'),
 ('Supplementary Data 3.', 'Identifier-proven UCSF molecular-profile linkage and cohort metadata used in the identity audit, including the record-level fields supporting the analyzed linkage.'),
]
W,H=landscape(A4)
styles=getSampleStyleSheet()
title_style=ParagraphStyle('TitleX',parent=styles['Title'],fontName='Helvetica-Bold',fontSize=18,leading=22,alignment=TA_CENTER,spaceAfter=9)
subtitle_style=ParagraphStyle('SubX',parent=styles['Normal'],fontName='Helvetica',fontSize=9,leading=12,alignment=TA_CENTER,textColor=colors.HexColor('#505050'))
head_style=ParagraphStyle('HeadX',parent=styles['Heading2'],fontName='Helvetica-Bold',fontSize=12,leading=15,spaceBefore=6,spaceAfter=5)
body_style=ParagraphStyle('BodyX',parent=styles['BodyText'],fontName='Helvetica',fontSize=9,leading=12,spaceAfter=7)
caption_style=ParagraphStyle('CaptionX',parent=styles['BodyText'],fontName='Helvetica',fontSize=10,leading=13,spaceAfter=0)

writer=PdfWriter()

# One overview page lists the separately supplied data files; each figure page below has its own legend.
buf=BytesIO(); c=canvas.Canvas(buf,pagesize=(W,H))
title=Paragraph('Supplementary Information',title_style)
tw,th=title.wrap(W-40*mm,H)
title.drawOn(c,20*mm,H-22*mm-th)
sub=Paragraph('Identity, training independence and selective prediction shape the transportability of a meningioma DNA methylation classifier',subtitle_style)
sw,sh=sub.wrap(W-40*mm,H)
sub.drawOn(c,20*mm,H-27*mm-th-sh)
y=H-42*mm-th-sh
head=Paragraph('Supplementary Data files',head_style); hw,hh=head.wrap(W-40*mm,H); head.drawOn(c,20*mm,y-hh); y-=hh+3*mm
for label,desc in data:
    para=Paragraph(f'<b>{label}</b> {desc}',body_style)
    pw,ph=para.wrap(W-42*mm,H); para.drawOn(c,21*mm,y-ph); y-=ph+3*mm
note=Paragraph('The four supplementary figures follow in numerical order. The data files are supplied separately with the manuscript.',body_style)
nw,nh=note.wrap(W-42*mm,H); note.drawOn(c,21*mm,y-nh)
c.showPage(); c.save(); buf.seek(0)
writer.append(PdfReader(buf))

for i,(label,desc) in enumerate(figs,1):
    source=PdfReader(str(OUT/f'Supplementary_Figure_S{i}.pdf'))
    if len(source.pages)!=1: raise ValueError(f'Supplementary Figure S{i} must be one page')
    sb=source.pages[0]
    sw=float(sb.mediabox.width); sh=float(sb.mediabox.height)
    buf=BytesIO(); c=canvas.Canvas(buf,pagesize=(W,H))
    cap=Paragraph(f'<b>{label}</b> {desc}',caption_style)
    cw,ch=cap.wrap(W-34*mm,H)
    cap.drawOn(c,17*mm,H-16*mm-ch)
    maxw=W-28*mm; maxh=H-45*mm-ch
    scale=min(maxw/sw,maxh/sh)
    tx=(W-sw*scale)/2
    ty=17*mm+(maxh-sh*scale)/2
    c.showPage(); c.save(); buf.seek(0)
    page=PdfReader(buf).pages[0]
    page.merge_transformed_page(sb,Transformation().scale(scale).translate(tx,ty),over=True)
    writer.add_page(page)

with FINAL_PDF.open('wb') as f: writer.write(f)
print(f'Created {FINAL_PDF} with {len(PdfReader(str(FINAL_PDF)).pages)} pages')
