#!/usr/bin/env python3
"""Create a short plain-language companion; requires reportlab in the local environment."""
from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib.pagesizes import letter

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'output/pdf/decision_review_simple_guide.pdf'
INK=colors.HexColor('#243746')
BLUE=colors.HexColor('#216b8e')
MUTED=colors.HexColor('#566775')
LIGHT=colors.HexColor('#edf4f7')
AMBER=colors.HexColor('#fff4df')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='GuideTitle',fontName='Helvetica-Bold',fontSize=23,leading=27,textColor=INK,spaceAfter=6))
styles.add(ParagraphStyle(name='GuideSub',fontName='Helvetica',fontSize=10.5,leading=15,textColor=MUTED,spaceAfter=10))
styles.add(ParagraphStyle(name='FigureTitle',fontName='Helvetica-Bold',fontSize=15,leading=20,textColor=BLUE,spaceBefore=4,spaceAfter=3))
styles.add(ParagraphStyle(name='FileLabel',fontName='Courier',fontSize=9,leading=13,textColor=MUTED,spaceAfter=6))
styles.add(ParagraphStyle(name='GuideBody',fontName='Helvetica',fontSize=11.2,leading=15.5,textColor=INK,spaceAfter=5))
styles.add(ParagraphStyle(name='Panel',fontName='Helvetica',fontSize=11.2,leading=15.5,textColor=INK,leftIndent=12,firstLineIndent=-12,spaceAfter=4))
styles.add(ParagraphStyle(name='SmallGuide',fontName='Helvetica',fontSize=9.5,leading=13.5,textColor=MUTED,spaceAfter=6))


def p(text,style='GuideBody'):
    return Paragraph(text,styles[style])


def box(text,background=LIGHT):
    t=Table([[p(text)]],colWidths=[516])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),background),('LEFTPADDING',(0,0),(-1,-1),11),
        ('RIGHTPADDING',(0,0),(-1,-1),11),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),3)]))
    return t


def section(title,filename,panels,takeaway,caution):
    story=[p(title,'FigureTitle'),p(filename,'FileLabel')]
    story.extend(p('<b>'+label+'</b> '+text,'Panel') for label,text in panels)
    story.extend([Spacer(1,3),p('<b>Main point:</b> '+takeaway),p('<b>Do not conclude:</b> '+caution,'SmallGuide'),Spacer(1,9)])
    return KeepTogether(story)


def footer(canvas,doc):
    canvas.setStrokeColor(colors.HexColor('#d8e2e7'))
    canvas.line(48,42,564,42)
    canvas.setFillColor(MUTED)
    canvas.setFont('Helvetica',8)
    canvas.drawString(48,28,'Glia enhancer design | Companion to reports/decision_review | October 15, 2026')
    canvas.drawRightString(564,28,str(doc.page))


def main():
    DEST.parent.mkdir(parents=True,exist_ok=True)
    doc=SimpleDocTemplate(str(DEST),pagesize=letter,rightMargin=48,leftMargin=48,topMargin=38,bottomMargin=54,
        title='Simple guide to the decision-review figures',author='Glia enhancer design project')
    story=[p('How to read the four figures','GuideTitle'),
        p('A simple guide for the discussion with Thanh. Open the matching PNG or HTML figure alongside this PDF.','GuideSub'),
        section('1. Does the MG model learn useful signal?','01_model_validation.png',[
            ('A - Counts:','Blue dots above the red bias-only line are better. All three model seeds improve peak-count prediction.'),
            ('B - Profile shape:','Here, lower is better. Blue dots below the red line mean the model predicts peak shape better than enzyme bias alone.'),
            ('C - Donors:','Blue bars exceed gray bars in all four donors. The model is useful across the donors that contributed data.'),
            ('D - Background:','Blue dots are below the red line, which is worse for this metric. Background-count prediction remains a limitation.')],
            'The model learns useful MG sequence signal, with stable overall performance across three training seeds.',
            'This does not prove performance in an unseen donor or that each designed edit is stable across seeds.'),
        section('2. Are the natural parent regions MG-enriched?','02_parent_specificity.png',[
            ('A - Signal:','Points above the diagonal have higher measured MG accessibility than their strongest off-target cell type.'),
            ('B - Cell types:','Each row is a parent region. Positive log2FC means greater MG accessibility than that off-target type.'),
            ('C - Donors:','Values above zero mean the MG-vs-rest effect stays positive even in the weakest of the four donors.'),
            ('D - Sensitivity:','At the red FDR = 0.05 line, 872 of 1,000 parents remain significant in every glial donor-omission comparison.')],
            'The selected natural parents have donor-supported MG-enriched accessibility.',
            'They are not confirmed enhancers, and their mutants are not yet shown to retain MG specificity.'),
        box('<b>Two useful terms:</b> <b>Pearson r</b> measures prediction-observation agreement; higher is better. <b>JSD</b> measures profile difference; lower is better.'),
        PageBreak(),p('Predicted designs and remaining gaps','GuideTitle'),
        p('The key distinction: natural accessibility is observed; mutation effects and reporter performance are not.','GuideSub'),
        section('3. Do the edits have consistent model predictions?','03_variant_predictions.png',[
            ('A - Filtering:','Many nominated edits are removed. The final set contains 938 predicted gains and 946 loss controls.'),
            ('B - Orientations:','Retained gains increase and loss controls decrease predictions in both the sequence and its reverse complement.'),
            ('C - Important bases:','Loss controls usually alter highly attributed parent bases, consistent with disrupting learned positive features.'),
            ('D - Effect size:','Attribution agrees more closely with loss-effect magnitude than gain-effect magnitude. Gains are less well explained.')],
            'The retained variants have screened, consistent effects in the MG model.',
            'This is not independent validation: the same model and attribution helped select the edits. Logcount change is not measured reporter activity.'),
        section('4. What is still missing before a large order?','04_open_questions.png',[
            ('A - Priority off-targets:','Astrocytes or microglia are the strongest observed off-target for 920 of the 1,000 parents.'),
            ('B - Data depth:','Those glial types have much less ATAC data than MG, making their model training less certain.'),
            ('C - Review flags:','Repeats, low reference mappability, orientation differences and attribution failures need review. Categories overlap.'),
            ('D - Missing tests:','We still need cross-seed edit predictions, off-target variant predictions and experimental reporter testing.')],
            'The evidence supports planning a controlled pilot, not claiming that the optimized variants are already MG-specific.',
            'Observed off-target BigWigs do not substitute for off-target sequence models or mutant experiments.'),
        box('<b>Decision for Thanh:</b> Plan a small matched parent/gain/loss-control pilot. Before a large specificity-focused order, check edit stability across model seeds, assess glial off-target modeling, and agree on the assay and controls.',AMBER),
        Spacer(1,8),p('<b>Current update:</b> Observed browser tracks for 13 cell types are now available. No off-target sequence model or reporter validation was added by that export.','SmallGuide')]
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    print(DEST)

if __name__=='__main__':main()
