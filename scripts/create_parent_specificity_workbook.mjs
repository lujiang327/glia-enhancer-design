// Run from a temporary directory linked to the bundled artifact-tool node_modules.
import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const [inputDir, outputDir] = process.argv.slice(2);
if (!inputDir || !outputDir) throw new Error('Usage: builder.mjs INPUT_DIR OUTPUT_DIR');
const data = JSON.parse(await fs.readFile(path.join(inputDir, 'workbook_inputs.json'), 'utf8'));
const wb = Workbook.create();
const sheets = Object.fromEntries(['Review','Parents','log2FC','FDR','Browser signal'].map(n=>[n,wb.worksheets.add(n)]));
const color = '#234B68';
const column = n => {let s=''; for(n++;n;n=Math.floor((n-1)/26)) s=String.fromCharCode(65+(n-1)%26)+s; return s;};
function base(sheet, end, rows, title, note) {
  sheet.showGridLines=false;
  sheet.getRange(`A1:${end}${rows}`).format.font={name:'Arial',size:10,color:'#24313C'};
  sheet.getRange(`A1:${end}${rows}`).format.rowHeight=20;
  sheet.getRange(`A1:${end}${rows}`).format.verticalAlignment='center';
  sheet.getRange(`A1:${end}${rows}`).format.columnWidth=17;
  sheet.getRange('A2').values=[[title]];
  sheet.getRange('A2').format.font={name:'Arial',size:14,bold:true,color:'#24313C'};
  if(note) {sheet.getRange('A3').values=[[note]];sheet.getRange('A3').format.font={name:'Arial',size:10,italic:true,color:'#56616B'};}
}
function table(sheet, headers, rows, name) {
  const end=column(headers.length-1), bottom=5+rows.length;
  sheet.getRange(`A5:${end}${bottom}`).values=[headers,...rows];
  sheet.tables.add(`A5:${end}${bottom}`,true,name);
  const h=sheet.getRange(`A5:${end}5`);
  h.format={fill:color,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:48};
  sheet.getRange(`A6:${end}${bottom}`).format.rowHeight=20;
  sheet.freezePanes.freezeRows(5);
  sheet.freezePanes.freezeColumns(2);
}
const rowMap = Object.fromEntries(data.tables.parents.columns.map((c,i)=>[c,i]));
const val=(row,key)=>row[rowMap[key]];
// Main review table: statistical values are retained at full precision.
const parentRows=data.tables.parents.rows.map(r=>[
  val(r,'candidate_rank'),val(r,'peak_id'),val(r,'IGV_parent_locus'),val(r,'genomic_class'),val(r,'nearest_gene'),
  val(r,'primary_logFC'),val(r,'primary_FDR'),null,null,null,null,
  val(r,'browser_strongest_offtarget'),val(r,'maximum_glial_lodo_FDR'),val(r,'primary_donors_positive'),
  val(r,'primary_minimum_donor_logFC'),val(r,'browser_donors_MG_above_both_glia_count')]);
base(sheets.Parents,'P',1005,'Review the 1,000 natural parent regions','Sorted by original rank. IGV coordinates are 1-based inclusive. Browser measurements describe natural loci.');
table(sheets.Parents,['Rank','Peak ID','IGV parent locus (hg38)','Genomic class','Nearest gene','MG vs rest log2FC','MG vs rest FDR','Minimum pairwise log2FC','Maximum pairwise FDR','Significant positive pairs (of 12)','Browser MG / max off-target','Strongest browser off-target','Worst glial LODO FDR','Positive primary donors (of 4)','Minimum donor primary log2FC','Browser donors MG > both glia (of 4)'],parentRows,'ParentsReview');
for(let i=0;i<1000;i++) {
  const r=i+6;
  sheets.Parents.getRange(`H${r}:K${r}`).formulas=[[
    `=MIN('log2FC'!H${r}:S${r})`, `=MAX('FDR'!H${r}:S${r})`,
    `=COUNTIFS('log2FC'!H${r}:S${r},">0",'FDR'!H${r}:S${r},"<0.05")`,
    `='Browser signal'!G${r}/MAX('Browser signal'!H${r}:S${r})`]];
}
sheets.Parents.getRange('A5:A1005').format.columnWidth=8;
sheets.Parents.getRange('B5:B1005').format.columnWidth=23;
sheets.Parents.getRange('C5:C1005').format.columnWidth=36;
sheets.Parents.getRange('D5:D1005').format.columnWidth=22;
sheets.Parents.getRange('F6:F1005').setNumberFormat('0.00');
for(const c of ['G','I','M']) sheets.Parents.getRange(`${c}6:${c}1005`).setNumberFormat('0.00E+00');
for(const c of ['H','K','O']) sheets.Parents.getRange(`${c}6:${c}1005`).setNumberFormat('0.00');
for(const c of ['J','N','P']) sheets.Parents.getRange(`${c}6:${c}1005`).setNumberFormat('0');
sheets.Parents.getRange('M6:M1005').conditionalFormats.add('cellIs',{operator:'greaterThanOrEqual',formula:.05,format:{fill:'#FCE6D5',font:{color:'#9C3C20'}}});
for(const [key,name,title,note] of [
  ['log2FC','log2FC','MG versus each comparator: donor-blocked log2 fold change','Positive values favor MG. 1 = twofold, 2 = fourfold, 3 = eightfold fitted enrichment.'],
  ['FDR','FDR','MG versus each comparator: original BH-adjusted FDR','Original tests across 147,359 peaks per contrast. FDR < 0.05 is the significance threshold.'],
  ['browser','Browser signal','Observed mean browser signal in each exact 500bp parent','IPM per base, 100bp mean bins. Global pooled normalization; not TMM or absolute accessibility per cell.']]) {
  const t=data.tables[key],sheet=sheets[name];
  base(sheet,'S',1005,title,note);
  const headers=['Rank','Peak ID','Chromosome','Start (0-based)','End (exclusive)','IGV parent locus (hg38)',...(key==='browser'?['MG',...data.comparisons.slice(1)]:['Rest',...data.comparisons.slice(1)])];
  table(sheet,headers,t.rows,name.replaceAll(' ','')+'Data');
  sheet.getRange('A5:A1005').format.columnWidth=8;
  sheet.getRange('B5:B1005').format.columnWidth=23;
  sheet.getRange('C5:C1005').format.columnWidth=13;
  sheet.getRange('D6:E1005').setNumberFormat('0');
  sheet.getRange('F5:F1005').format.columnWidth=36;
  sheet.getRange('G6:S1005').setNumberFormat(key==='FDR'?'0.00E+00':key==='browser'?'0.0000':'0.00');
  if(key==='log2FC') sheet.getRange('G6:S1005').conditionalFormats.add('colorScale',{colors:['#FFFFFF','#CCE1EF','#2878A5'],thresholds:[0,4,8]});
  if(key==='FDR') sheet.getRange('G6:S1005').conditionalFormats.add('cellIs',{operator:'greaterThanOrEqual',formula:.05,format:{fill:'#FCE6D5',font:{color:'#9C3C20'}}});
}
const review=sheets.Review;
base(review,'G',48,'MG parent differential-accessibility review','Four donors. Existing selected-region results, reviewed October 20, 2026.');
review.tabColor=color;
const summaryRows=data.comparisons.map(ct=>[ct==='Rest'?'Rest (equal type weights)':ct,null,null,null,null,null,null]);
table(review,['Comparator','Significant positive (of 1,000)','Minimum log2FC','Median log2FC','Median fitted fold enrichment','Median FDR','Maximum FDR'],summaryRows,'ComparisonSummary');
review.freezePanes.unfreeze();
review.getRange('A5:A18').format.columnWidth=30;
review.getRange('B5:G18').format.columnWidth=20;
for(let i=0;i<data.comparisons.length;i++) {
  const r=i+6,c=column(i+6),fc=`'log2FC'!${c}6:${c}1005`,fd=`'FDR'!${c}6:${c}1005`;
  review.getRange(`B${r}:G${r}`).formulas=[[
    `=COUNTIFS(${fc},">0",${fd},"<0.05")`,`=MIN(${fc})`,`=MEDIAN(${fc})`,`=2^D${r}`,`=MEDIAN(${fd})`,`=MAX(${fd})`]];
}
review.getRange('B6:B18').setNumberFormat('#,##0');
review.getRange('C6:E18').setNumberFormat('0.00');
review.getRange('F6:G18').setNumberFormat('0.00E+00');
const metrics=[
  ['Parents significant and positive in all 12 pairs',`=COUNTIFS('Parents'!J6:J1005,12)`],
  ['Minimum pooled browser MG / strongest off-target',`=MIN('Parents'!K6:K1005)`],
  ['Median pooled browser MG / strongest off-target',`=MEDIAN('Parents'!K6:K1005)`],
  ['Parents passing all eight glial LODO FDR tests',`=COUNTIFS('Parents'!M6:M1005,"<0.05")`],
  ['Parents with MG > both glia in all four donor tracks',`=COUNTIFS('Parents'!P6:P1005,4)`],
];
for(let i=0;i<metrics.length;i++) {
  const r=21+i;review.getRange(`A${r}`).values=[[metrics[i][0]]];review.getRange(`D${r}`).formulas=[[metrics[i][1]]];
  review.getRange(`D${r}`).setNumberFormat(i===1||i===2?'0.00':'#,##0');
}
const notes=[
  'How to interpret',
  'Enrichment allows off-target peaks. The closest glial contrasts are weaker than the neuronal contrasts.',
  'FDR measures differential evidence. It does not establish MG exclusivity, reporter activity or variant specificity.',
  'These 1,000 regions were selected using specificity metrics. This review is not independent validation.',
  'LODO = leave one donor out. Positive effects persist, but 128 parents lose FDR < 0.05 in at least one glial refit.',
  'In IGV use hg38, Group Autoscale, and a few-kb view. Compare the exact 500bp parent interval.',
  'The browser signal counts fragment endpoints. DA counts fragments overlapping each consensus interval.',
  'Pooled browser tracks weight donors by depth. DA uses biological donor blocking and sample-specific TMM factors.',
  'Sources and methods',
  'DA source: top_1000_candidates.tsv.gz, original Phase 3 candidate-ranking export (all 1,000 rows).',
  'Method: edgeR robust quasi-likelihood F tests, donor-blocked design, TMM normalization, BH FDR per contrast.',
  'Browser source: the delivered retina_hg38_review_100bp package. All 25 BigWig SHA256 hashes verified.',
  'Browser means were read at exact base values over each parent, including partial 100bp bins; no zoom summaries.',
  'Pseudobulks contain four contributing donors. Low glial depth limits inference from missing browser peaks.',
  'Detailed donor-track values and SHA256 provenance accompany the review in the audit report.'
];
notes.forEach((n,i)=>{const r=28+i;review.getRange(`A${r}`).values=[[n]];if(i===0||i===8)review.getRange(`A${r}`).format.font={name:'Arial',size:10,bold:true,color};});
// Allocate note text across otherwise unused cells without merging table columns.
review.getRange('A28:G42').format.wrapText=false;
await fs.mkdir(outputDir,{recursive:true});
wb.recalculate();
console.log((await wb.inspect({kind:'region',sheetId:'Review',range:'A6:G8',maxChars:1800,tableMaxRows:3,tableMaxCols:7})).ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Formula error scan'});
console.log(errors.ndjson);
// Check representative calculated values before export.
const expected=data.summary;
const close=(a,b)=>{if(Math.abs(a-b)>1e-10*Math.max(1,Math.abs(b)))throw new Error(`Mismatch ${a} vs ${b}`);};
close(review.getRange('D21').values[0][0],expected.all_twelve_positive_and_FDR_lt_0_05);
close(review.getRange('D22').values[0][0],expected.browser_ratio_min);
close(review.getRange('D23').values[0][0],expected.browser_ratio_median);
close(review.getRange('D24').values[0][0],expected.all_eight_glial_lodo_FDR_lt_0_05);
close(review.getRange('D25').values[0][0],expected.all_four_browser_donors_MG_above_both_glia);
for(const [name,range] of [['Review','A1:G25'],['Parents','A1:J11'],['log2FC','A1:J11'],['FDR','A1:J11'],['Browser signal','A1:J11']]) {
  const preview=await wb.render({sheetName:name,range,scale:1.3,format:'png'});
  await fs.writeFile(path.join(outputDir,name.replaceAll(' ','_')+'_preview.png'),new Uint8Array(await preview.arrayBuffer()));
}
const output=await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(outputDir,'MG_1000_parent_specificity.xlsx'));
console.log('Saved workbook with 1,000 parents, 13 original comparisons, and observed browser measurements.');
