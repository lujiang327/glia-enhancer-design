import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const root=process.cwd();
const out=path.join(root,'designs/CN4461_AAV_MPRA/pilot_order');
const data=JSON.parse(await fs.readFile(path.join(out,'order_sheet_inputs.json'),'utf8'));
const wb=Workbook.create();
const order=wb.worksheets.add('Pilot order');
const review=wb.worksheets.add('Assembly review');
for(const s of [order,review]){
  s.showGridLines=false;
  s.getRange('A1:E30').format.font={name:'Arial',size:10,color:'#243447'};
  s.getRange('A1:E30').format.verticalAlignment='center';
}
order.getRange('A2').values=[['CN4461 | initial three-fragment pilot']];
order.getRange('A2').format.font={name:'Arial',size:14,bold:true};
order.getRange('A3').values=[['Selected natural parents. Existing BarcodeS unchanged. Order preparation only; no order placed.']];
order.getRange('A3').format.font={name:'Arial',size:10,italic:true};
order.getRange('A5:E5').values=[['Fragment ID','Format','Length (bp)','Complete sequence (5′ to 3′)','Cloning note']];
order.getRange('A6:E8').values=data.orders.map(r=>[r.fragment_id,'Linear dsDNA',null,r.sequence_5to3,r.vendor_note]);
for(let i=6;i<=8;i++)order.getRange(`C${i}`).formulas=[[`=LEN(D${i})`]];
order.getRange('A5:E5').format={fill:'#30475E',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:28,wrapText:true};
for(const [col,width] of [['A',29],['B',18],['C',13],['D',95],['E',45]])order.getRange(`${col}1:${col}14`).format.columnWidth=width;
order.getRange('A6:E8').format.wrapText=true;
order.getRange('A6:E8').format.rowHeight=112;
order.getRange('A6:E8').format.verticalAlignment='top';
order.getRange('D6:D8').format.font={name:'Courier New',size:10};
order.getRange('C6:C8').setNumberFormat('0');
order.getRange('C6:C8').format.horizontalAlignment='right';
order.getRange('A10').values=[['Synthesis: exact supplied sequence; no additional adapters, substitutions or new barcode.']];
order.getRange('A11').values=[['Scale: quote the vendor’s standard minimum scale before purchase.']];
order.getRange('A13').values=[['Source: Addgene CN4461 #236600, sequence 472831; exact fragment and parent identities independently rechecked.']];
order.getRange('A14').values=[['Reference: https://www.addgene.org/236600/sequences/']];

review.getRange('A2').values=[['Assembly details and evidence limits']];
review.getRange('A2').format.font={name:'Arial',size:14,bold:true};
review.getRange('A4:C4').values=[['Detail','Value','Evidence status']];
const notes=[
['Left synthesis flank (5′ to 3′)',data.arms[0].sequence_5to3+' | 20 bp | reference 125–144','Sequence verified'],
['Right synthesis flank (5′ to 3′)',data.arms[1].sequence_5to3+' | 30 bp | reference 727–756','Sequence verified'],
['Actual retained overlaps','Left: 15 bp (125–139); right: 25 bp (732–756). Site-restoration bases supplied by insert.','Cut positions verified'],
['Fragment composition','20-bp left flank + exact 500-bp parent + 30-bp right flank = 550 bp. Both flanks in forward construct direction.','Sequence verified'],
['Backbone preparation','MluI/SacI digest of backbone only; remove original enhancer-containing fragment. Do not digest the synthetic inserts.','Vendor method review'],
['Assembly','Two-fragment homology assembly; complete MluI/SacI sites restored. Right synthesis flank is not reverse-complemented.','Digital simulation passed'],
['Expected plasmid','4,738 bp; reference 145–726 replaced by 500-bp parent, all other reference bases unchanged.','Sequence verified'],
['ITR-to-ITR span','2,141 bp including long-annotated repeats; 2,130 bp using alternative right-ITR endpoint.','Endpoint annotation review'],
['Parent 140 restriction site','One internal SacI at enhancer 392–397. Compatible with separate-backbone digestion; changes final diagnostic digest.','Sequence verified'],
['Left-arm ITR homology','Six terminal left-ITR bases occur in arm; longest exact tract elsewhere is 14 bp. Physical assembly efficiency not established.','Vendor feasibility review'],
['Protected cassette','All promoter/linker, SYFP2, WPRE/linker, barcode, polyA and ITR reference bases outside enhancer retained.','Sequence verified'],
['Named-feature limits','Exact minBG endpoints, WPRE3 terminal convention, BarcodeS name mapping and right-ITR endpoint remain partly inferred.','Annotation review'],
['Existing barcode','AGCAATTA retained identically in all three. Keep constructs separate for initial cloning QC.','No unique pooled IDs'],
['Laboratory validation','Verify complete insert, junctions, full plasmid, barcode region and ITR integrity after cloning.','Not yet performed'],
];
review.getRange(`A5:C${4+notes.length}`).values=notes;
review.getRange('A4:C4').format={fill:'#30475E',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:27};
review.getRange('A1:A28').format.columnWidth=32;
review.getRange('B1:B28').format.columnWidth=100;
review.getRange('C1:C28').format.columnWidth=30;
review.getRange(`A5:C${4+notes.length}`).format.wrapText=true;
review.getRange(`A5:C${4+notes.length}`).format.rowHeight=43;
review.getRange('A20:C20').values=[['Fragment ID','SHA256 of supplied 550-bp fragment','Nearest gene']];
review.getRange('A20:C20').format={fill:'#E6ECF1',font:{name:'Arial',size:10,bold:true},rowHeight:24};
review.getRange('A21:C23').values=data.orders.map(r=>[r.fragment_id,r.sequence_sha256,r.nearest_gene]);
review.getRange('A21:C23').format.rowHeight=25;
review.getRange('B21:B23').format.font={name:'Courier New',size:9};
review.getRange('A25').values=[['Reference: https://www.addgene.org/236600/sequences/ (sequence 472831).']];
review.getRange('A26').values=[['Enhancer-boundary evidence: Hunker et al., Neuron 2025, Table S1 row 30; DOI 10.1016/j.neuron.2025.04.035.']];

wb.recalculate();
const inspection=await wb.inspect({kind:'table',range:'Pilot order!A5:E8',include:'values,formulas',tableMaxRows:4,tableMaxCols:5,tableMaxCellChars:65,maxChars:2000});
await fs.writeFile(path.join(out,'spreadsheet_inspection.ndjson'),inspection.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Final formula error scan'});
await fs.writeFile(path.join(out,'spreadsheet_formula_scan.ndjson'),errors.ndjson);
for(const [sheet,range,name] of [['Pilot order','A1:E14','pilot_order_preview.png'],['Assembly review','A1:C26','assembly_review_preview.png']]){
 const preview=await wb.render({sheetName:sheet,range,scale:1.5,format:'png'});
 await fs.writeFile(path.join(out,name),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(out,'CN4461_three_construct_pilot_order.xlsx'));
console.log('Saved three-fragment pilot workbook and two sheet previews.');
