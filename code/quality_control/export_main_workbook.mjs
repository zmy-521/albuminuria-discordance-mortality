import fs from 'node:fs/promises';
import {Workbook,SpreadsheetFile} from '@oai/artifact-tool';
if (!process.env.ODKD_WORK_ROOT) throw new Error('Set ODKD_WORK_ROOT');
const S=process.env.ODKD_WORK_ROOT.replaceAll(String.fromCharCode(92),'/')+'/Stage3_finalization';
const ps=JSON.parse(await fs.readFile(S+'/internal/workbook_payload.json','utf8'));
const wb=Workbook.create();
function col(n){let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26);}return s;}
for(const p of ps){
 const s=wb.worksheets.add(p.name);s.showGridLines=false;const n=p.rows.length,c=p.headers.length;
 s.getCell(0,0).values=[[p.title]];s.getCell(1,0).values=[[p.note]];
 s.getRangeByIndexes(3,0,1,c).values=[p.headers];s.getRangeByIndexes(4,0,n,c).values=p.rows;
 const all=s.getRangeByIndexes(0,0,n+5,Math.max(c,7));all.format.font={name:'Times New Roman',size:11,color:'#000000'};all.format.rowHeight=25;all.format.verticalAlignment='center';
 s.getCell(0,0).format.font={name:'Times New Roman',size:14,bold:true,color:'#000000'};
 const h=s.getRangeByIndexes(3,0,1,c);h.format.font={name:'Times New Roman',size:11,bold:true,color:'#000000'};h.format.rowHeight=48;h.format.wrapText=true;h.format.borders={bottom:{style:'thin',color:'#000000'}};
 for(let j=0;j<c;j++){
  const vals=p.rows.map(r=>r[j]);const textLen=Math.max(p.headers[j].length,...vals.filter(x=>typeof x==='string').map(x=>x.length));
  let width=p.headers[j]==='URL'?100:Math.min(44,Math.max(15,textLen*.85+5));
  s.getRangeByIndexes(0,j,n+5,1).format.columnWidth=width;
  const rg=s.getRangeByIndexes(4,j,n,1);rg.format.wrapText=true;
  if(vals.every(x=>x===null||typeof x==='number'))rg.setNumberFormat(vals.every(x=>x===null||Number.isInteger(x))?'#,##0':p.headers[j]==='p'?'0.00E+00':'0.00000');
 }
 s.getRangeByIndexes(4,0,n,c).format.autofitRows();
 const note=s.getRangeByIndexes(1,0,1,Math.max(c,7));note.merge();note.format.wrapText=true;note.format.rowHeight=42;
 if(p.name==='Denominators')s.getRangeByIndexes(4,0,n,c).format.rowHeight=65;
 if(p.name==='Independent QC')s.getRangeByIndexes(4,0,n,c).format.rowHeight=45;
 s.freezePanes.freezeRows(4);
}
wb.recalculate();
await fs.writeFile(S+'/internal/workbook_errors.json',(await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!',options:{useRegex:true,maxResults:50},summary:'final error scan'})).ndjson);
await fs.writeFile(S+'/internal/workbook_inspection.json',(await wb.inspect({kind:'table',range:"'Stage3 Cox'!A4:I8",include:'values,formulas',tableMaxRows:5,tableMaxCols:11})).ndjson);
for(const p of ps){const blob=await wb.render({sheetName:p.name,range:`A1:${col(Math.max(p.headers.length,7))}${Math.min(p.rows.length+4,10)}`,scale:1,format:'png'});await fs.writeFile(S+'/internal/preview_'+p.name.replaceAll(' ','_')+'.png',new Uint8Array(await blob.arrayBuffer()));}
const file=await SpreadsheetFile.exportXlsx(wb);await file.save(S+'/19_final_clinical_results.xlsx');console.log('Export complete: 25 aggregate sheets');

