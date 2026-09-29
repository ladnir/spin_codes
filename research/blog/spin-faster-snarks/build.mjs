import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {createRequire} from 'node:module';

// Use the bundled runtime; no package installation is necessary.
const req=createRequire('C:/Users/peter/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/runtime-anchor.cjs');
const sharp=req('sharp');
const {marked}=req('marked');
const here=path.dirname(fileURLToPath(import.meta.url));
const figures=path.join(here,'figures');
await fs.mkdir(figures,{recursive:true});
const data=JSON.parse(await fs.readFile(path.join(here,'performance.json'),'utf8'));
const C={bg:'#FAF9F6',ink:'#202326',muted:'#5E656B',line:'#D7D8D4',zero:'#E8E8E2',orange:'#CD461A',pale:'#FADED1',purple:'#7259A6',green:'#267660',gray:'#606B75'};
const esc=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;').replaceAll('"','&quot;');
let parts=[];
function start(title,desc,h){parts=[`<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="${h}" viewBox="0 0 1100 ${h}" role="img" aria-labelledby="title desc"><title id="title">${esc(title)}</title><desc id="desc">${esc(desc)}</desc><defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0L8 4L0 8Z" fill="${C.orange}"/></marker></defs><rect width="1100" height="${h}" fill="${C.bg}"/>`];}
function text(x,y,s,size=23,color=C.ink,anchor='start',weight=400){parts.push(`<text x="${x}" y="${y}" font-family="Arial, sans-serif" font-size="${size}" font-weight="${weight}" fill="${color}" text-anchor="${anchor}">${esc(s)}</text>`);}
function lines(x,y,ss,size=23,color=C.ink,anchor='start',gap=size*1.35,weight=400){ss.forEach((s,i)=>text(x,y+i*gap,s,size,color,anchor,weight));}
function rect(x,y,w,h,fill='none',stroke='none',sw=1,rx=0){parts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${rx}" fill="${fill}" stroke="${stroke}" stroke-width="${sw}"/>`);}
function line(x1,y1,x2,y2,color=C.line,width=1,arrow=false,dash=''){parts.push(`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="${width}" ${arrow?'marker-end="url(#arrow)"':''} ${dash?`stroke-dasharray="${dash}"`:''}/>`);}
function arrow(x1,y1,x2,y2){line(x1,y1,x2,y2,C.orange,2,true);}
function cell(x,y,v,w=25,h=31,active=false){rect(x,y,w,h,active?(v?C.orange:C.pale):C.zero,active?C.orange:C.bg,1.4,2);if(w>=18)text(x+w/2,y+h/2+7,v,22,active&&v?'#fff':(active?C.orange:C.muted),'middle',active&&v?700:400);}
async function save(name){const svg=parts.join('\n')+'</svg>';await fs.writeFile(path.join(figures,name+'.svg'),svg);await sharp(Buffer.from(svg),{density:180}).png().toFile(path.join(figures,name+'.png'));}

start('Where encoding sits in a prover','The witness enters an encoding and hashing pipeline. Computation checks produce claims, and the PCS provides authenticated openings.',480);
text(42,53,'Where the prover spends work',30,C.ink,'start',600);
text(42,86,'Conceptual roles in a code-based proof system',19,C.muted);
rect(292,128,755,172,'none',C.line,1.5,8);
text(316,155,'COMMITMENT LAYER',16,C.muted,'start',600);
lines(123,218,['Large witness','table'],24,C.ink,'middle');
arrow(224,222,314,222);
rect(331,181,180,82,C.pale,C.orange,1.5,5);
text(421,216,'Encode',28,C.orange,'middle',600);
text(421,245,'touch every row',18,C.orange,'middle');
arrow(521,222,598,222);
lines(697,212,['Hash encoded','columns'],24,C.ink,'middle');
arrow(792,222,868,222);
text(953,231,'Commitment',22,C.ink,'middle',600);
lines(168,382,['Computation checks','reduce correctness to claims'],22,C.ink,'middle');
arrow(365,384,486,384);
lines(664,372,['Claims about the witness','answered by the PCS'],22,C.ink,'middle');
arrow(859,384,935,384);
lines(1000,374,['Sample +','authenticate'],19,C.ink,'middle');
text(43,449,'Fast algebraic checks still need an efficient commitment layer.',22,C.orange);
await save('01-prover-pipeline');

start('Distance makes sampling useful','Two illustrative pairs of words differ in two and eight positions out of 32. Four fixed sample locations miss the first and detect the second.',475);
text(42,53,'More disagreement, more chances to catch it',30,C.ink,'start',600);
text(42,86,'The same four sampled positions are outlined in green.',19,C.muted);
const a=Array.from({length:32},(_,i)=>[0,1,0,1,1,0,0,1][i%8]);
const small=[7,22],large=[1,5,9,13,17,21,25,29],queries=[2,9,18,25];
for(const [top,diff,label] of [[149,small,'2 / 32 differ'],[316,large,'8 / 32 differ']]){
 text(42,top+22,'Word A',20);text(42,top+65,'Word B',20);
 for(let j=0;j<32;j++){const hit=diff.includes(j);cell(160+j*27,top,a[j],25,32,hit);cell(160+j*27,top+43,a[j]^(hit?1:0),25,32,hit);if(queries.includes(j))rect(158+j*27,top-4,29,83,'none',C.green,2.5,4);}
 text(160,top+105,label,21,C.orange,'start',600);
 text(1024,top+105,diff===small?'All four samples miss.':'Two samples detect disagreement.',20,C.muted,'end');
}
await save('02-distance-and-sampling');

// A genuine small linear encoder example: repeated [8,4,4] outer blocks,
// a fixed permutation, then a zero-initialized recursive inner with memory 6.
const generators=['11110000','11001100','10101010','00001111'];
const toyWords=Array.from({length:16},(_,m)=>Array.from({length:8},(_,j)=>generators.reduce((v,g,k)=>v^(((m>>k)&1)*Number(g[j])),0)));
if(Math.min(...toyWords.slice(1).map(w=>w.reduce((a,b)=>a+b,0)))!==4)throw Error('Toy outer distance mismatch');
const input=Array(32).fill(0);input[4]=1;
const outer=Array(64).fill(0);for(let j=0;j<8;j++)outer[8+j]=Number(generators[0][j]);
let rngstate=137;function rand(){rngstate=(Math.imul(rngstate,1664525)+1013904223)>>>0;return rngstate/4294967296;}
const perm=Array.from({length:64},(_,i)=>i);for(let j=63;j>0;j--){const k=Math.floor(rand()*(j+1));[perm[j],perm[k]]=[perm[k],perm[j]];}
const shuffled=perm.map(i=>outer[i]);
const out=[],masks=[];let state=0;for(let i=0;i<64;i++){const mask=Math.floor(rand()*64);masks.push(mask);let z=state&mask,parity=0;while(z){parity^=z&1;z>>>=1;}const bit=shuffled[i]^parity;out.push(bit);state=((state<<1)|bit)&63;}
start('A sparse input through a toy SPIN encoder','One nonzero input bit becomes four outer ones, which are permuted and passed through a recursive encoder with six bits of memory.',610);
text(42,53,'Following one nonzero bit',30,C.ink,'start',600);
text(42,86,'A small execution of the outer–interleaver–inner architecture',19,C.muted);
const rows=[{y:143,label:'Sparse input',bits:input,sub:'32 input bits'},{y:260,label:'Block encoding',bits:outer,sub:'64 intermediate bits'},{y:377,label:'Interleaving',bits:shuffled,sub:'The same four ones'},{y:494,label:'Recursive mixing',bits:out,sub:'Influence persists'}];
rows.forEach((r,ri)=>{
 text(42,r.y+20,r.label,22,C.ink,'start',600);text(42,r.y+48,r.sub,17,C.muted);
 const step=12; r.bits.forEach((v,j)=>cell(282+j*step,r.y,v,10,36,v===1));
 if(ri<3){arrow(612,r.y+49,612,r.y+99);text(638,r.y+81,['Create several ones','Scatter activation points','Mix with the running state'][ri],19,C.muted);}
});
text(282,566,'One illustrative execution; distance is a guarantee over all nonzero inputs.',19,C.muted);
await save('03-sparse-input');
await fs.writeFile(path.join(figures,'sparse-example.json'),JSON.stringify({description:'Toy illustration, not production parameters or evidence of global minimum distance.',outerGenerators:generators,outerMinimumDistance:4,input,outer,permutationOutputToInput:perm,shuffled,memory:6,feedbackMasks:masks,output:out},null,2));

start('Structured routing spreads each block across regions','Four six-bit outer blocks are transposed into six four-bit regions. Three active coordinates from a single block land in three distinct regions.',475);
text(42,53,'Structured routing',30,C.ink,'start',600);
text(42,86,'One active outer block, tracked through the permutation',19,C.muted);
text(195,142,'Shuffle each block',23,C.ink,'middle',600);
text(553,142,'Transpose',23,C.ink,'middle',600);
text(907,142,'Shuffle each region',23,C.ink,'middle',600);
const cw=33,ch=33;const positions=[2,0,3,1,0,2];
for(let r=0;r<4;r++)for(let c=0;c<6;c++)cell(96+c*34,195+r*34,r===1&&c%2===0?1:0,cw,ch,r===1);
for(let r=0;r<6;r++)for(let c=0;c<4;c++){
cell(486+c*34,162+r*34,c===1&&r%2===0?1:0,cw,ch,c===1);
cell(842+c*34,162+r*34,c===positions[r]&&r%2===0?1:0,cw,ch,c===positions[r]);
}
arrow(329,263,449,263);arrow(650,263,787,263);
text(198,362,'One block of weight 3',20,C.orange,'middle');
text(910,396,'3 distinct active regions',20,C.orange,'middle',600);
text(550,445,'Tiled transpose + shuffles within smaller regions',24,C.ink,'middle');
await save('04-structured-routing');

const r=data.flock['14'];
const series=[['Ligerito',r['native-fast100']],['SPIN–Brakedown',r.spin]];
const breakdown=series.map(([name,x])=>({name,commit:x.commit_ms,open:x.open_ms,other:x.total_ms-x.commit_ms-x.open_ms,total:x.total_ms}));
start('FLOCK prover time at 16,384 BLAKE3 compressions',`Measured stacked bars separate commitment, opening, and all remaining prover work. Ligerito totals ${r['native-fast100'].total_ms.toFixed(2)} milliseconds and SPIN–Brakedown ${r.spin.total_ms.toFixed(2)} milliseconds.`,500);
text(42,53,`${(r['native-fast100'].total_ms/r.spin.total_ms).toFixed(2)}× prover throughput`,30,C.ink,'start',600);
text(42,87,'16,384 BLAKE3 compressions · one Ryzen 9 7950X core',19,C.muted);
const x0=241,scale=3.38,ybase=375;
for(let tick=0;tick<=200;tick+=50){const x=x0+tick*scale;line(x,155,x,ybase,C.line,1);text(x,ybase+30,tick,18,C.muted,'middle');}
text(x0+scale*200+5,ybase+61,'Prover time (ms)',18,C.muted,'end');
breakdown.forEach((s,i)=>{const y=183+i*120;let x=x0;text(218,y+32,s.name,21,C.ink,'end',600);
for(const [v,color] of [[s.commit,C.orange],[s.open,C.purple],[s.other,C.gray]]){rect(x,y,v*scale,55,color);text(x+v*scale/2,y+34,v.toFixed(1),18,'#fff','middle',600);x+=v*scale;}
text(x+15,y+33,Math.round(s.total)+' ms',23,C.ink,'start',600);
});
[['Commitment',C.orange,42],['Opening',C.purple,235],['Remaining prover work',C.gray,411]].forEach(([label,color,x])=>{rect(x,458,18,18,color);text(x+28,473,label,18,C.ink);});
await save('05-flock-performance');
await fs.writeFile(path.join(figures,'flock-breakdown.json'),JSON.stringify({source:'performance.json',workload:16384,units:'milliseconds',remainingDefinition:'total_ms - commit_ms - open_ms',breakdown},null,2));

const md=await fs.readFile(path.join(here,'post.md'),'utf8');
const body=marked.parse(md);
const html=`<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>SPIN: Faster Codes for Faster SNARKs</title><style>
:root{color-scheme:light}*{box-sizing:border-box}body{margin:0;background:#faf9f6;color:#25282b;font-family:Georgia,'Times New Roman',serif;font-size:20px;line-height:1.67}article{max-width:820px;margin:80px auto 100px;padding:0 26px}h1,h2,th{font-family:Arial,sans-serif}h1{font-size:52px;line-height:1.12;letter-spacing:-1.6px;margin:0 0 42px;font-weight:700}h2{font-size:30px;line-height:1.25;letter-spacing:-.4px;margin:65px 0 22px}p{margin:0 0 24px}a{color:#ad3817;text-decoration-thickness:1px;text-underline-offset:3px}strong{font-weight:700}p:has(>img){margin:36px -110px 12px}img{display:block;width:100%;height:auto;border:1px solid #e0dfd9;border-radius:4px}p:has(>img)+p{font-family:Arial,sans-serif;font-size:15px;line-height:1.55;color:#606569;margin:0 -40px 36px}p:has(>img)+p em{font-style:normal}table{border-collapse:collapse;width:100%;font-family:Arial,sans-serif;font-size:15px;margin:32px 0}th{font-weight:600;border-bottom:2px solid #cd461a}td,th{padding:12px 8px;vertical-align:top}td{border-bottom:1px solid #deded8}thead{color:#35393d}code{font-size:.83em}footer{font-family:Arial,sans-serif;font-size:14px;color:#727773;margin-top:65px;border-top:1px solid #d8d8d0;padding-top:20px}@media(max-width:1060px){p:has(>img){margin-left:0;margin-right:0}p:has(>img)+p{margin-left:0;margin-right:0}}@media(max-width:600px){body{font-size:18px}article{margin-top:36px;padding:0 20px}h1{font-size:38px}h2{font-size:27px}table{font-size:12px}td,th{padding:9px 5px}}
</style></head><body><article>${body}</article></body></html>`;
await fs.writeFile(path.join(here,'preview.html'),html);
const count=md.replace(/!\[[^\]]*\]\([^)]*\)/g,'').replace(/\]\([^)]*\)/g,']').split(/\s+/).length;
console.log(JSON.stringify({wordCount:count,figures:5,draft:path.join(here,'post.md'),preview:path.join(here,'preview.html'),breakdown},null,2));
