const fs=require('fs'),path=require('path');
const root=__dirname;process.chdir(root);
const numericPolicy=require('./numeric_policy.cjs');
fs.mkdirSync('outputs',{recursive:true});fs.mkdirSync('work',{recursive:true});fs.mkdirSync('kingdom_null',{recursive:true});
let splitDesign=null;
function run(name,files,aaOrder){
 const src=fs.readFileSync(`scripts/${name}.js`,'utf8');
 const ctx={files,numericPolicy,splitDesign,parseCsv:numericPolicy.parseCsv,
  readText:p=>{
   if(p==='POOL'){
    const raw=JSON.parse(fs.readFileSync('work/pooled-acc.json','utf8'));
    const base='ARNDCQEGHILKMFPSTWYV',aa=aaOrder||base;
    for(const c of Object.values(raw))for(const key of ['cnt','sum'])c[key]=Array.from(aa,a=>c[key][base.indexOf(a)]);
    return JSON.stringify({chains:raw,order:Object.keys(raw),aa_order:aa});
   }
   return fs.readFileSync(p,'utf8');
  },readBytes:()=>fs.readFileSync('inputs/selected_6908_organisms.csv'),
  writeFile:(p,text)=>{
   const target=(name==='pool'?'work/':'outputs/')+path.basename(p);fs.writeFileSync(target,text);
   if(name==='classify'&&path.basename(p)==='assign.json')fs.writeFileSync('kingdom_null/assign.json',text);
  }};
 new Function('ctx',src)(ctx);console.log('Completed: '+name);
}
const alpha='ACDEFGHIKLMNPQRSTVWY';
run('pool',['inputs/residue_RSA_6908_corrected.csv']);
if(process.argv.includes('--create-design')){
 const hash=numericPolicy.createSplitDesign(root,JSON.parse(fs.readFileSync('work/pooled-acc.json','utf8')));
 console.log('Frozen split design created: '+hash);
}
splitDesign=numericPolicy.loadSplitDesign(root);
run('classify',[]);
for(const name of ['decompose','s6_windows'])run(name,['POOL','outputs/assign.json'],alpha);
run('admissible_parameters',['POOL','inputs/selected_6908_organisms.csv']);
if(!process.argv.includes('--tables')){
 for(const name of ['pc3_stability','kingdom_composition','kingdom_joint_control'])run(name,['POOL','outputs/assign.json'],alpha);
 run('aaindex_ranks',['inputs/aaindex1.txt']);
 for(const name of ['parameter_correlations','length_slope'])run(name,['POOL','inputs/selected_6908_organisms.csv']);
}
