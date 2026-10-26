const fs=require('fs'),path=require('path');const root=__dirname;process.chdir(root);
fs.mkdirSync('outputs',{recursive:true});fs.mkdirSync('work',{recursive:true});
function run(name,files,aaOrder){
 const src=fs.readFileSync(`scripts/${name}.js`,'utf8');
 const ctx={files,readText:p=>{
  if(p==='POOL'){
   const raw=JSON.parse(fs.readFileSync('work/pooled-acc.json','utf8'));
   const base='ARNDCQEGHILKMFPSTWYV';const aa=aaOrder||base;
   for(const c of Object.values(raw))for(const key of ['cnt','sum'])c[key]=Array.from(aa,a=>c[key][base.indexOf(a)]);
   return JSON.stringify({chains:raw,order:Object.keys(raw),aa_order:aa});
  }
  return fs.readFileSync(p,'utf8');},
 readBytes:p=>fs.readFileSync('inputs/selected_6908_organisms.csv'),
 writeFile:(p,t)=>{const dest=(name==='pool'?'work/':'outputs/')+path.basename(p);fs.writeFileSync(dest,t);if(name==='classify'&&path.basename(p)==='assign.json')fs.writeFileSync('kingdom_null/assign.json',t);}}
 const result=new Function('ctx',src)(ctx);console.log(name,JSON.stringify(result));
}
const alpha='ACDEFGHIKLMNPQRSTVWY';run('pool',['inputs/residue_RSA_6908_corrected.csv']);run('classify',[]);
for(const name of ['decompose','pc3_stability','s6_windows','kingdom_composition','kingdom_joint_control'])run(name,['POOL','outputs/assign.json'],alpha);
run('aaindex_ranks',['inputs/aaindex1.txt']);
for(const name of ['admissible_parameters','parameter_correlations','length_slope'])run(name,['POOL','inputs/selected_6908_organisms.csv']);
