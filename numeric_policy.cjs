const fs=require('fs'),path=require('path'),crypto=require('crypto');
const intervalTolerance=1e-12;
function rsaInterval(v){
  if(!Number.isFinite(v)||v<0)throw new Error('Invalid RSA interval value');
  if(v===0)return 0;
  const nearest=Math.round(v*10);
  const adjusted=nearest>=1&&nearest<=9&&Math.abs(v-nearest/10)<=intervalTolerance?nearest/10:v;
  for(let i=1;i<10;i++)if(adjusted<=i/10)return i;
  return 10;
}
function parseCsv(text){
  const rows=[];let row=[],field='',quoted=false;
  for(let i=0;i<text.length;i++){
    const c=text[i];
    if(c==='"'){
      if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;
    }else if(c===','&&!quoted){row.push(field);field='';}
    else if((c==='\n'||c==='\r')&&!quoted){
      if(c==='\r'&&text[i+1]==='\n')i++;
      row.push(field);if(row.some(x=>x!==''))rows.push(row);row=[];field='';
    }else field+=c;
  }
  if(quoted)throw new Error('Unclosed CSV quote');
  if(field!==''||row.length){row.push(field);rows.push(row);}
  return rows;
}
function eligibleChains(root,pool){
  const rows=parseCsv(fs.readFileSync(path.join(root,'inputs/selected_6908_organisms.csv'),'utf8'));
  const sequences=new Map(rows.slice(1).map(r=>[(r[0]+r[1]).toLowerCase(),r[2]]));
  return Object.keys(pool).filter(id=>pool[id].L>=40&&sequences.get(id));
}
function createSplitDesign(root,pool){
  const aa='ARNDCQEGHILKMFPSTWYV',ids=eligibleChains(root,pool);
  const byType=Array.from(aa,(_,j)=>ids.filter(id=>pool[id].cnt[j]>0));
  const ranks=new Map();
  for(let j=0;j<20;j++)for(const id of byType[j])ranks.set(id+'|'+aa[j],new Array(20).fill(0));
  let seed=20260930>>>0;
  const rnd=()=>{seed=(seed*1664525+1013904223)>>>0;return seed/4294967296;};
  for(let sp=0;sp<20;sp++)for(let j=0;j<20;j++){
    const shuffled=byType[j].slice();
    for(let i=shuffled.length-1;i>0;i--){const q=Math.floor(rnd()*(i+1));[shuffled[i],shuffled[q]]=[shuffled[q],shuffled[i]];}
    shuffled.forEach((id,k)=>ranks.get(id+'|'+aa[j])[sp]=k+1);
  }
  const header=['chain_key','aa',...Array.from({length:20},(_,i)=>'permutation_rank_'+String(i+1).padStart(2,'0'))];
  const lines=[header.join(',')];
  for(const id of ids)for(const t of aa){const r=ranks.get(id+'|'+t);if(r)lines.push([id,t,...r].join(','));}
  fs.mkdirSync(path.join(root,'design'),{recursive:true});
  const target=path.join(root,'design/train_test_splits.csv');fs.writeFileSync(target,lines.join('\n')+'\n');
  const hash=crypto.createHash('sha256').update(fs.readFileSync(target)).digest('hex');
  fs.writeFileSync(path.join(root,'design/design_metadata.json'),JSON.stringify({seed:20260930,replicates:20,
    aminoAcidLoopOrder:aa,chainOrder:'first appearance in the frozen residue CSV',
    generator:'LCG: state = (1664525 * state + 1013904223) mod 2^32; Fisher-Yates shuffle',
    trainingRule:'rank <= floor(number of eligible chains containing the amino-acid type / 2)',
    testRule:'remaining ranks; original permutation order retained within both halves',
    nEligibleChains:ids.length,nRows:ranks.size,sha256:hash},null,2)+'\n');
  return hash;
}
function loadSplitDesign(root){
  const file=path.join(root,'design/train_test_splits.csv');
  if(!fs.existsSync(file))throw new Error('Missing frozen split design; restore the matching data archive.');
  const rows=parseCsv(fs.readFileSync(file,'utf8'));
  if(rows[0].length!==22)throw new Error('Unexpected split design header');
  const data=new Map();
  for(const row of rows.slice(1)){
    if(row.length!==22||data.has(row[0]+'|'+row[1]))throw new Error('Invalid or duplicate split design row');
    data.set(row[0]+'|'+row[1],row.slice(2).map(Number));
  }
  const cache=new Map();
  return {partition(aa,sp,ids){
    const key=aa+'|'+sp;
    if(!cache.has(key)){
      const sorted=ids.slice().sort((a,b)=>data.get(a+'|'+aa)?.[sp]-data.get(b+'|'+aa)?.[sp]);
      sorted.forEach((id,i)=>{if(data.get(id+'|'+aa)?.[sp]!==i+1)throw new Error('Split design does not match eligible chains: '+aa);});
      const half=Math.floor(sorted.length/2);cache.set(key,{train:sorted.slice(0,half),test:sorted.slice(half)});
    }
    const saved=cache.get(key);if(saved.train.length+saved.test.length!==ids.length)throw new Error('Inconsistent split eligibility');
    return saved;
  }};
}
module.exports={intervalTolerance,rsaInterval,parseCsv,createSplitDesign,loadSplitDesign};
