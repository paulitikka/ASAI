// kingdom_composition.js — is ΔPC1 a compositional effect?
// Reproduces v5b: "within composition tertiles the shift is 0.029 against 0.031 unmatched".
// Method: chains binned into tertiles of the first principal component of the twenty composition
// fractions, cut at the tertile boundaries of the pooled distribution; ΔPC1 recomputed inside
// each tertile. A surviving ΔPC1 implies a non-compositional origin.
// Inputs: [0] scratch/pool-full.json, [1] scratch/kg/assign.json
const AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
const Q = 99;
const pool = JSON.parse(ctx.readText(ctx.files[0]));
const asgRaw = JSON.parse(ctx.readText(ctx.files[1]));
const asg = {}; for (const k in asgRaw) asg[k.toLowerCase()] = asgRaw[k];
const ch = pool.chains, order = pool.order;

function quant(v){
  if(!v.length) return new Array(Q).fill(0);
  const a=v.slice().sort((x,y)=>x-y), o=[];
  for(let q=1;q<=Q;q++){ const pos=(q/(Q+1))*(a.length-1), lo=Math.floor(pos), hi=Math.min(lo+1,a.length-1);
    o.push(a[lo]+(a[hi]-a[lo])*(pos-lo)); }
  return o;
}
function jacobi(C,n){
  let A=C.map(r=>r.slice()), V=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  for(let sw=0;sw<800;sw++){
    let off=0; for(let i=0;i<n;i++) for(let j=i+1;j<n;j++) off+=A[i][j]*A[i][j];
    if(off<1e-24) break;
    for(let p=0;p<n;p++) for(let q=p+1;q<n;q++){
      if(Math.abs(A[p][q])<1e-16) continue;
      const th=0.5*Math.atan2(2*A[p][q],A[q][q]-A[p][p]), c=Math.cos(th), s=Math.sin(th);
      for(let i=0;i<n;i++){const a=A[i][p],b=A[i][q];A[i][p]=c*a-s*b;A[i][q]=s*a+c*b;}
      for(let i=0;i<n;i++){const a=A[p][i],b=A[q][i];A[p][i]=c*a-s*b;A[q][i]=s*a+c*b;}
      for(let i=0;i<n;i++){const a=V[i][p],b=V[i][q];V[i][p]=c*a-s*b;V[i][q]=s*a+c*b;}
    }
  }
  const ev=Array.from({length:n},(_,i)=>Math.max(A[i][i],0));
  const ord=Array.from({length:n},(_,i)=>i).sort((x,y)=>ev[y]-ev[x]);
  return {V,ord};
}
const matrix=(ids,j)=>quant(ids.filter(id=>ch[id].cnt[j]>0).map(id=>ch[id].sum[j]/ch[id].cnt[j]));
function decomp(rowsMat,n){
  const cm=new Array(Q).fill(0);
  for(let q=0;q<Q;q++){let t=0;for(let i=0;i<n;i++)t+=rowsMat[i][q];cm[q]=t/n;}
  const C=Array.from({length:n},()=>new Array(n).fill(0));
  for(let i=0;i<n;i++) for(let j=0;j<n;j++){let t=0;for(let q=0;q<Q;q++)t+=(rowsMat[i][q]-cm[q])*(rowsMat[j][q]-cm[q]);C[i][j]=t/Q;}
  return jacobi(C,n);
}
function kingdomShift(euk, prok){
  if (euk.length < 60 || prok.length < 60) return null;
  const rows=[];
  for(let j=0;j<20;j++) rows.push(matrix(euk,j));
  for(let j=0;j<20;j++) rows.push(matrix(prok,j));
  const d=decomp(rows,40);
  let p1=Array.from({length:40},(_,i)=>d.V[i][d.ord[0]]);
  const gm=AA.map((_,j)=>{let t=0,n=0;for(const id of [...euk,...prok]){const c=ch[id];t+=c.sum[j];n+=c.cnt[j];}return t/n;});
  { let a=0; for(let i=0;i<20;i++) a += p1[i]*gm[i] + p1[i+20]*gm[i]; if(a<0) p1=p1.map(x=>-x); }
  const dPC1=AA.map((_,j)=>p1[j]-p1[j+20]);
  const mean=dPC1.reduce((a,b)=>a+b,0)/20;
  const sd=Math.sqrt(dPC1.reduce((a,b)=>a+(b-mean)**2,0)/19);
  return { mean:+mean.toFixed(5), t:+(mean/(sd/Math.sqrt(20))).toFixed(2),
           nPositive:dPC1.filter(v=>v>0).length, nEuk:euk.length, nProk:prok.length };
}
const euk0=order.filter(id=>asg[id]&&asg[id].k==='euk'&&ch[id].L>=40);
const prok0=order.filter(id=>asg[id]&&asg[id].k==='prok'&&ch[id].L>=40);
const baseline=kingdomShift(euk0,prok0);

// composition PC1 over the pooled chains
const allIds=order.filter(id=>ch[id].L>=40);
const frac={};
for (const id of allIds){ const c=ch[id]; frac[id]=c.cnt.map(v=>v/c.L); }
const cm=new Array(20).fill(0);
for (const id of allIds) for (let j=0;j<20;j++) cm[j]+=frac[id][j];
for (let j=0;j<20;j++) cm[j]/=allIds.length;
const Cov=Array.from({length:20},()=>new Array(20).fill(0));
for (const id of allIds) for (let i=0;i<20;i++) for (let j=0;j<20;j++) Cov[i][j]+=(frac[id][i]-cm[i])*(frac[id][j]-cm[j]);
for (let i=0;i<20;i++) for (let j=0;j<20;j++) Cov[i][j]/=allIds.length;
const dcp=jacobi(Cov,20);
const cpc1=Array.from({length:20},(_,i)=>dcp.V[i][dcp.ord[0]]);
const cscore={};
for (const id of allIds){ let s=0; for(let j=0;j<20;j++) s+=(frac[id][j]-cm[j])*cpc1[j]; cscore[id]=s; }

const pooled=[...euk0,...prok0].filter(id=>cscore[id]!==undefined);
const vals=pooled.map(id=>cscore[id]).sort((a,b)=>a-b);
const cut=[vals[Math.floor(vals.length/3)], vals[Math.floor(2*vals.length/3)]];
const cbin=id=>cscore[id]<cut[0]?0:(cscore[id]<cut[1]?1:2);
const tertiles=[];
for (let b=0;b<3;b++){
  const e=euk0.filter(id=>cbin(id)===b), p=prok0.filter(id=>cbin(id)===b);
  const s=kingdomShift(e,p);
  tertiles.push({ tertile:b+1, nEuk:e.length, nProk:p.length,
    mean:s?s.mean:null, t:s?s.t:null, nPositive:s?s.nPositive:null });
}
const usable=tertiles.filter(o=>o.mean!==null);
const wsum=usable.reduce((a,o)=>a+Math.min(o.nEuk,o.nProk),0);
const pooledMean=usable.reduce((a,o)=>a+o.mean*Math.min(o.nEuk,o.nProk),0)/wsum;
const out = {
  chains: { euk: euk0.length, prok: prok0.length },
  baseline,
  compositionTertiles: tertiles,
  pooledTertileMean: +pooledMean.toFixed(5),
  ratioToBaseline: +(pooledMean/baseline.mean).toFixed(3),
  note: 'bins are tertiles of the first principal component of the twenty composition fractions, cut on the pooled distribution'
};
ctx.writeFile('scratch/v53c/results/kingdom_composition.json', JSON.stringify(out, null, 1));
return { baseline: out.baseline, pooled: out.pooledTertileMean, ratio: out.ratioToBaseline, tertiles };
