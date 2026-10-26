// kingdom_joint_control.js — length and composition controlled together.
// Reproduces v5b: "restricting the kingdoms to a common length range reduces the shift to 0.018...
// Matching on both simultaneously leaves 0.020."
// Design: 3 composition tertiles x 3 length tertiles = 9 joint cells, cut on the pooled
// distribution; ΔPC1 recomputed inside each cell; pooled by the smaller kingdom in each cell.
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
           nPositive:dPC1.filter(v=>v>0).length };
}
const euk0=order.filter(id=>asg[id]&&asg[id].k==='euk'&&ch[id].L>=40);
const prok0=order.filter(id=>asg[id]&&asg[id].k==='prok'&&ch[id].L>=40);
const baseline=kingdomShift(euk0,prok0);

// composition PC1
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

// tertile cut points on the pooled distribution
const pooled=[...euk0,...prok0];
const cv = pooled.map(id=>cscore[id]).sort((a,b)=>a-b);
const lv = pooled.map(id=>ch[id].L).sort((a,b)=>a-b);
const ccut=[cv[Math.floor(cv.length/3)], cv[Math.floor(2*cv.length/3)]];
const lcut=[lv[Math.floor(lv.length/3)], lv[Math.floor(2*lv.length/3)]];
const cbin=id=>cscore[id]<ccut[0]?0:(cscore[id]<ccut[1]?1:2);
const lbin=id=>ch[id].L<lcut[0]?0:(ch[id].L<lcut[1]?1:2);

// length-only reference (one length control)
const lengthOnly=[];
for (let b=0;b<3;b++){
  const e=euk0.filter(id=>lbin(id)===b), p=prok0.filter(id=>lbin(id)===b);
  const s=kingdomShift(e,p);
  if (s) lengthOnly.push({ tertile:b+1, nEuk:e.length, nProk:p.length, mean:s.mean, t:s.t });
}
const loW = lengthOnly.reduce((a,o)=>a+Math.min(o.nEuk,o.nProk),0);
const loPooled = lengthOnly.reduce((a,o)=>a+o.mean*Math.min(o.nEuk,o.nProk),0)/loW;

// joint cells
const cells=[];
for (let ci=0;ci<3;ci++) for (let li=0;li<3;li++){
  const e=euk0.filter(id=>cbin(id)===ci&&lbin(id)===li);
  const p=prok0.filter(id=>cbin(id)===ci&&lbin(id)===li);
  const s=kingdomShift(e,p);
  cells.push({ cell:'C'+(ci+1)+'L'+(li+1), compTertile:ci+1, lenTertile:li+1,
    nEuk:e.length, nProk:p.length, mean:s?s.mean:null, t:s?s.t:null, usable:!!s });
}
const usable=cells.filter(c=>c.usable);
const wsum=usable.reduce((a,c)=>a+Math.min(c.nEuk,c.nProk),0);
const jointPooled=usable.reduce((a,c)=>a+c.mean*Math.min(c.nEuk,c.nProk),0)/wsum;
const out = {
  chains: { euk: euk0.length, prok: prok0.length },
  baseline,
  lengthOnly: { tertiles: lengthOnly, pooled: +loPooled.toFixed(5),
                ratioToBaseline: +(loPooled/baseline.mean).toFixed(3) },
  joint: { cells, usableCells: usable.length, pooled: +jointPooled.toFixed(5),
           ratioToBaseline: +(jointPooled/baseline.mean).toFixed(3),
           tRange: [Math.min(...usable.map(c=>c.t)), Math.max(...usable.map(c=>c.t))] }
};
ctx.writeFile('scratch/v53c/results/kingdom_joint_control.json', JSON.stringify(out, null, 1));
return { baseline: out.baseline, lengthOnly: out.lengthOnly.pooled, joint: out.joint.pooled,
         ratios: { len: out.lengthOnly.ratioToBaseline, joint: out.joint.ratioToBaseline },
         reducedTo: { 'length 0.018': +(loPooled).toFixed(3), 'joint 0.020': +(jointPooled).toFixed(3) } };
