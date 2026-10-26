// s6_windows.js — Supplementary Table S6: the kingdom comparison under length restriction.
// DETERMINISTIC: a filter on chain length, no sampling, so no seed is involved.
// Designs: (a) all chains, (b) 150 <= L < 350, (c) 200 <= L < 300.
// For each design the eukaryotic and prokaryotic quantile matrices are stacked into a joint
// 40 x 99 matrix, decomposed by Jacobi, and PC1/PC2 oriented as in Table 2. The per-type
// difference between kingdoms gives a paired t across the twenty types.
const AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
const Q = 99;
const pool = JSON.parse(ctx.readText(ctx.files[0]));
const asgRaw = JSON.parse(ctx.readText(ctx.files[1]));
const asg = {}; for (const k in asgRaw) asg[k.toLowerCase()] = asgRaw[k];
const ch = pool.chains, order = pool.order;

function quant(v){
  if (!v.length) return new Array(Q).fill(0);
  const s = v.slice().sort((a,b)=>a-b), o = [];
  for (let q=1;q<=Q;q++){
    const pos=(q/(Q+1))*(s.length-1), lo=Math.floor(pos), hi=Math.min(lo+1,s.length-1);
    o.push(s[lo]+(s[hi]-s[lo])*(pos-lo));
  }
  return o;
}
function jacobi(C, n){
  let A=C.map(r=>r.slice()), V=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  for (let sw=0; sw<1000; sw++){
    let off=0; for(let i=0;i<n;i++) for(let j=i+1;j<n;j++) off+=A[i][j]*A[i][j];
    if (off<1e-24) break;
    for (let p=0;p<n;p++) for (let q=p+1;q<n;q++){
      if (Math.abs(A[p][q])<1e-16) continue;
      const th=0.5*Math.atan2(2*A[p][q],A[q][q]-A[p][p]), c=Math.cos(th), s=Math.sin(th);
      for (let i=0;i<n;i++){const a=A[i][p],b=A[i][q];A[i][p]=c*a-s*b;A[i][q]=s*a+c*b;}
      for (let i=0;i<n;i++){const a=A[p][i],b=A[q][i];A[p][i]=c*a-s*b;A[q][i]=s*a+c*b;}
      for (let i=0;i<n;i++){const a=V[i][p],b=V[i][q];V[i][p]=c*a-s*b;V[i][q]=s*a+c*b;}
    }
  }
  const ev=Array.from({length:n},(_,i)=>Math.max(A[i][i],0));
  const ord=Array.from({length:n},(_,i)=>i).sort((x,y)=>ev[y]-ev[x]);
  const tot=ev.reduce((a,b)=>a+b,0);
  return { V, ord, shares: ord.map(i=>+(100*ev[i]/tot).toFixed(2)) };
}
const matrix = (ids,j) => quant(ids.filter(id=>ch[id].cnt[j]>0).map(id=>ch[id].sum[j]/ch[id].cnt[j]));
const mean = a => a.reduce((x,y)=>x+y,0)/a.length;
const sd = a => { const m=mean(a); return Math.sqrt(a.reduce((s,v)=>s+(v-m)*(v-m),0)/(a.length-1)); };
const meanLen = g => g.reduce((s,id)=>s+ch[id].L,0)/g.length;

function run(label, gE, gP){
  const rows = [];
  for (let j=0;j<20;j++) rows.push(matrix(gE,j));
  for (let j=0;j<20;j++) rows.push(matrix(gP,j));
  const cm = new Array(Q).fill(0);
  for (let q=0;q<Q;q++){ let s=0; for (let i=0;i<40;i++) s+=rows[i][q]; cm[q]=s/40; }
  const C = Array.from({length:40},()=>new Array(40).fill(0));
  for (let i=0;i<40;i++) for (let j=0;j<40;j++){
    let s=0; for (let q=0;q<Q;q++) s+=(rows[i][q]-cm[q])*(rows[j][q]-cm[q]);
    C[i][j]=s/Q;
  }
  const d = jacobi(C, 40);
  let p1 = Array.from({length:40},(_,i)=>d.V[i][d.ord[0]]);
  // orient towards higher accessibility
  const gm = AA.map((_,j)=>{ let s=0,n=0; for (const id of [...gE,...gP]){ const c=ch[id]; s+=c.sum[j]; n+=c.cnt[j]; } return s/n; });
  { let acc=0; for (let i=0;i<20;i++) acc += p1[i]*gm[i] + p1[i+20]*gm[i]; if (acc<0) p1=p1.map(x=>-x); }
  let p2 = Array.from({length:40},(_,i)=>d.V[i][d.ord[1]]);
  const iK=AA.indexOf('K'), iG=AA.indexOf('G');
  { const k1=(p2[iK]+p2[iK+20])/2, g1=(p2[iG]+p2[iG+20])/2; if (k1<g1) p2=p2.map(x=>-x); }
  const d1 = AA.map((_,i)=>p1[i]-p1[i+20]);
  const d2 = AA.map((_,i)=>p2[i]-p2[i+20]);
  const tOf = a => mean(a)/(sd(a)/Math.sqrt(a.length));
  return { label, nEuk:gE.length, nProk:gP.length,
    meanLenEuk:+meanLen(gE).toFixed(1), meanLenProk:+meanLen(gP).toFixed(1),
    gap:+(meanLen(gP)-meanLen(gE)).toFixed(1),
    meanLenEukFullPrecision:meanLen(gE),meanLenProkFullPrecision:meanLen(gP),
    pc1_tFullPrecision:tOf(d1),pc2_tFullPrecision:tOf(d2),
    deltaPC1FullPrecision:d1,deltaPC2FullPrecision:d2,
    shares: d.shares.slice(0,3),
    pc1_t:+tOf(d1).toFixed(2), pc1_pos:d1.filter(v=>v>0).length,
    pc2_t:+tOf(d2).toFixed(2), pc2_pos:d2.filter(v=>v>0).length };
}
const win = (g,lo,hi) => g.filter(id => { const L=ch[id].L; return L>=lo && L<hi; });
const euk = order.filter(id=>asg[id.toLowerCase()] && asg[id.toLowerCase()].k==='euk');
const prok = order.filter(id=>asg[id.toLowerCase()] && asg[id.toLowerCase()].k==='prok');

const out = [
  run('All chains', euk, prok),
  run('150 <= L < 350', win(euk,150,350), win(prok,150,350)),
  run('200 <= L < 300', win(euk,200,300), win(prok,200,300)),
];
ctx.writeFile('scratch/bundle/results/s6_windows.json', JSON.stringify(out,null,1));
return out;