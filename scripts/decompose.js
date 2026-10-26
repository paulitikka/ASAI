// decompose.js — per-chain accumulators, Tables 1/2/S3 and Supplementary Table S5.
// Written to match the Methods exactly; run against the hashed inputs.
//
// Table 1 (Component analysis): per-chain per-type mean RSA; per-type distribution of those
//   per-chain means summarised by 1st-99th percentile quantiles at pos = q/(Q+1)*(n-1), Q = 99,
//   forming a 20 x 99 matrix; centred by quantile; covariance over Q; Jacobi decomposition.
//   PC1 oriented towards higher accessibility; PC2 so that lysine lies above glycine.
//   The mean RSA column is the unweighted average of the per-chain means over all chains
//   (no length filter, no weighting by residue count); verified against Table 1 to <= 0.0001.
// Table 2: the eukaryotic and prokaryotic quantile matrices stacked into a single 40 x 99 matrix,
//   each amino acid contributing one row per kingdom.
// S3: weighted least squares of per-chain mean RSA on log10 chain length, weighted by that
//   type's residue count, chains of 40 residues or more.
// S5 (Divergence across the accessibility range): each chain-type pair — one per residue type
//   present in a chain — is assigned the mean RSA of that type in that chain; the pairs are
//   binned into eleven intervals with exact zeros separate; in each interval the per-type
//   fraction of that type's pairs falling in it is computed, and the CV of those twenty
//   fractions taken across the types. Counts are numbers of chain-type pairs, 120,822 in total.
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
    const pos = (q/(Q+1))*(s.length-1), lo = Math.floor(pos), hi = Math.min(lo+1, s.length-1);
    o.push(s[lo] + (s[hi]-s[lo])*(pos-lo));
  }
  return o;
}
function jacobi(C, n){
  let A = C.map(r=>r.slice()), V = Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  for (let sw=0; sw<1000; sw++){
    let off=0; for (let i=0;i<n;i++) for (let j=i+1;j<n;j++) off += A[i][j]*A[i][j];
    if (off < 1e-24) break;
    for (let p=0;p<n;p++) for (let q=p+1;q<n;q++){
      if (Math.abs(A[p][q]) < 1e-16) continue;
      const th = 0.5*Math.atan2(2*A[p][q], A[q][q]-A[p][p]), c=Math.cos(th), s=Math.sin(th);
      for (let i=0;i<n;i++){ const a=A[i][p], b=A[i][q]; A[i][p]=c*a-s*b; A[i][q]=s*a+c*b; }
      for (let i=0;i<n;i++){ const a=A[p][i], b=A[q][i]; A[p][i]=c*a-s*b; A[q][i]=s*a+c*b; }
      for (let i=0;i<n;i++){ const a=V[i][p], b=V[i][q]; V[i][p]=c*a-s*b; V[i][q]=s*a+c*b; }
    }
  }
  const ev = Array.from({length:n},(_,i)=>Math.max(A[i][i],0));
  const ord = Array.from({length:n},(_,i)=>i).sort((x,y)=>ev[y]-ev[x]);
  const tot = ev.reduce((a,b)=>a+b,0);
  return { V, ord, shares: ord.map(i=>100*ev[i]/tot) };
}
const matrix = (ids, j) => quant(ids.filter(id=>ch[id].cnt[j]>0).map(id=>ch[id].sum[j]/ch[id].cnt[j]));
function decompFromRows(rowsMat, n){
  const cm = new Array(Q).fill(0);
  for (let q=0;q<Q;q++){ let s=0; for (let i=0;i<n;i++) s+=rowsMat[i][q]; cm[q]=s/n; }
  const C = Array.from({length:n},()=>new Array(n).fill(0));
  for (let i=0;i<n;i++) for (let j=0;j<n;j++){
    let s=0; for (let q=0;q<Q;q++) s += (rowsMat[i][q]-cm[q])*(rowsMat[j][q]-cm[q]);
    C[i][j] = s/Q;
  }
  return jacobi(C, n);
}

// ---------- Table 1 ----------
const allIds = order.slice();
const d1 = decompFromRows(AA.map((_,j)=>matrix(allIds, j)), 20);
let pc1 = Array.from({length:20},(_,i)=>d1.V[i][d1.ord[0]]);
// mean RSA column: UNWEIGHTED AVERAGE of the per-chain means, over all chains
const meanRSA = AA.map((_,j)=>{ const v=[]; for (const id of allIds){ const c=ch[id]; if (c.cnt[j]>0) v.push(c.sum[j]/c.cnt[j]); } return v.reduce((a,b)=>a+b,0)/v.length; });
{ let acc=0; for (let i=0;i<20;i++) acc += pc1[i]*meanRSA[i]; if (acc<0) pc1=pc1.map(x=>-x); }
let pc2 = Array.from({length:20},(_,i)=>d1.V[i][d1.ord[1]]);
const iK = AA.indexOf('K'), iG = AA.indexOf('G');
if (pc2[iK] < pc2[iG]) pc2 = pc2.map(x=>-x);
const pc3 = Array.from({length:20},(_,i)=>d1.V[i][d1.ord[2]]);
if (pc3[AA.indexOf('M')] < 0) for (let i=0;i<20;i++) pc3[i]*=-1;
const table1FullPrecision = AA.map((t,j)=>({t,meanRSA:meanRSA[j],PC1:pc1[j],PC2:pc2[j],PC3:pc3[j]}));
const table1 = AA.map((t,j)=>({ t, meanRSA:+meanRSA[j].toFixed(4), PC1:+pc1[j].toFixed(4), PC2:+pc2[j].toFixed(4), PC3:+pc3[j].toFixed(4) }));
const shares1 = d1.shares.slice(0,4).map(v=>+v.toFixed(2));

// ---------- Table 2 ----------
const euk = order.filter(id=>asg[id.toLowerCase()] && asg[id.toLowerCase()].k==='euk');
const prok = order.filter(id=>asg[id.toLowerCase()] && asg[id.toLowerCase()].k==='prok');
const m2 = []; for (let j=0;j<20;j++) m2.push(matrix(euk,j)); for (let j=0;j<20;j++) m2.push(matrix(prok,j));
const d2 = decompFromRows(m2, 40);
let p1a = Array.from({length:40},(_,i)=>d2.V[i][d2.ord[0]]);
const gm = AA.map((_,j)=>{ let s=0,n=0; for (const id of [...euk,...prok]){ const c=ch[id]; s+=c.sum[j]; n+=c.cnt[j]; } return s/n; });
{ let acc=0; for (let i=0;i<20;i++) acc += p1a[i]*gm[i] + p1a[i+20]*gm[i]; if (acc<0) p1a=p1a.map(x=>-x); }
let p2a = Array.from({length:40},(_,i)=>d2.V[i][d2.ord[1]]);
{ const k1=(p2a[iK]+p2a[iK+20])/2, g1=(p2a[iG]+p2a[iG+20])/2; if (k1<g1) p2a=p2a.map(x=>-x); }
const table2FullPrecision = AA.map((t,j)=>({t,PC1_euk:p1a[j],PC1_prok:p1a[j+20],dPC1:p1a[j]-p1a[j+20],PC2_euk:p2a[j],PC2_prok:p2a[j+20],dPC2:p2a[j]-p2a[j+20]}));
const table2 = AA.map((t,j)=>({ t, PC1_euk:+p1a[j].toFixed(4), PC1_prok:+p1a[j+20].toFixed(4),
  dPC1:+(p1a[j]-p1a[j+20]).toFixed(4), PC2_euk:+p2a[j].toFixed(4), PC2_prok:+p2a[j+20].toFixed(4),
  dPC2:+(p2a[j]-p2a[j+20]).toFixed(4) }));
const shares2 = d2.shares.slice(0,3).map(v=>+v.toFixed(2));

// ---------- S3 ----------
const ids40 = order.filter(id=>ch[id].L >= 40);
function wls(xs, ys, ws){
  let S=0,Sx=0,Sy=0,Sxx=0,Sxy=0;
  for (let i=0;i<xs.length;i++){ const w=ws[i], x=xs[i], y=ys[i];
    S+=w; Sx+=w*x; Sy+=w*y; Sxx+=w*x*x; Sxy+=w*x*y; }
  const den = S*Sxx - Sx*Sx;
  const b = (S*Sxy - Sx*Sy)/den, a = (Sy - b*Sx)/S;
  let ss=0; for (let i=0;i<xs.length;i++){ const w=ws[i], x=xs[i], y=ys[i]; const r=y-(a+b*x); ss+=w*r*r; }
  return { b, se: Math.sqrt((ss/(xs.length-2))*S/den), n: xs.length };
}
const S3 = [];
for (let j=0;j<20;j++){
  const xs=[], ys=[], ws=[];
  for (const id of ids40){ const c=ch[id]; if (c.cnt[j]===0) continue;
    xs.push(Math.log10(c.L)); ys.push(c.sum[j]/c.cnt[j]); ws.push(c.cnt[j]); }
  const f = wls(xs, ys, ws);
  S3.push({ t: AA[j], exponent:+f.b.toFixed(3), se:+f.se.toFixed(3), chains:f.n });
}
S3.sort((a,b)=>a.exponent-b.exponent);

// ---------- S5: chain-type pairs ----------
const IV = [{n:'exactly 0',kind:'zero'},{n:'0-10',lo:0,hi:0.10},{n:'10-20',lo:0.10,hi:0.20},
  {n:'20-30',lo:0.20,hi:0.30},{n:'30-40',lo:0.30,hi:0.40},{n:'40-50',lo:0.40,hi:0.50},
  {n:'50-60',lo:0.50,hi:0.60},{n:'60-70',lo:0.60,hi:0.70},{n:'70-80',lo:0.70,hi:0.80},
  {n:'80-90',lo:0.80,hi:0.90},{n:'90+',lo:0.90,hi:99}];
const cnt = IV.map(()=>new Array(20).fill(0));
const perTypeTotal = new Array(20).fill(0);
for (const id of order){
  const c = ch[id];
  for (let j=0;j<20;j++){
    if (c.cnt[j] === 0) continue;
    const pr = c.sum[j]/c.cnt[j];
    perTypeTotal[j]++;
    cnt[ctx.numericPolicy.rsaInterval(pr)][j]++;

  }
}
const totalPairs = perTypeTotal.reduce((a,b)=>a+b,0);
const S5FullPrecision = IV.map((iv,k) => {
  const frac = cnt[k].map((v,j)=> perTypeTotal[j] ? v/perTypeTotal[j] : 0);
  const m = frac.reduce((a,b)=>a+b,0)/20;
  const sd = Math.sqrt(frac.reduce((s,v)=>s+(v-m)*(v-m),0)/19);
  const n = cnt[k].reduce((a,b)=>a+b,0);
  return {interval:iv.n,pairs:n,pct:100*n/totalPairs,cv:m?sd/m:0};
});

const S5=S5FullPrecision.map(r=>({...r,pct:+r.pct.toFixed(2),cv:+r.cv.toFixed(3)}));
const out = { shares1, table1, shares2, table2, S3, S5, totalPairs,
  table1FullPrecision,table2FullPrecision,S5FullPrecision,intervalCountsByType:cnt,chainsWithType:perTypeTotal,aminoAcidOrder:AA,
  shares1FullPrecision:d1.shares,shares2FullPrecision:d2.shares,
  intervalPolicy:{tolerance:ctx.numericPolicy.intervalTolerance,closure:'(lower, upper]; exact zeros separate; last >0.9'},
  kingdomCounts: { euk: euk.length, prok: prok.length } };
ctx.writeFile('scratch/bundle/results/decompose.json', JSON.stringify(out,null,1));
return { shares1, shares2, totalPairs, S5,
  meanRSA_sample: AA.slice(0,5).map((a,j)=>a+':'+meanRSA[j].toFixed(4)),
  S3_lowest: S3[0], S3_highest: S3[S3.length-1] };
