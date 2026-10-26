// parameter_correlations.js — each context-setting parameter tested alone.
// Reproduces v5b's correlation paragraph: mean Pearson r, mean Spearman rho, permutation p,
// and the held-out R^2 of each parameter fitted alone. 500 permutation replicates, as Methods.
// Inputs: [0] scratch/pool-full.json, [1] scratch/bundle/inputs/selected_6908_organisms.csv
const AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];
const CHG = { D:-1, E:-1, K:1, R:1, H:0.1 };
const KD = { A:1.8,R:-4.5,N:-3.5,D:-3.5,C:2.5,Q:-3.5,E:-3.5,G:-0.4,H:-3.2,I:4.5,L:3.8,
             K:-3.9,M:1.9,F:2.8,P:-1.6,S:-0.8,T:-0.7,W:-0.9,Y:-1.3,V:4.2 };
const RIDGE = 1e-9, NSPLIT = 20, SEED = 20260930, NPERM = 500;
const pool = JSON.parse(ctx.readText(ctx.files[0]));
const ch = pool.chains, order = pool.order;
const seq = {};
for(const f of ctx.parseCsv(ctx.readText(ctx.files[1])).slice(1)) {
  if(f.length>=3)seq[(f[0].trim()+f[1].trim()).toLowerCase()]=f[2].trim();
}
const ids = order.filter(id => ch[id].L >= 40 && seq[id]);
const P = {};
for (const id of ids) {
  const c = ch[id], s = seq[id];
  let chg=0,pos=0,neg=0;
  for (let j=0;j<20;j++){ const q=CHG[AA[j]]||0; chg+=q*c.cnt[j]; if(q>0)pos+=c.cnt[j]; else if(q<0)neg+=c.cnt[j]; }
  const h = new Float64Array(s.length);
  for (let i=0;i<s.length;i++){ const v=KD[s[i]]; h[i]=v===undefined?0:v; }
  let hMean=0; for(let i=0;i<h.length;i++) hMean+=h[i]; hMean/=h.length;
  let hMax=-Infinity, run=0; const W=11;
  for(let i=0;i<h.length;i++){ run+=h[i]; if(i>=W) run-=h[i-W]; if(i>=W-1){const m=run/W; if(m>hMax)hMax=m;} }
  let inLow=0,tot=0;
  for(let i=0;i+20<=s.length;i+=20){ let mask=0; for(let k=0;k<20;k++) mask|=1<<(s.charCodeAt(i+k)%32);
    let bits=0; for(let b=0;b<32;b++) if(mask&(1<<b)) bits++; tot+=20; if(bits<=12) inLow+=20; }
  P[id] = { length: Math.log10(c.L), charge: chg/c.L, asymmetry: (pos-neg)/c.L,
            hydroMean: hMean, hydroWindow: hMax===-Infinity?0:hMax,
            lowComplexity: tot?inLow/tot:0, own: c.cnt.map(v=>v/c.L) };
}
const TERMS = ['length','charge','asymmetry','hydroMean','hydroWindow','lowComplexity'];
const byType = AA.map((_,j)=>ids.filter(id=>ch[id].cnt[j]>0));
const resp = byType.map((l,j)=>{ const m={}; for(const id of l) m[id]=ch[id].sum[j]/ch[id].cnt[j]; return m; });
// cache y and x per term per type, once
const yBy = byType.map((l,j)=>l.map(id=>resp[j][id]));
const xBy = {}; for (const tm of TERMS) xBy[tm] = byType.map(l=>l.map(id=>P[id][tm]));

function pearson(a,b){ const n=a.length; let ma=0,mb=0;
  for(let i=0;i<n;i++){ma+=a[i];mb+=b[i];} ma/=n; mb/=n;
  let sab=0,sa=0,sb=0;
  for(let i=0;i<n;i++){ const da=a[i]-ma,db=b[i]-mb; sab+=da*db; sa+=da*da; sb+=db*db; }
  return sab/Math.sqrt(sa*sb); }
function rank(v){ const ix=v.map((x,i)=>i).sort((p,q)=>v[p]-v[q]); const r=new Array(v.length);
  for(let i=0;i<ix.length;i++){ let j=i; while(j+1<ix.length && v[ix[j+1]]===v[ix[i]]) j++;
    const avg=(i+j)/2+1; for(let k=i;k<=j;k++) r[ix[k]]=avg; i=j; } return r; }
const spearman=(a,b)=>pearson(rank(a),rank(b));
function olsFit(X,y){ const p=X[0].length,n=X.length;
  const A=Array.from({length:p},()=>new Array(p).fill(0)), b=new Array(p).fill(0);
  for(let i=0;i<n;i++){ const xi=X[i];
    for(let a=0;a<p;a++){ b[a]+=xi[a]*y[i]; for(let c=a;c<p;c++) A[a][c]+=xi[a]*xi[c]; } }
  for(let a=0;a<p;a++) for(let c=0;c<a;c++) A[a][c]=A[c][a];
  let tr=0; for(let a=0;a<p;a++) tr+=A[a][a];
  for(let a=0;a<p;a++) A[a][a]+=RIDGE*(tr/p);
  const M=A.map((r,i)=>[...r,b[i]]);
  for(let col=0;col<p;col++){
    let piv=col; for(let r=col+1;r<p;r++) if(Math.abs(M[r][col])>Math.abs(M[piv][col])) piv=r;
    if(piv!==col){const t=M[col];M[col]=M[piv];M[piv]=t;}
    const d=M[col][col]; if(Math.abs(d)<1e-14) return null;
    for(let c=col;c<=p;c++) M[col][c]/=d;
    for(let r=0;r<p;r++){ if(r===col) continue; const f=M[r][col]; if(!f) continue;
      for(let c=col;c<=p;c++) M[r][c]-=f*M[col][c]; }
  }
  return M.map(r=>r[p]); }
function r2On(beta,X,y){ const n=y.length; let yb=0; for(const v of y) yb+=v; yb/=n;
  let ssr=0,sst=0;
  for(let i=0;i<n;i++){ let yh=0; for(let a=0;a<beta.length;a++) yh+=beta[a]*X[i][a];
    ssr+=(y[i]-yh)**2; sst+=(y[i]-yb)**2; }
  return sst===0?null:1-ssr/sst; }
const avg = a => a.reduce((x,y)=>x+y,0)/a.length;

// simple correlations
const simple = {}; for (const tm of TERMS) simple[tm] = { r: [], rho: [] };
for (let j=0;j<20;j++){
  if (byType[j].length < 200) continue;
  for (const tm of TERMS) {
    simple[tm].r.push(pearson(xBy[tm][j], yBy[j]));
    simple[tm].rho.push(spearman(xBy[tm][j], yBy[j]));
  }
}
// permutation on mean |r|; one shuffle order reused across terms within a replicate
let sg = 4242 >>> 0; const rnd = () => { sg=(sg*1664525+1013904223)>>>0; return sg/4294967296; };
const obsAbs = {}; for (const tm of TERMS) obsAbs[tm] = avg(simple[tm].r.map(Math.abs));
const ge = {}; for (const tm of TERMS) ge[tm] = 0;
const permOrder = AA.map((_,j)=>j).filter(j=>byType[j].length>=200);
for (let k=0;k<NPERM;k++){
  const ord = {};
  for (const j of permOrder) {
    const m = byType[j].length, ix = new Array(m);
    for (let i=0;i<m;i++) ix[i]=i;
    for (let i=m-1;i>0;i--){ const q=Math.floor(rnd()*(i+1)); const t=ix[i]; ix[i]=ix[q]; ix[q]=t; }
    ord[j] = ix;
  }
  for (const tm of TERMS) {
    let sum=0,cnt=0;
    for (const j of permOrder) {
      const ix = ord[j], x = xBy[tm][j], y = yBy[j], m = x.length;
      const xs = new Array(m); for (let i=0;i<m;i++) xs[i]=x[ix[i]];
      sum += Math.abs(pearson(xs,y)); cnt++;
    }
    if (sum/cnt >= obsAbs[tm]) ge[tm]++;
  }
}
// held-out R^2 of each parameter alone, and of the own-fraction predictor
const SETS = {
  length: id=>[1,P[id].length], charge: id=>[1,P[id].charge], asymmetry: id=>[1,P[id].asymmetry],
  hydroMean: id=>[1,P[id].hydroMean], hydroWindow: id=>[1,P[id].hydroWindow],
  lowComplexity: id=>[1,P[id].lowComplexity]
};
const Xs = {}; for (const nm of Object.keys(SETS)) Xs[nm] = {};
for (const id of ids) for (const nm of Object.keys(SETS)) Xs[nm][id] = SETS[nm](id);
let s1 = SEED >>> 0; const r1 = () => { s1=(s1*1664525+1013904223)>>>0; return s1/4294967296; };
const held = {}; for (const nm of Object.keys(SETS)) held[nm] = [];
const heldOwn = {}; for (const t of AA) heldOwn[t] = [];
for (let sp=0; sp<NSPLIT; sp++){
  for (let j=0;j<20;j++){
    const l = byType[j]; if (l.length < 200) continue;
    const assigned = ctx.splitDesign.partition(AA[j],sp,l);
    const tr=assigned.train,te=assigned.test;
    const yTr = tr.map(id=>resp[j][id]), yTe = te.map(id=>resp[j][id]);
    for (const nm of Object.keys(SETS)) {
      const beta = olsFit(tr.map(id=>Xs[nm][id]), yTr);
      const v = beta ? r2On(beta, te.map(id=>Xs[nm][id]), yTe) : null;
      if (v!==null) held[nm].push(v);
    }
    const Xtr = tr.map(id=>[1, P[id].own[j]]), Xte = te.map(id=>[1, P[id].own[j]]);
    const bo = olsFit(Xtr,yTr);
    if (bo) { const v = r2On(bo,Xte,yTe); if (v!==null) heldOwn[AA[j]].push(v); }
  }
}
const out = {
  chains: ids.length,
  permutationReplicates: NPERM,
  permutationStatistic:"mean absolute Pearson r across the twenty amino-acid types, identical for observed and null",
  observedPermutationStatistics:obsAbs,
  correlations: Object.fromEntries(TERMS.map(tm => [tm, {
    meanPearson: +avg(simple[tm].r).toFixed(3), meanSpearman: +avg(simple[tm].rho).toFixed(3),
    permutationP: +((ge[tm]+1)/(NPERM+1)).toFixed(3),
    negativeForTypes: simple[tm].r.filter(v=>v<0).length
  }])),
  heldOutSingleParameter: Object.fromEntries(Object.keys(SETS).map(nm => [nm, +avg(held[nm]).toFixed(4)])),
  compositionOwnFraction: {
    meanHeldOutR2: +avg(Object.values(heldOwn).flat()).toFixed(4),
    perType: AA.map(t=>({ type:t, R2:+avg(heldOwn[t]).toFixed(4) })).sort((a,b)=>b.R2-a.R2)
  }
};
ctx.writeFile('scratch/v53c/results/parameter_correlations.json', JSON.stringify(out, null, 1));
return out;
