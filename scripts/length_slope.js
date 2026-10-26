// length_slope.js — the slope of per-type mean RSA on log10 chain length.
// Reproduces v5b: "The mean slope on log10 chain length is -0.14 +- 0.04 across types."
// Reports the pooled slope with its standard error, and the mean and SD of the twenty per-type
// slopes, so the +/- in the paper can be read as the between-type spread.
// Inputs: [0] scratch/pool-full.json
const AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];
const pool = JSON.parse(ctx.readText(ctx.files[0]));
const ch = pool.chains, order = pool.order;
const ids = order.filter(id => ch[id].L >= 40);

function olsWithSE(X, y) {
  const p = X[0].length, n = X.length;
  const A = Array.from({length:p},()=>new Array(p).fill(0)), b = new Array(p).fill(0);
  for (let i=0;i<n;i++){ const xi=X[i];
    for (let a=0;a<p;a++){ b[a]+=xi[a]*y[i]; for (let c=a;c<p;c++) A[a][c]+=xi[a]*xi[c]; } }
  for (let a=0;a<p;a++) for (let c=0;c<a;c++) A[a][c]=A[c][a];
  const Ad = A.map(r=>r.slice());
  const M = A.map((r,i)=>[...r,b[i]]);
  for (let col=0;col<p;col++){
    let piv=col; for (let r=col+1;r<p;r++) if (Math.abs(M[r][col])>Math.abs(M[piv][col])) piv=r;
    if (piv!==col){ const t=M[col]; M[col]=M[piv]; M[piv]=t; }
    const d=M[col][col]; if (Math.abs(d)<1e-14) return null;
    for (let c=col;c<=p;c++) M[col][c]/=d;
    for (let r=0;r<p;r++){ if (r===col) continue; const f=M[r][col]; if (!f) continue;
      for (let c=col;c<=p;c++) M[r][c]-=f*M[col][c]; }
  }
  const beta = M.map(r=>r[p]);
  const G = Ad.map(r=>r.slice());
  const inv = Array.from({length:p},(_,i)=>Array.from({length:p},(_,j)=>i===j?1:0));
  for (let col=0;col<p;col++){
    let piv=col; for (let r=col+1;r<p;r++) if (Math.abs(G[r][col])>Math.abs(G[piv][col])) piv=r;
    if (piv!==col){ const t=G[col];G[col]=G[piv];G[piv]=t; const t2=inv[col];inv[col]=inv[piv];inv[piv]=t2; }
    const d=G[col][col]; if (Math.abs(d)<1e-14) return null;
    for (let c=0;c<p;c++){ G[col][c]/=d; inv[col][c]/=d; }
    for (let r=0;r<p;r++){ if (r===col) continue; const f=G[r][col];
      for (let c=0;c<p;c++){ G[r][c]-=f*G[col][c]; inv[r][c]-=f*inv[col][c]; } }
  }
  const yb = y.reduce((a,b2)=>a+b2,0)/n;
  let ssr=0,sst=0;
  for (let i=0;i<n;i++){ let yh=0; for (let a=0;a<p;a++) yh+=beta[a]*X[i][a];
    ssr+=(y[i]-yh)**2; sst+=(y[i]-yb)**2; }
  const sigma2 = ssr/(n-p);
  const se = beta.map((_,i)=>Math.sqrt(sigma2*inv[i][i]));
  return { beta, se, R2: 1-ssr/sst, n };
}
const perType = [];
for (let j=0;j<20;j++){
  const l = ids.filter(id => ch[id].cnt[j] > 0);
  if (l.length < 200) continue;
  const X = l.map(id => [1, Math.log10(ch[id].L)]);
  const y = l.map(id => ch[id].sum[j]/ch[id].cnt[j]);
  const f = olsWithSE(X, y);
  if (!f) continue;
  perType.push({ type: AA[j], n: l.length, slope: +f.beta[1].toFixed(4), se: +f.se[1].toFixed(4),
                 R2: +f.R2.toFixed(4) });
}
// pooled over all chain-type observations
const Xp=[], yp=[];
for (let j=0;j<20;j++){
  const l = ids.filter(id => ch[id].cnt[j] > 0);
  for (const id of l){ Xp.push([1, Math.log10(ch[id].L)]); yp.push(ch[id].sum[j]/ch[id].cnt[j]); }
}
const pooled = olsWithSE(Xp, yp);
const slopes = perType.map(o=>o.slope);
const m = slopes.reduce((a,b)=>a+b,0)/slopes.length;
const sdv = Math.sqrt(slopes.reduce((a,b)=>a+(b-m)**2,0)/(slopes.length-1));
const out = {
  chains: ids.length,
  pooled: { slope: +pooled.beta[1].toFixed(4), se: +pooled.se[1].toFixed(4),
            R2: +pooled.R2.toFixed(4), n: pooled.n },
  acrossTypes: { meanSlope: +m.toFixed(4), sdOfSlopes: +sdv.toFixed(4),
                 min: +Math.min(...slopes).toFixed(4), max: +Math.max(...slopes).toFixed(4) },
  perType
};
ctx.writeFile('scratch/v53c/results/length_slope.json', JSON.stringify(out, null, 1));
return { pooled: out.pooled, acrossTypes: out.acrossTypes };
