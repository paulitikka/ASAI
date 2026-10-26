// pc3_stability.js — loading-vector stability of the components.
// Reproduces the manuscript's stability figures: across 60 half-splits and 60 bootstrap
// resamples of the chain set, each replicate loading vector is correlated with the
// corresponding full-data vector after sign orientation.
// Parameters stated in the manuscript Methods: 60 half-splits, 60 bootstrap resamples.
// Inputs: [0] scratch/pool-full.json
const AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
const Q = 99;
const SPLITS = 60, BOOTS = 60, SEED = 53;
const pool = JSON.parse(ctx.readText(ctx.files[0]));
const ch = pool.chains, order = pool.order;

let s = SEED >>> 0;
const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };

function quant(v) {
  if (!v.length) return new Array(Q).fill(0);
  const a = v.slice().sort((x, y) => x - y), o = [];
  for (let q = 1; q <= Q; q++) {
    const pos = (q / (Q + 1)) * (a.length - 1), lo = Math.floor(pos), hi = Math.min(lo + 1, a.length - 1);
    o.push(a[lo] + (a[hi] - a[lo]) * (pos - lo));
  }
  return o;
}
function jacobi(C, n) {
  let A = C.map(r => r.slice()), V = Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  for (let sw = 0; sw < 1000; sw++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
    if (off < 1e-24) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(A[p][q]) < 1e-16) continue;
      const th = 0.5 * Math.atan2(2 * A[p][q], A[q][q] - A[p][p]), c = Math.cos(th), sn = Math.sin(th);
      for (let i = 0; i < n; i++) { const a = A[i][p], b = A[i][q]; A[i][p] = c*a - sn*b; A[i][q] = sn*a + c*b; }
      for (let i = 0; i < n; i++) { const a = A[p][i], b = A[q][i]; A[p][i] = c*a - sn*b; A[q][i] = sn*a + c*b; }
      for (let i = 0; i < n; i++) { const a = V[i][p], b = V[i][q]; V[i][p] = c*a - sn*b; V[i][q] = sn*a + c*b; }
    }
  }
  const ev = Array.from({length:n},(_,i)=>Math.max(A[i][i],0));
  const ord = Array.from({length:n},(_,i)=>i).sort((x,y)=>ev[y]-ev[x]);
  return { V, ord };
}
function loadings(ids) {
  const mat = AA.map((_, j) => quant(ids.filter(id => ch[id].cnt[j] > 0).map(id => ch[id].sum[j] / ch[id].cnt[j])));
  const cm = new Array(Q).fill(0);
  for (let q = 0; q < Q; q++) { let t = 0; for (let i = 0; i < 20; i++) t += mat[i][q]; cm[q] = t / 20; }
  const C = Array.from({length:20},()=>new Array(20).fill(0));
  for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) {
    let t = 0; for (let q = 0; q < Q; q++) t += (mat[i][q]-cm[q]) * (mat[j][q]-cm[q]);
    C[i][j] = t / Q;
  }
  const d = jacobi(C, 20);
  return [0,1,2].map(k => Array.from({length:20},(_,i)=>d.V[i][d.ord[k]]));
}
function corr(a, b) {
  const n = a.length; let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n;
  let sab = 0, sa = 0, sb = 0;
  for (let i = 0; i < n; i++) { const da = a[i]-ma, db = b[i]-mb; sab += da*db; sa += da*da; sb += db*db; }
  return sab / Math.sqrt(sa * sb);
}
const orient = (v, ref) => (corr(v, ref) < 0 ? v.map(x => -x) : v);

const full = loadings(order.slice());
const half = [[], [], []], boot = [[], [], []];
for (let r = 0; r < SPLITS; r++) {
  const sh = order.slice();
  for (let i = sh.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); const t = sh[i]; sh[i] = sh[j]; sh[j] = t; }
  const L = loadings(sh.slice(0, Math.floor(sh.length / 2)));
  for (let k = 0; k < 3; k++) half[k].push(corr(orient(L[k], full[k]), full[k]));
}
for (let r = 0; r < BOOTS; r++) {
  const ids = [];
  for (let i = 0; i < order.length; i++) ids.push(order[Math.floor(rnd() * order.length)]);
  const L = loadings(ids);
  for (let k = 0; k < 3; k++) boot[k].push(corr(orient(L[k], full[k]), full[k]));
}
const sum = a => ({ mean: +(a.reduce((x,y)=>x+y,0)/a.length).toFixed(3),
                    min: +Math.min(...a).toFixed(3), max: +Math.max(...a).toFixed(3) });
const out = {
  params: { halfSplits: SPLITS, bootstrapResamples: BOOTS, seed: SEED },
  halfSplit: { PC1: sum(half[0]), PC2: sum(half[1]), PC3: sum(half[2]) },
  bootstrap: { PC1: sum(boot[0]), PC2: sum(boot[1]), PC3: sum(boot[2]) }
};
ctx.writeFile('scratch/v53d/results/pc3_stability.json', JSON.stringify(out, null, 1));
return out;
