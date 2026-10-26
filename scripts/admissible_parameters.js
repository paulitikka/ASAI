// admissible_parameters.js — context-setting parameters and the joint models.
// Reproduces v5b: Table 2 (parameters and terms), Table 3 (joint models), 13.6% / 86.4%.
// ADMISSIBILITY RULE: a parameter is admitted only if computable for a chain without using that
// chain's per-residue RSA values. Burial level and RSA dispersion are excluded by the rule.
// Inputs: [0] scratch/pool-full.json, [1] scratch/bundle/inputs/selected_6908_organisms.csv
const AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];
const CHG = { D:-1, E:-1, K:1, R:1, H:0.1 };
const KD = { A:1.8,R:-4.5,N:-3.5,D:-3.5,C:2.5,Q:-3.5,E:-3.5,G:-0.4,H:-3.2,I:4.5,L:3.8,
             K:-3.9,M:1.9,F:2.8,P:-1.6,S:-0.8,T:-0.7,W:-0.9,Y:-1.3,V:4.2 };
const EP = 0.5;                 // 0.5-residue pseudocount, as in the composition models
const RIDGE = 1e-9;             // ridge relative to the diagonal
const NSPLIT = 20;              // held-out splits
const SEED = 20260930;

const pool = JSON.parse(ctx.readText(ctx.files[0]));
const ch = pool.chains, order = pool.order;
const seq = {};
for(const f of ctx.parseCsv(ctx.readText(ctx.files[1])).slice(1)) {
  if(f.length>=3)seq[(f[0].trim()+f[1].trim()).toLowerCase()]=f[2].trim();
}
const ids = order.filter(id => ch[id].L >= 40 && seq[id]);

// ---- per-chain context-setting parameters ----
const P = {};
for (const id of ids) {
  const c = ch[id], s = seq[id];
  let chg = 0, pos = 0, neg = 0;
  for (let j = 0; j < 20; j++) { const q = CHG[AA[j]] || 0; chg += q*c.cnt[j];
    if (q > 0) pos += c.cnt[j]; else if (q < 0) neg += c.cnt[j]; }
  const h = new Float64Array(s.length);
  for (let i = 0; i < s.length; i++) { const v = KD[s[i]]; h[i] = v === undefined ? 0 : v; }
  let hMean = 0; for (let i = 0; i < h.length; i++) hMean += h[i]; hMean /= h.length;
  let hMax = -Infinity, run = 0; const W = 11;
  for (let i = 0; i < h.length; i++) { run += h[i]; if (i >= W) run -= h[i-W];
    if (i >= W-1) { const m = run/W; if (m > hMax) hMax = m; } }
  let inLow = 0, tot = 0;
  for (let i = 0; i + 20 <= s.length; i += 20) {
    let mask = 0; for (let k = 0; k < 20; k++) mask |= 1 << (s.charCodeAt(i+k) % 32);
    let bits = 0; for (let b = 0; b < 32; b++) if (mask & (1 << b)) bits++;
    tot += 20; if (bits <= 12) inLow += 20;
  }
  const Lp = c.L + 20*EP;
  P[id] = {
    length: Math.log10(c.L),                       // 1 term
    comp: c.cnt.map(v => (v + EP) / Lp),           // 20 terms
    charge: chg / c.L,                             // 1 term
    asymmetry: (pos - neg) / c.L,                  // 1 term
    hydroMean: hMean,                              // 3 terms (with the two below)
    hydroWindow: hMax === -Infinity ? 0 : hMax,
    lowComplexity: tot ? inLow / tot : 0
  };
}
const SETS = {
  length:      id => [1, P[id].length],
  composition: id => [1, ...P[id].comp],
  charge:      id => [1, P[id].charge, P[id].asymmetry],
  hydrophobicity: id => [1, P[id].hydroMean, P[id].hydroWindow, P[id].lowComplexity],
  all:         id => [1, ...P[id].comp, P[id].charge, P[id].asymmetry,
                      P[id].hydroMean, P[id].hydroWindow, P[id].lowComplexity]
};
// length is part of every model except the length-only one, matching the paper's usage
const SETS_WITH_LENGTH = {
  length:      id => [1, P[id].length],
  composition: id => [1, ...P[id].comp, P[id].length],
  charge:      id => [1, P[id].charge, P[id].asymmetry, P[id].length],
  hydrophobicity: id => [1, P[id].hydroMean, P[id].hydroWindow, P[id].lowComplexity, P[id].length],
  all:         id => [1, ...P[id].comp, P[id].charge, P[id].asymmetry,
                      P[id].hydroMean, P[id].hydroWindow, P[id].lowComplexity, P[id].length]
};

function olsFit(X, y) {
  const p = X[0].length, n = X.length;
  const A = Array.from({length:p},()=>new Array(p).fill(0)), b = new Array(p).fill(0);
  for (let i=0;i<n;i++){ const xi=X[i];
    for (let a=0;a<p;a++){ b[a]+=xi[a]*y[i]; for (let c=a;c<p;c++) A[a][c]+=xi[a]*xi[c]; } }
  for (let a=0;a<p;a++) for (let c=0;c<a;c++) A[a][c]=A[c][a];
  let tr=0; for (let a=0;a<p;a++) tr+=A[a][a];
  for (let a=0;a<p;a++) A[a][a] += RIDGE*(tr/p);
  const M = A.map((r,i)=>[...r,b[i]]);
  for (let col=0;col<p;col++){
    let piv=col; for (let r=col+1;r<p;r++) if (Math.abs(M[r][col])>Math.abs(M[piv][col])) piv=r;
    if (piv!==col){ const t=M[col]; M[col]=M[piv]; M[piv]=t; }
    const d=M[col][col]; if (Math.abs(d)<1e-14) return null;
    for (let c=col;c<=p;c++) M[col][c]/=d;
    for (let r=0;r<p;r++){ if (r===col) continue; const f=M[r][col]; if (!f) continue;
      for (let c=col;c<=p;c++) M[r][c]-=f*M[col][c]; }
  }
  return M.map(r=>r[p]);
}
let lastEvaluation=null;
function r2On(beta, X, y) {
  const n = y.length; let yb = 0; for (const v of y) yb += v; yb /= n;
  let ssr = 0, sst = 0;
  for (let i=0;i<n;i++){ let yh=0; for (let a=0;a<beta.length;a++) yh += beta[a]*X[i][a];
    ssr += (y[i]-yh)**2; sst += (y[i]-yb)**2; }
  lastEvaluation={n_test:n,SSE:ssr,SST:sst};
  return sst === 0 ? null : 1 - ssr/sst;
}
const byType = AA.map((_, j) => ids.filter(id => ch[id].cnt[j] > 0));
const resp = byType.map((l, j) => { const m = {}; for (const id of l) m[id] = ch[id].sum[j]/ch[id].cnt[j]; return m; });
const X = {}; for (const id of ids) { X[id] = {}; for (const nm of Object.keys(SETS_WITH_LENGTH)) X[id][nm] = SETS_WITH_LENGTH[nm](id); }

let seed = SEED >>> 0;
const rnd = () => { seed = (seed*1664525 + 1013904223) >>> 0; return seed/4294967296; };
const cell = [];
const modelDiagnostics=[];
for (let sp = 0; sp < NSPLIT; sp++) {
  for (let j = 0; j < 20; j++) {
    const l = byType[j];
    if (l.length < 200) continue;
    const assigned = ctx.splitDesign.partition(AA[j],sp,l);
    const tr=assigned.train,te=assigned.test;
    const yTr = tr.map(id=>resp[j][id]), yTe = te.map(id=>resp[j][id]);
    const rec = { split: sp, type: AA[j] };
    for (const nm of Object.keys(SETS_WITH_LENGTH)) {
      const beta = olsFit(tr.map(id=>X[id][nm]), yTr);
      if(!beta)throw new Error("Model fit failed: "+AA[j]+" split "+(sp+1)+" "+nm);
      rec[nm] = r2On(beta, te.map(id=>X[id][nm]), yTe);
      modelDiagnostics.push({split:sp+1,aa:AA[j],model:nm,n_train:tr.length,...lastEvaluation,R2:rec[nm]});
    }
    cell.push(rec);
  }
}
const mean = a => a.reduce((x,y)=>x+y,0)/a.length;
const models = {};
for (const nm of Object.keys(SETS_WITH_LENGTH)) {
  const v = cell.map(r=>r[nm]).filter(x=>x!==null);
  const m = mean(v);
  models[nm] = { meanHeldOutR2: +m.toFixed(4), unexplainedPct: +(100*(1-m)).toFixed(1), n: v.length };
}
// paired test: all-composition, across splits, per type
function paired(u, v) {
  const ts = [];
  for (let sp = 0; sp < NSPLIT; sp++) {
    const rows = cell.filter(r=>r.split===sp && r[u]!==null && r[v]!==null);
    if (rows.length < 5) continue;
    const d = rows.map(r=>r[u]-r[v]);
    const m = mean(d), s = Math.sqrt(d.reduce((a,b)=>a+(b-m)**2,0)/Math.max(1,d.length-1));
    ts.push(s === 0 ? 0 : m/(s/Math.sqrt(d.length)));
  }
  const mt = mean(ts);
  return { meanDiff: +mean(cell.filter(r=>r[u]!==null&&r[v]!==null).map(r=>r[u]-r[v])).toFixed(4),
           meanPerSplitT: +mt.toFixed(2), splitsFavouring: ts.filter(t=>t>0).length, nSplits: ts.length };
}
ctx.writeFile("chain_predictors.json",JSON.stringify(P));
ctx.writeFile("model_scores_by_split.json",JSON.stringify(cell,null,1));
ctx.writeFile("model_diagnostics.json",JSON.stringify(modelDiagnostics,null,1));
const out = {
  chains: ids.length,
  parameters: [
    { n: 1, name: 'Chain length', terms: '1 (log10)', nTerms: 1 },
    { n: 2, name: 'Amino-acid composition', terms: '20 fractions', nTerms: 20 },
    { n: 3, name: 'Net charge', terms: '1', nTerms: 1 },
    { n: 4, name: 'Charge asymmetry', terms: '1', nTerms: 1 },
    { n: 5, name: 'Sequence-encoded hydrophobicity', terms: '3 (mean, 11-residue maximum, low-complexity fraction)', nTerms: 3 }
  ],
  models: models,
  paired: { 'all minus composition': paired('all','composition'),
            'composition minus length': paired('composition','length') },
  headline: { explainedPct: +(100*models.all.meanHeldOutR2).toFixed(1),
              intrinsicUpperPct: models.all.unexplainedPct }
};
ctx.writeFile('scratch/v53c/results/admissible_parameters.json', JSON.stringify(out, null, 1));
return out;
