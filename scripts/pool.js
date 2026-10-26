// Pool all six residue parts into one accumulator, keyed exactly as collect_a.js does:
//   chainKey = (pdb + chain).toLowerCase(), per-chain length, per-type count, per-type RSA sum.
// Plus whole-dataset tallies needed to check Supplementary Table S4 and reconcile.json.
const AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];
const acc = {};
let rows = 0, badRows = 0, unknownAA = 0, nonFinite = 0;
let zeroCount = 0;
// S4 buckets: exactly 0, then 0-10, 10-20, ... 80-90, 90+
const buckets = [0,0,0,0,0,0,0,0,0,0,0];
let sampleMaxRsa = -Infinity, sampleMinRsa = Infinity;
const chainSet = new Set();

for (const p of ctx.files) {
  const t = ctx.readText(p);
  const ls = t.split('\n');
  for (let i = 1; i < ls.length; i++) {
    const line = ls[i];
    if (!line || !line.trim()) continue;
    const f = line.split(',');
    if (f.length < 7) { badRows++; continue; }
    const j = AA.indexOf(f[5]);
    if (j < 0) { unknownAA++; continue; }
    const rsa = parseFloat(f[6]);
    if (!isFinite(rsa)) { nonFinite++; continue; }
    const k = (f[0] + f[1]).toLowerCase();
    const c = acc[k] || (acc[k] = { L:0, cnt:new Array(20).fill(0), sum:new Array(20).fill(0) });
    c.L++; c.cnt[j]++; c.sum[j] += rsa;
    chainSet.add(k);
    rows++;
    if (rsa === 0) zeroCount++;
    buckets[ctx.numericPolicy.rsaInterval(rsa)]++;
    if (rsa > sampleMaxRsa) sampleMaxRsa = rsa;
    if (rsa < sampleMinRsa) sampleMinRsa = rsa;
  }
}

// reconcile.js logic: all-twenty-types readings and the <40 chain count
let all20 = 0, all20res = 0, ge5 = 0, ge5res = 0, ge20 = 0, ge20res = 0;
let lt40 = 0, ge40 = 0;
for (const k of Object.keys(acc)) {
  const c = acc[k];
  let n20 = 0, n5 = 0, n20r = 0;
  for (let j = 0; j < 20; j++) { if (c.cnt[j] >= 1) n20++; if (c.cnt[j] >= 5) n5++; if (c.cnt[j] >= 20) n20r++; }
  if (n20 === 20) { all20++; all20res += c.L; }
  if (n5 === 20) { ge5++; ge5res += c.L; }
  if (n20r === 20) { ge20++; ge20res += c.L; }
  if (c.L < 40) lt40++; else ge40++;
}

// per-type totals and mean RSA (Table 1 Mean RSA column, and the ASAI ordering)
const perType = AA.map((a, j) => {
  let n = 0, s = 0;
  for (const k of Object.keys(acc)) { n += acc[k].cnt[j]; s += acc[k].sum[j]; }
  return { a, n, meanRSA: n ? s / n : null };
});

// persist the accumulator so later steps can reuse it without re-reading 51 MB
ctx.writeFile('scratch/pooled-acc.json', JSON.stringify(acc));
ctx.writeFile('scratch/pooled-summary.json', JSON.stringify({
  rows, badRows, unknownAA, nonFinite, zeroCount, buckets,
  chains: Object.keys(acc).length, all20, all20res, ge5, ge5res, ge20, ge20res,
  byLen: { lt40, ge40 }, perType, sampleMinRsa, sampleMaxRsa
}, null, 1));

return {
  rows, badRows, unknownAA, nonFinite, zeroCount,
  chains: Object.keys(acc).length,
  all20, all20res, ge5, ge5res, ge20, ge20res,
  byLen: { lt40, ge40 },
  buckets,
  minRsa: sampleMinRsa, maxRsa: sampleMaxRsa
};