// aaindex_ranks.js — rank enrichment of AAindex1 hydrophobicity scales against PC1/PC2/PC3.
//
// Method (Methods, "Rank-enrichment test against AAindex1"): all 553 complete entries of
// AAindex1 are ranked by the absolute Pearson correlation of their twenty value vectors with
// each component. A scale is labelled hydrophobicity-related when its DESCRIPTION contains any
// of these six terms:
//     hydrophobicity, hydropathy, hydrophilicity, partition energy,
//     transfer free energy, solvent accessibility
// The labelling is fixed before the correlations are computed.
//
// AAindex1 record layout (one record per "//" block):
//     H <accession>
//     D <description>            (may run over several D lines)
//     ...
//     I    A/L     R/K     N/M ...     <- label row
//          v1      v2      v3  ...     <- ten values
//          v11     v12     v13 ...     <- ten values
// The twenty numbers follow the I label row, over one or two continuation lines.
// The I value order is A R N D C Q E G H I L K M F P S T W Y V.
//
// Input: inputs/aaindex1.txt  (passed as ctx.files[0])

const AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'];
const ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'];
const TERMS = ['hydrophobicity','hydropathy','hydrophilicity','partition energy','transfer free energy','solvent accessibility'];

const text = ctx.readText(ctx.files[0]);
const records = [];

for (const block of text.split('//')) {
  const lines = block.split('\n');
  let acc = null;
  const descParts = [];
  let iIdx = -1;
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    if (iIdx < 0 && /^\s*H\s/.test(l)) acc = l.replace(/^\s*H\s+/, '').trim().split(/\s+/)[0];
    else if (iIdx < 0 && /^\s*D\s/.test(l)) descParts.push(l.replace(/^\s*D\s+/, '').trim());
    else if (iIdx < 0 && /^\s*I\s/.test(l)) iIdx = i;
  }
  if (!acc || iIdx < 0) continue;

  // collect numbers from the I line and the lines after it, stopping at the next field
  let nums = [];
  for (let j = iIdx; j < lines.length && nums.length < 20; j++) {
    const line = lines[j];
    if (j > iIdx && /^\s*[A-Z]\s/.test(line) && !/^\s*\d/.test(line)) break; // a new field letter
    const found = (line.match(/-?\d+(?:\.\d+)?/g) || []).map(Number);
    nums = nums.concat(found);
  }
  if (nums.length < 20) continue;
  const vals = nums.slice(0, 20);
  const desc = descParts.join(' ');
  records.push({ acc, desc, vals, hydro: TERMS.some(t => desc.toLowerCase().includes(t)) });
}

function toOurs(v) { const m = {}; ORDER.forEach((a, i) => m[a] = v[i]); return AA.map(a => m[a]); }
function pearson(a, b) {
  const n = a.length;
  const ma = a.reduce((x, y) => x + y, 0) / n;
  const mb = b.reduce((x, y) => x + y, 0) / n;
  let sab = 0, sa = 0, sb = 0;
  for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; sa += da * da; sb += db * db; }
  return sab / Math.sqrt(sa * sb);
}

const COORD = {
  PC1: {A:-0.0825,C:-0.3628,D:0.2533,E:0.3258,F:-0.2669,G:0.0593,H:0.0544,I:-0.2743,K:0.3873,L:-0.2426,M:-0.1562,N:0.1867,P:0.1955,Q:0.2194,R:0.2200,S:0.0645,T:0.0346,V:-0.2352,W:-0.2254,Y:-0.1550},
  PC2: {A:0.0746,C:0.0758,D:0.1655,E:0.2580,F:0.1068,G:-0.0400,H:-0.4627,I:0.2028,K:0.1976,L:0.2999,M:-0.5794,N:-0.1572,P:-0.1423,Q:-0.0376,R:-0.0172,S:-0.0053,T:0.0457,V:0.2524,W:-0.2216,Y:-0.0160},
  PC3: {A:-0.3895,C:0.2121,D:0.2498,E:0.2641,F:0.0368,G:0.2716,H:-0.2559,I:-0.0531,K:0.2718,L:-0.0130,M:0.3913,N:-0.2966,P:0.0527,Q:-0.2612,R:-0.1120,S:-0.2155,T:-0.1884,V:-0.0917,W:0.1873,Y:-0.0605}
};

const perPC = {};
for (const key of ['PC1','PC2','PC3']) {
  const comp = AA.map(a => COORD[key][a]);
  const ranked = records
    .map(r => ({ acc: r.acc, desc: r.desc.slice(0, 70), hydro: r.hydro, r: +Math.abs(pearson(toOurs(r.vals), comp)).toFixed(3) }))
    .sort((x, y) => y.r - x.r);
  ranked.forEach((s, i) => s.rank = i + 1);

  const hydroRanks = ranked.filter(s => s.hydro).map(s => s.rank).sort((a, b) => a - b);
  const nonHydroRanks = ranked.filter(s => !s.hydro).map(s => s.rank);
  const median = hydroRanks.length ? hydroRanks[Math.floor((hydroRanks.length - 1) / 2)] : null;
  const averageRank = +(ranked.reduce((s, x) => s + x.rank, 0) / ranked.length).toFixed(1);

  let U = 0;
  for (const r of hydroRanks) for (const q of nonHydroRanks) if (r > q) U++;
  const n1 = hydroRanks.length, n2 = nonHydroRanks.length;
  const mu = n1 * n2 / 2;
  const sd = Math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12);
  const z = sd ? +((U - mu) / sd).toFixed(2) : null;

  perPC[key] = { medianHydroRank: median, averageRank, nHydro: n1, mannWhitneyU: U, z,
    top10: ranked.slice(0, 10).map(s => ({ acc: s.acc, desc: s.desc, r: s.r, hydro: s.hydro })) };
}

const out = { entries: records.length, hydroCount: records.filter(r => r.hydro).length, terms: TERMS, perPC };
ctx.writeFile('scratch/bundle/results/aaindex_ranks.json', JSON.stringify(out, null, 1));
return {
  entries: out.entries, hydroCount: out.hydroCount,
  PC1: { median: perPC.PC1.medianHydroRank, avg: perPC.PC1.averageRank, z: perPC.PC1.z },
  PC2: { median: perPC.PC2.medianHydroRank, avg: perPC.PC2.averageRank, z: perPC.PC2.z, nHydro: perPC.PC2.nHydro,
         top3: perPC.PC2.top10.slice(0, 3).map(s => s.acc + '=' + s.r) },
  PC3: { median: perPC.PC3.medianHydroRank, avg: perPC.PC3.averageRank, z: perPC.PC3.z }
};
