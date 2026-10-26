// Kingdom classification, specification-compliant and reproducible.
// Rule order (first match wins):
//   1. synthetic construct            -> 'synthetic'
//   2. virus / phage / bacteriophage  -> 'virus'
//   3. exact species match (EUK_SPECIES)          -> 'euk'
//   4. genus in PROK_GENERA                       -> 'prok'
//   5. genus in EUK_GENERA                        -> 'euk'
//   6. regex fallbacks (prokaryote words, then eukaryote words) -> 'prok' / 'euk'
//   7. otherwise                                  -> 'unknown'
// Conventions: archaea are counted as prokaryotes; 'virus' and 'synthetic' are excluded from the
// euk/prok comparison; 'unknown' is excluded and its size reported.
//
// Genera are declared ONCE, in one list, with their kingdom. A genus cannot appear twice.
const PROK_GENERA = new Set(`
escherichia bacillus pseudomonas mycobacterium staphylococcus streptococcus thermus thermotoga
bacteroides vibrio pyrococcus enterococcus helicobacter salmonella klebsiella lactococcus
lactobacillus clostridium corynebacterium rhodobacter synechocystis synechococcus anabaena
nostoc thermococcus archaeoglobus methanocaldococcus methanobacterium halobacterium haloarcula
sulfolobus aeropyrum nanoarchaeum mycoplasma ureaplasma chlamydia rickettsia borrelia treponema
leptospira neisseria haemophilus pasteurella yersinia shigella proteus serratia acinetobacter
burkholderia ralstonia campylobacter deinococcus geobacillus enterobacter citrobacter
stenotrophomonas xanthomonas agrobacterium rhizobium sinorhizobium mesorhizobium bradyrhizobium
caulobacter sphingomonas flavobacterium cytophaga aquifex chlorobium chloroflexus
rhodopseudomonas rhodospirillum thiobacillus nitrosomonas desulfovibrio geobacter shewanella
psychrobacter moraxella halomonas chromohalobacter ferroplasma picrophilus thermoplasma
methanosarcina methanopyrus methanothermobacter methanococcus archaeon thermoproteus
pyrobaculum pyrococcus methanothermus natronomonas haloferax vibrionales streptomyces
`.trim().split(/\s+/));

const EUK_GENERA = new Set(`
homo mus rattus saccharomyces arabidopsis drosophila caenorhabditis bos oryctolagus sus danio
gallus xenopus schizosaccharomyces candida plasmodium trypanosoma leishmania aspergillus
neurospora kluyveromyces pichia komagataella oryza zea nicotiana solanum chlamydomonas
tetrahymena dictyostelium entamoeba giardia trichomonas toxoplasma cryptosporidium spodoptera
trichoplusia bombyx apis strongylocentrotus ciona branchiostoma petromyzon takifugu oryzias
macaca pan canis felis ovis equus cavia mesocricetus cricetulus chlorocebus saimiri mustela
loxodonta ornithorhynchus monodelphis taeniopygia anas meleagris columba physcomitrella
selaginella populus vitis glycine medicago lactuca helianthus sorghum brachypodium physarum
paramecium euglena naegleria cerevisiae pombe
`.trim().split(/\s+/));

const EUK_SPECIES = new Set([
  'homo sapiens','mus musculus','rattus norvegicus','saccharomyces cerevisiae',
  'arabidopsis thaliana','drosophila melanogaster','caenorhabditis elegans','bos taurus',
  'oryctolagus cuniculus','sus scrofa','danio rerio','gallus gallus','xenopus laevis',
  'schizosaccharomyces pombe','candida albicans','plasmodium falciparum','trypanosoma brucei',
  'leishmania major','aspergillus niger','neurospora crassa','kluyveromyces lactis',
  'pichia pastoris','komagataella phaffii','oryza sativa','zea mays','nicotiana tabacum',
  'solanum lycopersicum','chlamydomonas reinhardtii','tetrahymena thermophila',
  'dictyostelium discoideum','entamoeba histolytica','giardia lamblia','trichomonas vaginalis',
  'toxoplasma gondii','cryptosporidium parvum','spodoptera frugiperda','trichoplusia ni',
  'bombyx mori','apis mellifera','strongylocentrotus purpuratus','ciona intestinalis',
  'branchiostoma lanceolatum','petromyzon marinus','takifugu rubripes','oryzias latipes',
  'macaca mulatta','pan troglodytes','canis lupus familiaris','felis catus','ovis aries',
  'equus caballus','cavia porcellus','mesocricetus auratus','cricetulus griseus',
  'chlorocebus aethiops','saimiri sciureus','mustela putorius furo','loxodonta africana',
  'ornithorhynchus anatinus','monodelphis domestica','taeniopygia guttata',
  'anas platyrhynchos','meleagris gallopavo','columba livia',
]);

const PROK_RE = /\b(escherichia|bacillus|pseudomonas|mycobacterium|streptococcus|staphylococcus|clostridium|salmonella|vibrio|helicobacter|lactobacillus|lactococcus|mycoplasma|chlamydia|rickettsia|synechocystis|thermus|sulfolobus|halobacterium|methanococcus|archaeon|archaea|bacterium|bacterial)\b/;
const EUK_RE = /\b(homo|mus|rattus|saccharomyces|arabidopsis|drosophila|caenorhabditis|danio|xenopus|gallus|bos|sus|oryza|zea|human|mouse|rat|yeast|thale cress|fruit fly|zebrafish)\b/;

function classify(org) {
  const o = (org || '').toLowerCase().trim();
  if (!o) return 'unknown';
  if (o.includes('synthetic construct')) return 'synthetic';
  if (o.includes('virus') || o.includes('phage') || o.includes('bacteriophage')) return 'virus';
  if (EUK_SPECIES.has(o)) return 'euk';
  const first = o.split(/[\s,;()\/]+/)[0];
  if (PROK_GENERA.has(first)) return 'prok';
  if (EUK_GENERA.has(first)) return 'euk';
  if (PROK_RE.test(o)) return 'prok';
  if (EUK_RE.test(o)) return 'euk';
  return 'unknown';
}

// a genus may not appear in both lists
const overlap = [...PROK_GENERA].filter(g => EUK_GENERA.has(g));

// classify every row of the organism file, keyed by chain (entry+chain, case-folded)
const dec = new TextDecoder();
const t = dec.decode(ctx.readBytes("uploads/0rgkv2a-selected_6908_organisms.csv"));
const lines = ctx.parseCsv(t);
const chains = {};
const unknowns = [];
let rows = 0;
for (let i = 1; i < lines.length; i++) {
  const f = lines[i];
  if (f.length < 4) continue;
  const key = ((f[0] || '').trim() + (f[1] || '').trim()).toLowerCase();
  const org = (f[3] || '').trim();
  const k = classify(org);
  chains[key] = { k, org };
  if (k === 'unknown' && unknowns.length < 25) unknowns.push(org);
  rows++;
}
const counts = {};
for (const key in chains) counts[chains[key].k] = (counts[chains[key].k] || 0) + 1;

ctx.writeFile("scratch/kg/assign.json", JSON.stringify(chains));
ctx.writeFile("scratch/kg/spec.json", JSON.stringify({
  csvParsing:"quote-aware CSV; no external taxonomy refresh",
  prokaryoteFallbackPattern:PROK_RE.source,eukaryoteFallbackPattern:EUK_RE.source,
  ruleOrder: ["synthetic","virus","exact species","prokaryote genus","eukaryote genus","regex fallbacks","unknown"],
  prokGenera: [...PROK_GENERA].sort(),
  eukGenera: [...EUK_GENERA].sort(),
  eukSpecies: [...EUK_SPECIES].sort(),
  conventions: {
    archaea: "counted as prokaryotes",
    virus: "classified 'virus' and excluded from the euk/prok comparison",
    synthetic: "classified 'synthetic' and excluded",
    unknown: "unmatched, excluded, size reported",
  },
}, null, 1));

return { rows, chains: Object.keys(chains).length, counts,
  genusOverlap: overlap, unknownExamples: unknowns.slice(0, 12) };