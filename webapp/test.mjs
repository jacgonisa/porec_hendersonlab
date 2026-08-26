// Regression test for the core logic. Run: `node webapp/test.mjs`
// The heavy algorithm is validated against the real genome + ground-truth CSV
// (see webapp/README.md); this guards the browser-side helpers using the
// committed example genome. Helpers below MIRROR worker.js / app.js — keep in sync.
import fs from 'fs';
import assert from 'assert';

const DIR = new URL('./example/', import.meta.url).pathname;

// --- mirrors of worker.js -------------------------------------------------
const IUPAC = { A:'A',C:'C',G:'G',T:'T',R:'[AG]',Y:'[CT]',S:'[GC]',W:'[AT]',K:'[GT]',M:'[AC]',B:'[CGT]',D:'[AGT]',H:'[ACT]',V:'[ACG]',N:'[ACGT]' };
const COMP = { A:'T',T:'A',C:'G',G:'C',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N' };
const revComp = m => m.split('').reverse().map(c => COMP[c] || c).join('');
const toRegex = m => new RegExp(m.split('').map(c => IUPAC[c] || c).join(''), 'g');
const normChrom = s => String(s).toLowerCase().split(/[:\s]/)[0].replace(/^chr(omosome)?[_\-]?/, '');
function findSites(seq, motif) {
  const pal = revComp(motif) === motif, out = [];
  const scan = (r) => { r.lastIndex = 0; let m; while ((m = r.exec(seq))) { out.push(m.index); if (m.index === r.lastIndex) r.lastIndex++; } };
  scan(toRegex(motif)); if (!pal) scan(toRegex(revComp(motif)));
  return [...new Set(out)].sort((a, b) => a - b);
}
// --- mirror of app.js parseBed & score ------------------------------------
function parseBed(text) {
  const out = {};
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line[0] === '#' || /^(track|browser)\b/.test(line)) continue;
    const f = line.split(/\s+/); let chrom = null, start, end;
    for (let i = 1; i < f.length - 1; i++) if (/^\d+$/.test(f[i]) && /^\d+$/.test(f[i + 1])) { chrom = f[i - 1]; start = +f[i]; end = +f[i + 1]; break; }
    if (chrom === null || end <= start) continue;
    (out[chrom] = out[chrom] || []).push([start, end]);
  }
  return Object.entries(out).map(([name, intervals]) => ({ name, intervals }));
}
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
function score(s) {
  const sDD = clamp(1 - Math.abs(s.densityDiff) / 2, 0, 1);
  const ls = v => v <= 0 ? 0 : clamp(v < 800 ? (v - 300) / 500 : v <= 1400 ? 1 : (3000 - v) / 1600, 0, 1);
  const sLen = (ls(s.cenMean) + ls(s.armMean)) / 2;
  const sHom = clamp(Math.min(s.cenHom || 0, s.armHom || 0) / 0.8, 0, 1);
  const enough = (s.cenSites >= 200 && s.armSites >= 200) ? 1 : 0.3;
  return (0.5 * sDD + 0.3 * sLen + 0.2 * sHom) * enough;
}

// --- tests ----------------------------------------------------------------
// 1. BED parsing: standard 3-col AND the lab's id/name/start/end 4-col variant.
assert.deepStrictEqual(parseBed('Chr1\t100\t200').map(b => [b.name, b.intervals[0]]), [['Chr1', [100, 200]]]);
assert.deepStrictEqual(parseBed('1\tChr1:1-32640075\t14841147\t17216861')[0].intervals[0], [14841147, 17216861], 'lab 4-col BED');

// 2. Chromosome-name normalization: "Chr1:1-..", "chr1", "1" all match.
assert.strictEqual(normChrom('Chr1:1-32640075'), normChrom('1'));
assert.strictEqual(normChrom('chr1'), normChrom('Chr1'));
assert.notStrictEqual(normChrom('Chr10'), normChrom('Chr1'));

// 3. Palindrome dedupe: AATT counted once (both strands identical), GGATC twice.
assert.strictEqual(findSites('AATTAATT', 'AATT').length, 2, 'palindrome not double-counted');
assert.strictEqual(findSites('GGATCxGATCC', 'GGATC').length, 2, 'both strands of non-palindrome found');

// 4. End-to-end on committed demo genome: AT-cutter (MluCI, AATT) over-digests
//    the AT-rich centromere; balanced enzyme beats an over-digesting one on score.
const chroms = {}; let cur = null, parts = [];
for (const l of fs.readFileSync(DIR + 'demo.fa', 'utf8').split('\n')) {
  if (l[0] === '>') { if (cur) chroms[cur] = parts.join('').toUpperCase(); cur = l.slice(1).trim(); parts = []; }
  else if (cur) parts.push(l.trim());
}
if (cur) chroms[cur] = parts.join('').toUpperCase();
const bed = new Map(parseBed(fs.readFileSync(DIR + 'demo.bed', 'utf8')).map(b => [normChrom(b.name), b.intervals]));

function metrics(motif) {
  let cN = 0, aN = 0, cSize = 0, aSize = 0, cSum = 0, cSq = 0, cD = 0, aSum = 0, aSq = 0, aDn = 0;
  for (const [name, seq] of Object.entries(chroms)) {
    const iv = bed.get(normChrom(name)) || [];
    const cz = iv.reduce((s, [a, b]) => s + (b - a), 0); cSize += cz; aSize += seq.length - cz;
    const cen = [], arm = [];
    for (const p of findSites(seq, motif)) (iv.some(([a, b]) => p >= a && p <= b) ? cen : arm).push(p);
    cN += cen.length; aN += arm.length;
    for (let i = 1; i < cen.length; i++) { const d = cen[i] - cen[i-1]; cSum += d; cSq += d*d; cD++; }
    for (let i = 1; i < arm.length; i++) { const d = arm[i] - arm[i-1]; aSum += d; aSq += d*d; aDn++; }
  }
  const cMean = cD ? cSum/cD : 0, cStd = cD ? Math.sqrt(Math.max(0, cSq/cD - cMean*cMean)) : 0;
  const aMean = aDn ? aSum/aDn : 0, aStd = aDn ? Math.sqrt(Math.max(0, aSq/aDn - aMean*aMean)) : 0;
  return { cenSites: cN, armSites: aN, cenDensity: cN/cSize*1000, armDensity: aN/aSize*1000,
           densityDiff: cN/cSize*1000 - aN/aSize*1000, cenMean: cMean, armMean: aMean,
           cenHom: cMean/(cStd+1e-6), armHom: aMean/(aStd+1e-6) };
}
const mlucI = metrics('AATT'), dpnII = metrics('GATC');
assert.ok(mlucI.cenDensity > 2 * mlucI.armDensity, 'MluCI over-digests AT-rich centromere');
assert.ok(mlucI.densityDiff > dpnII.densityDiff, 'AATT more centromere-biased than GATC');
assert.ok(score(dpnII) > score(mlucI), 'less-biased enzyme scores higher than AT-overdigester');

console.log('All tests passed.');
