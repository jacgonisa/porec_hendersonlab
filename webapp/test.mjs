// Regression test for the core logic. Run: `node webapp/test.mjs`
// The heavy algorithm is validated against the real genome + ground-truth CSV
// (see README). This guards the browser-side helpers. Self-contained — no big files.
// Helpers below MIRROR worker.js / app.js — keep in sync.
import fs from 'fs';
import assert from 'assert';

const COMP = { A:'T',T:'A',C:'G',G:'C',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N' };
const IUPAC = { A:'A',C:'C',G:'G',T:'T',R:'[AG]',Y:'[CT]',S:'[GC]',W:'[AT]',K:'[GT]',M:'[AC]',B:'[CGT]',D:'[AGT]',H:'[ACT]',V:'[ACG]',N:'[ACGT]' };
const revComp = m => m.split('').reverse().map(c => COMP[c] || c).join('');
const toRegex = m => new RegExp(m.split('').map(c => IUPAC[c] || c).join(''), 'g');
const normChrom = s => String(s).toLowerCase().split(/[:\s]/)[0].replace(/^chr(omosome)?[_\-]?/, '');
function findSites(seq, motif) {
  const pal = revComp(motif) === motif, out = [];
  const scan = r => { r.lastIndex = 0; let m; while ((m = r.exec(seq))) { out.push(m.index); if (m.index === r.lastIndex) r.lastIndex++; } };
  scan(toRegex(motif)); if (!pal) scan(toRegex(revComp(motif)));
  return Int32Array.from(new Set(out)).sort();
}
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

// 1. BED parsing: standard 3-col, lab id/name/start/end, and the bundled example format.
assert.deepStrictEqual(parseBed('Chr1\t100\t200')[0].intervals[0], [100, 200]);
assert.deepStrictEqual(parseBed('1\tChr1:1-32640075\t14841147\t17216861')[0].intervals[0], [14841147, 17216861]);
const exBed = parseBed(fs.readFileSync(new URL('./example/arabidopsis.bed', import.meta.url).pathname, 'utf8'));
assert.strictEqual(exBed.length, 5, 'example bed has 5 centromeres');

// 2. Chromosome-name normalization: "Chr1:1-..", "chr1", "1" all match; Chr10 ≠ Chr1.
assert.strictEqual(normChrom('Chr1:1-32640075'), normChrom('1'));
assert.strictEqual(normChrom('chr1'), normChrom('Chr1'));
assert.notStrictEqual(normChrom('Chr10'), normChrom('Chr1'));
assert.strictEqual(normChrom(exBed[0].name), '1', 'example bed chr1 normalizes to "1"');

// 3. Palindrome dedupe: AATT (revcomp==self) counted once; GGATC (both strands) twice.
assert.strictEqual(findSites('AATTAATT', 'AATT').length, 2, 'palindrome not double-counted');
assert.strictEqual(findSites('GGATCxGATCC', 'GGATC').length, 2, 'both strands of non-palindrome found');

// 4. NEB catalog integrity.
const neb = JSON.parse(fs.readFileSync(new URL('./neb_enzymes.json', import.meta.url).pathname, 'utf8'));
assert.ok(neb.length > 150, `expected 200ish NEB sites, got ${neb.length}`);
assert.ok(neb.every(e => e.len >= 4 && /^[ACGTRYSWKMBDHVN]+$/.test(e.site)), 'all sites valid IUPAC, len>=4');
for (const b of ['AlwI', 'DpnII', 'NlaIII', 'BfaI', 'MluCI']) assert.ok(neb.some(e => e.enzymes.includes(b)), `${b} present`);

// 5. AT-cutter over-digests an AT-rich centromere (density diff logic on a built genome).
function seqIID(n, gc, seed) { let s = ''; const b = 'ACGT'; const p = [(1-gc)/2, gc/2, gc/2, (1-gc)/2]; let x = seed;
  for (let i = 0; i < n; i++) { x = (x * 1103515245 + 12345) & 0x7fffffff; let r = (x / 0x7fffffff), k = 0; while (r > p[k]) { r -= p[k]; k++; } s += b[k]; } return s; }
const genome = seqIID(30000, 0.38, 1) + seqIID(20000, 0.20, 2) + seqIID(30000, 0.38, 3); // AT-rich middle
const iv = [[30000, 50000]];
function density(motif) {
  const pos = findSites(genome, motif), cen = [], arm = [];
  for (const p of pos) (p >= 30000 && p < 50000 ? cen : arm).push(p);
  return { cD: cen.length / 20000 * 1000, aD: arm.length / 60000 * 1000 };
}
const aatt = density('AATT');
assert.ok(aatt.cD > 1.5 * aatt.aD, `AATT over-digests AT-rich centromere (CEN ${aatt.cD.toFixed(2)} vs arm ${aatt.aD.toFixed(2)})`);

console.log('All tests passed.');
