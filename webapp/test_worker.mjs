// Executes the REAL worker.js in Node via a web-streams shim, so runtime bugs
// in the worker (undefined vars, bad message shape) are caught without a browser.
// Run: `node webapp/test_worker.mjs`
import { ReadableStream, TextDecoderStream, DecompressionStream } from 'node:stream/web';
import assert from 'assert';

const posted = [];
globalThis.self = { onmessage: null, postMessage: (m) => posted.push(m) };
globalThis.DecompressionStream = DecompressionStream;
globalThis.TextDecoderStream = TextDecoderStream;

// tiny in-memory genome: chr1 with an AT-rich centromere in the middle, chr2 arm-only
function mk(n, at, seed) { let s = ''; const b = 'ACGT'; const p = [at/2, (1-at)/2, (1-at)/2, at/2]; let x = seed;
  for (let i = 0; i < n; i++) { x = (x*1103515245 + 12345) & 0x7fffffff; let r = x/0x7fffffff, k = 0; while (r > p[k]) { r -= p[k]; k++; } s += b[k]; } return s; }
const text = `>Chr1\n${mk(4000,0.5,1) + mk(3000,0.85,2) + mk(4000,0.5,3)}\n>Chr2\n${mk(6000,0.5,4)}\n`;
const bytes = new TextEncoder().encode(text);
const file = { name: 'mini.fa', size: bytes.length, stream() { return new ReadableStream({ start(c) { c.enqueue(bytes); c.close(); } }); } };

await import(new URL('./worker.js', import.meta.url));
await globalThis.self.onmessage({ data: {
  file,
  enzymes: [{ name: 'DpnII', motif: 'GATC' }, { name: 'AlwI', motif: 'GGATC' }, { name: 'MluCI', motif: 'AATT' }],
  bedList: [{ name: 'Chr1', intervals: [[4000, 7000]] }],
} });

const err = posted.find(m => m.type === 'error'); if (err) { console.error('WORKER ERROR:', err.message); process.exit(1); }
const done = posted.find(m => m.type === 'done'); assert.ok(done, 'worker posted a done message');
for (const s of done.summary) {
  assert.ok(Number.isFinite(s.cenMean) && Number.isFinite(s.densityDiff) && Number.isFinite(s.armMean), `finite metrics for ${s.name}`);
  assert.ok(s.totalSites > 0, `sites found for ${s.name}`);
}
assert.ok(done.positions.Chr1.DpnII.length > 0, 'per-enzyme positions emitted');
assert.strictEqual(done.chromMeta.length, 2, 'two chromosomes reported');
assert.strictEqual(done.perChromRows.length, 3 * 2, 'perChromRows = enzymes × chroms');
assert.ok(done.summary.find(s => s.name === 'MluCI').densityDiff > 0, 'AATT over-cuts the AT-rich centromere');
console.log('worker.js OK —', done.summary.map(s => `${s.name}:Δ${s.densityDiff.toFixed(2)}`).join('  '));
