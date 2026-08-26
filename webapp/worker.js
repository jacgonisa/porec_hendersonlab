// worker.js — streaming FASTA parse + in-silico restriction-enzyme metrics.
// Emits per-enzyme cut POSITIONS per chromosome (transferable Int32Array) so the
// UI can bin density tracks at any resolution and build per-chromosome monomer
// histograms without re-reading the genome. Also emits genome-wide + per-chrom
// metric rows for the table/scatter/CSV.
//
// Fix vs the Python scripts: palindromic motifs (GATC, CATG, CTAG, AATT ...) are
// counted ONCE (the Python scans motif + revcomp and concatenates, double-counting
// palindromes). Non-palindromic enzymes are unchanged (verified vs ground-truth CSV).
//
// ponytail: one chromosome's sequence is held in memory at a time (fine to ~200 Mb
// chromosomes). Plant genomes are well within this.

const IUPAC = { A:'A', C:'C', G:'G', T:'T', R:'[AG]', Y:'[CT]', S:'[GC]', W:'[AT]',
                K:'[GT]', M:'[AC]', B:'[CGT]', D:'[AGT]', H:'[ACT]', V:'[ACG]', N:'[ACGT]' };
const COMP  = { A:'T', T:'A', C:'G', G:'C', R:'Y', Y:'R', S:'S', W:'W', K:'M',
                M:'K', B:'V', V:'B', D:'H', H:'D', N:'N' };
const revComp = (m) => m.split('').reverse().map(c => COMP[c] || c).join('');
const toRegex = (m) => new RegExp(m.split('').map(c => IUPAC[c] || c).join(''), 'g');
const normChrom = (s) => String(s).toLowerCase().split(/[:\s]/)[0].replace(/^chr(omosome)?[_\-]?/, '');

function findPositions(seq, regex, out) {  // non-overlapping, matches Python re.finditer
  regex.lastIndex = 0; let m;
  while ((m = regex.exec(seq)) !== null) { out.push(m.index); if (m.index === regex.lastIndex) regex.lastIndex++; }
}
function diffStats(sorted) {
  let sum = 0, sumSq = 0, n = 0;
  for (let i = 1; i < sorted.length; i++) { const d = sorted[i] - sorted[i-1]; sum += d; sumSq += d*d; n++; }
  const mean = n > 0 ? sum / n : 0;
  return { sum, sumSq, n, mean, std: n > 0 ? Math.sqrt(Math.max(0, sumSq / n - mean * mean)) : 0 };
}

self.onmessage = async (e) => {
  const { file, enzymes, bedList } = e.data;
  try { await run(file, enzymes, bedList); }
  catch (err) { self.postMessage({ type: 'error', message: err.message + (err.stack ? '\n' + err.stack : '') }); }
};

async function run(file, enzymes, bedList) {
  const enz = enzymes.map(en => {
    const rc = revComp(en.motif);
    return { name: en.name, motif: en.motif, fwd: toRegex(en.motif), rev: toRegex(rc), palindrome: rc === en.motif };
  });
  const bedMap = new Map();
  for (const b of bedList) bedMap.set(normChrom(b.name), { name: b.name, intervals: b.intervals, used: false });
  const getIv = (fastaName) => { const e = bedMap.get(normChrom(fastaName)); if (e) { e.used = true; return e.intervals; } return []; };

  // genome-wide accumulators per enzyme
  const acc = enz.map(() => ({ totalSites: 0, cenSites: 0, armSites: 0, cenSize: 0, armSize: 0,
    cenSum: 0, cenSumSq: 0, cenN: 0, armSum: 0, armSumSq: 0, armN: 0 }));
  const perChromRows = [];
  const positions = {};              // positions[chrom][enzymeName] = Int32Array (cut positions)
  const chromMeta = [];              // [{name, len, intervals}]

  const total = file.size || 0;
  let readBytes = 0, lastPct = -1;
  let stream = file.stream();
  if (/\.gz$/i.test(file.name)) stream = stream.pipeThrough(new DecompressionStream('gzip'));
  const reader = stream.pipeThrough(new TextDecoderStream()).getReader();

  let curName = null, curParts = [], partial = '';
  const flush = () => { if (curName === null) return; processChrom(curName, curParts.join('').toUpperCase()); curParts = []; };

  function processChrom(name, seq) {
    const len = seq.length;
    const intervals = getIv(name);
    const cenSize = intervals.reduce((s, [a, b]) => s + (b - a), 0), armSize = len - cenSize;
    chromMeta.push({ name, len, intervals });
    positions[name] = {};
    const inCen = (p) => { for (const [a, b] of intervals) if (p >= a && p <= b) return true; return false; };

    for (let ei = 0; ei < enz.length; ei++) {
      const en = enz[ei], A = acc[ei];
      self.postMessage({ type: 'progress', phase: 'scan', chrom: name, enzI: ei, enzN: enz.length });
      const buf = [];
      findPositions(seq, en.fwd, buf);
      if (!en.palindrome) findPositions(seq, en.rev, buf);
      const sites = buf.length ? Int32Array.from(new Set(buf)).sort() : new Int32Array(0);
      positions[name][en.name] = sites;

      const cen = [], arm = [];
      for (const p of sites) (inCen(p) ? cen : arm).push(p);
      const cSp = diffStats(cen), aSp = diffStats(arm);
      A.totalSites += sites.length; A.cenSites += cen.length; A.armSites += arm.length;
      A.cenSize += cenSize; A.armSize += armSize;
      A.cenSum += cSp.sum; A.cenSumSq += cSp.sumSq; A.cenN += cSp.n;
      A.armSum += aSp.sum; A.armSumSq += aSp.sumSq; A.armN += aSp.n;

      const cD = cenSize > 0 ? cen.length / cenSize * 1000 : 0, aD = armSize > 0 ? arm.length / armSize * 1000 : 0;
      const eps = 1e-6;
      const cHom = cSp.mean > 0 ? cSp.mean / (cSp.std + eps) : NaN, aHom = aSp.mean > 0 ? aSp.mean / (aSp.std + eps) : NaN;
      perChromRows.push({
        Chromosome: name, Enzyme: en.name, 'Enzyme Recognition Sequence': en.motif,
        'Total Sites': sites.length, 'Centromere Sites': cen.length, 'Arm Sites': arm.length,
        'Centromere Density (cuts/kb)': cD, 'Arm Density (cuts/kb)': aD, 'Density Difference (Cen - Arm)': cD - aD,
        'Centromere/Arm Ratio': aD > 0 ? cD / aD : Infinity,
        'Centromere Mean Spacing (bp)': cSp.mean, 'Centromere Spacing Std Dev': cSp.std,
        'Arm Mean Spacing (bp)': aSp.mean, 'Arm Spacing Std Dev': aSp.std,
        'Centromere Homogeneity': cHom, 'Arm Homogeneity': aHom, 'Homogeneity Difference (Cen - Arm)': cHom - aHom,
      });
    }
  }

  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    readBytes += value.length;
    const pct = total ? Math.min(99, Math.floor(readBytes / total * 100)) : 0;
    if (pct !== lastPct) { lastPct = pct; self.postMessage({ type: 'progress', phase: 'read', pct }); }
    let text = partial + value, nl, start = 0;
    while ((nl = text.indexOf('\n', start)) !== -1) {
      let line = text.slice(start, nl); if (line.endsWith('\r')) line = line.slice(0, -1); start = nl + 1;
      if (line[0] === '>') { flush(); curName = line.slice(1).split(/\s+/)[0]; }
      else if (curName !== null && line) curParts.push(line);
    }
    partial = text.slice(start);
  }
  if (partial) { let line = partial.endsWith('\r') ? partial.slice(0, -1) : partial;
    if (line[0] === '>') { flush(); curName = line.slice(1).split(/\s+/)[0]; } else if (curName !== null && line) curParts.push(line); }
  flush();

  const eps = 1e-6;
  const summary = enz.map((en, i) => {
    const a = acc[i];
    const cD = a.cenSize > 0 ? a.cenSites / a.cenSize * 1000 : 0, aD = a.armSize > 0 ? a.armSites / a.armSize * 1000 : 0;
    const cMean = a.cenN > 0 ? a.cenSum / a.cenN : 0, cStd = a.cenN > 0 ? Math.sqrt(Math.max(0, a.cenSumSq / a.cenN - cMean*cMean)) : 0;
    const aMean = a.armN > 0 ? a.armSum / a.armN : 0, aStd = a.armN > 0 ? Math.sqrt(Math.max(0, a.armSumSq / a.armN - aMean*aMean)) : 0;
    return { name: en.name, motif: en.motif, baseLen: en.motif.length, palindrome: en.palindrome,
      totalSites: a.totalSites, cenSites: a.cenSites, armSites: a.armSites,
      cenDensity: cD, armDensity: aD, densityDiff: cD - aD, ratio: aD > 0 ? cD / aD : Infinity,
      cenMean, cenStd, armMean, armStd,
      cenHom: cMean > 0 ? cMean / (cStd + eps) : NaN, armHom: aMean > 0 ? aMean / (aStd + eps) : NaN,
      homDiff: (cMean > 0 ? cMean / (cStd + eps) : NaN) - (aMean > 0 ? aMean / (aStd + eps) : NaN) };
  });

  const unusedBed = [...bedMap.values()].filter(e => !e.used).map(e => e.name);
  const transfer = [];
  for (const c of Object.values(positions)) for (const arr of Object.values(c)) transfer.push(arr.buffer);
  self.postMessage({ type: 'done', summary, perChromRows, positions, chromMeta, unusedBed }, transfer);
}
