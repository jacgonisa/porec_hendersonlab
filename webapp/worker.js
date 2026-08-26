// worker.js — streaming FASTA parse + in-silico restriction-enzyme metrics.
// Runs off the main thread so the UI stays responsive on ~100 MB+ genomes.
//
// Ports insilico_choice/table_with_mean_and_sd_ratio.py, with one fix:
// palindromic motifs (GATC, CATG, CTAG, AATT ...) are counted ONCE, not twice.
// The Python scans motif + reverse-complement and concatenates, so a palindrome
// (whose revcomp == itself) is double-counted. We dedupe positions in a Set.
//
// ponytail: we hold ONE chromosome's sequence in memory at a time (not the whole
// genome). Fine up to ~200 Mb chromosomes; beyond that switch to a rolling-window
// scan. Plant genomes (largest chr ~35 Mb) are well within this.

const IUPAC = { A:'A', C:'C', G:'G', T:'T', R:'[AG]', Y:'[CT]', S:'[GC]', W:'[AT]',
                K:'[GT]', M:'[AC]', B:'[CGT]', D:'[AGT]', H:'[ACT]', V:'[ACG]', N:'[ACGT]' };
const COMP  = { A:'T', T:'A', C:'G', G:'C', R:'Y', Y:'R', S:'S', W:'W', K:'M',
                M:'K', B:'V', V:'B', D:'H', H:'D', N:'N' };

function revComp(motif) {
  return motif.split('').reverse().map(c => COMP[c] || c).join('');
}
function toRegex(motif) {
  return new RegExp(motif.split('').map(c => IUPAC[c] || c).join(''), 'g');
}

// Non-overlapping match positions (matches Python re.finditer).
function findPositions(seq, regex, out) {
  regex.lastIndex = 0;
  let m;
  while ((m = regex.exec(seq)) !== null) {
    out.push(m.index);
    if (m.index === regex.lastIndex) regex.lastIndex++; // guard zero-width
  }
}

// Monomer-length histogram: 0..MONO_MAX bp in MONO_STEP bins + one overflow bin.
const MONO_MAX = 5000, MONO_STEP = 50, MONO_NBINS = MONO_MAX / MONO_STEP + 1;
function histIndex(v) { return v >= MONO_MAX ? MONO_NBINS - 1 : Math.floor(v / MONO_STEP); }

// Normalize a chromosome name so "Chr1:1-32640075", "chr1" and "1" all match.
function normChrom(s) {
  s = String(s).toLowerCase().split(/[:\s]/)[0];
  return s.replace(/^chr(omosome)?[_\-]?/, '');
}

self.onmessage = async (e) => {
  const { file, enzymes, bedList, binSize } = e.data;
  try { await run(file, enzymes, bedList, binSize); }
  catch (err) { self.postMessage({ type: 'error', message: err.message + (err.stack ? '\n' + err.stack : '') }); }
};

async function run(file, enzymes, bedList, binSize) {
  // Build normalized centromere lookup; track which BED entries actually match.
  const bedMap = new Map();
  for (const b of bedList) bedMap.set(normChrom(b.name), { name: b.name, intervals: b.intervals, used: false });
  const bed = (fastaName) => { const e = bedMap.get(normChrom(fastaName)); if (e) { e.used = true; return e.intervals; } return []; };
  // Pre-compile per enzyme; skip the reverse scan for palindromes.
  const enz = enzymes.map(en => {
    const rc = revComp(en.motif);
    return { name: en.name, motif: en.motif, len: en.motif.length,
             fwd: toRegex(en.motif), rev: toRegex(rc), palindrome: rc === en.motif };
  });

  // Genome-wide accumulators, per enzyme.
  const acc = enz.map(() => ({
    totalSites: 0, cenSites: 0, armSites: 0, cenSize: 0, armSize: 0,
    cenSum: 0, cenSumSq: 0, cenN: 0, armSum: 0, armSumSq: 0, armN: 0,
    cenHist: new Float64Array(MONO_NBINS), armHist: new Float64Array(MONO_NBINS),
  }));
  const perChromRows = [];           // one row per (enzyme, chrom) — CSV parity
  const tracks = {};                 // tracks[chrom] = { len, binSize, byEnzyme: {name: Int32Array} }

  const total = file.size || 0;
  let readBytes = 0, lastPct = -1;

  let stream = file.stream();
  if (/\.gz$/i.test(file.name)) stream = stream.pipeThrough(new DecompressionStream('gzip'));
  const reader = stream.pipeThrough(new TextDecoderStream()).getReader();

  let curName = null, curParts = [], partial = '';

  const flush = () => {
    if (curName === null) return;
    processChrom(curName, curParts.join('').toUpperCase(), enz, acc, bed(curName), binSize, perChromRows, tracks);
    curParts = [];
  };

  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    readBytes += value.length; // approx (chars, pre-gzip differs — good enough for a bar)
    const pct = total ? Math.min(99, Math.floor(readBytes / total * 100)) : 0;
    if (pct !== lastPct) { lastPct = pct; self.postMessage({ type: 'progress', pct, chrom: curName }); }

    let text = partial + value;
    let nl, start = 0;
    while ((nl = text.indexOf('\n', start)) !== -1) {
      let line = text.slice(start, nl);
      if (line.endsWith('\r')) line = line.slice(0, -1);
      start = nl + 1;
      if (line[0] === '>') { flush(); curName = line.slice(1).split(/\s+/)[0]; }
      else if (curName !== null && line) curParts.push(line);
    }
    partial = text.slice(start);
  }
  // trailing line without newline
  if (partial) {
    let line = partial.endsWith('\r') ? partial.slice(0, -1) : partial;
    if (line[0] === '>') { flush(); curName = line.slice(1).split(/\s+/)[0]; }
    else if (curName !== null && line) curParts.push(line);
  }
  flush();

  // Genome-wide summary + recommendation-ready metrics.
  const eps = 1e-6;
  const summary = enz.map((en, i) => {
    const a = acc[i];
    const cenDensity = a.cenSize > 0 ? a.cenSites / a.cenSize * 1000 : 0;
    const armDensity = a.armSize > 0 ? a.armSites / a.armSize * 1000 : 0;
    const cenMean = a.cenN > 0 ? a.cenSum / a.cenN : 0;
    const cenStd  = a.cenN > 0 ? Math.sqrt(Math.max(0, a.cenSumSq / a.cenN - cenMean * cenMean)) : 0;
    const armMean = a.armN > 0 ? a.armSum / a.armN : 0;
    const armStd  = a.armN > 0 ? Math.sqrt(Math.max(0, a.armSumSq / a.armN - armMean * armMean)) : 0;
    const cenHom = cenMean > 0 ? cenMean / (cenStd + eps) : NaN;
    const armHom = armMean > 0 ? armMean / (armStd + eps) : NaN;
    return {
      name: en.name, motif: en.motif, baseLen: en.len, palindrome: en.palindrome,
      totalSites: a.totalSites, cenSites: a.cenSites, armSites: a.armSites,
      cenDensity, armDensity, densityDiff: cenDensity - armDensity,
      ratio: armDensity > 0 ? cenDensity / armDensity : Infinity,
      cenMean, cenStd, armMean, armStd,
      cenHom, armHom, homDiff: cenHom - armHom,
      cenHist: Array.from(a.cenHist), armHist: Array.from(a.armHist),
    };
  });

  const unusedBed = [...bedMap.values()].filter(e => !e.used).map(e => e.name);
  self.postMessage({
    type: 'done', summary, perChromRows, tracks, unusedBed,
    monoMeta: { max: MONO_MAX, step: MONO_STEP, nbins: MONO_NBINS },
  });
}

function processChrom(name, seq, enz, acc, intervals, binSize, perChromRows, tracks) {
  const len = seq.length;
  const cenSize = intervals.reduce((s, [a, b]) => s + (b - a), 0);
  const armSize = len - cenSize;
  const nbins = Math.max(1, Math.ceil(len / binSize));
  tracks[name] = { len, binSize, byEnzyme: {} };
  const inCen = (p) => { for (const [a, b] of intervals) if (p >= a && p <= b) return true; return false; };

  for (let ei = 0; ei < enz.length; ei++) {
    const en = enz[ei], A = acc[ei];
    const positions = [];
    findPositions(seq, en.fwd, positions);
    if (!en.palindrome) findPositions(seq, en.rev, positions);
    // dedupe + sort
    let sites = positions.length ? Array.from(new Set(positions)).sort((x, y) => x - y) : [];

    const bins = new Int32Array(nbins);
    const cenSites = [], armSites = [];
    for (const p of sites) {
      bins[Math.floor(p / binSize)]++;
      if (inCen(p)) cenSites.push(p); else armSites.push(p);
    }
    tracks[name].byEnzyme[en.name] = bins;

    // spacings within each region (np.diff on the region's sorted positions)
    const cenSp = diffStats(cenSites), armSp = diffStats(armSites);

    // accumulate genome-wide
    A.totalSites += sites.length; A.cenSites += cenSites.length; A.armSites += armSites.length;
    A.cenSize += cenSize; A.armSize += armSize;
    A.cenSum += cenSp.sum; A.cenSumSq += cenSp.sumSq; A.cenN += cenSp.n;
    A.armSum += armSp.sum; A.armSumSq += armSp.sumSq; A.armN += armSp.n;
    for (let k = 1; k < cenSites.length; k++) A.cenHist[histIndex(cenSites[k] - cenSites[k-1])]++;
    for (let k = 1; k < armSites.length; k++) A.armHist[histIndex(armSites[k] - armSites[k-1])]++;

    // per-chrom row (CSV parity, matches the Python column semantics)
    const cenDensity = cenSize > 0 ? cenSites.length / cenSize * 1000 : 0;
    const armDensity = armSize > 0 ? armSites.length / armSize * 1000 : 0;
    const eps = 1e-6;
    const cenHom = cenSp.mean > 0 ? cenSp.mean / (cenSp.std + eps) : NaN;
    const armHom = armSp.mean > 0 ? armSp.mean / (armSp.std + eps) : NaN;
    perChromRows.push({
      Chromosome: name, Enzyme: en.name, 'Enzyme Recognition Sequence': en.motif,
      'Total Sites': sites.length, 'Centromere Sites': cenSites.length, 'Arm Sites': armSites.length,
      'Centromere Density (cuts/kb)': cenDensity, 'Arm Density (cuts/kb)': armDensity,
      'Density Difference (Cen - Arm)': cenDensity - armDensity,
      'Centromere/Arm Ratio': armDensity > 0 ? cenDensity / armDensity : Infinity,
      'Centromere Mean Spacing (bp)': cenSp.mean, 'Centromere Spacing Std Dev': cenSp.std,
      'Arm Mean Spacing (bp)': armSp.mean, 'Arm Spacing Std Dev': armSp.std,
      'Centromere Homogeneity': cenHom, 'Arm Homogeneity': armHom,
      'Homogeneity Difference (Cen - Arm)': cenHom - armHom,
    });
  }
}

function diffStats(sorted) {
  let sum = 0, sumSq = 0, n = 0;
  for (let i = 1; i < sorted.length; i++) { const d = sorted[i] - sorted[i-1]; sum += d; sumSq += d*d; n++; }
  const mean = n > 0 ? sum / n : 0;
  const std = n > 0 ? Math.sqrt(Math.max(0, sumSq / n - mean * mean)) : 0; // population std (ddof=0), like np.std
  return { sum, sumSq, n, mean, std };
}
