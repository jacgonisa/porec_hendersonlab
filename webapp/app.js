'use strict';

// Benchmarked panel (Henderson-lab set). AlwI is the reference "balanced" enzyme.
const PRESETS = [
  { name: 'DpnII', motif: 'GATC' },   { name: 'NlaIII', motif: 'CATG' },  { name: 'AlwI', motif: 'GGATC' },
  { name: 'BfaI', motif: 'CTAG' },    { name: 'MluCI', motif: 'AATT' },   { name: 'BbsI', motif: 'GAAGAC' },
  { name: 'BglII', motif: 'AGATCT' }, { name: 'AflII', motif: 'CTTAAG' }, { name: 'MluI', motif: 'ACGCGT' },
];
const BIN_SIZE = 10000; // ponytail: fixed 10 kb density bins; expose a slider only if asked.

const $ = (id) => document.getElementById(id);
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const normChrom = (s) => String(s).toLowerCase().split(/[:\s]/)[0].replace(/^chr(omosome)?[_\-]?/, '');
const fmt = (v, d = 2) => (v == null || Number.isNaN(v)) ? '–' : (!Number.isFinite(v) ? '∞' : Number(v).toFixed(d));

let fastaFile = null, bedList = null, lastRun = null;

// ---- enzyme panel ----------------------------------------------------------
$('enzpanel').innerHTML = PRESETS.map((e, i) =>
  `<label class="enz"><input type="checkbox" data-i="${i}" checked> ${e.name} <code>${e.motif}</code></label>`).join('');

function selectedEnzymes() {
  const list = [...document.querySelectorAll('#enzpanel input:checked')].map(c => PRESETS[+c.dataset.i]);
  const custom = $('custom').value.split(',').map(s => s.trim().toUpperCase()).filter(Boolean);
  for (const m of custom) if (/^[ACGTRYSWKMBDHVN]+$/.test(m)) list.push({ name: m, motif: m });
  return list;
}

// ---- inputs ----------------------------------------------------------------
$('fasta').onchange = (e) => { fastaFile = e.target.files[0] || null; refresh(); };
$('bed').onchange = async (e) => { const f = e.target.files[0]; bedList = f ? parseBed(await f.text()) : null; refresh(); };

$('example').onclick = async () => {
  try {
    const [fa, bd] = await Promise.all([fetch('example/demo.fa'), fetch('example/demo.bed')]);
    fastaFile = new File([await fa.blob()], 'demo.fa');
    bedList = parseBed(await bd.text());
    refresh(); $('filestatus').textContent += '  (example loaded)';
  } catch (err) { alert('Could not load example: ' + err.message); }
};

function refresh() {
  const parts = [];
  if (fastaFile) parts.push(`FASTA: ${fastaFile.name} (${(fastaFile.size / 1e6).toFixed(1)} MB)`);
  if (bedList) parts.push(`BED: ${bedList.length} chromosome(s), ${bedList.reduce((n, b) => n + b.intervals.length, 0)} interval(s)`);
  $('filestatus').textContent = parts.join('   ·   ');
  $('run').disabled = !(fastaFile && bedList && bedList.length);
}

// Accepts standard BED (chrom start end) and the lab's id/name/start/end variant:
// use the field just before the first pair of adjacent integer columns as the chrom.
function parseBed(text) {
  const out = {};
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line[0] === '#' || /^(track|browser)\b/.test(line)) continue;
    const f = line.split(/\s+/);
    let chrom = null, start, end;
    for (let i = 1; i < f.length - 1; i++) {
      if (/^\d+$/.test(f[i]) && /^\d+$/.test(f[i + 1])) { chrom = f[i - 1]; start = +f[i]; end = +f[i + 1]; break; }
    }
    if (chrom === null || end <= start) continue;
    (out[chrom] = out[chrom] || []).push([start, end]);
  }
  return Object.entries(out).map(([name, intervals]) => ({ name, intervals }));
}

// ---- run -------------------------------------------------------------------
$('run').onclick = () => {
  const enzymes = selectedEnzymes();
  if (!enzymes.length) { alert('Pick at least one enzyme.'); return; }
  $('run').disabled = true; $('warn').classList.add('hidden');
  $('progress').classList.remove('hidden'); setBar(0, 'starting…');

  const worker = new Worker('worker.js');
  worker.onmessage = (ev) => {
    const d = ev.data;
    if (d.type === 'progress') setBar(d.pct, `reading ${d.chrom || ''}… ${d.pct}%`);
    else if (d.type === 'error') { setBar(0, ''); $('progress').classList.add('hidden'); $('run').disabled = false; showWarn('Error: ' + d.message); worker.terminate(); }
    else if (d.type === 'done') { setBar(100, 'done'); setTimeout(() => $('progress').classList.add('hidden'), 400); $('run').disabled = false; render(d, enzymes); worker.terminate(); }
  };
  worker.postMessage({ file: fastaFile, enzymes, bedList, binSize: BIN_SIZE });
};

function setBar(pct, msg) { $('bar').style.width = pct + '%'; $('pmsg').textContent = msg; }
function showWarn(msg) { $('warn').innerHTML = msg; $('warn').classList.remove('hidden'); }

// ---- scoring ---------------------------------------------------------------
function score(s) {
  const sDD = clamp(1 - Math.abs(s.densityDiff) / 2, 0, 1);          // near-0 density diff is best
  const lenScore = (v) => v <= 0 ? 0 : clamp(v < 800 ? (v - 300) / 500 : v <= 1400 ? 1 : (3000 - v) / 1600, 0, 1);
  const sLen = (lenScore(s.cenMean) + lenScore(s.armMean)) / 2;      // ~0.8–1.4 kb monomers ideal
  const hmin = Math.min(s.cenHom || 0, s.armHom || 0);
  const sHom = clamp(hmin / 0.8, 0, 1);                              // uniform spacing in both regions
  const enough = (s.cenSites >= 200 && s.armSites >= 200) ? 1 : 0.3; // too few cuts = unusable
  return { total: (0.5 * sDD + 0.3 * sLen + 0.2 * sHom) * enough, sDD, sLen, sHom, enough };
}

// ---- render ----------------------------------------------------------------
function render(d, enzymes) {
  lastRun = d;
  const ranked = d.summary.map(s => ({ ...s, s: score(s) })).sort((a, b) => b.s.total - a.s.total);
  d.ranked = ranked;

  if (d.unusedBed && d.unusedBed.length)
    showWarn(`Heads up: these BED chromosomes matched no FASTA sequence (naming mismatch — they were treated as absent): <b>${d.unusedBed.join(', ')}</b>`);
  if (ranked.every(r => r.totalSites === 0)) { showWarn('No cut sites found — check that the FASTA and enzyme motifs are valid.'); return; }

  $('results').classList.remove('hidden');
  renderReco(ranked[0]);
  renderTable(ranked);
  renderScatter(ranked);
  fillSelect($('trackChrom'), Object.keys(d.tracks));
  fillSelect($('trackEnz'), ranked.map(r => r.name));
  fillSelect($('monoEnz'), ranked.map(r => r.name));
  $('trackChrom').onchange = $('trackEnz').onchange = () => renderTrack(d);
  $('monoEnz').onchange = () => renderMono(d);
  renderTrack(d); renderMono(d);
}

function fillSelect(sel, opts) { sel.innerHTML = opts.map(o => `<option>${o}</option>`).join(''); }

function renderReco(r) {
  const cells = [
    ['Motif', r.motif + (r.palindrome ? ' (palindrome)' : '')],
    ['Suitability score', fmt(r.s.total, 2) + ' / 1.00'],
    ['Density diff (CEN−arm)', fmt(r.densityDiff, 2) + ' cuts/kb'],
    ['CEN monomer', fmt(r.cenMean, 0) + ' bp'],
    ['Arm monomer', fmt(r.armMean, 0) + ' bp'],
    ['CEN cut sites', r.cenSites.toLocaleString()],
  ];
  const why = Math.abs(r.densityDiff) < 0.6
    ? 'balanced cutting across centromere and arms'
    : (r.densityDiff > 0 ? 'best available, but leans toward over-cutting the centromere' : 'best available, but under-cuts the centromere');
  $('reco').innerHTML =
    `<div class="reco-name">${r.name}</div><p class="hint">Top pick — ${why}.</p>` +
    `<div class="reco-grid">${cells.map(([k, v]) => `<div><b>${k}</b>${v}</div>`).join('')}</div>`;
}

// [key, header, accessor, kind]  kind: 'str' | 'int' | number-of-decimals
const COLS = [
  ['name', 'Enzyme', s => s.name, 'str'], ['motif', 'Motif', s => s.motif, 'str'],
  ['score', 'Score', s => s.s.total, 2], ['cenDensity', 'CEN cuts/kb', s => s.cenDensity, 2],
  ['armDensity', 'Arm cuts/kb', s => s.armDensity, 2], ['densityDiff', 'Δ density', s => s.densityDiff, 2],
  ['ratio', 'CEN/arm', s => s.ratio, 2], ['cenMean', 'CEN monomer', s => s.cenMean, 0],
  ['armMean', 'Arm monomer', s => s.armMean, 0], ['cenHom', 'CEN homog.', s => s.cenHom, 2],
  ['armHom', 'Arm homog.', s => s.armHom, 2], ['totalSites', 'Total sites', s => s.totalSites, 'int'],
];
let sortKey = 'score', sortDir = -1;
function renderTable(ranked) {
  const rows = [...ranked].sort((a, b) => {
    const c = COLS.find(c => c[0] === sortKey)[2];
    const x = c(a), y = c(b);
    return (typeof x === 'string' ? String(x).localeCompare(String(y)) : (x - y)) * sortDir;
  });
  const cell = (v, kind) => kind === 'str' ? v : kind === 'int' ? v.toLocaleString() : fmt(v, kind);
  const head = COLS.map(c => `<th data-k="${c[0]}">${c[1]}${sortKey === c[0] ? (sortDir < 0 ? ' ▾' : ' ▴') : ''}</th>`).join('');
  const body = rows.map(s => {
    const best = s === ranked[0] ? ' class="best"' : '';
    return `<tr${best}>${COLS.map(c => `<td>${cell(c[2](s), c[3])}</td>`).join('')}</tr>`;
  }).join('');
  $('ranktable').innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
  $('ranktable').querySelectorAll('th').forEach(th => th.onclick = () => {
    const k = th.dataset.k; if (k === sortKey) sortDir *= -1; else { sortKey = k; sortDir = k === 'name' || k === 'motif' ? 1 : -1; }
    renderTable(ranked);
  });
}

function renderScatter(ranked) {
  const t = {
    x: ranked.map(r => r.densityDiff), y: ranked.map(r => r.homDiff),
    text: ranked.map(r => r.name), mode: 'markers+text', textposition: 'top center', type: 'scatter',
    marker: { size: ranked.map(r => 10 + 26 * r.s.total), color: ranked.map(r => r.s.total), colorscale: 'YlGnBu', showscale: true, colorbar: { title: 'score' }, line: { width: 1, color: '#345' } },
    hovertemplate: '%{text}<br>Δdensity %{x:.2f}<br>Δhomog %{y:.2f}<extra></extra>',
  };
  Plotly.newPlot('scatter', [t], {
    margin: { t: 10, r: 10 }, xaxis: { title: 'Density difference (CEN − arm, cuts/kb)', zeroline: true, zerolinecolor: '#2b6cb0', zerolinewidth: 2 },
    yaxis: { title: 'Homogeneity difference (CEN − arm)' }, shapes: [{ type: 'line', x0: 0, x1: 0, yref: 'paper', y0: 0, y1: 1, line: { color: '#2b6cb0', dash: 'dot' } }],
  }, PCFG);
}

function renderTrack(d) {
  const chrom = $('trackChrom').value, enz = $('trackEnz').value;
  const tr = d.tracks[chrom]; if (!tr) return;
  const bins = tr.byEnzyme[enz]; const x = [], y = [];
  for (let i = 0; i < bins.length; i++) { x.push(i * tr.binSize / 1e6); y.push(bins[i] / (tr.binSize / 1000)); } // cuts/kb
  const shapes = (bedIntervalsFor(chrom, d)).map(([a, b]) => ({ type: 'rect', xref: 'x', yref: 'paper', x0: a / 1e6, x1: b / 1e6, y0: 0, y1: 1, fillcolor: 'rgba(197,48,48,.12)', line: { width: 0 } }));
  Plotly.newPlot('track', [{ x, y, type: 'scattergl', mode: 'lines', line: { width: 1, color: '#2b6cb0' }, fill: 'tozeroy', name: enz }], {
    margin: { t: 10, r: 10 }, xaxis: { title: `${chrom} position (Mb)` }, yaxis: { title: 'cut density (cuts/kb)' },
    shapes, annotations: shapes.length ? [{ x: (shapes[0].x0 + shapes[0].x1) / 2, y: 1, yref: 'paper', text: 'centromere', showarrow: false, font: { color: '#c53030', size: 11 } }] : [],
  }, PCFG);
}
function bedIntervalsFor(fastaChrom, d) {
  const b = bedList.find(x => normChrom(x.name) === normChrom(fastaChrom));
  return b ? b.intervals : [];
}

function renderMono(d) {
  const enz = $('monoEnz').value; const s = d.summary.find(x => x.name === enz);
  const { step, nbins } = d.monoMeta;
  const centers = Array.from({ length: nbins }, (_, i) => i === nbins - 1 ? (nbins - 1) * step : i * step + step / 2);
  const bar = (h, name, color) => ({ x: centers, y: h, type: 'bar', name, marker: { color }, opacity: 0.6 });
  Plotly.newPlot('mono', [bar(s.armHist, 'arms', '#2b6cb0'), bar(s.cenHist, 'centromere', '#c53030')], {
    barmode: 'overlay', margin: { t: 10, r: 10 }, xaxis: { title: 'monomer length (bp), last bin = ≥5000' }, yaxis: { title: 'count' },
    shapes: [800, 1400].map(v => ({ type: 'line', x0: v, x1: v, yref: 'paper', y0: 0, y1: 1, line: { color: '#2f855a', dash: 'dot' } })),
  }, PCFG);
}

const PCFG = { responsive: true, displaylogo: false, toImageButtonOptions: { scale: 2 } };

// ---- CSV -------------------------------------------------------------------
$('csv').onclick = () => {
  if (!lastRun) return;
  const rows = lastRun.perChromRows; if (!rows.length) return;
  const cols = Object.keys(rows[0]);
  const esc = (v) => typeof v === 'number' ? (Number.isFinite(v) ? v : (Number.isNaN(v) ? '' : 'inf')) : `"${String(v).replace(/"/g, '""')}"`;
  const csv = [cols.join(','), ...rows.map(r => cols.map(c => esc(r[c])).join(','))].join('\n');
  const url = URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
  const a = document.createElement('a'); a.href = url; a.download = 'porec_enzyme_metrics.csv'; a.click(); URL.revokeObjectURL(url);
};

// Smoke test: open index.html#selftest to auto-load the example and run.
const _render = render;
render = function (d, e) { _render(d, e); document.body.dataset.done = d.ranked ? d.ranked[0].name : 'none'; };
if (location.hash === '#selftest') window.addEventListener('load', async () => {
  await $('example').onclick();
  $('run').click();
});
