'use strict';

const APP_VERSION = 4; // bump whenever worker.js changes, to bypass browser worker cache
const $ = (id) => document.getElementById(id);
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const normChrom = (s) => String(s).toLowerCase().split(/[:\s]/)[0].replace(/^chr(omosome)?[_\-]?/, '');
const fmt = (v, d = 2) => (v == null || Number.isNaN(v)) ? '–' : (!Number.isFinite(v) ? '∞' : Number(v).toFixed(d));
const PCFG = { responsive: true, displaylogo: false, toImageButtonOptions: { scale: 2 } };
const MONO_MAX = 5000, MONO_STEP = 50, MONO_NBINS = MONO_MAX / MONO_STEP + 1;

let ENZYMES = [];   // [{site, len, enzymes:[names], bench:[names], name(rep), search}]
let fastaFile = null, bedList = null, lastRun = null;

// ---- load NEB enzyme catalog + build picker --------------------------------
fetch('neb_enzymes.json').then(r => r.json()).then(list => {
  ENZYMES = list.map(e => {
    const name = (e.bench[0] || e.enzymes[0]);
    return { ...e, name, search: (e.site + ' ' + e.enzymes.join(' ') + ' ' + e.len + 'bp').toLowerCase() };
  });
  buildEnzList();
  selectSites(new Set(ENZYMES.filter(e => e.bench.length).map(e => e.site))); // benchmarked default
});

function buildEnzList() {
  $('enzlist').innerHTML = ENZYMES.map((e, i) => {
    const iso = e.enzymes.length > 1 ? ` <small>+${e.enzymes.length - 1}</small>` : '';
    const star = e.bench.length ? '<span class="star" title="benchmarked">★</span> ' : '';
    return `<label class="enz" data-i="${i}"><input type="checkbox" data-site="${e.site}"> ${star}${e.name} <code>${e.site}</code>${iso}</label>`;
  }).join('');
  $('enzlist').oninput = updateCount;
  $('enzsearch').oninput = filterEnzList;
  $('selBench').onclick = () => selectSites(new Set(ENZYMES.filter(e => e.bench.length).map(e => e.site)));
  $('selAll').onclick = () => document.querySelectorAll('#enzlist input').forEach(c => c.checked = true) || updateCount();
  $('selNone').onclick = () => document.querySelectorAll('#enzlist input').forEach(c => c.checked = false) || updateCount();
}
function filterEnzList() {
  const q = $('enzsearch').value.toLowerCase().trim();
  document.querySelectorAll('#enzlist .enz').forEach(el => {
    el.style.display = (!q || ENZYMES[+el.dataset.i].search.includes(q)) ? '' : 'none';
  });
}
function selectSites(set) {
  document.querySelectorAll('#enzlist input').forEach(c => c.checked = set.has(c.dataset.site));
  updateCount();
}
function selectedEnzymes() {
  const seen = new Set(), out = [];
  document.querySelectorAll('#enzlist input:checked').forEach(c => {
    const e = ENZYMES.find(x => x.site === c.dataset.site);
    if (e && !seen.has(e.site)) { seen.add(e.site); out.push({ name: e.name, motif: e.site, enzymes: e.enzymes }); }
  });
  for (const m of $('custom').value.split(',').map(s => s.trim().toUpperCase()).filter(Boolean))
    if (/^[ACGTRYSWKMBDHVN]+$/.test(m) && !seen.has(m)) { seen.add(m); out.push({ name: m, motif: m, enzymes: [m] }); }
  return out;
}
function updateCount() {
  const n = selectedEnzymes().length;
  $('enzcount').textContent = `${n} selected of ${ENZYMES.length} NEB sites`;
  $('runhint').textContent = n > 60 ? 'Screening many enzymes on a whole genome can take 1–2 min — progress shown below.' : '';
}

// ---- inputs ----------------------------------------------------------------
$('fasta').onchange = (e) => { fastaFile = e.target.files[0] || null; refresh(); };
$('bed').onchange = async (e) => { const f = e.target.files[0]; bedList = f ? parseBed(await f.text()) : null; refresh(); };
$('example').onclick = async () => {
  const btn = $('example'); btn.disabled = true; btn.textContent = 'Loading… (34 MB)';
  try {
    const [fa, bd] = await Promise.all([fetch('example/arabidopsis.fa.gz'), fetch('example/arabidopsis.bed')]);
    fastaFile = new File([await fa.blob()], 'arabidopsis.fa.gz');
    bedList = parseBed(await bd.text());
    refresh();
  } catch (err) { alert('Could not load example: ' + err.message); }
  btn.disabled = false; btn.textContent = 'Load Arabidopsis example';
};
function refresh() {
  const parts = [];
  if (fastaFile) parts.push(`FASTA: ${fastaFile.name} (${(fastaFile.size / 1e6).toFixed(1)} MB)`);
  if (bedList) parts.push(`BED: ${bedList.length} chromosome(s), ${bedList.reduce((n, b) => n + b.intervals.length, 0)} interval(s)`);
  $('filestatus').textContent = parts.join('   ·   ');
  $('run').disabled = !(fastaFile && bedList && bedList.length);
}
// Accepts standard BED and the lab's id/name/start/end variant.
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

// ---- run -------------------------------------------------------------------
$('run').onclick = () => {
  const enzymes = selectedEnzymes();
  if (!enzymes.length) { alert('Pick at least one enzyme.'); return; }
  $('run').disabled = true; $('warn').classList.add('hidden');
  $('progress').classList.remove('hidden'); setBar(0, 'starting…');
  const worker = new Worker('worker.js?v=' + APP_VERSION); // ?v busts stale worker cache
  worker.onmessage = (ev) => {
    const d = ev.data;
    if (d.type === 'progress') {
      if (d.phase === 'read') setBar(Math.floor(d.pct * 0.5), `reading genome… ${d.pct}%`);
      else setBar(50 + Math.floor((d.enzI / d.enzN) * 50), `scanning ${d.chrom}: enzyme ${d.enzI + 1}/${d.enzN}`);
    } else if (d.type === 'error') {
      $('progress').classList.add('hidden'); $('run').disabled = false; showWarn('Error: ' + d.message); worker.terminate();
    } else if (d.type === 'done') {
      setBar(100, 'done'); setTimeout(() => $('progress').classList.add('hidden'), 400); $('run').disabled = false;
      render(d, enzymes); worker.terminate();
    }
  };
  worker.postMessage({ file: fastaFile, enzymes, bedList });
};
function setBar(pct, msg) { $('bar').style.width = pct + '%'; $('pmsg').textContent = msg; }
function showWarn(msg) { $('warn').innerHTML = msg; $('warn').classList.remove('hidden'); }

// ---- render ----------------------------------------------------------------
function render(d, enzymes) {
  lastRun = d;
  d.enzByName = new Map(enzymes.map(e => [e.name, e]));
  d.ivByChrom = new Map(d.chromMeta.map(c => [c.name, c.intervals]));
  d.chroms = d.chromMeta.map(c => c.name);
  d.rowsByChrom = {};
  for (const r of d.perChromRows) (d.rowsByChrom[r.Chromosome] = d.rowsByChrom[r.Chromosome] || []).push(r);

  if (d.unusedBed && d.unusedBed.length)
    showWarn(`Heads up: these BED chromosomes matched no FASTA sequence (naming mismatch — treated as absent): <b>${d.unusedBed.join(', ')}</b>`);
  if (d.summary.every(r => r.totalSites === 0)) { showWarn('No cut sites found — check the FASTA and motifs.'); return; }

  $('results').classList.remove('hidden');
  renderTable(d);
  fillSelect($('scChrom'), ['Genome-wide', ...d.chroms]);
  fillSelect($('trackChrom'), d.chroms);
  fillSelect($('trackEnz'), d.summary.map(s => s.name));
  fillSelect($('monoChrom'), ['Genome-wide', ...d.chroms]);
  fillSelect($('monoEnz'), d.summary.map(s => s.name));
  $('scChrom').onchange = () => renderScatter(d);
  $('trackChrom').onchange = $('trackEnz').onchange = $('trackBin').onchange = () => renderTrack(d);
  $('monoChrom').onchange = $('monoEnz').onchange = () => renderMono(d);
  $('tablefilter').oninput = () => renderTable(d);
  renderScatter(d); renderTrack(d); renderMono(d);
}
function fillSelect(sel, opts) { const cur = sel.value; sel.innerHTML = opts.map(o => `<option>${o}</option>`).join(''); if (opts.includes(cur)) sel.value = cur; }

// metrics table (no recommendation / score) --------------------------------
const COLS = [
  ['name', 'Enzyme', s => s.name, 'str'], ['motif', 'Site', s => s.motif, 'str'],
  ['nIso', '#NEB', s => (lastRun.enzByName.get(s.name)?.enzymes.length || 1), 'int'],
  ['cenDensity', 'CEN cuts/kb', s => s.cenDensity, 2], ['armDensity', 'Arm cuts/kb', s => s.armDensity, 2],
  ['densityDiff', 'Δ density', s => s.densityDiff, 2], ['ratio', 'CEN/arm', s => s.ratio, 2],
  ['cenMean', 'CEN monomer', s => s.cenMean, 0], ['armMean', 'Arm monomer', s => s.armMean, 0],
  ['cenHom', 'CEN homog.', s => s.cenHom, 2], ['armHom', 'Arm homog.', s => s.armHom, 2],
  ['totalSites', 'Total sites', s => s.totalSites, 'int'],
];
let sortKey = 'totalSites', sortDir = -1;
function renderTable(d) {
  const q = $('tablefilter').value.toLowerCase().trim();
  let rows = d.summary;
  if (q) rows = rows.filter(s => (s.name + ' ' + s.motif + ' ' + (d.enzByName.get(s.name)?.enzymes.join(' ') || '')).toLowerCase().includes(q));
  const acc = COLS.find(c => c[0] === sortKey)[2];
  rows = [...rows].sort((a, b) => { const x = acc(a), y = acc(b); return (typeof x === 'string' ? String(x).localeCompare(String(y)) : (x - y)) * sortDir; });
  const cell = (v, kind) => kind === 'str' ? v : kind === 'int' ? Number(v).toLocaleString() : fmt(v, kind);
  const head = COLS.map(c => `<th data-k="${c[0]}">${c[1]}${sortKey === c[0] ? (sortDir < 0 ? ' ▾' : ' ▴') : ''}</th>`).join('');
  const body = rows.map(s => {
    const iso = d.enzByName.get(s.name)?.enzymes || [s.name];
    const title = iso.length > 1 ? ` title="${iso.join(', ')}"` : '';
    return `<tr${title}>${COLS.map(c => `<td>${cell(c[2](s), c[3])}</td>`).join('')}</tr>`;
  }).join('');
  $('ranktable').innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
  $('ranktable').querySelectorAll('th').forEach(th => th.onclick = () => {
    const k = th.dataset.k; if (k === sortKey) sortDir *= -1; else { sortKey = k; sortDir = (k === 'name' || k === 'motif') ? 1 : -1; }
    renderTable(d);
  });
}

// balance map (genome-wide or per chromosome) -------------------------------
function renderScatter(d) {
  const sel = $('scChrom').value;
  let pts;
  if (sel === 'Genome-wide') pts = d.summary.map(s => ({ name: s.name, dd: s.densityDiff, hd: s.homDiff, tot: s.totalSites }));
  else pts = (d.rowsByChrom[sel] || []).map(r => ({ name: r.Enzyme, dd: r['Density Difference (Cen - Arm)'], hd: r['Homogeneity Difference (Cen - Arm)'], tot: r['Total Sites'] }));
  pts = pts.filter(p => Number.isFinite(p.dd) && Number.isFinite(p.hd));
  const many = pts.length > 40;
  Plotly.newPlot('scatter', [{
    x: pts.map(p => p.dd), y: pts.map(p => p.hd), text: pts.map(p => p.name),
    mode: many ? 'markers' : 'markers+text', textposition: 'top center', type: 'scattergl',
    marker: { size: pts.map(p => 8 + 4 * Math.log10(1 + p.tot)), color: pts.map(p => p.dd), colorscale: 'RdBu', reversescale: true, cmid: 0, showscale: true, colorbar: { title: 'Δ dens' }, line: { width: 0.5, color: '#456' } },
    hovertemplate: '%{text}<br>Δdensity %{x:.2f} cuts/kb<br>Δhomog %{y:.2f}<extra></extra>',
  }], {
    margin: { t: 10, r: 10 }, xaxis: { title: 'Density difference (CEN − arm, cuts/kb)', zeroline: true, zerolinecolor: '#2b6cb0', zerolinewidth: 2 },
    yaxis: { title: 'Homogeneity difference (CEN − arm)' },
  }, PCFG);
}

// cut-density track (any bin size, per chromosome) --------------------------
function renderTrack(d) {
  const chrom = $('trackChrom').value, enz = $('trackEnz').value, bin = +$('trackBin').value;
  const pos = d.positions[chrom] && d.positions[chrom][enz]; if (!pos) return;
  const len = d.chromMeta.find(c => c.name === chrom).len;
  const nb = Math.max(1, Math.ceil(len / bin)); const counts = new Float64Array(nb);
  for (let i = 0; i < pos.length; i++) counts[Math.floor(pos[i] / bin)]++;
  const x = [], y = [], perKb = bin / 1000;
  for (let i = 0; i < nb; i++) { x.push(i * bin / 1e6); y.push(counts[i] / perKb); }
  const shapes = (d.ivByChrom.get(chrom) || []).map(([a, b]) => ({ type: 'rect', xref: 'x', yref: 'paper', x0: a / 1e6, x1: b / 1e6, y0: 0, y1: 1, fillcolor: 'rgba(197,48,48,.12)', line: { width: 0 } }));
  Plotly.newPlot('track', [{ x, y, type: 'scattergl', mode: 'lines', line: { width: 1, color: '#2b6cb0' }, fill: 'tozeroy', name: enz }], {
    margin: { t: 10, r: 10 }, xaxis: { title: `${chrom} position (Mb)` }, yaxis: { title: 'cut density (cuts/kb)' },
    shapes, annotations: shapes.length ? [{ x: (shapes[0].x0 + shapes[0].x1) / 2, y: 1, yref: 'paper', text: 'centromere', showarrow: false, font: { color: '#c53030', size: 11 } }] : [],
  }, PCFG);
}

// monomer-length histogram (genome-wide or per chromosome) ------------------
function monoHist(d, enz, chrom) {
  const cen = new Float64Array(MONO_NBINS), arm = new Float64Array(MONO_NBINS);
  const idx = v => v >= MONO_MAX ? MONO_NBINS - 1 : Math.floor(v / MONO_STEP);
  const chroms = chrom === 'Genome-wide' ? d.chroms : [chrom];
  for (const cn of chroms) {
    const pos = d.positions[cn] && d.positions[cn][enz]; if (!pos) continue;
    const iv = d.ivByChrom.get(cn) || [];
    const inCen = p => { for (const [a, b] of iv) if (p >= a && p <= b) return true; return false; };
    let lastC = null, lastA = null;
    for (let i = 0; i < pos.length; i++) {
      const p = pos[i];
      if (inCen(p)) { if (lastC !== null) cen[idx(p - lastC)]++; lastC = p; }
      else { if (lastA !== null) arm[idx(p - lastA)]++; lastA = p; }
    }
  }
  return { cen, arm };
}
function renderMono(d) {
  const enz = $('monoEnz').value, chrom = $('monoChrom').value;
  const h = monoHist(d, enz, chrom);
  const centers = Array.from({ length: MONO_NBINS }, (_, i) => i === MONO_NBINS - 1 ? (MONO_NBINS - 1) * MONO_STEP : i * MONO_STEP + MONO_STEP / 2);
  const bar = (y, name, color) => ({ x: centers, y: Array.from(y), type: 'bar', name, marker: { color }, opacity: 0.6 });
  Plotly.newPlot('mono', [bar(h.arm, 'arms', '#2b6cb0'), bar(h.cen, 'centromere', '#c53030')], {
    barmode: 'overlay', margin: { t: 10, r: 10 }, xaxis: { title: 'monomer length (bp), last bin = ≥5000' }, yaxis: { title: 'count' },
    shapes: [800, 1400].map(v => ({ type: 'line', x0: v, x1: v, yref: 'paper', y0: 0, y1: 1, line: { color: '#2f855a', dash: 'dot' } })),
  }, PCFG);
}

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

// Smoke test: index.html#selftest auto-loads the example and runs benchmarked enzymes.
if (location.hash === '#selftest') window.addEventListener('load', async () => {
  await $('example').onclick(); $('run').click();
});
