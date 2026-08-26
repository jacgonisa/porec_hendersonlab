# Pore-C Enzyme Chooser (web app)

A zero-backend browser app: drop in a **genome FASTA** + a **centromere BED**,
get an in-silico restriction-enzyme ranking, metrics, and plots for a Pore-C
experiment. Everything runs client-side — **no data leaves the browser**, so
unpublished genomes are safe.

It ports [`enzyme_choice/table_with_mean_and_sd_ratio.py`](../enzyme_choice/)
to JavaScript. The goal is the same one the lab benchmark proved (in-silico cut
spacing correlates with experimental monomer length at R²≈0.86): pick the enzyme
that digests centromeres and arms most **evenly**, avoiding DpnII-style
over-digestion of AT-rich centromeres.

## Use it
1. Open the hosted page (see *Deploy* below) or `index.html` via a local server.
2. Upload a genome FASTA (`.fa`, `.fasta`, `.fna`, or gzipped `.fa.gz`) and a
   centromere BED (`chrom  start  end`).
3. Tick the enzymes to test (add custom IUPAC motifs if you like) → **Analyse**.
4. Read the recommendation, sortable metrics table, balance map, per-chromosome
   cut-density track, and predicted monomer-length histograms. Download the full
   per-chromosome table as CSV.

Or click **Load example data** to run a bundled synthetic genome
(`example/demo.fa` + `example/demo.bed`) where an AT-cutter (MluCI) visibly
over-digests the AT-rich centromere.

## Inputs & robustness
- **FASTA** streamed and scanned in a Web Worker, so a ~130 MB plant genome
  stays responsive. gzip is handled natively (`DecompressionStream`).
- **BED** accepts standard `chrom start end` and the lab's
  `id  name  start  end` variant (start/end are auto-detected).
- **Chromosome names** are matched leniently: `Chr1:1-32640075`, `chr1` and `1`
  all resolve to the same chromosome. Any BED entry that matches no FASTA
  sequence is reported as a warning.

### What it fixes vs the Python scripts
Palindromic motifs (GATC, CATG, CTAG, AATT, …) are counted **once**, not twice.
The original scripts scan motif + reverse-complement and concatenate, so a
palindrome (whose revcomp is itself) is double-counted — which also inflated its
density relative to non-palindromic enzymes. The app dedupes positions, giving a
fair cross-enzyme comparison. (So palindrome counts here are ~½ of the old CSVs;
non-palindromic enzymes like AlwI/BbsI are unchanged — verified identical to the
ground-truth CSV to 3 decimals.)

### Suitability score (transparent)
`0.5·balance + 0.3·fragment-size + 0.2·uniformity`, all in [0,1]:
- **balance** — small |CEN−arm density difference| (near 0 = even cutting).
- **fragment-size** — mean monomer length in the nanopore sweet spot ~0.8–1.4 kb.
- **uniformity** — high, and similar, homogeneity (mean/std spacing) in both regions.
- enzymes with too few centromere/arm cuts are penalised.

Every input metric is shown in the table, so the ranking is auditable.

## Limits
- Practical for plant-sized genomes (largest chromosome ≲ 200 Mb). Mammalian
  3 GB genomes are impractical to upload in-browser — out of scope; use the
  Python scripts for those.

## Deploy (GitHub Pages)
The app is plain static files — no build step.
1. Push this `webapp/` folder to the repo.
2. Repo **Settings → Pages** → *Deploy from a branch* → branch `main`, folder
   `/webapp` isn't directly selectable, so either (a) set Pages source to
   `/ (root)` and move these files to the repo root, **or** (b) enable Pages on
   the branch root and link to `…/webapp/`, **or** (c) use a `docs/` folder.
   Simplest: enable Pages (root), then the app lives at
   `https://<user>.github.io/porec_hendersonlab/webapp/`.
3. Netlify alternative: drag the `webapp/` folder onto app.netlify.com — done.

## Dev
- Run locally: `python3 -m http.server` in this folder, open
  `http://localhost:8000/`. (A `file://` open won't work — Workers and `fetch`
  need http.)
- Logic test: `node webapp/test.mjs` (uses the committed example genome).
- Manual smoke test: open `index.html#selftest` — it auto-loads the example and
  runs, and sets `document.body.dataset.done` to the top enzyme when finished.

## Files
`index.html` · `app.js` (UI, scoring, plots, CSV) · `worker.js` (streaming scan
+ metrics) · `style.css` · `example/` · `test.mjs`.
