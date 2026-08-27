# Pore-C Enzyme Explorer (web app)

> **A joint contribution by [Jacob Gonzalez Isa](https://github.com/jacgonisa)** (jg2070@cam.ac.uk), **[Katie Jenike](https://github.com/kjenike)** (kj436@cam.ac.uk), **and Hanwen Cao** (hc617@cam.ac.uk) — with special thanks to Katie Jenike for the design of the cut-density track.

A zero-backend browser app for **exploring** restriction-enzyme choice for
Pore-C. Drop in a **genome FASTA** + a **centromere BED** and inspect, for any
NEB enzyme, its cut density, predicted fragment (monomer) lengths, and
centromere-vs-arm balance — genome-wide and **per chromosome**. Everything runs
client-side, so **no data leaves the browser**.

It's exploratory by design: no "best enzyme" is prescribed. You get the metrics
and plots; you decide.

Live: **https://jacgonisa.github.io/porec_hendersonlab/webapp/**

## Use it
1. Upload a genome FASTA (`.fa`/`.fasta`/`.fna`, plain or gzipped) and a
   centromere BED (`chrom  start  end`), or click **Load Arabidopsis example**
   (bundled Col-0 genome + centromeres).
2. Pick enzymes: the **benchmarked 9** are ticked by default; **Select all NEB**
   screens the whole catalog (236 enzymes / 204 unique recognition sites); or
   search/add custom IUPAC motifs.
3. **Analyse**, then explore:
   - **Metrics table** — sortable + filterable; every enzyme, no ranking imposed.
     Download the full per-chromosome table as CSV.
   - **Balance map** — density-difference vs homogeneity-difference per enzyme,
     genome-wide or for a chosen chromosome. `x=0` = equal cutting in centromere
     and arms.
   - **Cut-density track** — per chromosome, per enzyme, at a **selectable bin
     size (200 bp – 10 kb)**, centromere shaded.
   - **Monomer-length histogram** — predicted fragment sizes, centromere vs arm,
     genome-wide or per chromosome.

## Enzyme set
NEB-supplied enzymes pulled from REBASE via Biopython (`Bio.Restriction`),
grouped by recognition site (isoschizomers share a row; hover the table row to
see them). Sites shorter than 4 bp are excluded. Regenerate with
`python3 make_enzymes.py` if you want to refresh the catalog. Methylation-
dependent enzymes (DpnI, MspJI, …) are shown by recognition sequence only —
in-silico counts don't reflect their real methylation-gated cutting.

## What it computes / fixes
Per enzyme, both strands are scanned (non-overlapping, like Python
`re.finditer`) and **palindromic motifs are counted once** — the original lab
scripts scan motif + reverse-complement and concatenate, double-counting
palindromes (GATC, CATG, CTAG, AATT …) and inflating their density vs
non-palindromic enzymes. Verified: non-palindromic enzymes (AlwI, BbsI …) match
the ground-truth `enzyme_choice/validation/correlation_exp_insilico.csv` to 3
decimals; palindromes come out at exactly ½ the old (buggy) counts.

Metrics = cut density (cuts/kb) for centromere and arms, their difference and
ratio, mean/std spacing (≈ monomer length), and homogeneity (mean/std spacing).
In-silico spacing correlates with experimental monomer length at R²≈0.86.

## Performance & limits
- Practical for plant-sized genomes (largest chromosome ≲ 200 Mb). A full
  204-site NEB screen of the ~125 Mb Arabidopsis genome takes ~40–90 s in the
  browser (progress bar shown); the benchmarked-9 example run is a few seconds.
- Screening the full NEB set keeps ~150 MB of cut positions in memory (so tracks
  can re-bin at any resolution) — fine on a desktop; screen fewer enzymes on
  low-RAM machines.
- Mammalian 3 GB genomes are impractical to upload in-browser — out of scope;
  use the Python scripts in `enzyme_choice/`.

## Hosting
Plain static files, no build step — **GitHub Pages hosts it fine** (the bundled
34 MB gzipped genome is well under the 100 MB file limit; all compute is
client-side). Netlify works too (drag the `webapp/` folder onto app.netlify.com).
Pages is configured to deploy from `main` / root, so the app lives at
`…github.io/porec_hendersonlab/webapp/`.

## Dev
- Run locally: `python3 -m http.server` here, open `http://localhost:8000/`
  (a `file://` open won't work — Workers and `fetch` need http).
- Logic test: `node webapp/test.mjs` (self-contained).
- Manual smoke test: open `index.html#selftest` — auto-loads the example and runs.

## Files
`index.html` · `app.js` (UI, table, plots, CSV) · `worker.js` (streaming scan →
per-enzyme cut positions + metrics) · `style.css` · `neb_enzymes.json` ·
`make_enzymes.py` · `example/` (Arabidopsis genome + BED) · `test.mjs`.
