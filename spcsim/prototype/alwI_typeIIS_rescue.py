"""
AlwI is Type IIS: GGATC(4/5) cuts 4-5 nt 3' of the site, so each cut leaves ONE site-retaining end
and ONE site-less end. pore-c-py digests the READ by motif, so a junction where both ligated ends
are site-less has no GGATC to cut -> the two monomers stay fused. This quantifies the loss and shows
it is recoverable from the alignment (a fused monomer is a chimera -> minimap2 gives it a supplementary
alignment, whose breakpoint is the junction the digest missed).

Parts:
 A. digest real Col-0 Chr1 with AlwI -> fraction of fragment ends that are site-less, and the
    both-site-less junction fraction (= the digest-missed fraction).
 B. build synthetic fused monomers (two random AlwI fragments) -> minimap2 -> fraction recovered (SA).
 C. real-data grounding: SA-tag (chimera) fraction among monomers in the BfaI WT cs.bam.
Run: python alwI_typeIIS_rescue.py
"""
import re, os, random, subprocess, sys
import pysam

GENOME = '/home/jg2070/Desktop/PhD/PoreC/Col-0.ragtag_chrs.mito.chloro.fa'
CHR1 = 'Chr1:1-32640075'
MM2 = '/home/jg2070/minimap2/minimap2'
BFAI_BAM = '/home/jg2070/Desktop/PhD/PoreC/experiment_results/2025December/results_BfaI_WT/WT.cs.bam'
WORK = '/tmp/alwI_rescue'; os.makedirs(WORK, exist_ok=True)
random.seed(1)
REGION_BP = 8_000_000                      # analyse first 8 Mb of Chr1


def rc(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def alwI_cuts(seq):
    """return sorted list of (cut_coord, kind) where kind 'F'=forward GGATC (site on LEFT fragment),
    'R'=reverse GATCC (site on RIGHT fragment).  GGATC(4/5): cut ~4 nt 3' of the 5-bp site."""
    cuts = []
    for m in re.finditer('GGATC', seq):
        cuts.append((m.start() + 9, 'F'))          # cut 4 nt after the 5-bp site
    for m in re.finditer('GATCC', seq):            # reverse-strand site
        cuts.append((m.start() - 4, 'R'))
    cuts = sorted(set(cuts))
    return [(c, k) for c, k in cuts if 0 < c < len(seq)]


def part_A(seq):
    cuts = alwI_cuts(seq)
    # fragment i spans (cuts[i], cuts[i+1]); left end retained iff left cut is 'R', right end iff right cut is 'F'
    siteless = total = 0
    frags = []
    for i in range(len(cuts) - 1):
        (ca, ka), (cb, kb) = cuts[i], cuts[i + 1]
        left_retained = (ka == 'R'); right_retained = (kb == 'F')
        siteless += (not left_retained) + (not right_retained); total += 2
        frags.append((ca, cb))
    p = siteless / total
    print(f"A. AlwI digest of Chr1 (first {REGION_BP/1e6:.0f} Mb): {len(cuts)} cuts, {len(frags)} fragments")
    print(f"   site-less fragment ends : {100*p:.1f}%   (expected ~50%: one site-retaining + one site-less per cut)")
    print(f"   both-site-less junctions: {100*p*p:.1f}%   -> fraction the motif-digest CANNOT split\n")
    return frags, p


def part_B(seq, frags, n=2000):
    """build n fused monomers (two random fragments, random orientation), 300-2000 bp each, align, count SA."""
    usable = [(a, b) for a, b in frags if 300 <= b - a <= 2000]
    recs = []
    with open(f'{WORK}/ref.fa', 'w') as f:
        f.write(f'>chr1_region\n{seq}\n')
    with open(f'{WORK}/fused.fa', 'w') as f:
        for i in range(n):
            a1, b1 = random.choice(usable); a2, b2 = random.choice(usable)
            s1 = seq[a1:b1]; s2 = seq[a2:b2]
            if random.random() < 0.5: s1 = rc(s1)
            if random.random() < 0.5: s2 = rc(s2)
            f.write(f'>fused_{i}\n{s1}{s2}\n')          # a fused (under-digested) monomer
    out = subprocess.run([MM2, '-ax', 'map-ont', '--secondary=no', f'{WORK}/ref.fa', f'{WORK}/fused.fa'],
                         capture_output=True, text=True, cwd=WORK)
    with open(f'{WORK}/fused.sam', 'w') as f: f.write(out.stdout)
    aln = {}
    sf = pysam.AlignmentFile(f'{WORK}/fused.sam')
    for r in sf:
        if r.is_unmapped or r.is_secondary: continue
        aln.setdefault(r.query_name, 0)
        aln[r.query_name] += 1                          # primary + supplementary
    split = sum(1 for v in aln.values() if v >= 2)
    mapped = len(aln)
    print(f"B. recovery by alignment: of {mapped} fused monomers minimap2 mapped, "
          f"{100*split/mapped:.1f}% split into >=2 alignments (SA) -> junction recovered\n")
    return split / mapped if mapped else 0


def part_C(bam, n=200_000):
    if not os.path.exists(bam): print("C. BfaI BAM not found, skipping"); return
    sf = pysam.AlignmentFile(bam, 'rb')
    prim = sa = 0
    for r in sf:
        if r.is_secondary or r.is_supplementary: continue
        prim += 1
        if r.has_tag('SA'): sa += 1
        if prim >= n: break
    print(f"C. real-data grounding — BfaI WT cs.bam (first {prim:,} monomers): "
          f"{100*sa/prim:.2f}% carry an SA tag (chimeric / under-digested).")
    print("   (BfaI is a normal palindromic cutter, so this is the baseline incomplete-digestion rate,")
    print("    NOT the ~25% structural Type IIS loss — that extra loss is what AlwI would additionally suffer.)")


def main():
    fa = pysam.FastaFile(GENOME)
    seq = fa.fetch(CHR1, 0, REGION_BP).upper()
    frags, p = part_A(seq)
    recov = part_B(seq, frags)
    miss = p * p
    print(f"NET: motif-digest misses ~{100*miss:.0f}% of AlwI-AlwI junctions; alignment recovers "
          f"~{100*recov:.0f}% of those -> residual loss ~{100*miss*(1-recov):.1f}%.")
    # monomer/contact gain for a mean 3.5-monomer read (2.5 junctions)
    J = 2.5
    print(f"     for a 3.5-monomer read (~{J} junctions): recover ~{miss*recov*J:.2f} monomers/read "
          f"(~{100*miss*recov*J/3.5:.0f}% more), and disproportionately more multi-way contacts.\n")
    part_C(BFAI_BAM)


if __name__ == '__main__':
    main()
