"""
rescue_monomers.py — recover monomers lost to the digest-first step (Type IIS site-less junctions
for AlwI, and ordinary incomplete digestion for any enzyme).

A fused/under-digested monomer is a chimera: pore-c-py keeps it as ONE record, but minimap2 reports
the second fragment as a supplementary alignment (SA tag). Splitting that record at the SA breakpoints
recovers the extra monomer (and its contacts). Any large soft-clip WITHOUT an SA (minimap2 didn't
place the clipped piece) is re-aligned to recover it too.

Enzyme-agnostic; works on any pore-c-py monomer BAM (grouped by the MI concatemer tag).
Run (report):  python rescue_monomers.py --bam WT.cs.bam --limit 100000
"""
import argparse, re, subprocess, os
import pysam

MM2 = '/home/jg2070/minimap2/minimap2'
CIG = re.compile(r'(\d+)([MIDNSHP=X])')


def query_span(cigar, qlen=None):
    """read(query) start,end that a CIGAR consumes (soft/hard clips are the clipped ends)."""
    ops = CIG.findall(cigar)
    start = 0; i = 0
    while i < len(ops) and ops[i][1] in 'SH': start += int(ops[i][0]); i += 1
    qaln = sum(int(n) for n, o in ops if o in 'MI=X')
    return start, start + qaln


def sa_segments(sa):
    segs = []
    for rec in sa.strip(';').split(';'):
        if not rec: continue
        chrom, pos, strand, cig, mapq, nm = rec.split(',')
        s, e = query_span(cig)
        # if SA is reverse relative to primary, its query coords are on the reverse read; normalise length only
        segs.append(dict(chrom=chrom, pos=int(pos), strand=strand, mapq=int(mapq), qstart=s, qend=e, reflen=sum(int(n) for n,o in CIG.findall(cig) if o in 'MD=X')))
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bam', required=True)
    ap.add_argument('--limit', type=int, default=100000, help='primary monomers to scan')
    ap.add_argument('--realign-clips', action='store_true', help='also re-align large soft-clips lacking an SA')
    ap.add_argument('--ref', help='reference .mmi/.fa for --realign-clips')
    a = ap.parse_args()

    sf = pysam.AlignmentFile(a.bam, 'rb')
    n_prim = n_sa = extra_from_sa = 0
    sa_reflen_tot = 0
    n_clip_noSA = clip_bp_noSA = 0
    clip_seqs = []                                  # (name, seq) for optional re-alignment
    for r in sf:
        if r.is_secondary or r.is_supplementary: continue
        n_prim += 1
        if r.has_tag('SA'):
            n_sa += 1
            segs = sa_segments(r.get_tag('SA'))
            extra_from_sa += len(segs)              # each SA segment = a recovered monomer
            sa_reflen_tot += sum(s['reflen'] for s in segs)
        else:
            # large soft-clip with no SA -> minimap2 didn't place it; candidate for re-alignment
            clip = (r.query_length or 0) - r.query_alignment_length
            if clip >= 300 and r.query_sequence:
                n_clip_noSA += 1; clip_bp_noSA += clip
                if a.realign_clips and len(clip_seqs) < 5000:
                    qs, qe = r.query_alignment_start, r.query_alignment_end
                    for tag, seg in [('L', r.query_sequence[:qs]), ('R', r.query_sequence[qe:])]:
                        if len(seg) >= 300:
                            clip_seqs.append((f'{r.query_name}_{tag}', seg))
        if n_prim >= a.limit: break

    print(f"scanned {n_prim:,} primary monomers from {os.path.basename(a.bam)}")
    print(f"  SA-tagged (chimeric / under-digested): {n_sa:,}  ({100*n_sa/n_prim:.2f}%)")
    print(f"  -> monomers recovered by splitting at SA breakpoints: {extra_from_sa:,} "
          f"(+{100*extra_from_sa/n_prim:.2f}% monomers; mean recovered fragment {sa_reflen_tot/max(extra_from_sa,1):.0f} bp)")
    print(f"  large soft-clip (>=300 bp) with NO SA: {n_clip_noSA:,}  ({100*n_clip_noSA/n_prim:.2f}%), "
          f"{clip_bp_noSA/1e6:.1f} Mb clipped — recoverable by re-alignment")
    # contact gain ~ scales super-linearly: adding k monomers to an n-monomer read adds ~ (n+k choose 2)-(n choose 2)
    gain = (extra_from_sa) / n_prim
    print(f"  net monomer gain from SA split alone: +{100*gain:.1f}%  (more for multi-way contacts)")

    if a.realign_clips and clip_seqs:
        ref = a.ref
        with open('/tmp/clips.fa', 'w') as f:
            for name, seq in clip_seqs: f.write(f'>{name}\n{seq}\n')
        out = subprocess.run([MM2, '-ax', 'map-ont', '--secondary=no', ref, '/tmp/clips.fa'],
                             capture_output=True, text=True)
        mapped = sum(1 for l in out.stdout.splitlines() if l and l[0] != '@' and l.split('\t')[2] != '*')
        print(f"  re-aligned {len(clip_seqs)} no-SA clips -> {mapped} mapped "
              f"({100*mapped/len(clip_seqs):.0f}%) = additional recovered monomers")


if __name__ == '__main__':
    main()
