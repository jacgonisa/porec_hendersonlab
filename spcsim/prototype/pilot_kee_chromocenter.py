"""
Pilot 3: are KEE / chromocenter interactions intra-sister, inter-sister, or both?
And can we tell "sister-cohesion-driven" from "homolog/bulk-driven" clustering?

For an enriched hub interaction between two loci A,B we plant a TRUE mixture of contact types:
  - same-homolog sister contacts, split by 'sister_trans_frac' (trans vs cis-sister)
  - inter-homolog contacts (the confound: in inbred Col-0 they masquerade 50:50 as cis/trans)
We then read strand + BrdU (via quick_mc), classify cis/trans, and in F1 mode try to HAPLOTAG each
monomer (region-dependent success) to identify & remove inter-homolog pairs.

Readouts:
  1. recovered trans-sister fraction vs truth, inbred vs F1  -> the inter-homolog bias & F1 rescue
  2. sister-resolution efficiency eps at the hub (async 3 h pulse)
  3. contacts needed for +/-0.05 CI, and to distinguish cohesion(0.5) vs bulk(0.15) at 80% power
  4. reads = contacts / (p_hub * eps)   -> plug your MEASURED hub contact fraction p_hub
  5. multi-way (3-locus) sister-spanning hub: extra label cost

Regions: KEE (arm, low inter-homolog, good haplotag) vs pericentromere/chromocenter
(high inter-homolog, poor haplotag). Run: python pilot_kee_chromocenter.py [n]
"""
import math, random, sys
import numpy as np
import quick_mc as q

CHR = 0
OPP = {'W': 'C', 'C': 'W'}
REGIONS = {
    # A,B,C bins on Chr1; inter_homolog_frac; F1 haplotag success (SNP density/mappability)
    'KEE (arm)':            dict(A=520, B=1040, C=1560, p_hom=0.05, hsucc=0.70),
    'pericentromere':       dict(A=1450, B=1560, C=1505, p_hom=0.40, hsucc=0.25),
}
SCEN = {'sister-cohesion': 0.50, 'bulk/independent': 0.15}   # true trans-sister fraction (same-homolog)
POP = dict(q.BASE)                                            # async 3 h pulse (HU dropped)


def detect_probs(p, n=1_000_000, seed=0):
    rng = np.random.default_rng(seed)
    nT = rng.binomial(rng.lognormal(math.log(750), 0.45, n).astype(int) + 50, 0.32)
    lab = (rng.binomial(rng.binomial(nT, p['s']), p['sens']) + rng.binomial(nT, p['fp'])) >= p['thr']
    unlab = rng.binomial(nT, p['fp']) >= p['thr']
    return lab.mean(), unlab.mean()


def call_from_lab(lab, p_lab, p_unlab):
    """strand read + BrdU call given a chromatid's labelled-strand set. returns (ident, called)."""
    rs = random.choice('WC')
    called = random.random() < (p_lab if rs in lab else p_unlab)
    return (rs if called else OPP[rs]), called


def simulate(region, genotype, stf, n_reads, p_lab, p_unlab, seed=1):
    random.seed(seed)
    A, B, C = region['A'], region['B'], region['C']
    p_hom, hsucc = region['p_hom'], (region['hsucc'] if genotype == 'F1' else 0.0)
    kept = trans_calls = 0
    n3_resolved = 0                                   # fully-labelled 3-way (sister-spanning)
    for _ in range(n_reads):
        cell = q.sample_cell(POP)
        # choose homologs + true type FIRST, then read each locus' state exactly once (state() is stochastic)
        if random.random() < p_hom:                   # inter-homolog contact (the confound)
            hA, hB, interhom = 0, 1, True
        else:                                         # same-homolog sister contact
            hA = hB = random.randint(0, 1); interhom = False
            is_trans = random.random() < stf
        sA = cell.state(CHR, hA, A); sB = cell.state(CHR, hB, B)
        if not (sA[0] and sB[0]):                     # need two sisters present at both loci
            continue
        if interhom:
            cidA = random.choice(list(sA[1])); cidB = random.choice(list(sB[1]))
        elif is_trans:
            cidA, cidB = 'pW', 'pC'
        else:
            cidA = cidB = random.choice(('pW', 'pC'))
        identA, cA = call_from_lab(sA[1][cidA], p_lab, p_unlab)
        identB, cB = call_from_lab(sB[1][cidB], p_lab, p_unlab)
        if not (cA and cB):
            continue                                  # not both BrdU+ -> uninformative
        if genotype == 'F1':                          # haplotag: identified inter-homolog is removed
            if (random.random() < hsucc) and (random.random() < hsucc) and (hA != hB):
                continue
        kept += 1
        trans_calls += (identA != identB)
        if not interhom:                              # opportunistic sister-spanning 3-way at C
            sC = cell.state(CHR, hA, C)
            if sC[0]:
                cidC = 'pC' if cidA == 'pW' else 'pW'
                if call_from_lab(sC[1][cidC], p_lab, p_unlab)[1]:
                    n3_resolved += 1
    return dict(recovered=trans_calls / kept if kept else float('nan'),
                eps=kept / n_reads, eps3=n3_resolved / n_reads)


def contacts_precision(f, halfwidth=0.05, z=1.96):
    return (z / halfwidth) ** 2 * f * (1 - f)
def contacts_discriminate(f1, f2, za=1.96, zb=0.84):
    return (za + zb) ** 2 * (f1 * (1 - f1) + f2 * (1 - f2)) / (f1 - f2) ** 2


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150_000
    p_lab, p_unlab = detect_probs(POP)
    print(f"P(call|labelled)={p_lab:.3f}  (async 3 h pulse)\n", file=sys.stderr)

    res = {}   # (region, scen, geno) -> dict
    for region in REGIONS:
        for scen, stf in SCEN.items():
            for geno in ('inbred', 'F1'):
                res[(region, scen, geno)] = simulate(REGIONS[region], geno, stf, n, p_lab, p_unlab)

    print("=== 1. recovered trans-sister fraction vs TRUTH ===")
    print(f"{'region':16s} {'scenario':18s} {'true':>5s} {'inbred':>8s} {'F1':>7s}")
    for region in REGIONS:
        for scen, stf in SCEN.items():
            ib = res[(region, scen, 'inbred')]['recovered']; f1 = res[(region, scen, 'F1')]['recovered']
            print(f"{region:16s} {scen:18s} {stf:>5.2f} {ib:>8.2f} {f1:>7.2f}")
    print("\n  -> inter-homolog inflates the inbred estimate toward 0.5; worst at the pericentromere,")
    print("     where 'bulk' clustering can FAKE sister cohesion. F1 haplotagging recovers the truth.\n")

    print("=== 2. sister-resolution efficiency eps (informative hub contacts per hub contact) ===")
    for region in REGIONS:
        e_ib = res[(region, 'sister-cohesion', 'inbred')]['eps']
        e_f1 = res[(region, 'sister-cohesion', 'F1')]['eps']
        print(f"  {region:16s} inbred eps={e_ib:.3f}   F1 eps={e_f1:.3f} (haplotag drops some pairs)")

    print("\n=== 3. contacts needed (robust; independent of contact frequency) ===")
    for region in REGIONS:
        for geno in ('inbred', 'F1'):
            fc = res[(region, 'sister-cohesion', geno)]['recovered']
            fb = res[(region, 'bulk/independent', geno)]['recovered']
            n_prec = contacts_precision(fc)
            n_disc = contacts_discriminate(fc, fb)
            print(f"  {region:16s} {geno:6s}: +/-0.05 CI -> {n_prec:5.0f} contacts | "
                  f"cohesion-vs-bulk -> {n_disc:6.0f} contacts")

    print("\n=== 4. reads = contacts / (p_hub * eps).  Plug your MEASURED hub contact fraction p_hub ===")
    print(f"{'region':16s} {'geno':6s} {'contacts(disc)':>14s} " + " ".join(f"p_hub={p:>7.0e}" for p in (1e-3, 1e-4, 1e-5)))
    for region in REGIONS:
        for geno in ('inbred', 'F1'):
            fc = res[(region, 'sister-cohesion', geno)]['recovered']; fb = res[(region, 'bulk/independent', geno)]['recovered']
            nd = contacts_discriminate(fc, fb); eps = res[(region, 'sister-cohesion', geno)]['eps']
            reads = [nd / (ph * eps) for ph in (1e-3, 1e-4, 1e-5)]
            print(f"  {region:14s} {geno:6s} {nd:>14.0f} " + " ".join(f"{r:>12.1e}" for r in reads))

    print("\n=== 5. multi-way (3-locus sister-spanning) extra cost ===")
    for region in REGIONS:
        e2 = res[(region, 'sister-cohesion', 'inbred')]['eps']
        e3 = res[(region, 'sister-cohesion', 'inbred')]['eps3']
        print(f"  {region:16s} 2-way eps={e2:.3f}  3-way eps={e3:.4f}  -> {e2/e3 if e3 else float('inf'):.1f}x more reads for a 3-way")

    # ---- figure ----
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        # panel 1: recovered fraction inbred vs F1
        regs = list(REGIONS); scens = list(SCEN); x = np.arange(len(regs)); w = 0.2
        colors = {'inbred': '#c53030', 'F1': '#2f855a'}
        for si, scen in enumerate(scens):
            for gi, geno in enumerate(('inbred', 'F1')):
                vals = [res[(r, scen, geno)]['recovered'] for r in regs]
                ax[0].bar(x + (si*2+gi-1.5)*w, vals, w, color=colors[geno],
                          alpha=0.6 if scen == 'bulk/independent' else 1.0,
                          label=f"{geno}, {scen}" )
        for si, scen in enumerate(scens):
            ax[0].axhline(SCEN[scen], color='#2b6cb0', ls=':', lw=1)
            ax[0].text(len(regs)-0.4, SCEN[scen], f'true {SCEN[scen]}', color='#2b6cb0', fontsize=8, va='bottom')
        ax[0].set_xticks(x); ax[0].set_xticklabels(regs); ax[0].set_ylabel('recovered trans-sister fraction')
        ax[0].set_title('Inter-homolog bias & F1 rescue'); ax[0].legend(fontsize=7); ax[0].set_ylim(0, 0.6)
        # panel 2: reads to discriminate vs p_hub
        ph = np.logspace(-5, -2.5, 30)
        for region in regs:
            for geno in ('inbred', 'F1'):
                fc = res[(region, 'sister-cohesion', geno)]['recovered']; fb = res[(region, 'bulk/independent', geno)]['recovered']
                nd = contacts_discriminate(fc, fb); eps = res[(region, 'sister-cohesion', geno)]['eps']
                ax[1].plot(ph, nd / (ph * eps), label=f"{region.split(' ')[0]}, {geno}",
                           ls='-' if geno == 'F1' else '--')
        ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].axhline(18e6, color='k', ls=':', lw=0.8)
        ax[1].text(ph[-1], 18e6, ' 1 flow cell (50 Gb)', fontsize=7, ha='right', va='bottom')
        ax[1].set_xlabel('measured hub contact fraction p_hub'); ax[1].set_ylabel('reads to distinguish cohesion vs bulk')
        ax[1].set_title('Reads vs how enriched your hub is'); ax[1].legend(fontsize=7)
        # panel 3: eps 2-way vs 3-way
        e2 = [res[(r, 'sister-cohesion', 'inbred')]['eps'] for r in regs]
        e3 = [res[(r, 'sister-cohesion', 'inbred')]['eps3'] for r in regs]
        ax[2].bar(x - 0.2, e2, 0.4, color='#2b6cb0', label='2-way (pairwise)')
        ax[2].bar(x + 0.2, e3, 0.4, color='#bc8cff', label='3-way (sister-spanning)')
        ax[2].set_xticks(x); ax[2].set_xticklabels(regs); ax[2].set_ylabel('resolution efficiency eps')
        ax[2].set_title('Multi-way costs one more label'); ax[2].legend(fontsize=8)
        fig.suptitle('Sister-Pore-C for KEE / chromocenter interactions (Chr1, async 3 h pulse)', fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig('pilot_kee_chromocenter.png', dpi=130)
        print('\nwrote pilot_kee_chromocenter.png', file=sys.stderr)
    except Exception:
        import traceback; traceback.print_exc()


if __name__ == '__main__':
    main()
