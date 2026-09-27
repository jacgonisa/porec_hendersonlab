"""
Can sister-Pore-C detect the "conformational asymmetry of replicated chromosomes"?

Gerlich-lab result (Science 2025, biorxiv 2025.07.09.663929): sister chromatids are consistently
STAGGERED by a register shift delta along the 5'->3' direction of the inherited (parental) strand.
Trans-sister contacts between locus i (sister pW) and j (sister pC) peak at j - i = +delta, not 0.
The shift is registered to parental-strand orientation.

Why sister-Pore-C can see it: the BrdU strand read tells you which sister a labelled monomer is on
(nascent C -> pW, nascent W -> pC), so each trans-sister contact can be SIGNED by the parental-strand
axis. Bulk Hi-C/Pore-C cannot orient the pair, so +delta and -delta staggers superimpose into a
symmetric diagonal and the asymmetry is invisible.

This plants a staggered trans-sister map on a Chr1 arm window and asks, as a function of BrdU PULSE
LENGTH (1/3/6 h; ideal G2 = ceiling):
  1. do strand-oriented trans contacts recover the +delta stagger, while the bulk map stays symmetric?
  2. how many oriented trans-sister contacts -> how many genome reads to detect the shift at 80% power?

Reuses quick_mc.py for cells + labelling. Run: python contactmap_sister.py [n_reads]
Outputs: contactmap_sister.png, contactmap_sister_readsneeded.csv
"""
import math, random, sys
import numpy as np
import quick_mc as q

CHR = 0
R = 300                      # window bins (10 kb each) = 3 Mb
OFF = 500                    # window starts at Chr1 bin 500 (arm, ~5 Mb)
DELTA = 3                    # register shift (bins) = 30 kb  (parameter; real value from paper)
SIGMA = 8                    # trans-contact spread (bins) = 80 kb
# BrdU pulse sweep (HU dropped). ideal = theoretical ceiling (fully-labelled G2).
POPS = {'ideal G2 (ceiling)': dict(ideal=True), '1 h pulse': dict(P=1.0),
        '3 h pulse': dict(P=3.0), '6 h pulse': dict(P=6.0)}
MAP_PANELS = ['3 h pulse', '6 h pulse']
OPP = {'W': 'C', 'C': 'W'}


def detect_probs(p, n=1_000_000, seed=0):
    rng = np.random.default_rng(seed)
    nT = rng.binomial(rng.lognormal(math.log(750), 0.45, n).astype(int) + 50, 0.32)
    nB = rng.binomial(nT, p['s'])
    lab = (rng.binomial(nB, p['sens']) + rng.binomial(nT - nB, p['fp'])) >= p['thr']
    unlab = rng.binomial(nT, p['fp']) >= p['thr']
    return lab.mean(), unlab.mean()


def truth_trans_map(n=400_000, seed=3):
    rng = np.random.default_rng(seed)
    i = rng.integers(0, R, n)
    j = np.round(i + DELTA + rng.normal(0, SIGMA, n)).astype(int)
    ok = (j >= 0) & (j < R)
    return i[ok], j[ok]


def geometry_window(n_geo=400_000, seed=7):
    """P(a random read makes a cis same-homolog contact with BOTH loci inside the window)."""
    random.seed(seed); nprng = np.random.default_rng(seed)
    p = q.BASE; hit = 0
    for _ in range(n_geo):
        n = 1 + nprng.poisson(p['mean_mono'] - 1)
        c = random.choices(range(5), q.CHR_W)[0]; h = random.randint(0, 1); b = random.randrange(q.NB[c])
        mons = [(c, h, b)]
        for _k in range(n - 1):
            u = random.random()
            if u < p['p_inter']:
                c2 = random.choice([x for x in range(5) if x != c]); h2 = random.randint(0, 1); b2 = random.randrange(q.NB[c2])
            elif u < p['p_inter'] + p['p_hom']:
                c2, h2 = c, 1 - h; b2 = q.clipb(c, b + random.choice((-1, 1)) * q.pick_d(c))
            else:
                c2, h2 = c, h; b2 = q.clipb(c, b + random.choice((-1, 1)) * q.pick_d(c))
            mons.append((c2, h2, b2))
        win = [(hh, bb) for (cc, hh, bb) in mons if cc == CHR and OFF <= bb < OFF + R]
        found = False
        for i in range(len(win)):
            for j in range(i + 1, len(win)):
                if win[i][0] == win[j][0]:
                    found = True
        hit += found
    return hit / n_geo


def simulate(pop_over, n_reads, p_lab, p_unlab, seed=1):
    """Each read attempts one window trans-sister contact; classify by strand+label."""
    random.seed(seed); rng = np.random.default_rng(seed)
    p = dict(q.BASE); p.update(pop_over)
    ori_i, ori_j, bulk_signed = [], [], []
    for _ in range(n_reads):
        cell = q.sample_cell(p)
        i = random.randrange(R)
        j = int(round(i + DELTA + rng.normal(0, SIGMA)))
        if not (0 <= j < R):
            continue
        repA, cpsA, _ = cell.state(CHR, 0, OFF + i)
        repB, cpsB, _ = cell.state(CHR, 0, OFF + j)
        if not (repA and repB):
            continue                                  # no two sisters -> no trans-sister contact
        labA, labB = cpsA['pW'], cpsB['pC']
        rsA, rsB = random.choice('WC'), random.choice('WC')
        cA = random.random() < (p_lab if rsA in labA else p_unlab)
        cB = random.random() < (p_lab if rsB in labB else p_unlab)
        if cA and cB:
            ori_i.append(i); ori_j.append(j)
            bulk_signed.append(j - i if random.random() < 0.5 else i - j)
    n = n_reads
    return np.array(ori_i), np.array(ori_j), np.array(bulk_signed), len(ori_i) / n


def hist2d(i, j):
    H, _, _ = np.histogram2d(i, j, bins=[np.arange(R + 1), np.arange(R + 1)])
    return H


CONTACTS_80 = (1.96 + 0.84) ** 2 * (SIGMA / DELTA) ** 2   # oriented trans contacts for 80% power


def main():
    n_reads = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    p_lab, p_unlab = detect_probs(q.BASE)
    f_window = geometry_window()
    print(f"P(call|labelled)={p_lab:.3f}  delta={DELTA*10} kb sigma={SIGMA*10} kb  "
          f"f_window(read makes cis contact in 3Mb window)={f_window:.2e}", file=sys.stderr)
    print(f"oriented trans contacts needed for 80% power = {CONTACTS_80:.0f}", file=sys.stderr)

    gi, gj = truth_trans_map()
    results = {}
    for pop, over in POPS.items():
        oi, oj, bulk, y_attempt = simulate(over, n_reads, p_lab, p_unlab)
        # y_attempt = oriented trans per window-contact-ATTEMPT; genome yield = f_window * y_attempt
        y_genome = f_window * y_attempt
        reads_needed = CONTACTS_80 / y_genome if y_genome > 0 else float('inf')
        results[pop] = dict(oi=oi, oj=oj, bulk=bulk, y_attempt=y_attempt, y_genome=y_genome,
                            mean_signed=(oj - oi).mean() if len(oi) else float('nan'),
                            n_oriented=len(oi), reads_needed=reads_needed)
        print(f"{pop:20s} oriented/attempt={y_attempt:.3f} genome/read={y_genome:.2e} "
              f"mean offset={results[pop]['mean_signed']*10:5.1f} kb  reads80%={reads_needed:.2e}", file=sys.stderr)

    import csv
    with open('contactmap_sister_readsneeded.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['population', 'oriented_trans_per_read_genome', 'contacts_for_80pct_power',
                    'reads_for_80pct_power', 'flowcells_at_7M'])
        for pop, r in results.items():
            w.writerow([pop, f"{r['y_genome']:.3e}", f"{CONTACTS_80:.0f}",
                        f"{r['reads_needed']:.3e}", f"{r['reads_needed']/7e6:.3f}"])
    print("\n=== reads to DETECT the register shift (mean signed offset != 0, 80% power) ===")
    print(f"delta={DELTA*10} kb, sigma={SIGMA*10} kb  ->  {CONTACTS_80:.0f} oriented trans contacts needed")
    print(f"{'BrdU pulse':20s} {'oriented trans/read':>19s} {'reads needed':>14s} {'flow cells (7M)':>15s}")
    for pop, r in results.items():
        print(f"{pop:20s} {r['y_genome']:>19.2e} {r['reads_needed']:>14.2e} {r['reads_needed']/7e6:>15.3f}")

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        from scipy.ndimage import gaussian_filter
        fig, ax = plt.subplots(2, 3, figsize=(14, 9))
        def show(a, H, title):
            a.imshow(np.log1p(gaussian_filter(H, 2)).T, origin='lower', cmap='magma', aspect='equal')
            a.plot([0, R], [0, R], 'c--', lw=0.6, alpha=0.6)
            a.plot([0, R - DELTA], [DELTA, R], 'w:', lw=0.7, alpha=0.7)
            a.set_title(title, fontsize=10); a.set_xlabel('locus on pW (10 kb bins)'); a.set_ylabel('locus on pC')
        show(ax[0, 0], hist2d(gi, gj), f'ground-truth trans-sister\n(staggered +{DELTA*10} kb)')
        for k, pop in enumerate(MAP_PANELS):
            r = results[pop]
            show(ax[0, 1 + k], hist2d(r['oi'], r['oj']),
                 f"recovered oriented trans — {pop}\n(n={r['n_oriented']} @ {n_reads:,} attempts)")

        bins = np.arange(-40, 41)
        a = ax[1, 0]
        for pop in ['1 h pulse', '3 h pulse', '6 h pulse']:
            off = results[pop]['oj'] - results[pop]['oi']
            if len(off):
                a.hist(off, bins=bins, histtype='step', density=True, label=f'oriented: {pop}')
        a.hist(results['3 h pulse']['bulk'], bins=bins, histtype='step', density=True, color='gray', lw=1.5, label='bulk (no strand)')
        a.axvline(0, color='c', ls='--', lw=0.8); a.axvline(DELTA, color='k', ls=':', lw=0.8)
        a.set_xlabel('signed offset j-i (10 kb bins)'); a.set_ylabel('density')
        a.set_title('Strand-oriented recovers +delta stagger;\nbulk stays symmetric at 0', fontsize=10); a.legend(fontsize=7)

        a = ax[1, 1]
        deltas = np.array([1, 2, 3, 5, 8, 12])
        for pop in POPS:
            yg = results[pop]['y_genome']
            nr = [(1.96 + 0.84) ** 2 * (SIGMA / d) ** 2 / yg for d in deltas]
            a.plot(deltas * 10, nr, 'o-', label=pop)
        a.set_yscale('log'); a.set_xlabel('true register shift delta (kb)'); a.set_ylabel('reads for 80% power')
        a.axhline(7e6, color='k', ls=':', lw=0.8); a.text(deltas[-1]*10, 7e6, ' 1 flow cell', fontsize=7, va='bottom', ha='right')
        a.set_title('Reads to detect the shift vs its size & pulse', fontsize=10); a.legend(fontsize=7)

        a = ax[1, 2]
        names = list(POPS); means = [results[p]['mean_signed'] * 10 for p in names]
        errs = [(SIGMA * 10) / math.sqrt(max(results[p]['n_oriented'], 1)) for p in names]
        a.bar(range(len(names)), means, yerr=errs, capsize=4,
              color=['#2f855a', '#f6ad55', '#2b6cb0', '#c53030'])
        a.axhline(DELTA * 10, color='k', ls=':', lw=1); a.axhline(0, color='c', ls='--', lw=0.8)
        a.set_xticks(range(len(names))); a.set_xticklabels([n.replace(' pulse', '').replace(' (ceiling)', '') for n in names], fontsize=8)
        a.set_ylabel('recovered mean signed offset (kb)')
        a.set_title(f'Recovered stagger @ {n_reads:,} attempts\n(dotted = truth {DELTA*10} kb)', fontsize=10)

        fig.suptitle('Sister-Pore-C detection of replicated-chromosome conformational asymmetry — by BrdU pulse length', fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig('contactmap_sister.png', dpi=130)
        print('\nwrote contactmap_sister.png', file=sys.stderr)
    except Exception:
        import traceback; traceback.print_exc()


if __name__ == '__main__':
    main()
