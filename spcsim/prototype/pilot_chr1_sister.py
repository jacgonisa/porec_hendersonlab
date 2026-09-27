"""
Pilot: one Chr1 trans-(inter-)sister contact — how well is it captured, and how many
reads to detect it, under synchronised vs unsynchronised populations?

Built directly on quick_mc.py (same cell/strand/BrdU machinery). We PLANT a trans-sister
contact between two Chr1 loci A and B (truth = opposite sister chromatids, i.e. q_trans=1),
then measure:

  eps  = P(the A-B contact is captured AND correctly called trans-sister | a read spans A-B)
         -- this is the sister-resolution efficiency; sync >> async lives here.
  f    = P(a read physically spans A-B in cis)  (geometry only, from the contact model)
  N    = reads needed = target / (f * eps)

Populations: ideal fully-labelled G2 (upper bound), HU-sync 40%, and async 3 h pulse.
Two loci: an arm locus and a pericentromere locus (co-labelling bias differs by timing).

Run:  python pilot_chr1_sister.py [n_cells] [n_geo]
Outputs: pilot_chr1_sister.csv, pilot_chr1_sister.png, summary to stdout.
"""
import math, random, sys
import numpy as np
import pandas as pd
import quick_mc as q

CHR = 0                                   # Chr1 (0-indexed)
NB1 = q.NB[CHR]                           # number of 10-kb bins on Chr1
SEPS = [1, 5, 10, 50, 200]               # bin separations
SEP_LABEL = {1: '10 kb', 5: '50 kb', 10: '100 kb', 50: '500 kb', 200: '2 Mb'}
REGIONS = {'arm': 500, 'pericentromere': 1400}   # anchor bins: 5 Mb (arm), 14 Mb (peri, CEN~15.1 Mb)
POPS = {
    'ideal G2 (sync, fully labelled)': dict(ideal=True),
    'HU-sync 40%': dict(sync_frac=0.4),
    'async 3 h pulse': dict(),
}
OPP = {'W': 'C', 'C': 'W'}


def detect_probs(p, n=2_000_000, seed=0):
    """P(BrdU called | monomer read on labelled strand) and the false-positive equivalent.
    Depends only on the monomer length / T / substitution model, not on locus, so precompute once."""
    rng = np.random.default_rng(seed)
    L = rng.lognormal(math.log(750), 0.45, n).astype(int) + 50
    nT = rng.binomial(L, 0.32)
    nB = rng.binomial(nT, p['s'])
    calls_lab = rng.binomial(nB, p['sens']) + rng.binomial(nT - nB, p['fp'])
    calls_unlab = rng.binomial(nT, p['fp'])
    thr = p['thr']
    return (calls_lab >= thr).mean(), (calls_unlab >= thr).mean()


def eps_for(pop_over, region_bin, n_cells, p_lab, p_unlab, seed=1):
    """Monte-Carlo eps (capture-and-correct-call) for the planted trans contact, per separation."""
    random.seed(seed)
    p = dict(q.BASE); p.update(pop_over)
    a0 = region_bin
    n_sep = len(SEPS)
    captured = np.zeros(n_sep); correct = np.zeros(n_sep); ll_correct = np.zeros(n_sep)
    false_trans = np.zeros(n_sep); g2_both = np.zeros(n_sep); g1d = np.zeros(n_sep)
    for _ in range(n_cells):
        cell = q.sample_cell(p)
        repA, cpsA, g1dA = cell.state(CHR, 0, a0)
        for si, sep in enumerate(SEPS):
            b0 = min(NB1 - 1, a0 + sep)
            repB, cpsB, g1dB = cell.state(CHR, 0, b0)
            if repA and repB:                       # both sisters present -> plant trans (q_trans=1)
                cidA, cidB, truth_trans = 'pW', 'pC', True
                g2_both[si] += 1
            else:                                   # contact is not a resolvable sister pair
                cidA = random.choice(list(cpsA)); cidB = random.choice(list(cpsB)); truth_trans = False
            labA, labB = cpsA[cidA], cpsB[cidB]
            rsA, rsB = random.choice('WC'), random.choice('WC')
            calledA = random.random() < (p_lab if rsA in labA else p_unlab)
            calledB = random.random() < (p_lab if rsB in labB else p_unlab)
            identA = rsA if calledA else OPP[rsA]
            identB = rsB if calledB else OPP[rsB]
            both_L = calledA and calledB
            # rescue: within 100 kb (<=10 bins), an unlabelled monomer next to a labelled one is usable
            usable = both_L or (sep <= 10 and (calledA or calledB))
            call_trans = identA != identB
            if usable and call_trans:
                if truth_trans:
                    correct[si] += 1
                    if both_L:
                        ll_correct[si] += 1
                else:
                    false_trans[si] += 1
            if usable:
                captured[si] += 1
            if (g1dA or g1dB):
                g1d[si] += 1
    return dict(eps=correct / n_cells, eps_ll=ll_correct / n_cells,
                purity=correct / np.maximum(correct + false_trans, 1),
                p_g2both=g2_both / n_cells, g1d_frac=g1d / n_cells)


def geometry_f(n_geo, seed=7):
    """P(a read physically spans a specific Chr1 same-homolog cis pair at each separation)."""
    random.seed(seed)
    p = q.BASE
    counts = np.zeros(NB1, dtype=np.int64)   # counts[sep] = # same-homolog cis Chr1 pairs at that sep
    mean_mono, p_inter, p_hom = p['mean_mono'], p['p_inter'], p['p_hom']
    nprng = np.random.default_rng(seed)
    for _ in range(n_geo):
        n = 1 + nprng.poisson(mean_mono - 1)
        c = random.choices(range(5), q.CHR_W)[0]; h = random.randint(0, 1); b = random.randrange(q.NB[c])
        mons = [(c, h, b)]
        for _k in range(n - 1):
            u = random.random()
            if u < p_inter:
                c2 = random.choice([x for x in range(5) if x != c]); h2 = random.randint(0, 1); b2 = random.randrange(q.NB[c2])
            elif u < p_inter + p_hom:
                c2, h2 = c, 1 - h; b2 = q.clipb(c, b + random.choice((-1, 1)) * q.pick_d(c))
            else:
                c2, h2 = c, h; b2 = q.clipb(c, b + random.choice((-1, 1)) * q.pick_d(c))
            mons.append((c2, h2, b2))
        # same-homolog cis Chr1 pairs
        chr1 = [(hh, bb) for (cc, hh, bb) in mons if cc == CHR]
        for i in range(len(chr1)):
            for j in range(i + 1, len(chr1)):
                if chr1[i][0] == chr1[j][0]:
                    d = abs(chr1[i][1] - chr1[j][1])
                    if 0 < d < NB1:
                        counts[d] += 1
    f_band = counts / n_geo                       # per-read prob of ANY Chr1 same-hom pair at sep d
    f_pixel = f_band / np.maximum(NB1 - np.arange(NB1), 1)   # per specific pixel pair
    return f_band, f_pixel


def main():
    n_cells = int(sys.argv[1]) if len(sys.argv) > 1 else 150_000
    n_geo = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
    p_lab, p_unlab = detect_probs(q.BASE)
    print(f"per-monomer detection: P(call|labelled)={p_lab:.3f}  P(call|unlabelled,FP)={p_unlab:.4f}", file=sys.stderr)
    print("estimating contact geometry f(s)...", file=sys.stderr)
    f_band, f_pixel = geometry_f(n_geo)

    rows = []
    for region, a0 in REGIONS.items():
        for pop, over in POPS.items():
            print(f"eps: {region:14s} {pop}", file=sys.stderr, flush=True)
            r = eps_for(over, a0, n_cells, p_lab, p_unlab)
            for si, sep in enumerate(SEPS):
                eps = r['eps'][si]
                fb, fp = f_band[sep], f_pixel[sep]
                yield_pixel = fp * eps
                yield_band = fb * eps
                rows.append(dict(
                    region=region, population=pop, sep=SEP_LABEL[sep], sep_bins=sep,
                    P_contact_is_sister=round(r['p_g2both'][si], 4),
                    eps_capture=round(eps, 5), eps_LL_only=round(r['eps_ll'][si], 5),
                    purity=round(r['purity'][si], 3), g1daughter_frac=round(r['g1d_frac'][si], 4),
                    f_pixel=fp, f_band_chr1=fb,
                    reads_pixel_1=(1 / yield_pixel if yield_pixel > 0 else float('inf')),
                    reads_pixel_10=(10 / yield_pixel if yield_pixel > 0 else float('inf')),
                    reads_band_10=(10 / yield_band if yield_band > 0 else float('inf')),
                ))
    df = pd.DataFrame(rows)
    df.to_csv('pilot_chr1_sister.csv', index=False)

    # headline: async/sync read multiplier from eps ratio (region=pericentromere)
    piv = df[df.region == 'pericentromere'].pivot(index='sep', columns='population', values='eps_capture')
    print("\n=== eps (capture & correct trans-call), pericentromere Chr1 ===")
    print(piv.reindex([SEP_LABEL[s] for s in SEPS]).round(4).to_string())
    ideal_col = 'ideal G2 (sync, fully labelled)'; async_col = 'async 3 h pulse'
    mult = (piv[ideal_col] / piv[async_col].replace(0, np.nan)).reindex([SEP_LABEL[s] for s in SEPS])
    print("\nasync/sync read multiplier (eps_ideal / eps_async):")
    print(mult.round(1).to_string())

    print("\n=== reads to see the SPECIFIC Chr1 pixel >=10x (pericentromere) ===")
    pv = df[df.region == 'pericentromere'].pivot(index='sep', columns='population', values='reads_pixel_10')
    print(pv.reindex([SEP_LABEL[s] for s in SEPS]).map(lambda v: f"{v:.2e}").to_string())
    print("\n=== reads to accumulate 10 informative trans reads at that separation ANYWHERE on Chr1 ===")
    pb = df[df.region == 'pericentromere'].pivot(index='sep', columns='population', values='reads_band_10')
    print(pb.reindex([SEP_LABEL[s] for s in SEPS]).map(lambda v: f"{v:.2e}").to_string())

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for pop in POPS:
            sub = df[(df.region == 'pericentromere') & (df.population == pop)].set_index('sep_bins').reindex(SEPS)
            ax[0].plot(SEPS, sub['eps_capture'], 'o-', label=pop)
            ax[1].plot(SEPS, sub['reads_band_10'], 'o-', label=pop)
        for a in ax:
            a.set_xscale('log'); a.set_xticks(SEPS); a.set_xticklabels([SEP_LABEL[s] for s in SEPS]); a.set_xlabel('A-B separation')
        ax[0].set_yscale('log'); ax[0].set_ylabel('eps (capture & correct trans call)'); ax[0].set_title('Sister-resolution efficiency (Chr1 pericentromere)')
        ax[1].set_yscale('log'); ax[1].set_ylabel('reads for 10 informative trans pairs'); ax[1].set_title('Reads to detect (per separation, genome-wide Chr1)')
        ax[0].legend(fontsize=8); fig.tight_layout(); fig.savefig('pilot_chr1_sister.png', dpi=130)
        print("\nwrote pilot_chr1_sister.png", file=sys.stderr)
    except Exception as e:
        print("figure skipped:", e, file=sys.stderr)


if __name__ == '__main__':
    main()
