"""
REAL Pore-C contact maps from Jacob's AlwI and BfaI WT data (.mcool via cooler).
Genome-wide (nuclear Chr1-5) at 100 kb: shows the actual distance decay, pericentromere
blocks, chromocenter clustering and KNOT. Also computes a real chromocenter contact
fraction p_hub to make the KEE/chromocenter read estimates concrete.

Run: python real_contactmap.py
Outputs: real_contactmap_AlwI_BfaI.png, prints p_hub table.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cooler

BASE = '/home/jg2070/Desktop/PhD/PoreC/experiment_results'
MCOOLS = {
    'AlwI (WT, May 2025)':  f'{BASE}/2025May/results_all_formats_all_runs_AlwI/pairs/all.mcool',
    'BfaI (WT, Dec 2025)':  f'{BASE}/2025December/results_BfaI_WT/WT.mcool',
}
RES = 100_000
CENBED = '/home/jg2070/Desktop/PhD/PoreC/chr_cens.bed'
FLANK = 2_000_000   # pericentromere flank around the centromere interval


def load_cens():
    cen = {}
    for line in open(CENBED):
        f = line.split()
        # find the two integer columns = start,end; the field before is the chrom name
        for i in range(1, len(f) - 1):
            if f[i].isdigit() and f[i + 1].isdigit():
                cen[f[i - 1]] = (int(f[i]), int(f[i + 1])); break
    return cen


def get_matrix(mcool, res):
    c = cooler.Cooler(f'{mcool}::/resolutions/{res}')
    nuc = [ch for ch in c.chromnames if ch.lower().startswith('chr')]
    bins = c.bins()[:]
    mask = bins['chrom'].isin(nuc).values
    try:
        M = c.matrix(balance=True)[:][np.ix_(mask, mask)]
        if np.isnan(M).all(): raise ValueError
        M = np.nan_to_num(M)
    except Exception:
        M = c.matrix(balance=False)[:][np.ix_(mask, mask)].astype(float)
    bn = bins[mask].reset_index(drop=True)
    return M, bn, nuc, c


def main():
    cen = load_cens()
    fig, axes = plt.subplots(1, len(MCOOLS), figsize=(13.5, 6.6))
    print(f"{'enzyme':22s} {'contacts':>12s} {'inter%':>7s} {'p_hub (chromocenter / intra-peri)':>34s}")
    for ax, (name, mc) in zip(axes, MCOOLS.items()):
        M, bn, nuc, c = get_matrix(mc, RES)
        NB = M.shape[0]
        # chromosome boundaries + centromere bin positions (in the nuclear-only index)
        bounds = [0]; cen_bins = []; labels = []
        for ch in nuc:
            sub = bn[bn['chrom'] == ch]
            bounds.append(sub.index[-1] + 1); labels.append(ch.split(':')[0])
            if ch in cen:
                s, e = cen[ch]; mid = (s + e) // 2
                cen_bins.append(sub.index[0] + mid // RES)
        # display
        vmax = np.percentile(M[M > 0], 99.5) if (M > 0).any() else 1
        ax.imshow(np.log1p(M), cmap='Reds', vmax=np.log1p(vmax), interpolation='nearest', origin='upper')
        for b in bounds[1:-1]:
            ax.axhline(b - 0.5, color='#4169b0', lw=0.5, alpha=0.6); ax.axvline(b - 0.5, color='#4169b0', lw=0.5, alpha=0.6)
        for cb in cen_bins:
            ax.scatter([cb], [-4], marker='v', s=22, color='#2dd4bf', clip_on=False)
            ax.scatter([-4], [cb], marker='>', s=22, color='#2dd4bf', clip_on=False)
        centers = [(bounds[i] + bounds[i + 1]) / 2 for i in range(len(bounds) - 1)]
        ax.set_xticks(centers); ax.set_xticklabels(labels, fontsize=8)
        ax.set_yticks(centers); ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(name, fontsize=11)

        # real chromocenter p_hub: inter-chromosomal pericentromere-pericentromere contact fraction
        peri = np.zeros(NB, bool)
        chrom_of = bn['chrom'].astype(str).values; start_of = bn['start'].values
        for ch in nuc:
            if ch in cen:
                s, e = cen[ch]
                peri |= (chrom_of == ch) & (start_of >= s - FLANK) & (start_of <= e + FLANK)
        total = M.sum()
        interchrom = chrom_of[:, None] != chrom_of[None, :]
        pp = peri[:, None] & peri[None, :]
        p_cc = M[interchrom & pp].sum() / total
        p_intra_peri = M[(~interchrom) & pp].sum() / total       # intra-chrom pericentromere = the sister-relevant hub
        p_inter = M[interchrom].sum() / total
        print(f"{name:22s} {total:>12.3e} {100*p_inter:>11.1f}% inter-CC p_hub={p_cc:.2e}  intra-peri p_hub={p_intra_peri:.2e}")
    fig.suptitle('Real Pore-C contact maps (WT, 100 kb, nuclear) — teal ▾ = centromere', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig('real_contactmap_AlwI_BfaI.png', dpi=140)
    print('\nwrote real_contactmap_AlwI_BfaI.png')


if __name__ == '__main__':
    main()
