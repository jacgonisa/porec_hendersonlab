"""
REAL Pore-C contact maps from Jacob's AlwI and BfaI WT data (.mcool via cooler).
Two normalisations: ICE-balanced (top) and observed/expected log2 (bottom).
Run: python real_contactmap.py   ->  real_contactmap_AlwI_BfaI.png  (+ prints p_hub).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm
import cooler

BASE = '/home/jg2070/Desktop/PhD/PoreC/experiment_results'
MCOOLS = {
    'AlwI (WT, May 2025)':  f'{BASE}/2025May/results_all_formats_all_runs_AlwI/pairs/all.mcool',
    'BfaI (WT, Dec 2025)':  f'{BASE}/2025December/results_BfaI_WT/WT.mcool',
}
RES = 100_000
CENBED = '/home/jg2070/Desktop/PhD/PoreC/chr_cens.bed'
FLANK = 2_000_000
HIC = LinearSegmentedColormap.from_list('hic', ['#ffffff', '#fff0e0', '#fdcc8a', '#fc8d59', '#e34a33', '#b30000', '#6d0000'])
HIC.set_bad('#f2f2f2')
OEMAP = plt.get_cmap('RdBu_r').copy(); OEMAP.set_bad('#f2f2f2')


def load_cens():
    cen = {}
    for line in open(CENBED):
        f = line.split()
        for i in range(1, len(f) - 1):
            if f[i].isdigit() and f[i + 1].isdigit():
                cen[f[i - 1]] = (int(f[i]), int(f[i + 1])); break
    return cen


def ice(M, iters=60):
    M = M.astype(float).copy()
    cov = M.sum(1); good = cov > np.percentile(cov[cov > 0], 2)
    M[~good, :] = 0; M[:, ~good] = 0
    for _ in range(iters):
        marg = M.sum(1); mean = marg[good].mean()
        b = np.where(marg > 0, marg / mean, 1.0)
        M /= b[:, None]; M /= b[None, :]
    return M


def oe(M, chrom, nuc):
    """observed/expected (log2): divide each contact by the mean contact at its genomic distance."""
    OE = np.full(M.shape, np.nan)
    for ch in nuc:
        idx = np.where(chrom == ch)[0]; sub = M[np.ix_(idx, idx)]; n = len(idx)
        exp = np.array([(lambda d=np.diagonal(sub, k): d[d > 0].mean() if (d > 0).any() else np.nan)() for k in range(n)])
        D = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        with np.errstate(invalid='ignore', divide='ignore'):
            s = np.where(sub > 0, sub, np.nan)
            OE[np.ix_(idx, idx)] = s / exp[D]
    inter = chrom[:, None] != chrom[None, :]; iv = M[inter]
    ie = iv[iv > 0].mean() if (iv > 0).any() else np.nan
    if ie and not np.isnan(ie):
        OE[inter] = np.where(M[inter] > 0, M[inter] / ie, np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.log2(OE)


def annotate(ax, bn, nuc, cen, NB, marks=True):
    bounds = [0]; labels = []; cen_bins = []
    for ch in nuc:
        sub = bn[bn['chrom'] == ch]; bounds.append(sub.index[-1] + 1); labels.append(ch.split(':')[0])
        if ch in cen:
            s, e = cen[ch]; cen_bins.append(sub.index[0] + (s + e) // 2 // RES)
    for b in bounds[1:-1]:
        ax.axhline(b - .5, color='#9aa7b4', lw=.4); ax.axvline(b - .5, color='#9aa7b4', lw=.4)
    if marks:
        for cb in cen_bins:
            ax.scatter([cb], [NB * 1.012], marker='^', s=18, color='#2b6cb0', clip_on=False)
    cent = [(bounds[i] + bounds[i + 1]) / 2 for i in range(5)]
    ax.set_xticks(cent); ax.set_xticklabels(labels, fontsize=7)
    ax.set_yticks(cent); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(-0.5, NB - 0.5); ax.set_ylim(NB - 0.5, -0.5)


def main():
    cen = load_cens()
    fig, axes = plt.subplots(2, len(MCOOLS), figsize=(13, 12.6))
    print(f"{'enzyme':22s} {'contacts':>12s} {'inter%':>7s}  p_hub(intra-peri/inter-CC)")
    for k, (name, mc) in enumerate(MCOOLS.items()):
        c = cooler.Cooler(f'{mc}::/resolutions/{RES}')
        nuc = [ch for ch in c.chromnames if ch.lower().startswith('chr')]
        bins = c.bins()[:]; mask = bins['chrom'].isin(nuc).values
        Mraw = c.matrix(balance=False)[:][np.ix_(mask, mask)].astype(float)
        bn = bins[mask].reset_index(drop=True); chrom = bn['chrom'].astype(str).values; NB = Mraw.shape[0]
        Mb = ice(Mraw); OE = oe(Mb, chrom, nuc)
        # ICE panel
        a = axes[0, k]; pos = Mb[Mb > 0]
        im = a.imshow(np.where(Mb > 0, Mb, np.nan), cmap=HIC,
                      norm=LogNorm(vmin=np.percentile(pos, 45), vmax=np.percentile(pos, 99.7)),
                      interpolation='nearest', origin='upper')
        annotate(a, bn, nuc, cen, NB); a.set_title(f'{name} — ICE balanced', fontsize=11)
        fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
        # O/E panel
        a = axes[1, k]
        im = a.imshow(OE, cmap=OEMAP, norm=TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=2.5), interpolation='nearest', origin='upper')
        annotate(a, bn, nuc, cen, NB); a.set_title(f'{name} — observed/expected (log2)', fontsize=11)
        fig.colorbar(im, ax=a, fraction=0.046, pad=0.02)
        # p_hub
        start = bn['start'].values; peri = np.zeros(NB, bool)
        for ch in nuc:
            if ch in cen:
                s, e = cen[ch]; peri |= (chrom == ch) & (start >= s - FLANK) & (start <= e + FLANK)
        inter = chrom[:, None] != chrom[None, :]; pp = peri[:, None] & peri[None, :]; tot = Mraw.sum()
        print(f"{name:22s} {tot:>12.2e} {100*Mraw[inter].sum()/tot:>6.1f}%  "
              f"intra-peri={Mraw[(~inter)&pp].sum()/tot:.3f}  inter-CC={Mraw[inter&pp].sum()/tot:.2e}")
    fig.suptitle('Real Pore-C contact maps (WT, 100 kb) — ICE (top) and O/E log2 (bottom); ▲ centromere', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig('real_contactmap_AlwI_BfaI.png', dpi=150)
    print('\nwrote real_contactmap_AlwI_BfaI.png')


if __name__ == '__main__':
    main()
