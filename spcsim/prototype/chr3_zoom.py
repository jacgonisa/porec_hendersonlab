"""
Chr3 zoom of the real AlwI & BfaI Pore-C maps, with data-driven KEE detection.
A KEE (KNOT-engaged element) makes anomalous long-range / inter-chromosomal contacts,
so it shows up as a peak in an inter-chromosomal contact track (outside the pericentromere).
We plot the Chr3 intra map for each enzyme + the KEE-signature track, and annotate the
strongest p-arm peak.

Run: python chr3_zoom.py   ->  chr3_zoom_AlwI_BfaI.png
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cooler

BASE = '/home/jg2070/Desktop/PhD/PoreC/experiment_results'
MCOOLS = {'AlwI': f'{BASE}/2025May/results_all_formats_all_runs_AlwI/pairs/all.mcool',
          'BfaI': f'{BASE}/2025December/results_BfaI_WT/WT.mcool'}
CENBED = '/home/jg2070/Desktop/PhD/PoreC/chr_cens.bed'
CHR3 = 'Chr3:1-26150667'
RES_MAP = 50_000        # Chr3 intra map resolution
RES_TRK = 100_000       # inter-chrom track resolution
FLANK = 2_000_000


def cen_of(name):
    for line in open(CENBED):
        f = line.split()
        for i in range(1, len(f) - 1):
            if f[i].isdigit() and f[i + 1].isdigit() and f[i - 1] == name:
                return int(f[i]), int(f[i + 1])
    return None


def chr3_intra(mcool, res):
    c = cooler.Cooler(f'{mcool}::/resolutions/{res}')      # chrom names contain ':', so slice, don't fetch()
    bins = c.bins()[:]
    m = (bins['chrom'].astype(str) == CHR3).values
    return c.matrix(balance=False)[:][np.ix_(m, m)].astype(float)


def kee_track(mcool, res):
    """per-Chr3-bin inter-chromosomal contact fraction (KEE / KNOT signature)."""
    c = cooler.Cooler(f'{mcool}::/resolutions/{res}')
    bins = c.bins()[:]
    nuc = [ch for ch in c.chromnames if ch.lower().startswith('chr')]
    mask = bins['chrom'].isin(nuc).values
    M = c.matrix(balance=False)[:][np.ix_(mask, mask)].astype(float)
    bn = bins[mask].reset_index(drop=True)
    chrom = bn['chrom'].astype(str).values
    is3 = chrom == CHR3
    inter = chrom[:, None] != chrom[None, :]
    rows = np.where(is3)[0]
    cov = np.array([M[r].sum() for r in rows])
    frac = np.array([M[r, inter[r]].sum() / max(M[r].sum(), 1) for r in rows])
    start = bn['start'].values[is3]
    return start, frac, cov


def main():
    s, e = cen_of(CHR3)
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.1], hspace=0.28, wspace=0.15)
    kee_pos = {}
    for k, (name, mc) in enumerate(MCOOLS.items()):
        ax = fig.add_subplot(gs[0, k])
        M = chr3_intra(mc, RES_MAP)
        nb = M.shape[0]
        vmax = np.percentile(M[M > 0], 99.5)
        ax.imshow(np.log1p(M), cmap='Reds', vmax=np.log1p(vmax), origin='upper', interpolation='nearest')
        cs, ce = s // RES_MAP, e // RES_MAP
        for b in (cs, ce):
            ax.axhline(b, color='#2b6cb0', lw=0.7, ls='--'); ax.axvline(b, color='#2b6cb0', lw=0.7, ls='--')
        ax.set_title(f'{name}  —  Chr3 (50 kb)', fontsize=11)
        ticks = np.arange(0, nb, 100); ax.set_xticks(ticks); ax.set_xticklabels((ticks * RES_MAP // 1_000_000));
        ax.set_yticks(ticks); ax.set_yticklabels((ticks * RES_MAP // 1_000_000)); ax.set_xlabel('Mb'); ax.set_ylabel('Mb')
        ax.text(cs + (ce - cs) / 2, -6, 'CEN', color='#2b6cb0', ha='center', fontsize=8)
        ax.text(nb * 0.15, -6, 'p arm', color='#555', ha='center', fontsize=8)
        ax.text(nb * 0.9, -6, 'q arm', color='#555', ha='center', fontsize=8)
    # KEE track (bottom, spans both) — detect p-arm peak from AlwI+BfaI averaged
    axt = fig.add_subplot(gs[1, :])
    TELO = 1_500_000; parm_lo, parm_hi = TELO, s - FLANK      # p arm, minus telomere & pericentromere
    smooth = lambda x: np.convolve(x, np.ones(3) / 3, 'same')
    peaks = {}
    for name, mc in MCOOLS.items():
        start, frac, cov = kee_track(mc, RES_TRK); fs = smooth(frac)
        axt.plot(start / 1e6, fs, label=name, lw=1.5)
        valid = (start >= parm_lo) & (start <= parm_hi) & (cov >= np.median(cov[cov > 0]))
        peaks[name] = start[valid][np.argmax(fs[valid])]
    axt.axvspan(0, TELO / 1e6, color='#888', alpha=0.10)
    axt.axvspan((s - FLANK) / 1e6, (e + FLANK) / 1e6, color='#2b6cb0', alpha=0.12)
    ymax = axt.get_ylim()[1]
    axt.text(TELO / 2e6, ymax * 0.8, 'telomere', color='#888', ha='center', fontsize=7)
    axt.text((s + e) / 2e6, ymax * 0.85, 'pericentromere\n(chromocenter)', color='#2b6cb0', ha='center', fontsize=8)
    kee = peaks['BfaI']      # anchor on the deeper, cleaner-baseline dataset
    axt.axvline(kee / 1e6, color='#d97706', lw=1.8)
    axt.annotate(f'candidate KEE ~{kee/1e6:.1f} Mb', xy=(kee / 1e6, ymax * 0.5),
                 xytext=(kee / 1e6 + 2.5, ymax * 0.75), color='#d97706', fontsize=9,
                 arrowprops=dict(arrowstyle='->', color='#d97706'))
    axt.set_xlabel('Chr3 position (Mb)'); axt.set_ylabel('inter-chromosomal\ncontact fraction')
    axt.set_title('KEE / KNOT signature: bins that contact other chromosomes (peaks outside the pericentromere = KEEs)', fontsize=10)
    axt.legend(fontsize=8); axt.set_xlim(0, 26.15)
    # mark the KEE on the maps
    for k, name in enumerate(MCOOLS):
        ax = fig.axes[k]; kb = int(kee // RES_MAP)
        ax.scatter([kb], [-3], marker='v', s=40, color='#d97706', clip_on=False)
        ax.scatter([-3], [kb], marker='>', s=40, color='#d97706', clip_on=False)
    fig.suptitle('Chr3 zoom of the real Pore-C maps + data-driven KEE detection (p arm)', fontsize=13)
    fig.savefig('chr3_zoom_AlwI_BfaI.png', dpi=140, bbox_inches='tight')
    print('candidate KEE (p arm), Mb:', {k: round(v / 1e6, 2) for k, v in peaks.items()})
    print('wrote chr3_zoom_AlwI_BfaI.png')


if __name__ == '__main__':
    main()
