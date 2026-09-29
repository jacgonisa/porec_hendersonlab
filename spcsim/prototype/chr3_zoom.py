"""
Chr3 zoom of the real AlwI & BfaI Pore-C maps, with INTRACHROMOSOMAL KEE detection.
In Arabidopsis inter-chromosomal contacts are weak; KEEs (KNOT-engaged elements) interact
largely in cis (long-range within a chromosome — KEE-KEE and to heterochromatin). So we score
each arm bin by its long-range (>4 Mb) intrachromosomal contact fraction and flag the p-arm peak,
then show that bin's virtual-4C (its cis contact profile) to reveal its partners.

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
RES = 50_000
FLANK, TELO, LR = 2_000_000, 1_500_000, 4_000_000


def cen_of(name):
    for line in open(CENBED):
        f = line.split()
        for i in range(1, len(f) - 1):
            if f[i].isdigit() and f[i + 1].isdigit() and f[i - 1] == name:
                return int(f[i]), int(f[i + 1])


def chr3_intra(mcool, res):
    c = cooler.Cooler(f'{mcool}::/resolutions/{res}')     # names contain ':', so slice not fetch()
    bins = c.bins()[:]
    m = (bins['chrom'].astype(str) == CHR3).values
    return c.matrix(balance=False)[:][np.ix_(m, m)].astype(float)


def lr_score(M, cen):
    nb = M.shape[0]; pos = np.arange(nb) * RES
    peri = (pos >= cen[0] - FLANK) & (pos <= cen[1] + FLANK)
    telo = (pos < TELO) | (pos > pos[-1] - TELO)
    lr = int(LR / RES)
    dd = np.abs(np.arange(nb)[:, None] - np.arange(nb)[None, :])
    far = dd > lr
    cov = M.sum(1)
    score = np.array([M[i, far[i]].sum() / cov[i] if cov[i] > 0 else 0 for i in range(nb)])
    return pos, score, peri, telo, cov


def main():
    cen = cen_of(CHR3)
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.2], hspace=0.3, wspace=0.15)
    mats = {}
    for k, (name, mc) in enumerate(MCOOLS.items()):
        M = chr3_intra(mc, RES); mats[name] = M; nb = M.shape[0]
        ax = fig.add_subplot(gs[0, k])
        vmax = np.percentile(M[M > 0], 99.5)
        ax.imshow(np.log1p(M), cmap='Reds', vmax=np.log1p(vmax), origin='upper', interpolation='nearest')
        cs, ce = cen[0] // RES, cen[1] // RES
        for b in (cs, ce):
            ax.axhline(b, color='#2b6cb0', lw=0.7, ls='--'); ax.axvline(b, color='#2b6cb0', lw=0.7, ls='--')
        t = np.arange(0, nb, 100); ax.set_xticks(t); ax.set_xticklabels(t * RES // 10**6)
        ax.set_yticks(t); ax.set_yticklabels(t * RES // 10**6); ax.set_xlabel('Mb'); ax.set_ylabel('Mb')
        ax.set_title(f'{name}  —  Chr3 (50 kb)', fontsize=11)
        ax.text(nb * 0.16, -6, 'p arm', color='#555', ha='center', fontsize=8)
        ax.text((cs + ce) / 2, -6, 'CEN', color='#2b6cb0', ha='center', fontsize=8)
        ax.text(nb * 0.88, -6, 'q arm', color='#555', ha='center', fontsize=8)

    # intrachromosomal long-range track + KEE detection (anchor on BfaI, deeper)
    axt = fig.add_subplot(gs[1, :]); peaks = {}
    smooth = lambda x: np.convolve(x, np.ones(3) / 3, 'same')
    for name in MCOOLS:
        pos, score, peri, telo, cov = lr_score(mats[name], cen)
        axt.plot(pos / 1e6, smooth(score), label=name, lw=1.5)
        arm = (~peri) & (~telo) & (pos < cen[0]) & (cov >= np.median(cov[cov > 0]))
        peaks[name] = pos[arm][np.argmax(smooth(score)[arm])]
    kee = peaks['BfaI']
    axt.axvspan((cen[0] - FLANK) / 1e6, (cen[1] + FLANK) / 1e6, color='#2b6cb0', alpha=0.12)
    axt.axvspan(0, TELO / 1e6, color='#888', alpha=0.10)
    ym = axt.get_ylim()[1]
    axt.text((cen[0] + cen[1]) / 2e6, ym * 0.85, 'pericentromere', color='#2b6cb0', ha='center', fontsize=8)
    kb = int(kee // RES)
    axt.axvline(kee / 1e6, color='#d97706', lw=1.5, ls='--')
    axt.annotate(f'p-arm long-range-cis peak ~{kee/1e6:.1f} Mb\n(candidate KEE — confirm with published coords)',
                 xy=(kee / 1e6, ym * 0.5), xytext=(kee / 1e6 + 2.0, ym * 0.82), color='#d97706', fontsize=8,
                 arrowprops=dict(arrowstyle='->', color='#d97706'))
    axt.set_ylabel('long-range cis\ncontact fraction'); axt.set_xlim(0, 26.15); axt.set_xlabel('Chr3 position (Mb)'); axt.legend(fontsize=8)
    axt.set_title('Intrachromosomal structure: long-range (>4 Mb) cis contact fraction per bin '
                  '(rises into the pericentromere; euchromatic KEEs are subtle at this depth)', fontsize=9.5)
    for k, name in enumerate(MCOOLS):
        fig.axes[k].scatter([kb], [-3], marker='v', s=36, color='#d97706', clip_on=False)
        fig.axes[k].scatter([-3], [kb], marker='>', s=36, color='#d97706', clip_on=False)
    fig.suptitle('Chr3 zoom of the real AlwI & BfaI Pore-C maps (KEEs interact in cis; inter-chromosomal signal is weak in Arabidopsis)', fontsize=12)
    fig.savefig('chr3_zoom_AlwI_BfaI.png', dpi=140, bbox_inches='tight')
    print('candidate KEE (p arm, cis), Mb:', {k: round(v / 1e6, 2) for k, v in peaks.items()})
    print('wrote chr3_zoom_AlwI_BfaI.png')


if __name__ == '__main__':
    main()
