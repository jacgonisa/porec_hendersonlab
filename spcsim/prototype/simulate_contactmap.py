"""
Simulated Pore-C contact maps for the report.
  Fig A: genome-wide bulk map (5 Arabidopsis chromosomes) with distance decay,
         pericentromere heterochromatin blocks, chromocenter (inter-chromosomal
         pericentromere) clustering, and KNOT/KEE dots.
  Fig B: sister decomposition of Chr1 (bulk | cis-sister | trans-sister).

Phenomenological expected matrices + Poisson sampling (a "simulated" map at a given depth).
Coordinates TAIR10 (from quick_mc). KEE positions are illustrative. Run: python simulate_contactmap.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import quick_mc as q

CMAP = 'Reds'
CEN = q.CEN_MID
# illustrative KEE loci (chr index 0-4, position Mb) — the Arabidopsis KNOT engages ~10 such islands
KEES = [(0, 5.0), (0, 24.0), (1, 4.0), (1, 15.0), (2, 8.0), (2, 21.0),
        (3, 11.0), (3, 16.0), (4, 6.0), (4, 23.0)]


def genome_bins(binw):
    nb = [L // binw + 1 for L in q.CHR_LEN]
    off = np.concatenate([[0], np.cumsum(nb)])
    NB = off[-1]
    chrid = np.zeros(NB, int); posbp = np.zeros(NB)
    for c in range(5):
        idx = np.arange(nb[c])
        chrid[off[c]:off[c + 1]] = c
        posbp[off[c]:off[c + 1]] = (idx + 0.5) * binw
    return nb, off, NB, chrid, posbp


def expected_genome(binw=200_000):
    nb, off, NB, chrid, posbp = genome_bins(binw)
    same = chrid[:, None] == chrid[None, :]
    d = np.abs(posbp[:, None] - posbp[None, :])
    peri = np.abs(posbp - CEN_arr(chrid)) < 2.5e6
    # distance decay within chromosome (~ s^-1 with a short-range plateau)
    E = np.where(same, 1.0 / (1.0 + (d / 3e5) ** 1.1), 0.0)
    E += (~same) * 0.02                                   # inter-chromosomal background
    pp = peri[:, None] & peri[None, :]
    E += (same & pp) * 2.5                                # pericentromere heterochromatin block
    E += (~same & pp) * 1.1                               # chromocenter: inter-chr pericentromere clustering
    # KNOT / KEE dots (intra- and inter-chromosomal)
    kee = np.zeros(NB, bool)
    for c, mb in KEES:
        kee[off[c] + int(mb * 1e6 // binw)] = True
    kp = kee[:, None] & kee[None, :]
    E += kp * 3.0
    np.fill_diagonal(E, E.max())
    return E, nb, off, NB, chrid, kee, peri


def CEN_arr(chrid):
    return np.array([CEN[c] for c in chrid])


def poisson_sample(E, total, seed=0):
    rng = np.random.default_rng(seed)
    P = E / E.sum()
    return rng.poisson(P * total)


def draw_map(ax, M, title, boundaries=None, labels=None, kee=None):
    im = ax.imshow(np.log1p(M), cmap=CMAP, origin='upper', interpolation='nearest')
    ax.set_title(title, fontsize=11)
    if boundaries is not None:
        for b in boundaries[1:-1]:
            ax.axhline(b - 0.5, color='#4169b0', lw=0.6, alpha=0.7)
            ax.axvline(b - 0.5, color='#4169b0', lw=0.6, alpha=0.7)
    if labels is not None:
        cent = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(boundaries) - 1)]
        ax.set_xticks(cent); ax.set_xticklabels(labels, fontsize=8)
        ax.set_yticks(cent); ax.set_yticklabels(labels, fontsize=8)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    if kee is not None:
        k = np.where(kee)[0]
        ax.scatter(k, [-3] * len(k), marker='v', s=18, color='#2dd4bf', clip_on=False)
        ax.scatter([-3] * len(k), k, marker='>', s=18, color='#2dd4bf', clip_on=False)
    return im


# ---------------- Fig A: genome-wide bulk map ----------------
def fig_genome():
    E, nb, off, NB, chrid, kee, peri = expected_genome(200_000)
    M = poisson_sample(E, total=6_000_000)
    fig, ax = plt.subplots(figsize=(7.6, 7.0))
    draw_map(ax, M, 'Simulated bulk Pore-C map — Arabidopsis (200 kb bins)',
             boundaries=off, labels=[f'Chr{c+1}' for c in range(5)], kee=kee)
    # annotate features
    ax.text(0.99, 0.01, 'teal ▾ = KEE (KNOT)   blue lines = chromosome borders',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8, color='#333')
    ax.text(off[0] + nb[0] * 0.62, off[0] + nb[0] * 0.62, 'pericentromere\nblocks', fontsize=7,
            color='#111', ha='center')
    fig.tight_layout(); fig.savefig('contactmap_genome.png', dpi=140); plt.close(fig)
    print('wrote contactmap_genome.png')


# ---------------- Fig B: Chr1 sister decomposition ----------------
def expected_chr1(binw=100_000):
    L = q.CHR_LEN[0]; nb = L // binw + 1
    pos = (np.arange(nb) + 0.5) * binw
    d = np.abs(pos[:, None] - pos[None, :])
    peri = np.abs(pos - CEN[0]) < 2.5e6
    kee = np.zeros(nb, bool)
    for c, mb in KEES:
        if c == 0:
            kee[int(mb * 1e6 // binw)] = True
    # a couple of arm TAD-like domains + a loop dot, for texture
    decay = 1.0 / (1.0 + (d / 2e5) ** 1.1)
    cis = decay.copy()
    cis += (peri[:, None] & peri[None, :]) * 3.0                       # pericentromere block
    kp = kee[:, None] & kee[None, :]; cis += kp * 3.5                  # KEE dots
    for (a, b) in [(20, 45), (60, 90), (150, 175)]:            # arm domains
        cis[a:b, a:b] += 0.8
    np.fill_diagonal(cis, cis.max())
    # trans-sister: aligned sisters -> near-diagonal band (width ~sigma) + faint peri cohesion + faint KEE-trans
    delta = 0.3  # bins (30 kb) sub-bin stagger, essentially on-diagonal at 100 kb
    band = np.exp(-((d / binw - delta) ** 2) / (2 * (1.5 ** 2)))       # aligned band
    trans = 0.6 * band
    trans += (peri[:, None] & peri[None, :]) * 0.9                     # partial sister cohesion at CEN
    trans += kp * 0.7                                                  # some KEEs engage across sisters
    return cis, trans, peri, kee


def fig_sister():
    cis, trans, peri, kee = expected_chr1(100_000)
    bulk = cis + trans
    total = 3_000_000
    Mb = poisson_sample(bulk, total, 1)
    Mc = poisson_sample(cis, int(total * cis.sum() / bulk.sum()), 2)
    Mt = poisson_sample(trans, int(total * trans.sum() / bulk.sum()), 3)
    fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))
    draw_map(ax[0], Mb, 'Chr1 bulk (no sister resolution)', kee=kee)
    draw_map(ax[1], Mc, 'cis-sister (within one chromatid)', kee=kee)
    draw_map(ax[2], Mt, 'trans-sister (between the two sisters)', kee=kee)
    for a in ax:
        cb = int(CEN[0] // 100_000)
        a.add_patch(plt.Rectangle((cb - 25, cb - 25), 50, 50, fill=False, ec='#2b6cb0', lw=1.0, ls='--'))
    ax[2].text(0.5, -0.08, 'trans-sister is dominated by the aligned near-diagonal band; '
               'pericentromere cohesion & some KEEs bridge sisters',
               transform=ax[2].transAxes, ha='center', fontsize=8, color='#333')
    fig.suptitle('Sister-Pore-C decomposition of the Chr1 contact map (100 kb bins, simulated)', fontsize=12)
    fig.tight_layout(rect=[0, 0.02, 1, 0.96]); fig.savefig('contactmap_sister_decomp.png', dpi=140); plt.close(fig)
    print('wrote contactmap_sister_decomp.png')


if __name__ == '__main__':
    fig_genome()
    fig_sister()
    print('done')
