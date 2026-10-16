"""
Side-by-side contact maps: unrescued (digest-first, Type IIS site-less junctions lost)
vs rescued (alignment recovers the fused monomers). Simulated AlwI Pore-C on the real
Col-0 chromosome structure (no AlwI monomer BAM exists to rescue directly).

Model: generate concatemer walks (distance decay + inter-chromosomal + pericentromere pull).
At each junction, prob 0.25 it is site-less (both ends lack GGATC) -> the two monomers FUSE
-> the digest keeps one, the partner is soft-clipped away (lost). Rescued = keep all monomers.

Run: python typeIIS_map_rescue.py [n_reads]   ->  contactmap_typeIIS_rescue.png
"""
import sys, math, random
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import quick_mc as q

BINW = 300_000
P_SITELESS = 0.25            # Type IIS: fraction of junctions with both ends site-less
P_INTER, P_PERI = 0.12, 0.15
nb = [L // BINW + 1 for L in q.CHR_LEN]
off = np.concatenate([[0], np.cumsum(nb)]); NB = off[-1]
CEN_BIN = [int(c // BINW) for c in q.CEN_MID]


def gbin(c, b): return off[c] + min(b, nb[c] - 1)


def walk(nprng):
    n = 1 + nprng.poisson(2.5)
    c = random.choices(range(5), q.CHR_W)[0]; b = random.randrange(nb[c])
    mons = [(c, b)]
    for _ in range(n - 1):
        u = random.random()
        if u < P_INTER:
            c2 = random.choice([x for x in range(5) if x != c]); b2 = random.randrange(nb[c2])
        elif u < P_INTER + P_PERI:
            c2 = c; b2 = min(nb[c] - 1, max(0, int(random.gauss(CEN_BIN[c], 4))))   # pericentromere pull
        else:
            c2 = c; d = int(math.exp(random.uniform(0, math.log(max(nb[c] / 2, 2)))))
            b2 = min(nb[c] - 1, max(0, b + random.choice((-1, 1)) * d))
        mons.append((c2, b2))
    return mons


def main():
    n_reads = int(sys.argv[1]) if len(sys.argv) > 1 else 800_000
    random.seed(1); nprng = np.random.default_rng(1)
    Mr = np.zeros((NB, NB)); Mu = np.zeros((NB, NB))
    c_resc = c_unresc = mw_resc = mw_unresc = 0
    for _ in range(n_reads):
        mons = walk(nprng)
        # Type IIS loss: walk the junctions; a site-less junction fuses -> drop the later monomer
        lost = set()
        for i in range(len(mons) - 1):
            if (i not in lost) and (i + 1 not in lost) and random.random() < P_SITELESS:
                lost.add(i + 1)
        surv = [m for i, m in enumerate(mons) if i not in lost]
        for mons_set, M, isresc in [(mons, Mr, True), (surv, Mu, False)]:
            gb = [gbin(c, b) for c, b in mons_set]
            k = len(gb)
            for i in range(k):
                for j in range(i + 1, k):
                    M[gb[i], gb[j]] += 1; M[gb[j], gb[i]] += 1
            if isresc: c_resc += k * (k - 1) // 2; mw_resc += (k >= 3)
            else: c_unresc += k * (k - 1) // 2; mw_unresc += (k >= 3)

    def draw(ax, M, title):
        vmax = np.percentile(M[M > 0], 99.5) if (M > 0).any() else 1
        ax.imshow(np.log1p(M), cmap='Reds', vmax=np.log1p(vmax), origin='upper', interpolation='nearest')
        for b in off[1:-1]:
            ax.axhline(b - .5, color='#4169b0', lw=.4, alpha=.6); ax.axvline(b - .5, color='#4169b0', lw=.4, alpha=.6)
        ax.set_xticks([(off[i] + off[i+1]) / 2 for i in range(5)]); ax.set_xticklabels([f'Chr{i+1}' for i in range(5)], fontsize=7)
        ax.set_yticks([(off[i] + off[i+1]) / 2 for i in range(5)]); ax.set_yticklabels([f'Chr{i+1}' for i in range(5)], fontsize=7)
        ax.set_title(title, fontsize=11)

    fig, ax = plt.subplots(1, 3, figsize=(16, 5.4))
    draw(ax[0], Mu, f'UNRESCUED (digest-first)\n{c_unresc/1e6:.2f}M contacts')
    draw(ax[1], Mr, f'RESCUED (alignment recovers fused monomers)\n{c_resc/1e6:.2f}M contacts')
    diff = Mr - Mu
    vmax = np.percentile(diff[diff > 0], 99) if (diff > 0).any() else 1
    ax[2].imshow(np.log1p(diff), cmap='Greens', vmax=np.log1p(vmax), origin='upper', interpolation='nearest')
    for b in off[1:-1]:
        ax[2].axhline(b - .5, color='#4169b0', lw=.4, alpha=.5); ax[2].axvline(b - .5, color='#4169b0', lw=.4, alpha=.5)
    ax[2].set_xticks([]); ax[2].set_yticks([])
    gain = 100 * (c_resc - c_unresc) / c_resc
    ax[2].set_title(f'RECOVERED contacts (rescued − unrescued)\n+{gain:.0f}% total; multiway reads +{100*(mw_resc-mw_unresc)/mw_resc:.0f}%', fontsize=11)
    fig.suptitle('Type IIS (AlwI) monomer loss and its rescue — simulated Pore-C on Col-0 (300 kb bins)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig('contactmap_typeIIS_rescue.png', dpi=140)
    print(f"contacts: unrescued {c_unresc/1e6:.2f}M  rescued {c_resc/1e6:.2f}M  (+{gain:.1f}%)")
    print(f"multiway (>=3mono) reads: unrescued {mw_unresc:,}  rescued {mw_resc:,}  (+{100*(mw_resc-mw_unresc)/mw_resc:.1f}%)")
    print('wrote contactmap_typeIIS_rescue.png')


if __name__ == '__main__':
    main()
