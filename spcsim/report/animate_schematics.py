"""3Blue1Brown-style schematic GIFs for the sister-Pore-C report (matplotlib + pillow, no ffmpeg).
Produces: anim_labelling.gif, anim_read.gif, anim_asymmetry.gif
Run: python animate_schematics.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch, Circle

BG = '#0e1117'; FG = '#e6edf3'; BLUE = '#58a6ff'; TEAL = '#2dd4bf'; YEL = '#f2cc60'
RED = '#f85149'; GREY = '#6e7681'; GREEN = '#3fb950'; PURP = '#bc8cff'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'text.color': FG,
                     'axes.edgecolor': BG, 'figure.facecolor': BG, 'savefig.facecolor': BG})

def blank(ax, xlim=(0, 10), ylim=(0, 6)):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_facecolor(BG)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

def ease(t): return t * t * (3 - 2 * t)            # smoothstep
def fade(t, a, b): return max(0.0, min(1.0, (t - a) / (b - a)))


# ---------------------------------------------------------------- 1. labelling
def make_labelling(path, N=110, fps=20):
    fig, ax = plt.subplots(figsize=(7.2, 4.2)); blank(ax)
    def strand(y, x0, x1, color, lw=5, alpha=1, beads=False, bead_alpha=1):
        ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle='round', alpha=alpha)
        if beads:
            for x in np.arange(x0 + 0.4, x1, 0.7):
                ax.add_patch(Circle((x, y), 0.11, color=YEL, alpha=bead_alpha, zorder=5))
    def update(f):
        ax.clear(); blank(ax)
        t = f / N
        ax.text(5, 5.6, 'One S phase in BrdU  →  each sister labelled on a defined strand',
                ha='center', fontsize=13, color=FG)
        split = ease(fade(t, 0.15, 0.55))            # 0..1 vertical separation
        if t < 0.15:                                  # parental duplex
            strand(3.1, 1.5, 8.5, BLUE); strand(2.9, 1.5, 8.5, TEAL)
            ax.text(1.2, 3.1, "W (+)", color=BLUE, ha='right', va='center', fontsize=11)
            ax.text(1.2, 2.7, "C (−)", color=TEAL, ha='right', va='center', fontsize=11)
            ax.text(5, 1.5, 'parental duplex', ha='center', color=GREY, fontsize=11)
        else:
            yU, yL = 3.0 + 1.4 * split, 3.0 - 1.4 * split
            ba = ease(fade(t, 0.5, 0.85))
            # top sister pW: parental W + nascent C (BrdU)
            strand(yU + 0.12, 1.5, 8.5, BLUE)
            strand(yU - 0.12, 1.5, 8.5, TEAL, beads=True, bead_alpha=ba, alpha=0.9)
            # bottom sister pC: nascent W (BrdU) + parental C
            strand(yL + 0.12, 1.5, 8.5, BLUE, beads=True, bead_alpha=ba, alpha=0.9)
            strand(yL - 0.12, 1.5, 8.5, TEAL)
            if ba > 0.05:
                ax.text(8.8, yU, "pW\nBrdU on C (−)", color=YEL, va='center', fontsize=11, alpha=ba)
                ax.text(8.8, yL, "pC\nBrdU on W (+)", color=YEL, va='center', fontsize=11, alpha=ba)
                ax.text(1.2, yU, "sister 1", color=FG, ha='right', va='center', fontsize=10, alpha=ba)
                ax.text(1.2, yL, "sister 2", color=FG, ha='right', va='center', fontsize=10, alpha=ba)
        if t > 0.9:
            ax.text(5, 0.5, 'a BrdU+ monomer read on strand X  →  it is the sister with BrdU on X',
                    ha='center', color=GREEN, fontsize=11, alpha=ease(fade(t, 0.9, 1.0)))
        # legend
        ax.add_patch(Circle((0.5, 5.6), 0.09, color=YEL)); ax.text(0.7, 5.6, 'BrdU', color=YEL, va='center', fontsize=9)
        return []
    anim = FuncAnimation(fig, update, frames=N, interval=1000 / fps)
    anim.save(path, writer=PillowWriter(fps=fps)); plt.close(fig); print('wrote', path)


# ---------------------------------------------------------------- 2. read -> call
def make_read(path, fps=20):
    # cases: (name, sisterA, strandreadA, sisterB, strandreadB) ; nascent: pW->C, pC->W
    cases = [
        ('cis-sister', 'pW', 'C', 'pW', 'C', TEAL),      # both pW, both read nascent C -> same -> cis
        ('trans-sister', 'pW', 'C', 'pC', 'W', RED),     # opposite sisters, opposite strands -> trans
        ('uninformative', 'pW', 'W', 'pC', 'C', GREY),   # neither read the nascent strand -> no BrdU
    ]
    nascent = {'pW': 'C', 'pC': 'W'}
    per = 34; N = per * len(cases)
    fig, ax = plt.subplots(figsize=(7.2, 4.2)); blank(ax)
    def monomer(x, label, sister, readstrand, reveal, brdu_on):
        ax.add_patch(plt.Rectangle((x - 0.7, 2.6), 1.4, 0.8, color='#20262e', ec=BLUE, lw=1.5, alpha=reveal))
        ax.text(x, 3.0, label, ha='center', va='center', color=FG, fontsize=10, alpha=reveal)
        # strand-read arrow above (coin flip result)
        if reveal > 0.5:
            arr = '→ read '+readstrand if readstrand == 'W' else '← read '+readstrand
            ax.text(x, 4.0, arr, ha='center', color=YEL if brdu_on else GREY, fontsize=11)
            if brdu_on:
                ax.add_patch(Circle((x, 3.55), 0.12, color=YEL, zorder=6))
                ax.text(x, 4.4, 'BrdU+', ha='center', color=YEL, fontsize=9)
            else:
                ax.text(x, 4.4, 'BrdU−', ha='center', color=GREY, fontsize=9)
    def update(f):
        ax.clear(); blank(ax)
        ci = min(len(cases) - 1, f // per); lf = (f % per) / per
        name, sA, rA, sB, rB, col = cases[ci]
        ax.text(5, 5.5, 'Reading a Pore-C concatemer: strand + BrdU  →  sister call', ha='center', fontsize=13)
        ax.plot([1, 9], [3.0, 3.0], color='#30363d', lw=2, zorder=0)  # genome axis
        brA = (rA == nascent[sA]); brB = (rB == nascent[sB])
        monomer(3.2, 'monomer i', sA, rA, ease(fade(lf, 0.0, 0.3)), brA)
        monomer(6.8, 'monomer j', sB, rB, ease(fade(lf, 0.15, 0.45)), brB)
        if lf > 0.55:
            a = ease(fade(lf, 0.55, 0.8))
            ax.add_patch(FancyArrowPatch((3.2, 2.55), (6.8, 2.55), connectionstyle='arc3,rad=-0.35',
                         color=col, lw=2.5, alpha=a, arrowstyle='-'))
            if name == 'uninformative':
                verdict = 'not both BrdU+  →  cannot assign sisters'
            else:
                same = (rA == rB)
                verdict = f'both BrdU+, {"same" if same else "opposite"} strand  →  {name}'
            ax.text(5, 1.5, verdict, ha='center', color=col, fontsize=13, alpha=a)
        ax.text(0.6, 4.0, 'nanopore reads\none random strand', color=GREY, fontsize=8, va='center')
        return []
    anim = FuncAnimation(fig, update, frames=N, interval=1000 / fps)
    anim.save(path, writer=PillowWriter(fps=fps)); plt.close(fig); print('wrote', path)


# ---------------------------------------------------------------- 3. asymmetry / register shift
def make_asymmetry(path, N=120, fps=20, delta=3.0, sigma=8.0):
    rng = np.random.default_rng(1)
    Kmax = 4000
    oriented = rng.normal(delta, sigma, Kmax)                 # signed offsets j-i (pW->pC frame)
    bulk = oriented * rng.choice([-1, 1], Kmax)               # no strand -> symmetrised
    bins = np.arange(-40, 41, 2)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.4, 4.3), gridspec_kw={'width_ratios': [1, 1.25]})
    def update(f):
        for ax in (axL, axR): ax.clear()
        blank(axL, (0, 10), (0, 10))
        t = f / N; k = int(ease(t) * Kmax)
        # left: two staggered sisters + a contact
        axL.text(5, 9.3, 'sisters are staggered by δ', ha='center', fontsize=12)
        for y, col, lab in [(6.2, TEAL, 'sister pW'), (4.2, BLUE, 'sister pC')]:
            axL.plot([1.5, 8.5], [y, y], color=col, lw=6, solid_capstyle='round')
            axL.text(1.2, y, lab, ha='right', va='center', color=col, fontsize=9)
        sh = delta * 0.18
        i0 = 4.3
        axL.add_patch(Circle((i0, 6.2), 0.16, color=YEL, zorder=5))
        axL.add_patch(Circle((i0 + sh + 1.4, 4.2), 0.16, color=YEL, zorder=5))
        axL.annotate('', xy=(i0 + sh + 1.4, 4.2), xytext=(i0, 6.2),
                     arrowprops=dict(arrowstyle='-', color=RED, lw=2))
        axL.text(5, 2.6, 'trans-sister contact:\ni on pW  ↔  j on pC,  with j = i + δ',
                 ha='center', color=FG, fontsize=10)
        axL.text(5, 1.1, 'strand read tells which sister  →  we can SIGN j−i',
                 ha='center', color=GREEN, fontsize=10)
        # right: accumulating signed-offset histogram
        axR.set_facecolor(BG)
        for s in axR.spines.values(): s.set_color('#30363d')
        axR.tick_params(colors=GREY, labelsize=8)
        if k > 5:
            axR.hist(bulk[:k], bins=bins, density=True, histtype='stepfilled',
                     color=GREY, alpha=0.35, label='bulk (no strand)')
            axR.hist(oriented[:k], bins=bins, density=True, histtype='step',
                     color=BLUE, lw=2.2, label='strand-oriented')
            m = oriented[:k].mean()
            axR.axvline(0, color=GREY, ls='--', lw=1)
            axR.axvline(m, color=YEL, lw=2)
            axR.text(m, axR.get_ylim()[1]*0.92, f' mean = {m*10:.0f} kb', color=YEL, fontsize=10)
        axR.axvline(delta, color=GREEN, ls=':', lw=1.5)
        axR.set_xlim(-40, 40); axR.set_xlabel('signed offset  j − i  (10 kb bins)', color=FG, fontsize=10)
        axR.set_ylabel('density', color=FG, fontsize=10)
        axR.set_title(f'reads accumulated: {k:,}', color=FG, fontsize=11)
        axR.legend(fontsize=8, framealpha=0, labelcolor=FG, loc='upper left')
        return []
    anim = FuncAnimation(fig, update, frames=N, interval=1000 / fps)
    fig.tight_layout()
    anim.save(path, writer=PillowWriter(fps=fps)); plt.close(fig); print('wrote', path)


if __name__ == '__main__':
    make_labelling('anim_labelling.gif')
    make_read('anim_read.gif')
    make_asymmetry('anim_asymmetry.gif')
    print('done')
