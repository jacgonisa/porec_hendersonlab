"""Compute + render the monomer-type breakdown slide (async 3 h pulse) -> slide_read_4.png.
Answers: 'both BrdU- => same chromatid?'  No -- most BrdU- monomers are unreplicated or are the
parental strand of *either* sister. Style matches animate_schematics.py."""
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'prototype'))
import quick_mc as q

BG = '#0e1117'; FG = '#e6edf3'; YEL = '#f2cc60'; GREY = '#6e7681'; TEAL = '#2dd4bf'; RED = '#f85149'; BLUE = '#58a6ff'
random.seed(0)
p = dict(q.BASE)   # async 3 h pulse

cats = {'BrdU+  (labelled nascent) — informative': 0,
        'BrdU−  parental strand of a labelled sister — rescuable': 0,
        'BrdU−  unlabelled sister (replicated outside pulse)': 0,
        'BrdU−  unreplicated DNA — no sister': 0}
keys = list(cats)
N = 400_000
for _ in range(N):
    cell = q.sample_cell(p)
    c = random.randrange(5); h = random.randint(0, 1); b = random.randrange(q.NB[c])
    rep, cps, _ = cell.state(c, h, b)
    cid = random.choice(list(cps)); lab = cps[cid]; rs = random.choice('WC')
    brdu = rs in lab
    if not rep:                cats[keys[3]] += 1
    elif brdu:                 cats[keys[0]] += 1
    elif len(lab) > 0:         cats[keys[1]] += 1
    else:                      cats[keys[2]] += 1
pct = {k: 100 * v / N for k, v in cats.items()}

fig, ax = plt.subplots(figsize=(7.6, 4.2), facecolor=BG); ax.set_facecolor(BG)
ax.set_xlim(0, 100); ax.set_ylim(-0.5, 4.6); ax.axis('off')
ax.text(50, 4.4, 'Monomer types in unsynchronised root tips (3 h BrdU pulse)', ha='center', color=FG, fontsize=12)
colors = [YEL, TEAL, BLUE, GREY]
for i, k in enumerate(keys):
    y = 3.3 - i
    ax.barh(y, pct[k], height=0.55, color=colors[i])
    ax.text(0, y + 0.42, k, color=colors[i], fontsize=9, va='bottom')
    ax.text(pct[k] + 1, y, f'{pct[k]:.0f}%', color=FG, va='center', fontsize=10)
brdu_minus = pct[keys[1]] + pct[keys[2]] + pct[keys[3]]
ax.text(50, -0.35,
        f'{brdu_minus:.0f}% of monomers are BrdU− — but that does NOT mean "same chromatid":\n'
        f'a BrdU− read is a parental strand of EITHER sister, or unreplicated DNA (≈50:50 same vs different if replicated).',
        ha='center', color="#3fb950", fontsize=8.5)
fig.tight_layout()
fig.savefig('slide_read_4.png', dpi=110, facecolor=BG)
print('percentages:', {k: round(v, 1) for k, v in pct.items()})
print('wrote slide_read_4.png')
