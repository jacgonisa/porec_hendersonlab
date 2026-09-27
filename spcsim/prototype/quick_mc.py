"""
Back-of-envelope Monte Carlo: sister-Pore-C in unsynchronised Arabidopsis root tips.
Truth-level simulation (no signal): cell -> chromatids -> concatemer monomers -> strand read -> BrdU call.
Purpose: estimate yields / confounds to size the real simulation, not to be the final tool.

Run:   python quick_mc.py 150000        (reads per scenario; ~2 min for all 9 scenarios)
Output: table to stdout + quick_mc_results.csv

Strand convention (the whole method rests on this):
  chromatid 'pW' carries the parental Watson (+) strand; if replicated during the pulse its nascent
  strand, Crick (-), carries BrdU. Chromatid 'pC' is the reverse (BrdU on +). Identity is continuous
  along a chromosome (no SCE). Each monomer is read on a random strand (random ligation orientation).
  BrdU+ monomer aligned to strand X  -> chromatid with BrdU on X.
  Two BrdU+ monomers, same aligned strand -> cis-sister; opposite strands -> trans-sister.

Known simplifications (fixed in the full simulator, see PLAN.md): replication = timing profile + noise
(no explicit origins/forks); no endocycling cells; phenomenological contact model with made-up
q_trans / inter-homolog rates; no mappability or inverted-repeat mis-mapping; monomer lengths drawn
from a lognormal rather than an in-silico digest; per-read O(n^2) Python loops (fine for 1e5 reads only).
"""
import math, random, sys
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

BIN = 10_000
CHR_LEN = [30_427_671, 19_698_289, 23_459_830, 18_585_056, 26_975_502]   # TAIR10
CEN_MID = [15.1e6, 3.6e6, 13.8e6, 3.9e6, 11.8e6]                        # approx.
NB = [L // BIN for L in CHR_LEN]
CHR_W = list(np.array(CHR_LEN) / sum(CHR_LEN))
DOM = 50  # bins per 500-kb noise domain


def make_timing(seed=0):
    g = np.random.default_rng(seed)
    tim = []
    for c, nb in enumerate(NB):
        x = gaussian_filter1d(g.normal(size=nb), 30)
        r = (np.argsort(np.argsort(x)) + 0.5) / nb
        pos = (np.arange(nb) + 0.5) * BIN
        peri = np.abs(pos - CEN_MID[c]) < 2.5e6
        r = np.where(peri, 0.6 + 0.4 * r, 0.75 * r)   # arms early-mid, pericentromeres mid-late
        tim.append(r)
    return tim


TIM = make_timing()


def lognorm(mean, sd):
    s2 = math.log(1 + (sd / mean) ** 2)
    return random.lognormvariate(math.log(mean) - s2 / 2, math.sqrt(s2))


def draw_cycle(p):
    G2 = lognorm(p['G2'], p['G2sd'])
    S, M = p['S'], p['M']
    G1 = max(0.5, p['T'] - S - G2 - M)
    return G1, S, G2, M


class Cell:
    def __init__(self, p):
        self.p = p
        self.dom_noise = {}
        self.dom_noise_m = {}
        self.inherit = {}
        self.cyc = p.get('ideal', False) or (random.random() < p['f_prolif'])
        if not self.cyc:
            self.content = 1.0
            return
        self.G1, self.S, self.G2, self.M = draw_cycle(p)
        Tc = self.G1 + self.S + self.G2 + self.M
        if p.get('ideal'):
            self.a = self.G1 + self.S + random.random() * self.G2
        elif random.random() < p.get('sync_frac', 0.0):
            self.a = self.G1 + p['P'] + p['C']          # released into S at pulse start
        else:
            u = random.random()
            self.a = -Tc * math.log2(1 - u / 2)          # exponential-population age density
        self.mother = draw_cycle(p)
        fr = min(1.0, max(0.0, (self.a - self.G1) / self.S))
        self.content = 1.0 + fr

    def tau(self, c, b, mother=False):
        d = self.dom_noise_m if mother else self.dom_noise
        k = (c, b // DOM)
        if k not in d:
            d[k] = random.gauss(0, 0.07)
        return min(1.0, max(0.0, TIM[c][b] + d[k] + random.gauss(0, 0.02)))

    def state(self, c, h, b):
        """returns (replicated, {chromatid_id: frozenset(labelled strands)}, 'g1d' flag)"""
        if not self.cyc:
            return False, {'u': frozenset()}, False
        p = self.p
        tau = self.tau(c, b)
        if p.get('ideal'):
            return True, {'pW': frozenset('C'), 'pC': frozenset('W')}, False
        P, C = p['P'], p['C']
        t_rep = -self.a + self.G1 + tau * self.S
        rep = t_rep <= 0
        in_pulse = (-C - P) <= t_rep <= -C
        L0 = frozenset()
        if (-C - P) < -self.a:   # pulse began before this cell was born -> check mother
            G1m, Sm, G2m, Mm = self.mother
            t_rep_m = -self.a - Mm - G2m - Sm + self.tau(c, b, mother=True) * Sm
            if (-C - P) <= t_rep_m <= -C:
                k = (c, h)
                if k not in self.inherit:
                    self.inherit[k] = random.choice('WC')
                L0 = frozenset(self.inherit[k])
        if rep:
            pW = (L0 & {'W'}) | ({'C'} if in_pulse else set())
            pC = (L0 & {'C'}) | ({'W'} if in_pulse else set())
            return True, {'pW': frozenset(pW), 'pC': frozenset(pC)}, False
        return False, {'u': L0}, bool(L0)


def sample_cell(p):
    while True:
        cell = Cell(p)
        if p.get('gate') is not None and cell.content < p['gate']:
            continue
        if random.random() < cell.content / 2.0:     # reads drawn proportional to DNA mass
            return cell


def pick_d(c):
    return int(math.exp(random.uniform(0, math.log(NB[c] / 2)))) - 1


def clipb(c, b):
    if b < 0:
        b = -b
    if b >= NB[c]:
        b = 2 * (NB[c] - 1) - b
    return max(0, min(NB[c] - 1, b))


def simulate(p, n_reads, seed=1):
    random.seed(seed)
    nprng = np.random.default_rng(seed)
    acc = dict(reads=0, mono=0, mono_L=0, reads_L1=0, reads_LL=0,
               LL=0, LL_correct_samehom=0, LL_samehom=0, LL_hom=0, LL_g1d=0, LL_dbl=0,
               resc_cand=0, resc_ok=0, resc_pairs=0, resc_pairs_correct=0, resc_pairs_samehom=0, resc_pairs_hom=0,
               allreads_pairs=0, allreads_correct=0,
               pairs_dt_lo=0, LL_dt_lo=0, pairs_dt_hi=0, LL_dt_hi=0, intra_pairs=0)
    q_trans, p_inter, p_hom = p['q_trans'], p['p_inter'], p['p_hom']
    for _ in range(n_reads):
        cell = sample_cell(p)
        n = 1 + nprng.poisson(p['mean_mono'] - 1)
        # anchor, weighted by local copy number
        while True:
            c = random.choices(range(5), CHR_W)[0]
            h = random.randint(0, 1)
            b = random.randrange(NB[c])
            rep, cps, g1d = cell.state(c, h, b)
            if random.random() < len(cps) / 2:
                break
        cid = random.choice(list(cps))
        mons = [(c, h, b, cid, cps[cid], rep, g1d)]
        for _k in range(n - 1):
            u = random.random()
            if u < p_inter:
                c2 = random.choice([x for x in range(5) if x != c]); h2 = random.randint(0, 1)
                b2 = random.randrange(NB[c2])
                rep2, cps2, g1d2 = cell.state(c2, h2, b2)
                cid2 = random.choice(list(cps2))
            elif u < p_inter + p_hom:
                c2, h2 = c, 1 - h
                b2 = clipb(c, b + random.choice((-1, 1)) * pick_d(c))
                rep2, cps2, g1d2 = cell.state(c2, h2, b2)
                cid2 = random.choice(list(cps2))
            else:
                c2, h2 = c, h
                b2 = clipb(c, b + random.choice((-1, 1)) * pick_d(c))
                rep2, cps2, g1d2 = cell.state(c2, h2, b2)
                if rep2 and rep:
                    cid2 = ({'pW': 'pC', 'pC': 'pW'}[cid] if random.random() < q_trans else cid)
                else:
                    cid2 = random.choice(list(cps2))
            mons.append((c2, h2, b2, cid2, cps2[cid2], rep2, g1d2))
        # strand read + BrdU calling per monomer
        m = len(mons)
        L = nprng.lognormal(math.log(750), 0.45, m).astype(int) + 50
        nT = nprng.binomial(L, 0.32)
        rs = [random.choice('WC') for _ in range(m)]
        truelab = np.array([rs[i] in mons[i][4] for i in range(m)])
        nB = np.where(truelab, nprng.binomial(nT, p['s']), 0)
        calls = nprng.binomial(nB, p['sens']) + nprng.binomial(nT - nB, p['fp'])
        called = calls >= p['thr']
        acc['reads'] += 1; acc['mono'] += m; acc['mono_L'] += int(called.sum())
        acc['reads_L1'] += int(called.any())
        # sister identity implied by (strand, label): label on strand X => chromatid with BrdU on X
        ident = [rs[i] if called[i] else ('W' if rs[i] == 'C' else 'C') for i in range(m)]
        Lidx = [i for i in range(m) if called[i]]
        # rescue: unlabelled monomers within 100 kb (cis) of a labelled monomer in same read
        resc = set()
        for i in range(m):
            if called[i] or not Lidx:
                continue
            if any(mons[j][0] == mons[i][0] and abs(mons[j][2] - mons[i][2]) <= 10 for j in Lidx):
                resc.add(i)
                acc['resc_cand'] += 1
                lab = mons[i][4]
                acc['resc_ok'] += int(len(lab) == 1 and rs[i] not in lab and mons[i][5])
        has_LL = False
        for i in range(m):
            for j in range(i + 1, m):
                if mons[i][0] != mons[j][0]:
                    continue
                acc['intra_pairs'] += 1
                samehom = mons[i][1] == mons[j][1]
                if samehom:
                    truth = 'trans' if {mons[i][3], mons[j][3]} == {'pW', 'pC'} else 'cis'
                call = 'cis' if ident[i] == ident[j] else 'trans'
                both_L = called[i] and called[j]
                # timing co-labelling bias (same homolog, cycling cells)
                if samehom and cell.cyc:
                    dt = abs(TIM[mons[i][0]][mons[i][2]] - TIM[mons[j][0]][mons[j][2]])
                    if dt < 0.1:
                        acc['pairs_dt_lo'] += 1; acc['LL_dt_lo'] += int(both_L)
                    elif dt > 0.4:
                        acc['pairs_dt_hi'] += 1; acc['LL_dt_hi'] += int(both_L)
                # naive 'all_reads' rule applied to every intra pair
                if samehom:
                    acc['allreads_pairs'] += 1
                    acc['allreads_correct'] += int(call == truth)
                if both_L:
                    has_LL = True
                    acc['LL'] += 1
                    if samehom:
                        acc['LL_samehom'] += 1; acc['LL_correct_samehom'] += int(call == truth)
                    else:
                        acc['LL_hom'] += 1
                    acc['LL_g1d'] += int(mons[i][6] or mons[j][6])
                    acc['LL_dbl'] += int(len(mons[i][4]) == 2 or len(mons[j][4]) == 2)
                elif (called[i] or i in resc) and (called[j] or j in resc):
                    acc['resc_pairs'] += 1
                    if samehom:
                        acc['resc_pairs_samehom'] += 1; acc['resc_pairs_correct'] += int(call == truth)
                    else:
                        acc['resc_pairs_hom'] += 1
        acc['reads_LL'] += int(has_LL)
    R = acc['reads']
    sd = lambda a, b: (a / b) if b else float('nan')
    return {
        '% monomers BrdU+': 100 * acc['mono_L'] / acc['mono'],
        '% reads >=1 BrdU+ monomer': 100 * acc['reads_L1'] / R,
        '% reads with >=1 L-L cis-chrom pair': 100 * acc['reads_LL'] / R,
        'L-L sister-informative pairs / 1M reads': 1e6 * acc['LL'] / R,
        'L-L call accuracy (same homolog) %': 100 * sd(acc['LL_correct_samehom'], acc['LL_samehom']),
        'L-L pairs that are inter-homolog %': 100 * sd(acc['LL_hom'], acc['LL']),
        'L-L pairs from post-mitotic G1 daughters %': 100 * sd(acc['LL_g1d'], acc['LL']),
        'L-L pairs with doubly-labelled chromatid %': 100 * sd(acc['LL_dbl'], acc['LL']),
        'Rescue (<=100kb) precision %': 100 * sd(acc['resc_ok'], acc['resc_cand']),
        'Extra pairs from rescue / 1M reads': 1e6 * acc['resc_pairs'] / R,
        'Rescued-pair call accuracy (same homolog) %': 100 * sd(acc['resc_pairs_correct'], acc['resc_pairs_samehom']),
        'Co-labelling: LL rate early x late / same-timing': sd(sd(acc['LL_dt_hi'], acc['pairs_dt_hi']), sd(acc['LL_dt_lo'], acc['pairs_dt_lo'])),
        '_allreads_acc': 100 * sd(acc['allreads_correct'], acc['allreads_pairs']),
    }


BASE = dict(T=17.0, S=2.5, G2=3.5, G2sd=1.0, M=0.5, f_prolif=0.6, P=3.0, C=0.0,
            s=0.20, sens=0.6, fp=0.001, thr=4, mean_mono=3.5,
            q_trans=0.25, p_inter=0.12, p_hom=0.03)

SCEN = {
    'A  Ideal G2 (HeLa-like reference)': dict(ideal=True),
    'B  Unsync, 3 h pulse (your plan)': dict(),
    'C  Unsync, 1 h pulse': dict(P=1.0),
    'D  Unsync, 6 h pulse': dict(P=6.0),
    'E  Unsync, 3 h pulse + 3 h chase': dict(C=3.0),
    'F  Unsync, 3 h pulse, low BrdU (s=5%)': dict(s=0.05),
    'G  3 h pulse + sort >=3.7C nuclei': dict(gate=1.85),
    'H  3 h pulse + HU release (40% synced)': dict(sync_frac=0.4),
    'I  3 h pulse, faster cycle (T=10 h)': dict(T=10.0, S=2.0, G2=3.0),
}

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 60000
    rows = {}
    for name, over in SCEN.items():
        p = dict(BASE); p.update(over)
        r = simulate(p, n)
        r['Naive all_reads accuracy (same homolog) %'] = r.pop('_allreads_acc')
        rows[name] = r
        print(name, 'done', file=sys.stderr, flush=True)
    df = pd.DataFrame(rows).T
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    df.to_csv('quick_mc_results.csv')
    print(df.round(2).T.to_string())
