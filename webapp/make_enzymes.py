#!/usr/bin/env python3
"""Generate neb_enzymes.json: NEB-supplied enzymes from REBASE (via Biopython),
grouped by recognition site (isoschizomers share an entry). Sites <4 bp dropped.
Requires: pip install biopython.  Run: python3 make_enzymes.py"""
import json
from Bio import Restriction as R

BENCH = {'DpnII', 'NlaIII', 'AlwI', 'BfaI', 'MluCI', 'BbsI', 'BglII', 'AflII', 'MluI'}
IUPAC = set('ACGTRYSWKMBDHVN')

by_site = {}
for e in R.AllEnzymes:
    if 'New England Biolabs' not in e.supplier_list():
        continue
    site = str(e.site).upper() if e.site else ''
    if len(site) < 4 or any(c not in IUPAC for c in site):
        continue
    by_site.setdefault(site, set()).add(str(e))

entries = []
for site, names in by_site.items():
    names = sorted(names)
    entries.append({'site': site, 'len': len(site), 'enzymes': names,
                    'bench': [n for n in names if n in BENCH]})
entries.sort(key=lambda x: (x['len'], x['site']))

with open('neb_enzymes.json', 'w') as f:
    json.dump(entries, f, separators=(',', ':'))
print(f'{len(entries)} unique sites, '
      f'{sum(len(x["enzymes"]) for x in entries)} enzymes, '
      f'{sum(1 for x in entries if x["bench"])} benchmarked')
