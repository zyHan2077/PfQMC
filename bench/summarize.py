"""Summarize the measured wall times and paired speedups in a benchmark CSV."""
import csv
import statistics as st
import sys
from collections import defaultdict

rows=list(csv.DictReader(open(sys.argv[1])))
groups=defaultdict(list)
for r in rows:
    key=tuple(r[k] for k in ['model','sites','ndim','ltau','pairs','threads','seed'])
    groups[key].append(r)
print('| Model | sites | matrix dim | Ltau | pairs | threads | seed | baseline sweep s | optimized sweep s | speedup | paired range | total speedup |')
print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|')
for key,rs in groups.items():
    b={r['repeat']:r for r in rs if r['version']=='baseline'}
    o={r['repeat']:r for r in rs if r['version']=='blocks'}
    if b.keys()!=o.keys():
        print('Incomplete pairs:',key,file=sys.stderr); continue
    bm=st.median(float(r['sweeps_s']) for r in b.values())
    om=st.median(float(r['sweeps_s']) for r in o.values())
    paired=[float(b[i]['sweeps_s'])/float(o[i]['sweeps_s']) for i in b]
    total=st.median(float(r['total_s']) for r in b.values())/st.median(float(r['total_s']) for r in o.values())
    print('| '+' | '.join(key)+f' | {bm:.3f} | {om:.3f} | {bm/om:.2f}x | {min(paired):.2f}-{max(paired):.2f}x | {total:.2f}x |')
print('\nMaximum raw signed-Pfaffian discrepancy:',max(float(r['raw_error']) for r in rows))
print('Maximum rebuilt Green relative error:',max(float(r['rebuild_error']) for r in rows))
