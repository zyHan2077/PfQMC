"""Sequential, interleaved wall-clock comparisons of identical seeded sweeps."""
import argparse
import csv
import math
import statistics
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--versions', nargs='+', default=['baseline','blocks'])
parser.add_argument('--suite', choices=['pilot','large'], default='large')
parser.add_argument('--repeats', type=int, default=3)
args = parser.parse_args()
# MODEL, L, Ltau, timed right+left pairs, MKL threads, seed
cases = [('majorana',8,12,8,1,42), ('majorana',12,12,8,1,42),
         ('honeycomb',8,12,8,1,42)]
if args.suite == 'large':
    cases += [('majorana',16,12,4,1,42), ('square',12,12,4,1,42),
              ('chain0',128,12,4,1,42), ('chain1',128,12,4,1,42),
              ('majorana',12,12,4,4,923), ('majorana',8,40,4,1,123)]
header = 'model sites ndim ltau pairs threads seed init_s sweeps_s total_s sign_re sign_im checksum fields raw_error rebuild_error'.split()
rows=[]
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('w',newline='') as f:
    writer=csv.writer(f,lineterminator="\n"); writer.writerow(['version','repeat']+header)
    for repeat in range(args.repeats):
        for case in cases:
            reference=None
            for version in (args.versions if repeat%2==0 else args.versions[::-1]):
                cmd=['build/perf/pfqmc_bench_'+version,*map(str,case),'1']
                run=subprocess.run(cmd,text=True,capture_output=True)
                if run.returncode:
                    raise RuntimeError(f'{cmd}: {run.stdout}\n{run.stderr}')
                values=next(csv.reader([run.stdout.strip()]))
                row=dict(zip(header,values))
                for name in header[7:]:
                    if not math.isfinite(float(row[name])): raise RuntimeError(row)
                if reference:
                    assert row['fields']==reference['fields'], (reference,row)
                    for name in ['checksum','sign_re','sign_im']:
                        assert math.isclose(float(row[name]),float(reference[name]),abs_tol=1e-7,rel_tol=1e-9), (reference,row)
                reference=row
                writer.writerow([version,repeat]+values); f.flush()
                rows.append(dict(version=version,repeat=repeat,**row))
                print(version,repeat,case,'sweeps_s='+row['sweeps_s'],flush=True)
for case in cases:
    print('CASE',case)
    for version in args.versions:
        selected=[r for r in rows if r['version']==version and r['model']==case[0]
                  and int(r['ltau'])==case[2] and int(r['pairs'])==case[3]
                  and int(r['threads'])==case[4] and int(r['seed'])==case[5]
                  and int(r['sites'])==(2*case[1]**2 if case[0] in ('majorana','honeycomb') else case[1]**2 if case[0]=='square' else case[1])]
        print(version,{k:statistics.median(float(r[k]) for r in selected) for k in ['init_s','sweeps_s','total_s']})
