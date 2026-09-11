#!/usr/bin/env python3
"""Aggregate recorded rounds without mutating model counts or historical evidence."""
import csv,json,math
from pathlib import Path
c=Path(__file__).resolve().parent;o=c/'search_results'
keys=('rows','positive','zero','negative_controls','complex_controls')
first={k:0 for k in keys};first.update(max_relative_error=0.,max_positive_phase_error=0.)
for p in sorted(o.glob('round1_*.csv')):
 for row in csv.DictReader(p.open()):
  kind=row['kind'];key={'negative':'negative_controls','complex_control':'complex_controls','zero':'zero','search_fields':'positive','search_sweep':'positive'}[kind]
  z=complex(float(row['weight_real']),float(row['weight_imag']));ref=complex(float(row['reference_real']),float(row['reference_imag']))
  error=float(row['relative_error']);assert math.isfinite(error) and error<2e-7
  first['rows']+=1;first[key]+=1;first['max_relative_error']=max(first['max_relative_error'],error)
  if key=='positive':
   assert abs(z)>0 and ref.real>0 and abs(z/abs(z)-1)<2e-7
   first['max_positive_phase_error']=max(first['max_positive_phase_error'],abs(z/abs(z)-1))
assert first['rows']==3198 and first['positive']==3190,first
rounds=[dict(round=1,**first)]
for n in range(2,8):
 s=json.loads((o/f'round{n}_summary.json').read_text())['totals'];assert s['rows']==sum(s.get(k,0) for k in keys[1:])
 rounds.append(dict(round=n,**{k:s.get(k,0) for k in keys},max_relative_error=s['max_relative_error'],max_positive_phase_error=s['max_positive_phase_error']))
tot={k:sum(s[k] for s in rounds) for k in keys}
for k in ('max_relative_error','max_positive_phase_error'):tot[k]=max(s[k] for s in rounds)
registry=json.loads((c/'search_registry.json').read_text());count=sum(x['counts_as'] for x in registry['results'] if x['status']=='verified');assert count==registry['verified_count']
result=dict(count_scope=registry['scope'],distinct_hamiltonians=count,totals=tot,rounds=rounds)
(o/'all_rounds_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
