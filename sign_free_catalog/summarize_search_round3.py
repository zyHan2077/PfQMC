#!/usr/bin/env python3
"""Check and summarize the third round without rewriting earlier evidence."""
import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent;files=sorted((root/'search_results').glob('round3_*.csv'))
assert len(files)==4
stats=dict(rows=0,positive=0,zero=0,negative_controls=0,max_relative_error=0.,max_positive_phase_error=0.)
per=[]
for path in files:
 local=dict(file=path.name,rows=0,positive=0,zero=0,negative_controls=0)
 for row in csv.DictReader(path.open()):
  if 'expected_sign' in row:
   sign=int(row['expected_sign']);z=complex(float(row['real']),float(row['imag']))
  else:
   sign=-1 if row['kind']=='negative' else 1;z=complex(float(row['weight_real']),float(row['weight_imag']))
  error=float(row['relative_error']);assert math.isfinite(error) and error<2e-7
  key={1:'positive',0:'zero',-1:'negative_controls'}[sign]
  stats[key]+=1;local[key]+=1;stats['rows']+=1;local['rows']+=1
  stats['max_relative_error']=max(stats['max_relative_error'],error)
  if sign==1:
   assert abs(z)>0 and abs(z/abs(z)-1)<2e-7
   stats['max_positive_phase_error']=max(stats['max_positive_phase_error'],abs(z/abs(z)-1))
 per.append(local)
(root/'search_results/round3_summary.json').write_text(json.dumps(dict(totals=stats,files=per),indent=2)+'\n')
print(json.dumps(stats,indent=2))
