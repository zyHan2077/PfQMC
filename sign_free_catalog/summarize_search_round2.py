#!/usr/bin/env python3
"""Summarize immutable round-two evidence; never modify the search registry."""
import csv,json,math
from pathlib import Path
root=Path(__file__).resolve().parent
files=sorted((root/'search_results').glob('round2_*.csv'))
assert len(files)==24, f'Expected 24 evidence files, found {len(files)}'
stats=dict(files=len(files),rows=0,positive=0,zero=0,negative_controls=0,complex_controls=0,max_relative_error=0.,max_positive_phase_error=0.)
per_file=[]
for path in files:
 local=dict(file=path.name,rows=0,positive=0,zero=0,negative_controls=0,complex_controls=0)
 for row in csv.DictReader(path.open()):
  if 'expected_sign' in row:
   sign=int(row['expected_sign']);z=complex(float(row['real']),float(row['imag']));error=float(row['relative_error'])
  elif 'expected' in row:
   sign=int(row['expected']);z=complex(float(row['weight_real']),float(row['weight_imag']));error=float(row['scaled_error'])
  else:
   z=complex(float(row['weight_real']),float(row['weight_imag']));error=float(row['relative_error'])
   kind=row['kind'];sign=0 if kind=='zero' else -1 if kind=='negative' else 2 if kind=='complex_control' else 1
  assert math.isfinite(error) and error<3e-7,(path,row)
  key={1:'positive',0:'zero',-1:'negative_controls',2:'complex_controls'}[sign]
  stats[key]+=1;local[key]+=1;stats['rows']+=1;local['rows']+=1
  stats['max_relative_error']=max(stats['max_relative_error'],error)
  if sign==1:
   assert abs(z)>1e-11 and abs(z/abs(z)-1)<3e-7,(path,row)
   stats['max_positive_phase_error']=max(stats['max_positive_phase_error'],abs(z/abs(z)-1))
 per_file.append(local)
output=dict(totals=stats,files=per_file)
(root/'search_results/round2_summary.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(stats,indent=2))
