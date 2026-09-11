#!/usr/bin/env python3
"""Validate and summarize fifth-round CSVs without changing the registry."""
import csv,json,math
from pathlib import Path
r=Path(__file__).resolve().parent;o=r/'search_results';files=sorted(o.glob('round5_*.csv'));assert len(files)==8
stats=dict(rows=0,positive=0,zero=0,negative_controls=0,max_relative_error=0.,max_positive_phase_error=0.);per=[]
for p in files:
 local=dict(file=p.name,rows=0,positive=0,zero=0,negative_controls=0)
 for row in csv.DictReader(p.open()):
  sign=int(row.get('expected',row.get('expected_sign')))
  z=complex(float(row.get('weight_real',row.get('real'))),float(row.get('weight_imag',row.get('imag'))))
  ref=complex(float(row['reference_real']),float(row['reference_imag']))
  error=float(row.get('relative_error',row.get('scaled_error')))
  assert all(math.isfinite(x) for x in (z.real,z.imag,ref.real,ref.imag,error)) and error<2e-7
  key={1:'positive',0:'zero',-1:'negative_controls'}[sign]
  stats['rows']+=1;local['rows']+=1;stats[key]+=1;local[key]+=1;stats['max_relative_error']=max(stats['max_relative_error'],error)
  if sign:
   assert abs(z)>0 and abs(z/abs(z)-sign)<2e-7 and ref.real*sign>0
   if sign==1:stats['max_positive_phase_error']=max(stats['max_positive_phase_error'],abs(z/abs(z)-1))
  else:assert abs(ref)<2e-7 and abs(z)<2e-7
 per.append(local)
(o/'round5_summary.json').write_text(json.dumps(dict(totals=stats,files=per),indent=2)+'\n');print(json.dumps(stats,indent=2))
