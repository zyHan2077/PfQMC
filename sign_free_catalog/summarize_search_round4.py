#!/usr/bin/env python3
"""Summarize round four; preserve old evidence and the independent registry."""
import csv,json,math
from pathlib import Path
r=Path(__file__).resolve().parent;o=r/'search_results';files=sorted(o.glob('round4_*.csv'));assert len(files)==4
stats=dict(rows=0,positive=0,zero=0,negative_controls=0,max_relative_error=0.,max_positive_phase_error=0.);per=[]
for p in files:
 local=dict(file=p.name,rows=0,positive=0,zero=0,negative_controls=0)
 for row in csv.DictReader(p.open()):
  sign=int(row['expected']);z=complex(float(row['weight_real']),float(row['weight_imag']));error=float(row.get('relative_error',row.get('scaled_error')))
  assert math.isfinite(error) and error<2e-7
  key={1:'positive',0:'zero',-1:'negative_controls'}[sign]
  stats['rows']+=1;local['rows']+=1;stats[key]+=1;local[key]+=1;stats['max_relative_error']=max(stats['max_relative_error'],error)
  if sign==1:
   assert abs(z)>0 and abs(z/abs(z)-1)<2e-7
   stats['max_positive_phase_error']=max(stats['max_positive_phase_error'],abs(z/abs(z)-1))
 per.append(local)
(o/'round4_summary.json').write_text(json.dumps(dict(totals=stats,files=per),indent=2)+'\n');print(json.dumps(stats,indent=2))
