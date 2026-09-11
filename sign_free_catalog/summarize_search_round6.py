#!/usr/bin/env python3
"""Verify all six new validators' CSV evidence; do not mutate the registry."""
import csv,json,math
from pathlib import Path
r=Path(__file__).resolve().parent;o=r/'search_results';files=sorted(o.glob('round6_*.csv'));assert len(files)==12,(len(files),files)
tot=dict(rows=0,positive=0,zero=0,negative_controls=0,complex_controls=0,max_relative_error=0.,max_positive_phase_error=0.);per=[]
for p in files:
 local=dict(file=p.name,rows=0,positive=0,zero=0,negative_controls=0,complex_controls=0)
 for row in csv.DictReader(p.open()):
  assert None not in row,'Malformed CSV columns'
  z=complex(float(row.get('weight_real',row.get('real_pf'))),float(row.get('weight_imag',row.get('imag_pf'))))
  ref=complex(float(row.get('reference_real',row.get('real_reference'))),float(row.get('reference_imag',row.get('imag_reference'))))
  sign=int(row['expected']);error=float(row.get('relative_error',row.get('scaled_error')))
  assert all(math.isfinite(x) for x in (z.real,z.imag,ref.real,ref.imag,error)) and error<2e-7
  # Independent universal scaled comparison, in addition to validator-specific errors.
  assert abs(z-ref)/max(1,abs(ref))<2e-7
  key={1:'positive',0:'zero',-1:'negative_controls',2:'complex_controls'}[sign]
  tot['rows']+=1;local['rows']+=1;tot[key]+=1;local[key]+=1;tot['max_relative_error']=max(tot['max_relative_error'],error)
  if sign in (1,-1):
   assert abs(z)>0 and ref.real*sign>0 and abs(z/abs(z)-sign)<2e-7
   if sign==1:tot['max_positive_phase_error']=max(tot['max_positive_phase_error'],abs(z/abs(z)-1))
  elif sign==0:assert abs(z)<2e-7 and abs(ref)<2e-7
  else:assert abs(ref.imag)>0 and abs(z/abs(z)-ref/abs(ref))<2e-7
 per.append(local)
(o/'round6_summary.json').write_text(json.dumps(dict(totals=tot,files=per),indent=2)+'\n');print(json.dumps(tot,indent=2))
