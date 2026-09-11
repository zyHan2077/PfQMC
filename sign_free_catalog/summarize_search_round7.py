#!/usr/bin/env python3
"""Independently summarize the final PXP evidence without changing the registry."""
import csv,json,math
from pathlib import Path
o=Path(__file__).resolve().parent/'search_results';files=sorted(o.glob('round7_*.csv'));assert len(files)==2
s=dict(rows=0,positive=0,zero=0,negative_controls=0,complex_controls=0,max_relative_error=0.,max_positive_phase_error=0.);per=[]
for p in files:
 local=dict(file=p.name,rows=0,positive=0,zero=0,negative_controls=0)
 for row in csv.DictReader(p.open()):
  assert None not in row
  z=complex(float(row['real_pf']),float(row['imag_pf']));ref=complex(float(row['real_reference']),float(row['imag_reference']));e=float(row['relative_error']);sign=int(row['expected'])
  assert all(math.isfinite(x) for x in [z.real,z.imag,ref.real,ref.imag,e]) and e<2e-7 and abs(z-ref)/max(1,abs(ref))<2e-7
  key={1:'positive',0:'zero',-1:'negative_controls'}[sign];s[key]+=1;local[key]+=1;s['rows']+=1;local['rows']+=1;s['max_relative_error']=max(s['max_relative_error'],e)
  if sign:
   assert abs(z)>0 and sign*ref.real>0 and abs(z/abs(z)-sign)<2e-7
   if sign==1:s['max_positive_phase_error']=max(s['max_positive_phase_error'],abs(z/abs(z)-1))
  else:assert abs(z)<2e-7 and ref==0
 per.append(local)
(o/'round7_summary.json').write_text(json.dumps(dict(totals=s,files=per),indent=2)+'\n');print(json.dumps(s,indent=2))
