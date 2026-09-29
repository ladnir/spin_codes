from pathlib import Path
from fractions import Fraction as F
import json
R=Path(__file__).resolve().parents[1];P=R/'SpinCodes/Structured';base=R.parent/'workstreams/inner_design/imt_asymptotic'
D=json.loads((base/'FIXED_LIMIT.json').read_text());O=json.loads((base/'d11/OUTER_REFINED.json').read_text());rows=[]
for i in range(31):
 a,b=D['checks'][2*i:2*i+2];u=F(a['u']);lo=F(a['x']);hi=F(b['x'])
 if i<15:n=O['left_segments'][i];s,c=F(n['slope']),F(n['intercept'])
 elif i==15:s,c=F(0),F(O['central_constant_upper'])
 else:n=O['left_segments'][30-i];s,c=-F(n['slope']),F(n['slope'])+F(n['intercept'])
 # Old slope1/20 is the refined index16, not index14.
 if i==14:n=next(a for a in O['left_segments'] if F(a['slope'])==F(1,20));s,c=F(n['slope']),F(n['intercept'])
 if i==16:n=next(a for a in O['left_segments'] if F(a['slope'])==F(1,20));s,c=-F(n['slope']),F(n['slope'])+F(n['intercept'])
 v=max(1+u*F(127,250)*F(3,1600),u*F(1,524287)/F(3,1600)+F(127,250)*(1+u))
 rows.append(dict(lo=lo,hi=hi,s=s,c=c,u=u,v=v))
lines=['import SpinCodes.Structured.ConcreteFixedNumericDefs','','namespace Spin.Structured.ConcreteFixedNumeric','open DenseOccupationFixed','set_option maxHeartbeats 0','set_option maxRecDepth 100000']
for i,row in enumerate(rows):
 lines.append(f'namespace Segment{i:02d}')
 for name,q in row.items():lines.append(f'def {name} : QInput := ⟨{q.numerator}, {q.denominator}⟩')
 for side in ['lo','hi']:lines.append(f'theorem check_{side} : rowCheck {side} s c u v=true := by decide +kernel')
 lines.append(f'end Segment{i:02d}')
lines.append('end Spin.Structured.ConcreteFixedNumeric')
(P/'ConcreteFixedNumericData.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(R/'scripts/map_data/fixed_numeric_candidates.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATES','rows':[{k:str(v) for k,v in row.items()} for row in rows]},indent=2)+'\n')
print('generated31 segments,62 checks')
