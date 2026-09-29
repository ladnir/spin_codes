from pathlib import Path
from fractions import Fraction as Q
import json
r=Path('.')
outer=json.load(open('../workstreams/inner_design/imt_asymptotic/d11/OUTER_REFINED.json'))
template=(r/'SpinCodes/Structured/DenseOccupationFixedFirstBox.lean').read_text(encoding='utf-8-sig')
for i in range(1,9):
 tag=f'W{i:03d}';d=json.load(open(f'scripts/map_data/dense_occupation_fixed_{tag}_candidate.json'));b=d['source_box'];w=b['rational_witness'];seg=b['segment'];assert seg<19
 o=outer['left_segments'][seg]
 args={'m':o['slope'],'c':o['intercept'],'p':w['p'],'y':w['y'],'radius':str(Q(d['radius_scaled'],10**30)),'z':str(Q(d['z_scaled'],10**30)),'a0':b['alpha'][0],'a1':b['alpha'][1],'x0':b['row_density'][0],'x1':b['row_density'][1]}
 lines=['import SpinCodes.Structured.DenseOccupationFixedVertexDefs','',f'namespace Spin.Structured.DenseOccupationFixed.{tag}Box','set_option maxHeartbeats 0','set_option maxRecDepth 100000']
 for k,v in args.items():
  q=Q(v);lines.append(f'def {k} : QInput := ⟨{q.numerator}, {q.denominator}⟩')
 lines.append(f'def upper : Int := {-4*10**23}')
 for ai in [0,1]:
  for xi in [0,1]:lines.append(f'theorem check_{ai}{xi} : vertexCheck m c p y radius z a{ai} x{xi} 50 upper = true := by decide')
 lines.append(f'end Spin.Structured.DenseOccupationFixed.{tag}Box')
 (r/f'SpinCodes/Structured/DenseOccupationFixed{tag}BoxData.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 s=template.replace('DenseOccupationFixedFirstBoxData',f'DenseOccupationFixed{tag}BoxData').replace('DenseOccupationFixedFirst',f'DenseOccupationFixed{tag}').replace('DenseOccupationFixed.FirstBox',f'DenseOccupationFixed.{tag}Box').replace('First.',f'{tag}.').replace('first retained occupation box',f'retained occupation box {i}')
 (r/f'SpinCodes/Structured/DenseOccupationFixed{tag}Box.lean').write_text(s,encoding='utf-8')
 print(tag,'box files generated')
