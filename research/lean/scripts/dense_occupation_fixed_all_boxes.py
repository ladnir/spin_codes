"""Generate all 433 occupation rectangles; Lean alone certifies their claims."""
from pathlib import Path
from fractions import Fraction as Q
import json
r=Path(__file__).resolve().parents[1]
src=r.parent/'workstreams/inner_design/imt_asymptotic/d11'
outer=json.loads((src/'OUTER_REFINED.json').read_text())
boxes=[b for b in json.loads((src/'DENSE_REPLAY.json').read_text())['leaves'] if b['rational_witness']['family']=='occupation']
template=(r/'SpinCodes/Structured/DenseOccupationFixedFirstBox.lean').read_text(encoding='utf-8-sig')
seen={}; manifest=[]
for i,b in enumerate(boxes):
 key=json.dumps(b['rational_witness'],sort_keys=True)
 wi=seen.setdefault(key,len(seen)); witness='First' if wi==0 else f'W{wi:03d}'
 candidate='first' if wi==0 else witness
 d=json.loads((r/f'scripts/map_data/dense_occupation_fixed_{candidate}_candidate.json').read_text())
 seg=b['segment']
 if seg<19:
  o=outer['left_segments'][seg];m=Q(o['slope']);c=Q(o['intercept'])
 elif seg==19:m=Q(0);c=Q(outer['central_constant_upper'])
 else:
  o=outer['left_segments'][38-seg];m=-Q(o['slope']);c=Q(o['intercept'])+Q(o['slope'])
 tag=f'B{i:03d}';w=b['rational_witness']
 args={'m':m,'c':c,'p':w['p'],'y':w['y'],'radius':Q(d['radius_scaled'],10**30),'z':Q(d['z_scaled'],10**30),'a0':b['alpha'][0],'a1':b['alpha'][1],'x0':b['row_density'][0],'x1':b['row_density'][1]}
 lines=['import SpinCodes.Structured.DenseOccupationFixedVertexDefs','',f'namespace Spin.Structured.DenseOccupationFixed.{tag}','set_option maxHeartbeats 0','set_option maxRecDepth 100000']
 for name,value in args.items():
  q=Q(value);lines.append(f'def {name} : QInput := ⟨{q.numerator}, {q.denominator}⟩')
 lines.append(f'def upper : Int := {-4*10**23}')
 for ai in [0,1]:
  for xi in [0,1]:lines.append(f'theorem check_{ai}{xi} : vertexCheck m c p y radius z a{ai} x{xi} 50 upper = true := by decide')
 lines.append(f'end Spin.Structured.DenseOccupationFixed.{tag}')
 (r/f'SpinCodes/Structured/DenseOccupationFixed{tag}Data.lean').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 s=template.replace('DenseOccupationFixedFirstBoxData',f'DenseOccupationFixed{tag}Data').replace('DenseOccupationFixed.FirstBox',f'DenseOccupationFixed.{tag}')
 s=s.replace('DenseOccupationFixedFirst\n',f'DenseOccupationFixed{witness}\n').replace('First.',f'{witness}.').replace('first retained occupation box',f'occupation box {i}')
 (r/f'SpinCodes/Structured/DenseOccupationFixed{tag}.lean').write_text(s,encoding='utf-8')
 manifest.append({'box':tag,'source_occupation_index':i,'segment':seg,'witness':witness,'source_box':b})
(r/'scripts/map_data/dense_occupation_fixed_all_boxes_candidates.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATES','boxes':manifest},indent=2)+'\n')
print('Generated',len(boxes),'boxes associated with',len(seen),'witnesses')
