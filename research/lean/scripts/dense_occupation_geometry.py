"""Propose a guillotine cover of all 1023 original dense boxes. Lean checks every leaf."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
r=Path(__file__).resolve().parents[1]
source=r.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json'
d=json.loads(source.read_text())
boxes=[tuple(map(F,b['alpha']+b['row_density'])) for b in d['leaves']]
def q(x):return f'({x.numerator}/{x.denominator})'
nodes=leaves=depth=0
def split(region,ids,level=0):
 global nodes,leaves,depth
 nodes+=1;depth=max(depth,level)
 for i in ids:
  b=boxes[i]
  if b[0]<=region[0] and region[1]<=b[1] and b[2]<=region[2] and region[3]<=b[3]:
   leaves+=1;return f'(.leaf {i})'
 choices=[]
 for a in [0,2]:
  for v in sorted(set(boxes[i][a+1] for i in ids)):
   if not region[a]<v<region[a+1]:continue
   left=[i for i in ids if boxes[i][a]<v]
   right=[i for i in ids if v<boxes[i][a+1]]
   choices.append((len(left)+len(right),max(len(left),len(right)),a,v,left,right))
 _,_,a,v,left,right=min(choices)
 lo=list(region);lo[a+1]=v;hi=list(region);hi[a]=v
 kind='alpha' if a==0 else 'density'
 return f'(.{kind} {q(v)} {split(tuple(lo),left,level+1)} {split(tuple(hi),right,level+1)})'
root=(F('1/10000'),F(1),F('13/125'),F('112/125'))
tree=split(root,list(range(len(boxes))))
rect=lambda b:'⟨'+','.join(map(q,b))+'⟩'
text='import SpinCodes.Structured.DenseOccupationGeometryDefs\nnamespace Spin.Structured.DenseGeometry\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n'
text+='def boxes : List Rect := [\n'+',\n'.join(map(rect,boxes))+']\n'
text+=f'def whole : Rect := {rect(root)}\ndef tree : Tree := {tree}\n'
text+='theorem checked : check boxes tree whole = true := by decide +kernel\n'
text+='theorem boxes_length : boxes.length = 1023 := by decide +kernel\nend Spin.Structured.DenseGeometry\n'
(r/'SpinCodes/Structured/DenseOccupationGeometryData.lean').write_text(text,encoding='utf-8')
manifest={'status':'UNTRUSTED_CANDIDATE','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'nodes':nodes,'leaves':leaves,'depth':depth,'index_mapping':[{'index':i,'family':b['rational_witness']['family'],'segment':b['segment'],'alpha':b['alpha'],'row_density':b['row_density']} for i,b in enumerate(d['leaves'])]}
(r/'scripts/map_data/dense_occupation_geometry_candidate.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(nodes,'nodes;',leaves,'leaves; depth',depth)
