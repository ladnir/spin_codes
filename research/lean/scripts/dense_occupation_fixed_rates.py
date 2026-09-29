"""Link each occupation certificate to its exact original global rectangle."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,re
parser=argparse.ArgumentParser();parser.add_argument("--only",nargs="*");args=parser.parse_args()
r=Path(__file__).resolve().parents[1]
d=json.loads((r.parent/'workstreams/inner_design/imt_asymptotic/d11/DENSE_REPLAY.json').read_text())
local=json.loads((r/'scripts/map_data/dense_occupation_fixed_all_boxes_candidates.json').read_text())['boxes']
global_ids=[i for i,b in enumerate(d['leaves']) if b['rational_witness']['family']=='occupation']
template=(r/'SpinCodes/Structured/DenseOccupationFixedB000Rate.lean').read_text(encoding='utf-8')
manifest=[]
for entry,gi in zip(local,global_ids):
 tag=entry['box'];wi=entry['witness'];candidate='first' if wi=='First' else wi
 c=json.loads((r/f'scripts/map_data/dense_occupation_fixed_{candidate}_candidate.json').read_text())
 pref=c.get('prefactor',372);assert pref<=7000000
 s=re.sub(r'\b372\b',str(pref),template)
 s=s.replace('DenseOccupationFixedB000',f'DenseOccupationFixed{tag}').replace('DenseOccupationFixed.B000',f'DenseOccupationFixed.{tag}').replace('First.',f'{wi}.')
 s=s.replace('open Spin.Numeric DenseGeometry','open Spin.Numeric DenseGeometry\nset_option maxRecDepth 100000')
 s=s.replace('boxes.getD 89','boxes.getD '+str(gi)).replace('CertifiedBox 89','CertifiedBox '+str(gi))
 q=lambda x:f'{Q(x).numerator}/{Q(x).denominator}'
 bounds=d['leaves'][gi]['alpha']+d['leaves'][gi]['row_density']
 s=re.sub(r'change \(⟨.*?⟩ : Rect\) = _', 'change (⟨'+','.join(map(q,bounds))+'⟩ : Rect) = _',s)
 if 'theorem certified_uniform' not in s:
  s=s.replace('#print axioms certified',f'theorem certified_uniform : CertifiedBox {gi} (4/10000000) 7000000 :=\n  certified.mono_constant (by norm_num)\n\n#print axioms certified\n#print axioms certified_uniform')
 assert f'import SpinCodes.Structured.DenseOccupationFixed{tag}\n' in s
 assert f'namespace Spin.Structured.DenseOccupationFixed.{tag}\n' in s
 if args.only is None or tag in args.only:
  (r/f'SpinCodes/Structured/DenseOccupationFixed{tag}Rate.lean').write_text(s,encoding='utf-8')
 manifest.append({'box':tag,'global_index':gi,'witness':wi,'prefactor':pref})
(r/'scripts/map_data/dense_occupation_fixed_rates_candidates.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATES','entries':manifest},indent=2)+'\n')
print('Generated',sum(args.only is None or e['box'] in args.only for e in manifest),'indexed occupation rate wrappers; validated',len(manifest),'candidate identities')
