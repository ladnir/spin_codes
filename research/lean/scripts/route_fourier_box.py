"""Propose one Fourier exponent box and connect its original global rectangle."""
from pathlib import Path
from fractions import Fraction as Q
import json,argparse,os,sys
R=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('index',type=int);args=p.parse_args()
# These guards stop existing dispatchers at exact ownership boundaries.
lane=os.environ.get('FOURIER_BOX_LANE','')
lanes=json.loads((R/'scripts/map_data/dense_fourier_box_lanes.json').read_text())
stop=None
if 115<=args.index<154 and lane!='prefix_tail':
 stop=('dense_fourier_box_replay.json',115,list(range(115)))
if args.index in lanes['suffix_tail']['indices'] and lane!='suffix_tail':
 stop=('dense_fourier_box_second_replay.json',105,None)
if stop:
 path=R/'scripts/map_data'/stop[0];state=json.loads(path.read_text());assert state['checked']==stop[1]
 state.update(status='PASS',total=stop[1],ownership_stop=True)
 state['owned_indices']=[row['index'] for row in state['boxes']]
 path.write_text(json.dumps(state,indent=2)+'\n')
 print('LANE_COMPLETE: all owned indices checked; exit75 is ownership stop',flush=True);sys.exit(75)
if args.index>=154 and os.environ.get('FOURIER_SECOND_HALF')!='1':
 raise RuntimeError('Original prefix dispatcher crossed its ownership range')
S=R.parent/'workstreams/inner_design/imt_asymptotic/d11'
leaves=json.loads((S/'DENSE_REPLAY.json').read_text())['leaves']
selected=[(i,b) for i,b in enumerate(leaves) if b['rational_witness']['family']=='fourier']
unique={}
for gi,b in selected:
 key=json.dumps(b['rational_witness'],sort_keys=True)
 if key not in unique:unique[key]=len(unique)
gi,box=selected[args.index];wi=unique[json.dumps(box['rational_witness'],sort_keys=True)];wt=f'W{wi:03d}';tag=f'B{args.index:03d}'
cand=json.loads((R/f'scripts/map_data/dense_fourier_exact_{wt}_candidate.json').read_text());w=box['rational_witness']
outer=json.loads((S/'OUTER_REFINED.json').read_text());seg=box['segment']
if seg<19:a=outer['left_segments'][seg];m,c=Q(a['slope']),Q(a['intercept'])
elif seg==19:m,c=Q(0),Q(outer['central_constant_upper'])
else:a=outer['left_segments'][38-seg];m,c=-Q(a['slope']),Q(a['slope'])+Q(a['intercept'])
values=dict(m=m,c=c,p=Q(w['p']),y=Q(w['y']),radius=Q(cand['radius']),z=Q(cand['z']),a0=Q(box['alpha'][0]),a1=Q(box['alpha'][1]),x0=Q(box['row_density'][0]),x1=Q(box['row_density'][1]))
data=['import SpinCodes.Structured.DenseOccupationFixedVertexDefs','',f'namespace Spin.Structured.DenseFourierExact.{tag}','open Spin.Structured.DenseOccupationFixed','set_option maxHeartbeats 0','set_option maxRecDepth 100000']
for name,q in values.items():data.append(f'def {name} : QInput := ⟨{q.numerator}, {q.denominator}⟩')
data.append('def upper : Int := -400000000000000000000000')
for a in range(2):
 for x in range(2):data.append(f'theorem check_{a}{x} : vertexCheck m c p y radius z a{a} x{x} 50 upper=true := by decide')
data.append(f'end Spin.Structured.DenseFourierExact.{tag}')
(R/f'SpinCodes/Structured/DenseFourierExact{tag}Data.lean').write_text('\n'.join(data)+'\n',encoding='utf-8')
template=(R/'SpinCodes/Structured/DenseOccupationFixedFirstBox.lean').read_text(encoding='utf-8-sig')
body=template[template.index('lemma vertex00'):template.index('theorem parameters_match')].replace('first retained occupation box',f'Fourier box {args.index}')
sem=[f'import SpinCodes.Structured.DenseFourierExact{tag}Data',f'import SpinCodes.Structured.DenseFourierExact{wt}','import SpinCodes.Structured.DenseOccupationFixedVertex','import SpinCodes.Structured.DenseFourierExactRate','import SpinCodes.Structured.DenseGeometryRate','', 'noncomputable section', f'namespace Spin.Structured.DenseFourierExact.{tag}', 'open Spin.Numeric Spin.Imt DenseOccupationFixed DenseGeometry Set',body,
 f'theorem parameters_match : p.real*y.real = ({wt}.qn:ℝ)/{wt}.qd ∧ radius.real = ({wt}.rn:ℝ)/{wt}.rd ∧ z.real = ({wt}.zn:ℝ)/{wt}.zd := by',f'  norm_num [p,y,radius,z,QInput.real,{wt}.qn,{wt}.qd,{wt}.rn,{wt}.rd,{wt}.zn,{wt}.zd]',
 f'theorem collatz : ((ConcreteFourier.matrix (p.real*y.real) z.real).applyCol {wt}.v.real).le (Coords.smul radius.real {wt}.v.real) := by', '  rw [parameters_match.1,parameters_match.2.1,parameters_match.2.2]',f'  exact {wt}.collatz',
 'def supportLine : ℚ×ℚ := ((m.num:ℚ)/m.den,(c.num:ℚ)/c.den)', 'lemma supportLine_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel',
 f'lemma geometry_match : boxes.getD {gi} zeroRect = ⟨(a0.num:ℚ)/a0.den,(a1.num:ℚ)/a1.den,(x0.num:ℚ)/x0.den,(x1.num:ℚ)/x1.den⟩ := by decide +kernel',
 f'theorem certified : CertifiedBox {gi} (4/10000000) {cand["prefactor"]} := by',
 '  intro α x hg', '  rw [geometry_match] at hg','  have hα : α∈Icc a0.real a1.real := by','    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢','    exact ⟨hg.1,hg.2.1⟩','  have hx : x∈Icc x0.real x1.real := by','    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢','    exact ⟨hg.2.2.1,hg.2.2.2⟩',f'  have hf := {wt}.witness_floor',
 f'  have hr : PointRate α x m.real c.real (4/10000000) ({wt}.v.real.Z/{wt}.wmin) := fourier_pointRate',
 '    (by norm_num [a0,QInput.real] at hα; linarith [hα.1] : 0<α)',
 '    (by norm_num [a1,QInput.real] at hα; linarith [hα.2] : α≤1)',
 '    (by norm_num [x0,QInput.real] at hx; linarith [hx.1] : 0<x)',
 '    (by norm_num [x1,QInput.real] at hx; linarith [hx.2] : x≤1)',
 '    (by norm_num [p,QInput.real] : 0<p.real) (by norm_num [p,QInput.real] : p.real<1)',
 '    (by norm_num [y,QInput.real] : 0<y.real) (by norm_num [y,QInput.real] : y.real<1)',
 '    (by norm_num [radius,QInput.real] : 0<radius.real)',
 '    (by norm_num [z,QInput.real] : 0<z.real) (by norm_num [z,QInput.real] : z.real≤1)',
 '    hf.1 hf.2.1 hf.2.2.1 hf.2.2.2 collatz (exponent_bound hα hx)',
 f'  have hrFinal : PointRate α x m.real c.real (4/10000000) {wt}.prefactor := hr.mono_constant {wt}.witness_prefactor',
 '  refine ⟨supportLine,supportLine_mem,?_⟩',
 f'  simpa only [PointRate,supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,{wt}.prefactor] using @hrFinal',
 f'theorem certified_uniform : CertifiedBox {gi} (4/10000000) 24000000000000 := certified.mono_constant (by norm_num)', '#print axioms exponent_bound','#print axioms collatz','#print axioms certified','#print axioms certified_uniform',f'end Spin.Structured.DenseFourierExact.{tag}']
# Apostrophes in Python adjacent literals above are intentionally avoided in emitted identifiers.
sem='\n'.join(sem).replace('have hr := hr.mono_constant','have hrFinal := hr.mono_constant').replace('using hr\n#print','using @hrFinal\n#print')+'\n'
(R/f'SpinCodes/Structured/DenseFourierExact{tag}.lean').write_text(sem,encoding='utf-8')
(R/f'scripts/map_data/dense_fourier_box_{tag}_candidate.json').write_text(json.dumps({'status':'UNTRUSTED_CANDIDATE','index':args.index,'global_index':gi,'witness_index':wi,'prefactor':cand['prefactor'],'source_box':box},indent=2)+'\n')
print(tag,'global',gi,'witness',wt,'C',cand['prefactor'])




