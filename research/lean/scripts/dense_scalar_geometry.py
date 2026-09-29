"""Generate scalar box-to-global-geometry proofs; Lean checks every association."""
from fractions import Fraction as Q
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent/'workstreams/inner_design/imt_asymptotic/d11'

def rat(q):
    q = Q(q)
    return f'({q.numerator}/{q.denominator})'

def generate(index):
    leaves = json.loads((SOURCE/'DENSE_REPLAY.json').read_text())['leaves']
    entries = [(i,b) for i,b in enumerate(leaves) if b['rational_witness']['family']=='statefree']
    global_index, box = entries[index]
    outer = json.loads((SOURCE/'OUTER_REFINED.json').read_text())
    seg = box['segment']
    if seg < 19:
        a = outer['left_segments'][seg]; m,c = Q(a['slope']),Q(a['intercept'])
    elif seg == 19:
        m,c = Q(0),Q(outer['central_constant_upper'])
    else:
        a = outer['left_segments'][38-seg]; m,c = -Q(a['slope']),Q(a['slope'])+Q(a['intercept'])
    tag = f'B{index:03d}'
    rect = ','.join(rat(q) for q in box['alpha']+box['row_density'])
    text = f'''import SpinCodes.Structured.DenseScalarExact{tag}
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.{tag}
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox {global_index} (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD {global_index} zeroRect = ⟨{rect}⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨({rat(m)},{rat(c)}), Spin.Majorant.member{seg}, ?_⟩
  have hm : (({rat(m)}:ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (({rat(c)}:ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (({rat(m)}:ℚ):ℝ) (({rat(c)}:ℚ):ℝ) (4/10000000) 1
  rw [hm,hc]
  exact pointRate
    (lt_of_lt_of_le (by norm_num [a0,QInput.real]) hα.1)
    (le_trans hα.2 (by norm_num [a1,QInput.real]))
    (lt_of_lt_of_le (by norm_num [x0,QInput.real]) hx.1)
    (le_trans hx.2 (by norm_num [x1,QInput.real]))
    (by norm_num [p,QInput.real]) (by norm_num [p,QInput.real])
    (by norm_num [y,QInput.real]) (by norm_num [y,QInput.real])
    (by norm_num [radius,QInput.real]) (by norm_num [z,QInput.real])
    (by norm_num [z,QInput.real]) scalar_bound (exponent_bound hα hx)

#print axioms certified
end Spin.Structured.DenseScalarExact.{tag}
'''
    path = ROOT/f'SpinCodes/Structured/DenseScalarGeometry{tag}.lean'
    path.write_text(text,encoding='utf-8')
    print(tag,'global index',global_index,flush=True)
    return path

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('index',type=int)
    generate(parser.parse_args().index)
