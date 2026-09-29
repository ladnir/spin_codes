from pathlib import Path
from fractions import Fraction as F
import json
R=Path(__file__).resolve().parents[1];P=R/'SpinCodes/Structured'
rows=json.loads((R/'scripts/map_data/fixed_numeric_candidates.json').read_text())['rows']
L=['import SpinCodes.Structured.ConcreteFixedNumericData','import SpinCodes.Structured.ConcreteFixedNumericInterval','','noncomputable section','namespace Spin.Structured.ConcreteFixedNumeric','open Spin.Numeric DenseOccupationFixed Set','set_option maxHeartbeats 0','set_option maxRecDepth 100000']
for j,row in enumerate(rows):
 L.append(f'namespace Segment{j:02d}')
 for side in ['lo','hi']:
  L += [f'theorem bound_{side} : rowRate {side}.real s.real c.real u.real v.real ≤ -(8679/10000000) :=',f'  rowCheck_sound (by decide) (by decide) (by decide) check_{side}']
 L += ['def supportLine : ℚ×ℚ := ((s.num:ℚ)/s.den,(c.num:ℚ)/c.den)',
 'theorem support_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel',
 'theorem v_eq : v.real=Spin.mv (127/250) (3/1600) (1/524287) u.real := by norm_num [v,u,QInput.real,Spin.mv]',
 'theorem u_pos : 0<u.real := by norm_num [u,QInput.real]',
 'theorem interval_bound {x : ℝ} (hx : x∈Icc lo.real hi.real) : rate x u.real≤-(8679/10000000) := by',
 '  have hh := rowRate_interval (by norm_num [lo,QInput.real]) (by norm_num [hi,QInput.real]) hx bound_lo bound_hi',
 '  have hm := rate_le_rowRate support_mem x u.real v.real v_eq',
 '  have hm2 : rate x u.real ≤ rowRate x s.real c.real u.real v.real := by',
 '    simpa only [supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,Rat.cast_natCast] using hm',
 '  exact hm2.trans hh',
 f'#print axioms interval_bound',f'end Segment{j:02d}']
L.append('end Spin.Structured.ConcreteFixedNumeric')
(P/'ConcreteFixedNumericSegments.lean').write_text('\n'.join(L)+'\n',encoding='utf-8')
