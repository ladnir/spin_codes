import re
from fractions import Fraction as F
s=open('SpinCodes/Majorant/RefinedData.lean',encoding='utf-8').read().split('def refined')[0]
pairs=[tuple(F(t.replace(' ','')) for t in p) for p in re.findall(r'\(\(([-\d /]+)\), \(([-\d /]+)\)\)',s)]
xs=[F(13,125),F(3,20),F(1,5),F(1,4),F(3,10),F(7,20),F(2,5),F(9,20),F(1,2)]
ix=[0,5,7,8,9,10,11,13]
def frac(x):return f'({x.numerator}:ℝ)/{x.denominator}'
def fix(x):return f'(ofFrac ({x.numerator}) {x.denominator})'
out=['import SpinCodes.Structured.ConcreteOuterDefectNumeric','import SpinCodes.Majorant.RefinedData','','set_option maxRecDepth 100000','set_option maxHeartbeats 0','noncomputable section','namespace Spin.Structured.ConcreteOuter','open Spin.Numeric Spin.Numeric.Fix','']
for j,i in enumerate(ix):
 a,b=pairs[i]
 for side,x in [('lo',xs[j]),('hi',xs[j+1])]:
  out +=[f'lemma defect_endpoint_{j}_{side} :',f'    ({frac(a)})*({frac(x)})+({frac(b)})-hEnt ({frac(x)})+Real.log 2/2 ≤ 1281/100000 := by',f'  have hc : defectOK {fix(x)} {fix(a)} {fix(b)} = true := by decide +kernel',f'  have hh := defectOK_sound (ofFrac_mem (p:={x.numerator}) (q:={x.denominator}) (by norm_num))',f'    (ofFrac_mem (p:={a.numerator}) (q:={a.denominator}) (by norm_num))',f'    (ofFrac_mem (p:={b.numerator}) (q:={b.denominator}) (by norm_num)) hc','  norm_num at hh ⊢','  exact hh','']
 out +=[f'lemma defect_interval_{j} {{x : ℝ}} (hx : x ∈ Set.Icc ({frac(xs[j])}) ({frac(xs[j+1])})) :', '    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by',f'  have hl := defect_endpoint_{j}_lo',f'  have hh := defect_endpoint_{j}_hi',f'  have hbound := entropy_defect_interval (d := -Real.log 2/2+1281/100000) (s:=({frac(a)})) (i:=({frac(b)}))', '    (by norm_num) (by norm_num) hx (by linarith) (by linarith)',f'  have hm := Spin.Majorant.refined.toFun_le_support Spin.Majorant.member{i} x','  norm_num at hm','  linarith','']
out+=['theorem refined_defect_left {x : ℝ} (hx : x ∈ Set.Icc ((13:ℝ)/125) (1/2)) :','    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by']
for j in range(7):
 out +=[f'  by_cases h{j} : x ≤ ({frac(xs[j+1])})',f'  · exact defect_interval_{j} ⟨by linarith [hx.1], h{j}⟩']
out+=['  exact defect_interval_7 ⟨by linarith, hx.2⟩','','theorem refined_defect {x : ℝ} (hx : x ∈ Set.Icc ((13:ℝ)/125) (112/125)) :','    Spin.Majorant.refined.toFun x-hEnt x ≤ -Real.log 2/2+1281/100000 := by','  by_cases h : x ≤ (1:ℝ)/2','  · exact refined_defect_left ⟨hx.1,h⟩','  · have hh := refined_defect_left (x:=1-x) ⟨by linarith [hx.2], by linarith⟩','    rw [Spin.Majorant.refined_reflect] at hh','    have he : hEnt (1-x) = hEnt x := by rw [hEnt_eq_binEntropy, Real.binEntropy_one_sub, ←hEnt_eq_binEntropy]','    rwa [he] at hh','','end Spin.Structured.ConcreteOuter','']
open('SpinCodes/Structured/ConcreteOuterDefect.lean','w',encoding='utf-8').write('\n'.join(out))

