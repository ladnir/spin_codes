import pathlib,re
src=pathlib.Path('SpinCodes/Numeric/BAEval.lean').read_text(encoding='utf8')
a=src.index('lemma infeasible_vacuous '); z=src.index('/-! ## The direct instances',a)
body=src[a:z]
body=body.replace('DenseClaim','ClosedDenseClaim').replace('infeasible_vacuous','infeasible_closed_vacuous').replace('leafOK_sound','leafOK_closed_sound').replace('denseTail_of_checkTree','denseTail_closed_of_checkTree')
body=body.replace('have hb0 : 0 < b','have hb0 : 0 ≤ b')
body=body.replace('  have hb1 : b < 1 := by linarith\n','')
body=body.replace('fpiBestClampL_mem','fpiBestClampL_closed_mem')
body=body.replace('gBA_le ha0 ha1.le huR','gBA_le ha0 ha1 huR')
body=body.replace('piBA_eq_piEval ha0 ha1 f1 f2','piBA_eq_piEval_closed ha0 f1 f2')
body=body.replace('piBA_eq_piEval hb0.le hb1 f3 f4','piBA_eq_piEval_closed hb0 f3 f4')
head='''import SpinCodes.Numeric.MajorantEval

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Spin.Numeric.Fix

def ClosedDenseClaim (thr : Int) (a b w : ℝ) : Prop :=
  (a/2≤b ∧ b≤1-a/2 ∧ b/2≤w ∧ w≤1-b/2 ∧ 0≤a ∧ a≤1) →
    gBA a+piBA a b+piBA b w < (thr:ℝ)/(scale:ℝ)

'''
pathlib.Path('SpinCodes/Structured/ConcreteOuterTailDenseSound.lean').write_text(head+body+'\nend Spin.Structured.ConcreteOuter\n',encoding='utf8')
