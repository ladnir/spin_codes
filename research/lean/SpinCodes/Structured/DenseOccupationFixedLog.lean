import SpinCodes.Structured.DenseOccupationFixedLogDefs
import SpinCodes.Numeric.BAEval
import SpinCodes.Structured.DenseOccupationFixedVertices
import SpinCodes.Structured.DenseOccupationFixedBasic

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric

lemma klTermEval_mem {un ud vn vd : Int} (hud : 0<ud) (hvn : 0<vn) (hvd : 0<vd)
    {n : Nat} {res : Fix} (h : klTermEval un ud vn vd n = some res) :
    Fix.Mem res (((un:ℝ)/ud)*Real.log (((un:ℝ)/ud)/((vn:ℝ)/vd))) := by
  unfold klTermEval at h
  split_ifs at h with hu
  · cases Option.some.inj h
    subst un
    simpa using zero_mem
  · cases he : Fix.flogQ (un*vd) (ud*vn) n with
    | none => simp [he] at h
    | some lg =>
      simp only [he, Option.bind_eq_bind,Option.bind_some, Option.pure_def, Option.some.injEq] at h
      subst res
      have hl := Fix.flogQ_mem he
      have hh := Fix.mul_mem (Fix.ofFrac_mem (p := un) hud) hl
      convert hh using 1
      push_cast
      congr 2
      have hd : (ud:ℝ) ≠ 0 := by exact_mod_cast hud.ne'
      have hv : (vn:ℝ) ≠ 0 := by exact_mod_cast hvn.ne'
      have hvd' : (vd:ℝ) ≠ 0 := by exact_mod_cast hvd.ne'
      field_simp

/-- Certified interval evaluation of binary KL, including both first-argument endpoints. -/
theorem klEval_mem {un ud vn vd : Int} (hud : 0<ud) (hvn : 0<vn) (hvnd : vn<vd)
    {n : Nat} {res : Fix} (h : klEval un ud vn vd n = some res) :
    Fix.Mem res (binKL ((un:ℝ)/ud) ((vn:ℝ)/vd)) := by
  have hvd : 0<vd := lt_trans hvn hvnd
  unfold klEval at h
  cases ha : klTermEval un ud vn vd n with
  | none => simp [ha] at h
  | some a =>
    cases hb : klTermEval (ud-un) ud (vd-vn) vd n with
    | none => simp [ha,hb] at h
    | some b =>
      simp only [ha,hb,Option.bind_eq_bind,Option.bind_some,Option.pure_def,Option.some.injEq] at h
      subst res
      have hh := Fix.add_mem (klTermEval_mem hud hvn hvd ha)
        (klTermEval_mem hud (sub_pos.mpr hvnd) hvd hb)
      have hd : (ud:ℝ) ≠ 0 := by exact_mod_cast hud.ne'
      have hvd' : (vd:ℝ) ≠ 0 := by exact_mod_cast hvd.ne'
      have he : ((ud-un:ℤ):ℝ)/(ud:ℝ) = 1-(un:ℝ)/ud := by push_cast; field_simp
      have hf : ((vd-vn:ℤ):ℝ)/(vd:ℝ) = 1-(vn:ℝ)/vd := by push_cast; field_simp
      simpa only [he,hf,binKL] using hh

#print axioms klEval_mem
end Spin.Structured.DenseOccupationFixed

