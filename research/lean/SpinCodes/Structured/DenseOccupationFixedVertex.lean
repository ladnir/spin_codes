import SpinCodes.Structured.DenseOccupationFixedVertexDefs
import SpinCodes.Structured.DenseOccupationFixedLog

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric

def QInput.real (a : QInput) : ℝ := (a.num:ℝ)/a.den

lemma QInput.interval_mem (a : QInput) (h : 0 < a.den) : Fix.Mem a.interval a.real :=
  Fix.ofFrac_mem h

theorem vertexEval_mem (m c p y radius z α x : QInput)
    (hm : 0 < m.den) (hc : 0 < c.den) (hα : 0 < α.den) (hx : 0 < x.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den) (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    {n : Nat} {res : Fix} (h : vertexEval m c p y radius z α x n = some res) :
    Fix.Mem res (boxExponent m.real c.real p.real y.real radius.real z.real α.real x.real) := by
  unfold vertexEval at h
  cases ha : klEval α.num α.den p.num p.den n with
  | none => simp [ha] at h
  | some a =>
    cases hb : klEval x.num x.den y.num y.den n with
    | none => simp [ha,hb] at h
    | some b =>
      cases hr : Fix.flogQ radius.num radius.den n with
      | none => simp [ha,hb,hr] at h
      | some lr =>
        cases hz : Fix.flogQ z.num z.den n with
        | none => simp [ha,hb,hr,hz] at h
        | some lz =>
          simp only [ha,hb,hr,hz,Option.bind_eq_bind,Option.bind_some,Option.pure_def,Option.some.injEq] at h
          subst res
          exact Fix.sub_mem
            (Fix.add_mem (Fix.add_mem (Fix.add_mem
              (Fix.mul_mem (α.interval_mem hα) (Fix.add_mem
                (Fix.mul_mem (m.interval_mem hm) (x.interval_mem hx)) (c.interval_mem hc)))
              (klEval_mem hα hp0 hp1 ha))
              (Fix.mul_mem (α.interval_mem hα) (klEval_mem hx hy0 hy1 hb)))
              (Fix.divInt_mem (k := 128) (Fix.flogQ_mem hr) (by norm_num)))
            (Fix.mul_mem (by simpa using Fix.ofFrac_mem (p := 11) (q := 100) (by norm_num))
              (Fix.flogQ_mem hz))

/-- A successful integer vertex verdict implies the actual analytic exponent bound. -/
theorem vertexCheck_sound (m c p y radius z α x : QInput)
    (hm : 0 < m.den) (hc : 0 < c.den) (hα : 0 < α.den) (hx : 0 < x.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den) (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    {n : Nat} {upper : Int} (h : vertexCheck m c p y radius z α x n upper = true) :
    boxExponent m.real c.real p.real y.real radius.real z.real α.real x.real ≤
      (upper:ℝ)/scale := by
  unfold vertexCheck at h
  cases he : vertexEval m c p y radius z α x n with
  | none => simp [he] at h
  | some res =>
    have hi : res.hi ≤ upper := by simpa [he] using h
    exact (vertexEval_mem m c p y radius z α x hm hc hα hx hp0 hp1 hy0 hy1 he).2.trans
      (div_le_div_of_nonneg_right (by exact_mod_cast hi) scaleR_pos.le)

#print axioms vertexCheck_sound
end Spin.Structured.DenseOccupationFixed
