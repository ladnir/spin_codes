import SpinCodes.UniformMixer.ScalableDenseExponent
import SpinCodes.Structured.DenseOccupationFixedVertex

/-! Kernel-checked rational interval evaluation of a 64-bit dense box vertex. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

def scalableVertexEval (m c p y ρ z α x : QInput) (n : ℕ) : Option Fix := do
  let a ← klEval α.num α.den p.num p.den n
  let b ← klEval x.num x.den y.num y.den n
  let lr ← Fix.flogQ ρ.num ρ.den n
  let lz ← Fix.flogQ z.num z.den n
  pure (Fix.sub
    (Fix.add (Fix.add (Fix.add
      (Fix.mul α.interval
        (Fix.add (Fix.mul m.interval x.interval) c.interval)) a)
      (Fix.mul α.interval b)) (Fix.divInt lr 64))
    (Fix.mul (Fix.ofFrac 11 100) lz))

def scalableVertexCheck (m c p y ρ z α x : QInput)
    (n : ℕ) (upper : Int) : Bool :=
  match scalableVertexEval m c p y ρ z α x n with
  | none => false
  | some v => decide (v.hi ≤ upper)

theorem scalableVertexEval_mem (m c p y ρ z α x : QInput)
    (hm : 0 < m.den) (hc : 0 < c.den)
    (hα : 0 < α.den) (hx : 0 < x.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den)
    (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    {n : ℕ} {res : Fix}
    (h : scalableVertexEval m c p y ρ z α x n = some res) :
    Fix.Mem res (scalableBoxExponent m.real c.real p.real y.real
      ρ.real z.real α.real x.real) := by
  unfold scalableVertexEval at h
  cases ha : klEval α.num α.den p.num p.den n with
  | none => simp [ha] at h
  | some a =>
    cases hb : klEval x.num x.den y.num y.den n with
    | none => simp [ha, hb] at h
    | some b =>
      cases hr : Fix.flogQ ρ.num ρ.den n with
      | none => simp [ha, hb, hr] at h
      | some lr =>
        cases hz : Fix.flogQ z.num z.den n with
        | none => simp [ha, hb, hr, hz] at h
        | some lz =>
          simp only [ha, hb, hr, hz, Option.bind_eq_bind,
            Option.bind_some, Option.pure_def, Option.some.injEq] at h
          subst res
          exact Fix.sub_mem
            (Fix.add_mem (Fix.add_mem (Fix.add_mem
              (Fix.mul_mem (α.interval_mem hα) (Fix.add_mem
                (Fix.mul_mem (m.interval_mem hm) (x.interval_mem hx))
                (c.interval_mem hc)))
              (klEval_mem hα hp0 hp1 ha))
              (Fix.mul_mem (α.interval_mem hα)
                (klEval_mem hx hy0 hy1 hb)))
              (Fix.divInt_mem (k := 64) (Fix.flogQ_mem hr) (by norm_num)))
            (Fix.mul_mem (by simpa using Fix.ofFrac_mem (p := 11) (q := 100) (by norm_num))
              (Fix.flogQ_mem hz))

theorem scalableVertexCheck_sound (m c p y ρ z α x : QInput)
    (hm : 0 < m.den) (hc : 0 < c.den)
    (hα : 0 < α.den) (hx : 0 < x.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den)
    (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    {n : ℕ} {upper : Int}
    (h : scalableVertexCheck m c p y ρ z α x n upper = true) :
    scalableBoxExponent m.real c.real p.real y.real ρ.real z.real
      α.real x.real ≤ (upper : ℝ) / scale := by
  unfold scalableVertexCheck at h
  cases he : scalableVertexEval m c p y ρ z α x n with
  | none => simp [he] at h
  | some res =>
    have hi : res.hi ≤ upper := by simpa [he] using h
    exact (scalableVertexEval_mem m c p y ρ z α x
      hm hc hα hx hp0 hp1 hy0 hy1 he).2.trans
      (div_le_div_of_nonneg_right (by exact_mod_cast hi) scaleR_pos.le)

#print axioms scalableVertexCheck_sound

end Spin.UniformMixer
