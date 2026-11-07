import SpinCodes.UniformMixer.ScalableDenseVertex

/-! Reuse checked KL and logarithm intervals across dense-cover vertices. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

def scalableVertexCombine (m c α x : QInput) (a b lr lz : Fix) : Fix :=
  Fix.sub
    (Fix.add (Fix.add (Fix.add
      (Fix.mul α.interval
        (Fix.add (Fix.mul m.interval x.interval) c.interval)) a)
      (Fix.mul α.interval b)) (Fix.divInt lr 64))
    (Fix.mul (Fix.ofFrac 11 100) lz)

def scalableVertexCachedCheck (m c α x : QInput)
    (a b lr lz : Fix) (upper : Int) : Bool :=
  decide ((scalableVertexCombine m c α x a b lr lz).hi ≤ upper)

theorem scalableVertexCachedCheck_to_full
    (m c p y ρ z α x : QInput) (a b lr lz : Fix)
    {n : ℕ} {upper : Int}
    (ha : klEval α.num α.den p.num p.den n = some a)
    (hb : klEval x.num x.den y.num y.den n = some b)
    (hr : Fix.flogQ ρ.num ρ.den n = some lr)
    (hz : Fix.flogQ z.num z.den n = some lz)
    (h : scalableVertexCachedCheck m c α x a b lr lz upper = true) :
    scalableVertexCheck m c p y ρ z α x n upper = true := by
  have he : scalableVertexEval m c p y ρ z α x n =
      some (scalableVertexCombine m c α x a b lr lz) := by
    simp [scalableVertexEval, scalableVertexCombine, ha, hb, hr, hz]
  simpa [scalableVertexCheck, scalableVertexCachedCheck, he] using h

theorem scalableVertexCachedCheck_sound
    (m c p y ρ z α x : QInput) (a b lr lz : Fix)
    (hm : 0 < m.den) (hc : 0 < c.den)
    (hα : 0 < α.den) (hx : 0 < x.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den)
    (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    {n : ℕ} {upper : Int}
    (ha : klEval α.num α.den p.num p.den n = some a)
    (hb : klEval x.num x.den y.num y.den n = some b)
    (hr : Fix.flogQ ρ.num ρ.den n = some lr)
    (hz : Fix.flogQ z.num z.den n = some lz)
    (h : scalableVertexCachedCheck m c α x a b lr lz upper = true) :
    scalableBoxExponent m.real c.real p.real y.real ρ.real z.real
      α.real x.real ≤ (upper : ℝ) / scale := by
  apply scalableVertexCheck_sound m c p y ρ z α x
    hm hc hα hx hp0 hp1 hy0 hy1 (n := n) (upper := upper)
  exact scalableVertexCachedCheck_to_full m c p y ρ z α x a b lr lz
    ha hb hr hz h

#print axioms scalableVertexCachedCheck_sound

end Spin.UniformMixer
