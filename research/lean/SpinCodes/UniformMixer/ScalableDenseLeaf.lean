import SpinCodes.UniformMixer.ScalableDenseTransferExact
import SpinCodes.UniformMixer.ScalableDenseVertex
import SpinCodes.UniformMixer.ScalableDenseRateDefs

/-! Reusable soundness theorem for one checked dense rectangle. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

theorem scalable_dense_checked_leaf
    (s c p y ρ z a₀ a₁ x₀ x₁ : QInput)
    (v : Fin 3 → QInput) (C : ℝ)
    (hs : 0 < s.den) (hc : 0 < c.den)
    (ha₀ : 0 < a₀.den) (ha₁ : 0 < a₁.den)
    (hx₀ : 0 < x₀.den) (hx₁ : 0 < x₁.den)
    (hp0 : 0 < p.num) (hp1 : p.num < p.den)
    (hy0 : 0 < y.num) (hy1 : y.num < y.den)
    (hpr0 : 0 < p.real) (hpr1 : p.real < 1)
    (hyr0 : 0 < y.real) (hyr1 : y.real < 1)
    (hρ0 : 0 < ρ.real) (hz0 : 0 < z.real) (hz1 : z.real ≤ 1)
    (ha0 : 0 < a₀.real) (ha1 : a₁.real ≤ 1)
    (hx0 : 0 < x₀.real) (hx1 : x₁.real ≤ 1)
    (hC : 0 ≤ C) (hones : ∀ i, (1 : ℝ) ≤ C * (v i).real)
    (ht : transferRatCheck
      (transferQInputRat p * transferQInputRat y)
      (transferQInputRat z) (transferQInputRat ρ)
      (fun i => transferQInputRat (v i)) = true)
    (h00 : scalableVertexCheck s c p y ρ z a₀ x₀ 50
      (-100000000000000000000000) = true)
    (h01 : scalableVertexCheck s c p y ρ z a₀ x₁ 50
      (-100000000000000000000000) = true)
    (h10 : scalableVertexCheck s c p y ρ z a₁ x₀ 50
      (-100000000000000000000000) = true)
    (h11 : scalableVertexCheck s c p y ρ z a₁ x₁ 50
      (-100000000000000000000000) = true)
    {α x : ℝ}
    (hα : α ∈ Set.Icc a₀.real a₁.real)
    (hx : x ∈ Set.Icc x₀.real x₁.real) :
    ScalableDensePointRate α x s.real c.real
      (1 / 10000000) (C * (v 0).real) := by
  have htransfer : ∀ i,
      (paperSpectrumTransfer (p.real * y.real) z.real).mulVec
        (fun j => (v j).real) i ≤ ρ.real * (v i).real := by
    have hb0 : 0 ≤ p.real * y.real := mul_nonneg hpr0.le hyr0.le
    have hb1 : p.real * y.real ≤ 1 := by
      calc
        p.real * y.real ≤ 1 * y.real :=
          mul_le_mul_of_nonneg_right hpr1.le hyr0.le
        _ = y.real := one_mul _
        _ ≤ 1 := hyr1.le
    have h := transferRatCheck_sound
      (transferQInputRat p * transferQInputRat y)
      (transferQInputRat z) (transferQInputRat ρ)
      (fun i => transferQInputRat (v i))
      (by simpa [transferQInputRat_cast] using hb0)
      (by simpa [transferQInputRat_cast] using hb1)
      (by simpa [transferQInputRat_cast] using hz0.le) ht
    simpa [transferQInputRat_cast] using h
  have corner (a x : QInput) (had : 0 < a.den) (hxd : 0 < x.den)
      (h : scalableVertexCheck s c p y ρ z a x 50
        (-100000000000000000000000) = true) :
      scalableBoxExponent s.real c.real p.real y.real ρ.real z.real
        a.real x.real ≤ -(1 / 10000000 : ℝ) := by
    have hh := scalableVertexCheck_sound s c p y ρ z a x
      hs hc had hxd hp0 hp1 hy0 hy1 h
    norm_num [Spin.Numeric.scale] at hh
    exact hh
  have hbox : scalableBoxExponent s.real c.real p.real y.real
      ρ.real z.real α x ≤ -(1 / 10000000 : ℝ) := by
    have ha₀I : a₀.real ∈ Set.Icc (0 : ℝ) 1 :=
      ⟨ha0.le, (hα.1.trans hα.2).trans ha1⟩
    have ha₁I : a₁.real ∈ Set.Icc (0 : ℝ) 1 :=
      ⟨ha0.le.trans (hα.1.trans hα.2), ha1⟩
    have hx₀I : x₀.real ∈ Set.Icc (0 : ℝ) 1 :=
      ⟨hx0.le, (hx.1.trans hx.2).trans hx1⟩
    have hx₁I : x₁.real ∈ Set.Icc (0 : ℝ) 1 :=
      ⟨hx0.le.trans (hx.1.trans hx.2), hx1⟩
    have h := scalableBoxExponent_le_vertices
      s.real c.real p.real y.real ρ.real z.real
      hpr0 hpr1 hyr0 hyr1 ha₀I ha₁I hx₀I hx₁I hα hx
    apply h.trans
    apply max_le (max_le ?_ ?_) (max_le ?_ ?_)
    · exact corner a₀ x₀ ha₀ hx₀ h00
    · exact corner a₀ x₁ ha₀ hx₁ h01
    · exact corner a₁ x₀ ha₁ hx₀ h10
    · exact corner a₁ x₁ ha₁ hx₁ h11
  exact scalable_dense_witness_pointRate
    (lt_of_lt_of_le ha0 hα.1) (le_trans hα.2 ha1)
    (lt_of_lt_of_le hx0 hx.1) (le_trans hx.2 hx1)
    hpr0 hpr1 hyr0 hyr1
    hρ0 hz0 hz1 (fun i => (v i).real) hC hones htransfer hbox

#print axioms scalable_dense_checked_leaf

end Spin.UniformMixer
