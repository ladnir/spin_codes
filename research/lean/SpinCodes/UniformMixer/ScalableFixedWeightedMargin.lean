import SpinCodes.UniformMixer.ScalableFixedRegionLimit
import SpinCodes.Structured.ConcreteFixedProductMargin

/-! A fixed-occupation continuum margin transfers to all Scalable regions
and composes across their full-state boundaries. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.FiniteKernel Filter

private theorem finite_positive_sum {ι : Type*} [Fintype ι]
    [Nonempty ι] (w : ι → ℝ) (hw : ∀ i, 0 < w i) :
    0 < ∑ i, w i := by
  classical
  let i : ι := Classical.choice inferInstance
  have h := Finset.single_le_sum
    (fun j _ => (hw j).le)
    (Finset.mem_univ i)
  exact (hw i).trans_le h

theorem scalableAveragedRegion_nonneg (m : ℕ) (z : ℝ)
    (hz : 0 ≤ z) (hz1 : z ≤ 1) (S : Finset (Fin (Lsched m)))
    (q r : State22) :
    0 ≤ scalableAveragedRegion m z S q r := by
  unfold scalableAveragedRegion matrixExpect22
  apply Spin.FinPMF.expect_nonneg
  intro T
  change 0 ≤ chronologicalRegionKernel z
    (scalableRegionSupport m T) q r
  exact (chronologicalRegionKernel_unit_interval z hz
    hz1 (scalableRegionSupport m T) q r).1

theorem scalable_fixed_region_weighted_margin (a : ℕ)
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hδ : 0 < δ)
    (w : State22 → ℝ) (hw : ∀ r, 0 < w r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ w r)
    (hK : ∀ q : State22,
      (∑ r, continuumEndpoint θ a q r * w r) ≤ c * w q) :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = a →
        WeightedBound w (c + δ)
          (scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ)))) S) := by
  let W : ℝ := ∑ r : State22, w r
  have hW : 0 < W := finite_positive_sum w hw
  have hε : 0 < δ * m₀ / W :=
    div_pos (mul_pos hδ hm₀) hW
  filter_upwards [scalable_fixed_region_continuum_upper a hθ hε]
    with m hR
  intro S hS
  have hL : (0 : ℝ) < (Lsched m : ℝ) := by
    unfold Lsched
    positivity
  have hz1 : Real.exp (-(θ / (Lsched m : ℝ))) ≤ 1 :=
    Real.exp_le_one_iff.mpr (neg_nonpos.mpr (div_nonneg hθ hL.le))
  constructor
  · intro q r
    exact scalableAveragedRegion_nonneg m _ (Real.exp_pos _).le hz1 S q r
  · intro q
    have hsum :
        (∑ r, scalableAveragedRegion m
          (Real.exp (-(θ / (Lsched m : ℝ)))) S q r * w r) ≤
        (∑ r, (continuumEndpoint θ a q r + δ * m₀ / W) * w r) := by
      apply Finset.sum_le_sum
      intro r _
      exact mul_le_mul_of_nonneg_right (hR S hS q r).le (hw r).le
    calc
      _ ≤ (∑ r, (continuumEndpoint θ a q r + δ * m₀ / W) * w r) :=
        hsum
      _ = (∑ r, continuumEndpoint θ a q r * w r) + δ * m₀ := by
        simp only [add_mul, Finset.sum_add_distrib, ← Finset.mul_sum]
        change _ + (δ * m₀ / W) * W = _
        rw [div_mul_cancel₀ _ hW.ne']
      _ ≤ c * w q + δ * w q :=
        add_le_add (hK q) (mul_le_mul_of_nonneg_left (hmw q) hδ.le)
      _ = (c + δ) * w q := by ring

theorem scalable_fixed_product_weighted_margin (Q : ℕ)
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hδ : 0 < δ)
    (w : State22 → ℝ) (hw : ∀ r, 0 < w r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ w r)
    (hK : ∀ a ≤ Q, ∀ q : State22,
      (∑ r, continuumEndpoint θ a q r * w r) ≤ c * w q) :
    ∀ᶠ m : ℕ in atTop, ∀ b : ℕ,
      ∀ S : Fin b → Finset (Fin (Lsched m)),
        (∀ j, (S j).card ≤ Q) →
        WeightedBound w ((c + δ) ^ b)
          (List.ofFn (fun j => scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ)))) (S j))).prod := by
  have hAll : ∀ᶠ m : ℕ in atTop,
      ∀ a : Fin (Q + 1),
        ∀ S : Finset (Fin (Lsched m)), S.card = a.val →
          WeightedBound w (c + δ)
            (scalableAveragedRegion m
              (Real.exp (-(θ / (Lsched m : ℝ)))) S) := by
    apply Filter.eventually_all.mpr
    intro a
    exact scalable_fixed_region_weighted_margin a.val hθ hδ
      w hw hm₀ hmw (hK a.val (by omega))
  filter_upwards [hAll] with m hm
  intro b S hS
  apply WeightedBound.ofFn (by linarith)
  intro j
  have hj := hS j
  exact hm ⟨(S j).card, by omega⟩ (S j) rfl

theorem scalable_fixed_product_total (Q : ℕ)
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hδ : 0 < δ)
    (w : State22 → ℝ) (hw : ∀ r, 0 < w r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ w r)
    (hK : ∀ a ≤ Q, ∀ q : State22,
      (∑ r, continuumEndpoint θ a q r * w r) ≤ c * w q) :
    ∀ᶠ m : ℕ in atTop, ∀ b : ℕ,
      ∀ S : Fin b → Finset (Fin (Lsched m)),
        (∀ j, (S j).card ≤ Q) →
        (∑ r : State22,
          (List.ofFn (fun j => scalableAveragedRegion m
            (Real.exp (-(θ / (Lsched m : ℝ)))) (S j))).prod 0 r) ≤
          (c + δ) ^ b * w 0 / m₀ := by
  filter_upwards [scalable_fixed_product_weighted_margin Q hθ hc hδ
    w hw hm₀ hmw hK] with m hm
  intro b S hS
  exact (hm b S hS).total hm₀ hmw 0

end Spin.UniformMixer
