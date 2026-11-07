import SpinCodes.UniformMixer.ScalableFixedWeightedMargin
import SpinCodes.Structured.ConcreteFixedFugacity

/-! Fixed-row fugacity kernel for the actual Scalable regional shuffle.
Each marked row may contribute its bit or remain absent from a region. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.Placement Filter

def scalableFugacityKernel {m Q : ℕ} (θ : ℝ)
    (marks : Fin Q ↪ Fin (Lsched m)) (u : Fin Q → ℝ) :
    Matrix State22 State22 ℝ :=
  fun q r => ∑ A : Finset (Fin Q),
    markCoefficient u A *
      scalableAveragedRegion m
        (Real.exp (-(θ / (Lsched m : ℝ)))) (A.map marks) q r

def scalableFugacityContinuum {Q : ℕ} (θ : ℝ)
    (u : Fin Q → ℝ) : Matrix State22 State22 ℝ :=
  fun q r => ∑ A : Finset (Fin Q),
    markCoefficient u A * continuumEndpoint θ A.card q r

theorem scalableFugacityKernel_nonneg {m Q : ℕ}
    {θ : ℝ} (hθ : 0 ≤ θ) {u : Fin Q → ℝ}
    (hu : ∀ i, 0 ≤ u i) (marks : Fin Q ↪ Fin (Lsched m))
    (q r : State22) :
    0 ≤ scalableFugacityKernel θ marks u q r := by
  have hL : (0 : ℝ) < (Lsched m : ℝ) := by
    unfold Lsched
    positivity
  have hz1 : Real.exp (-(θ / (Lsched m : ℝ))) ≤ 1 :=
    Real.exp_le_one_iff.mpr (neg_nonpos.mpr (div_nonneg hθ hL.le))
  unfold scalableFugacityKernel
  apply Finset.sum_nonneg
  intro A _
  exact mul_nonneg (markCoefficient_nonneg hu A)
    (scalableAveragedRegion_nonneg m _ (Real.exp_pos _).le hz1
      (A.map marks) q r)

theorem scalableFugacityKernel_uniform_upper {Q : ℕ}
    {θ : ℝ} (hθ : 0 ≤ θ) (u : Fin Q → ℝ)
    (hu : ∀ i, 0 ≤ u i) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin Q ↪ Fin (Lsched m), ∀ q r : State22,
        scalableFugacityKernel θ marks u q r <
          scalableFugacityContinuum θ u q r + ε := by
  let C : ℝ := ∑ A : Finset (Fin Q), markCoefficient u A
  have hC : 0 < C := markCoefficient_sum_pos hu
  have hsmall : 0 < ε / (2 * C) :=
    div_pos hε (mul_pos (by norm_num) hC)
  have hAll : ∀ᶠ m : ℕ in atTop,
      ∀ A : Finset (Fin Q),
        ∀ S : Finset (Fin (Lsched m)), S.card = A.card →
          ∀ q r : State22,
            scalableAveragedRegion m
              (Real.exp (-(θ / (Lsched m : ℝ)))) S q r <
            continuumEndpoint θ A.card q r + ε / (2 * C) := by
    apply Filter.eventually_all.mpr
    intro A
    exact scalable_fixed_region_continuum_upper A.card hθ hsmall
  filter_upwards [hAll] with m hm
  intro marks q r
  have hterm (A : Finset (Fin Q)) :
      markCoefficient u A *
        scalableAveragedRegion m
          (Real.exp (-(θ / (Lsched m : ℝ)))) (A.map marks) q r ≤
      markCoefficient u A *
        (continuumEndpoint θ A.card q r + ε / (2 * C)) := by
    exact mul_le_mul_of_nonneg_left
      (hm A (A.map marks) (by simp) q r).le
      (markCoefficient_nonneg hu A)
  have hsum := Finset.sum_le_sum
    (fun A (_ : A ∈ (Finset.univ : Finset (Finset (Fin Q)))) => hterm A)
  change scalableFugacityKernel θ marks u q r ≤
    ∑ A : Finset (Fin Q),
      markCoefficient u A * (continuumEndpoint θ A.card q r + ε / (2 * C)) at hsum
  have heq :
      (∑ A : Finset (Fin Q),
        markCoefficient u A * (continuumEndpoint θ A.card q r + ε / (2 * C))) =
      scalableFugacityContinuum θ u q r + ε / 2 := by
    unfold scalableFugacityContinuum
    simp only [mul_add, Finset.sum_add_distrib, ← Finset.sum_mul]
    change _ + C * (ε / (2 * C)) = _
    field_simp
  calc
    scalableFugacityKernel θ marks u q r ≤
        ∑ A : Finset (Fin Q),
          markCoefficient u A * (continuumEndpoint θ A.card q r + ε / (2 * C)) := hsum
    _ = scalableFugacityContinuum θ u q r + ε / 2 := heq
    _ < scalableFugacityContinuum θ u q r + ε := by linarith

theorem scalable_fugacity_weighted_margin {Q : ℕ}
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (w : State22 → ℝ) (hw : ∀ r, 0 < w r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ w r)
    (hK : ∀ q : State22,
      (∑ r, scalableFugacityContinuum θ u q r * w r) ≤
        c * w q) :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin Q ↪ Fin (Lsched m),
        Spin.Structured.FiniteKernel.WeightedBound w (c + δ)
          (scalableFugacityKernel θ marks u) := by
  let W : ℝ := ∑ r : State22, w r
  have hW : 0 < W := finite_positive_sum w hw
  have hε : 0 < δ * m₀ / W :=
    div_pos (mul_pos hδ hm₀) hW
  filter_upwards [scalableFugacityKernel_uniform_upper hθ u hu hε]
    with m hm
  intro marks
  constructor
  · intro q r
    exact scalableFugacityKernel_nonneg hθ hu marks q r
  · intro q
    have hsum :
        (∑ r, scalableFugacityKernel θ marks u q r * w r) ≤
        (∑ r, (scalableFugacityContinuum θ u q r + δ * m₀ / W) * w r) := by
      apply Finset.sum_le_sum
      intro r _
      exact mul_le_mul_of_nonneg_right (hm marks q r).le (hw r).le
    calc
      _ ≤ (∑ r, (scalableFugacityContinuum θ u q r + δ * m₀ / W) * w r) :=
        hsum
      _ = (∑ r, scalableFugacityContinuum θ u q r * w r) + δ * m₀ := by
        simp only [add_mul, Finset.sum_add_distrib, ← Finset.mul_sum]
        change _ + (δ * m₀ / W) * W = _
        rw [div_mul_cancel₀ _ hW.ne']
      _ ≤ c * w q + δ * w q :=
        add_le_add (hK q) (mul_le_mul_of_nonneg_left (hmw q) hδ.le)
      _ = (c + δ) * w q := by ring

theorem scalable_fugacity_product_total {Q : ℕ}
    {θ c δ m₀ : ℝ} (hθ : 0 ≤ θ) (hc : 0 ≤ c) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (w : State22 → ℝ) (hw : ∀ r, 0 < w r)
    (hm₀ : 0 < m₀) (hmw : ∀ r, m₀ ≤ w r)
    (hK : ∀ q : State22,
      (∑ r, scalableFugacityContinuum θ u q r * w r) ≤
        c * w q) :
    ∀ᶠ m : ℕ in atTop, ∀ b : ℕ,
      ∀ marks : Fin b → (Fin Q ↪ Fin (Lsched m)),
        (∑ r : State22,
          (List.ofFn (fun j => scalableFugacityKernel θ (marks j) u)).prod 0 r) ≤
          (c + δ) ^ b * w 0 / m₀ := by
  filter_upwards [scalable_fugacity_weighted_margin hθ hδ
    u hu w hw hm₀ hmw hK] with m hm
  intro b marks
  have hprod : Spin.Structured.FiniteKernel.WeightedBound w
      ((c + δ) ^ b)
      (List.ofFn (fun j => scalableFugacityKernel θ (marks j) u)).prod := by
    apply Spin.Structured.FiniteKernel.WeightedBound.ofFn (by linarith)
    exact fun j => hm (marks j)
  exact hprod.total hm₀ hmw 0

end Spin.UniformMixer
