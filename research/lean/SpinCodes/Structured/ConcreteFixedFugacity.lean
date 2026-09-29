import SpinCodes.Structured.ConcreteFixedWeightedMargin

/-! Actual fixed-mark fugacity kernels and their continuum limit. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder ConcreteMaps FiniteKernel Filter
open scoped Topology
attribute [local instance] Classical.propDecidable

def markCoefficient {Q : Nat} (u : Fin Q → ℝ) (A : Finset (Fin Q)) : ℝ := ∏ i ∈ A, u i

def fugacityKernel {R Q : Nat} (θ : ℝ) (marks : Fin Q ↪ Fin (128 * R)) (u : Fin Q → ℝ) :
    Matrix State State ℝ := fun q r => ∑ A : Finset (Fin Q),
      markCoefficient u A * shuffledRegionKernel θ (A.map marks) q r

def fugacityContinuum {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  ∑ A : Finset (Fin Q), markCoefficient u A • continuumRegionKernel θ A.card

theorem markCoefficient_nonneg {Q : Nat} {u : Fin Q → ℝ} (hu : ∀ i, 0 ≤ u i) (A : Finset (Fin Q)) :
    0 ≤ markCoefficient u A := Finset.prod_nonneg (fun i _ => hu i)

theorem markCoefficient_sum_pos {Q : Nat} {u : Fin Q → ℝ} (hu : ∀ i, 0 ≤ u i) :
    0 < ∑ A : Finset (Fin Q), markCoefficient u A := by
  have h := Finset.single_le_sum (fun A _ => markCoefficient_nonneg hu A) (mem_univ (∅ : Finset (Fin Q)))
  simp only [markCoefficient, prod_empty] at h
  unfold markCoefficient
  linarith

theorem fugacityContinuum_lift {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (q r : State) :
    (liveLift * fugacityContinuum θ u * liveProjection) q r =
      ∑ A : Finset (Fin Q), markCoefficient u A * (liveLift * continuumRegionKernel θ A.card * liveProjection) q r := by
  simp only [lifted_entry, fugacityContinuum, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]
  by_cases hr : r = ∅
  · simp only [hr, ite_true]
  · simp only [hr, ite_false, sum_div, mul_div_assoc]

theorem fugacityKernel_uniform_limit {Q : Nat} {θ : ℝ} (hθ : 0 ≤ θ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i) {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ R : Nat in atTop, ∀ marks : Fin Q ↪ Fin (128 * R), ∀ q r : State,
      |fugacityKernel θ marks u q r - (liveLift * fugacityContinuum θ u * liveProjection) q r| < ε := by
  let C := ∑ A : Finset (Fin Q), markCoefficient u A
  have hC : 0 < C := markCoefficient_sum_pos hu
  have hh : ∀ᶠ R : Nat in atTop, ∀ A : Finset (Fin Q),
      ∀ S : Finset (Fin (128 * R)), S.card = A.card → ∀ q r : State,
        |shuffledRegionKernel θ S q r - (liveLift * continuumRegionKernel θ A.card * liveProjection) q r| < ε / (2 * C) := by
    apply Filter.eventually_all.mpr
    intro A
    exact shuffledRegionKernel_uniform_limit A.card hθ (div_pos hε (by positivity))
  filter_upwards [hh] with R hR
  intro marks q r
  rw [fugacityKernel, fugacityContinuum_lift, ← sum_sub_distrib]
  calc
    _ ≤ ∑ A : Finset (Fin Q), |markCoefficient u A * shuffledRegionKernel θ (A.map marks) q r -
        markCoefficient u A * (liveLift * continuumRegionKernel θ A.card * liveProjection) q r| :=
      Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ A : Finset (Fin Q), markCoefficient u A * (ε / (2 * C)) := by
      apply sum_le_sum
      intro A _
      rw [← mul_sub, abs_mul, abs_of_nonneg (markCoefficient_nonneg hu A)]
      exact mul_le_mul_of_nonneg_left (hR A (A.map marks) (card_map _) q r).le (markCoefficient_nonneg hu A)
    _ = ε / 2 := by rw [← sum_mul]; change C * (ε / (2 * C)) = _; field_simp
    _ < ε := half_lt_self hε

theorem fugacityKernel_nonneg {R Q : Nat} (θ : ℝ) (marks : Fin Q ↪ Fin (128 * R))
    (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i) (q r : State) : 0 ≤ fugacityKernel θ marks u q r := by
  apply sum_nonneg
  intro A _
  apply mul_nonneg (markCoefficient_nonneg hu A)
  apply Spin.FinPMF.expect_nonneg
  intro T
  exact endpointKernel_nonneg (Real.exp_pos _).le _ q r

theorem actual_fugacity_weighted_margin {Q : Nat} {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hδ : 0 < δ) (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : Nat in atTop, ∀ marks : Fin Q ↪ Fin (128 * R), ∀ q : State,
      (∑ r, fugacityKernel θ marks u q r * stateWeight v r) ≤ (c + δ) * stateWeight v q := by
  let m := min 1 v
  let W := 1 + 524287 * v
  have hm : 0 < m := lt_min (by norm_num) hv
  have hW : 0 < W := by dsimp [W]; positivity
  have hε : 0 < δ * m / W := div_pos (mul_pos hδ hm) hW
  filter_upwards [fugacityKernel_uniform_limit hθ u hu hε] with R hR
  intro marks q
  have hh : (∑ r, fugacityKernel θ marks u q r * stateWeight v r) ≤
      (∑ r, (liveLift * fugacityContinuum θ u * liveProjection) q r * stateWeight v r) + δ * m := by
    calc
      _ ≤ ∑ r, ((liveLift * fugacityContinuum θ u * liveProjection) q r + δ * m / W) * stateWeight v r := by
        apply sum_le_sum
        intro r _
        apply mul_le_mul_of_nonneg_right _ (stateWeight_pos hv r).le
        have he := (abs_lt.mp (hR marks q r)).2
        linarith
      _ = _ := by
        simp only [add_mul, sum_add_distrib, ← mul_sum, stateWeight_sum]
        change _ + δ * m / W * W = _
        rw [div_mul_cancel₀ _ hW.ne']
  have hmq : m ≤ stateWeight v q := by
    by_cases hq : q = ∅
    · simpa only [stateWeight, hq, ite_true] using min_le_left (1 : ℝ) v
    · simpa only [stateWeight, hq, ite_false] using min_le_right (1 : ℝ) v
  have hb := lifted_weight_bound hv.le hK q
  have hmδ := mul_le_mul_of_nonneg_left hmq hδ.le
  nlinarith

end Spin.Structured.Placement

