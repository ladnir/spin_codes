import SpinCodes.Structured.ConcreteFixedCoefficient

/-! Coefficient extraction at prescribed marked-row weights. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel Filter
attribute [local instance] Classical.propDecidable

def rowCount {Q b : Nat} (A : Fin b → Finset (Fin Q)) (i : Fin Q) : Nat :=
  (univ.filter (fun j => i ∈ A j)).card

def profileCoefficient {Q : Nat} (u : Fin Q → ℝ) (w : Fin Q → Nat) : ℝ := ∏ i, u i ^ w i

theorem markCoefficient_eq_prod_ite {Q : Nat} (u : Fin Q → ℝ) (A : Finset (Fin Q)) :
    markCoefficient u A = ∏ i, if i ∈ A then u i else 1 := by
  simp [markCoefficient, ← prod_filter]

theorem markCoefficient_product {Q b : Nat} (u : Fin Q → ℝ) (A : Fin b → Finset (Fin Q)) :
    (∏ j, markCoefficient u (A j)) = profileCoefficient u (rowCount A) := by
  simp only [markCoefficient_eq_prod_ite]
  rw [Finset.prod_comm]
  unfold profileCoefficient rowCount
  apply prod_congr rfl
  intro i _
  simp [← prod_filter]

theorem profileCoefficient_pos {Q : Nat} {u : Fin Q → ℝ} (hu : ∀ i, 0 < u i) (w : Fin Q → Nat) :
    0 < profileCoefficient u w := prod_pos (fun i _ => pow_pos (hu i) _)

theorem profile_coefficient_bound {Q b : Nat} (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (w : Fin Q → Nat) (F : (Fin b → Finset (Fin Q)) → ℝ) (hF : ∀ A, 0 ≤ F A) :
    profileCoefficient u w * (∑ A : Fin b → Finset (Fin Q), if rowCount A = w then F A else 0) ≤
      ∑ A : Fin b → Finset (Fin Q), (∏ j, markCoefficient u (A j)) * F A := by
  rw [mul_sum]
  apply sum_le_sum
  intro A _
  by_cases hA : rowCount A = w
  · rw [if_pos hA, markCoefficient_product, hA]
  · rw [if_neg hA, mul_zero]
    exact mul_nonneg (prod_nonneg (fun j _ => markCoefficient_nonneg hu (A j))) (hF A)

theorem actual_profile_regionStream_margin {Q : Nat} {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ) (u : Fin Q → ℝ) (hu : ∀ i, 0 < u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : Nat in atTop, ∀ b : Nat, ∀ marks : Fin b → (Fin Q ↪ Fin (128 * R)), ∀ w : Fin Q → Nat,
      (∑ A : Fin b → Finset (Fin Q), if rowCount A = w then
        (Spin.piPMF (fun j => shuffleLaw ((A j).map (marks j)))).expect
          (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get ∅) else 0) ≤
        ((c + δ) ^ b / min 1 v) / profileCoefficient u w := by
  filter_upwards [actual_fugacity_product_total hθ hv hc hδ u (fun i => (hu i).le) hK] with R hR
  intro b marks w
  apply (le_div_iff₀ (profileCoefficient_pos hu w)).mpr
  rw [mul_comm]
  calc
    _ ≤ ∑ A : Fin b → Finset (Fin Q), (∏ j, markCoefficient u (A j)) *
        (Spin.piPMF (fun j => shuffleLaw ((A j).map (marks j)))).expect
          (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get ∅) := by
      apply profile_coefficient_bound u (fun i => (hu i).le)
      intro A
      apply Spin.FinPMF.expect_nonneg
      intro regions
      exact inputMoment_nonneg (Real.exp_pos _).le _ _
    _ = ∑ r, (List.ofFn (fun j => fugacityKernel θ (marks j) u)).prod ∅ r := fugacity_regionStream_moment θ marks u ∅
    _ ≤ _ := hR b marks

end Spin.Structured.Placement

