import SpinCodes.Structured.ConcreteFixedProductMargin
import SpinCodes.Structured.ConcreteMarkedUniformSubset
import SpinCodes.Structured.ConcreteMarkedShuffle

/-! Actual fair marked inputs: exact region kernel and multiplicative margin. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder ConcreteMarked FiniteKernel Filter
attribute [local instance] Classical.propDecidable

def fairRegionKernel {R : Nat} (θ : ℝ) (S : Finset (Fin (128 * R))) : Matrix State State ℝ :=
  matrixExpect (regionLaw S) (fun T => endpointKernel (Real.exp (-(θ / (128 * R)))) (serializedInput T))

theorem fairRegionKernel_eq_fugacity {R Q : Nat} (θ : ℝ) (marks : Fin Q ↪ Fin (128 * R)) (q r : State) :
    fairRegionKernel θ (univ.map marks) q r = fugacityKernel θ marks (fun _ => 1) q r / (2:ℝ)^Q := by
  unfold fairRegionKernel matrixExpect
  rw [← shuffled_fairOnMarks, shuffledLaw, Spin.FinPMF.expect_bind, fairOnMarks_expect_uniform_sum]
  simp only [fugacityKernel, markCoefficient, prod_const_one, one_mul, shuffledRegionKernel]

theorem fair_regionStream_endpoint {R b : Nat} (θ : ℝ) (S : Fin b → Finset (Fin (128 * R))) :
    matrixExpect (Spin.piPMF (fun j => regionLaw (S j)))
      (fun regions => endpointKernel (Real.exp (-(θ / (128 * R)))) (regionStream regions).get) =
        (List.ofFn (fun j => fairRegionKernel θ (S j))).prod := by
  simp only [endpointKernel_regionStream]
  exact matrixExpect_pi_product (fun j => regionLaw (S j))
    (fun j T => endpointKernel (Real.exp (-(θ / (128 * R)))) (serializedInput T))

theorem fair_regionStream_moment {R b : Nat} (θ : ℝ) (S : Fin b → Finset (Fin (128 * R))) (q : State) :
    (Spin.piPMF (fun j => regionLaw (S j))).expect
      (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get q) =
        ∑ r, (List.ofFn (fun j => fairRegionKernel θ (S j))).prod q r := by
  simp only [← endpointKernel_total, Spin.FinPMF.expect_sum]
  apply sum_congr rfl
  intro r _
  exact congrFun (congrFun (fair_regionStream_endpoint θ S) q) r

theorem actual_fair_regionStream_margin {Q : Nat} {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ (fun _ : Fin Q => 1))) :
    ∀ᶠ R : Nat in atTop, ∀ b : Nat, ∀ marks : Fin b → (Fin Q ↪ Fin (128 * R)),
      (Spin.piPMF (fun j => regionLaw (univ.map (marks j)))).expect
        (fun regions => inputMoment (Real.exp (-(θ / (128 * R)))) (regionStream regions).get ∅) ≤
          ((c + δ) / (2:ℝ)^Q) ^ b / min 1 v := by
  filter_upwards [actual_fugacity_weighted_margin hθ hv hδ (fun _ : Fin Q => 1) (by intro i; norm_num) hK] with R hR
  intro b marks
  rw [fair_regionStream_moment]
  have hk (j : Fin b) : WeightedBound (stateWeight v) ((c + δ) / (2:ℝ)^Q)
      (fairRegionKernel θ (univ.map (marks j))) := by
    constructor
    · intro q r
      rw [fairRegionKernel_eq_fugacity]
      exact div_nonneg (fugacityKernel_nonneg θ (marks j) _ (by intro i; norm_num) q r) (by positivity)
    · intro q
      simp only [fairRegionKernel_eq_fugacity, div_mul_eq_mul_div, ← sum_div]
      simpa only [div_mul_eq_mul_div] using div_le_div_of_nonneg_right (hR (marks j) q) (by positivity : 0 ≤ (2:ℝ)^Q)
  have hp := WeightedBound.ofFn (by positivity : 0 ≤ (c + δ) / (2:ℝ)^Q)
    (fun j => fairRegionKernel θ (univ.map (marks j))) hk
  have hw (q : State) : min 1 v ≤ stateWeight v q := by
    unfold stateWeight
    split_ifs
    · exact min_le_left _ _
    · exact min_le_right _ _
  simpa [stateWeight] using hp.total (lt_min (by norm_num) hv) hw ∅

end Spin.Structured.Placement
