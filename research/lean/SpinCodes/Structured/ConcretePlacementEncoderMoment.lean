import SpinCodes.Structured.ConcretePlacementEncoderStream

/-! Actual serialized region moments and their deterministic coarse-product approximation. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder FiniteKernel
attribute [local instance] Classical.propDecidable

def regionMoment {R : Nat} (z : ℝ) (T : Finset (Fin (128 * R))) (q : State) : ℝ :=
  inputMoment z (serializedInput T) q

theorem regionMoment_unit {R : Nat} {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (T : Finset (Fin (128 * R))) (q : State) : 0 ≤ regionMoment z T q ∧ regionMoment z T q ≤ 1 := by
  constructor
  · exact inputMoment_nonneg hz _ q
  · exact (Spin.FinPMF.expect_mono _ (fun _ => pow_le_one₀ hz hz1)).trans_eq (Spin.FinPMF.expect_const _ 1)

theorem endpointKernel_serialized_singletons {R a : Nat} (z : ℝ) (B : BlockSubset R a)
    (coords : Fin a → Fin 128) :
    endpointKernel z (serializedInput (singletonSupport (orderedBlocks B) coords)) =
      impulseProduct z (fun i => emptyGaps B i.castSucc) coords (emptyGaps B (Fin.last a)) := by
  rw [endpointKernel_eq_path, serializedInput_list_eq_impulseInputs, pathKernel_impulseInputs]

theorem regionMoment_coordinates {R a : Nat} (z : ℝ) (B : BlockSubset R a) (q : State) :
    (Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
      (fun coords => regionMoment z (singletonSupport (orderedBlocks B) coords) q) =
        ∑ r, uniformImpulseProduct z (fun i => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a)) q r := by
  simp only [regionMoment, ← endpointKernel_total, endpointKernel_serialized_singletons,
    Spin.FinPMF.expect_sum]
  apply sum_congr rfl
  intro r _
  exact congrFun (congrFun (impulseProduct_average z (fun i => emptyGaps B i.castSucc)
    (emptyGaps B (Fin.last a))) q) r

theorem rowError_total {α : Type*} [Fintype α] [DecidableEq α]
    {K L : Matrix α α ℝ} {ε : ℝ} (h : RowError K L ε) (q : α) :
    |(∑ r, K q r) - ∑ r, L q r| ≤ ε := by
  rw [← sum_sub_distrib]
  exact (Finset.abs_sum_le_sum_abs _ _).trans (h q)

def coarseRegionMoment {R a : Nat} (θ : ℝ) (B : BlockSubset R a) (q : State) : ℝ :=
  ∑ r, (liveLift * placementProduct θ B * liveProjection) q r

theorem regionMoment_coordinates_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (B : BlockSubset R a) (q : State) :
    |(Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
        (fun coords => regionMoment (Real.exp (-(θ / (128 * R))))
          (singletonSupport (orderedBlocks B) coords) q) - coarseRegionMoment θ B q| ≤
      impulseErrorBudget (θ / (128 * R)) (fun i => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a)) := by
  rw [regionMoment_coordinates]
  have hc : 0 ≤ θ / (128 * (R : ℝ)) := by positivity
  have h := uniformImpulseProduct_rowError hc (fun i => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a))
  rw [approximateImpulseProduct_eq_lift] at h
  exact rowError_total h q

def regionKernelError (a R H : Nat) (θ : ℝ) : ℝ :=
  ((a : ℝ) + 1) * 524288 * ((θ / (128 * R)) * (128 * Real.sqrt (3 * R)) + (1 / 2 : ℝ)^H) +
    (a : ℝ) * (128 * (θ / (128 * R)))

theorem regionMoment_good_coordinates_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (B : BlockSubset R a) (H : Nat) (hHR : H ≤ R) (hg : ¬ badBlockIndices (orderedBlocks B) H) (q : State) :
    |(Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
        (fun coords => regionMoment (Real.exp (-(θ / (128 * R))))
          (singletonSupport (orderedBlocks B) coords) q) - coarseRegionMoment θ B q| ≤
      regionKernelError a R H θ := by
  apply (regionMoment_coordinates_error hR hθ B q).trans
  simpa only [regionKernelError, mul_assoc] using impulseErrorBudget_gap_bound
    (by positivity : 0 ≤ θ / (128 * (R : ℝ)))
    (fun i : Fin a => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a)) H R
    (fun i => good_emptyGaps_lower B H hHR hg i.castSucc)
    (fun i => emptyGaps_le B i.castSucc)
    (good_emptyGaps_lower B H hHR hg (Fin.last a)) (emptyGaps_le B (Fin.last a))

end Spin.Structured.Placement

