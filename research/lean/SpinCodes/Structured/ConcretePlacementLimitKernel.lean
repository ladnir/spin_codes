import SpinCodes.Structured.ConcretePlacementLimitBadSets

/-! Componentwise actual shuffled endpoint kernels versus their continuum Riemann sums. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel Filter MeasureTheory
open scoped Topology
attribute [local instance] Classical.propDecidable

def shuffledRegionKernel {R : Nat} (θ : ℝ) (S : Finset (Fin (128 * R))) (q r : State) : ℝ :=
  (shuffleLaw S).expect (fun T => endpointKernel (Real.exp (-(θ / (128 * R)))) (serializedInput T) q r)

theorem endpoint_coordinates {R a : Nat} (z : ℝ) (B : BlockSubset R a) (q r : State) :
    (Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
      (fun coords => endpointKernel z (serializedInput (singletonSupport (orderedBlocks B) coords)) q r) =
        uniformImpulseProduct z (fun i => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a)) q r := by
  simp only [endpointKernel_serialized_singletons]
  exact congrFun (congrFun (impulseProduct_average z (fun i => emptyGaps B i.castSucc)
    (emptyGaps B (Fin.last a))) q) r

theorem endpoint_good_coordinates_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (B : BlockSubset R a) (H : Nat) (hHR : H ≤ R) (hg : ¬ badBlockIndices (orderedBlocks B) H) (q r : State) :
    |(Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
        (fun coords => endpointKernel (Real.exp (-(θ / (128 * R))))
          (serializedInput (singletonSupport (orderedBlocks B) coords)) q r) -
        (liveLift * placementProduct θ B * liveProjection) q r| ≤ regionKernelError a R H θ := by
  rw [endpoint_coordinates]
  have he := uniformImpulseProduct_rowError (by positivity : 0 ≤ θ / (128 * (R : ℝ)))
    (fun i : Fin a => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a))
  rw [approximateImpulseProduct_eq_lift] at he
  apply (rowError_entry_general he q r).trans
  simpa only [regionKernelError, mul_assoc] using impulseErrorBudget_gap_bound
    (by positivity : 0 ≤ θ / (128 * (R : ℝ)))
    (fun i : Fin a => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a)) H R
    (fun i => good_emptyGaps_lower B H hHR hg i.castSucc) (fun i => emptyGaps_le B i.castSucc)
    (good_emptyGaps_lower B H hHR hg (Fin.last a)) (emptyGaps_le B (Fin.last a))

theorem shuffledRegionKernel_coarse_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat) (hHR : H ≤ R) (q r : State) :
    |shuffledRegionKernel θ S q r - finiteGoodKernel θ R a H q r| ≤
      placementOmissionError a R H + regionKernelError a R H θ := by
  let z := Real.exp (-(θ / (128 * (R : ℝ))))
  let f : Finset (Fin (128 * R)) → ℝ := fun T => endpointKernel z (serializedInput T) q r
  let F : BlockSubset R a → ℝ := fun B =>
    (Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
      (fun coords => f (singletonSupport (orderedBlocks B) coords))
  let w : ℝ := (128 : ℝ)^a / ((128 * R).choose a : ℝ)
  let M : ℝ := w * ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then F B else 0
  have hz0 : 0 ≤ z := (Real.exp_pos _).le
  have hz1 : z ≤ 1 := Real.exp_le_one_iff.mpr (neg_nonpos.mpr (by positivity))
  have hfirst : |shuffledRegionKernel θ S q r - M| ≤ placementOmissionError a R H := by
    exact shuffle_good_blocks_sum_error hR S hS H f
      (fun T => substochastic_entry_unit (endpointKernel_substochastic hz0 hz1 (serializedInput T)) q r)
  have hsecond : |M - finiteGoodKernel θ R a H q r| ≤ regionKernelError a R H θ := by
    exact weighted_restricted_sum_error
      (fun B : BlockSubset R a => ¬ badBlockIndices (orderedBlocks B) H)
      F (fun B => (liveLift * placementProduct θ B * liveProjection) q r)
      (w := w) (ε := regionKernelError a R H θ)
      (by positivity) (by unfold regionKernelError; positivity)
      (good_placement_mass_le_one S hS H)
      (fun B hg => endpoint_good_coordinates_error hR hθ B H hHR hg q r)
  exact (abs_sub_le _ M _).trans (add_le_add hfirst hsecond)

theorem shuffledRegionKernel_site_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat) (hHR : H ≤ R) (q r : State) :
    |shuffledRegionKernel θ S q r - finiteSiteKernel θ R a q r| ≤
      (placementOmissionError a R H + regionKernelError a R H θ) +
        ((θ * epochMean / 128) * ((a : ℝ) / R) + placementOmissionError a R H) :=
  (abs_sub_le _ (finiteGoodKernel θ R a H q r) _).trans (add_le_add
    (shuffledRegionKernel_coarse_error hR hθ S hS H hHR q r)
    (finiteGoodKernel_site_error hR hθ S hS H q r))

theorem finiteSiteKernel_tendsto (θ : ℝ) (a : Nat) (q r : State) :
    Tendsto (fun R : Nat => finiteSiteKernel θ R a q r) atTop
      (𝓝 ((liveLift * continuumRegionKernel θ a * liveProjection) q r)) := by
  rw [lifted_entry]
  by_cases hr : r = ∅
  · simpa only [finiteSiteKernel, lifted_entry, hr, ite_true] using siteProduct_sum_tendsto θ a (stateClass q) 0
  · have h := (siteProduct_sum_tendsto θ a (stateClass q) 1).div_const 524287
    simpa only [finiteSiteKernel, lifted_entry, hr, ite_false, ← sum_div, mul_div_assoc] using h

end Spin.Structured.Placement
