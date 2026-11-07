import SpinCodes.UniformMixer.ScalableRoutedMoment

/-!
# Route domination for the realized Scalable encoder moment

The previously checked route comparison applies to the new encoder only
after its conditional weighted moment is shown nonnegative. That
nonnegativity follows from the exact sampled-map moment identity.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute

def scalableRoutedMoment (m : ℕ) (z : ℝ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) : ℝ :=
  ∑ r : State22,
    inputListKernel z (List.ofFn (scalableFlatInput m regions)) 0 r

theorem scalableRoutedMoment_nonneg (m : ℕ) {z : ℝ} (hz : 0 ≤ z)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    0 ≤ scalableRoutedMoment m z regions := by
  unfold scalableRoutedMoment
  rw [← scalable_fixed_route_moment]
  exact FinPMF.expect_nonneg _ (fun _ => pow_nonneg hz _)

theorem scalable_route_iid_bound (m : ℕ) {z : ℝ} (hz : 0 ≤ z)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m)))
    (hw0 : 0 < Spin.Structured.ConcreteRoute.totalWeight rows)
    (hw1 : Spin.Structured.ConcreteRoute.totalWeight rows < Lsched m * bsched m) :
    (Spin.Structured.ConcreteRoute.law rows).expect
      (scalableRoutedMoment m z) ≤
    (((bsched m : ℝ) + 1) ^ Spin.Structured.ConcreteRoute.activeRows rows *
      ((Lsched m : ℝ) + 1) ^ bsched m) *
      (Spin.Structured.ConcreteRoute.iidLaw (Lsched m) (bsched m)
        (Spin.Structured.ConcreteRoute.density rows)
        (Spin.Structured.ConcreteRoute.density_nonneg rows)
        (Spin.Structured.ConcreteRoute.density_le_one
          (Lsched_pos m) (bsched_pos m) rows)).expect
        (scalableRoutedMoment m z) := by
  have h := Spin.Structured.ConcreteRoute.permutation_expect_le_of_weight
    (Lsched_pos m) (bsched_pos m) rows hw0 hw1
    (scalableRoutedMoment m z)
    (scalableRoutedMoment_nonneg m hz)
  calc
    _ = ((Spin.Structured.ConcreteRoute.seedLaw (Lsched m) (bsched m)).map
        (Spin.Structured.ConcreteRoute.routeEval rows)).expect
        (scalableRoutedMoment m z) := by
        rw [Spin.Structured.ConcreteRoute.permutation_law_eq]
    _ = (Spin.Structured.ConcreteRoute.seedLaw (Lsched m) (bsched m)).expect
        (fun seed => scalableRoutedMoment m z
          (Spin.Structured.ConcreteRoute.routeEval rows seed)) :=
        FinPMF.expect_map _ _ _
    _ ≤ _ := h

/-- A selected shared outer gives the nontrivial weight range required by
the route comparison for every nonzero message. -/
theorem scalable_good_outer_weight (m : ℕ)
    {out : Spin.Structured.ConcreteOuter.NativeSeed m}
    (hg : Spin.Structured.ConcreteNativeFamily.nativeGood m out)
    {message : Fin (Nsched m / 2) → Bool} (hm : message ≠ 0) :
    0 < totalWeight (nativeRows m out message) ∧
      totalWeight (nativeRows m out message) < Lsched m * bsched m := by
  let e := Spin.Structured.ConcreteBinaryEncoding.wordEquiv (Nsched m / 2)
  have hm' : e message ≠ 0 := by
    intro h
    exact hm (e.injective (h.trans (by simp [e])))
  have h := Spin.Structured.ConcreteNativeFamily.native_good_total_weight
    m hg hm'
  change 0 < totalWeight (nativeRows m out (e.symm (e message))) ∧
    totalWeight (nativeRows m out (e.symm (e message))) <
      Lsched m * bsched m at h
  rw [e.symm_apply_apply] at h
  exact h

theorem scalable_good_outer_route_iid_bound (m : ℕ)
    {z : ℝ} (hz : 0 ≤ z)
    {out : Spin.Structured.ConcreteOuter.NativeSeed m}
    (hg : Spin.Structured.ConcreteNativeFamily.nativeGood m out)
    {message : Fin (Nsched m / 2) → Bool} (hm : message ≠ 0) :
    (Spin.Structured.ConcreteRoute.law (nativeRows m out message)).expect
      (scalableRoutedMoment m z) ≤
    (((bsched m : ℝ) + 1) ^
        activeRows (nativeRows m out message) *
      ((Lsched m : ℝ) + 1) ^ bsched m) *
      (Spin.Structured.ConcreteRoute.iidLaw (Lsched m) (bsched m)
        (density (nativeRows m out message))
        (density_nonneg (nativeRows m out message))
        (density_le_one (Lsched_pos m) (bsched_pos m)
          (nativeRows m out message))).expect
        (scalableRoutedMoment m z) := by
  obtain ⟨hw0, hw1⟩ := scalable_good_outer_weight m hg hm
  exact scalable_route_iid_bound m hz
    (nativeRows m out message) hw0 hw1

end Spin.UniformMixer
