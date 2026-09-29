import SpinCodes.Structured.ConcreteNativeSparseSelected

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter

def nativeTail (m : ℕ) : ℝ :=
  ∑ w ∈ badWeights (bsched m) (weightWindow (bsched m)),
    Abar (nativeSeedLaw m) spectrum w

theorem raw_selection_failure_tendsto (htail : Tendsto nativeTail atTop (nhds 0)) :
    Tendsto (fun m => (nativeSeedLaw m).prob (fun seed => ¬ nativeGood m seed)) atTop (nhds 0) := by
  have h := prob_not_good_tendsto_zero nativeSeedLaw (fun m => spectrum (k := nativeBlocks m))
    bsched (fun m => weightWindow (bsched m)) bsched_pos
    (fun _ => Finset.filter_subset _ _) bsched_tendsto htail
  apply h.congr
  intro m
  apply FinPMF.prob_congr
  intro seed
  simp only [nativeGood, constituentGood, nativeSeedLaw, native_width]

theorem raw_selection_positive_eventually (htail : Tendsto nativeTail atTop (nhds 0)) :
    ∀ᶠ m in atTop, 0 < (nativeSeedLaw m).prob (nativeGood m) := by
  have h := (raw_selection_failure_tendsto htail).eventually
    (gt_mem_nhds (show (0 : ℝ) < 1 by norm_num))
  filter_upwards [h] with m hm
  rw [FinPMF.prob_compl] at hm
  linarith

theorem concrete_selection_failure_tendsto (htail : Tendsto nativeTail atTop (nhds 0)) :
    Tendsto concreteFamily.probNotGood atTop (nhds 0) :=
  concrete_probNotGood_tendsto (raw_selection_failure_tendsto htail)

def nativeCut (m : ℕ) : ℕ := Lsched m / 10000

theorem nativeCut_spec :
    ∀ᶠ m in atTop, 4095 ≤ nativeCut m ∧ nativeCut m ≤ Lsched m := by
  apply eventually_atTop.mpr
  refine ⟨320000, fun m hm => ?_⟩
  unfold nativeCut Lsched
  omega

theorem nativeCut_density (m Q : ℕ) (hQ : Q ≤ nativeCut m) :
    (Q : ℝ)/(Lsched m : ℝ) ≤ 1/10000 := by
  have hL : (0 : ℝ) < Lsched m := by exact_mod_cast Lsched_pos m
  have hQ' : 10000*Q ≤ Lsched m := by unfold nativeCut at hQ; omega
  have hR : (10000 : ℝ)*(Q : ℝ) ≤ Lsched m := by exact_mod_cast hQ'
  apply (div_le_iff₀ hL).mpr
  linarith

/-- The selected native family has the uniform sparse regime, conditional only
on the two remaining outer asymptotic facts: vanishing tail and shell envelope. -/
theorem concrete_sparse_eventually (htail : Tendsto nativeTail atTop (nhds 0))
    {r : ℕ → ℝ} (hr : Tendsto r atTop (nhds 0))
    (hbar : ∀ᶠ m in atTop, ∀ w ∈ weightWindow (bsched m),
      Abar (nativeSeedLaw m) spectrum w ≤
        spectrumRatioBound (bsched m) (r m) * ((bsched m).choose w : ℝ)) :
    ∀ᶠ m in atTop, ∀ Q ∈ Ico 4096 (nativeCut m + 1),
      concreteFamily.EZ m Q ≤ Real.exp (-(0.006 : ℝ)*(Q : ℝ)*bsched m) := by
  filter_upwards [raw_selection_positive_eventually htail,
    SparseRate.native_thresholds (selectedRemainder_tendsto hr), hbar] with m hpos hth hbar
  intro Q hQ
  have hQlo := (mem_Ico.mp hQ).1
  have hQcut : Q ≤ nativeCut m := by have := (mem_Ico.mp hQ).2; omega
  have hQL : Q ≤ Lsched m := le_trans hQcut (Nat.div_le_self _ _)
  simpa only [show (3/500 : ℝ) = 0.006 by norm_num] using
    selected_EZ_sparse_rate m Q hpos hQlo hQL hth.1 (Real.exp_pos _).le hbar
      (selected_row_cost m (r m)).le (nativeCut_density m Q hQcut) hth.2

end Spin.Structured.ConcreteNativeFamily
