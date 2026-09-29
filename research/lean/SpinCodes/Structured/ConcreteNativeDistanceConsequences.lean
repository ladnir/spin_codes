import SpinCodes.Structured.ConcreteNativeLinearCodeDistance

/-! Success and existence consequences, independent of the pending numerical assembly. -/
noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter Finset ConcreteOuter
attribute [local instance] Classical.propDecidable

/-- The exact floor threshold is equivalent to strict relative distance above 11 percent. -/
theorem relative_distance_gt_iff_threshold (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (11/100 : ℝ) < (minimumDistance m out inner : ℝ)/(Nsched m : ℝ) ↔
      threshold m < minimumDistance m out inner := by
  have hN : (0:ℝ) < Nsched m := by exact_mod_cast Nsched_pos m
  rw [lt_div_iff₀ hN, threshold_floor, Nat.floor_lt (by positivity : (0:ℝ) ≤ 0.11*Nsched m)]
  norm_num

/-- Complementation refers to the same actual product distribution of all setup seeds. -/
theorem probability_good_distance_eq_one_sub_bad (m : ℕ) :
    ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => (11/100 : ℝ) < (minimumDistance m ω.1 ω.2 : ℝ)/(Nsched m : ℝ)) =
    1 - ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊) := by
  rw [← FinPMF.prob_compl]
  apply FinPMF.prob_congr
  intro ω
  rw [relative_distance_gt_iff_threshold, threshold_floor, not_le]

/-- An explicit vanishing failure probability gives success probability tending to one. -/
theorem success_probability_tendsto_of_failure
    (hfailure : Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊)) atTop (nhds 0)) :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => (11/100:ℝ) < (minimumDistance m ω.1 ω.2 : ℝ)/(Nsched m : ℝ)))
      atTop (nhds 1) := by
  have h := (tendsto_const_nhds (x := (1:ℝ))).sub hfailure
  simp only [sub_zero] at h
  exact h.congr (fun m => (probability_good_distance_eq_one_sub_bad m).symm)

theorem exists_seed_of_probability_pos {Ω : Type*} [Fintype Ω]
    (P : Spin.FinPMF Ω) (E : Ω → Prop) (h : 0 < P.prob E) :
    ∃ ω, E ω ∧ 0 < P.p ω := by
  unfold FinPMF.prob at h
  obtain ⟨ω, hω, hp⟩ := (sum_pos_iff_of_nonneg (fun ω _ => P.nonneg ω)).mp h
  exact ⟨ω, (mem_filter.mp hω).2, hp⟩

/-- Eventually the actual setup has a positive-mass realization with strict distance and exact rate. -/
theorem eventually_exists_good_code_of_failure
    (hfailure : Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊)) atTop (nhds 0)) :
    ∀ᶠ m in atTop, ∃ out : NativeSeed m, ∃ inner : InnerSeed m,
      (11/100:ℝ) < (minimumDistance m out inner : ℝ)/(Nsched m : ℝ) ∧
      (Module.finrank (ZMod 2) (realizedCode m out inner) : ℝ)/(Nsched m : ℝ) = 1/2 ∧
      0 < ((nativeSeedLaw m).prod
        (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).p (out, inner) := by
  have h := (success_probability_tendsto_of_failure hfailure).eventually
    (lt_mem_nhds (show (0:ℝ) < 1 by norm_num))
  filter_upwards [h] with m hm
  obtain ⟨ω, hω, hp⟩ := exists_seed_of_probability_pos _ _ hm
  exact ⟨ω.1, ω.2, hω, realizedCode_rate m ω.1 ω.2, hp⟩

end Spin.Structured.ConcreteNativeFamily
