import SpinCodes.Structured.Schedule

/-! The actual least-admissible-multiple schedule has rounding error less
than 24 eventually. This strengthens the existing growth and gap results;
it does not introduce a requested-length wrapper. -/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured
open Filter

/-- Minimality of `Nat.find` bounds the logarithmic schedule rounding error. -/
theorem bsched_logb_error_lt (m : ℕ) (hb : 24 < bsched m) :
    0 ≤ (bsched m : ℝ) - (39 / 4 : ℝ) * Real.logb 2 (Nsched m : ℝ) ∧
    (bsched m : ℝ) - (39 / 4 : ℝ) * Real.logb 2 (Nsched m : ℝ) < 24 := by
  have hspec : (39 / 4 : ℝ) * Real.logb 2 (Nsched m : ℝ) ≤ (bsched m : ℝ) := by
    simpa only [SchedOK, Nsched_eq, Nat.cast_mul] using bsched_spec m
  refine ⟨sub_nonneg.mpr hspec, ?_⟩
  have hn : 0 < Nat.find (exists_sched m) := by
    unfold bsched at hb
    omega
  obtain ⟨k, hk⟩ : ∃ k, Nat.find (exists_sched m) = k + 1 := by
    exact ⟨Nat.find (exists_sched m) - 1, by omega⟩
  have hfail : ¬ SchedOK m (24 * (k + 1)) :=
    Nat.find_min (exists_sched m) (by omega)
  have hprev : (24 * (k + 1) : ℝ) <
      (39 / 4 : ℝ) * Real.logb 2 ((Lsched m : ℝ) * (24 * (k + 1))) := by
    simpa only [SchedOK, Nat.cast_mul, Nat.cast_add, Nat.cast_ofNat, Nat.cast_one, not_le] using hfail
  have hmono : Real.logb 2 ((Lsched m : ℝ) * (24 * (k + 1))) ≤
      Real.logb 2 (Nsched m : ℝ) := by
    apply Real.logb_le_logb_of_le (by norm_num)
    · have hp : (0 : ℝ) < Lsched m := by exact_mod_cast Lsched_pos m
      positivity
    · have he : bsched m = 24 * (k + 2) := by simp only [bsched, hk]
      rw [Nsched_eq, Nat.cast_mul, he]
      push_cast
      have hp : (0 : ℝ) ≤ Lsched m := Nat.cast_nonneg _
      nlinarith
  have he : (bsched m : ℝ) = 24 * (k + 2 : ℝ) := by
    simp only [bsched, hk, Nat.cast_mul, Nat.cast_add, Nat.cast_ofNat]
    ring
  have hm := mul_le_mul_of_nonneg_left hmono (by norm_num : (0 : ℝ) ≤ 39 / 4)
  rw [he]
  linarith

/-- The paper's native schedule satisfies `b = (39/4) log₂ N + O(1)`,
with explicit eventual error in `[0,24)`. -/
theorem bsched_logb_error_eventually :
    ∀ᶠ m : ℕ in atTop,
      0 ≤ (bsched m : ℝ) - (39 / 4 : ℝ) * Real.logb 2 (Nsched m : ℝ) ∧
      (bsched m : ℝ) - (39 / 4 : ℝ) * Real.logb 2 (Nsched m : ℝ) < 24 := by
  filter_upwards [(tendsto_atTop.mp bsched_tendsto) 25] with m hm
  exact bsched_logb_error_lt m (by omega)

#print axioms bsched_logb_error_lt
#print axioms bsched_logb_error_eventually
end Spin.Structured
