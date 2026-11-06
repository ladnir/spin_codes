import SpinCodes.UniformMixer.RegionalSiteProduct

/-!
# Passing two-state row error through full-state boundaries

The entry maps select the zero/live class of the entering state. The exit
map has entries in `[0,1]`, so a two-state row error also bounds every
full-state endpoint entry after the lift and projection.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.FiniteKernel (RowError rowAbs)

theorem projectM_unit_interval (b : Bool) (q : State22) :
    0 ≤ projectM b q ∧ projectM b q ≤ 1 := by
  have hD : (0 : ℝ) < liveCount State22 := by
    exact_mod_cast (liveCount_pos (S := State22))
  have hD1 : (1 : ℝ) ≤ liveCount State22 := by
    exact_mod_cast (liveCount_pos (S := State22))
  have hfrac0 : (0 : ℝ) ≤ 1 / liveCount State22 := by positivity
  have hfrac1 : (1 : ℝ) / liveCount State22 ≤ 1 := by
    apply (div_le_iff₀ hD).mpr
    nlinarith
  cases b
  · change 0 ≤ (if q = 0 then (1 : ℝ) else 0) ∧
      (if q = 0 then (1 : ℝ) else 0) ≤ 1
    split_ifs <;> norm_num
  · change 0 ≤ (if q = 0 then (0 : ℝ) else 1 / liveCount State22) ∧
      (if q = 0 then (0 : ℝ) else 1 / liveCount State22) ≤ 1
    split_ifs
    · norm_num
    · exact ⟨hfrac0, hfrac1⟩

theorem liftM_mul_apply (K : Matrix Bool Bool ℝ)
    (q : State22) (b : Bool) :
    (liftM * K) q b = K (if q = 0 then false else true) b := by
  change (∑ s : Bool, classLift q s * K s b) = _
  rw [Fintype.sum_bool]
  by_cases hq : q = 0
  · simp [classLift, hq]
  · simp [classLift, hq]

theorem boundary_rowError_entry {K L : Matrix Bool Bool ℝ} {ε : ℝ}
    (h : RowError K L ε) (q r : State22) :
    (liftM * K * projectM) q r ≤
      (liftM * L * projectM) q r + ε := by
  let c : Bool := if q = 0 then false else true
  have hrow :
      (∑ b : Bool, |K c b - L c b|) ≤ ε := by
    simpa only [RowError, rowAbs, Matrix.sub_apply] using h c
  have hsum :
      (∑ b : Bool, K c b * projectM b r) ≤
        (∑ b : Bool, L c b * projectM b r) + ε := by
    calc
      _ = (∑ b : Bool, L c b * projectM b r) +
            ∑ b : Bool, (K c b - L c b) * projectM b r := by
          rw [← Finset.sum_add_distrib]
          apply Finset.sum_congr rfl
          intro b _
          ring
      _ ≤ (∑ b : Bool, L c b * projectM b r) +
            ∑ b : Bool, |K c b - L c b| := by
          apply add_le_add_right
          apply Finset.sum_le_sum
          intro b _
          obtain ⟨hp0, hp1⟩ := projectM_unit_interval b r
          have h1 := mul_le_mul_of_nonneg_right
            (le_abs_self (K c b - L c b)) hp0
          have h2 := mul_le_mul_of_nonneg_left hp1
            (abs_nonneg (K c b - L c b))
          nlinarith
      _ ≤ (∑ b : Bool, L c b * projectM b r) + ε :=
          add_le_add_right hrow _
  change (∑ b : Bool, (liftM * K) q b * projectM b r) ≤
    (∑ b : Bool, (liftM * L) q b * projectM b r) + ε
  simp only [liftM_mul_apply]
  exact hsum

theorem good_support_boundary_error {R a : ℕ} (hR : 0 < R)
    (B : Spin.Structured.Placement.BlockSubset (64 * R) a)
    (hgood : ¬ badRegionalSites B.val) {γ : ℝ} (hγ : 0 ≤ γ)
    (q r : State22) :
    (liftM * coarseRegionProduct
        (fun g => Real.exp (-γ * ((g : ℝ) / R)))
        (supportMarks R (B.val.sort (· ≤ ·)) 0)
        (supportLastGap R (B.val.sort (· ≤ ·)) 0) * projectM) q r ≤
      (liftM * continuumSiteProduct γ
        (Spin.Structured.Placement.normalizedSites B) * projectM) q r +
        γ * (3 * (a + 1 : ℕ) / (R : ℝ)) := by
  have h := good_support_coarse_ideal_rowError hR B.val hgood hγ
  rw [B.property] at h
  have hh := boundary_rowError_entry h q r
  simpa only [idealSupportProduct_sorted] using hh

end Spin.UniformMixer
