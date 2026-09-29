import SpinCodes.Structured.ConcreteFourierProduct
import SpinCodes.Structured.ConcreteFourier

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteMaps
open scoped symmDiff

theorem orth_univ_iff {n : ℕ} (x : Finset (Fin n)) : Spin.OrthTo univ x ↔ x = ∅ := by
  constructor
  · intro h
    apply eq_empty_iff_forall_notMem.mpr
    intro i hi
    have he := h {i} (mem_univ _)
    simpa [hi] using he
  · rintro rfl
    simp [Spin.OrthTo]

theorem character_orthogonality {n : ℕ} (x y : Finset (Fin n)) :
    (∑ q : Finset (Fin n), (Spin.chi q x : ℝ)*(Spin.chi q y : ℝ)) =
      if x = y then (2 : ℝ)^n else 0 := by
  have h := Spin.sum_chi (univ : Finset (Finset (Fin n))) (fun _ _ _ _ => mem_univ _) (x ∆ y)
  simp only [orth_univ_iff, Finset.symmDiff_eq_empty] at h
  simp only [Spin.chi_symmDiff_right, card_univ, Fintype.card_finset, Fintype.card_fin] at h
  have hc := congrArg (fun a : ℤ => (a : ℝ)) h
  push_cast at hc
  exact hc

/-- Exact Fourier identity for any input law and the concrete syndrome map. -/
theorem syndrome_probability_fourier (P : FinPMF (Finset (Fin 128))) (s : Finset (Fin 19)) :
    (524288 : ℝ) * P.prob (fun X => Cset X = s) =
      ∑ q : Finset (Fin 19), (Spin.chi q s : ℝ) *
        P.expect (fun X => (Spin.chi (CtransposeSet q) X : ℝ)) := by
  have hi (X : Finset (Fin 128)) :
      (524288 : ℝ)*(if Cset X = s then 1 else 0) =
        ∑ q : Finset (Fin 19), (Spin.chi q s : ℝ)*(Spin.chi (CtransposeSet q) X : ℝ) := by
    simp_rw [C_character_pairing]
    rw [character_orthogonality]
    by_cases hx : Cset X = s
    · norm_num [hx]
    · simp [hx, Ne.symm hx]
  rw [FinPMF.prob_eq_expect_indicator]
  calc
    _ = P.expect (fun X => ∑ q : Finset (Fin 19),
        (Spin.chi q s : ℝ)*(Spin.chi (CtransposeSet q) X : ℝ)) := by
      simp only [FinPMF.expect, ← hi, Finset.mul_sum]
      apply sum_congr rfl
      intro X _
      ring
    _ = _ := by
      simp only [FinPMF.expect, Finset.mul_sum]
      rw [Finset.sum_comm]
      apply sum_congr rfl
      intro q _
      apply sum_congr rfl
      intro X _
      ring

theorem abs_real_chi {n : ℕ} (q x : Finset (Fin n)) : |(Spin.chi q x : ℝ)| = 1 := by
  rw [← Int.cast_abs, Spin.abs_chi, Int.cast_one]

/-- Every syndrome atom of an independent input law is bounded by its absolute Fourier sum. -/
theorem syndrome_probability_le (p : Fin 128 → ℝ) (hp0 : ∀ i, 0 ≤ p i)
    (hp1 : ∀ i, p i ≤ 1) (s : Finset (Fin 19)) :
    (poissonBinom p hp0 hp1).prob (fun X => Cset X = s) ≤
      (∑ q : Finset (Fin 19), ∏ i ∈ CtransposeSet q, |1-2*p i|) / 524288 := by
  apply (le_div_iff₀ (by norm_num : (0 : ℝ) < 524288)).mpr
  rw [mul_comm, syndrome_probability_fourier]
  simp_rw [poissonBinom_character]
  apply sum_le_sum
  intro q _
  calc
    _ ≤ |(Spin.chi q s : ℝ)*(∏ i ∈ CtransposeSet q, (1-2*p i))| := le_abs_self _
    _ = _ := by rw [abs_mul, abs_real_chi, one_mul, Finset.abs_prod]

end Spin.Structured.ConcreteFourier

