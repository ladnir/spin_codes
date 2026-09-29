import SpinCodes.Structured.ConcreteOuterEntropyTransition

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

lemma xlogx_unit_bounds {x : ℝ} (hx : 0 ≤ x) (hx1 : x ≤ 1) :
    -1 ≤ xlogx x ∧ xlogx x ≤ 0 := by
  constructor
  · have ht := xlogx_tangent hx (by norm_num : (0 : ℝ) < 1)
    simp only [Real.log_one, mul_zero, zero_add] at ht
    linarith
  · unfold xlogx
    exact mul_nonpos_of_nonneg_of_nonpos hx (Real.log_nonpos hx hx1)

lemma xlogx_step_ordered {B x y : ℝ} (hB : 1 ≤ B)
    (hy : 0 ≤ y) (hyx : y ≤ x) (hxB : x ≤ B) (hd : x - y ≤ 1) :
    |xlogx x - xlogx y| ≤ Real.log B + 1 := by
  have hx : 0 ≤ x := hy.trans hyx
  have hlB : 0 ≤ Real.log B := Real.log_nonneg hB
  by_cases hx1 : x ≤ 1
  · obtain ⟨hxlo, hxhi⟩ := xlogx_unit_bounds hx hx1
    obtain ⟨hylo, hyhi⟩ := xlogx_unit_bounds hy (hyx.trans hx1)
    exact abs_le.mpr ⟨by linarith, by linarith⟩
  have hx1' : 1 ≤ x := by linarith
  have hxpos : 0 < x := by linarith
  have hnonneg : 0 ≤ xlogx x - xlogx y := by
    by_cases hy1 : y ≤ 1
    · have := (xlogx_unit_bounds hy hy1).2
      have : 0 ≤ xlogx x := mul_nonneg hx (Real.log_nonneg hx1')
      linarith
    · have hypos : 0 < y := by linarith
      have hlog : Real.log y ≤ Real.log x := Real.log_le_log hypos hyx
      have hmul := mul_le_mul hyx hlog (Real.log_nonneg (by linarith)) hx
      exact sub_nonneg.mpr hmul
  rw [abs_of_nonneg hnonneg]
  by_cases hy0 : y = 0
  · subst y
    have : x = 1 := by linarith
    simp [this, xlogx]
    linarith
  have hypos : 0 < y := lt_of_le_of_ne hy (Ne.symm hy0)
  have ht := Real.log_le_sub_one_of_pos (show 0 < x / y by positivity)
  rw [Real.log_div (ne_of_gt hxpos) hy0] at ht
  have hm := mul_le_mul_of_nonneg_left ht hy
  have he : y * (x / y - 1) = x - y := by field_simp
  rw [he] at hm
  have hlx : Real.log x ≤ Real.log B := Real.log_le_log hxpos hxB
  have hd0 : 0 ≤ x - y := sub_nonneg.mpr hyx
  have hml := mul_le_mul_of_nonneg_left hlx hd0
  have hmB := mul_le_mul_of_nonneg_right hd (show 0 ≤ Real.log B + 1 by linarith)
  unfold xlogx
  nlinarith

/-- A count-scale continuity bound which remains valid at zero. -/
theorem xlogx_step {B x y : ℝ} (hB : 1 ≤ B)
    (hx : 0 ≤ x) (hxB : x ≤ B) (hy : 0 ≤ y) (hyB : y ≤ B)
    (hd : |x - y| ≤ 1) :
    |xlogx x - xlogx y| ≤ Real.log B + 1 := by
  rcases le_total y x with h | h
  · exact xlogx_step_ordered hB hy h hxB ((abs_le.mp hd).2)
  · rw [abs_sub_comm]
    exact xlogx_step_ordered hB hx h hyB (by linarith [(abs_le.mp hd).1])

/-- Entropy with real, unnormalized counts. -/
def realLayerEntropy (n k : ℝ) : ℝ := xlogx n - xlogx k - xlogx (n-k)

lemma realLayerEntropy_eq {n k : ℝ} (hk : 0 ≤ k) (hkn : k ≤ n) :
    realLayerEntropy n k = n * hEnt (k/n) := by
  have hn : 0 ≤ n := hk.trans hkn
  by_cases hn0 : n = 0
  · have hk0 : k = 0 := by linarith
    simp [hn0, hk0, realLayerEntropy]
  have hnpos : 0 < n := lt_of_le_of_ne hn (Ne.symm hn0)
  by_cases hk0 : k = 0
  · simp [hk0, realLayerEntropy]
  by_cases hkn0 : k = n
  · simp [hkn0, hn0, realLayerEntropy]
  have hnk : n-k ≠ 0 := by intro h; apply hkn0; linarith
  have he : 1-k/n = (n-k)/n := by field_simp
  unfold realLayerEntropy hEnt xlogx
  rw [he, Real.log_div hk0 hn0, Real.log_div hnk hn0]
  field_simp
  ring

lemma layerEntropy_eq_real (n k : ℕ) (hk : k ≤ n) :
    layerEntropy n k = realLayerEntropy n k := by
  symm
  exact realLayerEntropy_eq (by positivity) (by exact_mod_cast hk)

/-- Three endpoint-uniform scalar bounds control each finite binomial entropy. -/
theorem realLayerEntropy_step {B n k n' k' : ℝ} (hB : 1 ≤ B)
    (hk : 0 ≤ k) (hkn : k ≤ n) (hnB : n ≤ B)
    (hk' : 0 ≤ k') (hkn' : k' ≤ n') (hnB' : n' ≤ B)
    (hn : |n-n'| ≤ 1) (hkk : |k-k'| ≤ 1)
    (hc : |(n-k)-(n'-k')| ≤ 1) :
    |realLayerEntropy n k - realLayerEntropy n' k'| ≤ 3 * (Real.log B + 1) := by
  have h₁ := xlogx_step hB (hk.trans hkn) hnB (hk'.trans hkn') hnB' hn
  have h₂ := xlogx_step hB hk (hkn.trans hnB) hk' (hkn'.trans hnB') hkk
  have h₃ := xlogx_step hB (by linarith : 0 ≤ n-k) (by linarith)
    (by linarith : 0 ≤ n'-k') (by linarith) hc
  unfold realLayerEntropy
  apply abs_le.mpr
  constructor <;> linarith [(abs_le.mp h₁).1, (abs_le.mp h₁).2,
    (abs_le.mp h₂).1, (abs_le.mp h₂).2, (abs_le.mp h₃).1, (abs_le.mp h₃).2]

end Spin.Structured.ConcreteOuter

