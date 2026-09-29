import SpinCodes.Numeric.BAClosed
import SpinCodes.Numeric.FixedLog

namespace Spin.Numeric

set_option maxRecDepth 100000

lemma hEnt_concave : ConcaveOn ℝ (Set.Icc (0 : ℝ) 1) hEnt := by
  refine ⟨convex_Icc _ _, ?_⟩
  intro x hx y hy p q hp hq hpq
  have h1 := Real.convexOn_mul_log.2 hx.1 hy.1 hp hq hpq
  have h2 := Real.convexOn_mul_log.2
    (show 0 ≤ 1 - x by linarith [hx.2])
    (show 0 ≤ 1 - y by linarith [hy.2]) hp hq hpq
  simp only [smul_eq_mul] at h1 h2 ⊢
  have he : p * (1 - x) + q * (1 - y) = 1 - (p * x + q * y) := by
    nlinarith
  rw [he] at h2
  unfold hEnt
  nlinarith [h1, h2]

lemma piBA_nonpos {a c : ℝ} (ha : 0 ≤ a)
    (hlo : a / 2 ≤ c) (hhi : c ≤ 1 - a / 2) : piBA a c ≤ 0 := by
  rcases eq_or_lt_of_le ha with he | ha
  · rw [← he]
    simp [piBA]
  have hc : 0 < c := by linarith
  have hd : 0 < 1 - c := by linarith
  have hleft : a / (2 * c) ∈ Set.Icc (0 : ℝ) 1 := by
    constructor
    · positivity
    · rw [div_le_one (by positivity)]
      linarith
  have hright : a / (2 * (1 - c)) ∈ Set.Icc (0 : ℝ) 1 := by
    constructor
    · positivity
    · rw [div_le_one (by positivity)]
      linarith
  have hj := hEnt_concave.2 hleft hright hc.le hd.le (by ring : c + (1-c) = 1)
  simp only [smul_eq_mul] at hj
  have he : c * (a / (2*c)) + (1-c) * (a / (2*(1-c))) = a := by
    field_simp
    ring
  rw [he] at hj
  unfold piBA
  linarith

lemma gBA_le_half_log_two {a : ℝ} (ha : a ∈ Set.Icc (0 : ℝ) 1) :
    gBA a ≤ Real.log 2 / 2 := by
  have hg := gBA_le ha.1 ha.2 (by norm_num : (0 : ℝ) < 1)
  have he : gObj 1 a = Real.log 2 / 2 := by
    unfold gObj
    rw [show GolayG 1 = (2 : ℝ)^12 by norm_num [GolayG], Real.log_pow]
    simp
    ring
  rwa [he] at hg

lemma ba_central {a b w : ℝ} (ha : a ∈ Set.Icc (0 : ℝ) 1)
    (f1 : a/2 ≤ b) (f2 : b ≤ 1-a/2) (f3 : b/2 ≤ w) (f4 : w ≤ 1-b/2) :
    gBA a + piBA a b + piBA b w ≤ Real.log 2 / 2 := by
  have hb : 0 ≤ b := by linarith [ha.1]
  linarith [gBA_le_half_log_two ha, piBA_nonpos ha.1 f1 f2, piBA_nonpos hb f3 f4]

lemma piBA_reflect (a c : ℝ) : piBA a (1-c) = piBA a c := by
  unfold piBA
  rw [show (1 : ℝ) - (1-c) = c by ring]
  ring

lemma central_rounding : Real.log 2 / 2 ≤ (3465735902799727 : ℝ) / 10^16 := by
  have h := Fix.log2_mem.2
  have hn : Fix.log2.hi ≤ 693147180559945400000000000000 := by decide
  have hnR : ((Fix.log2.hi : Int) : ℝ) ≤ 693147180559945400000000000000 := by
    exact_mod_cast hn
  have hscale : ((scale : Int) : ℝ) = 10^30 := by norm_num [scale]
  rw [hscale] at h
  linarith

end Spin.Numeric

