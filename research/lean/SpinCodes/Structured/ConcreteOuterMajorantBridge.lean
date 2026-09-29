import SpinCodes.Structured.ConcreteOuterEnvelopeShiftSupport
import SpinCodes.Structured.ConcreteOuterMajorantRegularity

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset

lemma normalized_feasible {b a c : ℝ} (hb : 0 < b)
    (hlo : a/2 ≤ c) (hhi : c ≤ b-a/2) :
    (a/b)/2 ≤ c/b ∧ c/b ≤ 1-(a/b)/2 := by
  have h₁ : (a/b)/2 = (a/2)/b := by ring
  have h₂ : 1-(a/b)/2 = (b-a/2)/b := by field_simp
  rw [h₂, h₁]
  exact ⟨(div_le_div_iff_of_pos_right hb).mpr hlo, (div_le_div_iff_of_pos_right hb).mpr hhi⟩

lemma refined_count_step {b t w : ℝ} (hb : 0 < b) (htw : t ≤ w) (hwt : w-t ≤ 1/2) :
    b * Spin.Majorant.refined.toFun (t/b) ≤ b * Spin.Majorant.refined.toFun (w/b) + 1 := by
  have hh := mul_le_mul_of_nonneg_left (refined_step (t/b) (w/b)) hb.le
  have he : b*(2*|t/b-w/b|) ≤ 1 := by
    rw [← sub_div, abs_div, abs_of_pos hb, abs_of_nonpos (by linarith : t-w ≤ 0)]
    have hcalc : b*(2*(-(t-w)/b)) = 2*(w-t) := by field_simp; ring
    rw [hcalc]
    linarith
  nlinarith

lemma clippedOutput_window {b a c : ℕ} (hb : 0 < b) (hab : a ≤ b)
    (hc : (c:ℝ)/b ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    clippedOutput b a c / b ∈ Set.Icc ((13:ℝ)/125) (112/125) := by
  have hbR : (0:ℝ)<b := by exact_mod_cast hb
  have habR : (a:ℝ)≤b := by exact_mod_cast hab
  unfold clippedOutput
  rcases le_total (c:ℝ) ((b:ℝ)-(a:ℝ)/2) with h | h
  · simpa only [min_eq_left h] using hc
  · rw [min_eq_right h]
    constructor
    · apply (le_div_iff₀ hbR).mpr
      nlinarith
    · exact (div_le_div_of_nonneg_right h hbR.le).trans hc.2

/-- Every genuinely supported finite path lies below the certified continuous
majorant, with a uniform explicit logarithmic remainder. -/
theorem pathEntropy_le_refined {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125))
    {p : ℕ×ℕ} (hp : p ∈ supportedPaths k w) :
    pathEntropy k w p ≤ ((k*24:ℕ):ℝ)*Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) +
      12*(Real.log (k*24:ℕ)+1)+1 := by
  obtain ⟨hmem,hcoef,ht₁,ht₂⟩ := mem_filter.mp hp
  obtain ⟨ha,hc⟩ := mem_product.mp hmem
  have ha' := mem_Icc.mp ha
  have hc' : p.2 ≤ k*24 := by simpa only [mem_range, Nat.lt_succ_iff] using hc
  have hb : 0 < k*24 := by omega
  have hbR : (0:ℝ)<(k*24:ℕ) := by exact_mod_cast hb
  have ha0 : 0 < p.1 := by omega
  have hc0 : 0 < p.2 := (accT_supported (by omega) ht₁).1
  obtain ⟨hlo₁,hhi₁⟩ := accT_even_feasible ha0 (golayCoefficient_even hcoef) ht₁
  obtain ⟨hlo₂,hhi₂,htc,hct⟩ := clippedOutput_spec hc0 hc' ht₂
  have h₁ := transitionEntropy_le_shift hb ha0 ha'.2 ht₁ hlo₁ hhi₁ (le_refl _) (by norm_num)
  have h₂ := transitionEntropy_le_shift hb hc0 hc' ht₂ hlo₂ hhi₂ htc hct
  have hm := Spin.Majorant.majorant_pointwise
    (a := (p.1:ℝ)/(k*24:ℕ)) (b := (p.2:ℝ)/(k*24:ℕ))
    ⟨by positivity, (div_le_one hbR).mpr (by exact_mod_cast ha'.2)⟩
    ⟨by positivity, (div_le_one hbR).mpr (by exact_mod_cast hc')⟩
    (clippedOutput_window hb hc' hw)
    ⟨(normalized_feasible hbR hlo₁ hhi₁).1, (normalized_feasible hbR hlo₁ hhi₁).2,
      (normalized_feasible hbR hlo₂ hhi₂).1, (normalized_feasible hbR hlo₂ hhi₂).2⟩
  have hm' := mul_le_mul_of_nonneg_left hm hbR.le
  have hs := refined_count_step hbR htc hct
  unfold pathEntropy
  linarith

/-- The exact finite maximum has a logarithmic remainder uniformly over the
closed certified output window. -/
theorem finiteEnvelope_le_refined {k w : ℕ} (hk : 0 < k)
    (hw : (w:ℝ)/(k*24:ℕ) ∈ Set.Icc ((13:ℝ)/125) (112/125)) :
    finiteEnvelope k w ≤ ((k*24:ℕ):ℝ)*Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) +
      12*(Real.log (k*24:ℕ)+1)+1 := by
  unfold finiteEnvelope
  split_ifs with h
  · exact Finset.sup'_le h _ (fun p hp => pathEntropy_le_refined hk hw hp)
  · have hm := refined_nonneg hw
    have hl : 0 ≤ Real.log ((k*24:ℕ):ℝ) := Real.log_nonneg (by exact_mod_cast (by omega : 1 ≤ k*24))
    have hp : 0 ≤ ((k*24:ℕ):ℝ)*Spin.Majorant.refined.toFun ((w:ℝ)/(k*24:ℕ)) := mul_nonneg (by positivity) hm
    linarith

end Spin.Structured.ConcreteOuter

