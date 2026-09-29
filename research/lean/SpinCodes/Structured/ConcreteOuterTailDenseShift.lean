import SpinCodes.Structured.ConcreteOuterEnvelopeShiftTransition
import SpinCodes.Numeric.BACentral

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

lemma xlogx_step_two {B x y : ℝ} (hB : 1 ≤ B)
    (hx : 0 ≤ x) (hxB : x ≤ B) (hy : 0 ≤ y) (hyB : y ≤ B)
    (hd : |x-y| ≤ 2) :
    |xlogx x-xlogx y| ≤ 2*(Real.log B+1) := by
  have h₁ := xlogx_step hB hx hxB
    (show 0 ≤ (x+y)/2 by linarith) (show (x+y)/2 ≤ B by linarith)
    (show |x-(x+y)/2| ≤ 1 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp hd).1,(abs_le.mp hd).2])
  have h₂ := xlogx_step hB
    (show 0 ≤ (x+y)/2 by linarith) (show (x+y)/2 ≤ B by linarith) hy hyB
    (show |(x+y)/2-y| ≤ 1 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp hd).1,(abs_le.mp hd).2])
  apply abs_le.mpr
  constructor <;> linarith [(abs_le.mp h₁).1,(abs_le.mp h₁).2,(abs_le.mp h₂).1,(abs_le.mp h₂).2]

lemma realLayerEntropy_step_two {B n k n' k' : ℝ} (hB : 1 ≤ B)
    (hk : 0 ≤ k) (hkn : k ≤ n) (hnB : n ≤ B)
    (hk' : 0 ≤ k') (hkn' : k' ≤ n') (hnB' : n' ≤ B)
    (hn : |n-n'| ≤ 2) (hkk : |k-k'| ≤ 2)
    (hc : |(n-k)-(n'-k')| ≤ 2) :
    |realLayerEntropy n k-realLayerEntropy n' k'| ≤ 6*(Real.log B+1) := by
  have h₁ := xlogx_step_two hB (hk.trans hkn) hnB (hk'.trans hkn') hnB' hn
  have h₂ := xlogx_step_two hB hk (hkn.trans hnB) hk' (hkn'.trans hnB') hkk
  have h₃ := xlogx_step_two hB (by linarith : 0 ≤ n-k) (by linarith)
    (by linarith : 0 ≤ n'-k') (by linarith) hc
  unfold realLayerEntropy
  apply abs_le.mpr
  constructor <;> linarith [(abs_le.mp h₁).1,(abs_le.mp h₁).2,
    (abs_le.mp h₂).1,(abs_le.mp h₂).2,(abs_le.mp h₃).1,(abs_le.mp h₃).2]

lemma piBA_count_step {B a c a' c' : ℝ} (hB : 1 ≤ B)
    (ha : 0 ≤ a) (hf₁ : a/2 ≤ c) (hf₂ : c ≤ B-a/2)
    (ha' : 0 ≤ a') (hf₁' : a'/2 ≤ c') (hf₂' : c' ≤ B-a'/2)
    (haa : |a-a'| ≤ 2) (hcc : |c-c'| ≤ 1) :
    |B*piBA (a/B) (c/B)-B*piBA (a'/B) (c'/B)| ≤ 18*(Real.log B+1) := by
  have hb : 0 < B := by linarith
  have h₁ := realLayerEntropy_step_two hB
    (show 0 ≤ a/2 by linarith) hf₁ (by linarith : c ≤ B)
    (show 0 ≤ a'/2 by linarith) hf₁' (by linarith : c' ≤ B)
    (show |c-c'| ≤ 2 by linarith)
    (show |a/2-a'/2| ≤ 2 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp haa).1,(abs_le.mp haa).2])
    (show |(c-a/2)-(c'-a'/2)| ≤ 2 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp haa).1,(abs_le.mp haa).2,(abs_le.mp hcc).1,(abs_le.mp hcc).2])
  have h₂ := realLayerEntropy_step_two hB
    (show 0 ≤ a/2 by linarith) (by linarith : a/2 ≤ B-c) (by linarith : B-c ≤ B)
    (show 0 ≤ a'/2 by linarith) (by linarith : a'/2 ≤ B-c') (by linarith : B-c' ≤ B)
    (show |(B-c)-(B-c')| ≤ 2 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp hcc).1,(abs_le.mp hcc).2])
    (show |a/2-a'/2| ≤ 2 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp haa).1,(abs_le.mp haa).2])
    (show |(B-c-a/2)-(B-c'-a'/2)| ≤ 2 by apply abs_le.mpr; constructor <;> linarith [(abs_le.mp haa).1,(abs_le.mp haa).2,(abs_le.mp hcc).1,(abs_le.mp hcc).2])
  have h₃ := realLayerEntropy_step_two hB ha (by linarith : a ≤ B) le_rfl
    ha' (by linarith : a' ≤ B) le_rfl (by simp) haa (by simpa [abs_sub_comm] using haa)
  rw [←realLayerEntropy_pi hb ha hf₁ hf₂, ←realLayerEntropy_pi hb ha' hf₁' hf₂']
  apply abs_le.mpr
  constructor <;> linarith [(abs_le.mp h₁).1,(abs_le.mp h₁).2,
    (abs_le.mp h₂).1,(abs_le.mp h₂).2,(abs_le.mp h₃).1,(abs_le.mp h₃).2]

lemma gBA_radial_lower {a a' : ℝ} (ha : 0 < a) (ha1 : a ≤ 1)
    (ha' : 0 ≤ a') (hle : a' ≤ a) : a'/a*gBA a ≤ gBA a' := by
  change a'/a*gBA a ≤ sInf ((fun u => gObj u a') '' Set.Ioi 0)
  refine le_csInf ⟨gObj 1 a', ⟨1, by norm_num, rfl⟩⟩ ?_
  rintro y ⟨u,hu,rfl⟩
  have hg := mul_le_mul_of_nonneg_left (gBA_le ha.le ha1 hu) (show 0 ≤ a'/a by positivity)
  have hz := mul_nonneg (show 0 ≤ 1-a'/a by rw [sub_nonneg, div_le_one ha]; exact hle)
    (gObj_nonneg hu (by norm_num : (0:ℝ) ≤ 0) (by norm_num : (0:ℝ) ≤ 1))
  have he : a'/a*gObj u a+(1-a'/a)*gObj u 0=gObj u a' := by
    unfold gObj
    field_simp
    ring
  linarith

/-- A one-sided modulus away from input density zero. -/
lemma gBA_count_step {B a a' : ℝ} (hB : 0 < B)
    (ha : B/125 ≤ a) (haB : a ≤ B) (ha' : 0 ≤ a') (hle : a' ≤ a)
    (hd : a-a' ≤ 2) : B*gBA (a/B) ≤ B*gBA (a'/B)+250 := by
  have ha0 : 0 < a := by linarith
  have hx0 : 0 < a/B := by positivity
  have hx1 : a/B ≤ 1 := (div_le_one hB).mpr haB
  have hy0 : 0 ≤ a'/B := by positivity
  have hyx : a'/B ≤ a/B := (div_le_div_iff_of_pos_right hB).mpr hle
  have hr := gBA_radial_lower hx0 hx1 hy0 hyx
  have hg := gBA_le_half_log_two ⟨hx0.le,hx1⟩
  have hlog := Real.log_le_sub_one_of_pos (by norm_num : (0:ℝ)<2)
  have hg1 : gBA (a/B) ≤ 1 := by linarith
  have hratio : (a'/B)/(a/B)=a'/a := by field_simp
  rw [hratio] at hr
  have hscale := mul_le_mul_of_nonneg_left hr (show 0 ≤ a by positivity)
  have he : a*(a'/a*gBA (a/B))=a'*gBA (a/B) := by field_simp
  rw [he] at hscale
  have hdiff := mul_le_mul_of_nonneg_left hg1 (show 0 ≤ a-a' by linarith)
  have hbound : a*(gBA (a/B)-gBA (a'/B)) ≤ 2 := by nlinarith
  have hnonneg := gBA_nonneg hx0.le hx1
  by_cases hh : gBA (a/B)-gBA (a'/B) ≤ 0
  · nlinarith
  · have := mul_le_mul_of_nonneg_right ha (le_of_lt (lt_of_not_ge hh))
    nlinarith

end Spin.Structured.ConcreteOuter


