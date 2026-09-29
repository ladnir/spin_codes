import SpinCodes.Structured.ConcreteOuterEnvelopeShiftEntropy

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

lemma accT_supported {b a c : ℕ} (ha : a ≠ 0) (ht : accT b a c ≠ 0) :
    0 < c ∧ c ≤ b ∧ (a+1)/2 ≤ c ∧ a-(a+1)/2 ≤ b-c := by
  have hc0 : c ≠ 0 := by intro h; subst c; exact ht (accT_zero_right b ha)
  have hcb : c ≤ b := by
    by_contra h
    exact ht (accT_of_lt b a c ha (by omega))
  have he : accT b a c = (c-1).choose ((a+1)/2-1) *
      (b-c).choose (a-(a+1)/2) := by
    simp only [accT, if_neg ha, if_neg (show ¬(c=0 ∨ b<c) by omega)]
  have h₁ : (a+1)/2-1 ≤ c-1 := by
    by_contra h
    exact ht (by rw [he, Nat.choose_eq_zero_of_lt (by omega), zero_mul])
  have h₂ : a-(a+1)/2 ≤ b-c := by
    by_contra h
    exact ht (by rw [he, Nat.choose_eq_zero_of_lt (n:=b-c) (k:=a-(a+1)/2) (by omega), mul_zero])
  omega

lemma realLayerEntropy_pi {b a c : ℝ} (hb : 0 < b) (ha : 0 ≤ a)
    (hlo : a/2 ≤ c) (hhi : c ≤ b-a/2) :
    realLayerEntropy c (a/2) + realLayerEntropy (b-c) (a/2) - realLayerEntropy b a =
      b * piBA (a/b) (c/b) := by
  have hab : a ≤ b := by linarith
  rw [realLayerEntropy_eq (by positivity) hlo,
    realLayerEntropy_eq (by positivity) (by linarith), realLayerEntropy_eq ha hab]
  unfold piBA
  have he : 1-c/b = (b-c)/b := by field_simp
  rw [he]
  have h₁ : a/b/(2*(c/b)) = (a/2)/c := by field_simp
  have h₂ : a/b/(2*((b-c)/b)) = (a/2)/(b-c) := by field_simp
  rw [h₁, h₂]
  field_simp

/-- Exact parity arguments differ from feasible real counts by at most one.
The resulting loss is uniform even on either accumulator boundary. -/
theorem transitionEntropy_le_shift {b a c : ℕ} (hb : 0 < b) (ha : 0 < a)
    (hab : a ≤ b) (ht : accT b a c ≠ 0) {t : ℝ}
    (hlo : (a:ℝ)/2 ≤ t) (hhi : t ≤ (b:ℝ)-(a:ℝ)/2)
    (htc : t ≤ c) (hct : (c:ℝ)-t ≤ 1/2) :
    transitionEntropy b a c ≤ (b:ℝ) * piBA ((a:ℝ)/b) (t/b) +
      6 * (Real.log b + 1) := by
  obtain ⟨hc, hcb, hr, hs⟩ := accT_supported (by omega) ht
  have hr1 : 1 ≤ (a+1)/2 := by omega
  have hra : (a+1)/2 ≤ a := by omega
  have hrlo : a ≤ 2*((a+1)/2) := by omega
  have hrhi : 2*((a+1)/2) ≤ a+1 := by omega
  have hrloR : (a:ℝ) ≤ 2*(((a+1)/2:ℕ):ℝ) := by exact_mod_cast hrlo
  have hrhiR : 2*(((a+1)/2:ℕ):ℝ) ≤ (a:ℝ)+1 := by exact_mod_cast hrhi
  have hcbR : (c:ℝ) ≤ b := by exact_mod_cast hcb
  have habR : (a:ℝ) ≤ b := by exact_mod_cast hab
  have hcR : (1:ℝ) ≤ c := by exact_mod_cast hc
  have hrR : (((a+1)/2:ℕ):ℝ) ≤ c := by exact_mod_cast hr
  have haR : (0:ℝ) ≤ a := by positivity
  have hbR : (1:ℝ) ≤ b := by exact_mod_cast hb
  have h₁ := realLayerEntropy_step hbR
    (show (0:ℝ) ≤ ((a+1)/2-1:ℕ) by positivity)
    (show (((a+1)/2-1:ℕ):ℝ) ≤ (c-1:ℕ) by exact_mod_cast (show (a+1)/2-1 ≤ c-1 by omega))
    (show ((c-1:ℕ):ℝ) ≤ b by exact_mod_cast (show c-1 ≤ b by omega))
    (show (0:ℝ) ≤ (a:ℝ)/2 by positivity) hlo (by linarith : t ≤ b)
  have h₂ := realLayerEntropy_step hbR
    (show (0:ℝ) ≤ (a-(a+1)/2:ℕ) by positivity)
    (show ((a-(a+1)/2:ℕ):ℝ) ≤ (b-c:ℕ) by exact_mod_cast hs)
    (show ((b-c:ℕ):ℝ) ≤ b by exact_mod_cast Nat.sub_le b c)
    (show (0:ℝ) ≤ (a:ℝ)/2 by positivity) (by linarith : (a:ℝ)/2 ≤ (b:ℝ)-t)
    (by linarith : (b:ℝ)-t ≤ b)
  rw [Nat.cast_sub (by omega : 1 ≤ c), Nat.cast_one,
    Nat.cast_sub hr1, Nat.cast_one] at h₁
  rw [Nat.cast_sub hcb, Nat.cast_sub hra] at h₂
  have h₁' := h₁ (abs_le.mpr ⟨by linarith, by linarith⟩)
    (abs_le.mpr ⟨by linarith, by linarith⟩) (abs_le.mpr ⟨by linarith, by linarith⟩)
  have h₂' := h₂ (abs_le.mpr ⟨by linarith, by linarith⟩)
    (abs_le.mpr ⟨by linarith, by linarith⟩) (abs_le.mpr ⟨by linarith, by linarith⟩)
  rw [transitionEntropy, if_neg (by omega : a ≠ 0),
    layerEntropy_eq_real _ _ (by omega), layerEntropy_eq_real _ _ hs,
    layerEntropy_eq_real _ _ hab]
  rw [Nat.cast_sub (by omega : 1 ≤ c), Nat.cast_one,
    Nat.cast_sub hr1, Nat.cast_one, Nat.cast_sub hcb, Nat.cast_sub hra]
  have he := realLayerEntropy_pi (by exact_mod_cast hb : (0:ℝ)<b) haR hlo hhi
  linarith [(abs_le.mp h₁').2, (abs_le.mp h₂').2]

end Spin.Structured.ConcreteOuter
