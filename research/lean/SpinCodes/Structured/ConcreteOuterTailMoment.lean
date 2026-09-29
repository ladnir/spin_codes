import SpinCodes.Structured.ConcreteOuterTailTransition

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

/-- A finite negative-binomial generating-function bound. -/
theorem choose_geometric_sum_le {r : ℝ} (hr0 : 0 ≤ r) (hr1 : r < 1) (b k : ℕ) :
    (∑ i ∈ range b, (i.choose k:ℝ)*r^(i+1)) ≤ (r/(1-r))^(k+1) := by
  have hr : ‖r‖<1 := by simpa only [Real.norm_eq_abs, abs_of_nonneg hr0] using hr1
  have hs := hasSum_choose_mul_geometric_of_norm_lt_one k hr
  have hsum : (∑ n ∈ range b, ((n+k).choose k:ℝ)*r^n) ≤ 1/(1-r)^(k+1) := by
    rw [←hs.tsum_eq]
    exact hs.summable.sum_le_tsum (range b) (fun _ _ => by positivity)
  have hzero : (∑ i ∈ range k, (i.choose k:ℝ)*r^(i+1))=0 := by
    apply Finset.sum_eq_zero
    intro i hi
    rw [Nat.choose_eq_zero_of_lt (mem_range.mp hi), Nat.cast_zero, zero_mul]
  calc
    (∑ i ∈ range b, (i.choose k:ℝ)*r^(i+1)) ≤
        ∑ i ∈ range (k+b), (i.choose k:ℝ)*r^(i+1) := by
      exact sum_le_sum_of_subset_of_nonneg (range_mono (by omega)) (fun _ _ _ => by positivity)
    _ = r^(k+1)*(∑ n ∈ range b, ((n+k).choose k:ℝ)*r^n) := by
      rw [sum_range_add, hzero, zero_add, mul_sum]
      apply sum_congr rfl
      intro n hn
      rw [show k+n+1=(k+1)+n by omega, pow_add, Nat.add_comm k n]
      ring
    _ ≤ r^(k+1)*(1/(1-r)^(k+1)) := mul_le_mul_of_nonneg_left hsum (by positivity)
    _ = _ := by rw [div_pow]; ring

/-- Exact first-accumulator geometric moment, bounded uniformly over length. -/
theorem even_accumulator_moment_le {b ℓ : ℕ} (hℓ : 0 < ℓ) (hb : 2*ℓ ≤ b) :
    (∑ h ∈ range (b+1), Pt b (2*ℓ) h*zeta^h) ≤
      (kappa*(ℓ:ℝ)/(b+1-2*ℓ:ℕ))^ℓ := by
  have hc : (0:ℝ)<b.choose (2*ℓ) := by exact_mod_cast Nat.choose_pos hb
  have hd : (0:ℝ)<(b+1-2*ℓ:ℕ) := by exact_mod_cast (by omega : 0<b+1-2*ℓ)
  have hratio : (b.choose ℓ:ℝ)/(b.choose (2*ℓ):ℝ) ≤
      ((4:ℝ)*(ℓ:ℝ)/(b+1-2*ℓ:ℕ))^ℓ := by
    have hh : (b.choose ℓ:ℝ)*((b+1-2*ℓ:ℕ):ℝ)^ℓ ≤
        (4:ℝ)^ℓ*(ℓ:ℝ)^ℓ*(b.choose (2*ℓ):ℝ) := by exact_mod_cast choose_ratio_le hℓ hb
    rw [div_pow, mul_pow]
    apply (div_le_iff₀ hc).mpr
    rw [div_mul_eq_mul_div]
    exact (le_div_iff₀ (pow_pos hd _)).mpr hh
  have hterm (i : ℕ) (hi : i∈range b) :
      Pt b (2*ℓ) (i+1)*zeta^(i+1) ≤
        ((b.choose ℓ:ℝ)/(b.choose (2*ℓ):ℝ))*((i.choose (ℓ-1):ℝ)*zeta^(i+1)) := by
    have hiB : i+1 ≤ b := by have := mem_range.mp hi; omega
    rw [Pt, accT_even hℓ (by omega) hiB]
    have he : i+1-1=i := by omega
    rw [he, Nat.cast_mul]
    have hch : ((b-(i+1)).choose ℓ:ℝ) ≤ b.choose ℓ := by
      exact_mod_cast Nat.choose_le_choose ℓ (Nat.sub_le b (i+1))
    have hh := mul_le_mul_of_nonneg_left hch (show (0:ℝ) ≤ (i.choose (ℓ-1):ℝ) by positivity)
    have hh' := mul_le_mul_of_nonneg_right (div_le_div_of_nonneg_right hh hc.le)
      (show 0≤zeta^(i+1) by exact pow_nonneg zeta_pos.le _)
    convert hh' using 1 <;> ring
  have hsum := sum_le_sum hterm
  rw [←mul_sum] at hsum
  have hgeom := choose_geometric_sum_le zeta_pos.le zeta_lt_one b (ℓ-1)
  rw [show ℓ-1+1=ℓ by omega] at hgeom
  have hz := zeta_pos
  have h1z := one_sub_zeta_pos
  have hfinal := (hsum.trans (mul_le_mul_of_nonneg_left hgeom (by positivity))).trans
    (mul_le_mul_of_nonneg_right hratio (by positivity : 0≤(zeta/(1-zeta))^ℓ))
  have he : ((4:ℝ)*(ℓ:ℝ)/(b+1-2*ℓ:ℕ))^ℓ*(zeta/(1-zeta))^ℓ =
      (kappa*(ℓ:ℝ)/(b+1-2*ℓ:ℕ))^ℓ := by
    rw [←mul_pow, kappa]
    congr 1
    ring
  rw [he] at hfinal
  rw [sum_range_succ', Pt, accT_zero_right b (by omega), Nat.cast_zero, zero_div, zero_mul, add_zero]
  exact hfinal

end Spin.Structured.ConcreteOuter

