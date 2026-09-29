import SpinCodes.Structured.SparseConditioningMassMode
import SpinCodes.Structured.SparseConditioningMassLikelihood

noncomputable section
namespace Spin.Structured.ConcreteMarked

/-- Sharp conditioning probability, with the paper's binary relative entropy. -/
theorem markMass_lower {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    Real.exp (-(L:ℝ)*Spin.binKL ((Q:ℝ)/L) p) / (8*Real.sqrt (Q:ℝ)) ≤
      markMass L Q p := by
  rw [markMass_likelihood_ratio hQ hQL hp0 hp1]
  have hh := mul_le_mul_of_nonneg_right (markMass_mode_lower hQ hQL)
    (Real.exp_pos (-(L:ℝ)*Spin.binKL ((Q:ℝ)/L) p)).le
  simpa only [one_div_mul_eq_div, mul_comm] using hh

theorem markMass_inverse_le {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    1 / markMass L Q p ≤ 8*Real.sqrt (Q:ℝ) *
      Real.exp ((L:ℝ)*Spin.binKL ((Q:ℝ)/L) p) := by
  have hs : 0 < 8*Real.sqrt (Q:ℝ) := by
    have hQr : (0:ℝ) < Q := by exact_mod_cast hQ
    positivity
  have hm : 0 < markMass L Q p := by
    unfold markMass
    have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
    have hpc : 0 < 1-p := by linarith
    positivity
  have hh := (div_le_iff₀ hs).mp (markMass_lower hQ hQL hp0 hp1)
  have hmul := mul_le_mul_of_nonneg_right hh
    (Real.exp_pos ((L:ℝ)*Spin.binKL ((Q:ℝ)/L) p)).le
  have he : Real.exp (-(L:ℝ)*Spin.binKL ((Q:ℝ)/L) p) *
      Real.exp ((L:ℝ)*Spin.binKL ((Q:ℝ)/L) p) = 1 := by
    rw [← Real.exp_add]
    convert Real.exp_zero using 1 <;> ring
  rw [he] at hmul
  apply (div_le_iff₀ hm).mpr
  nlinarith

theorem markMass_inverse_pow_le {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (b : ℕ) :
    (1 / markMass L Q p)^b ≤ (8*Real.sqrt (Q:ℝ))^b *
      Real.exp (((L:ℝ)*b)*Spin.binKL ((Q:ℝ)/L) p) := by
  have hm : 0 ≤ markMass L Q p := by unfold markMass; positivity
  have hh := pow_le_pow_left₀ (div_nonneg zero_le_one hm)
    (markMass_inverse_le hQ hQL hp0 hp1) b
  rw [mul_pow, ← Real.exp_nat_mul] at hh
  convert hh using 1 <;> ring

end Spin.Structured.ConcreteMarked
