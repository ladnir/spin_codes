import SpinCodes.Structured.ConcreteMarkedConditioning
import SpinCodes.Structured.KLChain

noncomputable section
namespace Spin.Structured.ConcreteMarked

theorem markMass_log {L Q : ℕ} (hQL : Q ≤ L) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    Real.log (markMass L Q p) = Real.log (L.choose Q : ℝ) +
      (Q:ℝ)*Real.log p + ((L-Q:ℕ):ℝ)*Real.log (1-p) := by
  have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
  have hpc : 0 < 1-p := by linarith
  unfold markMass
  rw [Real.log_mul (by positivity) (by positivity),
    Real.log_mul (by positivity) (by positivity), Real.log_pow, Real.log_pow]

/-- Exact binomial likelihood ratio; the endpoint Q=L is included. -/
theorem markMass_likelihood_ratio {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    markMass L Q p = markMass L Q ((Q:ℝ)/L) *
      Real.exp (-(L:ℝ)*Spin.binKL ((Q:ℝ)/L) p) := by
  have hL : 0 < L := lt_of_lt_of_le hQ hQL
  have hLp : (0:ℝ) < L := by exact_mod_cast hL
  have hpc : 0 < 1-p := by linarith
  rcases eq_or_lt_of_le hQL with he | hlt
  · subst Q
    have hkl : Spin.binKL 1 p = -Real.log p := by
      simp [Spin.binKL, Real.log_div one_ne_zero hp0.ne']
    simp only [div_self hLp.ne', markMass, Nat.choose_self, Nat.cast_one,
      one_mul, Nat.sub_self, pow_zero, mul_one, one_pow, hkl]
    rw [show -(L:ℝ)*(-Real.log p) = (L:ℝ)*Real.log p by ring,
      Real.exp_nat_mul, Real.exp_log hp0]
  · have hα0 : (0:ℝ) < (Q:ℝ)/L := by positivity
    have hα1 : (Q:ℝ)/L < 1 := (div_lt_one hLp).mpr (by exact_mod_cast hlt)
    have hαc : 0 < 1-(Q:ℝ)/L := by linarith
    have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
    have hmp : 0 < markMass L Q p := by unfold markMass; positivity
    have hma : 0 < markMass L Q ((Q:ℝ)/L) := by unfold markMass; positivity
    have he : Real.log (markMass L Q p) =
        Real.log (markMass L Q ((Q:ℝ)/L)) - (L:ℝ)*Spin.binKL ((Q:ℝ)/L) p := by
      rw [markMass_log hQL hp0 hp1, markMass_log hQL hα0 hα1]
      unfold Spin.binKL
      rw [Real.log_div hα0.ne' hp0.ne', Real.log_div hαc.ne' hpc.ne', Nat.cast_sub hQL]
      field_simp
      ring
    rw [← Real.exp_log hmp, he, sub_eq_add_neg, Real.exp_add, Real.exp_log hma]
    congr 1
    congr 1
    ring

end Spin.Structured.ConcreteMarked
