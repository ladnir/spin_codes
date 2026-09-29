import SpinCodes.Structured.SparseFiniteExponent

noncomputable section
namespace Spin.Structured.SparseRate

def finiteBound (L Q b R d : ℕ) (α ε : ℝ) : ℝ :=
  2048 * (L.choose Q:ℝ) * Real.exp ((Q:ℝ)*b*(Real.log 2/2+1281/100000+ε)) *
    (8*Real.sqrt (Q:ℝ))^b * Real.exp (((L:ℝ)*b)*Spin.binKL α ((8/5)*α)) *
      (1-96*α)^R / (1-(8/5)*α)^d

theorem finiteBound_pos {L Q b R d : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    0 < finiteBound L Q b R d α ε := by
  have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
  have hQr : (0:ℝ) < Q := by exact_mod_cast hQ
  have hr : 0 < 1-96*α := by linarith
  have hz : 0 < 1-(8/5:ℝ)*α := by linarith
  unfold finiteBound
  positivity

theorem finiteBound_log {L Q b R d : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    Real.log (finiteBound L Q b R d α ε) =
      Real.log 2048 + Real.log (L.choose Q:ℝ) + (Q:ℝ)*b*(Real.log 2/2+1281/100000+ε) +
        (b:ℝ)*Real.log (8*Real.sqrt (Q:ℝ)) + ((L:ℝ)*b)*Spin.binKL α ((8/5)*α) +
        (R:ℝ)*Real.log (1-96*α) - (d:ℝ)*Real.log (1-(8/5)*α) := by
  have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
  have hQr : (0:ℝ) < Q := by exact_mod_cast hQ
  have hr : 0 < 1-96*α := by linarith
  have hz : 0 < 1-(8/5:ℝ)*α := by linarith
  unfold finiteBound
  rw [Real.log_div (by positivity) (by positivity)]
  repeat' rw [Real.log_mul (by positivity) (by positivity)]
  simp only [Real.log_pow, Real.log_exp]
  rw [Real.log_mul (by norm_num : (8:ℝ) ≠ 0) (by positivity : Real.sqrt (Q:ℝ) ≠ 0)]

theorem finiteBound_log_le_exponent {L Q b R d : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000) :
    Real.log (finiteBound L Q b R d α ε) ≤ exponent L Q b R d α ε := by
  have hc : (0:ℝ) < (L.choose Q:ℝ) := by exact_mod_cast Nat.choose_pos hQL
  have hn : (L.choose Q:ℝ) ≤ (L:ℝ)^Q := by exact_mod_cast Nat.choose_le_pow L Q
  have hh := Real.log_le_log hc hn
  rw [Real.log_pow] at hh
  rw [finiteBound_log hQ hQL hα0 hα1]
  unfold exponent
  have hQnn : (0:ℝ) ≤ Q := by positivity
  linarith

/-- Finite sparse estimate with explicit b and ε thresholds. -/
theorem finiteBound_le {L Q b R d : ℕ} (hQ : 4096 ≤ Q) (hQL : Q ≤ L)
    (hb : 1000 ≤ b) {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000)
    (hα : (L:ℝ)*α = Q) (hrounds : 128*R = L*b)
    (hd : (d:ℝ) ≤ (11/100)*(L:ℝ)*b)
    (hschedule : Real.log (L:ℝ) ≤ (4*Real.log 2/39)*(b:ℝ)) (hε : ε ≤ 1/500) :
    finiteBound L Q b R d α ε ≤ Real.exp (-(3/500)*(Q:ℝ)*b) := by
  have hQ0 : 0 < Q := by omega
  have hL : 0 < L := lt_of_lt_of_le hQ0 hQL
  have he := exponent_le (by exact_mod_cast hL) (by exact_mod_cast hQ)
    (by exact_mod_cast hb) (Nat.cast_nonneg R) (Nat.cast_nonneg d) hα0 hα1 hα
    (by exact_mod_cast hrounds) hd hschedule hε
  have hh := (finiteBound_log_le_exponent hQ0 hQL hα0 hα1).trans he
  have hh2 := Real.exp_le_exp.mpr hh
  rwa [Real.exp_log (finiteBound_pos hQ0 hQL hα0 hα1)] at hh2

end Spin.Structured.SparseRate

