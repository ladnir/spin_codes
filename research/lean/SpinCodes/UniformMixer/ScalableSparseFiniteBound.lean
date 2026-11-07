import SpinCodes.UniformMixer.ScalableSparseExponent
import SpinCodes.Structured.SparseConditioningMass

/-! Convert the conditioned sparse count to the certified exponent. -/

noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.SparseRate
  Spin.Structured.ConcreteMarked

def scalableSparseFiniteBound (L Q b R d : ℕ) (α ε : ℝ) : ℝ :=
  1024 * (L.choose Q : ℝ) *
    Real.exp ((Q : ℝ) * b *
      (Real.log 2 / 2 + 1281 / 100000 + ε)) *
    (8 * Real.sqrt (Q : ℝ)) ^ b *
    Real.exp (((L : ℝ) * b) * Spin.binKL α ((8 / 5) * α)) *
    (1 - 48 * α) ^ R / (1 - (8 / 5) * α) ^ d

theorem scalableSparseFiniteBound_pos {L Q b R d : ℕ}
    (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    0 < scalableSparseFiniteBound L Q b R d α ε := by
  have hc : (0 : ℝ) < (L.choose Q : ℝ) := by
    exact_mod_cast Nat.choose_pos hQL
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have hr : 0 < 1 - 48 * α := by linarith
  have hz : 0 < 1 - (8 / 5 : ℝ) * α := by linarith
  unfold scalableSparseFiniteBound
  positivity

theorem scalableSparseFiniteBound_log {L Q b R d : ℕ}
    (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    Real.log (scalableSparseFiniteBound L Q b R d α ε) =
      Real.log 1024 + Real.log (L.choose Q : ℝ) +
        (Q : ℝ) * b *
          (Real.log 2 / 2 + 1281 / 100000 + ε) +
        (b : ℝ) * Real.log (8 * Real.sqrt (Q : ℝ)) +
        ((L : ℝ) * b) * Spin.binKL α ((8 / 5) * α) +
        (R : ℝ) * Real.log (1 - 48 * α) -
        (d : ℝ) * Real.log (1 - (8 / 5) * α) := by
  have hc : (0 : ℝ) < (L.choose Q : ℝ) := by
    exact_mod_cast Nat.choose_pos hQL
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  have hr : 0 < 1 - 48 * α := by linarith
  have hz : 0 < 1 - (8 / 5 : ℝ) * α := by linarith
  unfold scalableSparseFiniteBound
  rw [Real.log_div (by positivity) (by positivity)]
  repeat' rw [Real.log_mul (by positivity) (by positivity)]
  simp only [Real.log_pow, Real.log_exp]
  rw [Real.log_mul (by norm_num : (8 : ℝ) ≠ 0)
    (by positivity : Real.sqrt (Q : ℝ) ≠ 0)]

theorem scalableSparseFiniteBound_log_le {L Q b R d : ℕ}
    (hQ : 0 < Q) (hQL : Q ≤ L)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    Real.log (scalableSparseFiniteBound L Q b R d α ε) ≤
      scalableSparseExponent L Q b R d α ε := by
  have hc : (0 : ℝ) < (L.choose Q : ℝ) := by
    exact_mod_cast Nat.choose_pos hQL
  have hn : (L.choose Q : ℝ) ≤ (L : ℝ) ^ Q := by
    exact_mod_cast Nat.choose_le_pow L Q
  have hh := Real.log_le_log hc hn
  rw [Real.log_pow] at hh
  rw [scalableSparseFiniteBound_log hQ hQL hα0 hα1]
  unfold scalableSparseExponent
  have hQnn : (0 : ℝ) ≤ Q := by positivity
  linarith

theorem scalableSparseFiniteBound_le {L Q b R d : ℕ}
    (hQ : 4096 ≤ Q) (hQL : Q ≤ L) (hb : 1000 ≤ b)
    {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000)
    (hα : (L : ℝ) * α = Q) (hsteps : 64 * R = L * b)
    (hd : (d : ℝ) ≤ (11 / 100) * (L : ℝ) * b)
    (hschedule : Real.log (L : ℝ) ≤
      (4 * Real.log 2 / 39) * (b : ℝ))
    (hε : ε ≤ 1 / 500) :
    scalableSparseFiniteBound L Q b R d α ε ≤
      Real.exp (-(3 / 500) * (Q : ℝ) * b) := by
  have hQ0 : 0 < Q := by omega
  have hL : 0 < L := lt_of_lt_of_le hQ0 hQL
  have he := scalableSparseExponent_le
    (by exact_mod_cast hL) (by exact_mod_cast hQ)
    (by exact_mod_cast hb) (Nat.cast_nonneg R)
    (Nat.cast_nonneg d) hα0 hα1 hα
    (by exact_mod_cast hsteps) hd hschedule hε
  have hh := (scalableSparseFiniteBound_log_le
    hQ0 hQL hα0 hα1).trans he
  have hh2 := Real.exp_le_exp.mpr hh
  rwa [Real.exp_log
    (scalableSparseFiniteBound_pos hQ0 hQL hα0 hα1)] at hh2

#print axioms scalableSparseFiniteBound_le
end Spin.UniformMixer
