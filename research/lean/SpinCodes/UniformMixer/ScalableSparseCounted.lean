import SpinCodes.UniformMixer.ScalableSparseFiniteBound

/-! Exact inverse-binomial cost for Scalable sparse occupation layers. -/

noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.SparseRate
  Spin.Structured.ConcreteMarked

theorem scalable_sparse_counted_bound {L Q b R d : ℕ}
    (hQ : 4096 ≤ Q) (hQL : Q ≤ L) (hb : 1000 ≤ b)
    {α ε rowCost : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000)
    (hα : (L : ℝ) * α = Q) (hsteps : 64 * R = L * b)
    (hd : (d : ℝ) ≤ (11 / 100) * (L : ℝ) * b)
    (hschedule : Real.log (L : ℝ) ≤
      (4 * Real.log 2 / 39) * (b : ℝ))
    (hε : ε ≤ 1 / 500)
    (hrow : rowCost ≤ Real.exp ((Q : ℝ) * b *
      (Real.log 2 / 2 + 1281 / 100000 + ε))) :
    (L.choose Q : ℝ) * (rowCost *
      (((8 * Real.sqrt (Q : ℝ)) ^ b *
          Real.exp (((L : ℝ) * b) *
            Spin.binKL α ((8 / 5) * α))) *
        (1024 * (1 - 48 * α) ^ R) /
          (1 - (8 / 5) * α) ^ d)) ≤
      Real.exp (-(3 / 500) * (Q : ℝ) * b) := by
  have hr : 0 ≤ 1 - 48 * α := by linarith
  have hz : 0 ≤ 1 - (8 / 5 : ℝ) * α := by linarith
  have he : (L.choose Q : ℝ) * (rowCost *
      (((8 * Real.sqrt (Q : ℝ)) ^ b *
          Real.exp (((L : ℝ) * b) *
            Spin.binKL α ((8 / 5) * α))) *
        (1024 * (1 - 48 * α) ^ R) /
          (1 - (8 / 5) * α) ^ d)) ≤
      scalableSparseFiniteBound L Q b R d α ε := by
    unfold scalableSparseFiniteBound
    have hh := mul_le_mul_of_nonneg_left hrow
      (by positivity : 0 ≤ (L.choose Q : ℝ) *
        (((8 * Real.sqrt (Q : ℝ)) ^ b *
          Real.exp (((L : ℝ) * b) *
            Spin.binKL α ((8 / 5) * α))) *
          (1024 * (1 - 48 * α) ^ R) /
            (1 - (8 / 5) * α) ^ d))
    convert hh using 1 <;> ring
  exact he.trans (scalableSparseFiniteBound_le
    hQ hQL hb hα0 hα1 hα hsteps hd hschedule hε)

theorem scalable_sparse_counted_mass {L Q b R d : ℕ}
    (hQ : 4096 ≤ Q) (hQL : Q ≤ L) (hb : 1000 ≤ b)
    {α ε rowCost : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000)
    (hα : (L : ℝ) * α = Q) (hsteps : 64 * R = L * b)
    (hd : (d : ℝ) ≤ (11 / 100) * (L : ℝ) * b)
    (hschedule : Real.log (L : ℝ) ≤
      (4 * Real.log 2 / 39) * (b : ℝ))
    (hε : ε ≤ 1 / 500) (hrow0 : 0 ≤ rowCost)
    (hrow : rowCost ≤ Real.exp ((Q : ℝ) * b *
      (Real.log 2 / 2 + 1281 / 100000 + ε))) :
    (L.choose Q : ℝ) * (rowCost *
      (((1 / markMass L Q ((8 / 5) * α)) ^ b *
        (1024 * (1 - 48 * α) ^ R)) /
          (1 - (8 / 5) * α) ^ d)) ≤
      Real.exp (-(3 / 500) * (Q : ℝ) * b) := by
  have hQ0 : 0 < Q := by omega
  have hL : (0 : ℝ) < L := by
    exact_mod_cast lt_of_lt_of_le hQ0 hQL
  have hαeq : (Q : ℝ) / L = α :=
    (div_eq_iff hL.ne').mpr (by nlinarith [hα])
  have hi := markMass_inverse_pow_le hQ0 hQL
    (by positivity : 0 < (8 / 5 : ℝ) * α)
    (by linarith : (8 / 5 : ℝ) * α < 1) b
  rw [hαeq] at hi
  have hr : 0 ≤ 1 - 48 * α := by linarith
  have hz : 0 ≤ 1 - (8 / 5 : ℝ) * α := by linarith
  refine le_trans ?_ (scalable_sparse_counted_bound
    hQ hQL hb hα0 hα1 hα hsteps hd hschedule hε hrow)
  apply mul_le_mul_of_nonneg_left _ (Nat.cast_nonneg _)
  apply mul_le_mul_of_nonneg_left _ hrow0
  apply div_le_div_of_nonneg_right _ (pow_nonneg hz d)
  exact mul_le_mul_of_nonneg_right hi (by positivity)

#print axioms scalable_sparse_counted_mass
end Spin.UniformMixer
