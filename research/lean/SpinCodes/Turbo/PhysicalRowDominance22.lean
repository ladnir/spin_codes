import SpinCodes.Turbo.PhysicalSymbolicRow22

/-!
# Entrywise domination of a complete symbolic physical row

An outward numeric matrix may replace the symbolic operator only after
every entry of the complete chosen row is checked. This theorem does
not permit mixing entries from alternative birth-row decompositions.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

theorem physicalPiece22_mono {i : Fin 11} {b B : ℝ}
    {μ : State22 → ℝ} (h : physicalPiece22 i b μ)
    (hbound : b ≤ B) : physicalPiece22 i B μ := by
  refine ⟨h.1, h.2.1, h.2.2.1.trans hbound, ?_⟩
  intro hi q
  exact (h.2.2.2 hi q).trans
    (div_le_div_of_nonneg_right hbound (by norm_num))

theorem physicalRowSound22_mono
    {K : Matrix State22 State22 ℝ}
    {B C : Matrix (Fin 11) (Fin 11) ℝ}
    (hB : physicalRowSound22 K B)
    (hBC : ∀ i k, B i k ≤ C i k) :
    physicalRowSound22 K C := by
  intro i b μ hμ
  obtain ⟨parts, hsum, hparts⟩ := hB i b μ hμ
  refine ⟨parts, hsum, ?_⟩
  intro k
  have hb : 0 ≤ b :=
    (Finset.sum_nonneg (fun q _ => hμ.1 q)).trans hμ.2.2.1
  exact physicalPiece22_mono (hparts k)
    (mul_le_mul_of_nonneg_left (hBC i k) hb)

theorem turboSelectedKernel22_outward_row_sound
    (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (B : Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ i k, turboSymbolicRow22 j z i k ≤ B i k) :
    physicalRowSound22
      (Matrix.of (turboSelectedKernel22 j z)) B :=
  physicalRowSound22_mono (turboSymbolicRow22_sound j z hz) hB

#print axioms physicalPiece22_mono
#print axioms physicalRowSound22_mono
#print axioms turboSelectedKernel22_outward_row_sound

end Spin.Turbo
