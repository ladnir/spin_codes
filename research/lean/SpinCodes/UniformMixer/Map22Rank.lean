import SpinCodes.UniformMixer.Map22

/-!
# Full rank of the production 64-by-22 expansion

Twenty-two output coordinates form an invertible square submatrix.
The inverse masks below are a compact, kernel-checked witness. This
establishes injectivity of the actual expansion, independently of its
minimum nonzero output weight.
-/

noncomputable section

namespace Spin.UniformMixer

def pivotOutput : Fin 22 → Fin 64 := ![
  0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24,
  32, 33, 34, 36, 40, 48]

def pivotInverseMask : Fin 22 → ℕ := ![
  0x1, 0x3, 0x5, 0x11, 0x81, 0x801, 0x10001, 0x2dca00,
  0x1d4c40, 0x5c245, 0xc0291, 0xcc0d1, 0x94a85, 0x390891,
  0x3040c4, 0x54841, 0x210cca, 0x39ccf6, 0x54d53, 0x1418d2,
  0x38ec3, 0x142881]

def pivotInverseCoeff (i p : Fin 22) : ZMod 2 :=
  if (pivotInverseMask i).testBit p.val then 1 else 0

set_option maxRecDepth 8192 in
theorem pivot_inverse_certificate :
    ∀ i k : Fin 22,
      (∑ p : Fin 22, pivotInverseCoeff i p *
        coordCoeff (pivotOutput p) k) =
        if i = k then 1 else 0 := by
  decide

theorem pivot_reconstruct (q : State22) (i : Fin 22) :
    (∑ p : Fin 22, pivotInverseCoeff i p *
      coordMap (pivotOutput p) q) = q i := by
  calc
    _ = ∑ k : Fin 22,
        (∑ p : Fin 22, pivotInverseCoeff i p *
          coordCoeff (pivotOutput p) k) * q k := by
      simp only [coordMap, LinearMap.coe_mk, AddHom.coe_mk,
        Finset.mul_sum, mul_assoc, Finset.sum_mul]
      rw [Finset.sum_comm]
    _ = ∑ k : Fin 22, (if i = k then 1 else 0) * q k := by
      apply Finset.sum_congr rfl
      intro k _
      rw [pivot_inverse_certificate]
    _ = q i := by simp

theorem expansion_injective :
    Function.Injective (fun q : State22 =>
      fun j : Fin 64 => coordMap j q) := by
  intro q r h
  funext i
  rw [← pivot_reconstruct q i, ← pivot_reconstruct r i]
  apply Finset.sum_congr rfl
  intro p _
  have hp : coordMap (pivotOutput p) q =
      coordMap (pivotOutput p) r := congrFun h (pivotOutput p)
  rw [hp]

end Spin.UniformMixer
