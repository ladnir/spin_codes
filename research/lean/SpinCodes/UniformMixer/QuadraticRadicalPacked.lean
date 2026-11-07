import SpinCodes.UniformMixer.QuadraticRadicalBasis
import SpinCodes.UniformMixer.Map22PackedBridge

/-!
# Packed quadratic coefficients

The finite radical certificate ranges over the fifteen quadratic
coefficients. This module proves that fifteen-bit integers enumerate
those coefficient vectors exactly once, and that the executable six-basis
test measures their mathematical radicals.
-/

namespace Spin.UniformMixer

def packedQuadratic (q : Fin (2 ^ 15)) : Fin 15 → ZMod 2 :=
  fun i => bitScalar (q.val.testBit i.val)

theorem packedQuadratic_injective : Function.Injective packedQuadratic := by
  intro q r h
  apply Fin.ext
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < 15
  · exact bitScalar_injective (congrFun h ⟨i, hi⟩)
  · have hpow : 2 ^ 15 ≤ 2 ^ i :=
      Nat.pow_le_pow_right (by decide) (by omega)
    rw [Nat.testBit_lt_two_pow (q.isLt.trans_le hpow),
      Nat.testBit_lt_two_pow (r.isLt.trans_le hpow)]

theorem packedQuadratic_bijective : Function.Bijective packedQuadratic := by
  have hcard : Fintype.card (Fin (2 ^ 15)) =
      Fintype.card (Fin 15 → ZMod 2) := by
    simp [ZMod.card]
  exact (Fintype.bijective_iff_injective_and_card _).mpr
    ⟨packedQuadratic_injective, hcard⟩

noncomputable def packedQuadraticEquiv : Fin (2 ^ 15) ≃ (Fin 15 → ZMod 2) :=
  Equiv.ofBijective packedQuadratic packedQuadratic_bijective

def packedRadicalSize (q : Fin (2 ^ 15)) : Nat :=
  (basisRadical (packedQuadratic q)).card

def packedRadicalSizeNat (q : Nat) : Nat :=
  (basisRadical (fun i : Fin 15 => bitScalar (q.testBit i.val))).card

theorem packedRadicalSize_eq_nat (q : Fin (2 ^ 15)) :
    packedRadicalSize q = packedRadicalSizeNat q.val := rfl

def radicalBlockCount (start length radicalSize : Nat) : Nat :=
  ((List.range length).filter (fun j =>
    packedRadicalSizeNat (start + j) = radicalSize)).length

theorem packedRadicalSize_eq (q : Fin (2 ^ 15)) :
    packedRadicalSize q = (quadraticRadical (packedQuadratic q)).card := by
  rw [packedRadicalSize, basisRadical_eq_quadraticRadical]

def packedRadicalCount (radicalSize : Nat) : Nat :=
  (Finset.univ.filter (fun q : Fin (2 ^ 15) =>
    packedRadicalSize q = radicalSize)).card

theorem packedRadicalCount_eq_actual (radicalSize : Nat) :
    packedRadicalCount radicalSize =
      (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
        (quadraticRadical c).card = radicalSize)).card := by
  classical
  unfold packedRadicalCount
  rw [Finset.card_filter, Finset.card_filter]
  apply Fintype.sum_equiv packedQuadraticEquiv
  intro q
  change (if packedRadicalSize q = radicalSize then 1 else 0) =
    if (quadraticRadical (packedQuadratic q)).card = radicalSize then 1 else 0
  rw [packedRadicalSize_eq]

end Spin.UniformMixer
