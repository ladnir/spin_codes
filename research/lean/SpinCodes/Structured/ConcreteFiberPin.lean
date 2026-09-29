import SpinCodes.Structured.ConcreteFiberBounds

/-! Statement pins for the actual feedback map and its counting bounds. -/

noncomputable section
namespace Spin.Structured.ConcreteFiberPin
open ConcreteMaps

example (q : Finset (Fin 19)) (x : Finset (Fin 128)) :
    Spin.chi (CtransposeSet q) x = Spin.chi q (Cset x) := C_character_pairing q x

example : dual.card = 2 ^ 19 := dual_card

example (x : Finset (Fin 128)) : Spin.OrthTo dual x ↔ Cset x = ∅ := orthTo_iff_kernel x

example {x : Finset (Fin 128)} (hx : Cset x = ∅) (hne : x ≠ ∅) :
    4 ≤ x.card := kernel_distance hx hne

example (j : ℕ) : (2 ^ 19 : ℤ) * ((ConcreteMaps.kernelLayer j).card : ℤ) =
    ∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) * Spin.krawtchouk 128 j w :=
  kernel_fourier_weights j

example (j : ℕ) : (2 ^ 19 : ℤ) * ((equalSyndromePairs j).card : ℤ) =
    ∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) * (Spin.krawtchouk 128 j w) ^ 2 :=
  feedback_parseval_weights j

example (j : ℕ) : ∑ q : Finset (Fin 19), (syndromeFiber j q).card ^ 2 =
    (equalSyndromePairs j).card := sum_syndromeFiber_sq j

example (j : ℕ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (2 ^ 19 - 1) * (syndromeFiber j q).card ^ 2 +
        (Nat.choose 128 j - (ConcreteMaps.kernelLayer j).card) ^ 2 ≤
      2 * (Nat.choose 128 j - (ConcreteMaps.kernelLayer j).card) * (syndromeFiber j q).card +
        (2 ^ 19 - 2) * ((equalSyndromePairs j).card - (ConcreteMaps.kernelLayer j).card ^ 2) :=
  concrete_variance_bound j hq

example (j : ℕ) (q : Finset (Fin 19)) :
    (syndromeFiber j q).card * j.choose (j - 1) ≤ Nat.choose 128 (j - 1) :=
  fiber_packing j q

example (j : ℕ) (q : Finset (Fin 19)) :
    (syndromeFiber j q).card * (128 - j).choose ((128 - j) - 1) ≤
      Nat.choose 128 ((128 - j) - 1) := fiber_packing_compl j q

end Spin.Structured.ConcreteFiberPin
