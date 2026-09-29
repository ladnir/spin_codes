import SpinCodes.Structured.ConcreteKernel
import SpinCodes.Structured.ConcreteVariance
import SpinCodes.Structured.ConcreteWeightSums
import SpinCodes.Structured.FiberCertificate

/-! Sufficient integer checks for each of the paper's four fiber bounds.
All counting terms refer to the actual feedback map. Spectrum substitution
and finite arithmetic certificates can therefore be checked separately. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset

theorem fiber_cap_complement (j cap : ℕ) {q : Finset (Fin 19)} (hq : q ≠ ∅)
    (hc : Nat.choose 128 j - (kernelLayer j).card ≤ cap) :
    (syndromeFiber j q).card ≤ cap := by
  have h := nonzero_fiber_le_complement j hq
  omega

theorem fiber_cap_fourier (j cap : ℕ) (q : Finset (Fin 19))
    (hc : (∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) *
        |Spin.krawtchouk 128 j w|) < (2 ^ 19 : ℤ) * (cap + 1)) :
    (syndromeFiber j q).card ≤ cap := by
  have h := fiber_le_fourier_abs_weights j q
  omega

theorem fiber_cap_packing (j cap : ℕ) (q : Finset (Fin 19))
    (hc : Nat.choose 128 (j - 1) < (cap + 1) * j.choose (j - 1)) :
    (syndromeFiber j q).card ≤ cap := by
  have h := fiber_packing j q
  by_contra ht
  have ht' : cap + 1 ≤ (syndromeFiber j q).card := by omega
  have hm := Nat.mul_le_mul_right (j.choose (j - 1)) ht'
  omega

theorem fiber_cap_packing_compl (j cap : ℕ) (q : Finset (Fin 19))
    (hc : Nat.choose 128 ((128 - j) - 1) <
      (cap + 1) * (128 - j).choose ((128 - j) - 1)) :
    (syndromeFiber j q).card ≤ cap := by
  have h := fiber_packing_compl j q
  by_contra ht
  have ht' : cap + 1 ≤ (syndromeFiber j q).card := by omega
  have hm := Nat.mul_le_mul_right ((128 - j).choose ((128 - j) - 1)) ht'
  omega

theorem fiber_cap_variance (j cap : ℕ) {q : Finset (Fin 19)} (hq : q ≠ ∅)
    (hm : Nat.choose 128 j - (kernelLayer j).card ≤ (2 ^ 19 - 1) * (cap + 1))
    (hc : 2 * (Nat.choose 128 j - (kernelLayer j).card) * (cap + 1) +
        (2 ^ 19 - 2) * ((equalSyndromePairs j).card - (kernelLayer j).card ^ 2) <
      (2 ^ 19 - 1) * (cap + 1) ^ 2 + (Nat.choose 128 j - (kernelLayer j).card) ^ 2) :
    (syndromeFiber j q).card ≤ cap :=
  Spin.variance_cap_of_certificate (concrete_variance_bound j hq) hm hc

end Spin.Structured.ConcreteMaps
