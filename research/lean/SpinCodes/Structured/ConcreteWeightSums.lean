import SpinCodes.Structured.ConcreteFourier

/-! Group character sums by the actual output weight. The bridge from packed
histograms to these support counts uses an equivalence, so no enumeration of
the noncomputable support representation is required. -/

noncomputable section
namespace Spin.Structured
open Finset

def weightCounts {s t : ℕ} (F : Finset (Fin s) → Finset (Fin t)) (i : ℕ) : ℕ :=
  (univ.filter fun q => (F q).card = i).card

theorem sum_by_weight {s t : ℕ} (F : Finset (Fin s) → Finset (Fin t)) (f : ℕ → ℤ) :
    ∑ q : Finset (Fin s), f (F q).card =
      ∑ w : Fin (t + 1), (weightCounts F w : ℤ) * f w := by
  let g : Finset (Fin s) → Fin (t + 1) := fun q =>
    ⟨(F q).card, Nat.lt_succ_of_le (by simpa using Finset.card_le_univ (F q))⟩
  have h := Finset.sum_fiberwise' (univ : Finset (Finset (Fin s))) g
    (fun w => f w.val)
  change (∑ w : Fin (t + 1), ∑ q ∈ univ.filter (fun q => g q = w), f w.val) =
    ∑ q : Finset (Fin s), f (F q).card at h
  rw [← h]
  apply Finset.sum_congr rfl
  intro w _
  have he : univ.filter (fun q => g q = w) = univ.filter (fun q => (F q).card = w.val) := by
    ext q
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Fin.ext_iff, g]
  rw [he, Finset.sum_const]
  simp only [nsmul_eq_mul, weightCounts]

theorem weightCounts_packed {s t : ℕ} (F : Fin (2 ^ s) → Fin (2 ^ t)) (i : ℕ) :
    weightCounts (fun q => PackedMap.supportEquiv t (F ((PackedMap.supportEquiv s).symm q))) i =
      (univ.filter fun q => PackedMap.weight t (F q) = i).card := by
  apply Finset.card_equiv (PackedMap.supportEquiv s).symm
  intro q
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  rw [PackedMap.weight_eq_card_support]
  rfl

namespace ConcreteMaps

theorem Aset_weightCounts_packed (i : ℕ) :
    weightCounts Aset i = (univ.filter fun q => PackedMap.weight 128 (A q) = i).card :=
  weightCounts_packed A i

theorem CtransposeSet_weightCounts_packed (i : ℕ) :
    weightCounts CtransposeSet i =
      (univ.filter fun q => PackedMap.weight 128 (Ctranspose q) = i).card :=
  weightCounts_packed Ctranspose i

theorem kernel_fourier_weights (j : ℕ) :
    (2 ^ 19 : ℤ) * ((kernelLayer j).card : ℤ) =
      ∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) * Spin.krawtchouk 128 j w := by
  rw [kernel_fourier, sum_by_weight]

theorem feedback_parseval_weights (j : ℕ) :
    (2 ^ 19 : ℤ) * ((equalSyndromePairs j).card : ℤ) =
      ∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) * (Spin.krawtchouk 128 j w) ^ 2 := by
  rw [feedback_parseval, sum_by_weight CtransposeSet (fun w => (Spin.krawtchouk 128 j w) ^ 2)]

theorem fiber_le_fourier_abs_weights (j : ℕ) (q : Finset (Fin 19)) :
    (2 ^ 19 : ℤ) * ((syndromeFiber j q).card : ℤ) ≤
      ∑ w : Fin 129, (weightCounts CtransposeSet w : ℤ) * |Spin.krawtchouk 128 j w| := by
  have h := fiber_le_fourier_abs j q
  rwa [sum_by_weight CtransposeSet (fun w => |Spin.krawtchouk 128 j w|)] at h

end ConcreteMaps
end Spin.Structured
