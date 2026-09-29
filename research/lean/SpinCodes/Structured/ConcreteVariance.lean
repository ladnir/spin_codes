import SpinCodes.Structured.ConcreteFourier
import SpinCodes.Structured.VarianceBound

/-! The variance inequality for the actual feedback map, indexed by all
nonzero syndromes. Empty fibers are included, so the number of terms is
exactly `2^19 - 1` without an assumption on which syndromes a layer meets. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset

def nonzeroSyndromes : Finset (Finset (Fin 19)) := univ.erase ∅

theorem nonzeroSyndromes_card : nonzeroSyndromes.card = 2 ^ 19 - 1 := by
  simp [nonzeroSyndromes]

theorem syndromeFiber_zero (j : ℕ) : syndromeFiber j ∅ = kernelLayer j := rfl

theorem sum_syndromeFiber_card (j : ℕ) :
    ∑ q : Finset (Fin 19), (syndromeFiber j q).card = Nat.choose 128 j := by
  rw [← Spin.card_layer]
  exact (Finset.card_eq_sum_card_fiberwise (f := Cset)
    (s := Spin.layer 128 j) (t := univ) (fun _ _ => Finset.mem_univ _)).symm

theorem sum_syndromeFiber_sq (j : ℕ) :
    ∑ q : Finset (Fin 19), (syndromeFiber j q).card ^ 2 = (equalSyndromePairs j).card := by
  have hrow : (equalSyndromePairs j).card =
      ∑ x ∈ Spin.layer 128 j, (syndromeFiber j (Cset x)).card := by
    rw [equalSyndromePairs, Finset.card_filter, Finset.sum_product]
    apply Finset.sum_congr rfl
    intro x _
    simp only [syndromeFiber, Finset.card_filter, eq_comm]
  rw [hrow]
  have h := Finset.sum_fiberwise_of_maps_to
    (s := Spin.layer 128 j) (t := (univ : Finset (Finset (Fin 19))))
    (g := Cset) (fun _ _ => Finset.mem_univ _)
    (fun x => (syndromeFiber j (Cset x)).card)
  rw [← h]
  apply Finset.sum_congr rfl
  intro q _
  have he : ∑ x ∈ (Spin.layer 128 j).filter (fun x => Cset x = q),
      (syndromeFiber j (Cset x)).card = (syndromeFiber j q).card ^ 2 := by
    rw [Finset.sum_congr rfl (fun x hx => by rw [(Finset.mem_filter.mp hx).2])]
    simp only [Finset.sum_const, smul_eq_mul, syndromeFiber, pow_two]
  exact he.symm

theorem sum_nonzero_fiber_card (j : ℕ) :
    ∑ q ∈ nonzeroSyndromes, (syndromeFiber j q).card =
      Nat.choose 128 j - (kernelLayer j).card := by
  have h := Finset.sum_erase_add (univ : Finset (Finset (Fin 19)))
    (fun q => (syndromeFiber j q).card) (Finset.mem_univ ∅)
  rw [sum_syndromeFiber_card, syndromeFiber_zero] at h
  exact Nat.eq_sub_of_add_eq h

theorem sum_nonzero_fiber_sq (j : ℕ) :
    ∑ q ∈ nonzeroSyndromes, (syndromeFiber j q).card ^ 2 =
      (equalSyndromePairs j).card - (kernelLayer j).card ^ 2 := by
  have h := Finset.sum_erase_add (univ : Finset (Finset (Fin 19)))
    (fun q => (syndromeFiber j q).card ^ 2) (Finset.mem_univ ∅)
  rw [sum_syndromeFiber_sq, syndromeFiber_zero] at h
  exact Nat.eq_sub_of_add_eq h

theorem concrete_variance_bound (j : ℕ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (2 ^ 19 - 1) * (syndromeFiber j q).card ^ 2 +
        (Nat.choose 128 j - (kernelLayer j).card) ^ 2 ≤
      2 * (Nat.choose 128 j - (kernelLayer j).card) * (syndromeFiber j q).card +
        (2 ^ 19 - 2) * ((equalSyndromePairs j).card - (kernelLayer j).card ^ 2) := by
  exact Spin.sq_max_le nonzeroSyndromes (fun p => (syndromeFiber j p).card)
    (by simp [nonzeroSyndromes, hq])
    (by rw [nonzeroSyndromes_card]; norm_num)
    (sum_nonzero_fiber_card j) (sum_nonzero_fiber_sq j)

end Spin.Structured.ConcreteMaps
