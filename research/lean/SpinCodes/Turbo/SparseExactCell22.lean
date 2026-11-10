import SpinCodes.Turbo.SparseExactRegion22

/-! Exact sparse support laws are products of complete nonzero packet labels.
This reuses `RS8ReferenceRegion22`'s transpose proof after replacing each
selected-cell uniform law by its exact nonzero version. -/
noncomputable section
namespace Spin.Turbo
open Finset

def sparseNonzeroPacketLaw22 : FinPMF Packet4 :=
  (FinPMF.uniform {p : Packet4 // p ≠ 0}).map Subtype.val

theorem sparseNonzeroPacketLaw22_apply (x : Packet4) :
    sparseNonzeroPacketLaw22.p x = if x ≠ 0 then 1 / (15 : ℝ) else 0 := by
  rw [sparseNonzeroPacketLaw22, FinPMF.map_p]
  by_cases hx : x ≠ 0
  · rw [if_pos hx]
    let xx : {p : Packet4 // p ≠ 0} := ⟨x, hx⟩
    have he (p : {p : Packet4 // p ≠ 0}) : p.val = x ↔ p = xx := by
      constructor
      · intro h
        apply Subtype.ext
        exact h
      · intro hp
        exact congrArg Subtype.val hp
    simp only [he, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    change 1 / (Fintype.card {p : Packet4 // p ≠ 0} : ℝ) = _
    rw [packet4_nonzero_card]
    norm_num
  · rw [if_neg hx]
    have hz : x = 0 := not_not.mp hx
    apply Finset.sum_eq_zero
    intro p _
    have hp : p.val ≠ x := by
      intro h
      exact p.property (h.trans hz)
    simp only [if_neg hp]

def sparseExactCellLaw22 {n : ℕ} (T : Finset (Fin n)) (i : Fin n) : FinPMF Packet4 :=
  if i ∈ T then sparseNonzeroPacketLaw22 else zeroPacketLaw22

theorem regionalNonzeroGiven22_eq_cell_product {n : ℕ} (T : Finset (Fin n)) :
    regionalNonzeroGiven22 T = Spin.piPMF (sparseExactCellLaw22 T) := by
  apply FinPMF.ext
  funext x
  rw [regionalNonzeroGiven22_apply, Spin.piPMF_apply]
  by_cases hT : regionalPacketSupport22 x = T
  · rw [if_pos hT]
    have hterm (i : Fin n) : (sparseExactCellLaw22 T i).p (x i) =
        if i ∈ T then 1 / (15 : ℝ) else 1 := by
      have hi : x i ≠ 0 ↔ i ∈ T := by
        have hh := Finset.ext_iff.mp hT i
        simpa only [regionalPacketSupport22, Finset.mem_filter, Finset.mem_univ, true_and] using hh
      by_cases ht : i ∈ T
      · simp only [sparseExactCellLaw22, ht, ite_true]
        rw [sparseNonzeroPacketLaw22_apply, if_pos (hi.mpr ht)]
      · have hx : x i = 0 := by simpa using (not_congr hi).mpr ht
        simp only [sparseExactCellLaw22, ht, ite_false, zeroPacketLaw22_apply, hx, ite_true]
    simp_rw [hterm]
    rw [Finset.prod_ite]
    simp [div_eq_mul_inv]
  · rw [if_neg hT]
    have hw : ∃ i : Fin n, ¬ (x i ≠ 0 ↔ i ∈ T) := by
      by_contra hn
      apply hT
      ext i
      have hi := not_not.mp (not_exists.mp hn i)
      simpa only [regionalPacketSupport22, Finset.mem_filter, Finset.mem_univ, true_and] using hi
    obtain ⟨i, hi⟩ := hw
    symm
    apply Finset.prod_eq_zero (Finset.mem_univ i)
    by_cases ht : i ∈ T
    · have hz : x i = 0 := by
        by_contra hne
        exact hi (iff_of_true hne ht)
      simp only [sparseExactCellLaw22, ht, ite_true, sparseNonzeroPacketLaw22_apply,
        hz, ne_eq, not_true_eq_false, ite_false]
    · have hz : x i ≠ 0 := by
        intro he
        exact hi (iff_of_false (by simpa [he]) ht)
      simp only [sparseExactCellLaw22, ht, ite_false, zeroPacketLaw22_apply,
        hz, ite_false]

#print axioms regionalNonzeroGiven22_eq_cell_product
end Spin.Turbo

