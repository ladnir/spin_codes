import SpinCodes.Structured.LowCancellationRows
import SpinCodes.Structured.LowCancellationData.Weight1
import SpinCodes.Structured.LowCancellationData.Weight2

/-! The exceptional weight-one and weight-two cancellation envelopes for the
actual hexadecimal maps. All 8,256 inputs are covered by checked certificates. -/
noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt LowCancellation LowCancellation.Data

theorem cancellationMoment_le_patterns_one (z : ℝ) (q) :
    cancellationMoment 1 z q ≤
      ((SparsePolynomial.Data.weight1.lowPatterns.getD []).map
        fun p => Occupation.Sparse.pattern 1 p z).foldr max 0 :=
  cancellation_le_patterns groups1_valid groups1_inputs groups1_length groups1_keys groups1_patterns z q

theorem cancellationMoment_le_patterns_two (z : ℝ) (q) :
    cancellationMoment 2 z q ≤
      ((SparsePolynomial.Data.weight2.lowPatterns.getD []).map
        fun p => Occupation.Sparse.pattern 2 p z).foldr max 0 :=
  cancellation_le_patterns groups2_valid groups2_inputs groups2_length groups2_keys groups2_patterns z q

theorem shellCancellationMoment_one (z : ℝ) (i : Fin 5) :
    shellCancellationMoment 1 z i = Occupation.Sparse.pattern 1
      (SparsePolynomial.Data.weight1.lowShells.getD i []) z / shellCount i := by
  apply shell_cancellation_pattern groups1_valid groups1_inputs groups1_length i _ z
  fin_cases i
  · exact groups1_shell0
  · exact groups1_shell1
  · exact groups1_shell2
  · exact groups1_shell3
  · exact groups1_shell4

theorem shellCancellationMoment_two (z : ℝ) (i : Fin 5) :
    shellCancellationMoment 2 z i = Occupation.Sparse.pattern 2
      (SparsePolynomial.Data.weight2.lowShells.getD i []) z / shellCount i := by
  apply shell_cancellation_pattern groups2_valid groups2_inputs groups2_length i _ z
  fin_cases i
  · exact groups2_shell0
  · exact groups2_shell1
  · exact groups2_shell2
  · exact groups2_shell3
  · exact groups2_shell4

theorem cancellationMoment_le_rowD_all {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelD := by
  by_cases hj1 : j.val = 1
  · have hj : j = ⟨1, by decide⟩ := Fin.ext hj1
    subst j
    exact cancellationMoment_le_rowD_of_patterns hz hz1 _ _ rfl hq
      (cancellationMoment_le_patterns_one z q)
  by_cases hj2 : j.val = 2
  · have hj : j = ⟨2, by decide⟩ := Fin.ext hj2
    subst j
    exact cancellationMoment_le_rowD_of_patterns hz hz1 _ _ rfl hq
      (cancellationMoment_le_patterns_two z q)
  exact cancellationMoment_le_rowD hz hz1 j hj1 hj2 hq

theorem shellCancellationMoment_le_rowS_all {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (i : Fin 5) :
    shellCancellationMoment j z i ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i := by
  by_cases hj1 : j.val = 1
  · have hj : j = ⟨1, by decide⟩ := Fin.ext hj1
    subst j
    exact shellCancellationMoment_le_rowS_of_patterns hz hz1 _ _ rfl i
      (shellCancellationMoment_one z i).le
  by_cases hj2 : j.val = 2
  · have hj : j = ⟨2, by decide⟩ := Fin.ext hj2
    subst j
    exact shellCancellationMoment_le_rowS_of_patterns hz hz1 _ _ rfl i
      (shellCancellationMoment_two z i).le
  exact shellCancellationMoment_le_rowS hz hz1 j hj1 hj2 i

end Spin.Structured.ConcreteMaps
