import SpinCodes.Structured.ConcreteRows

/-! Identification of the generic cancellation envelopes with the numerical
matrix for every weight except the two special low-input tables. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem lowPatterns_none_checked : ∀ j : Fin 129, j.val ≠ 1 → j.val ≠ 2 →
    (SparsePolynomial.weightN j).lowPatterns = none := by decide +kernel

theorem shellCount_getD (i : Fin 5) : shellCount i = SparsePolynomial.counts.getD i 1 := by
  fin_cases i <;> rfl

theorem cancellationMoment_le_rowD {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2)
    {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelD := by
  simp only [Occupation.Sparse.row, Occupation.Sparse.cancellationBoundsD,
    lowPatterns_none_checked j hj1 hj2, List.append_nil, List.foldr_cons, List.foldr_nil]
  exact le_min (cancellationMoment_le_maximum hz j hq)
    (le_min (cancellationMoment_le_live hz hz1 j hq)
      (le_min (cancellationMoment_le_minDistance hz hz1 j hq)
        (cancellationMoment_le_maximum hz j hq)))

theorem shellCancellationMoment_le_rowS {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2) (i : Fin 5) :
    shellCancellationMoment j z i ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i := by
  have hm := shellCancellationMoment_le_moment hz j i
  have hc := shellCancellationMoment_le_cap hz hz1 j i
  simp only [shellWeight, shellCount_getD] at hm hc
  simp only [Occupation.Sparse.row, Occupation.Sparse.cancellationBoundsS,
    lowPatterns_none_checked j hj1 hj2, List.append_nil, List.foldr_cons, List.foldr_nil]
  exact le_min hm (le_min hc hm)

theorem stepRow_dominates_D_generic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2)
    {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    LiveDominates actualShellSystem
      (Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowD (stepRow j z q) :=
  stepRow_dominates_D hz hz1 j hq (cancellationMoment_le_rowD hz hz1 j hj1 hj2 hq)

end Spin.Structured.ConcreteMaps
