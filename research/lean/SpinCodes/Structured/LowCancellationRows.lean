import SpinCodes.Structured.LowCancellationShells
import SpinCodes.Structured.ConcreteCancellationRows

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

theorem cancellationMoment_le_rowD_of_patterns {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (ps : List (List (Nat × Nat)))
    (hp : (SparsePolynomial.weightN j).lowPatterns = some ps)
    {q : Finset (Fin 19)} (hq : q ≠ ∅)
    (hc : cancellationMoment j z q ≤
      (ps.map fun p => Occupation.Sparse.pattern j p z).foldr max 0) :
    cancellationMoment j z q ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelD := by
  simp only [Occupation.Sparse.row, Occupation.Sparse.cancellationBoundsD, hp,
    List.cons_append, List.nil_append, List.foldr_cons, List.foldr_nil]
  exact le_min (cancellationMoment_le_maximum hz j hq)
    (le_min (cancellationMoment_le_live hz hz1 j hq)
      (le_min (cancellationMoment_le_minDistance hz hz1 j hq)
        (le_min hc (cancellationMoment_le_maximum hz j hq))))

theorem shellCancellationMoment_le_rowS_of_patterns {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (ps : List (List (Nat × Nat)))
    (hp : (SparsePolynomial.weightN j).lowPatterns = some ps) (i : Fin 5)
    (hc : shellCancellationMoment j z i ≤
      Occupation.Sparse.pattern j ((SparsePolynomial.weightN j).lowShells.getD i []) z / shellCount i) :
    shellCancellationMoment j z i ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i := by
  have hm := shellCancellationMoment_le_moment hz j i
  have hf := shellCancellationMoment_le_cap hz hz1 j i
  simp only [shellWeight, shellCount_getD] at hm hf hc
  simp only [Occupation.Sparse.row, Occupation.Sparse.cancellationBoundsS, hp,
    List.cons_append, List.nil_append, List.foldr_cons, List.foldr_nil]
  exact le_min hm (le_min hf (le_min hc hm))

end Spin.Structured.ConcreteMaps
