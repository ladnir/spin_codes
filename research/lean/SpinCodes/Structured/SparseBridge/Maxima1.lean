import SpinCodes.Structured.SparseSelection
import SpinCodes.Structured.SparseBridge.Dominance1
noncomputable section
namespace Spin.Structured.SparsePolynomial
set_option maxRecDepth 100000
set_option maxHeartbeats 0
theorem maximum_1 (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    Spin.Imt.Occupation.Sparse.maximumMoment 1 (1 - x / 6250) ≤
      (arbitrary 1 Data.weight1).eval x := by
  have hs : (arbitrary 1 Data.weight1).eval x = Dominance1.m0.eval x := by
    change (moment 48 1).eval x = _
    exact RatPoly.checkEq_sound _ _ Dominance1.m0_checked x
  apply maximumMoment_le_polynomial 1 x h1 (arbitrary 1 Data.weight1)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.m0_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.m0_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.m1_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.m1_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.m2_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.m2_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.m3_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.m3_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.m4_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.m4_le_checked h0 h1).trans_eq hs.symm)
theorem lowMaximum_1 (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    lowMaximum 1 Data.weight1 (1 - x / 6250) ≤
      (selectedLowPolynomial 1 Data.weight1).eval x := by
  have he : selectedLowPolynomial 1 Data.weight1 = pattern 1 [(49, 1)] := rfl
  have hs : (selectedLowPolynomial 1 Data.weight1).eval x = Dominance1.low0.eval x := by
    rw [he]
    exact RatPoly.checkEq_sound _ _ Dominance1.low0_checked x
  have hn : 0 ≤ (selectedLowPolynomial 1 Data.weight1).eval x := by
    rw [he, eval_pattern]
    exact Spin.Imt.Occupation.Sparse.pattern_nonneg _ _ (by linarith)
  change ([[(49, 1)], [(55, 1)], [(57, 1)], [(63, 1)], [(65, 1)], [(71, 1)], [(73, 1)], [(79, 1)], [(81, 1)]].map fun p => Spin.Imt.Occupation.Sparse.pattern 1 p (1 - x / 6250)).foldr max 0 ≤ _
  apply maximumPattern_le_polynomial 1 _ x _ hn
  intro p hp
  simp only [List.mem_cons, List.not_mem_nil, or_false] at hp
  rcases hp with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low0_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low0_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low1_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low1_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low2_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low2_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low3_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low3_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low4_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low4_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low5_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low5_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low6_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low6_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low7_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low7_le_checked h0 h1).trans_eq hs.symm)
  · exact (RatPoly.checkEq_sound _ _ Dominance1.low8_checked x).le.trans
      ((RatPoly.checkLE_sound _ _ Dominance1.low8_le_checked h0 h1).trans_eq hs.symm)
end Spin.Structured.SparsePolynomial
