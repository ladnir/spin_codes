import SpinCodes.Turbo.SparseListSound22

/-! Soundness of the scaled squaring chain: integer mantissa matrices with an
explicit binary exponent per stage, so that powers whose entries grow or
decay far beyond any absolute grid stay bounded in size. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix SparseFixedSound22 SparseListSound22

namespace SparseListSquare22

def chainExp (e : ℤ) : List (Nat × Nat × List (List Int)) → ℤ
  | [] => e
  | st :: rest => chainExp (2 * e + st.2.1 - st.1) rest

theorem matScale_spec (c : Int) (A : List (List Int)) (hA : MatShape A) :
    MatShape (sparseMatScaleL22 c A) ∧ ∀ i j, i < 11 → j < 11 →
      sparseLMat22 (sparseMatScaleL22 c A) i j = c * sparseLMat22 A i j := by
  refine ⟨⟨by simp [sparseMatScaleL22, hA.1], fun r hr => ?_⟩, fun i j hi hj => ?_⟩
  · obtain ⟨row, hrow, rfl⟩ := List.mem_map.mp hr
    simp [hA.2 row hrow]
  · have hiA : i < A.length := by rw [hA.1]; exact hi
    have hrow : (A[i]'hiA).length = 11 := hA.2 _ (List.getElem_mem hiA)
    show ((A.map fun row => row.map (c * ·)).getD i []).getD j 0 = c * (A.getD i []).getD j 0
    rw [gd (A.map fun row => row.map (c * ·)) [] (by simpa using hiA), List.getElem_map,
      gd A [] hiA, map_getD (by omega)]

theorem stage_sound (s : ℤ) (hs : 0 < s) (A : Matrix (Fin 11) (Fin 11) ℝ)
    (hA0 : ∀ i j, 0 ≤ A i j) (U : List (List Int)) (hU : MatShape U) (E : ℤ)
    (hAU : ∀ i j : Fin 11, A i j ≤ (2 : ℝ) ^ E * (sparseLMat22 U i j : ℝ) / s)
    (st : Nat × Nat × List (List Int)) (hc : sparseSquareStageCheckL22 s U st = true) :
    MatShape st.2.2 ∧ ∀ i j : Fin 11,
      (A * A) i j ≤ (2 : ℝ) ^ (2 * E + st.2.1 - st.1) * (sparseLMat22 st.2.2 i j : ℝ) / s := by
  obtain ⟨down, up, U'⟩ := st
  simp only [sparseSquareStageCheckL22, Bool.and_eq_true] at hc
  obtain ⟨hshape, hle⟩ := hc
  have hU' := matShape_of_check U' hshape
  refine ⟨hU', fun i j => ?_⟩
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  have hmm := matMul_spec U U hU hU
  have hsc := matScale_spec (2 ^ down) _ hmm.1
  have hc := matLeScaled_sound _ _ _ hsc.1 hU' hle i j i.isLt j.isLt
  rw [hsc.2 i j i.isLt j.isLt, hmm.2 i j i.isLt j.isLt] at hc
  have hcR : (2 : ℝ) ^ down * (∑ k : Fin 11, (sparseLMat22 U i k : ℝ) * sparseLMat22 U k j) ≤
      (s : ℝ) * 2 ^ up * sparseLMat22 U' i j := by exact_mod_cast hc
  have h2E : (0 : ℝ) < (2 : ℝ) ^ E := zpow_pos (by norm_num) E
  have hdot : (A * A) i j ≤ (∑ k : Fin 11, ((2 : ℝ) ^ E * sparseLMat22 U i k) *
      ((2 : ℝ) ^ E * sparseLMat22 U k j)) / (s * s) := by
    rw [Matrix.mul_apply]
    exact dot_le _ hsR _ _ _ _ (fun k => hA0 i k) (fun k => hA0 k j)
      (fun k => by simpa [mul_div_assoc] using hAU i k)
      (fun k => by simpa [mul_div_assoc] using hAU k j)
  refine hdot.trans ?_
  have hsum : (∑ k : Fin 11, ((2 : ℝ) ^ E * sparseLMat22 U i k) *
      ((2 : ℝ) ^ E * sparseLMat22 U k j)) =
      ((2 : ℝ) ^ E * (2 : ℝ) ^ E) *
        ∑ k : Fin 11, (sparseLMat22 U i k : ℝ) * sparseLMat22 U k j := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [hsum]
  have hexp : (2 : ℝ) ^ (2 * E + (up : ℤ) - (down : ℤ)) =
      (2 : ℝ) ^ E * (2 : ℝ) ^ E * (2 : ℝ) ^ up / (2 : ℝ) ^ down := by
    rw [zpow_sub₀ (by norm_num), zpow_add₀ (by norm_num), two_mul,
      zpow_add₀ (by norm_num), zpow_natCast, zpow_natCast]
  rw [hexp]
  have h2d : (0 : ℝ) < (2 : ℝ) ^ down := by positivity
  rw [div_le_div_iff₀ (mul_pos hsR hsR) hsR]
  have key : (2 : ℝ) ^ E * (2 : ℝ) ^ E *
      (∑ k : Fin 11, (sparseLMat22 U i k : ℝ) * sparseLMat22 U k j) ≤
      (2 : ℝ) ^ E * (2 : ℝ) ^ E * ((s : ℝ) * 2 ^ up * sparseLMat22 U' i j / 2 ^ down) := by
    apply mul_le_mul_of_nonneg_left _ (by positivity)
    rw [le_div_iff₀ h2d]
    linarith
  calc (2 : ℝ) ^ E * (2 : ℝ) ^ E *
        (∑ k : Fin 11, (sparseLMat22 U i k : ℝ) * sparseLMat22 U k j) * s
      ≤ (2 : ℝ) ^ E * (2 : ℝ) ^ E * ((s : ℝ) * 2 ^ up * sparseLMat22 U' i j / 2 ^ down) * s :=
        mul_le_mul_of_nonneg_right key hsR.le
    _ = _ := by field_simp

theorem chain_sound (s : ℤ) (hs : 0 < s) :
    ∀ (chain : List (Nat × Nat × List (List Int))) (U : List (List Int)) (E : ℤ)
      (A : Matrix (Fin 11) (Fin 11) ℝ), (∀ i j, 0 ≤ A i j) → MatShape U →
      (∀ i j : Fin 11, A i j ≤ (2 : ℝ) ^ E * (sparseLMat22 U i j : ℝ) / s) →
      sparseSquareChainCheckL22 s U chain = true →
      ∀ i j : Fin 11, (A ^ (2 ^ chain.length)) i j ≤
        (2 : ℝ) ^ (chainExp E chain) * (sparseLMat22 (sparseSquareChainLastL22 U chain) i j : ℝ) / s
  | [], U, E, A, _, _, hAU, _ => by simpa [chainExp, sparseSquareChainLastL22] using hAU
  | st :: rest, U, E, A, hA0, hU, hAU, hc => by
    simp only [sparseSquareChainCheckL22, Bool.and_eq_true] at hc
    obtain ⟨hst, hrest⟩ := hc
    have h1 := stage_sound s hs A hA0 U hU E hAU st hst
    have hAA0 : ∀ i j, 0 ≤ (A * A) i j := fun i j => by
      rw [Matrix.mul_apply]; exact Finset.sum_nonneg fun k _ => mul_nonneg (hA0 i k) (hA0 k j)
    have h2 := chain_sound s hs rest st.2.2 (2 * E + st.2.1 - st.1) (A * A) hAA0 h1.1 h1.2 hrest
    intro i j
    have hpow : A ^ (2 ^ (st :: rest).length) = (A * A) ^ (2 ^ rest.length) := by
      rw [List.length_cons, pow_succ', pow_mul, sq]
    rw [hpow]
    simpa [chainExp, sparseSquareChainLastL22] using h2 i j

end SparseListSquare22

open SparseListSquare22 in
#print axioms chain_sound
end Spin.Turbo
