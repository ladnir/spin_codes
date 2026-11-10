import SpinCodes.Turbo.SparseFixedDefs22
import SpinCodes.Turbo.MatrixMarker22
import Mathlib.Algebra.BigOperators.Field

/-! Soundness of sparse integer checks against real matrices. The exact
integer checker remains in the small companion Definitions module. The
Collatz argument is the matrix-index version of audited
`Structured.Collatz.collatz_total`, preserving the initial row and all
terminal coordinates. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix

private theorem sparseAll11Bool22_sound (f : Nat → Bool)
    (h : sparseAll11Bool22 f = true) (i : Fin 11) : f i.val = true := by
  simp only [sparseAll11Bool22, Bool.and_eq_true] at h
  fin_cases i <;> tauto

private theorem sparseSum11Int22_eq_finSum (f : Nat → Int) :
    sparseSum11Int22 f = ∑ i : Fin 11, f i.val := by
  simp [sparseSum11Int22, Fin.sum_univ_succ]
  ring
private theorem sparse_cross_upper22 {a b c d : Int} (hb : 0 < b) (hd : 0 < d)
    (h : a * d ≤ c * b) : (a : ℝ) / b ≤ (c : ℝ) / d := by
  apply (div_le_div_iff₀ (by exact_mod_cast hb : (0 : ℝ) < b)
    (by exact_mod_cast hd : (0 : ℝ) < d)).mpr
  exact_mod_cast h

theorem sparseIntMatrixMulUpperCheck22_sound
    (scale : Int) (hs : 0 < scale) (a b u : Array Int)
    (hcheck : sparseIntMatrixMulUpperCheck22 scale a b u = true)
    (M N : Matrix (Fin 11) (Fin 11) ℝ)
    (hM : ∀ i j, 0 ≤ M i j) (hN : ∀ i j, 0 ≤ N i j)
    (hMa : ∀ i j, M i j ≤ (sparseIntMatrix22 a i j : ℝ) / scale)
    (hNb : ∀ i j, N i j ≤ (sparseIntMatrix22 b i j : ℝ) / scale) :
    ∀ i j, (M * N) i j ≤ (sparseIntMatrix22 u i j : ℝ) / scale := by
  have hc : ∀ i j : Fin 11,
      (∑ k, sparseIntMatrix22 a i k * sparseIntMatrix22 b k j) ≤
        scale * sparseIntMatrix22 u i j := by
    intro i j
    have hi := sparseAll11Bool22_sound _ hcheck i
    have hij := sparseAll11Bool22_sound _ hi j
    have hx := of_decide_eq_true hij
    simpa only [sparseSum11Int22_eq_finSum, sparseIntMatrix22] using hx
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  intro i j
  calc
    _ ≤ ∑ k : Fin 11,
        ((sparseIntMatrix22 a i k : ℝ) / scale) * ((sparseIntMatrix22 b k j : ℝ) / scale) := by
      rw [Matrix.mul_apply]
      exact Finset.sum_le_sum fun k _ => mul_le_mul (hMa i k) (hNb k j) (hN k j)
        ((hM i k).trans (hMa i k))
    _ = ((∑ k : Fin 11, sparseIntMatrix22 a i k * sparseIntMatrix22 b k j : Int) : ℝ) /
        ((scale : ℝ) * scale) := by
      simp only [Int.cast_sum, Int.cast_mul, Finset.sum_div]
      apply Finset.sum_congr rfl
      intro k _
      ring
    _ ≤ ((scale * sparseIntMatrix22 u i j : Int) : ℝ) / ((scale : ℝ) * scale) := by
      apply div_le_div_of_nonneg_right _ (mul_nonneg hsR.le hsR.le)
      exact_mod_cast hc i j
    _ = (sparseIntMatrix22 u i j : ℝ) / scale := by
      push_cast
      field_simp

/-- One positive vector bounds the whole chronological matrix power. -/
theorem sparseMatrixCollatz_total22 {I : Type*} [Fintype I] [DecidableEq I]
    (M : Matrix I I ℝ) (hM : ∀ i j, 0 ≤ M i j)
    (v : I → ℝ) (rho : ℝ) (hrho : 0 ≤ rho)
    (hrow : ∀ i, (M *ᵥ v) i ≤ rho * v i)
    (vmin : ℝ) (hvmin : 0 < vmin) (hv : ∀ i, vmin ≤ v i) (n : ℕ) (i : I) :
    (∑ j, (M ^ n) i j) ≤ rho ^ n * v i / vmin := by
  have hpow : ∀ n : ℕ, ∀ i, ((M ^ n) *ᵥ v) i ≤ rho ^ n * v i := by
    intro n
    induction n with
    | zero => intro i; simp
    | succ n ih =>
      intro i
      rw [pow_succ', ← Matrix.mulVec_mulVec]
      change (∑ k, M i k * ((M ^ n) *ᵥ v) k) ≤ _
      calc
        _ ≤ ∑ k, M i k * (rho ^ n * v k) :=
          Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left (ih k) (hM i k)
        _ = rho ^ n * ((M *ᵥ v) i) := by
          simp only [Matrix.mulVec, dotProduct, Finset.mul_sum]
          apply Finset.sum_congr rfl
          intro k _
          ring
        _ ≤ rho ^ n * (rho * v i) :=
          mul_le_mul_of_nonneg_left (hrow i) (pow_nonneg hrho n)
        _ = _ := by rw [pow_succ]; ring
  apply (le_div_iff₀ hvmin).mpr
  calc
    _ = ∑ j, (M ^ n) i j * vmin := Finset.sum_mul _ _ _
    _ ≤ ∑ j, (M ^ n) i j * v j :=
      Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_left (hv j) (markerPowNonneg22 M hM n i j)
    _ ≤ rho ^ n * v i := hpow n i

/-- The proposed integer vector and rational radius certify the real marker
matrix after its entries have independently been bounded by the data. -/
theorem sparseIntCollatzCheck22_sound
    (scale rhoNum rhoDen : Int) (m vector : Array Int)
    (hcheck : sparseIntCollatzCheck22 scale rhoNum rhoDen m vector = true)
    (M : Matrix (Fin 11) (Fin 11) ℝ) (hM : ∀ i j, 0 ≤ M i j)
    (hMu : ∀ i j, M i j ≤ (sparseIntMatrix22 m i j : ℝ) / scale)
    (minNum : Int) (hmin : 0 < minNum)
    (hvector : ∀ i : Fin 11, minNum ≤ vector.getD i.val 0) (n : ℕ) (i : Fin 11) :
    (∑ j : Fin 11, (M ^ n) i j) ≤
      ((rhoNum : ℝ) / rhoDen) ^ n * (vector.getD i.val 0 : ℝ) / minNum := by
  have hc : 0 < scale ∧ 0 < rhoDen ∧ 0 ≤ rhoNum ∧
      (∀ i : Fin 11, 0 < vector.getD i.val 0) ∧
      ∀ i : Fin 11, rhoDen *
        (∑ j : Fin 11, sparseIntMatrix22 m i j * vector.getD j.val 0) ≤
        scale * rhoNum * vector.getD i.val 0 := by
    have hb : decide (0 < scale ∧ 0 < rhoDen ∧ 0 ≤ rhoNum) = true ∧
        ((sparseAll11Bool22 fun i => decide (0 < vector.getD i 0)) = true ∧
        (sparseAll11Bool22 fun i => decide
          (rhoDen * sparseSum11Int22 (fun j => sparseIntMatrixGet22 m i j * vector.getD j 0) ≤
            scale * rhoNum * vector.getD i 0)) = true) := by
      simpa only [sparseIntCollatzCheck22, Bool.and_eq_true] using hcheck
    have hg : 0 < scale ∧ 0 < rhoDen ∧ 0 ≤ rhoNum := of_decide_eq_true hb.1
    refine ⟨hg.1, hg.2.1, hg.2.2, ?_, ?_⟩
    · intro i
      exact of_decide_eq_true (sparseAll11Bool22_sound _ hb.2.1 i)
    · intro i
      have hi := of_decide_eq_true (sparseAll11Bool22_sound _ hb.2.2 i)
      simpa only [sparseSum11Int22_eq_finSum, sparseIntMatrix22] using hi
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hc.1
  have hdR : (0 : ℝ) < rhoDen := by exact_mod_cast hc.2.1
  have hrR : (0 : ℝ) ≤ rhoNum := by exact_mod_cast hc.2.2.1
  have hv (j : Fin 11) : (0 : ℝ) ≤ vector.getD j.val 0 := by
    exact_mod_cast (hc.2.2.2.1 j).le
  apply sparseMatrixCollatz_total22 M hM (fun j => (vector.getD j.val 0 : ℝ))
    ((rhoNum : ℝ) / rhoDen) (div_nonneg hrR hdR.le) _ (minNum : ℝ)
    (by exact_mod_cast hmin) (fun j => by exact_mod_cast hvector j) n i
  intro r
  change (∑ j : Fin 11, M r j * (vector.getD j.val 0 : ℝ)) ≤ _
  calc
    _ ≤ ∑ j : Fin 11,
        ((sparseIntMatrix22 m r j : ℝ) / scale) * (vector.getD j.val 0 : ℝ) :=
      Finset.sum_le_sum fun j _ => mul_le_mul_of_nonneg_right (hMu r j) (hv j)
    _ = ((∑ j : Fin 11, sparseIntMatrix22 m r j * vector.getD j.val 0 : Int) : ℝ) / scale := by
      simp only [Int.cast_sum, Int.cast_mul, Finset.sum_div]
      apply Finset.sum_congr rfl
      intro j _
      ring
    _ ≤ ((rhoNum : ℝ) / rhoDen) * (vector.getD r.val 0 : ℝ) := by
      have hcross : (∑ j : Fin 11, sparseIntMatrix22 m r j * vector.getD j.val 0) * rhoDen ≤
          (rhoNum * vector.getD r.val 0) * scale := by
        convert hc.2.2.2.2 r using 1 <;> ring
      simpa only [Int.cast_mul, div_mul_eq_mul_div] using
        sparse_cross_upper22 hc.1 hc.2.1 hcross
#print axioms sparseIntCollatzCheck22_sound
#print axioms sparseIntMatrixMulUpperCheck22_sound
#print axioms sparseMatrixCollatz_total22
end Spin.Turbo


