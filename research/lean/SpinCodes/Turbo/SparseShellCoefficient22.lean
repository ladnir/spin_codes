import SpinCodes.Turbo.MatrixPolynomial22

/-! The two sparse support degrees reuse the ordered univariate checker.
Base 65 has no carry because each of the 64 regions contributes at most
one packet per selected group. This is a semantic encoding; producers may
retain their efficient two-dimensional arrays. -/
noncomputable section
namespace Spin.Turbo
open Finset

def sparseBoolDegree22 (b : Bool) : ℕ := if b then 1 else 0

def sparseTwoRegionDegree22 (b : Bool × Bool) : ℕ :=
  sparseBoolDegree22 b.1 + 65 * sparseBoolDegree22 b.2

theorem sparseBoolDegree22_sum_le (n : ℕ) (b : Fin n → Bool) :
    (∑ i, sparseBoolDegree22 (b i)) ≤ n := by
  calc
    _ ≤ ∑ _i : Fin n, (1 : ℕ) := by
      apply Finset.sum_le_sum
      intro i _
      cases b i <;> decide
    _ = n := by simp

theorem sparseTwoRegionDegree22_sum (n : ℕ) (b : Fin n → Bool × Bool) :
    (∑ i, sparseTwoRegionDegree22 (b i)) =
      (∑ i, sparseBoolDegree22 (b i).1) +
        65 * (∑ i, sparseBoolDegree22 (b i).2) := by
  simp only [sparseTwoRegionDegree22, Finset.sum_add_distrib, Finset.mul_sum]

theorem sparseTwoDegree22_no_carry (a b u v : ℕ) (ha : a ≤ 64) (hu : u ≤ 64) :
    a + 65 * b = u + 65 * v ↔ a = u ∧ b = v := by
  omega

variable {I : Type*} [Fintype I] [DecidableEq I]

/-- Ordered complete regional choices, constrained separately for both
selected groups. W is a joint regional kernel, never a product of two
independently permuted group kernels. -/
def sparseTwoShellCoefficient22 (W : Bool × Bool → Matrix I I ℝ) (u v : ℕ) :
    Matrix I I ℝ :=
  ∑ bits : Fin 64 → Bool × Bool,
    if (∑ i, sparseBoolDegree22 (bits i).1) = u ∧
        (∑ i, sparseBoolDegree22 (bits i).2) = v then
      (List.ofFn (fun i => W (bits i))).prod else 0

theorem sparseTwoShellCoefficient22_eq_ordered
    (W : Bool × Bool → Matrix I I ℝ) (u v : ℕ) (hu : u ≤ 64) :
    sparseTwoShellCoefficient22 W u v =
      orderedMatrixCoefficient sparseTwoRegionDegree22 W 64 (u + 65 * v) := by
  unfold sparseTwoShellCoefficient22 orderedMatrixCoefficient
  apply Finset.sum_congr rfl
  intro bits _
  rw [sparseTwoRegionDegree22_sum]
  have ha := sparseBoolDegree22_sum_le 64 (fun i => (bits i).1)
  simp only [sparseTwoDegree22_no_carry _ _ _ _ ha hu]

theorem sparseTwoShellCoefficient22_eq_polynomial
    (W : Bool × Bool → Matrix I I ℝ) (u v : ℕ) (hu : u ≤ 64) :
    sparseTwoShellCoefficient22 W u v =
      ((orderedMatrixPolynomial sparseTwoRegionDegree22 W) ^ 64).coeff (u + 65 * v) := by
  rw [sparseTwoShellCoefficient22_eq_ordered W u v hu,
    orderedMatrixPolynomial_pow_coeff, orderedMatrixCoefficient_eq_dp]

#print axioms sparseTwoDegree22_no_carry
#print axioms sparseTwoShellCoefficient22_eq_polynomial
end Spin.Turbo

