import SpinCodes.Turbo.SparseSupportMeasure22
import SpinCodes.Turbo.SparseRegionKernel22

/-! The exact Cv support-class measure is connected to the original shared
outer setup and the continuous regional kernel, before finite shell checks. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

def sparseSupportProfileCount22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (supports : Fin (R * 16) → Finset (Fin 64)) : ℝ :=
  ∏ i, sparseSupportCount22 points (i ∈ S) (supports i)

theorem sparseSupportCount22_nonneg (points : Fin 8 → Field16)
    (active : Prop) [Decidable active] (T : Finset (Fin 64)) :
    0 ≤ sparseSupportCount22 points active T := by
  unfold sparseSupportCount22
  split
  · apply div_nonneg
    · exact Finset.sum_nonneg fun _ _ => FinPMF.prob_nonneg _ _
    · positivity
  · split <;> norm_num

theorem sparseSupportProfileCount22_nonneg (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (supports : Fin (R * 16) → Finset (Fin 64)) :
    0 ≤ sparseSupportProfileCount22 R points S supports :=
  Finset.prod_nonneg fun _ _ => sparseSupportCount22_nonneg _ _ _

/-- Exact original-law first moment, disintegrated into complete-valued
support classes. Selected messages still share the original outer seed. -/
theorem turboSetup22_activeSet_expect_eq_support_measure
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) (d : ℕ) :
    (turboOuterLaw22 R).expect (fun seed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) =
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        sparseSupportProfileCount22 R points S supports *
          (sparseSupportArrayLaw22 R supports).expect
            (fun words => turboRemainingFailure22 words d) := by
  rw [turboSetup22_activeSet_expect_eq_shuffled_count]
  simp_rw [turboShuffledTupleCount22_eq_support_measure, Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro supports _
  unfold sparseSupportProfileCount22 FinPMF.expect
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro words _
  ring

/-- Actual sparse contribution bounded by exact shared-permutation regional
kernels; this retains all support cardinalities and complete packet labels. -/
theorem turboSetup22_activeSet_le_sparse_kernel
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (turboOuterLaw22 R).expect (fun seed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        sparseSupportProfileCount22 R points S supports *
          ((∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R
            (sparseRegionalSupport22 supports region) z)).prod 0 q) / z ^ d) := by
  rw [turboSetup22_activeSet_expect_eq_support_measure]
  apply Finset.sum_le_sum
  intro supports _
  exact mul_le_mul_of_nonneg_left
    (sparseSupportFailure22_le_kernel R supports d z hz0 hz1)
    (sparseSupportProfileCount22_nonneg R points S supports)

/-- A separate positive Chernoff fugacity may be fixed for each complete
support profile before summing the unconditional expected measure. -/
theorem turboSetup22_activeSet_le_sparse_kernel_family
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : (Fin (R * 16) → Finset (Fin 64)) → ℝ)
    (hz0 : ∀ supports, 0 < z supports) (hz1 : ∀ supports, z supports ≤ 1) :
    (turboOuterLaw22 R).expect (fun seed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        sparseSupportProfileCount22 R points S supports *
          ((∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R
            (sparseRegionalSupport22 supports region) (z supports))).prod 0 q) /
              (z supports) ^ d) := by
  rw [turboSetup22_activeSet_expect_eq_support_measure]
  apply Finset.sum_le_sum
  intro supports _
  exact mul_le_mul_of_nonneg_left
    (sparseSupportFailure22_le_kernel R supports d (z supports) (hz0 supports) (hz1 supports))
    (sparseSupportProfileCount22_nonneg R points S supports)
#print axioms turboSetup22_activeSet_expect_eq_support_measure
#print axioms turboSetup22_activeSet_le_sparse_kernel_family
#print axioms turboSetup22_activeSet_le_sparse_kernel
end Spin.Turbo


