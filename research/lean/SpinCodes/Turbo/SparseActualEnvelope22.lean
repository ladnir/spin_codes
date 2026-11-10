import SpinCodes.Turbo.SparseActualMoment22
import SpinCodes.Turbo.SparseRegionalEnvelope22
import SpinCodes.Turbo.SparseActiveKernel22

/-! Original active-set first moments reach the normalized eleven-coordinate
physical DP, and then the exact one/two-group support choices. This uses
unconditional expected counts throughout. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem turboSetup22_activeSet_le_sparse_physicalDP
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        sparseSupportProfileCount22 R points S supports *
          ((∑ k : Fin 11, (List.ofFn (fun region => sparseRegionalDP22 B R
            (sparseRegionalSupport22 supports region).card)).prod 0 k) / z ^ d) := by
  apply (turboSetup22_activeSet_le_sparse_kernel R points S d z hz0 hz1).trans
  apply Finset.sum_le_sum
  intro supports _
  apply mul_le_mul_of_nonneg_left _ (sparseSupportProfileCount22_nonneg R points S supports)
  apply div_le_div_of_nonneg_right _ (pow_nonneg hz0.le d)
  have h := sparseExactRegional_product_mass_le R z
    (fun region => sparseRegionalSupport22 supports region) B hB
  simpa only [sparsePhysicalRegionMatrix22_eq_occupancy] using h

theorem turboSetup22_activeSet_le_sparse_active_profiles
    {R q : ℕ} (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S}) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ A : Fin q → Finset (Fin 64),
        (∏ a : Fin q, rs8ExpectedPacketCount points (A a).card /
          (Nat.choose 64 (A a).card : ℝ)) *
          ((∑ k : Fin 11, (List.ofFn (fun region => sparseRegionalDP22 B R
            (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region).card)).prod 0 k) / z ^ d) := by
  let F := fun supports : Fin (R * 16) → Finset (Fin 64) =>
    (∑ k : Fin 11, (List.ofFn (fun region => sparseRegionalDP22 B R
      (sparseRegionalSupport22 supports region).card)).prod 0 k) / z ^ d
  calc
    _ ≤ ∑ supports : Fin (R * 16) → Finset (Fin 64),
        (∏ i, sparseSupportCount22 points (i ∈ S) (supports i)) * F supports := by
      simpa only [sparseSupportProfileCount22, F] using
        turboSetup22_activeSet_le_sparse_physicalDP R points S d z hz0 hz1 B hB
    _ = _ := sparseSupportProfileCount22_reindex points S e F

theorem turboSetup22_activeSet_le_sparse_physicalDP_family
    (R : ℕ) (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : (Fin (R * 16) → Finset (Fin 64)) → ℝ)
    (hz0 : ∀ supports, 0 < z supports) (hz1 : ∀ supports, z supports ≤ 1)
    (B : (Fin (R * 16) → Finset (Fin 64)) → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ supports j, physicalRowSound22
      (Matrix.of (turboPhysicalKernel22 j (z supports))) (B supports j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ supports : Fin (R * 16) → Finset (Fin 64),
        sparseSupportProfileCount22 R points S supports *
          ((∑ k : Fin 11, (List.ofFn (fun region => sparseRegionalDP22 (B supports) R
            (sparseRegionalSupport22 supports region).card)).prod 0 k) / (z supports) ^ d) := by
  apply (turboSetup22_activeSet_le_sparse_kernel_family R points S d z hz0 hz1).trans
  apply Finset.sum_le_sum
  intro supports _
  apply mul_le_mul_of_nonneg_left _ (sparseSupportProfileCount22_nonneg R points S supports)
  apply div_le_div_of_nonneg_right _ (pow_nonneg (hz0 supports).le d)
  have h := sparseExactRegional_product_mass_le R (z supports)
    (fun region => sparseRegionalSupport22 supports region) (B supports) (hB supports)
  simpa only [sparsePhysicalRegionMatrix22_eq_occupancy] using h
theorem turboSetup22_activeSet_le_sparse_active_profiles_family
    {R q : ℕ} (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S}) (d : ℕ)
    (z : (Fin (R * 16) → Finset (Fin 64)) → ℝ)
    (hz0 : ∀ supports, 0 < z supports) (hz1 : ∀ supports, z supports ≤ 1)
    (B : (Fin (R * 16) → Finset (Fin 64)) → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ supports j, physicalRowSound22
      (Matrix.of (turboPhysicalKernel22 j (z supports))) (B supports j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      ∑ A : Fin q → Finset (Fin 64),
        (∏ a : Fin q, rs8ExpectedPacketCount points (A a).card /
          (Nat.choose 64 (A a).card : ℝ)) *
          ((∑ k : Fin 11, (List.ofFn (fun region =>
            sparseRegionalDP22 (B (sparseEnumeratedSupport22 S e A)) R
              (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region).card)).prod 0 k) /
                (z (sparseEnumeratedSupport22 S e A)) ^ d) := by
  let F := fun supports : Fin (R * 16) → Finset (Fin 64) =>
    (∑ k : Fin 11, (List.ofFn (fun region => sparseRegionalDP22 (B supports) R
      (sparseRegionalSupport22 supports region).card)).prod 0 k) / (z supports) ^ d
  calc
    _ ≤ ∑ supports : Fin (R * 16) → Finset (Fin 64),
        (∏ i, sparseSupportCount22 points (i ∈ S) (supports i)) * F supports := by
      simpa only [sparseSupportProfileCount22, F] using
        turboSetup22_activeSet_le_sparse_physicalDP_family R points S d z hz0 hz1 B hB
    _ = _ := sparseSupportProfileCount22_reindex points S e F
#print axioms turboSetup22_activeSet_le_sparse_physicalDP
#print axioms turboSetup22_activeSet_le_sparse_active_profiles
end Spin.Turbo


