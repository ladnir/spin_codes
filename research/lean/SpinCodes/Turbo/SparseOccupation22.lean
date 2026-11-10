import SpinCodes.Turbo.SparseActualShell22
import SpinCodes.Turbo.Occupation22
import SpinCodes.Turbo.SparseMarkerFold22
import Mathlib.Data.Finset.Sort

/-! Sum the actual fixed-active-set first moments over occupation classes.
The sampled setup is unchanged; the factors are L and choose(L,2). -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem turboSetup22_occupation_le_of_activeSet_bound
    (R : ℕ) (points : Fin 8 → Field16) (d q : ℕ) (hq : 0 < q) (U : ℝ)
    (hU : ∀ S : Finset (Fin (R * 16)), S.card = q →
      (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤ U) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) q ω : ℝ)) ≤
          (Nat.choose (R * 16) q : ℝ) * U := by
  rw [turboSetup22_occupation_first_moment_activeSets R points d q hq,
    FinPMF.expect_sum]
  calc
    _ ≤ ∑ S ∈ (Finset.univ : Finset (Fin (R * 16))).powersetCard q, U := by
      apply Finset.sum_le_sum
      intro S hS
      exact hU S (Finset.mem_powersetCard.mp hS).2
    _ = _ := by simp

theorem turboSetup22_occupation_one_le_sparse_shell_family
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (z : ℕ → ℝ) (hz0 : ∀ v, 0 < z v) (hz1 : ∀ v, z v ≤ 1)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ v j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j (z v))) (B v j)) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) 1 ω : ℝ)) ≤
          ((R * 16 : ℕ) : ℝ) * sparseOneShellFamilyBound22 points B R d z := by
  simpa only [Nat.choose_one_right] using
    turboSetup22_occupation_le_of_activeSet_bound R points d 1 (by decide)
      (sparseOneShellFamilyBound22 points B R d z)
      (fun S hS => turboSetup22_one_activeSet_le_sparse_shell_family points S
        (S.orderIsoOfFin hS).toEquiv d z hz0 hz1 B hB)

theorem turboSetup22_occupation_two_le_sparse_shells
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j)) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) 2 ω : ℝ)) ≤
          (Nat.choose (R * 16) 2 : ℝ) * sparseTwoShellBound22 points B R d z := by
  exact turboSetup22_occupation_le_of_activeSet_bound R points d 2 (by decide)
    (sparseTwoShellBound22 points B R d z)
    (fun S hS => turboSetup22_two_activeSet_le_sparse_shells points S
      (S.orderIsoOfFin hS).toEquiv d z hz0 hz1 B hB)

def sparseTwoMarkerBound22 (points : Fin 8 → Field16)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (R d : ℕ) (z t : ℝ) : ℝ :=
  ((∑ k : Fin 11, (orderedMatrixMarker22 sparseTotalRegionDegree22
      (sparseTwoRegionChoice22 B R) t ^ 64) 0 k) *
    (∑ v ∈ Finset.range 65,
      (rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ)) / t ^ v) ^ 2) / z ^ d

theorem sparseTwoShellBound22_le_common_marker
    (points : Fin 8 → Field16) (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (R d : ℕ) (z t : ℝ) (hz : 0 ≤ z) (ht : 0 < t) :
    sparseTwoShellBound22 points B R d z ≤ sparseTwoMarkerBound22 points B R d z t := by
  unfold sparseTwoShellBound22 sparseTwoMarkerBound22
  apply div_le_div_of_nonneg_right _ (pow_nonneg hz d)
  apply sparseTwoWeighted_common_marker_fold_le
    (sparseTwoRegionChoice22 B R) (sparseTwoRegionChoice22_nonneg B hB R) t ht
  intro v _
  apply div_nonneg
  · exact Finset.sum_nonneg fun _ _ => FinPMF.prob_nonneg _ _
  · positivity

theorem turboSetup22_occupation_two_le_common_marker
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) (t : ℝ) (ht : 0 < t)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j))
    (hBn : ∀ a i j, 0 ≤ B a i j) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) 2 ω : ℝ)) ≤
          (Nat.choose (R * 16) 2 : ℝ) * sparseTwoMarkerBound22 points B R d z t := by
  apply (turboSetup22_occupation_two_le_sparse_shells R points d z hz0 hz1 B hB).trans
  exact mul_le_mul_of_nonneg_left
    (sparseTwoShellBound22_le_common_marker points B hBn R d z t hz0.le ht)
    (Nat.cast_nonneg _)
#print axioms turboSetup22_occupation_one_le_sparse_shell_family
#print axioms turboSetup22_occupation_two_le_sparse_shells
end Spin.Turbo


