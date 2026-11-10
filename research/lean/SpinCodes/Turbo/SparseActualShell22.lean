import SpinCodes.Turbo.SparseActualEnvelope22
import SpinCodes.Turbo.SparseShellFold22

/-! The actual shared setup's one-group first moment reaches every current
RS8 shell. Regional kernels are normalized before the ordered coefficient;
the outer support binomial and Chernoff denominator occur only once. -/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

def sparseOneSupportEquiv22 :
    Finset (Fin 64) ≃ (Fin 1 → Finset (Fin 64)) where
  toFun T := fun _ => T
  invFun A := A 0
  left_inv _ := rfl
  right_inv A := by funext i; fin_cases i; rfl

def sparseOneShellFamilyBound22 (points : Fin 8 → Field16)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (R d : ℕ) (z : ℕ → ℝ) : ℝ :=
  ∑ v ∈ Finset.range 65,
    ((rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ)) / (z v) ^ d) *
      ∑ k : Fin 11, orderedMatrixCoefficient sparseBoolDegree22
        (sparseOneRegionChoice22 (B v) R) 64 v 0 k

/-- The shell selects its fugacity before the unconditional first-moment
sum. This is the interface required by the retained q1 tilt family. -/
theorem turboSetup22_one_activeSet_le_sparse_shell_family
    {R : ℕ} (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin 1 ≃ {i : Fin (R * 16) // i ∈ S}) (d : ℕ)
    (z : ℕ → ℝ) (hz0 : ∀ v, 0 < z v) (hz1 : ∀ v, z v ≤ 1)
    (B : ℕ → Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ v j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j (z v))) (B v j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      sparseOneShellFamilyBound22 points B R d z := by
  let Z := fun supports : Fin (R * 16) → Finset (Fin 64) => z (supports (e 0).val).card
  let BB := fun supports : Fin (R * 16) → Finset (Fin 64) => B (supports (e 0).val).card
  have h := turboSetup22_activeSet_le_sparse_active_profiles_family points S e d Z
    (fun supports => hz0 _) (fun supports => hz1 _) BB (fun supports => hB _)
  dsimp only [Z, BB] at h
  simp_rw [sparseEnumeratedSupport22_apply, sparseEnumeratedSupport22_one_DP] at h
  apply h.trans_eq
  rw [← Equiv.sum_comp sparseOneSupportEquiv22]
  simp only [sparseOneSupportEquiv22, Equiv.coe_fn_mk, Fin.prod_univ_one]
  simp_rw [← mul_div_assoc, ← div_mul_eq_mul_div]
  exact sparseOneWeightedSupportFold22_family_mass (0 : Fin 11)
    (fun v => (rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ)) / (z v) ^ d)
    (fun v => sparseOneRegionChoice22 (B v) R)
def sparseTwoShellBound22 (points : Fin 8 → Field16)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (R d : ℕ) (z : ℝ) : ℝ :=
  (∑ u ∈ Finset.range 65, ∑ v ∈ Finset.range 65,
    ((rs8ExpectedPacketCount points u / (Nat.choose 64 u : ℝ)) *
      (rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ))) *
      ∑ k : Fin 11, sparseTwoShellCoefficient22
        (sparseTwoRegionChoice22 B R) u v 0 k) / z ^ d

/-- Both selected groups retain one shared regional permutation. The full
ordered shell sum includes diagonal and both off-diagonal orientations. -/
theorem turboSetup22_two_activeSet_le_sparse_shells
    {R : ℕ} (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin 2 ≃ {i : Fin (R * 16) // i ∈ S}) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ j, physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) (B j)) :
    (turboOuterLaw22 R).expect (fun seed => ∑ msg ∈ turboActiveMessages22 S,
      (turboSetup22 R points).qd d (turboOuterWord22 points seed msg)) ≤
      sparseTwoShellBound22 points B R d z := by
  have h := turboSetup22_activeSet_le_sparse_active_profiles points S e d z hz0 hz1 B hB
  simp_rw [sparseEnumeratedSupport22_two_DP] at h
  apply h.trans_eq
  rw [← Equiv.sum_comp (finTwoArrowEquiv (Finset (Fin 64))).symm]
  simp only [Fin.prod_univ_two, finTwoArrowEquiv_symm_apply,
    Matrix.cons_val_zero, Matrix.cons_val_one]
  -- Reassociate only the outer Chernoff division; the product weight
  -- (a/C₁)*(b/C₂) must keep the exact shape of the shell-fold lemma.
  refine (Finset.sum_congr rfl fun AB _ => (mul_div_assoc _ _ _).symm).trans ?_
  rw [← Finset.sum_div]
  unfold sparseTwoShellBound22
  exact congrArg (fun value : ℝ => value / z ^ d)
    (sparseTwoWeightedSupportFold22_mass (0 : Fin 11)
      (fun v => rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ))
      (sparseTwoRegionChoice22 B R))
#print axioms turboSetup22_one_activeSet_le_sparse_shell_family
#print axioms turboSetup22_two_activeSet_le_sparse_shells
end Spin.Turbo




