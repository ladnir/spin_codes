import SpinCodes.Turbo.SparseActiveProfiles22
import SpinCodes.Turbo.SparsePhysicalKernel22
import SpinCodes.Turbo.SparseMarker22

/-! For one/two actual active groups, every region uses only occupancy0,1,2
under its shared permutation. These are the exact normalized physical DPs;
no canonical slot choices or independent regional permutations are needed. -/
noncomputable section
namespace Spin.Turbo
open Finset

def sparseRegionalDP22 {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R selected : ℕ) : Matrix I I ℝ :=
  (1 / (Nat.choose (R * 16) selected : ℝ)) •
    orderedMatrixDP Fin.val (fun j : Fin 17 => (Nat.choose 16 j.val : ℝ) • B j) R selected

def sparseOneRegionChoice22 {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) (bit : Bool) : Matrix I I ℝ :=
  if bit then sparseRegionalDP22 B R 1 else sparseRegionalDP22 B R 0

def sparseTwoRegionChoice22 {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) : Bool × Bool → Matrix I I ℝ :=
  sparseJointRegionChoice22 (sparseRegionalDP22 B R 0)
    (sparseRegionalDP22 B R 1) (sparseRegionalDP22 B R 2)

theorem sparsePhysicalRegionMatrix22_eq_occupancy {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (R : ℕ) (T : Finset (Fin (R * 16))) :
    sparsePhysicalRegionMatrix22 B R T = sparseRegionalDP22 B R T.card :=
  sparsePhysicalRegionMatrix22_eq_dp B R T

theorem sparseEnumeratedSupport22_one_region {I : Type*} [Fintype I] [DecidableEq I]
    {R : ℕ} (B : Fin 17 → Matrix I I ℝ) (S : Finset (Fin (R * 16)))
    (e : Fin 1 ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin 1 → Finset (Fin 64)) (region : Fin 64) :
    sparsePhysicalRegionMatrix22 B R
      (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region) =
      sparseOneRegionChoice22 B R (decide (region ∈ A 0)) := by
  rw [sparsePhysicalRegionMatrix22_eq_occupancy]
  have hc := sparseEnumeratedSupport22_region_card S e A region
  simp only [Fin.sum_univ_one] at hc
  rw [hc]
  by_cases h : region ∈ A 0 <;> simp [sparseOneRegionChoice22, h]

theorem sparseEnumeratedSupport22_two_region {I : Type*} [Fintype I] [DecidableEq I]
    {R : ℕ} (B : Fin 17 → Matrix I I ℝ) (S : Finset (Fin (R * 16)))
    (e : Fin 2 ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin 2 → Finset (Fin 64)) (region : Fin 64) :
    sparsePhysicalRegionMatrix22 B R
      (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region) =
      sparseTwoRegionChoice22 B R (decide (region ∈ A 0), decide (region ∈ A 1)) := by
  rw [sparsePhysicalRegionMatrix22_eq_occupancy]
  have hc := sparseEnumeratedSupport22_region_card S e A region
  simp only [Fin.sum_univ_two] at hc
  rw [hc]
  by_cases h0 : region ∈ A 0 <;> by_cases h1 : region ∈ A 1 <;>
    simp [sparseTwoRegionChoice22, sparseJointRegionChoice22, h0, h1]

theorem sparseEnumeratedSupport22_one_DP {I : Type*} [Fintype I] [DecidableEq I]
    {R : ℕ} (B : Fin 17 → Matrix I I ℝ) (S : Finset (Fin (R * 16)))
    (e : Fin 1 ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin 1 → Finset (Fin 64)) (region : Fin 64) :
    sparseRegionalDP22 B R
      (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region).card =
      sparseOneRegionChoice22 B R (decide (region ∈ A 0)) := by
  simpa only [sparsePhysicalRegionMatrix22_eq_occupancy] using
    sparseEnumeratedSupport22_one_region B S e A region

theorem sparseEnumeratedSupport22_two_DP {I : Type*} [Fintype I] [DecidableEq I]
    {R : ℕ} (B : Fin 17 → Matrix I I ℝ) (S : Finset (Fin (R * 16)))
    (e : Fin 2 ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin 2 → Finset (Fin 64)) (region : Fin 64) :
    sparseRegionalDP22 B R
      (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region).card =
      sparseTwoRegionChoice22 B R (decide (region ∈ A 0), decide (region ∈ A 1)) := by
  simpa only [sparsePhysicalRegionMatrix22_eq_occupancy] using
    sparseEnumeratedSupport22_two_region B S e A region
private theorem sparseActiveListNonneg22 {I : Type*} [Fintype I] [DecidableEq I]
    (ms : List (Matrix I I ℝ)) (hm : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j) :
    ∀ i j, 0 ≤ ms.prod i j := by
  induction ms with
  | nil => intro i j; simp only [List.prod_nil, Matrix.one_apply]; split_ifs <;> norm_num
  | cons M ms ih =>
    intro i j
    rw [List.prod_cons, Matrix.mul_apply]
    exact Finset.sum_nonneg fun k _ => mul_nonneg (hm M List.mem_cons_self i k)
      (ih (fun N hN => hm N (List.mem_cons_of_mem M hN)) k j)

theorem sparseRegionalDP22_nonneg {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (R selected : ℕ) (i j : I) : 0 ≤ sparseRegionalDP22 B R selected i j := by
  unfold sparseRegionalDP22
  rw [← orderedMatrixCoefficient_eq_dp]
  change 0 ≤ (1 / (Nat.choose (R * 16) selected : ℝ)) * _
  apply mul_nonneg (one_div_nonneg.mpr (Nat.cast_nonneg _))
  unfold orderedMatrixCoefficient
  rw [Matrix.sum_apply]
  apply Finset.sum_nonneg
  intro counts _
  split_ifs
  · apply sparseActiveListNonneg22
    intro M hM
    obtain ⟨step, rfl⟩ := List.mem_ofFn.mp hM
    intro a b
    exact mul_nonneg (Nat.cast_nonneg _) (hB _ _ _)
  · exact le_rfl
theorem sparseTwoRegionChoice22_nonneg {I : Type*} [Fintype I] [DecidableEq I]
    (B : Fin 17 → Matrix I I ℝ) (hB : ∀ a i j, 0 ≤ B a i j)
    (R : ℕ) (bits : Bool × Bool) (i j : I) :
    0 ≤ sparseTwoRegionChoice22 B R bits i j := by
  rcases bits with ⟨a, b⟩
  cases a <;> cases b <;>
    exact sparseRegionalDP22_nonneg B hB R _ i j
#print axioms sparseRegionalDP22_nonneg
#print axioms sparseEnumeratedSupport22_one_region
#print axioms sparseEnumeratedSupport22_two_region
end Spin.Turbo



