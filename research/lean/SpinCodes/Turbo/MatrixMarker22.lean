import SpinCodes.Turbo.MatrixPolynomial22
import Mathlib.Tactic.Ring

/-!
# Positive marker bounds for ordered coefficients

Audit: ConcreteFixedCoefficient moves scalar weights through ordered matrix
products; MatrixCoefficient22 and MatrixPolynomial22 preserve factor order.
The marker weights each choice by its degree, then discards an exact-degree
constraint using entrywise nonnegativity. No matrix factors are permuted.
-/
noncomputable section
namespace Spin.Turbo
open Finset
variable {A I : Type*} [Fintype A] [Fintype I] [DecidableEq I]

private theorem markerListNonneg22 (ms : List (Matrix I I ℝ))
    (hm : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j) : ∀ i j, 0 ≤ ms.prod i j := by
  induction ms with
  | nil =>
    intro i j
    simp only [List.prod_nil, Matrix.one_apply]
    split_ifs
    · exact zero_le_one
    · exact le_refl 0
  | cons M ms ih =>
    intro i j
    rw [List.prod_cons, Matrix.mul_apply]
    exact sum_nonneg (fun k _ => mul_nonneg (hm M (List.mem_cons_self) i k)
      (ih (fun N hN => hm N (List.mem_cons_of_mem _ hN)) k j))

private theorem markerScalarProduct22 {n : Nat} (c : Fin n → ℝ) (W : Fin n → Matrix I I ℝ) :
    (List.ofFn (fun i => c i • W i)).prod = (∏ i, c i) • (List.ofFn W).prod := by
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [List.ofFn_succ, List.prod_cons, Fin.prod_univ_succ]
    rw [ih, smul_mul_smul_comm]

private theorem markerSumProducts22 {n : Nat} (W : Fin n → A → Matrix I I ℝ) :
    (∑ xs : Fin n → A, (List.ofFn (fun i => W i (xs i))).prod) =
      (List.ofFn (fun i => ∑ a, W i a)).prod := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (n+1) => A)), Fintype.sum_prod_type]
    change (∑ a : A, ∑ xs : Fin n → A,
      (List.ofFn (fun i => W i ((Fin.cons a xs : Fin (n+1) → A) i))).prod) = _
    simp only [List.ofFn_succ, List.prod_cons, Fin.cons_zero, Fin.cons_succ]
    simp_rw [← Finset.mul_sum, ih]
    rw [← Finset.sum_mul]

def orderedMatrixMarker22 (degree : A → Nat) (W : A → Matrix I I ℝ) (t : ℝ) : Matrix I I ℝ :=
  ∑ a, t ^ degree a • W a

theorem orderedMatrixMarker22_pow (degree : A → Nat) (W : A → Matrix I I ℝ)
    (t : ℝ) (n : Nat) :
    (orderedMatrixMarker22 degree W t) ^ n =
      ∑ xs : Fin n → A, t ^ (∑ i, degree (xs i)) • (List.ofFn (fun i => W (xs i))).prod := by
  unfold orderedMatrixMarker22
  have h := markerSumProducts22 (fun (_ : Fin n) a => t ^ degree a • W a)
  simp only [List.ofFn_const, List.prod_replicate] at h
  rw [← h]
  apply sum_congr rfl
  intro xs _
  rw [markerScalarProduct22, Finset.prod_pow_eq_pow_sum]

theorem orderedMatrixCoefficient22_marker_le (degree : A → Nat) (W : A → Matrix I I ℝ)
    (hW : ∀ a i j, 0 ≤ W a i j) (t : ℝ) (ht : 0 ≤ t) (n q : Nat) (i j : I) :
    t ^ q * orderedMatrixCoefficient degree W n q i j ≤
      (orderedMatrixMarker22 degree W t ^ n) i j := by
  rw [orderedMatrixCoefficient, Matrix.sum_apply, mul_sum, orderedMatrixMarker22_pow,
    Matrix.sum_apply]
  apply sum_le_sum
  intro xs _
  simp only [Matrix.smul_apply, smul_eq_mul]
  split_ifs with he
  · rw [he]
  · simp only [Matrix.zero_apply, mul_zero]
    apply mul_nonneg (pow_nonneg ht _)
    apply markerListNonneg22
    intro M hM
    obtain ⟨a, rfl⟩ := List.mem_ofFn.mp hM
    exact hW _

private theorem markerCoefficientNonneg22 (degree : A → Nat) (W : A → Matrix I I ℝ)
    (hW : ∀ a i j, 0 ≤ W a i j) (n q : Nat) (i j : I) :
    0 ≤ orderedMatrixCoefficient degree W n q i j := by
  rw [orderedMatrixCoefficient, Matrix.sum_apply]
  apply sum_nonneg
  intro xs _
  split_ifs
  · apply markerListNonneg22
    intro M hM
    obtain ⟨a, rfl⟩ := List.mem_ofFn.mp hM
    exact hW _
  · exact le_refl 0

theorem markerPowNonneg22 (M : Matrix I I ℝ) (hM : ∀ i j, 0 ≤ M i j) (n : Nat) :
    ∀ i j, 0 ≤ (M^n) i j := by
  have h := markerListNonneg22 (List.replicate n M) (by
    intro N hN
    have he : N = M := List.eq_of_mem_replicate hN
    subst N
    exact hM)
  simpa only [List.prod_replicate] using h

theorem markerPowMono22 (M N : Matrix I I ℝ)
    (hM : ∀ i j, 0 ≤ M i j) (hN : ∀ i j, 0 ≤ N i j)
    (hMN : ∀ i j, M i j ≤ N i j) (n : Nat) : ∀ i j, (M^n) i j ≤ (N^n) i j := by
  induction n with
  | zero => intro i j; simp
  | succ n ih =>
    intro i j
    rw [pow_succ, pow_succ, Matrix.mul_apply, Matrix.mul_apply]
    exact sum_le_sum (fun k _ => mul_le_mul (ih i k) (hMN k j) (hM k j)
      (markerPowNonneg22 N hN n i k))

theorem markerSmulPow22 (c : ℝ) (M : Matrix I I ℝ) (n : Nat) :
    (c • M)^n = c^n • M^n := by
  induction n with
  | zero => simp
  | succ n ih => rw [pow_succ, ih, smul_mul_smul_comm, ← pow_succ, ← pow_succ]

private theorem pairedProfileMatrix22_marker_nonneg
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) : ∀ ab i j, 0 ≤ pairedProfileMatrix22 B α ab i j := by
  intro ab i j
  change 0 ≤ ((Nat.choose 16 ab.1.val : ℝ) * (Nat.choose 16 ab.2.val : ℝ)) *
    (((B ab.1 * B ab.2) i j) ^ α)
  apply mul_nonneg (mul_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _))
  apply Real.rpow_nonneg
  rw [Matrix.mul_apply]
  exact sum_nonneg (fun k _ => mul_nonneg (hB _ _ _) (hB _ _ _))

/-- One exact regional coefficient is bounded by any positive marker. -/
theorem pairedRegionDP22_marker_le
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (n q : Nat) (t : ℝ) (ht : 0 < t)
    (i j : Fin 11) :
    pairedRegionDP22 B α n q i j ≤
      (1 / (Nat.choose (32*n) q : ℝ) / t^q) *
        (orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B α) t ^ n) i j := by
  rw [pairedRegionDP22_eq_profile_sum]
  simp only [Matrix.smul_apply, smul_eq_mul]
  have h := orderedMatrixCoefficient22_marker_le pairedProfileDegree22
    (pairedProfileMatrix22 B α) (pairedProfileMatrix22_marker_nonneg B α hB) t ht.le n q i j
  have hd : orderedMatrixCoefficient pairedProfileDegree22 (pairedProfileMatrix22 B α) n q i j ≤
      (orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B α) t ^ n) i j / t^q :=
    (le_div_iff₀ (pow_pos ht q)).mpr (by simpa only [mul_comm] using h)
  have hm := mul_le_mul_of_nonneg_left hd
    (one_div_nonneg.mpr (Nat.cast_nonneg (Nat.choose (32*n) q) : (0:ℝ) ≤ _))
  convert hm using 1 <;> ring

/-- The complete boundary state persists through all 64 regions. The scalar
normalizer may be factored out only after the ordered matrix power. -/
theorem pairedRegionDP22_marker_pow64_le
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (n q : Nat) (t : ℝ) (ht : 0 < t)
    (i j : Fin 11) :
    ((pairedRegionDP22 B α n q)^64) i j ≤
      (1 / (Nat.choose (32*n) q : ℝ) / t^q)^64 *
        (orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B α) t ^ (n*64)) i j := by
  let M := orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B α) t
  let c := 1 / (Nat.choose (32*n) q : ℝ) / t^q
  have hW := pairedProfileMatrix22_marker_nonneg B α hB
  have hc : 0 ≤ c := div_nonneg (one_div_nonneg.mpr (Nat.cast_nonneg _)) (pow_nonneg ht.le _)
  have hM : ∀ i j, 0 ≤ M i j := by
    intro i j
    unfold M orderedMatrixMarker22
    rw [Matrix.sum_apply]
    exact sum_nonneg (fun ab _ => mul_nonneg (pow_nonneg ht.le _) (hW ab i j))
  have hDP : ∀ i j, 0 ≤ pairedRegionDP22 B α n q i j := by
    intro i j
    rw [pairedRegionDP22_eq_profile_sum]
    exact mul_nonneg (one_div_nonneg.mpr (Nat.cast_nonneg _))
      (markerCoefficientNonneg22 _ _ hW n q i j)
  have hE : ∀ i j, 0 ≤ (c • M^n) i j := by
    intro i j
    exact mul_nonneg hc (markerPowNonneg22 M hM n i j)
  have hle : ∀ i j, pairedRegionDP22 B α n q i j ≤ (c • M^n) i j :=
    pairedRegionDP22_marker_le B α hB n q t ht
  have h := markerPowMono22 _ _ hDP hE hle 64 i j
  simpa only [markerSmulPow22, ← pow_mul, Matrix.smul_apply, smul_eq_mul, c, M] using h
theorem pairedRegionDP22_marker_terminal_le
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (α : ℝ)
    (hB : ∀ a i j, 0 ≤ B a i j) (n q : Nat) (t : ℝ) (ht : 0 < t)
    (i : Fin 11) :
    (∑ j : Fin 11, ((pairedRegionDP22 B α n q)^64) i j) ≤
      (1 / (Nat.choose (32*n) q : ℝ) / t^q)^64 *
        ∑ j : Fin 11,
          (orderedMatrixMarker22 pairedProfileDegree22 (pairedProfileMatrix22 B α) t ^ (n*64)) i j := by
  rw [mul_sum]
  exact sum_le_sum (fun j _ => pairedRegionDP22_marker_pow64_le B α hB n q t ht i j)
#print axioms orderedMatrixCoefficient22_marker_le
#print axioms pairedRegionDP22_marker_pow64_le
end Spin.Turbo










