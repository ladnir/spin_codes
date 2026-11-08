import Mathlib.Analysis.MeanInequalitiesPow
import Mathlib.Data.Matrix.Mul

/-!
# Fractional powers of chronological nonnegative matrix paths

The retained K18 argument first multiplies each chronological pair, then
takes entrywise powers. Expanding matrix paths and using finite-sum
subadditivity bounds the fractional power of the scalar endpoint. The
coordinates need not describe disjoint events. No routing average or
active-set union is moved through a fractional power in this module.
-/
noncomputable section
namespace Spin.Turbo
open Finset

/-- Finite-sum subadditivity, from Mathlib's two-term real inequality. -/
theorem fractional_sum_le {I : Type*} (s : Finset I) (f : I → ℝ)
    (hf : ∀ i ∈ s, 0 ≤ f i) {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) :
    (∑ i ∈ s, f i) ^ α ≤ ∑ i ∈ s, (f i) ^ α := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [Real.zero_rpow hα0.ne']
  | @insert i s hi ih =>
    rw [sum_insert hi, sum_insert hi]
    exact (Real.rpow_add_le_add_rpow (hf i (mem_insert_self i s))
      (sum_nonneg (fun j hj => hf j (mem_insert_of_mem hj))) hα0.le hα1).trans
      (add_le_add (le_refl _) (ih (fun j hj => hf j (mem_insert_of_mem hj))))

def matrixEntrywiseRpow {I J : Type*} (M : Matrix I J ℝ) (α : ℝ) : Matrix I J ℝ :=
  fun i j => (M i j) ^ α

/-- Fractional subadditivity over one intermediate matrix-path coordinate. -/
theorem matrix_mul_fractional_le {I J K : Type*} [Fintype J]
    (A : Matrix I J ℝ) (B : Matrix J K ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hB : ∀ j k, 0 ≤ B j k)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) (i : I) (k : K) :
    matrixEntrywiseRpow (A * B) α i k ≤
      (matrixEntrywiseRpow A α * matrixEntrywiseRpow B α) i k := by
  change (∑ j, A i j * B j k) ^ α ≤ ∑ j, (A i j) ^ α * (B j k) ^ α
  calc
    _ ≤ ∑ j, (A i j * B j k) ^ α :=
      fractional_sum_le univ _ (fun j _ => mul_nonneg (hA i j) (hB j k)) hα0 hα1
    _ = _ := sum_congr rfl (fun j _ => Real.mul_rpow (hA i j) (hB j k))

private theorem matrix_list_prod_nonneg {I : Type*} [Fintype I] [DecidableEq I]
    (ms : List (Matrix I I ℝ)) (hms : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j) :
    ∀ i j, 0 ≤ ms.prod i j := by
  induction ms with
  | nil => intro i j; simp [Matrix.one_apply]; split_ifs <;> norm_num
  | cons M ms ih =>
    intro i j
    rw [List.prod_cons, Matrix.mul_apply]
    exact sum_nonneg (fun k _ => mul_nonneg
      (hms M (List.mem_cons_self) i k)
      (ih (fun N hN => hms N (List.mem_cons_of_mem M hN)) k j))

/-- The powered chronological scalar paths are bounded by a chronological
product of entrywise-powered blocks. -/
theorem matrix_list_prod_fractional_le {I : Type*} [Fintype I] [DecidableEq I]
    (ms : List (Matrix I I ℝ)) (hms : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) (i j : I) :
    (ms.prod i j) ^ α ≤
      (ms.map (fun M => matrixEntrywiseRpow M α)).prod i j := by
  induction ms generalizing i j with
  | nil => simp [Matrix.one_apply, Real.zero_rpow hα0.ne']
  | cons M ms ih =>
    have hM := hms M List.mem_cons_self
    have htail : ∀ N ∈ ms, ∀ a b, 0 ≤ N a b :=
      fun N hN => hms N (List.mem_cons_of_mem M hN)
    simp only [List.prod_cons, List.map_cons]
    calc
      _ ≤ (matrixEntrywiseRpow M α * matrixEntrywiseRpow ms.prod α) i j :=
        matrix_mul_fractional_le M ms.prod hM
          (matrix_list_prod_nonneg ms htail) hα0 hα1 i j
      _ ≤ _ := by
        rw [Matrix.mul_apply, Matrix.mul_apply]
        exact sum_le_sum (fun k _ => mul_le_mul_of_nonneg_left
          (ih htail k j) (Real.rpow_nonneg (hM i k) α))

/-- Initial coordinate `i`, summed terminal coordinates, and state carried
through every block. This is the `e₀ product(Q) 1` endpoint. -/
theorem matrix_list_row_sum_fractional_le {I : Type*} [Fintype I] [DecidableEq I]
    (ms : List (Matrix I I ℝ)) (hms : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) (i : I) :
    (∑ j, ms.prod i j) ^ α ≤
      ∑ j, (ms.map (fun M => matrixEntrywiseRpow M α)).prod i j := by
  exact (fractional_sum_le univ _
    (fun j _ => matrix_list_prod_nonneg ms hms i j) hα0 hα1).trans
    (sum_le_sum (fun j _ => matrix_list_prod_fractional_le ms hms hα0 hα1 i j))

private theorem chronological_pair_product {I : Type*} [Fintype I] [DecidableEq I]
    (pairs : List (Matrix I I ℝ × Matrix I I ℝ)) :
    (pairs.flatMap (fun p => [p.1, p.2])).prod =
      (pairs.map (fun p => p.1 * p.2)).prod := by
  induction pairs with
  | nil => rfl
  | cons p pairs ih =>
    simp only [List.flatMap_cons, List.prod_append, List.prod_cons, List.prod_nil,
      mul_one, List.map_cons, ih, mul_assoc]

/-- K18's two-step order: multiply each actual chronological pair first,
then take entrywise powers, and retain the full matrix product. -/
theorem chronological_pairs_row_fractional_le
    {I : Type*} [Fintype I] [DecidableEq I]
    (pairs : List (Matrix I I ℝ × Matrix I I ℝ))
    (hpairs : ∀ p ∈ pairs, (∀ i j, 0 ≤ p.1 i j) ∧ (∀ i j, 0 ≤ p.2 i j))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) (i : I) :
    (∑ j, (pairs.flatMap (fun p => [p.1, p.2])).prod i j) ^ α ≤
      ∑ j, (pairs.map (fun p => matrixEntrywiseRpow (p.1 * p.2) α)).prod i j := by
  rw [chronological_pair_product]
  have h := matrix_list_row_sum_fractional_le
    (pairs.map (fun p => p.1 * p.2)) (by
      intro M hM a b
      obtain ⟨p, hp, rfl⟩ := List.mem_map.mp hM
      rw [Matrix.mul_apply]
      exact sum_nonneg (fun k _ => mul_nonneg
        ((hpairs p hp).1 a k) ((hpairs p hp).2 k b))) hα0 hα1 i
  simpa only [List.map_map, Function.comp_def] using h

#print axioms matrix_list_row_sum_fractional_le
#print axioms chronological_pairs_row_fractional_le
end Spin.Turbo
