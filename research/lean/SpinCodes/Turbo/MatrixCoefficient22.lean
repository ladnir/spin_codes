import SpinCodes.Turbo.ConditionalFractionalMatrix
import Mathlib.Algebra.BigOperators.Fin

/-!
# Exact coefficient recurrence for ordered physical pairs

The historical `ConcreteFixedCoefficient.matrix_sum_ofFn_prod` expands
ordered matrix products over finite input choices. We retain that order
and add the exact total-degree constraint needed by Turbo's regional
hypergeometric calculation. The recurrence is a finite arithmetic checker
interface; numerical matrix entries need their own soundness proofs.
-/
noncomputable section
namespace Spin.Turbo
open Finset

variable {A I : Type*} [Fintype A] [Fintype I] [DecidableEq I]

def orderedMatrixCoefficient (degree : A → ℕ) (W : A → Matrix I I ℝ)
    (n q : ℕ) : Matrix I I ℝ :=
  ∑ xs : Fin n → A,
    if (∑ i, degree (xs i)) = q then (List.ofFn (fun i => W (xs i))).prod else 0

def orderedMatrixDP (degree : A → ℕ) (W : A → Matrix I I ℝ) :
    ℕ → ℕ → Matrix I I ℝ
  | 0, q => if q = 0 then 1 else 0
  | n + 1, q => ∑ a : A,
      if degree a ≤ q then W a * orderedMatrixDP degree W n (q - degree a) else 0

theorem orderedMatrixCoefficient_zero (degree : A → ℕ) (W : A → Matrix I I ℝ)
    (q : ℕ) : orderedMatrixCoefficient degree W 0 q = if q = 0 then 1 else 0 := by
  simp [orderedMatrixCoefficient, eq_comm]

theorem orderedMatrixCoefficient_succ (degree : A → ℕ)
    (W : A → Matrix I I ℝ) (n q : ℕ) :
    orderedMatrixCoefficient degree W (n + 1) q = ∑ a : A,
      if degree a ≤ q then W a * orderedMatrixCoefficient degree W n (q - degree a)
      else 0 := by
  unfold orderedMatrixCoefficient
  rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (n + 1) => A))]
  rw [Fintype.sum_prod_type]
  apply sum_congr rfl
  intro a _
  change (∑ xs : Fin n → A, if (∑ i : Fin (n + 1), degree ((Fin.cons a xs : Fin (n+1) → A) i)) = q
    then (List.ofFn (fun i : Fin (n + 1) => W ((Fin.cons a xs : Fin (n+1) → A) i))).prod else 0) = _
  simp only [Fin.sum_univ_succ, List.ofFn_succ, List.prod_cons,
    Fin.cons_zero, Fin.cons_succ]
  by_cases ha : degree a ≤ q
  · rw [if_pos ha, Finset.mul_sum]
    apply sum_congr rfl
    intro xs _
    have he : degree a + (∑ i, degree (xs i)) = q ↔
        (∑ i, degree (xs i)) = q - degree a := by omega
    simp only [he, mul_ite, mul_zero]
  · rw [if_neg ha]
    apply sum_eq_zero
    intro xs _
    have he : degree a + (∑ i, degree (xs i)) ≠ q := by omega
    simp only [he, ite_false]

theorem orderedMatrixCoefficient_eq_dp (degree : A → ℕ)
    (W : A → Matrix I I ℝ) (n q : ℕ) :
    orderedMatrixCoefficient degree W n q = orderedMatrixDP degree W n q := by
  induction n generalizing q with
  | zero => exact orderedMatrixCoefficient_zero degree W q
  | succ n ih =>
    rw [orderedMatrixCoefficient_succ]
    simp only [orderedMatrixDP, ih]

/-- Two physical steps stay in chronological order before entrywise powers.
The binomial factors are outside the fractional power. -/
def pairedProfileMatrix22 (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (α : ℝ) (ab : Fin 17 × Fin 17) : Matrix (Fin 11) (Fin 11) ℝ :=
  ((Nat.choose 16 ab.1.val : ℝ) * (Nat.choose 16 ab.2.val : ℝ)) •
    matrixEntrywiseRpow (B ab.1 * B ab.2) α

def pairedProfileDegree22 (ab : Fin 17 × Fin 17) : ℕ := ab.1.val + ab.2.val

def pairedRegionDP22 (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (α : ℝ) (pairs selected : ℕ) : Matrix (Fin 11) (Fin 11) ℝ :=
  (1 / (Nat.choose (32 * pairs) selected : ℝ)) •
    orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B α) pairs selected

/-- Exact unnormalized coefficient, then the single regional normalization.
No replacement by independent selections is made. -/
theorem pairedRegionDP22_eq_profile_sum
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (α : ℝ) (pairs selected : ℕ) :
    pairedRegionDP22 B α pairs selected =
      (1 / (Nat.choose (32 * pairs) selected : ℝ)) •
        orderedMatrixCoefficient pairedProfileDegree22
          (pairedProfileMatrix22 B α) pairs selected := by
  rw [orderedMatrixCoefficient_eq_dp]
  rfl

#print axioms orderedMatrixCoefficient_eq_dp
#print axioms pairedRegionDP22_eq_profile_sum
end Spin.Turbo
