import SpinCodes.Turbo.MatrixScaledDefs22
import SpinCodes.Numeric.Fixed
import Mathlib.Data.Matrix.Mul

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix Finset

def MatrixFixMem22 (a : List (List Fix)) (M : Matrix (Fin 11) (Fin 11) ℝ) : Prop :=
  ∀ i j, Mem (matrixFixEntry22 a i.val j.val) (M i j)

private theorem matrixFixGetOfFn22 {n : Nat} {α : Type*} (f : Fin n → α)
    (d : α) (i : Fin n) : (List.ofFn f).getD i.val d = f i := by
  rw [List.getD_eq_getElem _ _ (n := i.val) (by simpa only [List.length_ofFn] using i.isLt),
    List.getElem_ofFn]

private theorem matrixFixSumMem22 {n : Nat} (a : Fin n → Fix) (x : Fin n → ℝ)
    (h : ∀ i, Mem (a i) (x i)) :
    Mem (matrixFixSum22 (List.ofFn a)) (∑ i, x i) := by
  induction n with
  | zero => simp [matrixFixSum22, Mem, Fix.zero]
  | succ n ih =>
    rw [List.ofFn_succ, matrixFixSum22, Fin.sum_univ_succ]
    exact add_mem (h 0) (ih (fun i => a i.succ) (fun i => x i.succ) (fun i => h i.succ))

theorem matrixFixMul22_mem {a b : List (List Fix)}
    {A B : Matrix (Fin 11) (Fin 11) ℝ} (ha : MatrixFixMem22 a A) (hb : MatrixFixMem22 b B) :
    MatrixFixMem22 (matrixFixMul22 a b) (A * B) := by
  intro i j
  simp only [matrixFixEntry22, matrixFixMul22, matrixFixGetOfFn22]
  rw [Matrix.mul_apply]
  exact matrixFixSumMem22 _ _ (fun k => mul_mem (ha i k) (hb k j))

theorem matrixFixScale22_mem {a : List (List Fix)} {A : Matrix (Fin 11) (Fin 11) ℝ}
    (ha : MatrixFixMem22 a A) (num den : Int) (hd : 0 < den) :
    MatrixFixMem22 (matrixFixScale22 num den a) (((num : ℝ)/(den : ℝ)) • A) := by
  intro i j
  simp only [matrixFixEntry22, matrixFixScale22, matrixFixGetOfFn22,
    Matrix.smul_apply, smul_eq_mul]
  exact mul_mem (ofFrac_mem hd) (ha i j)

theorem matrixFixSquareStep22_mem {a : List (List Fix)}
    {A : Matrix (Fin 11) (Fin 11) ℝ} (ha : MatrixFixMem22 a A)
    (num den : Int) (hd : 0 < den) :
    MatrixFixMem22 (matrixFixSquareStep22 num den a)
      (((num : ℝ)/(den : ℝ)) • (A * A)) :=
  matrixFixScale22_mem (matrixFixMul22_mem ha ha) num den hd

def matrixScaledRealChain22 (A : Matrix (Fin 11) (Fin 11) ℝ) :
    List (Int × Int) → Matrix (Fin 11) (Fin 11) ℝ
  | [] => A
  | (num,den) :: ss => matrixScaledRealChain22 (((num : ℝ)/(den : ℝ)) • (A*A)) ss

theorem matrixFixSquareChain22_mem {a : List (List Fix)}
    {A : Matrix (Fin 11) (Fin 11) ℝ} (ha : MatrixFixMem22 a A)
    (ss : List (Int × Int)) (hd : ∀ s ∈ ss, 0 < s.2) :
    MatrixFixMem22 (matrixFixSquareChain22 a ss) (matrixScaledRealChain22 A ss) := by
  induction ss generalizing a A with
  | nil => exact ha
  | cons s ss ih =>
    exact ih (matrixFixSquareStep22_mem ha s.1 s.2 (hd s (by simp)))
      (fun t ht => hd t (by simp [ht]))

theorem matrixFixTerminal22_mem {a : List (List Fix)}
    {A : Matrix (Fin 11) (Fin 11) ℝ} (ha : MatrixFixMem22 a A) :
    Mem (matrixFixTerminal22 a) (∑ j, A 0 j) :=
  matrixFixSumMem22 _ _ (ha 0)

#print axioms matrixFixSquareChain22_mem
end Spin.Turbo

namespace Spin.Turbo

def matrixScaledFactor22 : List (Int × Int) → ℝ
  | [] => 1
  | (num,den) :: ss => matrixScaledFactor22 ss *
      ((num : ℝ)/(den : ℝ))^(2^ss.length)

private theorem matrixScaledSmulPow22 (c : ℝ) (A : Matrix (Fin 11) (Fin 11) ℝ) (n : Nat) :
    (c • A)^n = c^n • A^n := by
  induction n with
  | zero => simp
  | succ n ih => rw [pow_succ, ih, smul_mul_smul_comm, ← pow_succ, ← pow_succ]

theorem matrixScaledRealChain22_eq (A : Matrix (Fin 11) (Fin 11) ℝ)
    (ss : List (Int × Int)) :
    matrixScaledRealChain22 A ss = matrixScaledFactor22 ss • A^(2^ss.length) := by
  induction ss generalizing A with
  | nil => simp [matrixScaledRealChain22, matrixScaledFactor22]
  | cons s ss ih =>
    rw [matrixScaledRealChain22, ih, matrixScaledSmulPow22, smul_smul,
      ← pow_two, ← pow_mul]
    simp only [matrixScaledFactor22, List.length_cons, pow_succ]
    rw [Nat.mul_comm 2]

 theorem matrixScaledFactor22_pos (ss : List (Int × Int))
    (h : ∀ s ∈ ss, 0 < s.1 ∧ 0 < s.2) : 0 < matrixScaledFactor22 ss := by
  induction ss with
  | nil => exact zero_lt_one
  | cons s ss ih =>
    have hs := h s (by simp)
    exact mul_pos (ih (fun t ht => h t (by simp [ht])))
      (pow_pos (div_pos (by exact_mod_cast hs.1) (by exact_mod_cast hs.2)) _)

#print axioms matrixScaledRealChain22_eq
end Spin.Turbo

namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

theorem matrixFixEnclosesCheck22_sound {outer inner : List (List Fix)}
    (h : matrixFixEnclosesCheck22 outer inner = true)
    {M : Matrix (Fin 11) (Fin 11) ℝ} (hm : MatrixFixMem22 inner M) :
    MatrixFixMem22 outer M := by
  intro i j
  have hi := List.all_eq_true.mp h i.val (List.mem_range.mpr i.isLt)
  have hij := List.all_eq_true.mp hi j.val (List.mem_range.mpr j.isLt)
  have hc : (matrixFixEntry22 outer i.val j.val).lo ≤
      (matrixFixEntry22 inner i.val j.val).lo ∧
      (matrixFixEntry22 inner i.val j.val).hi ≤
      (matrixFixEntry22 outer i.val j.val).hi := of_decide_eq_true hij
  have hs : (0 : ℝ) ≤ (Spin.Numeric.scale : ℝ) := scaleR_pos.le
  exact ⟨(div_le_div_of_nonneg_right (by exact_mod_cast hc.1) hs).trans (hm i j).1,
    (hm i j).2.trans (div_le_div_of_nonneg_right (by exact_mod_cast hc.2) hs)⟩

theorem matrixFixSquareStageCheck22_sound {previous proposed : List (List Fix)}
    {num den : Int} (h : matrixFixSquareStageCheck22 previous proposed num den = true)
    {M : Matrix (Fin 11) (Fin 11) ℝ} (hm : MatrixFixMem22 previous M) :
    MatrixFixMem22 proposed (((num : ℝ)/(den : ℝ)) • (M*M)) := by
  simp only [matrixFixSquareStageCheck22, Bool.and_eq_true, decide_eq_true_eq] at h
  exact matrixFixEnclosesCheck22_sound h.2 (matrixFixSquareStep22_mem hm num den h.1.2)

#print axioms matrixFixSquareStageCheck22_sound
end Spin.Turbo

namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

def matrixScaledRealSteps22 (A : Matrix (Fin 11) (Fin 11) ℝ)
    (scales : Nat → Int × Int) : Nat → Matrix (Fin 11) (Fin 11) ℝ
  | 0 => A
  | k+1 => let X := matrixScaledRealSteps22 A scales k
    (((scales k).1 : ℝ)/((scales k).2 : ℝ)) • (X*X)

theorem matrixFixStages22_mem (matrices : Nat → List (List Fix))
    (scales : Nat → Int × Int) (A : Matrix (Fin 11) (Fin 11) ℝ)
    (h0 : MatrixFixMem22 (matrices 0) A) (n : Nat)
    (hstage : ∀ k < n, matrixFixSquareStageCheck22
      (matrices k) (matrices (k+1)) (scales k).1 (scales k).2 = true) :
    MatrixFixMem22 (matrices n) (matrixScaledRealSteps22 A scales n) := by
  induction n with
  | zero => exact h0
  | succ n ih =>
    exact matrixFixSquareStageCheck22_sound (hstage n (by omega))
      (ih (fun k hk => hstage k (by omega)))

def matrixScaledFactorSteps22 (scales : Nat → Int × Int) : Nat → ℝ
  | 0 => 1
  | k+1 => (((scales k).1 : ℝ)/((scales k).2 : ℝ)) *
      (matrixScaledFactorSteps22 scales k)^2

theorem matrixScaledRealSteps22_eq (A : Matrix (Fin 11) (Fin 11) ℝ)
    (scales : Nat → Int × Int) (n : Nat) :
    matrixScaledRealSteps22 A scales n = matrixScaledFactorSteps22 scales n • A^(2^n) := by
  induction n with
  | zero => simp [matrixScaledRealSteps22, matrixScaledFactorSteps22]
  | succ n ih =>
    rw [matrixScaledRealSteps22, ih, smul_mul_smul_comm, ← pow_two,
      ← pow_two, ← pow_mul, smul_smul]
    simp only [matrixScaledFactorSteps22, pow_succ]

#print axioms matrixFixStages22_mem
#print axioms matrixScaledRealSteps22_eq
end Spin.Turbo
