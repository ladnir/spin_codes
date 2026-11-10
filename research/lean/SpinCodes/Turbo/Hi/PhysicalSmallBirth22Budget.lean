import SpinCodes.Turbo.Hi.PhysicalSmallBirth22BudgetDefs
import SpinCodes.Turbo.PhysicalSmallBirth22Polynomial
import SpinCodes.Turbo.Hi.PhysicalProfileNumeric22

/-!
# Checked integer birth budgets for the actual small physical layers

Audit: Numeric.Fixed proves outward add/mul/divInt enclosures against reals;
FixedDefs keeps computation on Int and uses Horner to avoid repeated powers.
Here only the already authenticated class-count polynomial is evaluated.
Endpoint comparisons reuse profileFixUpperCheck22 and its soundness theorem.
-/
noncomputable section
namespace Spin.Turbo.Hi
open Finset Spin.NumericHi Spin.NumericHi.Fix

private theorem birthBudgetZeroMem22 : Mem Fix.zero 0 := by simp [Mem,Fix.zero]

private theorem smallBirthOneBudgetCoefficients22 : ∀ i : Fin 11,
    smallBirthFixClassCoefficients22 smallBirthOneHistogram22 i.val =
      List.ofFn (fun t : Fin 9 => smallBirthOneHistogram22.getD (9*i.val+t.val) 0) := by
  decide +kernel

private theorem smallBirthTwoBudgetCoefficients22 : ∀ i : Fin 11,
    smallBirthFixClassCoefficients22 smallBirthTwoHistogram22 i.val =
      List.ofFn (fun t : Fin 9 => smallBirthTwoHistogram22.getD (9*i.val+t.val) 0) := by
  decide +kernel

theorem smallBirthFixHorner22_mem {n : Nat} (cs : Fin n → Nat)
    {box : Fix} {z : ℝ} (hz : Mem box z) :
    Mem (smallBirthFixHorner22 (List.ofFn cs) box)
      (∑ t : Fin n, (cs t : ℝ) * z ^ t.val) := by
  induction n with
  | zero => simpa [smallBirthFixHorner22] using birthBudgetZeroMem22
  | succ n ih =>
    have hs : (∑ t : Fin (n+1), (cs t : ℝ) * z ^ t.val) =
        (cs 0 : ℝ) + z * ∑ t : Fin n, (cs t.succ : ℝ) * z ^ t.val := by
      rw [Fin.sum_univ_succ]
      simp only [Fin.val_zero, pow_zero, mul_one, Fin.val_succ, pow_succ]
      rw [Finset.mul_sum]
      congr 1
      apply Finset.sum_congr rfl
      intro t _
      ring
    rw [hs, List.ofFn_succ, smallBirthFixHorner22]
    exact add_mem (by simpa only [Int.cast_natCast] using ofInt_mem (cs 0 : Int))
      (mul_mem hz (ih (fun t => cs t.succ)))

theorem smallBirthFixOneBudget22_mem {box : Fix} {z : ℝ}
    (hz : Mem box z) (i : Fin 11) :
    Mem (smallBirthFixOneBudget22 box i.val) (physicalSmallBirthOneBudget22 z i) := by
  unfold smallBirthFixOneBudget22 smallBirthFixBudget22
  rw [smallBirthOneBudgetCoefficients22 i]
  have hh := divInt_mem (smallBirthFixHorner22_mem
    (fun t : Fin 9 => smallBirthOneHistogram22.getD (9*i.val+t.val) 0) hz)
    (show (0 : Int) < 240 by decide)
  simpa only [physicalSmallBirthOneBudget22, Int.cast_ofNat, div_eq_mul_inv,
    one_div, one_mul, mul_comm] using hh

theorem smallBirthFixTwoBudget22_mem {box : Fix} {z : ℝ}
    (hz : Mem box z) (i : Fin 11) :
    Mem (smallBirthFixTwoBudget22 box i.val) (physicalSmallBirthTwoBudget22 z i) := by
  unfold smallBirthFixTwoBudget22 smallBirthFixBudget22
  rw [smallBirthTwoBudgetCoefficients22 i]
  have hh := divInt_mem (smallBirthFixHorner22_mem
    (fun t : Fin 9 => smallBirthTwoHistogram22.getD (9*i.val+t.val) 0) hz)
    (show (0 : Int) < 27000 by decide)
  simpa only [physicalSmallBirthTwoBudget22, Int.cast_ofNat, div_eq_mul_inv,
    one_div, one_mul, mul_comm] using hh

theorem smallBirthFixZeroBudget22_mem (i : Fin 11) :
    Mem (smallBirthFixZeroBudget22 i.val) (physicalInitialEnvelope22 i) := by
  by_cases hi : i.val = 0
  · simp only [smallBirthFixZeroBudget22, physicalInitialEnvelope22, hi, ite_true]
    simpa only [Int.cast_one] using ofInt_mem 1
  · simp only [smallBirthFixZeroBudget22, physicalInitialEnvelope22, hi, ite_false]
    exact birthBudgetZeroMem22

theorem smallBirthFixVectorUpperCheck22_sound (values : Nat → Fix) (realValues : Fin 11 → ℝ)
    (hm : ∀ i, Mem (values i.val) (realValues i)) {nums : List Int} {den : Int}
    (h : smallBirthFixVectorUpperCheck22 values nums den = true) :
    ∀ i : Fin 11, realValues i ≤ (nums.getD i.val 0 : ℝ) / (den : ℝ) := by
  simp only [smallBirthFixVectorUpperCheck22, Bool.and_eq_true] at h
  have hp := h
  have hd : 0 < den := of_decide_eq_true hp.1
  intro i
  have hc := List.all_eq_true.mp hp.2 i.val (List.mem_range.mpr i.isLt)
  exact profileFixUpperCheck22_sound (hm i) hd hc

private theorem birthBudgetNonnegOfBox22 {box : Fix} {z : ℝ}
    (hz : Mem box z) (hlo : 0 ≤ box.lo) : 0 ≤ z := by
  have hl : (0 : ℝ) ≤ (box.lo : ℝ) := by exact_mod_cast hlo
  exact (div_nonneg hl scaleR_pos.le).trans hz.1

theorem smallBirthFixOneBudgetCheck22_cone {box : Fix} {z : ℝ}
    (hz : Mem box z) {nums : List Int} {den : Int}
    (h : smallBirthFixOneBudgetCheck22 box nums den = true) :
    physicalCone22 (fun i => (nums.getD i.val 0 : ℝ) / (den : ℝ))
      (turboExactBirthMass22 1 z) := by
  simp only [smallBirthFixOneBudgetCheck22, Bool.and_eq_true] at h
  have hp := h
  have hlo : 0 ≤ box.lo := of_decide_eq_true hp.1
  apply physicalSmallBirth22_one_cone_of_le z (birthBudgetNonnegOfBox22 hz hlo)
  exact smallBirthFixVectorUpperCheck22_sound _ _ (smallBirthFixOneBudget22_mem hz) hp.2

theorem smallBirthFixTwoBudgetCheck22_cone {box : Fix} {z : ℝ}
    (hz : Mem box z) {nums : List Int} {den : Int}
    (h : smallBirthFixTwoBudgetCheck22 box nums den = true) :
    physicalCone22 (fun i => (nums.getD i.val 0 : ℝ) / (den : ℝ))
      (turboExactBirthMass22 2 z) := by
  simp only [smallBirthFixTwoBudgetCheck22, Bool.and_eq_true] at h
  have hp := h
  have hlo : 0 ≤ box.lo := of_decide_eq_true hp.1
  apply physicalSmallBirth22_two_cone_of_le z (birthBudgetNonnegOfBox22 hz hlo)
  exact smallBirthFixVectorUpperCheck22_sound _ _ (smallBirthFixTwoBudget22_mem hz) hp.2

theorem smallBirthFixZeroBudgetCheck22_cone (z : ℝ) (hz : 0 ≤ z)
    {nums : List Int} {den : Int} (h : smallBirthFixZeroBudgetCheck22 nums den = true) :
    physicalCone22 (fun i => (nums.getD i.val 0 : ℝ) / (den : ℝ))
      (turboExactBirthMass22 0 z) := by
  apply physicalSmallBirth22_cone_of_class_bounds 0 z hz
  intro i
  rw [physicalSmallBirthClassMass22_empty]
  exact smallBirthFixVectorUpperCheck22_sound _ _ smallBirthFixZeroBudget22_mem h i

#print axioms smallBirthFixOneBudgetCheck22_cone
#print axioms smallBirthFixTwoBudgetCheck22_cone
#print axioms smallBirthFixZeroBudgetCheck22_cone
end Spin.Turbo.Hi

