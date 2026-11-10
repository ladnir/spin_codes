import SpinCodes.Turbo.PhysicalCharacterPolynomial22
import SpinCodes.Turbo.Hi.PhysicalProfileNumericTable22
import SpinCodes.Turbo.PhysicalCauchy22

/-!
# Profile formulas for the actual physical moments

The local character argument and exact occupation normalization identify
the physical moments with the polynomials authenticated by the fixed-point
checker. The profile census is applied only after this actual-map identity.
-/
noncomputable section
namespace Spin.Turbo.Hi
open Finset Polynomial Spin.UniformMixer Spin.NumericHi Spin.NumericHi.Fix

theorem expandedPacket22_weight_eq_slot (q : State22) (slot : Fin 16) :
    packetWeight22 (expandedPacket22 q slot) = physicalSlotWeight22 q slot := rfl

attribute [local irreducible] packetWeight22 expandedPacket22 physicalSlotWeight22
  physicalProfile22 physicalProfileProduct22

private theorem profileLinearFactor22 (a b : ℝ) :
    Polynomial.monomial 0 a + Polynomial.monomial 1 b = C a + X*C b := by
  rw [monomial_zero_left, X_mul_C, C_mul_X_eq_monomial]

theorem physicalOccupationPolynomial22_eq_bank (q : State22) (z : ℝ) :
    physicalOccupationPolynomial22 q z =
      profileBankPolynomial22 (List.ofFn (physicalProfile22 q)) (profileRealMomentBank22 z) := by
  rw [physicalOccupationPolynomial22_product, profileBankPolynomial22_moment]
  calc
    _ = ∏ slot : Fin 16, (C (z^physicalSlotWeight22 q slot) +
        X*C ((1+z)^4-z^physicalSlotWeight22 q slot)) := by
      apply Finset.prod_congr rfl
      intro slot _
      simp only [expandedPacket22_weight_eq_slot]
      exact profileLinearFactor22 _ _
    _ = _ := physicalProfile22_product_nat q (fun r => C (z^r) + X*C ((1+z)^4-z^r))

theorem physicalCharacterOccupationPolynomial22_eq_bank (q : State22) (z : ℝ) :
    physicalCharacterOccupationPolynomial22 q z =
      profileBankPolynomial22 (List.ofFn (physicalProfile22 q)) (profileRealFourierBank22 z) := by
  rw [physicalCharacterOccupationPolynomial22_product, profileBankPolynomial22_fourier]
  calc
    _ = ∏ slot : Fin 16, (C (1 : ℝ) +
        X*C ((1+z)^(4-physicalSlotWeight22 q slot)*(1-z)^physicalSlotWeight22 q slot-1)) := by
      apply Finset.prod_congr rfl
      intro slot _
      simp only [expandedPacket22_weight_eq_slot]
      exact profileLinearFactor22 _ _
    _ = _ := physicalProfile22_product_nat q
      (fun r => C (1 : ℝ) + X*C ((1+z)^(4-r)*(1-z)^r-1))

theorem physicalInputMoment22_eq_profile (j : Fin 17) (q : State22) (z : ℝ) :
    physicalInputMoment22 (packetStepLaw22 j) z q =
      profileNumericMoment22 (List.ofFn (physicalProfile22 q)) z j := by
  rw [physicalInputMoment22_eq_coefficient, physicalOccupationPolynomial22_eq_bank]
  rfl

theorem physicalCharacterMoment22_eq_profile (j : Fin 17) (q : State22) (z : ℝ) :
    (packetStepLaw22 j).expect (fun x => z ^ stepWeight (packetStepInput22 x) 0 *
      (binarySign (statePairing q (stepFeedback (packetStepInput22 x))) : ℝ)) =
      profileNumericFourier22 (List.ofFn (physicalProfile22 q)) z j := by
  rw [physicalCharacterMoment22_eq_coefficient, physicalCharacterOccupationPolynomial22_eq_bank]
  rfl

theorem physicalInputMoment22_fixed_mem {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) (q : State22) :
    Mem (profileFixMoment22 (List.ofFn (physicalProfile22 q)) box j.val)
      (physicalInputMoment22 (packetStepLaw22 j) z q) := by
  rw [physicalInputMoment22_eq_profile]
  exact profileFixMoment22_sound _ hz j

theorem physicalInputSecondMoment22_fixed_mem {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) (q : State22) :
    Mem (profileFixSecondMoment22 (List.ofFn (physicalProfile22 q)) box j.val)
      (physicalInputSecondMoment22 (packetStepLaw22 j) z q) := by
  rw [physicalInputSecondMoment22_eq, physicalInputMoment22_eq_profile]
  exact profileFixSecondMoment22_sound _ hz j

theorem physicalZeroHold22_eq_profile_sum (j : Fin 17) (z : ℝ) :
    physicalZeroHold22 (packetStepLaw22 j) z =
      (∑ q : State22, (profileBankPolynomial22 (List.ofFn (physicalProfile22 q))
        (profileRealFourierBank22 z)).coeff j.val) /
          ((4194304 : ℝ) * ((Nat.choose 16 j.val : ℝ)*15^j.val)) := by
  have h := physicalZeroHold22_fourier_sum (packetStepLaw22 j) z
  conv at h =>
    rhs
    simp only [physicalCharacterMoment22_eq_profile, profileNumericFourier22]
    rw [← Finset.sum_div]
  have hp : (2^22 : ℝ) = 4194304 := by norm_num
  rw [hp] at h
  calc
    _ = ((∑ q : State22, (profileBankPolynomial22 (List.ofFn (physicalProfile22 q))
        (profileRealFourierBank22 z)).coeff j.val) /
          ((Nat.choose 16 j.val : ℝ)*15^j.val)) / 4194304 := by
      apply (eq_div_iff (by norm_num : (4194304 : ℝ) ≠ 0)).mpr
      simpa only [mul_comm] using h
    _ = _ := by rw [div_div, mul_comm _ (4194304 : ℝ)]

theorem physicalZeroHold22_fixed_mem {box : Fix} {z : ℝ} (hz : Mem box z)
    (j : Fin 17) :
    Mem (profileFixAllFourier22 profileNumericEntries22 box j.val)
      (physicalZeroHold22 (packetStepLaw22 j) z) := by
  rw [physicalZeroHold22_eq_profile_sum]
  exact profileNumericAllFourier22_actual_mem hz j

#print axioms physicalInputMoment22_eq_profile
#print axioms physicalZeroHold22_fixed_mem
end Spin.Turbo.Hi
