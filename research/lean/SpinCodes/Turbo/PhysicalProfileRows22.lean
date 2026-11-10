import SpinCodes.Turbo.PhysicalProfileMoment22
import SpinCodes.Turbo.PhysicalProfileWeighted22
import SpinCodes.Turbo.PhysicalExactRow22
import SpinCodes.Turbo.PhysicalCoarseBirth22

/-!
# Actual physical rows from authenticated profile observables

This combines the old maximum/mean distinction with the exact packet
moment identity. The birth cone and numerical upper bounds remain explicit
premises; this file does not itself certify any finite endpoint.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem physicalProfile22_first_weighted (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (i : Fin 11) (hi : i ≠ 0) (M : ℝ) (hM : 0 ≤ M)
    (hcap : i ≠ 2 → ∀ q, physicalClassAllows22 i q →
      profileNumericMoment22 (List.ofFn (physicalProfile22 q)) z j ≤ M)
    (hmean : i = 2 → physicalProfileClassMean22 2
      (fun p => profileNumericMoment22 (List.ofFn p) z j) ≤ M)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ) :
    (∑ q, μ q * physicalInputMoment22 (packetStepLaw22 j) z q) ≤ b*M := by
  have hn (q : State22) :
      0 ≤ profileNumericMoment22 (List.ofFn (physicalProfile22 q)) z j := by
    rw [← physicalInputMoment22_eq_profile]
    exact FinPMF.expect_nonneg _ (fun _ => pow_nonneg hz _)
  have h := physicalProfile22_weighted_bound i hi b μ hμ
    (fun p => profileNumericMoment22 (List.ofFn p) z j) hn M hM hcap hmean
  simpa only [← physicalInputMoment22_eq_profile] using h

theorem physicalProfile22_second_weighted (j : Fin 17) (z : ℝ)
    (i : Fin 11) (hi : i ≠ 0) (S : ℝ) (hS : 0 ≤ S)
    (hcap : i ≠ 2 → ∀ q, physicalClassAllows22 i q →
      profileNumericMoment22 (List.ofFn (physicalProfile22 q)) (z*z) j ≤ S)
    (hmean : i = 2 → physicalProfileClassMean22 2
      (fun p => profileNumericMoment22 (List.ofFn p) (z*z) j) ≤ S)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ) :
    (∑ q, μ q * physicalInputSecondMoment22 (packetStepLaw22 j) z q) ≤ b*S := by
  simpa only [physicalInputSecondMoment22_eq] using
    physicalProfile22_first_weighted j (z*z) (mul_self_nonneg z) i hi S hS hcap hmean b μ hμ

private theorem turboPhysicalKernel22_eq_input (j : Fin 17) (z : ℝ) :
    Matrix.of (turboPhysicalKernel22 j z) = physicalInputKernel22 (packetStepLaw22 j) z := rfl

private theorem turboExactBirthMass22_eq_input_row (j : Fin 17) (z : ℝ) :
    turboExactBirthMass22 j z = physicalInputKernel22 (packetStepLaw22 j) z 0 := by
  calc
    _ = (Matrix.of (turboPhysicalKernel22 j z)) 0 := by
      funext r
      exact (turboPhysicalKernel22_zero_entry j z r).symm
    _ = _ := congrArg (fun K : Matrix State22 State22 ℝ => K 0)
      (turboPhysicalKernel22_eq_input j z)

theorem turboExactBirthMass22_coarse_cone (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (Z T : ℝ) (hZ : physicalZeroHold22 (packetStepLaw22 j) z ≤ Z)
    (hT : physicalInputMoment22 (packetStepLaw22 j) z 0 ≤ T) :
    physicalCone22 (physicalCoarseBirthRow22 Z T) (turboExactBirthMass22 j z) := by
  have hc : physicalCone22 (physicalCoarseBirthRow22 Z T)
      (physicalInputKernel22 (packetStepLaw22 j) z 0) := physicalCoarseBirth22_parts _ Z T
    (physicalInputKernel22_nonneg _ z hz 0)
    (by rw [physicalInputKernel22_zero_zero]; exact hZ)
    (by rw [physicalInputKernel22_row_sum]; exact hT)
  exact (congrArg (physicalCone22 (physicalCoarseBirthRow22 Z T))
    (turboExactBirthMass22_eq_input_row j z)).mpr hc

private theorem probability_decidable_eq22 {Ω : Type*} [Fintype Ω]
    (P : FinPMF Ω) (E : Ω → Prop) (a b : DecidablePred E) :
    @FinPMF.prob Ω _ P E a = @FinPMF.prob Ω _ P E b :=
  congrArg (fun d => @FinPMF.prob Ω _ P E d) (Subsingleton.elim a b)

attribute [local irreducible] FinPMF.prob

theorem physicalFeedbackProbability22_zero_hold (P : FinPMF PacketStep16) :
    physicalFeedbackProbability22 P = 1 - physicalZeroHold22 P 1 := by
  calc
    _ = P.prob (fun x => stepFeedback (packetStepInput22 x) ≠ 0) := by
      unfold physicalFeedbackProbability22
      exact probability_decidable_eq22 P _ _ _
    _ = _ := physicalFeedbackActive22_complement P

#print axioms physicalProfile22_first_weighted
#print axioms physicalProfile22_second_weighted
#print axioms turboExactBirthMass22_coarse_cone
#print axioms physicalFeedbackProbability22_zero_hold
end Spin.Turbo
