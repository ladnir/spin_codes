import SpinCodes.Turbo.PhysicalCone22

/-!
# Exact birth measure for the Turbo physical step

An entering zero state is not refreshed by the sampled invertible map.
Its next state is exactly the feedback of the complete four-bit input,
including selected slots whose value is zero. This is the source row
that the retained 11-coordinate producer represents by zero-state,
uniform-density, and expansion-weight-class masses.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

def turboBirthMass22 (j : Fin 17) (z : ℝ)
    (r : State22) : ℝ :=
  (selectedStepLaw22 j).expect (fun x =>
    if stepFeedback (packetStepInput22 x) = r then
      z ^ stepWeight (packetStepInput22 x) 0 else 0)

theorem turboSelectedKernel22_zero_entry (j : Fin 17)
    (z : ℝ) (r : State22) :
    turboSelectedKernel22 j z 0 r = turboBirthMass22 j z r := by
  unfold turboSelectedKernel22 turboBirthMass22
  apply congrArg
  funext x
  simp [stepKernel, update_zero, eq_comm]

theorem turboBirthMass22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (r : State22) :
    0 ≤ turboBirthMass22 j z r := by
  unfold turboBirthMass22
  apply FinPMF.expect_nonneg
  intro x
  by_cases h : stepFeedback (packetStepInput22 x) = r
  · simp only [h, ite_true]
    exact pow_nonneg hz _
  · simp only [h, ite_false]
    exact le_refl 0

theorem turboBirthMass22_total (j : Fin 17) (z : ℝ) :
    (∑ r : State22, turboBirthMass22 j z r) =
      (selectedStepLaw22 j).expect (fun x =>
        z ^ stepWeight (packetStepInput22 x) 0) := by
  simp only [turboBirthMass22]
  rw [← FinPMF.expect_sum]
  apply congrArg
  funext x
  simp [Finset.mem_univ]

#print axioms turboSelectedKernel22_zero_entry
#print axioms turboBirthMass22_nonneg
#print axioms turboBirthMass22_total

end Spin.Turbo
