import SpinCodes.Turbo.PhysicalCauchy22
import SpinCodes.Turbo.PhysicalCappedBirth22
import SpinCodes.Turbo.PhysicalThinningRows22

/-!
# Complete exact-nonzero physical rows

This follows `PhysicalSymbolicRow22`'s complete-row assembly, replacing its
selected-slot class-only birth measure by the actual exact-nonzero law and
one jointly capped birth decomposition. Numerical census bounds remain
explicit hypotheses; no numerical endpoint is accepted as an axiom.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer
attribute [local instance] Classical.propDecidable
local instance exactRowState22DecidableEq : DecidableEq State22 :=
  fun a b => Classical.propDecidable (a = b)

private theorem cappedBirthClassPart22_zero_index (μ : State22 → ℝ) (d : ℝ)
    (r : State22) : cappedBirthClassPart22 μ d 0 r = 0 := by
  unfold cappedBirthClassPart22
  split_ifs with hindex
  · have ha := physicalBirthIndex22_allows r
    rw [hindex] at ha
    have hr : r = 0 := by simpa only [physicalClassAllows22, Fin.val_zero, ite_true] using ha
    subst r
    simp [cappedBirthResidual22]
  · rfl

private theorem cappedBirthClassPart22_uniform_index (μ : State22 → ℝ) (d : ℝ)
    (r : State22) : cappedBirthClassPart22 μ d 2 r = 0 := by
  simp only [cappedBirthClassPart22, physicalBirthIndex22_ne_two, ite_false]

def cappedBirthBudget22 (μ U : State22 → ℝ) (d : ℝ) (i : Fin 11) : ℝ :=
  if i = 0 then μ 0 else if i = 2 then (liveCount State22 : ℝ)*d
  else cappedBirthClassBound22 μ U d i

/-- One complete cone witness for the same cap in every output coordinate. -/
theorem cappedBirth22_cone (μ U : State22 → ℝ) (d : ℝ)
    (hμ : ∀ r, 0 ≤ μ r) (hU : ∀ r, μ r ≤ U r) (hd : 0 ≤ d) :
    physicalCone22 (cappedBirthBudget22 μ U d) μ := by
  let zeroPart : State22 → ℝ := fun r => μ 0 * physicalInitialMass22 r
  let parts : Fin 11 → State22 → ℝ := fun i r =>
    if i = 0 then zeroPart r else if i = 2 then cappedBirthUniform22 μ d r
      else cappedBirthClassPart22 μ d i r
  have hzero : physicalPiece22 0 (μ 0) zeroPart := by
    simpa only [zeroPart, mul_one] using
      physicalPiece22_scale physicalInitialMass22_piece (hμ 0)
  have hpart (i : Fin 11) (r : State22) : parts i r =
      (if i = 0 then zeroPart r else 0) +
        (if i = 2 then cappedBirthUniform22 μ d r else 0) + cappedBirthClassPart22 μ d i r := by
    by_cases hi0 : i = 0
    · subst i
      simp [parts, cappedBirthClassPart22_zero_index]
    · by_cases hi2 : i = 2
      · subst i
        simp [parts, cappedBirthClassPart22_uniform_index]
      · simp only [parts, hi0, hi2, ite_false, zero_add]
  refine ⟨parts, ?_, ?_⟩
  · intro r
    symm
    calc
      (∑ i, parts i r) = zeroPart r + cappedBirthUniform22 μ d r +
          ∑ i, cappedBirthClassPart22 μ d i r := by
        simp_rw [hpart]
        rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
        simp only [Finset.sum_ite_eq', Finset.mem_univ, ite_true]
      _ = μ r := by
        have hzeroEq : zeroPart r = if r = 0 then μ 0 else 0 := by
          dsimp [zeroPart, physicalInitialMass22]
          split_ifs <;> ring
        rw [hzeroEq]
        have hdec : Fintype.decidablePiFintype r (0 : State22) =
            exactRowState22DecidableEq r 0 := Subsingleton.elim _ _
        rw [← hdec]
        exact (cappedBirth22_decomposition μ d r).symm
  · intro i
    by_cases hi0 : i = 0
    · subst i
      simpa only [parts, cappedBirthBudget22, ite_true] using hzero
    · by_cases hi2 : i = 2
      · subst i
        simpa only [parts, cappedBirthBudget22, hi0, ite_false, ite_true] using
          cappedBirthUniform22_piece μ d hμ hd
      · simpa only [parts, cappedBirthBudget22, hi0, hi2, ite_false] using
          cappedBirthClassPart22_piece μ U d hμ hU hd i

/-- Any authenticated birth cone transfers an arbitrary zero-source piece. -/
theorem turboPhysicalKernel22_zero_row_of_cone (j : Fin 17) (z : ℝ)
    (v : Fin 11 → ℝ) (hv : physicalCone22 v (turboExactBirthMass22 j z))
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 0 b μ) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k (b*v k) (parts k) := by
  obtain ⟨birth, hbirth, hpieces⟩ := hv
  have hc : 0 ≤ μ 0 := hμ.1 0
  have hcb : μ 0 ≤ b := physicalZeroPiece22_mass hμ
  refine ⟨fun k r => μ 0 * birth k r, ?_, ?_⟩
  · intro r
    have hk : physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r =
        μ 0 * turboExactBirthMass22 j z r := by
      unfold physicalKernelAction22
      rw [Finset.sum_eq_single 0]
      · exact congrArg (μ 0 * ·) (turboPhysicalKernel22_zero_entry j z r)
      · intro q _ hq
        rw [physicalZeroPiece22_supported hμ q hq]
        exact zero_mul _
      · intro h; exact (h (Finset.mem_univ 0)).elim
    rw [hk, hbirth r, Finset.mul_sum]
  · intro k
    have hv0 : 0 ≤ v k :=
      (Finset.sum_nonneg (fun r _ => (hpieces k).1 r)).trans (hpieces k).2.2.1
    exact physicalPiece22_mono (physicalPiece22_scale (hpieces k) hc)
      (mul_le_mul_of_nonneg_right hcb hv0)

def physicalExactCauchyRows22 (v M S : Fin 11 → ℝ) (F : ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ := Matrix.of fun i k =>
  if i = 0 then v k else physicalLiveCapRow22 (M i) (min (M i) (Real.sqrt (S i*F))) k

/-- Complete exact-nonzero physical matrix, before binomial thinning.
All numerical hypotheses concern the same selected birth representation and
actual moment functions of `packetStepLaw22`. -/
theorem turboPhysicalKernel22_cauchy_rows_sound (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (v M S : Fin 11 → ℝ) (F : ℝ)
    (hbirth : physicalCone22 v (turboExactBirthMass22 j z))
    (hS : ∀ i, i ≠ 0 → 0 ≤ S i)
    (hfirst : ∀ i, i ≠ 0 → ∀ b μ, physicalPiece22 i b μ →
      (∑ q, μ q * physicalInputMoment22 (packetStepLaw22 j) z q) ≤ b*M i)
    (hsecond : ∀ i, i ≠ 0 → ∀ b μ, physicalPiece22 i b μ →
      (∑ q, μ q * physicalInputSecondMoment22 (packetStepLaw22 j) z q) ≤ b*S i)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤ F) :
    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z))
      (physicalExactCauchyRows22 v M S F) := by
  intro i b μ hμ
  by_cases hi : i = 0
  · subst i
    simpa only [physicalExactCauchyRows22, Matrix.of_apply, ite_true] using
      turboPhysicalKernel22_zero_row_of_cone j z v hbirth b μ hμ
  · simpa only [physicalExactCauchyRows22, Matrix.of_apply, hi, ite_false] using
      turboPhysicalKernel22_live_row_of_weighted_cauchy j z hz i hi b μ hμ
        (M i) (S i) F (hS i hi) (hfirst i hi b μ hμ) (hsecond i hi b μ hμ) hfeedback

/-- The retained dyadic matrix can be used only after complete-row domination. -/
theorem turboPhysicalKernel22_outward_cauchy_rows_sound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (v M S : Fin 11 → ℝ) (F : ℝ)
    (hbirth : physicalCone22 v (turboExactBirthMass22 j z))
    (hS : ∀ i, i ≠ 0 → 0 ≤ S i)
    (hfirst : ∀ i, i ≠ 0 → ∀ b μ, physicalPiece22 i b μ →
      (∑ q, μ q * physicalInputMoment22 (packetStepLaw22 j) z q) ≤ b*M i)
    (hsecond : ∀ i, i ≠ 0 → ∀ b μ, physicalPiece22 i b μ →
      (∑ q, μ q * physicalInputSecondMoment22 (packetStepLaw22 j) z q) ≤ b*S i)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤ F)
    (B : Matrix (Fin 11) (Fin 11) ℝ)
    (hB : ∀ i k, physicalExactCauchyRows22 v M S F i k ≤ B i k) :
    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z)) B :=
  physicalRowSound22_mono
    (turboPhysicalKernel22_cauchy_rows_sound j z hz v M S F hbirth hS hfirst hsecond hfeedback) hB

#print axioms cappedBirth22_cone
#print axioms turboPhysicalKernel22_zero_row_of_cone
#print axioms turboPhysicalKernel22_cauchy_rows_sound
#print axioms turboPhysicalKernel22_outward_cauchy_rows_sound
end Spin.Turbo

