import SpinCodes.UniformMixer.ScalableDenseVertex
import SpinCodes.UniformMixer.ScalableDenseRateDefs
import SpinCodes.UniformMixer.ScalableSpectrumPolynomial

/-! A representative rational dense witness for testing the Lean certificate path. -/

noncomputable section
set_option maxHeartbeats 1000000
set_option maxRecDepth 8192

namespace Spin.UniformMixer

private def pilotVector : Fin 3 → ℝ := ![1, 1 / 10000, 1 / 1000000]

private theorem pilot_H_upper :
    spectrumH (2 / 5) (1 / 5) ≤ (7 / 1000000000000000 : ℝ) := by
  rw [spectrumH_polynomial]
  norm_num [productionSpectrumPolynomial, feedbackWeightMoment]

private theorem pilot_H_nonneg :
    0 ≤ spectrumH (2 / 5) (1 / 5) := by
  rw [spectrumH_polynomial]
  norm_num [productionSpectrumPolynomial, feedbackWeightMoment]

private theorem pilot_G_upper :
    spectrumG (2 / 5) (1 / 5) ≤ (2 / 100000000000 : ℝ) := by
  have hpow : ((17 / 25 : ℝ) ^ 64) ≤ 2 / 100000000000 := by norm_num
  unfold spectrumG
  norm_num only at *
  linarith [pilot_H_nonneg]

private theorem pilot_F_upper :
    spectrumF (2 / 5) (1 / 5) ≤ (3 / 10000000000000 : ℝ) := by
  unfold spectrumF scalableSpectrumMaxMoment
  apply Finset.sup'_le
  intro w hw
  simp only [liveSpectrumWeights, spectrumWeights, Finset.mem_erase,
    Finset.mem_insert, Finset.mem_singleton] at hw
  rcases hw.2 with h | h | h | h | h | h | h | h | h
  · exact (hw.1 h).elim
  all_goals subst w
  all_goals norm_num [productionWeightMoment]

private theorem pilot_U_upper :
    spectrumU (2 / 5) (1 / 5) ≤ (7 / 1000000000000000 : ℝ) := by
  rw [spectrumU_polynomial]
  rw [state22_liveCount]
  norm_num [productionSpectrumPolynomial, productionWeightMoment]

private theorem pilot_transfer :
    ∀ i, (paperSpectrumTransfer (2 / 5) (1 / 5)).mulVec
      pilotVector i ≤ (1 / 100000000000000 : ℝ) * pilotVector i := by
  intro i
  fin_cases i
  all_goals dsimp [paperSpectrumTransfer, spectrumTransfer,
    weightedTransfer, Matrix.mulVec, dotProduct,
    Fin.sum_univ_three, pilotVector]
  all_goals norm_num only [Fin.sum_univ_three]
  all_goals simp only [show (0 : Fin 3) ≠ 1 by decide,
    show (1 : Fin 3) ≠ 0 by decide,
    show (2 : Fin 3) ≠ 0 by decide,
    show (2 : Fin 3) ≠ 1 by decide,
    ite_true, ite_false, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, mul_one, zero_mul, add_zero]
  · change spectrumH (2 / 5) (1 / 5) +
      spectrumG (2 / 5) (1 / 5) * (1 / 10000) ≤
      (1 / 100000000000000 : ℝ)
    linarith [pilot_H_upper, pilot_G_upper]
  · change min (spectrumF (2 / 5) (1 / 5))
        (spectrumu (2 / 5)) / (liveCount State22 : ℝ) +
      spectrumF (2 / 5) (1 / 5) * (1 / 1000000) ≤
        (1 / 1000000000000000000 : ℝ)
    rw [state22_liveCount]
    have hv : min (spectrumF (2 / 5) (1 / 5))
        (spectrumu (2 / 5)) ≤ spectrumF (2 / 5) (1 / 5) := min_le_left _ _
    nlinarith [pilot_F_upper]
  · change min (spectrumU (2 / 5) (1 / 5))
        (spectrumu (2 / 5)) / (liveCount State22 : ℝ) +
      spectrumU (2 / 5) (1 / 5) * (1 / 1000000) ≤
        (1 / 100000000000000000000 : ℝ)
    rw [state22_liveCount]
    have hv : min (spectrumU (2 / 5) (1 / 5))
        (spectrumu (2 / 5)) ≤ spectrumU (2 / 5) (1 / 5) := min_le_left _ _
    nlinarith [pilot_U_upper]

#print axioms pilot_H_upper
#print axioms pilot_F_upper
#print axioms pilot_U_upper
#print axioms pilot_transfer

private def pilotSlope : Spin.Structured.DenseOccupationFixed.QInput := ⟨0, 1⟩
private def pilotIntercept : Spin.Structured.DenseOccupationFixed.QInput :=
  ⟨3465735902799727, 10000000000000000⟩
private def pilotP : Spin.Structured.DenseOccupationFixed.QInput := ⟨3, 5⟩
private def pilotY : Spin.Structured.DenseOccupationFixed.QInput := ⟨2, 3⟩
private def pilotRadius : Spin.Structured.DenseOccupationFixed.QInput :=
  ⟨1, 100000000000000⟩
private def pilotZ : Spin.Structured.DenseOccupationFixed.QInput := ⟨1, 5⟩
private def pilotA0 : Spin.Structured.DenseOccupationFixed.QInput := ⟨49, 100⟩
private def pilotA1 : Spin.Structured.DenseOccupationFixed.QInput := ⟨51, 100⟩
private def pilotX0 : Spin.Structured.DenseOccupationFixed.QInput := ⟨499, 1000⟩
private def pilotX1 : Spin.Structured.DenseOccupationFixed.QInput := ⟨501, 1000⟩

private theorem pilot_check00 :
    scalableVertexCheck pilotSlope pilotIntercept pilotP pilotY
      pilotRadius pilotZ pilotA0 pilotX0 50
      (-1000000000000000000000000000) = true := by decide +kernel
private theorem pilot_check01 :
    scalableVertexCheck pilotSlope pilotIntercept pilotP pilotY
      pilotRadius pilotZ pilotA0 pilotX1 50
      (-1000000000000000000000000000) = true := by decide +kernel
private theorem pilot_check10 :
    scalableVertexCheck pilotSlope pilotIntercept pilotP pilotY
      pilotRadius pilotZ pilotA1 pilotX0 50
      (-1000000000000000000000000000) = true := by decide +kernel
private theorem pilot_check11 :
    scalableVertexCheck pilotSlope pilotIntercept pilotP pilotY
      pilotRadius pilotZ pilotA1 pilotX1 50
      (-1000000000000000000000000000) = true := by decide +kernel

#print axioms pilot_check00
#print axioms pilot_check01
#print axioms pilot_check10
#print axioms pilot_check11

private theorem pilot_corner (a x : Spin.Structured.DenseOccupationFixed.QInput)
    (ha : 0 < a.den) (hx : 0 < x.den)
    (h : scalableVertexCheck pilotSlope pilotIntercept pilotP pilotY
      pilotRadius pilotZ a x 50
      (-1000000000000000000000000000) = true) :
    scalableBoxExponent pilotSlope.real pilotIntercept.real pilotP.real
      pilotY.real pilotRadius.real pilotZ.real a.real x.real ≤
      -(1 / 1000 : ℝ) := by
  have hh := scalableVertexCheck_sound pilotSlope pilotIntercept
    pilotP pilotY pilotRadius pilotZ a x
    (by decide) (by decide) ha hx
    (by decide) (by decide) (by decide) (by decide) h
  norm_num [Spin.Numeric.scale] at hh
  exact hh

private theorem pilot_box {α x : ℝ}
    (hα : α ∈ Set.Icc pilotA0.real pilotA1.real)
    (hx : x ∈ Set.Icc pilotX0.real pilotX1.real) :
    scalableBoxExponent pilotSlope.real pilotIntercept.real pilotP.real
      pilotY.real pilotRadius.real pilotZ.real α x ≤
      -(1 / 1000 : ℝ) := by
  have h := scalableBoxExponent_le_vertices
    pilotSlope.real pilotIntercept.real pilotP.real pilotY.real
    pilotRadius.real pilotZ.real
    (by norm_num [pilotP, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotP, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotY, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotY, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotA0, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotA1, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotX0, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotX1, Spin.Structured.DenseOccupationFixed.QInput.real])
    hα hx
  apply h.trans
  apply max_le (max_le ?_ ?_) (max_le ?_ ?_)
  · exact pilot_corner pilotA0 pilotX0 (by decide) (by decide) pilot_check00
  · exact pilot_corner pilotA0 pilotX1 (by decide) (by decide) pilot_check01
  · exact pilot_corner pilotA1 pilotX0 (by decide) (by decide) pilot_check10
  · exact pilot_corner pilotA1 pilotX1 (by decide) (by decide) pilot_check11

private theorem pilot_pointRate {α x : ℝ}
    (hα : α ∈ Set.Icc pilotA0.real pilotA1.real)
    (hx : x ∈ Set.Icc pilotX0.real pilotX1.real) :
    ScalableDensePointRate α x pilotSlope.real pilotIntercept.real
      (1 / 1000) (1000000 : ℝ) := by
  have hα0 : 0 < α := lt_of_lt_of_le
    (by norm_num [pilotA0, Spin.Structured.DenseOccupationFixed.QInput.real]) hα.1
  have hα1 : α ≤ 1 := le_trans hα.2
    (by norm_num [pilotA1, Spin.Structured.DenseOccupationFixed.QInput.real])
  have hx0 : 0 < x := lt_of_lt_of_le
    (by norm_num [pilotX0, Spin.Structured.DenseOccupationFixed.QInput.real]) hx.1
  have hx1 : x ≤ 1 := le_trans hx.2
    (by norm_num [pilotX1, Spin.Structured.DenseOccupationFixed.QInput.real])
  have hones : ∀ i, (1 : ℝ) ≤ 1000000 * pilotVector i := by
    intro i
    fin_cases i <;> norm_num [pilotVector]
  have hv : ∀ i,
      (paperSpectrumTransfer (pilotP.real * pilotY.real)
        pilotZ.real).mulVec pilotVector i ≤
        pilotRadius.real * pilotVector i := by
    intro i
    convert pilot_transfer i using 1 <;>
      norm_num [pilotP, pilotY, pilotZ, pilotRadius,
        Spin.Structured.DenseOccupationFixed.QInput.real]
  have h := scalable_dense_witness_pointRate
    hα0 hα1 hx0 hx1
    (by norm_num [pilotP, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotP, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotY, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotY, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotRadius, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotZ, Spin.Structured.DenseOccupationFixed.QInput.real])
    (by norm_num [pilotZ, Spin.Structured.DenseOccupationFixed.QInput.real])
    pilotVector (by norm_num : (0 : ℝ) ≤ 1000000)
    hones hv (pilot_box hα hx)
  simpa only [pilotVector, Matrix.cons_val_zero, mul_one] using h

#print axioms pilot_corner
#print axioms pilot_box
#print axioms pilot_pointRate

/-- One actual dense rectangle reaches the selected-outer rate interface. -/
theorem scalable_dense_pilot_coverage {α x : ℝ}
    (hα : α ∈ Set.Icc (49 / 100) (51 / 100))
    (hx : x ∈ Set.Icc (499 / 1000) (501 / 1000)) :
    ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧
      ScalableDensePointRate α x p.1 p.2 (1 / 1000) 1000000 := by
  let p : ℚ × ℚ := (0, 3465735902799727 / 10000000000000000)
  have hp : p ∈ Spin.Majorant.refined.supports := by
    change p ∈ Spin.Majorant.supports
    simpa [p] using Spin.Majorant.member19
  refine ⟨p, hp, ?_⟩
  have hα' : α ∈ Set.Icc pilotA0.real pilotA1.real := by
    simpa [pilotA0, pilotA1,
      Spin.Structured.DenseOccupationFixed.QInput.real] using hα
  have hx' : x ∈ Set.Icc pilotX0.real pilotX1.real := by
    simpa [pilotX0, pilotX1,
      Spin.Structured.DenseOccupationFixed.QInput.real] using hx
  convert (pilot_pointRate hα' hx') using 1 <;>
    norm_num [p, pilotSlope, pilotIntercept,
      Spin.Structured.DenseOccupationFixed.QInput.real]

#print axioms scalable_dense_pilot_coverage

end Spin.UniformMixer
