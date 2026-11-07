import SpinCodes.UniformMixer.ScalableDenseLeaf
import SpinCodes.UniformMixer.ScalableDenseTransferExactDefs

/-! A second dense witness, found by local fallback rather than the initial bank. -/

noncomputable section
set_option maxRecDepth 8192
set_option maxHeartbeats 1000000

namespace Spin.UniformMixer

open Spin.Numeric Spin.Structured.DenseOccupationFixed

private def mixedSlope : QInput := ⟨-833, 500⟩
private def mixedIntercept : QInput := ⟨373189, 250000⟩
private def mixedP : QInput := ⟨23104065205033067, 25000000000000000⟩
private def mixedY : QInput := ⟨78511931972952687, 100000000000000000⟩
private def mixedZ : QInput :=
  ⟨27906833776273592503728404371, 125000000000000000000000000000⟩
private def mixedRadius : QInput :=
  ⟨77209539290942994, 1000000000000000000000000000000⟩
private def mixedVector : Fin 3 → QInput :=
  ![⟨1, 1⟩, ⟨8597328264764681, 1000000⟩,
    ⟨1528038021384231, 500000000000⟩]

private theorem mixed_transfer_exact :
    transferRatCheck
      (transferQInputRat mixedP * transferQInputRat mixedY)
      (transferQInputRat mixedZ) (transferQInputRat mixedRadius)
      (fun i => transferQInputRat (mixedVector i)) = true := by
  decide +kernel

private def mixedA0 : QInput := ⟨504260003, 524288000⟩
private def mixedA1 : QInput := ⟨5042610029, 5242880000⟩
private def mixedX0 : QInput := ⟨14092855623973879, 15728640000000000⟩
private def mixedX1 : QInput := ⟨112, 125⟩

private theorem mixed00 :
    scalableVertexCheck mixedSlope mixedIntercept mixedP mixedY
      mixedRadius mixedZ mixedA0 mixedX0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem mixed01 :
    scalableVertexCheck mixedSlope mixedIntercept mixedP mixedY
      mixedRadius mixedZ mixedA0 mixedX1 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem mixed10 :
    scalableVertexCheck mixedSlope mixedIntercept mixedP mixedY
      mixedRadius mixedZ mixedA1 mixedX0 50
      (-100000000000000000000000) = true := by decide +kernel
private theorem mixed11 :
    scalableVertexCheck mixedSlope mixedIntercept mixedP mixedY
      mixedRadius mixedZ mixedA1 mixedX1 50
      (-100000000000000000000000) = true := by decide +kernel

private theorem mixed_pointRate {α x : ℝ}
    (hα : α ∈ Set.Icc mixedA0.real mixedA1.real)
    (hx : x ∈ Set.Icc mixedX0.real mixedX1.real) :
    ScalableDensePointRate α x mixedSlope.real mixedIntercept.real
      (1 / 10000000) 1 := by
  have h := scalable_dense_checked_leaf mixedSlope mixedIntercept
    mixedP mixedY mixedRadius mixedZ mixedA0 mixedA1 mixedX0 mixedX1
    mixedVector 1
    (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide)
    (by norm_num [mixedP, QInput.real])
    (by norm_num [mixedP, QInput.real])
    (by norm_num [mixedY, QInput.real])
    (by norm_num [mixedY, QInput.real])
    (by norm_num [mixedRadius, QInput.real])
    (by norm_num [mixedZ, QInput.real])
    (by norm_num [mixedZ, QInput.real])
    (by norm_num [mixedA0, QInput.real])
    (by norm_num [mixedA1, QInput.real])
    (by norm_num [mixedX0, QInput.real])
    (by norm_num [mixedX1, QInput.real])
    (by norm_num : (0 : ℝ) ≤ 1)
    (by intro i; fin_cases i <;> norm_num [mixedVector, QInput.real])
    mixed_transfer_exact mixed00 mixed01 mixed10 mixed11 hα hx
  simpa [mixedVector, QInput.real] using h

theorem scalable_dense_mixed_pilot_coverage {α x : ℝ}
    (hα : α ∈ Set.Icc mixedA0.real mixedA1.real)
    (hx : x ∈ Set.Icc mixedX0.real mixedX1.real) :
    ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧
      ScalableDensePointRate α x p.1 p.2 (1 / 10000000) 1 := by
  refine ⟨((-833 / 500 : ℚ), (373189 / 250000 : ℚ)),
    Spin.Majorant.member38, ?_⟩
  simpa [mixedSlope, mixedIntercept, QInput.real] using
    (mixed_pointRate hα hx)

#print axioms mixed_transfer_exact
#print axioms mixed00
#print axioms mixed11
#print axioms scalable_dense_mixed_pilot_coverage

end Spin.UniformMixer
