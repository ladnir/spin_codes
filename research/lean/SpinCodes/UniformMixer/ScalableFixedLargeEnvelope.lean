import SpinCodes.UniformMixer.ScalableFixedLargeInsertion
import SpinCodes.Structured.ConcreteFixedLargeNorm

/-! Scalable's all-row insertion kernel lies below the old enlarged kernel. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Spin.Structured.ConcreteEncoder
  Finset MeasureTheory

attribute [local instance] Classical.propDecidable

theorem scalableFugacityTwo_le_enlarged {Q : ℕ}
    (θ θOld : ℝ)
    (hγ : regionGamma θ = θOld * epochMean / 128)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k)
    (i j : Fin 2) :
    toFinTwo (scalableFugacityTwo θ u) i j ≤
      enlargedInsertionKernel θOld u i j := by
  rw [scalableFugacityTwo_eq_insertion]
  unfold scalableAveragedInsertionKernel enlargedInsertionKernel
  rw [← hγ]
  apply sum_le_sum
  intro π hπ
  have hi (P : Matrix (Fin 2) (Fin 2) ℝ) : IntegrableOn
      (fun x => fullSiteProduct (regionGamma θ) P x
        (fun k => u (π.symm k)) i j)
        (orderedSiteDomain Q) := by
    have hc := (continuous_apply j).comp
      ((continuous_apply i).comp
        (continuous_fullSiteProduct
          (regionGamma θ) P (fun k => u (π.symm k))))
    apply hc.integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun k => (hx.1 k).1,
      fun k => (hx.1 k).2.le⟩
  apply setIntegral_mono_on (hi _) (hi _)
    (orderedSiteDomain_measurable Q)
  intro x hx
  exact fullGapProduct_mono (regionGamma θ)
    scalableImpulseFin (Spin.Pplus (1 / 524287))
    (fun k => positionGaps x k.castSucc)
    (fun k => u (π.symm k))
    (positionGaps x (Fin.last Q))
    (fun k => hu (π.symm k))
    scalable_impulse_nonneg plus_impulse_nonneg
    scalableImpulseFin_le_plus i j

def scalableLargeTheta (Q : ℕ) : ℝ :=
  (Q : ℝ) * (133 / 125) * scalableThetaOne

theorem scalableLargeTheta_gamma (Q : ℕ) :
    regionGamma (scalableLargeTheta Q) =
      paperLargeTheta Q * epochMean / 128 := by
  unfold scalableLargeTheta paperLargeTheta
  norm_num [regionGamma, production_liveWeightAverage,
    scalableThetaOne, epochMean]
  ring

theorem scalableLargeFugacity_le_enlarged {Q : ℕ}
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k)
    (i j : Fin 2) :
    toFinTwo (scalableFugacityTwo
      (scalableLargeTheta Q) u) i j ≤
        enlargedInsertionKernel (paperLargeTheta Q) u i j :=
  scalableFugacityTwo_le_enlarged
    (scalableLargeTheta Q) (paperLargeTheta Q)
    (scalableLargeTheta_gamma Q) u hu i j

end Spin.UniformMixer
