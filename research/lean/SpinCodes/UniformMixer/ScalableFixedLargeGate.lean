import SpinCodes.UniformMixer.ScalableFixedLargeEnvelope
import SpinCodes.Structured.ConcreteNativeFixedLargeDefs

/-! The old spacing estimate bounds Scalable's fixed-Q continuum kernel. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Spin.Structured.ConcreteEncoder
  Spin.Structured.ConcreteNativeFamily
  Spin.Structured.ConcreteFixedNumeric Finset MeasureTheory

theorem scalableLargeNorm_proved {Q : ℕ} (hQ : 3 ≤ Q)
    (u : FixedFugacityChoice Q) :
    scalableTwoRowNormLe (3 / 1600) (fixedNormCost u)
      (scalableFugacityTwo (scalableLargeTheta Q)
        (fixedFugacity u)) := by
  let v : ℝ := 3 / 1600
  let c : ℝ := Real.exp ((Q : ℝ) * paperLargeExponent) *
    ∏ i, Spin.mv (127 / 250) v (1 / 524287)
      (fixedFugacity u i)
  have hu : ∀ i, 0 ≤ fixedFugacity u i :=
    fun i => (fixedFugacity_pos u i).le
  have hscalar (S : Finset (Fin (Q + 1))) :
      (Q.factorial : ℝ) *
        (∫ x in orderedSiteDomain Q,
          Real.exp (-(paperLargeTheta Q * epochMean / 128) *
            liveDuration S x)) ≤
      Real.exp ((Q : ℝ) * paperLargeExponent) *
        (127 / 250 : ℝ) ^ S.card := by
    simpa only [paperLargeTheta_epochMean, neg_mul] using
      paper_liveDuration_laplace_le hQ S
  have h0 := enlargedInsertionKernel_action_bound
    (by norm_num : (0 : ℝ) ≤ 127 / 250)
    (by norm_num : (127 / 250 : ℝ) ≤ 1)
    (by norm_num : (0 : ℝ) < v)
    (Real.exp_nonneg _) hscalar (fixedFugacity u) hu 0
  have h1 := enlargedInsertionKernel_action_bound
    (by norm_num : (0 : ℝ) ≤ 127 / 250)
    (by norm_num : (127 / 250 : ℝ) ≤ 1)
    (by norm_num : (0 : ℝ) < v)
    (Real.exp_nonneg _) hscalar (fixedFugacity u) hu 1
  have h0' :
      enlargedInsertionKernel (paperLargeTheta Q)
          (fixedFugacity u) 0 0 +
        enlargedInsertionKernel (paperLargeTheta Q)
          (fixedFugacity u) 0 1 * v ≤ c := by
    simpa [Matrix.mulVec, dotProduct,
      Fin.sum_univ_two, twoWeight, c] using h0
  have h1' :
      enlargedInsertionKernel (paperLargeTheta Q)
          (fixedFugacity u) 1 0 +
        enlargedInsertionKernel (paperLargeTheta Q)
          (fixedFugacity u) 1 1 * v ≤ c * v := by
    simpa [Matrix.mulVec, dotProduct,
      Fin.sum_univ_two, twoWeight, c] using h1
  have he (i j : Fin 2) :=
    scalableLargeFugacity_le_enlarged (fixedFugacity u) hu i j
  have hv : 0 ≤ v := by norm_num [v]
  have hz := add_le_add (he 0 0)
    (mul_le_mul_of_nonneg_right (he 0 1) hv)
  have hl := add_le_add (he 1 0)
    (mul_le_mul_of_nonneg_right (he 1 1) hv)
  unfold scalableTwoRowNormLe
  constructor
  · have hh :
        scalableFugacityTwo (scalableLargeTheta Q)
            (fixedFugacity u) false false +
          scalableFugacityTwo (scalableLargeTheta Q)
            (fixedFugacity u) false true * v ≤ c := by
      have hrow := hz.trans h0'
      simpa [toFinTwo_apply, Matrix.mulVec,
        dotProduct, Fin.sum_univ_two, twoWeight,
        finTwoEquiv, v, mul_one, add_zero] using hrow
    simpa only [fixedNormCost, gThree_eq, c, v] using hh
  · have hh :
        scalableFugacityTwo (scalableLargeTheta Q)
            (fixedFugacity u) true false +
          scalableFugacityTwo (scalableLargeTheta Q)
            (fixedFugacity u) true true * v ≤ c * v := by
      have hrow := hl.trans h1'
      simpa [toFinTwo_apply, Matrix.mulVec,
        dotProduct, Fin.sum_univ_two, twoWeight,
        finTwoEquiv, v, mul_one, add_zero] using hrow
    simpa only [fixedNormCost, gThree_eq, c, v] using hh

end Spin.UniformMixer
