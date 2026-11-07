import SpinCodes.UniformMixer.ScalableFixedRoutedProfile
import SpinCodes.UniformMixer.ScalableRouteAverage
import SpinCodes.Structured.ConcreteFixedFairProduct

/-! Fair active outer rows under the actual Scalable regional kernel. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteMarked
  Spin.Structured.ConcreteRoute Spin.Structured.Placement
  Spin.Structured.Routing Filter

def scalableFairRegionKernel (m : ℕ) (θ : ℝ)
    (S : Finset (Fin (Lsched m))) : Matrix State22 State22 ℝ :=
  matrixExpect22 (regionLaw S) (fun T =>
    Matrix.of (chronologicalRegionKernel
      (Real.exp (-(θ / (Lsched m : ℝ))))
      (scalableRegionSupport m T)))

theorem scalableFairRegionKernel_eq_fugacity {Q : ℕ}
    (m : ℕ) (θ : ℝ) (marks : Fin Q ↪ Fin (Lsched m))
    (q r : State22) :
    scalableFairRegionKernel m θ (Finset.univ.map marks) q r =
      scalableFugacityKernel θ marks (fun _ => 1) q r /
        (2 : ℝ) ^ Q := by
  unfold scalableFairRegionKernel matrixExpect22
  rw [← shuffled_fairOnMarks, shuffledLaw,
    FinPMF.expect_bind, fairOnMarks_expect_uniform_sum]
  simp only [scalableFugacityKernel, markCoefficient,
    Finset.prod_const_one, one_mul, scalableAveragedRegion,
    matrixExpect22]

theorem scalable_fair_region_product (m : ℕ) (θ : ℝ)
    (S : Finset (Fin (Lsched m))) :
    matrixExpect22 (regionsLaw S (bsched m)) (fun regions =>
      Matrix.of (inputListKernel
        (Real.exp (-(θ / (Lsched m : ℝ))))
        (List.ofFn (scalableFlatInput m regions)))) =
      (List.ofFn (fun _j : Fin (bsched m) =>
        scalableFairRegionKernel m θ S)).prod := by
  have hpoint (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
      Matrix.of (inputListKernel
        (Real.exp (-(θ / (Lsched m : ℝ))))
        (List.ofFn (scalableFlatInput m regions))) =
      (List.ofFn (fun j : Fin (bsched m) =>
        Matrix.of (chronologicalRegionKernel
          (Real.exp (-(θ / (Lsched m : ℝ))))
          (scalableRegionSupport m (regions j))))).prod := by
    rw [scalable_kernel_region_product, kernelListProduct_matrix]
    simp only [List.map_ofFn]
    rfl
  simp only [regionsLaw, hpoint, scalableFairRegionKernel]
  exact matrixExpect22_pi_product
    (fun _j : Fin (bsched m) => regionLaw S)
    (fun _j T => Matrix.of (chronologicalRegionKernel
      (Real.exp (-(θ / (Lsched m : ℝ))))
      (scalableRegionSupport m T)))

theorem scalable_fair_region_moment (m : ℕ) (θ : ℝ)
    (S : Finset (Fin (Lsched m))) :
    (regionsLaw S (bsched m)).expect
      (scalableRoutedMoment m
        (Real.exp (-(θ / (Lsched m : ℝ))))) =
      ∑ r : State22,
        (List.ofFn (fun _j : Fin (bsched m) =>
          scalableFairRegionKernel m θ S)).prod 0 r := by
  change (regionsLaw S (bsched m)).expect (fun regions =>
      ∑ r : State22,
        inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
          (List.ofFn (scalableFlatInput m regions)) 0 r) = _
  rw [FinPMF.expect_sum]
  apply Finset.sum_congr rfl
  intro r _
  exact congrFun (congrFun
    (scalable_fair_region_product m θ S) 0) r

theorem scalable_fair_region_margin {Q : ℕ}
    {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ)
    (hK : scalableTwoRowNormLe v c
      (scalableFugacityTwo θ (fun _ : Fin Q => 1))) :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin Q ↪ Fin (Lsched m),
        (regionsLaw (Finset.univ.map marks) (bsched m)).expect
          (scalableRoutedMoment m
            (Real.exp (-(θ / (Lsched m : ℝ))))) ≤
          (((c + δ) / (2 : ℝ) ^ Q) ^ bsched m /
            min 1 v) := by
  have hm₀ : 0 < min 1 v := lt_min (by norm_num) hv
  filter_upwards [scalable_fugacity_weighted_margin
    hθ hδ (fun _ : Fin Q => 1) (by intro i; norm_num)
    (scalableStateWeight v) (scalableStateWeight_pos hv)
    hm₀ scalableStateWeight_min
    (scalableFugacityContinuum_weighted_bound θ _ hK)]
    with m hm
  intro marks
  rw [scalable_fair_region_moment]
  have hk (j : Fin (bsched m)) :
      Spin.Structured.FiniteKernel.WeightedBound
        (scalableStateWeight v) ((c + δ) / (2 : ℝ) ^ Q)
        (scalableFairRegionKernel m θ (Finset.univ.map marks)) := by
    constructor
    · intro q r
      rw [scalableFairRegionKernel_eq_fugacity]
      exact div_nonneg
        (scalableFugacityKernel_nonneg hθ
          (by intro i; norm_num)
          marks q r) (by positivity)
    · intro q
      simp only [scalableFairRegionKernel_eq_fugacity,
        div_mul_eq_mul_div, ← Finset.sum_div]
      simpa only [div_mul_eq_mul_div] using
        div_le_div_of_nonneg_right ((hm marks).2 q)
          (by positivity : 0 ≤ (2 : ℝ) ^ Q)
  have hp := Spin.Structured.FiniteKernel.WeightedBound.ofFn
    (by positivity : 0 ≤ (c + δ) / (2 : ℝ) ^ Q)
    (fun _j : Fin (bsched m) =>
      scalableFairRegionKernel m θ (Finset.univ.map marks)) hk
  simpa [scalableStateWeight] using
    hp.total hm₀ scalableStateWeight_min 0

end Spin.UniformMixer
