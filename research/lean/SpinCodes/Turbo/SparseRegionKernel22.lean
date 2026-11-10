import SpinCodes.Turbo.SparseSupportArray22
import SpinCodes.Turbo.RS8ReferenceKernel22

/-!
# Continuous exact sparse regional kernel

Audit: `RS8ReferenceKernel22.turboReferenceRoutedKernel22` and its moment
and lower-tail lemmas multiply averaged full-state kernels chronologically.
Here each region uses its exact nonzero support-class law. Reusing that
argument preserves both shared regional permutations and every boundary
state; no region is reset and no failure probabilities are multiplied.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

def sparseExactRegionKernel22 (R : ℕ) (T : Finset (Fin (R * 16))) (z : ℝ) :
    Matrix State22 State22 ℝ :=
  matrixExpect22 (sparseExactRegionLaw22 R T) (fun row =>
    Matrix.of (inputListKernel z (turboRegionStepStream22 ((regionalStepWordEquiv22 R) row))))

theorem sparseExactRegionKernel22_eq_of_card (R : ℕ)
    (T U : Finset (Fin (R * 16))) (h : T.card = U.card) (z : ℝ) :
    sparseExactRegionKernel22 R T z = sparseExactRegionKernel22 R U z := by
  unfold sparseExactRegionKernel22
  rw [sparseExactRegionLaw22_eq_of_card R T U h]

theorem sparseSupportRoutedKernel22 (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) (z : ℝ) :
    matrixExpect22 (sparseSupportRoutedLaw22 R supports)
      (fun regions => Matrix.of (inputListKernel z (List.ofFn (flatInput regions)))) =
      (List.ofFn (fun region => sparseExactRegionKernel22 R
        (sparseRegionalSupport22 supports region) z)).prod := by
  rw [sparseSupportRoutedLaw22_product]
  simp_rw [inputListKernel_turbo_regions22]
  exact matrixExpect22_pi_product
    (fun region : Fin 64 => sparseExactRegionLaw22 R (sparseRegionalSupport22 supports region))
    (fun _ row => Matrix.of (inputListKernel z
      (turboRegionStepStream22 ((regionalStepWordEquiv22 R) row))))

theorem sparseSupportMoment22 (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) (z : ℝ) :
    (sparseSupportArrayLaw22 R supports).expect (fun words =>
      (turboRemainingLaw22 R).expect (fun seed =>
        z ^ sampledOutputWeight
          (flatInput (turboRouteRegions22 words seed.1)) seed.2 0)) =
      ∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R
        (sparseRegionalSupport22 supports region) z)).prod 0 q := by
  simp_rw [turboRemainingMoment22]
  calc
    _ = ((sparseSupportArrayLaw22 R supports).prod (turboRegionPermLaw22 R)).expect
        (fun pair => ∑ q : State22, inputListKernel z
          (List.ofFn (flatInput (turboRouteRegions22 pair.1 pair.2))) 0 q) := by
      rw [FinPMF.expect_prod]
    _ = (sparseSupportRoutedLaw22 R supports).expect (fun regions =>
        ∑ q : State22, inputListKernel z (List.ofFn (flatInput regions)) 0 q) := by
      rw [sparseSupportRoutedLaw22, FinPMF.expect_map]
    _ = ∑ q : State22, matrixExpect22 (sparseSupportRoutedLaw22 R supports)
        (fun regions => Matrix.of (inputListKernel z (List.ofFn (flatInput regions)))) 0 q := by
      rw [FinPMF.expect_sum]
      rfl
    _ = _ := by rw [sparseSupportRoutedKernel22]

theorem sparseSupportFailure22_le_kernel (R : ℕ)
    (supports : Fin (R * 16) → Finset (Fin 64)) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (sparseSupportArrayLaw22 R supports).expect (fun words => turboRemainingFailure22 words d) ≤
      (∑ q : State22, (List.ofFn (fun region => sparseExactRegionKernel22 R
        (sparseRegionalSupport22 supports region) z)).prod 0 q) / z ^ d := by
  have hjoint : (sparseSupportArrayLaw22 R supports).expect
      (fun words => turboRemainingFailure22 words d) =
      ((sparseSupportArrayLaw22 R supports).prod (turboRemainingLaw22 R)).prob
        (fun pair => sampledOutputWeight
          (flatInput (turboRouteRegions22 pair.1 pair.2.1)) pair.2.2 0 ≤ d) := by
    rw [FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod]
    unfold turboRemainingFailure22
    simp_rw [FinPMF.prob_eq_expect_indicator]
  rw [hjoint]
  apply FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
  rw [FinPMF.expect_prod]
  exact (sparseSupportMoment22 R supports z).le

#print axioms sparseSupportMoment22
#print axioms sparseSupportFailure22_le_kernel
end Spin.Turbo
