import SpinCodes.Structured.ConcreteFixedSmallKernel

/-! Closed forms for the actual zero- and one-impulse continuum kernels. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory
attribute [local instance] Classical.propDecidable

theorem orderedSiteDomain_zero : orderedSiteDomain 0 = Set.univ := by
  ext x
  simp only [Set.mem_univ, iff_true]
  constructor
  · intro i; exact Fin.elim0 i
  · intro i; exact Fin.elim0 i

theorem continuumRegionKernel_zero (θ : ℝ) :
    continuumRegionKernel θ 0 = timeEmpty (θ * epochMean / 128) 1 := by
  ext i j
  simp only [continuumRegionKernel, Nat.factorial_zero, Nat.cast_one, one_mul, orderedSiteDomain_zero,
    Measure.restrict_univ]
  rw [Measure.volume_pi_eq_dirac (fun i : Fin 0 => Fin.elim0 i), integral_dirac]
  simp [siteProduct, simplexProduct, timeProduct, positionGaps, Fin.snoc]

theorem actual_integrand_one (γ t : ℝ) :
    timeEmpty γ t * impulseMatrix * timeEmpty γ (1-t) =
      !![0, Real.exp (-γ * (1-t));
        (1/524287) * Real.exp (-γ*t), (1-1/524287) * Real.exp (-γ)] := by
  have he := Spin.exp_prod_one γ t
  simp only [neg_mul] at he
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [timeEmpty, impulseMatrix, Matrix.mul_apply, Fin.sum_univ_two] <;> nlinarith

theorem continuumRegionKernel_one {θ : ℝ} (hγ : θ * epochMean / 128 ≠ 0) :
    continuumRegionKernel θ 1 =
      !![0, (1 - Real.exp (-(θ * epochMean / 128))) / (θ * epochMean / 128);
        (1/524287) * ((1 - Real.exp (-(θ * epochMean / 128))) / (θ * epochMean / 128)),
        (1-1/524287) * Real.exp (-(θ * epochMean / 128))] := by
  ext i j
  rw [continuumRegionKernel_one_integral]
  simp only [actual_integrand_one]
  fin_cases i <;> fin_cases j
  · simp
  · simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
    exact Spin.integral_expNeg_rev hγ
  · change (∫ t in (0:ℝ)..1, (1/524287:ℝ) * Real.exp (-(θ * epochMean / 128) * t)) =
      (1/524287:ℝ) * ((1-Real.exp (-(θ * epochMean / 128))) / (θ * epochMean / 128))
    rw [intervalIntegral.integral_const_mul, Spin.integral_expNeg hγ]
  · simp [intervalIntegral.integral_const]

end Spin.Structured.Placement


