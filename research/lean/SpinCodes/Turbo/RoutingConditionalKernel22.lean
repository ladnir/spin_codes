import SpinCodes.Turbo.RoutingConditionalSteps22
import SpinCodes.Turbo.RS8ReferenceKernel22

/-!
# Actual routing-fiber first moment as a chronological physical product

The conditional input law is proved from the original setup, including
selected zero values. Every kernel is multiplied in region-major order;
the state is never reset or summed at an intermediate region boundary.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.UniformMixer

/-- Full-state conditional kernel, before any eleven-coordinate envelope. -/
def turboCountKernel22 {R : ℕ} (counts : TurboPotentialCounts22 R) (z : ℝ) :
    Matrix State22 State22 ℝ :=
  (List.ofFn (fun region : Fin 64 =>
    (List.ofFn (fun step : Fin R =>
      Matrix.of (turboSelectedKernel22 (counts region step) z))).prod)).prod

theorem turboStepArrayInput22_region_kernel {R : ℕ}
    (words : TurboStepArray22 R) (z : ℝ) :
    Matrix.of (inputListKernel z (List.ofFn (turboStepArrayInput22 words))) =
      (List.ofFn (fun region : Fin 64 =>
        Matrix.of (inputListKernel z (turboRegionStepStream22 (words region))))).prod := by
  let regions : Fin 64 → Fin (R * 16) → Packet4 :=
    fun region => (regionalStepWordEquiv22 R).symm (words region)
  have hre : regionalStepWords22 R regions = words := by
    funext region
    exact (regionalStepWordEquiv22 R).apply_symm_apply (words region)
  have hi : flatInput regions = turboStepArrayInput22 words := by
    funext index
    rw [flatInput_eq_regionalStepWords22]
    change packetStepInput22 ((regionalStepWords22 R regions) _ _) = _
    rw [hre]
    rfl
  rw [← hi, inputListKernel_turbo_regions22, hre]

/-- Independent physical inputs are averaged only after the actual routing
counts are fixed. Matrix multiplication retains the complete boundary state. -/
theorem turboCountStepArrayKernel22 {R : ℕ}
    (counts : TurboPotentialCounts22 R) (z : ℝ) :
    matrixExpect22 (turboCountStepArrayLaw22 counts)
      (fun words => Matrix.of (inputListKernel z
        (List.ofFn (turboStepArrayInput22 words)))) = turboCountKernel22 counts z := by
  simp_rw [turboStepArrayInput22_region_kernel]
  rw [turboCountStepArrayLaw22]
  rw [matrixExpect22_pi_product
    (fun region : Fin 64 => Spin.piPMF (fun step : Fin R =>
      selectedStepLaw22 (counts region step)))
    (fun (_ : Fin 64) (word : Fin R → PacketStep16) =>
      Matrix.of (inputListKernel z (turboRegionStepStream22 word)))]
  simp only [turboReferenceStepProduct22, turboCountKernel22]

/-- Exact conditional generating moment of the selected-step experiment. -/
theorem turboCountStepArrayMoment22 {R : ℕ}
    (counts : TurboPotentialCounts22 R) (z : ℝ) :
    ((turboCountStepArrayLaw22 counts).prod (turboStateLaw22 R)).expect
      (fun pair => z ^ sampledOutputWeight (turboStepArrayInput22 pair.1) pair.2 0) =
      ∑ q : State22, turboCountKernel22 counts z 0 q := by
  rw [FinPMF.expect_prod]
  have he (words : TurboStepArray22 R) :
      (turboStateLaw22 R).expect
        (fun maps => z ^ sampledOutputWeight (turboStepArrayInput22 words) maps 0) =
      ∑ q : State22, inputListKernel z
        (List.ofFn (turboStepArrayInput22 words)) 0 q :=
    sampledInputMoment_kernel z (turboStepArrayInput22 words) 0
  simp_rw [he]
  rw [FinPMF.expect_sum]
  change (∑ q : State22, matrixExpect22 (turboCountStepArrayLaw22 counts)
    (fun words => Matrix.of (inputListKernel z
      (List.ofFn (turboStepArrayInput22 words)))) 0 q) = _
  rw [turboCountStepArrayKernel22]

/-- Original conditional active-set mean, dominated by the actual selected
physical kernels. The remaining numerical task is to bound this product. -/
theorem turboActiveBadCount22_fiber_le_kernel
    (R : ℕ) (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R)
    (hp : 0 < (turboPotentialCountLaw22 R S).p counts)
    (d : ℕ) (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (turboSetupCountFiber22 R points S counts).expect
      (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) ≤
      turboOldUniformFactor22 ^ S.card *
        ((∑ q : State22, turboCountKernel22 counts z 0 q) / z ^ d) := by
  apply (turboActiveBadCount22_fiber_le_reference R points hpoints S counts d).trans
  apply mul_le_mul_of_nonneg_left _ (pow_nonneg turboOldUniformFactor22_pos.le _)
  rw [turboReferenceFiberFailure22_eq_selected R points S counts hp d]
  apply FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
  exact (turboCountStepArrayMoment22 counts z).le

#print axioms turboCountStepArrayKernel22
#print axioms turboCountStepArrayMoment22
#print axioms turboActiveBadCount22_fiber_le_kernel
end Spin.Turbo

