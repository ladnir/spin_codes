import SpinCodes.Turbo.RS8ReferenceStepArray22
import SpinCodes.UniformMixer.ScalableRouteAverage

/-!
# Exact full-state kernel for the Turbo reference route

This module carries the checked physical-step input law into the
chronological 22-state weighted kernel. The independent regional input
laws can be averaged separately, but their kernels are multiplied in
region order so that the state persists across every boundary.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

def turboRegionStepStream22 {R : ℕ}
    (words : Fin R → PacketStep16) : List Input64 :=
  List.ofFn (fun step : Fin R => packetStepInput22 (words step))

theorem inputListKernel_kernelListProduct22 (z : ℝ)
    (xs : List Input64) :
    inputListKernel z xs =
      kernelListProduct (xs.map (stepKernel z)) := by
  induction xs with
  | nil => rfl
  | cons x tail ih =>
      simp only [List.map_cons, inputListKernel, kernelListProduct]
      rw [ih]

theorem inputListKernel_turbo_steps22 {R : ℕ}
    (z : ℝ) (words : Fin R → PacketStep16) :
    Matrix.of (inputListKernel z (turboRegionStepStream22 words)) =
      (List.ofFn (fun step : Fin R =>
        Matrix.of (stepKernel z (packetStepInput22 (words step))))).prod := by
  rw [inputListKernel_kernelListProduct22,
    kernelListProduct_matrix]
  simp only [turboRegionStepStream22, List.map_ofFn,
    Function.comp_def]
  rfl

theorem flatInput_step_stream22 {R : ℕ}
    (regions : Fin 64 → Fin (R * 16) → Packet4) :
    List.ofFn (flatInput regions) =
      (List.ofFn (fun region : Fin 64 =>
        turboRegionStepStream22
          ((regionalStepWords22 R regions) region))).flatten := by
  rw [List.ofFn_mul]
  apply congrArg List.flatten
  apply congrArg List.ofFn
  funext region
  apply congrArg List.ofFn
  funext step
  have heq : (⟨region.val * R + step.val,
      by
        have hfirst : region.val * R + step.val <
            (region.val + 1) * R := by
          rw [add_mul, one_mul]
          omega
        exact hfirst.trans_le
          (Nat.mul_le_mul_right _ (Nat.succ_le_of_lt region.isLt))⟩ :
        Fin (64 * R)) =
      (finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R))
        (region, step) := by
    ext
    ac_rfl
  rw [heq]
  simpa only [Equiv.symm_apply_apply] using
    flatInput_eq_regionalStepWords22 regions
      ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R))
        (region, step))

theorem inputListKernel_turbo_regions22 {R : ℕ}
    (z : ℝ) (regions : Fin 64 → Fin (R * 16) → Packet4) :
    Matrix.of (inputListKernel z (List.ofFn (flatInput regions))) =
      (List.ofFn (fun region : Fin 64 =>
        Matrix.of (inputListKernel z
          (turboRegionStepStream22
            ((regionalStepWords22 R regions) region))))).prod := by
  rw [flatInput_step_stream22, inputListKernel_flatten,
    kernelListProduct_matrix]
  rfl

theorem turboSelectedKernel22_matrix (j : Fin 17) (z : ℝ) :
    matrixExpect22 (selectedStepLaw22 j)
      (fun x => Matrix.of (stepKernel z (packetStepInput22 x))) =
        Matrix.of (turboSelectedKernel22 j z) := by
  rfl

theorem turboReferenceStepProduct22 {R : ℕ}
    (counts : Fin R → Fin 17) (z : ℝ) :
    matrixExpect22 (Spin.piPMF (fun step : Fin R =>
        selectedStepLaw22 (counts step)))
      (fun words => Matrix.of
        (inputListKernel z (turboRegionStepStream22 words))) =
      (List.ofFn (fun step : Fin R =>
        Matrix.of (turboSelectedKernel22 (counts step) z))).prod := by
  have hpoint (words : Fin R → PacketStep16) :
      Matrix.of (inputListKernel z
        (turboRegionStepStream22 words)) =
      (List.ofFn (fun step : Fin R =>
        Matrix.of (stepKernel z
          (packetStepInput22 (words step))))).prod :=
    inputListKernel_turbo_steps22 z words
  simp_rw [hpoint]
  calc
    _ = (List.ofFn (fun step : Fin R =>
          matrixExpect22 (selectedStepLaw22 (counts step))
            (fun x => Matrix.of
              (stepKernel z (packetStepInput22 x))))).prod := by
        exact matrixExpect22_pi_product
          (fun step : Fin R => selectedStepLaw22 (counts step))
          (fun _ x => Matrix.of (stepKernel z (packetStepInput22 x)))
    _ = _ := by simp only [turboSelectedKernel22_matrix]

theorem matrixExpect22_map_turbo {α β : Type*}
    [Fintype α] [Fintype β]
    (P : Spin.FinPMF α) (f : α → β)
    (K : β → Matrix State22 State22 ℝ) :
    matrixExpect22 (P.map f) K =
      matrixExpect22 P (fun x => K (f x)) := by
  ext q r
  exact Spin.FinPMF.expect_map P f (fun x => K x q r)

theorem matrixExpect22_bind_turbo {α β : Type*}
    [Fintype α] [Fintype β]
    (P : Spin.FinPMF α) (F : α → Spin.FinPMF β)
    (K : β → Matrix State22 State22 ℝ) :
    matrixExpect22 (P.bind F) K =
      matrixExpect22 P (fun x => matrixExpect22 (F x) K) := by
  ext q r
  exact Spin.FinPMF.expect_bind P F (fun x => K x q r)

def turboReferenceAveragedRegion22 (R : ℕ)
    (S : Finset (Fin (R * 16))) (z : ℝ) :
    Matrix State22 State22 ℝ :=
  matrixExpect22 (regionalStepProfileLaw22 R S)
    (fun counts =>
      (List.ofFn (fun step : Fin R =>
        Matrix.of (turboSelectedKernel22 (counts step) z))).prod)

theorem turboReferenceOneRegionKernel22 (R : ℕ)
    (S : Finset (Fin (R * 16))) (z : ℝ) :
    matrixExpect22 (turboReferenceOneRegionLaw22 R S)
      (fun row => Matrix.of (inputListKernel z
        (turboRegionStepStream22
          ((regionalStepWordEquiv22 R) row)))) =
      turboReferenceAveragedRegion22 R S z := by
  calc
    _ = matrixExpect22
          ((turboReferenceOneRegionLaw22 R S).map
            (regionalStepWordEquiv22 R))
          (fun words => Matrix.of (inputListKernel z
            (turboRegionStepStream22 words))) := by
          rw [matrixExpect22_map_turbo]
    _ = matrixExpect22
          ((regionalStepProfileLaw22 R S).bind
            (fun counts => Spin.piPMF (fun step : Fin R =>
              selectedStepLaw22 (counts step))))
          (fun words => Matrix.of (inputListKernel z
            (turboRegionStepStream22 words))) := by
          rw [turboReferenceOneRegionLaw22_step_mixture]
    _ = _ := by
      rw [matrixExpect22_bind_turbo]
      unfold turboReferenceAveragedRegion22
      apply congrArg
      funext counts
      exact turboReferenceStepProduct22 counts z

/-- Exact averaged weighted kernel of the actual reference route.
Each factor is a full `2^22` by `2^22` state kernel; matrix
multiplication, not scalar multiplication, carries the state across
the 64 region boundaries. -/
theorem turboReferenceRoutedKernel22 (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (z : ℝ) :
    matrixExpect22 (turboReferenceRoutedLaw22 R points S)
      (fun regions => Matrix.of (inputListKernel z
        (List.ofFn (flatInput regions)))) =
      (List.ofFn (fun _ : Fin 64 =>
        turboReferenceAveragedRegion22 R S z)).prod := by
  rw [turboReferenceRoutedLaw22_product]
  have hpoint (regions : Fin 64 → Fin (R * 16) → Packet4) :
      Matrix.of (inputListKernel z
        (List.ofFn (flatInput regions))) =
      (List.ofFn (fun region : Fin 64 =>
        Matrix.of (inputListKernel z
          (turboRegionStepStream22
            ((regionalStepWords22 R regions) region))))).prod :=
    inputListKernel_turbo_regions22 z regions
  simp_rw [hpoint]
  calc
    _ = (List.ofFn (fun region : Fin 64 =>
          matrixExpect22 (turboReferenceOneRegionLaw22 R S)
            (fun row => Matrix.of (inputListKernel z
              (turboRegionStepStream22
                ((regionalStepWordEquiv22 R) row)))))).prod := by
          exact matrixExpect22_pi_product
            (fun _ : Fin 64 => turboReferenceOneRegionLaw22 R S)
            (fun _ row => Matrix.of (inputListKernel z
              (turboRegionStepStream22
                ((regionalStepWordEquiv22 R) row))))
    _ = _ := by simp only [turboReferenceOneRegionKernel22]

theorem turboRemainingMoment22 (R : ℕ) (z : ℝ)
    (words : TurboOuterWord22 R) :
    (turboRemainingLaw22 R).expect (fun seed =>
      z ^ sampledOutputWeight
        (flatInput (turboRouteRegions22 words seed.1)) seed.2 0) =
      (turboRegionPermLaw22 R).expect (fun perms =>
        ∑ q : State22,
          inputListKernel z
            (List.ofFn (flatInput
              (turboRouteRegions22 words perms))) 0 q) := by
  rw [turboRemainingLaw22, FinPMF.expect_prod]
  apply congrArg
  funext perms
  exact sampledInputMoment_kernel z
    (flatInput (turboRouteRegions22 words perms)) 0

/-- The reference's exact moment under the actual regional permutations
and independent sampled state maps. This is the full-state semantic
target for finite Turbo numerical certificates. -/
theorem turboReferenceMoment22 (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (z : ℝ) :
    (turboReferenceArrayLaw22 R points S).expect (fun words =>
      (turboRemainingLaw22 R).expect (fun seed =>
        z ^ sampledOutputWeight
          (flatInput (turboRouteRegions22 words seed.1)) seed.2 0)) =
      ∑ q : State22,
        (List.ofFn (fun _ : Fin 64 =>
          turboReferenceAveragedRegion22 R S z)).prod 0 q := by
  simp_rw [turboRemainingMoment22]
  calc
    _ = ((turboReferenceArrayLaw22 R points S).prod
          (turboRegionPermLaw22 R)).expect (fun pair =>
            ∑ q : State22,
              inputListKernel z
                (List.ofFn (flatInput
                  (turboRouteRegions22 pair.1 pair.2))) 0 q) := by
          rw [FinPMF.expect_prod]
    _ = (turboReferenceRoutedLaw22 R points S).expect
          (fun regions => ∑ q : State22,
            inputListKernel z
              (List.ofFn (flatInput regions)) 0 q) := by
          rw [turboReferenceRoutedLaw22, FinPMF.expect_map]
    _ = ∑ q : State22,
          matrixExpect22 (turboReferenceRoutedLaw22 R points S)
            (fun regions => Matrix.of (inputListKernel z
              (List.ofFn (flatInput regions)))) 0 q := by
          rw [FinPMF.expect_sum]
          rfl
    _ = _ := by rw [turboReferenceRoutedKernel22]

theorem turboReferenceFailure22_joint (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (d : ℕ) :
    (turboReferenceArrayLaw22 R points S).expect
      (fun words => turboRemainingFailure22 words d) =
      ((turboReferenceArrayLaw22 R points S).prod
        (turboRemainingLaw22 R)).prob (fun pair =>
          sampledOutputWeight
            (flatInput
              (turboRouteRegions22 pair.1 pair.2.1))
            pair.2.2 0 ≤ d) := by
  rw [FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod]
  unfold turboRemainingFailure22
  simp_rw [FinPMF.prob_eq_expect_indicator]

theorem turboReferenceFailure22_le_kernel (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (turboReferenceArrayLaw22 R points S).expect
      (fun words => turboRemainingFailure22 words d) ≤
      (∑ q : State22,
        (List.ofFn (fun _ : Fin 64 =>
          turboReferenceAveragedRegion22 R S z)).prod 0 q) /
        z ^ d := by
  rw [turboReferenceFailure22_joint]
  apply FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
  rw [FinPMF.expect_prod]
  exact le_of_eq (turboReferenceMoment22 R points S z)

theorem turboSetup22_activeSet_le_reference_kernel
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (turboOuterLaw22 R).expect (fun outerSeed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d
          (turboOuterWord22 points outerSeed msg)) ≤
      turboOldUniformFactor22 ^ S.card *
        ((∑ q : State22,
          (List.ofFn (fun _ : Fin 64 =>
            turboReferenceAveragedRegion22 R S z)).prod 0 q) /
          z ^ d) := by
  calc
    _ ≤ turboOldUniformFactor22 ^ S.card *
          (turboReferenceArrayLaw22 R points S).expect
            (fun words => turboRemainingFailure22 words d) :=
      turboSetup22_activeSet_moment_le_reference_expect
        R points hpoints S d
    _ ≤ _ := mul_le_mul_of_nonneg_left
      (turboReferenceFailure22_le_kernel R points S d z hz0 hz1)
      (pow_nonneg turboOldUniformFactor22_pos.le _)

#print axioms flatInput_step_stream22
#print axioms inputListKernel_turbo_regions22
#print axioms inputListKernel_turbo_steps22
#print axioms turboReferenceStepProduct22
#print axioms turboReferenceOneRegionKernel22
#print axioms turboReferenceRoutedKernel22
#print axioms turboReferenceMoment22
#print axioms turboReferenceFailure22_le_kernel
#print axioms turboSetup22_activeSet_le_reference_kernel

end Spin.Turbo
