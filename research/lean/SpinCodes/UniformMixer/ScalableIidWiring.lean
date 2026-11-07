import SpinCodes.UniformMixer.ScalableBernoulliStep
import SpinCodes.Structured.ConcreteReshape

/-!
# The iid route reference in chronological step order

The route comparison gives independent Bernoulli bits in region-major
coordinates. This module rewires those same bits into 64-bit encoder steps.
The wiring is a bijection, so it does not add independence assumptions.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute

def scalableBitWiring (m : ℕ) :
    (Fin (bsched m) × Fin (Lsched m)) ≃
      (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) :=
  finProdFinEquiv.trans
    ((finCongr (by rw [scalable_region_length]; ring)).trans
      finProdFinEquiv.symm)

theorem scalableBitWiring_regionBit (m : ℕ)
    (j : Fin (bsched m)) (s : Fin (scalableStepsPerRegion m))
    (c : Fin 64) :
    scalableBitWiring m
      (j, (finCongr (scalable_region_length m)).symm (regionBit s c)) =
        (finProdFinEquiv (j, s), c) := by
  apply finProdFinEquiv.injective
  rw [show finProdFinEquiv
      (scalableBitWiring m
        (j, (finCongr (scalable_region_length m)).symm (regionBit s c))) =
      finCongr (by rw [scalable_region_length]; ring)
        (finProdFinEquiv
          (j, (finCongr (scalable_region_length m)).symm
            (regionBit s c))) from
      finProdFinEquiv.apply_symm_apply _]
  apply Fin.ext
  change 64 * s.val + c.val + (Lsched m) * j.val =
    c.val + 64 * (s.val + scalableStepsPerRegion m * j.val)
  rw [scalable_region_length]
  ring

theorem scalableFlatInput_reshape (m : ℕ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    (fun i => supportInput (reshape (scalableBitWiring m) regions i)) =
      scalableFlatInput m regions := by
  funext i
  obtain ⟨⟨j, s⟩, rfl⟩ := finProdFinEquiv.surjective i
  funext c
  have he : (scalableBitWiring m).symm (finProdFinEquiv (j, s), c) =
      (j, (finCongr (scalable_region_length m)).symm (regionBit s c)) := by
    rw [← scalableBitWiring_regionBit, Equiv.symm_apply_apply]
  simp only [supportInput, mem_reshape, he, scalableFlatInput,
    Equiv.symm_apply_apply, scalableRegionInput, regionInput]
  simp [scalableRegionSupport]

theorem scalableFlatInput_iidLaw (m : ℕ) (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    (iidLaw (Lsched m) (bsched m) β hβ0 hβ1).map
      (scalableFlatInput m) =
    piPMF (fun _ : Fin (bsched m * scalableStepsPerRegion m) =>
      bernoulliInputLaw β hβ0 hβ1) := by
  let e := scalableBitWiring m
  let f : (Fin (bsched m * scalableStepsPerRegion m) → Finset (Fin 64)) →
      (Fin (bsched m * scalableStepsPerRegion m) → Input64) :=
    fun x i => supportInput (x i)
  calc
    _ = ((iidLaw (Lsched m) (bsched m) β hβ0 hβ1).map
        (reshape e)).map f := by
      rw [FinPMF.map_comp]
      congr 1
      funext regions
      exact (scalableFlatInput_reshape m regions).symm
    _ = (iidLaw 64 (bsched m * scalableStepsPerRegion m)
        β hβ0 hβ1).map f := by
      rw [reshape_iid]
    _ = _ := by
      change (piPMF (fun _ : Fin (bsched m * scalableStepsPerRegion m) =>
        poissonBinom (fun _ : Fin 64 => β) (fun _ => hβ0)
          (fun _ => hβ1))).map
        (fun x i => supportInput (x i)) = _
      rw [piPMF_map]
      rfl

theorem scalable_iid_moment_power (m : ℕ) (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    (iidLaw (Lsched m) (bsched m) β hβ0 hβ1).expect
      (scalableRoutedMoment m z) =
    ∑ r : State22,
      ((bernoulliStepKernel β z hβ0 hβ1) ^
        (bsched m * scalableStepsPerRegion m)) 0 r := by
  change (iidLaw (Lsched m) (bsched m) β hβ0 hβ1).expect
      (fun regions => ∑ r : State22,
        inputListKernel z (List.ofFn (scalableFlatInput m regions)) 0 r) = _
  calc
    _ = ((iidLaw (Lsched m) (bsched m) β hβ0 hβ1).map
        (scalableFlatInput m)).expect
        (fun xs => ∑ r : State22,
          inputListKernel z (List.ofFn xs) 0 r) := by
      exact (FinPMF.expect_map
        (iidLaw (Lsched m) (bsched m) β hβ0 hβ1)
        (scalableFlatInput m)
        (fun xs => ∑ r : State22,
          inputListKernel z (List.ofFn xs) 0 r)).symm
    _ = (piPMF (fun _ : Fin (bsched m * scalableStepsPerRegion m) =>
          bernoulliInputLaw β hβ0 hβ1)).expect
        (fun xs => ∑ r : State22,
          inputListKernel z (List.ofFn xs) 0 r) := by
      rw [scalableFlatInput_iidLaw]
    _ = _ := by
      rw [FinPMF.expect_sum]
      change ∑ r : State22,
          matrixExpect22
            (piPMF (fun _ : Fin (bsched m * scalableStepsPerRegion m) =>
              bernoulliInputLaw β hβ0 hβ1))
            (fun xs => Matrix.of (inputListKernel z (List.ofFn xs))) 0 r = _
      rw [bernoulli_input_kernel_power]

end Spin.UniformMixer
