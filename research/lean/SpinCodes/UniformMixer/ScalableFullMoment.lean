import SpinCodes.UniformMixer.ScalableMapLaw

/-!
# Exact first moment under the Scalable SPIN setup

For a fixed routed outer word, integrating the sampled maps gives the
chronological regional kernel. The second theorem integrates over the
actual outer and route laws as well. No independent code is resampled
per message.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily
  Spin.Structured.ConcreteOuter Spin.Structured.ConcreteRoute

theorem scalable_fixed_route_moment (m : ℕ) (z : ℝ)
    (regions : Fin (bsched m) → Finset (Fin (Lsched m))) :
    (scalableStateMapLaw m).expect (fun maps =>
      z ^ sampledOutputWeight (scalableFlatInput m regions)
        (scalableFlatMaps m maps) 0) =
    ∑ r : State22,
      inputListKernel z (List.ofFn (scalableFlatInput m regions)) 0 r := by
  calc
    _ = (piPMF (fun _ : Fin (bsched m * scalableStepsPerRegion m) =>
          FinPMF.uniform StateMap)).expect (fun maps =>
          z ^ sampledOutputWeight (scalableFlatInput m regions) maps 0) :=
        scalableFlatMaps_expect m _
    _ = _ := sampledInputMoment_kernel z (scalableFlatInput m regions) 0

theorem scalable_setup_moment (m : ℕ) (z : ℝ)
    (message : Fin (Nsched m / 2) → Bool) :
    (scalableSetupLaw m).expect (fun seed =>
      z ^ scalableOutputWeight m message seed) =
    (nativeSeedLaw m).expect (fun out =>
      (Spin.Structured.ConcreteRoute.law (nativeRows m out message)).expect
        (fun regions =>
          ∑ r : State22,
            inputListKernel z
              (List.ofFn (scalableFlatInput m regions)) 0 r)) := by
  rw [scalableSetupLaw, FinPMF.expect_prod]
  congr 1
  funext out
  rw [FinPMF.expect_prod]
  rw [← Spin.Structured.ConcreteRoute.permutation_law_eq
    (nativeRows m out message), FinPMF.expect_map]
  congr 1
  funext route
  simp only [scalableOutputWeight]
  rw [scalable_fixed_route_moment]

end Spin.UniformMixer
