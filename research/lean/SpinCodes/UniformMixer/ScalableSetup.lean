import SpinCodes.UniformMixer.RegionalLimit
import SpinCodes.Structured.ConcreteOuterTailClosure
import SpinCodes.Structured.ConcreteRoutePermutation

/-!
# The sampled setup for Scalable SPIN

The outer is the paper's one shared Golay–BAA constituent. Its two
permutations are sampled once. The one-bit route has independent group
and region permutations. Every 64-bit step samples an independent
uniform invertible linear map on the 22-bit state. This module records
the joint finite law and its outer and route marginals.
-/

noncomputable section

namespace Spin.UniformMixer

attribute [local instance] Classical.propDecidable

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute Spin.Structured.ConcreteNativeFamily Filter
open scoped Topology

def scalableStepsPerRegion (m : ℕ) : ℕ := 2 * (m + 1)

theorem scalable_region_length (m : ℕ) :
    Lsched m = 64 * scalableStepsPerRegion m := by
  unfold Lsched scalableStepsPerRegion
  ring

abbrev ScalableRouteSeed (m : ℕ) :=
  Spin.Structured.ConcreteRoute.Seeds (Lsched m) (bsched m)

abbrev StateMap := State22 ≃ₗ[ZMod 2] State22

noncomputable instance : Finite StateMap := state22GL_finite

noncomputable instance : Fintype StateMap :=
  Fintype.ofFinite StateMap

abbrev ScalableStateMaps (m : ℕ) :=
  Fin (bsched m) → Fin (scalableStepsPerRegion m) → StateMap

def scalableStateMapLaw (m : ℕ) : FinPMF (ScalableStateMaps m) :=
  piPMF (fun _ : Fin (bsched m) =>
    piPMF (fun _ : Fin (scalableStepsPerRegion m) =>
      FinPMF.uniform StateMap))

abbrev ScalableSetupSeed (m : ℕ) :=
  NativeSeed m × (ScalableRouteSeed m × ScalableStateMaps m)

def scalableSetupLaw (m : ℕ) : FinPMF (ScalableSetupSeed m) :=
  (nativeSeedLaw m).prod
    ((Spin.Structured.ConcreteRoute.seedLaw (Lsched m) (bsched m)).prod
      (scalableStateMapLaw m))

theorem scalable_outer_marginal (m : ℕ) :
    (scalableSetupLaw m).map (fun seed => seed.1) = nativeSeedLaw m := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [FinPMF.expect_map, scalableSetupLaw, FinPMF.expect_prod]
  simp only [FinPMF.expect_const]

theorem scalable_route_marginal (m : ℕ) :
    (scalableSetupLaw m).map (fun seed => seed.2.1) =
      Spin.Structured.ConcreteRoute.seedLaw (Lsched m) (bsched m) := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [FinPMF.expect_map, scalableSetupLaw, FinPMF.expect_prod]
  simp only [FinPMF.expect_prod, FinPMF.expect_const]

theorem scalable_outer_selection_failure_tendsto :
    Tendsto (fun m : ℕ =>
      (scalableSetupLaw m).prob (fun seed => ¬ nativeGood m seed.1))
      atTop (𝓝 0) := by
  have h := Spin.Structured.ConcreteNativeFamily.raw_selection_failure_tendsto
    Spin.Structured.ConcreteNativeFamily.nativeTail_tendsto
  apply h.congr
  intro m
  rw [← scalable_outer_marginal m, FinPMF.map_prob]

theorem scalable_route_law (m : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))) :
    (scalableSetupLaw m).map
      (fun seed => Spin.Structured.ConcreteRoute.routeEval rows seed.2.1) =
    Spin.Structured.ConcreteRoute.law rows := by
  calc
    _ = ((scalableSetupLaw m).map (fun seed => seed.2.1)).map
          (Spin.Structured.ConcreteRoute.routeEval rows) := by
        rw [FinPMF.map_comp]
        rfl
    _ = _ := by
      rw [scalable_route_marginal]
      exact Spin.Structured.ConcreteRoute.permutation_law_eq rows

theorem scalable_route_expect (m : ℕ)
    (message : Fin (Nsched m / 2) → Bool)
    (F : (Fin (bsched m) → Finset (Fin (Lsched m))) → ℝ) :
    (scalableSetupLaw m).expect (fun seed =>
      F (Spin.Structured.ConcreteRoute.routeEval
        (nativeRows m seed.1 message) seed.2.1)) =
    (nativeSeedLaw m).expect (fun out =>
      (Spin.Structured.ConcreteRoute.law (nativeRows m out message)).expect F) := by
  rw [scalableSetupLaw, FinPMF.expect_prod]
  congr 1
  funext out
  rw [FinPMF.expect_prod]
  simp only [FinPMF.expect_const]
  rw [← FinPMF.expect_map,
    Spin.Structured.ConcreteRoute.permutation_law_eq]

end Spin.UniformMixer
