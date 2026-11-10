import SpinCodes.Turbo.RS8ReferenceStepArray22

/-!
# Complete reference values under a fixed regional route

The old `ConcreteRouteTranspose.transpose_density` and the checked Turbo
reference transpose separate fixed wiring from randomness. Here the
regional permutations are fixed. Uniform reference values remain
independent on the selected positions, including zero values. A later
conditioning argument can therefore average only the selected sets.
-/
noncomputable section
namespace Spin.Turbo
open Spin.Structured.Routing

theorem turboReferenceArrayLaw22_fixed_route (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (perms : TurboRegionPerms22 R) :
    (turboReferenceArrayLaw22 R points S).map
        (fun words => turboRouteRegions22 words perms) =
      Spin.piPMF (fun region : Fin 64 =>
        regionalSelectedGiven22 (shuffleSupport (perms region) S)) := by
  calc
    _ = ((turboReferenceArrayLaw22 R points S).map
          (packetTransposeEquiv22 R)).map
          (fun regions region => regionalValueShuffleEquiv22
            (perms region) (regions region)) := by
      rw [FinPMF.map_comp]
      rfl
    _ = (turboReferenceRegionPreLaw22 R S).map
          (fun regions region => regionalValueShuffleEquiv22
            (perms region) (regions region)) := by
      rw [turboReferenceArrayLaw22_transpose]
    _ = Spin.piPMF (fun region : Fin 64 =>
          (Spin.piPMF (fun i : Fin (R * 16) =>
            turboReferenceCellLaw22 S i)).map
            (regionalValueShuffleEquiv22 (perms region))) :=
      Spin.piPMF_map _ _
    _ = _ := by simp only [turboReferenceGiven_shuffle]

theorem turboReferenceArrayLaw22_fixed_route_steps (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (perms : TurboRegionPerms22 R) :
    (turboReferenceArrayLaw22 R points S).map
        (fun words => regionalStepWords22 R (turboRouteRegions22 words perms)) =
      Spin.piPMF (fun region : Fin 64 => Spin.piPMF (fun step : Fin R =>
        regionalSelectedGiven22
          (regionalStepSelectedSet22 (shuffleSupport (perms region) S) step))) := by
  calc
    _ = ((turboReferenceArrayLaw22 R points S).map
          (fun words => turboRouteRegions22 words perms)).map
          (regionalStepWords22 R) := by rw [FinPMF.map_comp]; rfl
    _ = (Spin.piPMF (fun region : Fin 64 =>
          regionalSelectedGiven22 (shuffleSupport (perms region) S))).map
          (regionalStepWords22 R) := by rw [turboReferenceArrayLaw22_fixed_route]
    _ = Spin.piPMF (fun region : Fin 64 =>
          (regionalSelectedGiven22 (shuffleSupport (perms region) S)).map
            (regionalStepWordEquiv22 R)) := Spin.piPMF_map _ _
    _ = _ := by simp only [regionalSelectedGiven22_step_product]

#print axioms turboReferenceArrayLaw22_fixed_route
#print axioms turboReferenceArrayLaw22_fixed_route_steps
end Spin.Turbo
