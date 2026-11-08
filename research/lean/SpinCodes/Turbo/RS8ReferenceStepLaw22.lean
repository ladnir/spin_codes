import SpinCodes.Turbo.RS8StepValueFactor22

/-!
# Physical-step law under the actual Turbo regional route

The normalized reference first chooses a uniform regional selected set,
then independent uniform four-bit values on that set. Its physical-step
words are a hypergeometric mixture of independent `selectedStepLaw22`
inputs. The count profile is shared across the steps; their values are
independent only after conditioning on it.
-/

noncomputable section

namespace Spin.Turbo

open Finset

theorem uniformStepSelectedSet22_bind_values (j : Fin 17) :
    (uniformStepSelectedSet22 j).bind regionalSelectedGiven22 =
      selectedStepLaw22 j := by
  rw [uniformStepSelectedSet22, FinPMF.bind_map]
  rfl

theorem turboReferenceOneRegionLaw22_step_mixture (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (turboReferenceOneRegionLaw22 R S).map
        (regionalStepWordEquiv22 R) =
      (regionalStepProfileLaw22 R S).bind
        (fun counts => Spin.piPMF (fun step : Fin R =>
          selectedStepLaw22 (counts step))) := by
  rw [turboReferenceOneRegionLaw22_eq_selected_mixture,
    FinPMF.map_bind]
  simp_rw [regionalSelectedGiven22_step_product]
  have hsets :
      (Spin.Structured.Routing.shuffleLaw S).bind
        (fun T => Spin.piPMF (fun step : Fin R =>
          regionalSelectedGiven22 (regionalStepSelectedSet22 T step))) =
        (regionalSelectedStepSetsLaw22 R S).bind
          (fun sets => Spin.piPMF (fun step : Fin R =>
            regionalSelectedGiven22 (sets step))) := by
    rw [regionalSelectedStepSetsLaw22, FinPMF.bind_map]
    rfl
  rw [hsets, regionalSelectedStepSetsLaw22_eq_profile_mixture,
    FinPMF.bind_assoc]
  congr 1
  funext counts
  change (Spin.piPMF (fun step : Fin R =>
      uniformStepSelectedSet22 (counts step))).bind
        (fun sets => Spin.piPMF (fun step : Fin R =>
          regionalSelectedGiven22 (sets step))) = _
  rw [Spin.piPMF_bind]
  congr 1
  funext step
  exact uniformStepSelectedSet22_bind_values (counts step)

#print axioms uniformStepSelectedSet22_bind_values
#print axioms turboReferenceOneRegionLaw22_step_mixture

end Spin.Turbo
