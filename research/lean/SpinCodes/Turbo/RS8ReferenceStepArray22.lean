import SpinCodes.Turbo.RS8ReferenceStepLaw22

/-!
# Joint physical-step law across the 64 Turbo regions

The reference route has independent regional permutations. Each region
has its own hypergeometric count profile, and conditional on that
profile its ordered complete four-bit step words follow the existing
physical selected-step law. This is a law of inputs, not a state reset:
the inner will consume the resulting region-major stream with one state.
-/

noncomputable section

namespace Spin.Turbo

def regionalStepWords22 (R : ℕ)
    (regions : Fin 64 → Fin (R * 16) → Packet4) :
    Fin 64 → Fin R → PacketStep16 :=
  fun region => regionalStepWordEquiv22 R (regions region)

theorem stepInput_eq_regionalStepWords22 {R : ℕ}
    (regions : Fin 64 → Fin (R * 16) → Packet4)
    (region : Fin 64) (step : Fin R) :
    stepInput regions region step =
      packetStepInput22 ((regionalStepWords22 R regions) region step) := by
  rw [stepInput_eq_packetStepInput22]
  rfl

theorem flatInput_eq_regionalStepWords22 {R : ℕ}
    (regions : Fin 64 → Fin (R * 16) → Packet4)
    (index : Fin (64 * R)) :
    flatInput regions index =
      packetStepInput22
        ((regionalStepWords22 R regions)
          ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index).1
          ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index).2) := by
  change stepInput regions
    ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index).1
    ((finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index).2 = _
  exact stepInput_eq_regionalStepWords22 regions _ _

theorem turboReferenceRoutedLaw22_step_array (R : ℕ)
    (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) :
    (turboReferenceRoutedLaw22 R points S).map
        (regionalStepWords22 R) =
      Spin.piPMF (fun _ : Fin 64 =>
        (regionalStepProfileLaw22 R S).bind
          (fun counts => Spin.piPMF (fun step : Fin R =>
            selectedStepLaw22 (counts step)))) := by
  rw [turboReferenceRoutedLaw22_product]
  change (Spin.piPMF (fun _ : Fin 64 =>
    turboReferenceOneRegionLaw22 R S)).map
      (fun regions region => (regionalStepWordEquiv22 R) (regions region)) = _
  rw [Spin.piPMF_map]
  simp only [turboReferenceOneRegionLaw22_step_mixture]

#print axioms turboReferenceRoutedLaw22_step_array
#print axioms flatInput_eq_regionalStepWords22

end Spin.Turbo
