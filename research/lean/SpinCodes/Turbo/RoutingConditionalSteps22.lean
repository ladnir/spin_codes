import SpinCodes.Turbo.RoutingSelectedFiber22
import SpinCodes.Turbo.RoutingFixedValues22

/-!
# Complete physical-step input law inside actual routing-count fibers

Both the latent selected subsets and all selected four-bit values remain
explicit. A selected slot may be zero. The actual count-conditioned route
therefore induces independent `selectedStepLaw22` inputs; this is not inferred
from an unconditional first moment.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.UniformMixer

abbrev TurboStepArray22 (R : ℕ) := Fin 64 → Fin R → PacketStep16

def turboCountStepArrayLaw22 {R : ℕ} (counts : TurboPotentialCounts22 R) :
    FinPMF (TurboStepArray22 R) :=
  Spin.piPMF (fun region : Fin 64 => Spin.piPMF (fun step : Fin R =>
    selectedStepLaw22 (counts region step)))

/-- The reference array is independent of the regional route fiber. -/
def turboReferenceCountStepLaw22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R) :
    FinPMF (TurboStepArray22 R) :=
  ((turboReferenceArrayLaw22 R points S).prod
    (turboRegionCountFiber22 R S counts)).map
      (fun pair => regionalStepWords22 R (turboRouteRegions22 pair.1 pair.2))

/-- Exact conditional complete-value law in every region and physical step. -/
theorem turboReferenceCountStepLaw22_eq_selected (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (counts : TurboPotentialCounts22 R)
    (hp : 0 < (turboPotentialCountLaw22 R S).p counts) :
    turboReferenceCountStepLaw22 R points S counts = turboCountStepArrayLaw22 counts := by
  let K : TurboSelectedStepSets22 R → FinPMF (TurboStepArray22 R) :=
    fun sets => Spin.piPMF (fun region : Fin 64 => Spin.piPMF (fun step : Fin R =>
      regionalSelectedGiven22 (sets region step)))
  have hswap : turboReferenceCountStepLaw22 R points S counts =
      (turboRegionCountFiber22 R S counts).bind (fun perms =>
        (turboReferenceArrayLaw22 R points S).map
          (fun words => regionalStepWords22 R (turboRouteRegions22 words perms))) := by
    apply FinPMF.eq_of_expect_eq
    intro f
    simp only [turboReferenceCountStepLaw22, FinPMF.expect_map,
      FinPMF.expect_prod, FinPMF.expect_bind]
    unfold FinPMF.expect
    simp only [mul_sum]
    rw [sum_comm]
    apply sum_congr rfl
    intro perms _
    apply sum_congr rfl
    intro words _
    ring
  calc
    _ = (turboRegionCountFiber22 R S counts).bind
        (fun perms => K (turboSelectedStepSets22 S perms)) := by
      rw [hswap]
      apply congrArg (FinPMF.bind (turboRegionCountFiber22 R S counts))
      funext perms
      exact turboReferenceArrayLaw22_fixed_route_steps R points S perms
    _ = ((turboRegionCountFiber22 R S counts).map (turboSelectedStepSets22 S)).bind K :=
      (FinPMF.bind_map _ _ _).symm
    _ = (turboUniformSelectedStepSets22 counts).bind K := by
      rw [turboRegionCountFiber22_selected_sets R S counts hp]
    _ = Spin.piPMF (fun region : Fin 64 =>
        (uniformStepSelectedSets22 (counts region)).bind
          (fun sets => Spin.piPMF (fun step : Fin R => regionalSelectedGiven22 (sets step)))) := by
      change (Spin.piPMF (fun region : Fin 64 =>
        uniformStepSelectedSets22 (counts region))).bind
          (fun sets => Spin.piPMF (fun region : Fin 64 =>
            Spin.piPMF (fun step : Fin R => regionalSelectedGiven22 (sets region step)))) = _
      exact Spin.piPMF_bind
        (fun region : Fin 64 => uniformStepSelectedSets22 (counts region))
        (fun (_ : Fin 64) (sets : Fin R → Finset (Fin 16)) =>
          Spin.piPMF (fun step : Fin R => regionalSelectedGiven22 (sets step)))
    _ = Spin.piPMF (fun region : Fin 64 => Spin.piPMF (fun step : Fin R =>
        (uniformStepSelectedSet22 (counts region step)).bind regionalSelectedGiven22)) := by
      apply congrArg Spin.piPMF
      funext region
      exact Spin.piPMF_bind
        (fun step : Fin R => uniformStepSelectedSet22 (counts region step))
        (fun (_ : Fin R) (T : Finset (Fin 16)) => regionalSelectedGiven22 T)
    _ = _ := by simp only [uniformStepSelectedSet22_bind_values, turboCountStepArrayLaw22]

/-- Region-major chronology, with one state carried through all regions. -/
def turboStepArrayInput22 {R : ℕ} (words : TurboStepArray22 R) :
    Fin (64 * R) → Input64 :=
  fun index =>
    let pos := (finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index
    packetStepInput22 (words pos.1 pos.2)

/-- Conditional reference failure is exactly the complete selected-step
experiment with the original independent GL schedule. -/
theorem turboReferenceFiberFailure22_eq_selected
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R)
    (hp : 0 < (turboPotentialCountLaw22 R S).p counts) (d : ℕ) :
    (turboReferenceArrayLaw22 R points S).expect
      (fun words => turboFiberRemainingFailure22 R S counts words d) =
    ((turboCountStepArrayLaw22 counts).prod (turboStateLaw22 R)).prob
      (fun pair => sampledOutputWeight (turboStepArrayInput22 pair.1) pair.2 0 ≤ d) := by
  rw [← turboReferenceCountStepLaw22_eq_selected R points S counts hp,
    FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod,
    turboReferenceCountStepLaw22, FinPMF.expect_map, FinPMF.expect_prod]
  unfold turboFiberRemainingFailure22
  simp only [FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod]
  apply congrArg
  funext words
  apply congrArg
  funext perms
  apply congrArg
  funext maps
  have he : turboStepArrayInput22
      (regionalStepWords22 R (turboRouteRegions22 words perms)) =
      flatInput (turboRouteRegions22 words perms) := by
    funext index
    exact (flatInput_eq_regionalStepWords22 _ index).symm
  rw [he]

#print axioms turboReferenceCountStepLaw22_eq_selected
#print axioms turboReferenceFiberFailure22_eq_selected
end Spin.Turbo


