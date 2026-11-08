import SpinCodes.LinearCircuit
import SpinCodes.UniformMixer.ScalableInputEmbedding

/-! The actual sampled one-bit route and chronological reshape are wires. -/
noncomputable section
namespace Spin.UniformMixer
open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteBinaryEncoding Spin.Structured.ConcreteRoute
open Spin.LinearCircuit

/-- Row shuffles, array transposition, then regional shuffles. -/
def routeWiringEquiv {L b : ℕ} (seed : ConcreteRoute.Seeds L b) :
    (Fin L × Fin b) ≃ (Fin b × Fin L) where
  toFun p :=
    let j := (seed.1 p.1).symm p.2
    (j, (seed.2 j).symm p.1)
  invFun p :=
    let i := seed.2 p.1 p.2
    (i, seed.1 i p.1)
  left_inv p := by simp
  right_inv p := by simp

def regionBlockEquiv (R : ℕ) : (Fin R × Fin 64) ≃ Fin (64 * R) where
  toFun p := regionBit p.1 p.2
  invFun p := (regionStep p, regionLocal p)
  left_inv p := by simp only [regionBit_step, regionBit_local]
  right_inv := regionBit_reconstruct

/-- Reshape routed regional bits into the complete chronological step stream. -/
def scalableRouteReshape (m : ℕ) :
    (Fin (bsched m) × Fin (Lsched m)) ≃
      (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) :=
  (Equiv.prodCongr (Equiv.refl _) ((finCongr (scalable_region_length m)).trans
    (regionBlockEquiv (scalableStepsPerRegion m)).symm)).trans
    ((Equiv.prodAssoc _ _ _).symm.trans
      (Equiv.prodCongr finProdFinEquiv (Equiv.refl _)))

def scalableRouteCircuit (m : ℕ) (seed : ScalableRouteSeed m) :
    Circuit (Fin (Lsched m) × Fin (bsched m))
      (Fin (bsched m * scalableStepsPerRegion m) × Fin 64) :=
  .reindex ((routeWiringEquiv seed).trans (scalableRouteReshape m))

@[simp] theorem scalableRouteCircuit_size (m : ℕ) (seed : ScalableRouteSeed m) :
    (scalableRouteCircuit m seed).size = 0 := rfl

/-- This circuit uses the actual support convention and sampled pullbacks. -/
theorem scalableRouteCircuit_eval (m : ℕ) (seed : ScalableRouteSeed m)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))) :
    (scalableRouteCircuit m seed).eval
      (fun p => bitEquiv (decide (p.2 ∈ rows p.1))) =
      fun p => scalableFlatInput m (routeEval rows seed) p.1 p.2 := by
  classical
  funext p
  simp only [scalableRouteCircuit, LinearCircuit.Circuit.eval,
    Equiv.trans_apply, scalableRouteReshape, Equiv.prodCongr_symm,
    Equiv.prodCongr_apply, Equiv.refl_apply, Equiv.symm_symm,
    routeWiringEquiv, Equiv.coe_fn_mk, Equiv.prodAssoc_apply,
    scalableFlatInput, scalableRegionInput, regionInput, scalableRegionSupport]
  simp only [Finset.mem_map_equiv, routeEval, Routing.mem_shuffleSupport,
    mem_transpose]
  simp [regionBlockEquiv, bitEquiv]

#print axioms scalableRouteCircuit_eval
end Spin.UniformMixer
