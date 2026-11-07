import SpinCodes.Turbo.PacketRoute
import SpinCodes.UniformMixer.ScalableEncoderKernel

/-!
# Turbo packet stream with the shared 22-state inner

Each of 64 regions contains `R` physical 64-bit steps. Each step reads
16 consecutive four-bit packets. The same sampled state persists over the
entire region-major stream; there is no per-region reset or final flush.
The packet route and the state maps are independent in this ideal law.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured Spin.Structured.ConcreteRoute Spin.UniformMixer

abbrev Packet4 := Fin 4 → ZMod 2

/-- The 16 packet slots belonging to one physical 64-bit step. -/
def stepInput {R : ℕ}
    (packets : Fin 64 → Fin (R * 16) → Packet4)
    (region : Fin 64) (step : Fin R) : Input64 :=
  fun bit =>
    let slot : Fin 16 × Fin 4 :=
      (finProdFinEquiv : Fin 16 × Fin 4 ≃ Fin 64).symm bit
    packets region (finProdFinEquiv (step, slot.1)) slot.2

/-- Region-major serialization of all physical inner inputs. -/
def flatInput {R : ℕ}
    (packets : Fin 64 → Fin (R * 16) → Packet4) :
    Fin (64 * R) → Input64 :=
  fun index =>
    let position : Fin 64 × Fin R :=
      (finProdFinEquiv : Fin 64 × Fin R ≃ Fin (64 * R)).symm index
    stepInput packets position.1 position.2

abbrev RouteInnerSeed22 (R : ℕ) :=
  Seeds (R * 16) 64 × (Fin (64 * R) → StateMap)

/-- Uniform packet permutations and independent uniform GL(22,2) updates. -/
def routeInnerLaw22 (R : ℕ) : FinPMF (RouteInnerSeed22 R) :=
  (seedLaw (R * 16) 64).prod
    (piPMF (fun _ : Fin (64 * R) => FinPMF.uniform StateMap))

/-- Output weight for one fixed routed outer packet array and setup. -/
def outputWeight22 {R : ℕ}
    (rows : Fin (R * 16) → Fin 64 → Packet4)
    (seed : RouteInnerSeed22 R) : ℕ :=
  sampledOutputWeight (flatInput (routePackets rows seed.1)) seed.2 0

/-- The kernel is the exact moment of the sampled route and inner encoder. -/
theorem routeInnerMoment22 {R : ℕ}
    (rows : Fin (R * 16) → Fin 64 → Packet4) (z : ℝ) :
    (routeInnerLaw22 R).expect (fun seed => z ^ outputWeight22 rows seed) =
      (seedLaw (R * 16) 64).expect (fun route =>
        ∑ q : State22,
          inputListKernel z (List.ofFn (flatInput (routePackets rows route)))
            0 q) := by
  rw [routeInnerLaw22, FinPMF.expect_prod]
  apply congrArg
  funext route
  exact sampledInputMoment_kernel z
    (flatInput (routePackets rows route)) 0

#print axioms routeInnerMoment22

end Spin.Turbo
