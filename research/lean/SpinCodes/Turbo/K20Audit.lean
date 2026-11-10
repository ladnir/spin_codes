import SpinCodes.Turbo.Generated.K20.Endpoint
import SpinCodes.Turbo.Generated.K18.Endpoint

/-! Local audit of both Turbo endpoints (oleans built on Peach from
byte-identical sources, modulo line endings). -/
open Spin.Turbo

#check @K20Final.k20_endpoint
#print axioms K20Final.k20_endpoint
#print axioms K18Final.k18_endpoint

example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 512 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 512 points).Z 209715 ω) ≤ (2 : ℝ) ^ (-(56.25709 : ℝ)) :=
  K20Final.k20_endpoint points hp

example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 128 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 128 points).Z 52428 ω) ≤ (2 : ℝ) ^ (-(58.24083 : ℝ)) :=
  K18Final.k18_endpoint points hp

-- K20: N = 512 * 16 * 256 = 2097152, cutoff floor(N/10) = 209715.
example : 512 * 16 * 256 = 2097152 ∧ 2097152 / 10 = 209715 := by decide
