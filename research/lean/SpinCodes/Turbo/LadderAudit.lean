import SpinCodes.Turbo.Generated.G14.Endpoint
import SpinCodes.Turbo.Generated.G16.Endpoint
import SpinCodes.Turbo.Generated.G18.Endpoint
import SpinCodes.Turbo.Generated.G20.Endpoint

/-! Local audit of the s22 ladder rungs (oleans built on Peach from
byte-identical sources, modulo line endings). -/
open Spin.Turbo

#print axioms G14Final.g14_endpoint
#print axioms G16Final.g16_endpoint
#print axioms G18Final.g18_endpoint
#print axioms G20Final.g20_endpoint

-- K = 2^14: N = 8 * 16 * 256 = 32768, cutoff floor(N/10) = 3276.
example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 8 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 8 points).Z 3276 ω) ≤ (2 : ℝ) ^ (-(64.25152 : ℝ)) :=
  G14Final.g14_endpoint points hp

-- K = 2^16: N = 32 * 16 * 256 = 131072, cutoff floor(N/10) = 13107.
example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 32 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 32 points).Z 13107 ω) ≤ (2 : ℝ) ^ (-(62.34042 : ℝ)) :=
  G16Final.g16_endpoint points hp

example : 8 * 16 * 256 = 32768 ∧ 32768 / 10 = 3276 := by decide
example : 32 * 16 * 256 = 131072 ∧ 131072 / 10 = 13107 := by decide

-- K = 2^18: N = 128 * 16 * 256 = 524288, cutoff floor(N/10) = 52428.
example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 128 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 128 points).Z 52428 ω) ≤ (2 : ℝ) ^ (-(60.36266 : ℝ)) :=
  G18Final.g18_endpoint points hp

-- K = 2^20: N = 512 * 16 * 256 = 2097152, cutoff floor(N/10) = 209715.
example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 512 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 512 points).Z 209715 ω) ≤ (2 : ℝ) ^ (-(58.36785 : ℝ)) :=
  G20Final.g20_endpoint points hp

example : 128 * 16 * 256 = 524288 ∧ 524288 / 10 = 52428 := by decide
example : 512 * 16 * 256 = 2097152 ∧ 2097152 / 10 = 209715 := by decide

