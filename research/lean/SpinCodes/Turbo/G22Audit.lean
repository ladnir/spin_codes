import SpinCodes.Turbo.Generated.G22.Endpoint

/-! Local audit of the K = 2^22 s22 rung (oleans built on Peach from byte-identical
sources, modulo line endings). Separate from `LadderAudit` because importing every
rung at once exceeds this machine's memory. -/
open Spin.Turbo

#print axioms G22Final.g22_endpoint

-- K = 2^22: N = 2048 * 16 * 256 = 8388608, cutoff floor(N/10) = 838860.
example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 2048 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 2048 points).Z 838860 ω) ≤ (2 : ℝ) ^ (-(56.3691 : ℝ)) :=
  G22Final.g22_endpoint points hp

example : 2048 * 16 * 256 = 8388608 ∧ 8388608 / 10 = 838860 := by decide
