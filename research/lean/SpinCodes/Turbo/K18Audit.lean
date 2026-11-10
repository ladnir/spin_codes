import SpinCodes.Turbo.Generated.K18.Endpoint

/-! Local audit of the K18 Turbo endpoint (oleans built on Peach from
byte-identical sources). -/
open Spin.Turbo

#check @K18Final.k18_endpoint
#print axioms K18Final.k18_endpoint

example (points : Fin 8 → Field16) (hp : Function.Injective points) :
    (turboSetup22 128 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 128 points).Z 52428 ω) ≤ (2 : ℝ) ^ (-(58.24083 : ℝ)) :=
  K18Final.k18_endpoint points hp

-- N = 128 * 16 * 256 = 524288 output bits, cutoff floor(N/10).
example : 128 * 16 * 256 = 524288 ∧ 524288 / 10 = 52428 := by decide
