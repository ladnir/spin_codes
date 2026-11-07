import Mathlib.Algebra.Field.Rat
import SpinCodes.UniformMixer.ScalableDenseTransferEvalDefs
import SpinCodes.Structured.DenseOccupationFixedVertexDefs

/-! Exact rational checker for the production three-state transfer. -/

namespace Spin.UniformMixer

def transferRatG0 (b z : ℚ) : ℚ := 1 - b + b * z
def transferRatG1 (b z : ℚ) : ℚ := b + (1 - b) * z
def transferRatH1 (b z : ℚ) : ℚ := 1 - b - b * z

def transferRatWeight (b z : ℚ) (w : Nat) : ℚ :=
  transferRatG0 b z ^ (64 - w) * transferRatG1 b z ^ w

def transferRatFeedback (b z : ℚ) (w : Nat) : ℚ :=
  transferRatG0 b z ^ (64 - w) * transferRatH1 b z ^ w

def transferRatNeutral (b : ℚ) (w : Nat) : ℚ :=
  (1 - 2 * b) ^ w

def transferRatSpectrum (f : Nat → ℚ) : ℚ :=
  transferSpectrumTerms.foldl
    (fun a term => a + term.2 * f term.1) 0

def transferRatH (b z : ℚ) : ℚ :=
  transferRatSpectrum (transferRatFeedback b z) / 4194304

def transferRatG (b z : ℚ) : ℚ :=
  transferRatG0 b z ^ 64 - transferRatH b z

def transferRatU (b z : ℚ) : ℚ :=
  (transferRatSpectrum (transferRatWeight b z) -
    transferRatWeight b z 0) / 4194303

def transferRatu (b : ℚ) : ℚ :=
  1 - transferRatSpectrum (transferRatNeutral b) / 4194304

def transferRatF (b z : ℚ) : ℚ :=
  (transferLiveWeights.map (transferRatWeight b z)).foldl max 0

def transferRatRows (b z : ℚ) (v : Fin 3 → ℚ) : Fin 3 → ℚ :=
  let H := transferRatH b z
  let G := transferRatG b z
  let F := transferRatF b z
  let U := transferRatU b z
  let u := transferRatu b
  fun i =>
    if i = 0 then H * v 0 + G * v 1
    else if i = 1 then min F u / 4194303 * v 0 + F * v 2
    else min U u / 4194303 * v 0 + U * v 2

def transferRatCheck (b z ρ : ℚ) (v : Fin 3 → ℚ) : Bool :=
  decide (transferRatRows b z v 0 ≤ ρ * v 0) &&
  decide (transferRatRows b z v 1 ≤ ρ * v 1) &&
  decide (transferRatRows b z v 2 ≤ ρ * v 2)

def transferQInputRat
    (a : Spin.Structured.DenseOccupationFixed.QInput) : ℚ :=
  a.num / a.den

end Spin.UniformMixer
