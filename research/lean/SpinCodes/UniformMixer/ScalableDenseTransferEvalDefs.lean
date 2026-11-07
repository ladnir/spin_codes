import SpinCodes.Numeric.FixedDefs

/-! Reducible fixed-point evaluator for the actual three-state dense transfer. -/

namespace Spin.UniformMixer

open Spin.Numeric

def transferOne : Fix := Fix.ofInt 1

def transferG0 (β z : Fix) : Fix :=
  Fix.add (Fix.sub transferOne β) (Fix.mul β z)

def transferG1 (β z : Fix) : Fix :=
  Fix.add β (Fix.mul (Fix.sub transferOne β) z)

def transferH1 (β z : Fix) : Fix :=
  Fix.sub (Fix.sub transferOne β) (Fix.mul β z)

def transferWeightMoment (β z : Fix) (w : Nat) : Fix :=
  Fix.mul (Fix.pow (transferG0 β z) (64 - w))
    (Fix.pow (transferG1 β z) w)

def transferFeedbackMoment (β z : Fix) (w : Nat) : Fix :=
  Fix.mul (Fix.pow (transferG0 β z) (64 - w))
    (Fix.pow (transferH1 β z) w)

def transferNeutralMoment (β : Fix) (w : Nat) : Fix :=
  Fix.pow (Fix.sub transferOne (Fix.mul (Fix.ofInt 2) β)) w

def transferSpectrumTerms : List (Nat × Int) :=
  [(0, 1), (16, 2604), (24, 291648), (28, 888832),
    (32, 1828134), (36, 888832), (40, 291648),
    (48, 2604), (64, 1)]

def transferSpectrumPoly (f : Nat → Fix) : Fix :=
  transferSpectrumTerms.foldl
    (fun a term => Fix.add a (Fix.mul (Fix.ofInt term.2) (f term.1)))
    Fix.zero

def transferEvalH (β z : Fix) : Fix :=
  Fix.divInt (transferSpectrumPoly (transferFeedbackMoment β z)) 4194304

def transferEvalG (β z : Fix) : Fix :=
  Fix.sub (Fix.pow (transferG0 β z) 64) (transferEvalH β z)

def transferEvalU (β z : Fix) : Fix :=
  Fix.divInt
    (Fix.sub (transferSpectrumPoly (transferWeightMoment β z))
      (transferWeightMoment β z 0)) 4194303

def transferEvalu (β : Fix) : Fix :=
  Fix.sub transferOne
    (Fix.divInt (transferSpectrumPoly (transferNeutralMoment β)) 4194304)

def transferLiveWeights : List Nat := [16, 24, 28, 32, 36, 40, 48, 64]

def transferEvalF (β z : Fix) : Fix :=
  ⟨0, (transferLiveWeights.map
    (fun w => (transferWeightMoment β z w).hi)).foldl max 0⟩

def transferMin (a b : Fix) : Fix :=
  ⟨min a.lo b.lo, min a.hi b.hi⟩

def transferEvalRows (β z : Fix) (v : Fin 3 → Fix) : Fin 3 → Fix :=
  let H := transferEvalH β z
  let G := transferEvalG β z
  let F := transferEvalF β z
  let U := transferEvalU β z
  let u := transferEvalu β
  fun i =>
    if i = 0 then
      Fix.add (Fix.mul H (v 0)) (Fix.mul G (v 1))
    else if i = 1 then
      Fix.add (Fix.mul (Fix.divInt (transferMin F u) 4194303) (v 0))
        (Fix.mul F (v 2))
    else
      Fix.add (Fix.mul (Fix.divInt (transferMin U u) 4194303) (v 0))
        (Fix.mul U (v 2))

def transferCheck (β z ρ : Fix) (v : Fin 3 → Fix) : Bool :=
  decide ((transferEvalRows β z v 0).hi ≤ (Fix.mul ρ (v 0)).lo) &&
  decide ((transferEvalRows β z v 1).hi ≤ (Fix.mul ρ (v 1)).lo) &&
  decide ((transferEvalRows β z v 2).hi ≤ (Fix.mul ρ (v 2)).lo)

end Spin.UniformMixer
