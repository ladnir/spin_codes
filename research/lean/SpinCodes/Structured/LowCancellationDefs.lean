import SpinCodes.Structured.ConcreteMapData

/-! Small, executable certificates for the two exceptional input layers.
The generator is untrusted; all fields, coverage, and histograms are checked. -/
namespace Spin.Structured.LowCancellation
open PackedMap

structure Group where
  syndrome : Nat
  shell : Nat
  entries : List (Nat × Nat)
  deriving DecidableEq, Repr

def inputs (gs : List Group) : List Nat :=
  gs.flatMap fun g => g.entries.map Prod.fst

def expand (p : List (Nat × Nat)) : List Nat :=
  p.flatMap fun ec => List.replicate ec.2 ec.1

def valid (j : Nat) (g : Group) : Prop :=
  g.syndrome < 2 ^ 19 ∧ g.syndrome ≠ 0 ∧
  weight 128 (eval ConcreteMaps.aRows g.syndrome) = g.shell ∧
  ∀ r ∈ g.entries, r.1 < 2 ^ 128 ∧ weight 128 r.1 = j ∧
    eval ConcreteMaps.cRows r.1 = g.syndrome ∧
    weight 128 (r.1 ^^^ eval ConcreteMaps.aRows g.syndrome) = r.2

instance (j : Nat) (g : Group) : Decidable (valid j g) := by
  unfold valid
  infer_instance

def shellExponents (gs : List Group) (w : Nat) : List Nat :=
  gs.flatMap fun g => if g.shell = w then g.entries.map Prod.snd else []

end Spin.Structured.LowCancellation
