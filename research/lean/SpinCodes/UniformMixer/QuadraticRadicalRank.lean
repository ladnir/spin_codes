import SpinCodes.UniformMixer.QuadraticRadicalPacked

/-!
# Fast six-row radical calculation

This is a compact finite evaluator for the rank distribution of the
quadratic polar forms. It is not yet used as a theorem about the
production map: the row-mask and elimination correctness bridges to
`quadraticRadical` are separate obligations.
-/

namespace Spin.UniformMixer

def polarRowContribution (q : Nat) (i : Fin 6) (p : Fin 15) : Nat :=
  if q.testBit p.val then
    Nat.xor
      (if pairLeft p = i then 2 ^ (pairRight p).val else 0)
      (if pairRight p = i then 2 ^ (pairLeft p).val else 0)
  else 0

def polarRowMask (q : Nat) (i : Fin 6) : Nat :=
  (List.ofFn (fun p : Fin 15 => polarRowContribution q i p)).foldl
    Nat.xor 0

def polarRows (q : Nat) : List Nat :=
  List.ofFn (fun i : Fin 6 => polarRowMask q i)

def eliminatePolarColumn (bit : Nat) (rows : List Nat) : List Nat :=
  match rows.find? (fun row => row.testBit bit) with
  | none => rows
  | some pivot =>
      (rows.erase pivot).map (fun row =>
        if row.testBit bit then Nat.xor row pivot else row)

def polarRank (q : Nat) : Nat :=
  6 - ((List.range 6).reverse.foldl
    (fun rows bit => eliminatePolarColumn bit rows) (polarRows q)).length

def fastRadicalSize (q : Nat) : Nat := 2 ^ (6 - polarRank q)

def fastRadicalBlockCount (start length radicalSize : Nat) : Nat :=
  ((List.range length).filter (fun j =>
    fastRadicalSize (start + j) = radicalSize)).length

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem first256_fast_radical_count :
    fastRadicalBlockCount 0 256 16 = 59 := by
  decide +kernel

end Spin.UniformMixer
