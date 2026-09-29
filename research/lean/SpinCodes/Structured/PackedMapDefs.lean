/-! Executable binary linear maps. Row zero corresponds to the least
significant input bit, as in the paper's hexadecimal map table. -/

namespace Spin.Structured.PackedMap

def eval : List Nat → Nat → Nat
  | [], _ => 0
  | r :: rs, q =>
      let tail := eval rs (q >>> 1)
      if q.testBit 0 then r ^^^ tail else tail

def weight : Nat → Nat → Nat
  | 0, _ => 0
  | n + 1, q => (q.testBit 0).toNat + weight n (q >>> 1)

def pack : List Bool → Nat
  | [] => 0
  | b :: bs => b.toNat + 2 * pack bs

def transpose (width : Nat) (rows : List Nat) : List Nat :=
  (List.range width).map fun i => pack (rows.map fun r => r.testBit i)

def identityRows : Nat → List Nat
  | 0 => []
  | n + 1 => 1 :: (identityRows n).map (fun r => r <<< 1)

end Spin.Structured.PackedMap
