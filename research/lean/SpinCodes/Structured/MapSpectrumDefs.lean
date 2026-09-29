import SpinCodes.Structured.PackedMapDefs

/-! Executable spectra of the concrete packed maps, split into consecutive
input blocks so every finite check can be replayed independently. -/

namespace Spin.Structured.MapSpectrum

def fastWeight : Nat → Nat → Nat
  | 0, _ => 0
  | n + 1, q => q % 2 + fastWeight n (q / 2)

theorem fastWeight_eq (n q : Nat) : fastWeight n q = PackedMap.weight n q := by
  induction n generalizing q with
  | zero => rfl
  | succ n ih =>
    simp only [fastWeight, PackedMap.weight, ih, Nat.shiftRight_eq_div_pow]
    have h := Nat.mod_two_eq_zero_or_one q
    rcases h with h | h <;> simp [Nat.testBit_zero, h]

def bump : Nat → List Nat → List Nat
  | _, [] => []
  | 0, a :: as => (a + 1) :: as
  | i + 1, a :: as => a :: bump i as

def histogram (rows : List Nat) (inputs : List Nat) : List Nat :=
  inputs.foldl (fun h q => bump (PackedMap.weight 128 (PackedMap.eval rows q)) h)
    (List.replicate 129 0)

def block (rows : List Nat) (start size : Nat) : List Nat :=
  histogram rows (List.range' start size)

theorem block_of_weights (rows : List Nat) (start size : Nat) (ws : List Nat)
    (h : (List.range' start size).map (fun q => PackedMap.weight 128 (PackedMap.eval rows q)) = ws) :
    block rows start size = ws.foldl (fun h w => bump w h) (List.replicate 129 0) := by
  have hh := congrArg (fun xs : List Nat => xs.foldl (fun h w => bump w h)
    (List.replicate 129 0)) h
  simpa only [List.foldl_map, block, histogram] using hh

end Spin.Structured.MapSpectrum
