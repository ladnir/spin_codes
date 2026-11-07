import SpinCodes.UniformMixer.Map22Packed

/-!
# Relating packed words to the production linear expansion

The spectrum checker operates on `Nat` bit masks. The encoder proof uses
`ZMod 2` coordinates. This bridge identifies every output coordinate for
every packed input, so a full checked packed histogram will count the
actual production expansion.
-/

namespace Spin.UniformMixer

def bitScalar (b : Bool) : ZMod 2 := if b then 1 else 0

theorem bitScalar_xor (a b : Bool) :
    bitScalar (Bool.xor a b) = bitScalar a + bitScalar b := by
  cases a <;> cases b <;> decide

theorem packed_eval_coordinate (rows : List Nat) (q j : Nat) :
    bitScalar ((Spin.Structured.PackedMap.eval rows q).testBit j) =
      ∑ i : Fin rows.length,
        bitScalar ((rows.getD i.val 0).testBit j) *
          bitScalar (q.testBit i.val) := by
  induction rows generalizing q with
  | nil => simp [Spin.Structured.PackedMap.eval, bitScalar]
  | cons r rs ih =>
      simp only [List.length_cons]
      rw [Fin.sum_univ_succ]
      change bitScalar ((Spin.Structured.PackedMap.eval (r :: rs) q).testBit j) =
        bitScalar (r.testBit j) * bitScalar (q.testBit 0) +
          ∑ i : Fin rs.length,
            bitScalar ((rs.getD i.val 0).testBit j) *
              bitScalar (q.testBit (i.val + 1))
      have htail :
          (∑ i : Fin rs.length,
            bitScalar ((rs.getD i.val 0).testBit j) *
              bitScalar (q.testBit (i.val + 1))) =
            bitScalar ((Spin.Structured.PackedMap.eval rs (q >>> 1)).testBit j) := by
        rw [ih]
        apply Finset.sum_congr rfl
        intro i _
        simp [Nat.testBit_shiftRight, Nat.add_comm]
      rw [htail]
      cases hq : q.testBit 0 with
      | false => simp [Spin.Structured.PackedMap.eval, hq, bitScalar]
      | true =>
          simp only [Spin.Structured.PackedMap.eval, hq, ite_true]
          rw [Nat.testBit_xor, bitScalar_xor]
          simp [bitScalar]

def packedState (q : Nat) : State22 :=
  fun i => bitScalar (q.testBit i.val)

theorem productionPacked_coordinate (q : Nat) (j : Fin 64) :
    bitScalar ((Spin.Structured.PackedMap.eval productionPackedRows q).testBit
      j.val) = coordMap j (packedState q) := by
  rw [packed_eval_coordinate]
  change (∑ i : Fin 22,
    bitScalar ((productionPackedRows.getD i.val 0).testBit j.val) *
      bitScalar (q.testBit i.val)) =
    ∑ i : Fin 22, coordCoeff j i * packedState q i
  apply Finset.sum_congr rfl
  intro i _
  rw [← Spin.Structured.PackedMap.eval_onehot,
    productionPacked_basis]
  rfl

theorem bitScalar_injective : Function.Injective bitScalar := by
  intro a b h
  cases a <;> cases b <;> simp [bitScalar] at h ⊢

theorem packedState_injective :
    Function.Injective (fun q : Fin (2 ^ 22) => packedState q.val) := by
  intro q r h
  apply Fin.ext
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < 22
  · exact bitScalar_injective (congrFun h ⟨i, hi⟩)
  · have hpow : 2 ^ 22 ≤ 2 ^ i := Nat.pow_le_pow_right (by decide) (by omega)
    rw [Nat.testBit_lt_two_pow (q.isLt.trans_le hpow),
      Nat.testBit_lt_two_pow (r.isLt.trans_le hpow)]

theorem packedState_bijective :
    Function.Bijective (fun q : Fin (2 ^ 22) => packedState q.val) := by
  have hcard : Fintype.card (Fin (2 ^ 22)) = Fintype.card State22 := by
    simp [State22, ZMod.card]
  exact (Fintype.bijective_iff_injective_and_card _).mpr
    ⟨packedState_injective, hcard⟩

theorem packedBit_weight (b : Bool) :
    b.toNat = if bitScalar b = 1 then 1 else 0 := by
  cases b <;> decide

theorem productionPacked_weight (q : Nat) :
    Spin.Structured.PackedMap.weight 64
      (Spin.Structured.PackedMap.eval productionPackedRows q) =
        coordinateWeight coordMap (packedState q) := by
  rw [Spin.Structured.PackedMap.weight_eq_sum_bits]
  unfold coordinateWeight
  apply Finset.sum_congr rfl
  intro j _
  rw [packedBit_weight, productionPacked_coordinate]

theorem packedWeight_128_eq_64 (x : Nat) (hx : x < 2 ^ 64) :
    Spin.Structured.PackedMap.weight 128 x =
      Spin.Structured.PackedMap.weight 64 x := by
  rw [Spin.Structured.PackedMap.weight_eq_sum_bits,
    Spin.Structured.PackedMap.weight_eq_sum_bits]
  change (∑ i : Fin 128, (fun j : Nat => (x.testBit j).toNat) i.val) =
    ∑ i : Fin 64, (fun j : Nat => (x.testBit j).toNat) i.val
  rw [Fin.sum_univ_eq_sum_range (fun j : Nat => (x.testBit j).toNat),
    Fin.sum_univ_eq_sum_range (fun j : Nat => (x.testBit j).toNat)]
  conv_lhs =>
    rw [show 128 = 64 + 64 by decide, Finset.sum_range_add]
  have hhigh (i : Nat) : x.testBit (i + 64) = false := by
    apply Nat.testBit_lt_two_pow
    exact hx.trans_le (Nat.pow_le_pow_right (by decide) (by omega))
  have hhigh' (i : Nat) : x.testBit (64 + i) = false := by
    simpa only [Nat.add_comm] using hhigh i
  simp [hhigh']

theorem productionPacked_weight128 (q : Nat) :
    Spin.Structured.PackedMap.weight 128
      (Spin.Structured.PackedMap.eval productionPackedRows q) =
        coordinateWeight coordMap (packedState q) := by
  rw [packedWeight_128_eq_64, productionPacked_weight]
  exact Spin.Structured.PackedMap.eval_lt productionPackedRows
    (by decide) q

theorem productionHistogram_card (w : Nat) (hw : w < 129) :
    (Spin.Structured.MapSpectrum.histogram productionPackedRows
      (List.range (2 ^ 22))).getD w 0 =
        (Finset.univ.filter fun q : State22 =>
          coordinateWeight coordMap q = w).card := by
  classical
  rw [Spin.Structured.MapSpectrum.histogram_range_card productionPackedRows
    (2 ^ 22) hw]
  let e : Fin (2 ^ 22) ≃ State22 :=
    Equiv.ofBijective (fun q => packedState q.val) packedState_bijective
  rw [Finset.card_filter, Finset.card_filter]
  apply Fintype.sum_equiv e
  intro q
  change (if Spin.Structured.PackedMap.weight 128
      (Spin.Structured.PackedMap.eval productionPackedRows q.val) = w
    then 1 else 0) =
      if coordinateWeight coordMap (packedState q.val) = w then 1 else 0
  rw [productionPacked_weight128]

end Spin.UniformMixer
