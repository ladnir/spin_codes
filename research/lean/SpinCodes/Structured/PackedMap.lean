import SpinCodes.Structured.PackedMapDefs
import Mathlib.Data.Nat.Bitwise
import Mathlib.Tactic

/-! Soundness of the packed representation used to check the concrete maps.
The evaluator preserves XOR, and checked inverse basis images therefore
establish injectivity and surjectivity on the whole binary vector spaces. -/

namespace Spin.Structured.PackedMap

theorem eval_zero (rs : List ℕ) : eval rs 0 = 0 := by
  induction rs with
  | nil => rfl
  | cons r rs ih => simp [eval, ih]

theorem eval_onehot (rs : List ℕ) (i : ℕ) : eval rs (2 ^ i) = rs.getD i 0 := by
  induction rs generalizing i with
  | nil => simp [eval]
  | cons r rs ih =>
    cases i with
    | zero => simp [eval, eval_zero]
    | succ i =>
      simp [eval, Nat.testBit_two_pow, Nat.shiftRight_eq_div_pow, Nat.pow_succ, ih]

theorem weight_zero (n : ℕ) : weight n 0 = 0 := by
  induction n with
  | zero => rfl
  | succ n ih => simp [weight, ih]

theorem weight_le (n q : ℕ) : weight n q ≤ n := by
  induction n generalizing q with
  | zero => rfl
  | succ n ih =>
    simp only [weight]
    have h := ih (q >>> 1)
    cases q.testBit 0 <;> simp <;> omega

theorem eval_xor (rs : List ℕ) (q p : ℕ) :
    eval rs (q ^^^ p) = eval rs q ^^^ eval rs p := by
  induction rs generalizing q p with
  | nil => simp [eval]
  | cons r rs ih =>
    simp only [eval, Nat.testBit_xor, Nat.shiftRight_xor_distrib, ih]
    cases hq : q.testBit 0 <;> cases hp : p.testBit 0 <;>
      simp [Nat.xor_assoc, Nat.xor_left_comm, Nat.xor_comm]

theorem eval_lt {t : ℕ} (rs : List ℕ) (hr : ∀ r ∈ rs, r < 2 ^ t) (q : ℕ) :
    eval rs q < 2 ^ t := by
  induction rs generalizing q with
  | nil => simpa [eval] using (Nat.two_pow_pos t)
  | cons r rs ih =>
    have ht := ih (fun a ha => hr a (by simp [ha])) (q >>> 1)
    simp only [eval]
    split
    · exact Nat.xor_lt_two_pow (hr r (by simp)) ht
    · exact ht

theorem eval_comp (rs ss : List ℕ) (q : ℕ) :
    eval ss (eval rs q) = eval (rs.map (eval ss)) q := by
  induction rs generalizing q with
  | nil => simp [eval, eval_zero]
  | cons r rs ih =>
    simp only [eval, List.map_cons]
    split <;> simp [eval_xor, ih]

theorem eval_congr_bits (rs : List ℕ) (q p : ℕ)
    (h : ∀ i < rs.length, q.testBit i = p.testBit i) : eval rs q = eval rs p := by
  induction rs generalizing q p with
  | nil => rfl
  | cons r rs ih =>
    have h0 := h 0 (by simp)
    have ht := ih (q >>> 1) (p >>> 1) (by
      intro i hi
      simpa only [Nat.testBit_shiftRight] using h (1 + i) (by simp only [List.length_cons]; omega))
    simp only [eval, h0, ht]

theorem eval_mod (rs : List ℕ) (q : ℕ) : eval rs (q % 2 ^ rs.length) = eval rs q := by
  apply eval_congr_bits
  intro i hi
  simp [Nat.testBit_mod_two_pow, hi]

theorem eval_append (rs ss : List ℕ) (q : ℕ) :
    eval (rs ++ ss) q = eval rs q ^^^ eval ss (q >>> rs.length) := by
  induction rs generalizing q with
  | nil => simp [eval]
  | cons r rs ih =>
    simp only [List.cons_append, eval, List.length_cons, ih]
    split <;> simp [← Nat.shiftRight_add, Nat.add_comm, Nat.xor_assoc]

theorem eval_split (rs : List ℕ) {k low : ℕ} (hk : k ≤ rs.length)
    (hlow : low < 2 ^ k) (high : ℕ) :
    eval rs (high * 2 ^ k + low) =
      eval (rs.take k) low ^^^ eval (rs.drop k) high := by
  have hlen : (rs.take k).length = k := by simp [List.length_take, hk]
  have hlow' : (high * 2 ^ k + low) % 2 ^ k = low := by
    simp [Nat.add_mod, Nat.mod_eq_of_lt hlow]
  have hhigh : (high * 2 ^ k + low) >>> k = high := by
    rw [Nat.shiftRight_eq_div_pow]
    simp [Nat.add_div, Nat.div_eq_of_lt hlow, Nat.mod_eq_of_lt hlow, hlow]
  conv_lhs => rw [← List.take_append_drop k rs, eval_append, hlen, hhigh,
    ← eval_mod (rs.take k), hlen, hlow']

theorem eval_shiftLeft (rs : List ℕ) (q n : ℕ) :
    eval (rs.map fun r => r <<< n) q = (eval rs q) <<< n := by
  induction rs generalizing q with
  | nil => simp [eval]
  | cons r rs ih =>
    simp only [List.map_cons, eval]
    split <;> simp [ih, Nat.shiftLeft_xor_distrib]

theorem eval_identity_testBit (n q i : ℕ) :
    (eval (identityRows n) q).testBit i = (decide (i < n) && q.testBit i) := by
  induction n generalizing q i with
  | zero => simp [identityRows, eval]
  | succ n ih =>
    simp only [identityRows, eval, eval_shiftLeft]
    cases i with
    | zero =>
      cases h : q.testBit 0 <;>
        simp only [Bool.false_eq_true, ite_false, ite_true,
          Nat.testBit_xor, Nat.testBit_shiftLeft, Nat.testBit_one_zero] <;> simp
    | succ i =>
      have hone : (1 : ℕ).testBit (i + 1) = false := by
        apply Nat.testBit_two_pow_of_ne (n := 0)
        omega
      cases h : q.testBit 0 <;>
        simp only [Bool.false_eq_true, ite_false, ite_true,
          Nat.testBit_xor, Nat.testBit_shiftLeft, hone] <;>
        simp [ih, Nat.testBit_shiftRight, Nat.add_comm]

theorem eval_identity (n q : ℕ) : eval (identityRows n) q = q % 2 ^ n := by
  apply Nat.eq_of_testBit_eq
  intro i
  simp [eval_identity_testBit, Nat.testBit_mod_two_pow, Bool.and_comm]

theorem eval_inverse {s : ℕ} (rs ss : List ℕ)
    (h : rs.map (eval ss) = identityRows s) {q : ℕ} (hq : q < 2 ^ s) :
    eval ss (eval rs q) = q := by
  rw [eval_comp, h, eval_identity, Nat.mod_eq_of_lt hq]

theorem eval_injective {s : ℕ} (rs ss : List ℕ)
    (h : rs.map (eval ss) = identityRows s) :
    Function.Injective (fun q : Fin (2 ^ s) => eval rs q) := by
  intro q p he
  apply Fin.ext
  have hh := congrArg (eval ss) he
  simpa [eval_inverse rs ss h q.isLt, eval_inverse rs ss h p.isLt] using hh

theorem weight_eq_sum_bits (n q : ℕ) :
    weight n q = ∑ i : Fin n, (q.testBit i).toNat := by
  induction n generalizing q with
  | zero => simp [weight]
  | succ n ih =>
    simp [weight, Fin.sum_univ_succ, ih, Nat.testBit_shiftRight, Nat.add_comm]

open Finset
open scoped symmDiff

def support (n q : ℕ) : Finset (Fin n) := univ.filter fun i => q.testBit i

theorem weight_eq_card_support (n q : ℕ) : weight n q = (support n q).card := by
  rw [weight_eq_sum_bits, support, Finset.card_filter]
  apply Finset.sum_congr rfl
  intro i _
  cases q.testBit i <;> rfl

theorem support_zero (n : ℕ) : support n 0 = ∅ := by
  simp [support]

theorem support_xor (n q p : ℕ) : support n (q ^^^ p) = support n q ∆ support n p := by
  ext i
  simp only [support, Finset.mem_filter, Finset.mem_univ, true_and,
    Finset.mem_symmDiff, Nat.testBit_xor]
  cases q.testBit i <;> cases p.testBit i <;> simp

theorem support_injective (n : ℕ) :
    Function.Injective (fun q : Fin (2 ^ n) => support n q) := by
  intro q p h
  change support n q = support n p at h
  apply Fin.ext
  apply Nat.eq_of_testBit_eq
  intro i
  by_cases hi : i < n
  · have hh : (⟨i, hi⟩ : Fin n) ∈ support n q ↔ ⟨i, hi⟩ ∈ support n p := by rw [h]
    simpa only [support, Finset.mem_filter, Finset.mem_univ, true_and,
      Bool.coe_iff_coe] using hh
  · have hpow : 2 ^ n ≤ 2 ^ i := Nat.pow_le_pow_right (by decide) (by omega)
    rw [Nat.testBit_lt_two_pow (q.isLt.trans_le hpow),
      Nat.testBit_lt_two_pow (p.isLt.trans_le hpow)]

noncomputable def supportEquiv (n : ℕ) : Fin (2 ^ n) ≃ Finset (Fin n) :=
  Equiv.ofBijective (fun q => support n q)
    ((Fintype.bijective_iff_injective_and_card _).mpr
      ⟨support_injective n, by simp⟩)

end Spin.Structured.PackedMap
