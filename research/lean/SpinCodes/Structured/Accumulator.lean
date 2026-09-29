/-
The exact accumulator weight transition (`eq:structured-acc-transition`).

    T_b(a,c) = C(c-1, r_a - 1) · C(b-c, a - r_a),   r_a = ⌈a/2⌉

is claimed to count the binary words of length `b` and weight `a` whose
accumulation (running parity) has weight `c`; `P_b(a,c) = T_b(a,c)/C(b,a)` is
then the weight-transition probability of a uniformly interleaved accumulator.

The proof here avoids run-length combinatorics.  Splitting on the *last*
letter gives

    N_{b+1}(a, c + a mod 2) = N_b(a, c) + N_b(a-1, c),

because appending a letter extends the accumulation by the parity of the whole
word.  The closed form satisfies the same recursion, in both parities, by a
single application of Pascal's rule.

One trap worth recording: the naive translation of `T_b` to `ℕ` is *wrong* at
`c = 0`.  Truncated subtraction turns `C(-1, r-1)` into `C(0, r-1)`, which is
`1` when `r = 1`, whereas the paper's convention ("binomial coefficients
outside their natural range are zero") makes it `0`.  `accT` guards `c = 0`
explicitly.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open List

/-! ## Words, weight, parity, accumulated weight -/

/-- Number of ones. -/
def wtL : List Bool → ℕ
  | [] => 0
  | x :: t => (if x then 1 else 0) + wtL t

/-- Parity of the number of ones. -/
def parityL : List Bool → Bool
  | [] => false
  | x :: t => xor x (parityL t)

/-- Weight of the accumulation, from a given entering parity. -/
def accWtAux : Bool → List Bool → ℕ
  | _, [] => 0
  | p, x :: t => (if xor p x then 1 else 0) + accWtAux (xor p x) t

/-- Weight of the accumulation of a word. -/
def accWtL (v : List Bool) : ℕ := accWtAux false v

@[simp] lemma wtL_nil : wtL [] = 0 := rfl
@[simp] lemma accWtL_nil : accWtL [] = 0 := rfl

lemma wtL_le_length (v : List Bool) : wtL v ≤ v.length := by
  induction v with
  | nil => simp [wtL]
  | cons x t ih => cases x <;> simp [wtL] <;> omega

lemma wtL_append (v : List Bool) (x : Bool) :
    wtL (v ++ [x]) = wtL v + (if x then 1 else 0) := by
  induction v with
  | nil => simp [wtL]
  | cons y t ih => simp [wtL, ih]; omega

lemma parityL_append (v : List Bool) (x : Bool) :
    parityL (v ++ [x]) = xor (parityL v) x := by
  induction v with
  | nil => simp [parityL]
  | cons y t ih => simp [parityL, ih, Bool.xor_assoc]

/-- Appending a letter extends the accumulation by the parity of the whole
resulting word. -/
lemma accWtAux_append (p : Bool) (v : List Bool) (x : Bool) :
    accWtAux p (v ++ [x]) = accWtAux p v + (if xor (xor p (parityL v)) x then 1 else 0) := by
  induction v generalizing p with
  | nil => simp [accWtAux, parityL]
  | cons y t ih =>
      have hassoc : ∀ q r s w : Bool,
          xor (xor (xor q r) s) w = xor (xor q (xor r s)) w := by decide
      simp only [cons_append, accWtAux, parityL, ih (xor p y),
        hassoc p y (parityL t) x]
      exact (Nat.add_assoc _ _ _).symm

lemma accWtL_append (v : List Bool) (x : Bool) :
    accWtL (v ++ [x]) = accWtL v + (if parityL (v ++ [x]) then 1 else 0) := by
  unfold accWtL
  rw [accWtAux_append, parityL_append]
  simp

/-- The parity bit is the weight mod 2.  This is what turns the accumulation
increment into `a mod 2`, where `a` is the weight of the whole word. -/
lemma parityL_toNat (v : List Bool) : (if parityL v then 1 else 0) = wtL v % 2 := by
  induction v with
  | nil => simp [parityL, wtL]
  | cons x t ih =>
      cases x <;> cases hp : parityL t <;> simp_all [parityL, wtL] <;> omega

/-! ## The accumulation as a map

Composing two stages needs the accumulated *word*, not only its weight: the
second permutation acts on the output of the first accumulator. -/

/-- The accumulated word, from a given entering parity. -/
def accLAux : Bool -> List Bool -> List Bool
  | _, [] => []
  | p, x :: t => (xor p x) :: accLAux (xor p x) t

/-- The accumulated word. -/
def accL (v : List Bool) : List Bool := accLAux false v

@[simp] lemma length_accLAux (p : Bool) (v : List Bool) :
    (accLAux p v).length = v.length := by
  induction v generalizing p with
  | nil => rfl
  | cons x t ih => simp [accLAux, ih]

@[simp] lemma length_accL (v : List Bool) : (accL v).length = v.length :=
  length_accLAux false v

/-- The weight of the accumulated word is the accumulated weight: the two
definitions agree term by term. -/
lemma wtL_accLAux (p : Bool) (v : List Bool) : wtL (accLAux p v) = accWtAux p v := by
  induction v generalizing p with
  | nil => rfl
  | cons x t ih => simp [accLAux, wtL, accWtAux, ih]

lemma wtL_accL (v : List Bool) : wtL (accL v) = accWtL v := wtL_accLAux false v

/-! ## Flipping the first bit complements the accumulation

`Acc(V) + 1 = Acc(V + e₁)`: adding the first standard basis vector flips every
accumulated bit, because it flips the entering parity for every position.  This
is what turns the *upper* tail of the accumulator into a *lower* tail, which is
how `lem:structured-ba-tails` reuses `eq:structured-ba-one-acc-tail` for the
sparse upper tail. -/

/-- Add `e₁`: flip the first bit. -/
def flipHead : List Bool → List Bool
  | [] => []
  | x :: t => (!x) :: t

@[simp] lemma length_flipHead (v : List Bool) : (flipHead v).length = v.length := by
  cases v <;> simp [flipHead]

/-- Complementing the entering parity complements the whole accumulation. -/
lemma accLAux_not (p : Bool) (v : List Bool) :
    accLAux (!p) v = (accLAux p v).map not := by
  induction v generalizing p with
  | nil => simp [accLAux]
  | cons x t ih =>
      have hx : xor (!p) x = !(xor p x) := by cases p <;> cases x <;> rfl
      simp only [accLAux, hx, ih (xor p x), List.map_cons]

/-- `Acc(V + e₁) = ¬ Acc(V)`. -/
lemma accL_flipHead (v : List Bool) : accL (flipHead v) = (accL v).map not := by
  cases v with
  | nil => simp [accL, accLAux, flipHead]
  | cons x t =>
      simp only [flipHead, accL, accLAux, Bool.false_xor, List.map_cons]
      have : accLAux (!x) t = (accLAux x t).map not := by
        have := accLAux_not x t
        exact this
      rw [this]

lemma wtL_map_not (l : List Bool) : wtL (l.map not) = l.length - wtL l := by
  induction l with
  | nil => simp [wtL]
  | cons x t ih =>
      have hle : wtL t ≤ t.length := wtL_le_length t
      cases x <;> simp [wtL, ih] <;> omega

/-- **The upper tail is a lower tail.**  The accumulated weight of `V + e₁` is
the complement of that of `V`. -/
theorem accWtL_flipHead (v : List Bool) :
    accWtL (flipHead v) = v.length - accWtL v := by
  rw [← wtL_accL, ← wtL_accL, accL_flipHead, wtL_map_not, length_accL]

/-- Flipping the first bit moves the input weight by one, in the direction set
by that bit. -/
lemma wtL_flipHead_cons (x : Bool) (t : List Bool) :
    wtL (flipHead (x :: t)) = if x then wtL t else wtL t + 1 := by
  cases x <;> simp [flipHead, wtL] <;> omega

/-- The tail-swap in the form the upper-tail estimate uses. -/
theorem accWtL_ge_iff (v : List Bool) (D : ℕ) (hD : D ≤ v.length) :
    v.length - D ≤ accWtL v ↔ accWtL (flipHead v) ≤ D := by
  have hle : accWtL v ≤ v.length := by
    rw [← wtL_accL]
    simpa using wtL_le_length (accL v)
  rw [accWtL_flipHead]
  omega

/-! ## The word list -/

/-- All binary words of a given length, built by appending on the right. -/
def allWords : ℕ → List (List Bool)
  | 0 => [[]]
  | n + 1 => (allWords n).flatMap (fun v => [v ++ [false], v ++ [true]])

lemma mem_allWords_length {b : ℕ} {v : List Bool} (h : v ∈ allWords b) : v.length = b := by
  induction b generalizing v with
  | zero => simp [allWords] at h; simp [h]
  | succ n ih =>
      simp only [allWords, List.mem_flatMap] at h
      obtain ⟨u, hu, hv⟩ := h
      have hlen := ih hu
      simp only [List.mem_cons, List.not_mem_nil, or_false] at hv
      rcases hv with rfl | rfl <;> simp [hlen]

/-- The number of length-`b` words of weight `a` whose accumulation has
weight `c`.  This is `T_b(a,c)` as a *definition by counting*; the closed form
is proved equal to it below. -/
def accCount (b a c : ℕ) : ℕ :=
  (allWords b).countP (fun v => decide (wtL v = a ∧ accWtL v = c))

/-! ## The counting recursion -/

private lemma countP_eq_sum_map (p : List Bool → Bool) (l : List (List Bool)) :
    l.countP p = (l.map (fun v => if p v then 1 else 0)).sum := by
  induction l with
  | nil => simp
  | cons x t ih => by_cases h : p x <;> simp [List.countP_cons, h, ih] <;> omega

private lemma sum_map_add {α : Type*} (l : List α) (f g : α → ℕ) :
    (l.map (fun a => f a + g a)).sum = (l.map f).sum + (l.map g).sum := by
  induction l with
  | nil => simp
  | cons x t ih => simp [ih]; omega

/-- Appending a letter extends the accumulation by the parity of the *whole*
resulting word — which is why the recursion's shift depends only on `a`. -/
lemma accWtL_append_mod (v : List Bool) (x : Bool) :
    accWtL (v ++ [x]) = accWtL v + wtL (v ++ [x]) % 2 := by
  rw [accWtL_append, parityL_toNat]

/-- Splitting the word list on the last letter. -/
lemma accCount_succ_split (b : ℕ) (P : List Bool → Bool) :
    (allWords (b + 1)).countP P
      = (allWords b).countP (fun v => P (v ++ [false]))
        + (allWords b).countP (fun v => P (v ++ [true])) := by
  show (( allWords b).flatMap (fun v => [v ++ [false], v ++ [true]])).countP P = _
  rw [List.countP_flatMap]
  have hstep : ∀ v : List Bool,
      (List.countP P ∘ fun v => [v ++ [false], v ++ [true]]) v
        = (if P (v ++ [false]) then 1 else 0) + (if P (v ++ [true]) then 1 else 0) := by
    intro v
    by_cases h0 : P (v ++ [false]) <;> by_cases h1 : P (v ++ [true]) <;>
      simp [List.countP_cons, h0, h1]
  rw [List.map_congr_left (fun v _ => hstep v), sum_map_add,
    countP_eq_sum_map, countP_eq_sum_map]

lemma accCount_succ_zero (b c : ℕ) : accCount (b + 1) 0 c = accCount b 0 c := by
  unfold accCount
  rw [accCount_succ_split]
  have h0 : ∀ v : List Bool,
      decide (wtL (v ++ [false]) = 0 ∧ accWtL (v ++ [false]) = c)
        = decide (wtL v = 0 ∧ accWtL v = c) := by
    intro v
    rw [wtL_append, accWtL_append_mod, wtL_append]
    by_cases h : wtL v = 0 <;> simp [h]
  have h1 : ∀ v : List Bool,
      decide (wtL (v ++ [true]) = 0 ∧ accWtL (v ++ [true]) = c) = false := by
    intro v
    rw [wtL_append]
    simp
  simp only [h0, h1]
  simp

/-- **The counting recursion.**  Conditioning on the last letter:

  `N_{b+1}(a+1, c + (a+1) mod 2) = N_b(a+1, c) + N_b(a, c)`. -/
lemma accCount_succ (b a c : ℕ) :
    accCount (b + 1) (a + 1) (c + (a + 1) % 2)
      = accCount b (a + 1) c + accCount b a c := by
  unfold accCount
  rw [accCount_succ_split]
  congr 1
  · refine List.countP_congr fun v _ => ?_
    rw [wtL_append, accWtL_append_mod, wtL_append]
    by_cases h : wtL v = a + 1
    · simp [h]
    · simp [h]
  · refine List.countP_congr fun v _ => ?_
    rw [wtL_append, accWtL_append_mod, wtL_append]
    by_cases h : wtL v = a
    · simp [h]
    · simp [h]

/-! ## Nonzero words have nonzero accumulation -/

lemma wtL_eq_zero_of_accWtL (v : List Bool) (h : accWtL v = 0) : wtL v = 0 := by
  have key : ∀ u : List Bool, accWtAux false u = 0 → wtL u = 0 := by
    intro u
    induction u with
    | nil => intro _; simp [wtL]
    | cons x t ih =>
        intro hu
        simp only [accWtAux, Bool.false_xor] at hu
        have hx : x = false := by
          cases x
          · rfl
          · simp at hu
        subst hx
        simp only [Bool.false_eq_true, if_false, zero_add] at hu
        simp [wtL, ih hu]
  exact key v h

lemma accCount_zero_right (b : ℕ) {a : ℕ} (ha : a ≠ 0) : accCount b a 0 = 0 := by
  unfold accCount
  refine List.countP_eq_zero.mpr fun v _ => ?_
  simp only [decide_eq_true_eq, Bool.not_eq_true, decide_eq_false_iff_not, not_and]
  intro hw hc
  exact ha (by rw [← hw, wtL_eq_zero_of_accWtL v hc])

/-! ## The closed form -/

/-- `T_b(a,c) = C(c-1, r_a - 1) * C(b-c, a - r_a)` with `r_a = ⌈a/2⌉`, and the
paper's "binomial coefficients outside their natural range are zero"
convention made explicit.  Both guards are load-bearing under truncated
subtraction on `ℕ`: without the `c = 0` guard the first factor would become
`C(0, r-1)`, and without `b < c` the second would become `C(0, a - r)`. -/
def accT (b a c : ℕ) : ℕ :=
  if a = 0 then (if c = 0 then 1 else 0)
  else if c = 0 ∨ b < c then 0
  else (c - 1).choose ((a + 1) / 2 - 1) * (b - c).choose (a - (a + 1) / 2)

lemma accT_zero_left (b c : ℕ) : accT b 0 c = if c = 0 then 1 else 0 := by
  unfold accT; simp

lemma accT_succ_zero (b c : ℕ) : accT (b + 1) 0 c = accT b 0 c := by
  rw [accT_zero_left, accT_zero_left]

lemma accT_zero_right (b : ℕ) {a : ℕ} (ha : a ≠ 0) : accT b a 0 = 0 := by
  unfold accT; simp [ha]

lemma accT_of_lt (b a c : ℕ) (ha : a ≠ 0) (h : b < c) : accT b a c = 0 := by
  unfold accT
  rw [if_neg ha, if_pos (Or.inr h)]

lemma accT_zero_length {a : ℕ} (ha : a ≠ 0) (c : ℕ) : accT 0 a c = 0 := by
  rcases Nat.eq_zero_or_pos c with rfl | hc
  · exact accT_zero_right 0 ha
  · exact accT_of_lt 0 a c ha (by omega)

/-- Evaluation in the only regime where the closed form is not forced to zero:
`1 ≤ c ≤ b`, written as `c = e + 1` and `b = c + d` so that every truncated
subtraction disappears. -/
lemma accT_eval {a : ℕ} (ha : a ≠ 0) (e d : ℕ) :
    accT (e + 1 + d) a (e + 1)
      = e.choose ((a + 1) / 2 - 1) * d.choose (a - (a + 1) / 2) := by
  unfold accT
  rw [if_neg ha, if_neg (by omega)]
  have h1 : e + 1 - 1 = e := by omega
  have h2 : e + 1 + d - (e + 1) = d := by omega
  rw [h1, h2]

/-- The closed form satisfies the counting recursion, in both parities, by a
single application of Pascal's rule each time. -/
lemma accT_succ (b a c : ℕ) :
    accT (b + 1) (a + 1) (c + (a + 1) % 2) = accT b (a + 1) c + accT b a c := by
  obtain ⟨k, hk | hk⟩ := Nat.even_or_odd' a
  · ------------------------------------------------------------------
    -- a = 2k, so a+1 = 2k+1 is odd and the accumulation shifts by 1
    subst hk
    rw [show (2 * k + 1) % 2 = 1 from by omega]
    have hAidx1 : (2 * k + 1 + 1) / 2 - 1 = k := by omega
    have hAidx2 : 2 * k + 1 - (2 * k + 1 + 1) / 2 = k := by omega
    rcases Nat.eq_zero_or_pos c with rfl | hcpos
    · -- c = 0
      rw [accT_zero_right _ (by omega : (2 * k + 1) ≠ 0)]
      have hL : accT (b + 1) (2 * k + 1) (0 + 1)
          = (0 : ℕ).choose k * b.choose k := by
        rw [show (0 : ℕ) + 1 = 0 + 1 from rfl, show b + 1 = 0 + 1 + b from by omega,
          accT_eval (by omega : (2 * k + 1) ≠ 0) 0 b, hAidx1, hAidx2]
      rw [hL]
      rcases Nat.eq_zero_or_pos k with rfl | hk0
      · rw [show 2 * 0 = 0 from by omega, accT_zero_left]
        simp
      · rw [Nat.choose_eq_zero_of_lt (by omega), accT_zero_right _ (by omega)]
        simp
    · -- c ≥ 1
      obtain ⟨e, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : c ≠ 0)
      rcases Nat.lt_or_ge b (e + 1) with hbc | hbc
      · rw [accT_of_lt _ _ _ (by omega : (2 * k + 1) ≠ 0) (by omega),
          accT_of_lt _ _ _ (by omega : (2 * k + 1) ≠ 0) hbc]
        rcases Nat.eq_zero_or_pos k with rfl | hk0
        · rw [show 2 * 0 = 0 from by omega, accT_zero_left, if_neg (by omega)]
        · rw [accT_of_lt _ _ _ (by omega) hbc]
      · obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hbc
        have hL : accT (e + 1 + d + 1) (2 * k + 1) (e + 1 + 1)
            = (e + 1).choose k * d.choose k := by
          rw [show e + 1 + d + 1 = (e + 1) + 1 + d from by omega,
            accT_eval (by omega : (2 * k + 1) ≠ 0) (e + 1) d, hAidx1, hAidx2]
        rw [hL, accT_eval (by omega : (2 * k + 1) ≠ 0) e d, hAidx1, hAidx2]
        rcases Nat.eq_zero_or_pos k with rfl | hk0
        · rw [show 2 * 0 = 0 from by omega, accT_zero_left, if_neg (by omega)]
          simp
        · obtain ⟨j, rfl⟩ := Nat.exists_eq_add_of_le hk0
          have hidx1 : (2 * (1 + j) + 1) / 2 - 1 = j := by omega
          have hidx2 : 2 * (1 + j) - (2 * (1 + j) + 1) / 2 = 1 + j := by omega
          rw [accT_eval (by omega : (2 * (1 + j)) ≠ 0) e d, hidx1, hidx2]
          have hpascal : (e + 1).choose (1 + j) = e.choose j + e.choose (1 + j) := by
            rw [show (1 : ℕ) + j = j + 1 from by omega, Nat.choose_succ_succ]
          rw [hpascal]
          ring
  · ------------------------------------------------------------------
    -- a = 2k+1, so a+1 = 2k+2 is even and there is no shift
    subst hk
    rw [show (2 * k + 1 + 1) % 2 = 0 from by omega, Nat.add_zero]
    have hAidx1 : (2 * k + 1 + 1 + 1) / 2 - 1 = k := by omega
    have hAidx2 : 2 * k + 1 + 1 - (2 * k + 1 + 1 + 1) / 2 = k + 1 := by omega
    have haidx1 : (2 * k + 1 + 1) / 2 - 1 = k := by omega
    have haidx2 : 2 * k + 1 - (2 * k + 1 + 1) / 2 = k := by omega
    rcases Nat.eq_zero_or_pos c with rfl | hcpos
    · rw [accT_zero_right _ (by omega : (2 * k + 1 + 1) ≠ 0),
        accT_zero_right _ (by omega : (2 * k + 1 + 1) ≠ 0),
        accT_zero_right _ (by omega : (2 * k + 1) ≠ 0)]
    · obtain ⟨e, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : c ≠ 0)
      rcases Nat.lt_or_ge b (e + 1) with hbc | hbc
      · rw [accT_of_lt _ _ _ (by omega : (2 * k + 1 + 1) ≠ 0) hbc,
          accT_of_lt _ _ _ (by omega : (2 * k + 1) ≠ 0) hbc]
        rcases Nat.lt_or_ge (b + 1) (e + 1) with hb1 | hb1
        · rw [accT_of_lt _ _ _ (by omega : (2 * k + 1 + 1) ≠ 0) hb1]
        · -- b + 1 = e + 1, so the second factor is `C(0, k+1) = 0`
          have hbe : b = e := by omega
          subst hbe
          rw [show b + 1 = b + 1 + 0 from by omega,
            accT_eval (by omega : (2 * k + 1 + 1) ≠ 0) b 0, hAidx1, hAidx2,
            Nat.choose_eq_zero_of_lt (n := 0) (k := k + 1) (by omega)]
          simp
      · obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hbc
        rw [show e + 1 + d + 1 = e + 1 + (d + 1) from by omega,
          accT_eval (by omega : (2 * k + 1 + 1) ≠ 0) e (d + 1), hAidx1, hAidx2,
          accT_eval (by omega : (2 * k + 1 + 1) ≠ 0) e d, hAidx1, hAidx2,
          accT_eval (by omega : (2 * k + 1) ≠ 0) e d, haidx1, haidx2,
          Nat.choose_succ_succ]
        ring

/-- **`lem:structured-exact-ba-spectrum`, accumulator half.**  The closed form
counts exactly the words it is claimed to count, so `P_b(a,c) = T_b(a,c)/C(b,a)`
is the exact weight-transition law of a uniformly interleaved accumulator. -/
theorem accCount_eq_accT (b a c : ℕ) : accCount b a c = accT b a c := by
  induction b generalizing a c with
  | zero =>
      have hall : allWords 0 = [[]] := rfl
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [accT_zero_left]
        unfold accCount
        rw [hall]
        rcases Nat.eq_zero_or_pos c with rfl | hc
        · norm_num [List.countP_cons, wtL, accWtL, accWtAux]
        · rw [if_neg (by omega : ¬ c = 0)]
          have hne : ¬ ((0 : ℕ) = c) := by omega
          simp [List.countP_cons, wtL, accWtL, accWtAux, hne]
      · rw [accT_zero_length (by omega) c]
        unfold accCount
        rw [hall]
        have hne : ¬ ((0 : ℕ) = a) := by omega
        simp [List.countP_cons, wtL, accWtL, accWtAux, hne]
  | succ n ih =>
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [accCount_succ_zero, accT_succ_zero, ih]
      · obtain ⟨a', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : a ≠ 0)
        have e1 := accCount_succ n a' c
        have e2 := accT_succ n a' c
        rcases Nat.eq_zero_or_pos ((a' + 1) % 2) with hpar | hpar
        · simp only [hpar, Nat.add_zero] at e1 e2
          rw [e1, e2, ih, ih]
        · rcases Nat.eq_zero_or_pos c with rfl | hcpos
          · rw [accCount_zero_right _ (by omega), accT_zero_right _ (by omega)]
          · obtain ⟨c', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : c ≠ 0)
            have hmod : (a' + 1) % 2 = 1 := by omega
            have e1' := accCount_succ n a' c'
            have e2' := accT_succ n a' c'
            rw [hmod] at e1' e2'
            rw [e1', e2', ih, ih]


/-! ## Regression checks

These pin the *semantics*, not the theorem: if `accWtL` or `allWords` were
subtly wrong, `accCount_eq_accT` would still hold (both sides would move
together), so a few hand-computed values are checked directly.

Length 3, weight 1: `100 ↦ 111` (weight 3), `010 ↦ 011` (weight 2),
`001 ↦ 001` (weight 1). -/
example : accCount 3 1 3 = 1 := by decide
example : accCount 3 1 2 = 1 := by decide
example : accCount 3 1 1 = 1 := by decide
example : accCount 3 1 0 = 0 := by decide

/-- Every word is counted exactly once: the rows sum to `C(b,a)`. -/
example : (List.range 5).map (fun c => accCount 4 2 c) = [0, 3, 2, 1, 0] := by decide
example : ((List.range 5).map (fun c => accCount 4 2 c)).sum = Nat.choose 4 2 := by decide

/-- The accumulation is a bijection, so only the zero word has zero
accumulated weight. -/
example : accCount 4 0 0 = 1 := by decide


/-- The complement identity on a concrete word: `1011 ↦ Acc = 1101` (weight 3),
and flipping the first bit gives `0011 ↦ Acc = 0010` (weight 1 = 4 - 3). -/
example : accWtL [true, false, true, true] = 3 := by decide
example : accWtL (flipHead [true, false, true, true]) = 1 := by decide

end Spin.Structured
