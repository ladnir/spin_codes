/-
The Golay weight enumerator, on the tuple side.

`Golay.lean` computes the weight distribution over the list enumeration, which
is what the kernel can afford.  This file carries it to `Fin 12 → Bool`, where
`enumerator` and hence `card_tuples_weight` live, and assembles

    G(u) = 1 + 759 u^8 + 2576 u^12 + 759 u^16 + u^24.

Kept separate from `Golay.lean` so that editing the assembly does not re-run
the 4096-codeword kernel computation.
-/
import SpinCodes.Structured.Golay

set_option linter.unusedSectionVars false

namespace Spin.Structured.Golay

open Finset Polynomial Spin.Structured

/-- The Golay encoder on tuples. -/
def enc (m : Fin 12 → Bool) : Fin 24 → Bool :=
  fun j => (encList (List.ofFn m)).getD (j : ℕ) false

lemma ofFn_enc (m : Fin 12 → Bool) : List.ofFn (enc m) = encList (List.ofFn m) := by
  refine List.ext_getElem (by simp [length_encList]) fun n h1 h2 => ?_
  simp only [List.getElem_ofFn, enc]
  rw [List.getD_eq_getElem _ _ (by simpa [length_encList] using h2)]

/-- The block weight function of the Golay outer: the Hamming weight of the
encoded message. -/
def W (m : Fin 12 → Bool) : ℕ := wtF (enc m)

lemma W_eq (m : Fin 12 → Bool) : W m = wtL (encList (List.ofFn m)) := by
  unfold W wtF
  rw [ofFn_enc]

lemma W_le (m : Fin 12 → Bool) : W m ≤ 24 := by
  rw [W_eq]
  simpa [length_encList] using wtL_le_length (encList (List.ofFn m))

/-- The Golay weight distribution as a function of the weight. -/
def A (i : ℕ) : ℕ :=
  if i = 0 then 1 else if i = 8 then 759 else if i = 12 then 2576
  else if i = 16 then 759 else if i = 24 then 1 else 0

lemma hist_getD (i : ℕ) (hi : i < 25) : hist.getD i 0 = A i := by
  rw [hist_eq]
  interval_cases i <;> rfl

/-- **The Golay weight distribution, on the tuple side.** -/
theorem card_weight (i : ℕ) :
    (univ.filter (fun m : Fin 12 → Bool => W m = i)).card = A i := by
  classical
  rcases Nat.lt_or_ge i 25 with hi | hi
  · have hfilter : (univ.filter (fun m : Fin 12 → Bool => W m = i))
        = univ.filter (fun m : Fin 12 → Bool => wtL (encList (List.ofFn m)) = i) := by
      exact Finset.filter_congr fun m _ => by rw [W_eq]
    rw [hfilter, card_filter_ofFn_eq_countP (fun l => wtL (encList l) = i)]
    have hfold : hist.getD i 0
        = (allWords 12).countP (fun m => wtL (encList m) == i) := by
      show ((allWords 12).foldl (fun h m => bumpAt (wtL (encList m)) h)
        (List.replicate 25 0)).getD i 0 = _
      have h := getD_foldl_bumpAt (fun m : List Bool => wtL (encList m))
        (allWords 12) (List.replicate 25 0) i (by simpa using hi)
      have hz : (List.replicate 25 (0:ℕ)).getD i 0 = 0 := by
        interval_cases i <;> rfl
      rw [hz, Nat.zero_add] at h
      exact h
    have hcp : (allWords 12).countP (fun l => decide (wtL (encList l) = i))
        = (allWords 12).countP (fun m => wtL (encList m) == i) := by
      refine List.countP_congr fun l _ => ?_
      simp
    rw [hcp, ← hfold]
    exact hist_getD i hi
  · have hA : A i = 0 := by
      unfold A
      rw [if_neg (by omega), if_neg (by omega), if_neg (by omega), if_neg (by omega),
        if_neg (by omega)]
    rw [hA]
    refine Finset.card_eq_zero.mpr (Finset.filter_eq_empty_iff.mpr fun m _ => ?_)
    have := W_le m
    omega

/-- **`G(u) = 1 + 759 u^8 + 2576 u^12 + 759 u^16 + u^24`.**

Computed, not assumed: the coefficients come from a kernel evaluation of all
4096 codewords. -/
theorem enumerator_eq :
    enumerator W = 1 + C 759 * X ^ 8 + C 2576 * X ^ 12 + C 759 * X ^ 16 + X ^ 24 := by
  refine Polynomial.ext fun n => ?_
  rw [coeff_enumerator, card_weight]
  simp only [Polynomial.coeff_add, Polynomial.coeff_one, Polynomial.coeff_C_mul,
    Polynomial.coeff_X_pow, A]
  split_ifs <;> omega

/-- The code is `[24, 12, 8]`: minimum distance 8, since no nonzero codeword
has weight below 8. -/
theorem min_distance (m : Fin 12 → Bool) (hm : W m ≠ 0) : 8 ≤ W m := by
  by_contra hlt
  push_neg at hlt
  have hcard : (univ.filter (fun m' : Fin 12 → Bool => W m' = W m)).card = A (W m) :=
    card_weight (W m)
  have hne : (univ.filter (fun m' : Fin 12 → Bool => W m' = W m)).Nonempty :=
    ⟨m, by simp⟩
  have hpos : 0 < A (W m) := by
    rw [← hcard]
    exact Finset.card_pos.mpr hne
  have hA : A (W m) = 0 := by
    unfold A
    rw [if_neg hm, if_neg (by omega), if_neg (by omega), if_neg (by omega),
      if_neg (by omega)]
  omega

end Spin.Structured.Golay
