/-
The accumulator transition over tuples.

`Accumulator.lean` proves the count over `List Bool`, because appending on the
right is structural there and that is what makes the recursion work.  But a
*uniform interleaver* is a permutation of positions, which acts on
`Fin b → Bool`, not on lists.  This file carries the count across, so the
transition law is available where the permutation argument lives.

`List.ofFn (Fin.snoc v x) = List.ofFn v ++ [x]` is the whole bridge: the
append lemmas transport verbatim, and only the cardinality split needs new
work.
-/
import SpinCodes.Structured.Accumulator

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-- Weight of a tuple. -/
def wtF {b : ℕ} (v : Fin b → Bool) : ℕ := wtL (List.ofFn v)

/-- Accumulated weight of a tuple. -/
def accWtF {b : ℕ} (v : Fin b → Bool) : ℕ := accWtL (List.ofFn v)

lemma ofFn_snoc {b : ℕ} (v : Fin b → Bool) (x : Bool) :
    List.ofFn (Fin.snoc v x) = List.ofFn v ++ [x] := by
  rw [List.ofFn_succ']
  simp [List.concat_eq_append]

lemma wtF_snoc {b : ℕ} (v : Fin b → Bool) (x : Bool) :
    wtF (Fin.snoc v x) = wtF v + (if x then 1 else 0) := by
  unfold wtF
  rw [ofFn_snoc, wtL_append]

lemma accWtF_snoc {b : ℕ} (v : Fin b → Bool) (x : Bool) :
    accWtF (Fin.snoc v x) = accWtF v + wtF (Fin.snoc v x) % 2 := by
  unfold accWtF wtF
  rw [ofFn_snoc, accWtL_append_mod]

/-- The number of tuples of weight `a` whose accumulation has weight `c`. -/
def accCountF (b a c : ℕ) : ℕ :=
  (univ.filter (fun v : Fin b → Bool => wtF v = a ∧ accWtF v = c)).card

/-- Splitting the tuple space on the last coordinate. -/
lemma card_filter_succ_split (b : ℕ) (P : (Fin (b + 1) → Bool) → Prop)
    [DecidablePred P] :
    (univ.filter P).card
      = (univ.filter (fun v : Fin b → Bool => P (Fin.snoc v false))).card
        + (univ.filter (fun v : Fin b → Bool => P (Fin.snoc v true))).card := by
  classical
  simp only [Finset.card_filter]
  rw [← Fintype.sum_equiv (Fin.snocEquiv (fun _ : Fin (b + 1) => Bool))
    (fun p : Bool × (Fin b → Bool) => if P (Fin.snoc p.2 p.1) then 1 else 0)
    (fun f => if P f then 1 else 0) (fun p => rfl)]
  rw [Fintype.sum_prod_type, Fintype.sum_bool]
  simp only [Finset.card_filter]
  omega

lemma accCountF_succ_zero (b c : ℕ) : accCountF (b + 1) 0 c = accCountF b 0 c := by
  unfold accCountF
  rw [card_filter_succ_split]
  have h0 : ∀ v : Fin b → Bool,
      (wtF (Fin.snoc v false) = 0 ∧ accWtF (Fin.snoc v false) = c)
        ↔ (wtF v = 0 ∧ accWtF v = c) := by
    intro v
    rw [accWtF_snoc, wtF_snoc]
    by_cases h : wtF v = 0 <;> simp [h]
  have h1 : ∀ v : Fin b → Bool,
      ¬ (wtF (Fin.snoc v true) = 0 ∧ accWtF (Fin.snoc v true) = c) := by
    intro v
    rw [wtF_snoc]
    simp
  rw [Finset.filter_congr (fun v _ => by simpa using h0 v),
    Finset.filter_false_of_mem (fun v _ => h1 v)]
  simp

/-- The counting recursion, over tuples. -/
lemma accCountF_succ (b a c : ℕ) :
    accCountF (b + 1) (a + 1) (c + (a + 1) % 2)
      = accCountF b (a + 1) c + accCountF b a c := by
  unfold accCountF
  rw [card_filter_succ_split]
  congr 1
  · refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
    rw [accWtF_snoc, wtF_snoc]
    by_cases h : wtF v = a + 1 <;> simp [h]
  · refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
    rw [accWtF_snoc, wtF_snoc]
    by_cases h : wtF v = a <;> simp [h]

lemma accCountF_zero_right (b : ℕ) {a : ℕ} (ha : a ≠ 0) : accCountF b a 0 = 0 := by
  unfold accCountF
  refine Finset.card_eq_zero.mpr (Finset.filter_eq_empty_iff.mpr fun v _ => ?_)
  rintro ⟨hw, hc⟩
  exact ha (by rw [← hw]; exact wtL_eq_zero_of_accWtL _ hc)

/-- **The transition law, over tuples.**  `T_b(a,c)` counts the tuples of
weight `a` whose accumulation has weight `c`. -/
theorem accCountF_eq_accT (b a c : ℕ) : accCountF b a c = accT b a c := by
  induction b generalizing a c with
  | zero =>
      have hzero : ∀ v : Fin 0 → Bool, wtF v = 0 ∧ accWtF v = 0 := by
        intro v
        have : List.ofFn v = [] := by simp
        unfold wtF accWtF
        rw [this]
        exact ⟨rfl, rfl⟩
      unfold accCountF
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [accT_zero_left]
        rcases Nat.eq_zero_or_pos c with rfl | hc
        · rw [if_pos rfl, Finset.filter_true_of_mem (fun v _ => hzero v)]
          simp
        · rw [if_neg (by omega)]
          refine Finset.card_eq_zero.mpr (Finset.filter_eq_empty_iff.mpr fun v _ => ?_)
          rintro ⟨-, h2⟩
          rw [(hzero v).2] at h2
          omega
      · rw [accT_zero_length (by omega) c]
        refine Finset.card_eq_zero.mpr (Finset.filter_eq_empty_iff.mpr fun v _ => ?_)
        rintro ⟨h1, -⟩
        rw [(hzero v).1] at h1
        omega
  | succ n ih =>
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [accCountF_succ_zero, accT_succ_zero, ih]
      · obtain ⟨a', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : a ≠ 0)
        have e1 := accCountF_succ n a' c
        have e2 := accT_succ n a' c
        rcases Nat.eq_zero_or_pos ((a' + 1) % 2) with hpar | hpar
        · simp only [hpar, Nat.add_zero] at e1 e2
          rw [e1, e2, ih, ih]
        · rcases Nat.eq_zero_or_pos c with rfl | hcpos
          · rw [accCountF_zero_right _ (by omega), accT_zero_right _ (by omega)]
          · obtain ⟨c', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : c ≠ 0)
            have hmod : (a' + 1) % 2 = 1 := by omega
            have e1' := accCountF_succ n a' c'
            have e2' := accT_succ n a' c'
            rw [hmod] at e1' e2'
            rw [e1', e2', ih, ih]

/-- The list and tuple counts agree, as they must. -/
theorem accCountF_eq_accCount (b a c : ℕ) : accCountF b a c = accCount b a c := by
  rw [accCountF_eq_accT, accCount_eq_accT]


/-! ## The accumulation as a map on tuples -/

/-- The accumulated tuple. -/
def accF {b : ℕ} (v : Fin b → Bool) : Fin b → Bool :=
  fun i => (accL (List.ofFn v)).getD (i : ℕ) false

lemma ofFn_accF {b : ℕ} (v : Fin b → Bool) :
    List.ofFn (accF v) = accL (List.ofFn v) := by
  refine List.ext_getElem (by simp) fun n h1 h2 => ?_
  simp only [List.getElem_ofFn, accF]
  rw [List.getD_eq_getElem _ _ (by simpa using h2)]

/-- The accumulated tuple has the accumulated weight. -/
@[simp] lemma wtF_accF {b : ℕ} (v : Fin b → Bool) : wtF (accF v) = accWtF v := by
  unfold wtF accWtF
  rw [ofFn_accF, wtL_accL]

lemma wtF_le {b : ℕ} (v : Fin b → Bool) : wtF v ≤ b := by
  unfold wtF
  simpa using wtL_le_length (List.ofFn v)


/-! ## Counting over tuples versus over the word list

The list enumeration `allWords n` is cheap for the kernel to evaluate; the
`Fintype` enumeration of `Fin n → Bool` is not (`Fintype.piFinset` reduction is
orders of magnitude slower).  This lemma lets a count be *computed* over lists
and then *used* over tuples. -/

lemma card_filter_ofFn_eq_countP :
    ∀ {n : ℕ} (P : List Bool → Prop) [DecidablePred P],
      (univ.filter (fun v : Fin n → Bool => P (List.ofFn v))).card
        = (allWords n).countP (fun l => decide (P l)) := by
  intro n
  induction n with
  | zero =>
      intro P _
      classical
      have hnil : ∀ v : Fin 0 → Bool, List.ofFn v = [] := by intro v; simp
      show _ = (allWords 0).countP _
      unfold allWords
      by_cases h : P []
      · rw [Finset.filter_true_of_mem (fun v _ => by rw [hnil v]; exact h)]
        simp [h]
      · rw [Finset.filter_false_of_mem (fun v _ => by rw [hnil v]; exact h)]
        simp [h]
  | succ k ih =>
      intro P _
      classical
      rw [card_filter_succ_split, accCount_succ_split]
      have e1 : (univ.filter
            (fun v : Fin k → Bool => P (List.ofFn (Fin.snoc v false)))).card
          = (allWords k).countP (fun l => decide (P (l ++ [false]))) := by
        rw [← ih (fun l => P (l ++ [false]))]
        refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
        rw [ofFn_snoc]
      have e2 : (univ.filter
            (fun v : Fin k → Bool => P (List.ofFn (Fin.snoc v true)))).card
          = (allWords k).countP (fun l => decide (P (l ++ [true]))) := by
        rw [← ih (fun l => P (l ++ [true]))]
        refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
        rw [ofFn_snoc]
      rw [e1, e2]

end Spin.Structured
