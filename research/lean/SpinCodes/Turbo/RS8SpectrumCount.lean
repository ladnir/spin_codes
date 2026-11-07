import SpinCodes.Turbo.RS8PositionSets

/-!
# Intersections in the RS8 inclusion-exclusion sum

The inclusion-exclusion term for a subset `t` of candidate active columns
forces zeros on `sᶜ ∪ t`. This file identifies its count with the already
proved arbitrary-set shortening count.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F] [Fintype F] [DecidableEq F]

omit [Fintype F] in
theorem rs8SupportedIn_iff_zero_compl
    (points : Fin 8 → F) (s : Finset (Fin 8)) (rows : RS8Group F) :
    rs8SupportedIn points s rows ↔ rs8ZeroOnSet points sᶜ rows := by
  classical
  constructor
  · intro hs j hj
    by_contra hn
    have hmem : j ∈ rs8Support points rows :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩
    exact (Finset.mem_compl.mp hj) (hs hmem)
  · intro hz j hj
    by_contra hn
    have hcomp : j ∈ sᶜ := Finset.mem_compl.mpr hn
    have hzero := hz j hcomp
    exact (Finset.mem_filter.mp hj).2 hzero

/-- Each intersection in the exact-support formula is an ordinary
shortening on `sᶜ ∪ t`. -/
theorem rs8_missing_intersection_card
    (points : Fin 8 → F) (s t : Finset (Fin 8)) :
    (t.inf (rs8MissingAt points s)).card =
      Fintype.card {rows : RS8Group F //
        rs8ZeroOnSet points (sᶜ ∪ t) rows} := by
  classical
  let α := {rows : RS8Group F // rs8SupportedIn points s rows}
  let st : Finset α := t.inf (rs8MissingAt points s)
  have hiff (x : α) :
      x ∈ st ↔ rs8ZeroOnSet points (sᶜ ∪ t) x.val := by
    constructor
    · intro hx j hj
      rcases Finset.mem_union.mp hj with hc | ht
      · exact (rs8SupportedIn_iff_zero_compl points s x.val).mp x.property j hc
      · have hb : x ∈ rs8MissingAt points s j :=
          (Finset.mem_inf.mp hx) j ht
        exact (Finset.mem_filter.mp hb).2
    · intro hz
      apply Finset.mem_inf.mpr
      intro j hj
      apply Finset.mem_filter.mpr
      exact ⟨Finset.mem_univ _, hz j (Finset.mem_union_right _ hj)⟩
  let e : ↥st ≃ {rows : RS8Group F //
      rs8ZeroOnSet points (sᶜ ∪ t) rows} :=
    { toFun := fun x => ⟨x.val.val, (hiff x.val).mp x.property⟩
      invFun := fun x =>
        ⟨⟨x.val, (rs8SupportedIn_iff_zero_compl points s x.val).mpr
          (fun j hj => x.property j (Finset.mem_union_left _ hj))⟩,
          (hiff ⟨x.val, (rs8SupportedIn_iff_zero_compl points s x.val).mpr
            (fun j hj => x.property j (Finset.mem_union_left _ hj))⟩).mpr x.property⟩
      left_inv := fun _ => Subtype.ext (Subtype.ext rfl)
      right_inv := fun _ => Subtype.ext rfl }
  calc
    (t.inf (rs8MissingAt points s)).card = Fintype.card ↥st :=
      (Fintype.card_coe st).symm
    _ = _ := Fintype.card_congr e

omit [Fintype F] in
/-- More than four forced zero columns leave only the zero message. -/
theorem rs8ZeroOnSet_eq_zero_of_card_gt_four
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (z : Finset (Fin 8)) (hz : 4 < z.card)
    (rows : RS8Group F) (hzero : rs8ZeroOnSet points z rows) :
    rows = 0 := by
  classical
  by_contra hn
  let zeros : Finset (Fin 8) :=
    Finset.univ.filter (fun j => rs8Column points rows j = 0)
  have hsub : z ⊆ zeros := by
    intro j hj
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hzero j hj⟩
  have hzcard : 5 ≤ zeros.card := by
    have := Finset.card_le_card hsub
    omega
  have hpart : zeros.card + (rs8Support points rows).card = 8 := by
    simpa [zeros, rs8Support] using
      (Finset.card_filter_add_card_filter_not
        (s := (Finset.univ : Finset (Fin 8)))
        (p := fun j => rs8Column points rows j = 0))
  have hfloor := rs8Support_ge_five points hpoints rows hn
  omega

theorem rs8ZeroOnSet_card_of_gt_four
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (z : Finset (Fin 8)) (hz : 4 < z.card) :
    Fintype.card {rows : RS8Group F // rs8ZeroOnSet points z rows} = 1 := by
  classical
  let e : {rows : RS8Group F // rs8ZeroOnSet points z rows} ≃ Unit :=
    { toFun := fun _ => ()
      invFun := fun _ => ⟨0, by
        intro j hj
        funext r
        simp [rs8Column]⟩
      left_inv := fun x => Subtype.ext
        (rs8ZeroOnSet_eq_zero_of_card_gt_four points hpoints z hz x.val x.property).symm
      right_inv := fun _ => rfl }
  simpa using Fintype.card_congr e

/-- The exact cardinality of every intersection in inclusion-exclusion. -/
theorem rs8_missing_intersection_card_formula
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (s t : Finset (Fin 8)) :
    (t.inf (rs8MissingAt points s)).card =
      if (sᶜ ∪ t).card ≤ 4 then
        (Fintype.card F) ^ (8 * (4-(sᶜ ∪ t).card))
      else 1 := by
  rw [rs8_missing_intersection_card]
  split_ifs with h
  · exact rs8ZeroOnSet_card points hpoints _ h
  · exact rs8ZeroOnSet_card_of_gt_four points hpoints _ (by omega)

/-- A checked exact-support enumerator for the RS8 group, written as a
finite set sum. The paper's closed `A_h` expression is a simplification of
this formula using `|F|=16` and the cardinality of each subset. -/
theorem rs8_exact_support_formula
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (s : Finset (Fin 8)) :
    ((Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card : ℤ) =
    ∑ t ∈ s.powerset, (-1 : ℤ) ^ t.card *
      ((if (sᶜ ∪ t).card ≤ 4 then
          (Fintype.card F) ^ (8 * (4-(sᶜ ∪ t).card))
        else 1 : ℕ) : ℤ) := by
  rw [rs8_exact_support_inclusion_exclusion]
  apply Finset.sum_congr rfl
  intro t ht
  rw [rs8_missing_intersection_card_formula points hpoints s t]

#print axioms rs8_missing_intersection_card
#print axioms rs8ZeroOnSet_card_of_gt_four
#print axioms rs8_missing_intersection_card_formula
#print axioms rs8_exact_support_formula

end Spin.Turbo
