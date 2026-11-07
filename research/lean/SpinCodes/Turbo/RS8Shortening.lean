import SpinCodes.Turbo.RS8Code

/-!
# Shortening the eight-row RS outer

The old Structured outer proof counted each base-code weight class before
feeding the first-moment identity. For Turbo, the first count is the
three-zero-column shortening of the RS[8,4] group. It has one free
32-bit column, hence `|F|^8` messages over the eight-row message space.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F] [Fintype F]

def rs8ZeroOnThree (points : Fin 8 → F)
    (positions : Fin 4 → Fin 8) (rows : RS8Group F) : Prop :=
  ∀ i : Fin 3, rs8Column points rows (positions ⟨i.val, by omega⟩) = 0

noncomputable instance (points : Fin 8 → F) (positions : Fin 4 → Fin 8) :
    Fintype {rows : RS8Group F // rs8ZeroOnThree points positions rows} :=
  Fintype.ofFinite _

private def extendFourth (v : Fin 8 → F) : Fin 4 → Fin 8 → F :=
  fun i => if i.val = 3 then v else 0

/-- A shortened group is specified by its single free column. -/
noncomputable def rs8ThreeZeroEquiv
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    {rows : RS8Group F // rs8ZeroOnThree points positions rows} ≃
      (Fin 8 → F) := by
  let e := rs8_four_eval_equiv points hpoints positions hpositions
  let fourth : Fin 4 := ⟨3, by omega⟩
  refine ⟨(fun x => e x.val fourth),
    (fun v => ⟨e.symm (extendFourth v), ?_⟩), ?_, ?_⟩
  · intro i
    have hzero : e (e.symm (extendFourth v)) ⟨i.val, by omega⟩ = 0 := by
      rw [e.apply_symm_apply]
      have hi : i.val ≠ 3 := by omega
      simp [extendFourth, hi]
    funext r
    exact (rs8_four_eval_equiv_apply points hpoints positions
      hpositions _ _ r).symm.trans (congrFun hzero r)
  · intro x
    apply Subtype.ext
    apply e.injective
    funext i
    by_cases hi : i.val = 3
    · have hieq : i = fourth := Fin.ext hi
      subst i
      change e (e.symm (extendFourth (e x.val fourth))) fourth =
        e x.val fourth
      simp [extendFourth, fourth]
    · have hlt : i.val < 3 := by omega
      let j : Fin 3 := ⟨i.val, hlt⟩
      have hzero := x.property j
      have hieq : i = ⟨j.val, by omega⟩ := Fin.ext rfl
      have hzero' : e x.val i = 0 := by
        rw [hieq]
        funext r
        exact (rs8_four_eval_equiv_apply points hpoints positions
          hpositions x.val _ r).trans (congrFun hzero r)
      change e (e.symm (extendFourth (e x.val fourth))) i = e x.val i
      rw [e.apply_symm_apply]
      simp [extendFourth, hi, hzero']
  · intro v
    change e (e.symm (extendFourth v)) fourth = v
    simp [extendFourth, fourth]

/-- Exactly `|F|^8` group messages vanish in any fixed three columns. -/
theorem rs8ThreeZero_card
    [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    Fintype.card {rows : RS8Group F // rs8ZeroOnThree points positions rows} =
      (Fintype.card F) ^ 8 := by
  classical
  rw [Fintype.card_congr (rs8ThreeZeroEquiv points hpoints positions hpositions)]
  simp

/-- Removing the all-zero group leaves `|F|^8-1` nonzero messages. -/
theorem rs8ThreeZero_nonzero_card
    [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    (Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8ZeroOnThree points positions rows} =>
        x.val ≠ 0)).card = (Fintype.card F) ^ 8 - 1 := by
  classical
  let α := {rows : RS8Group F // rs8ZeroOnThree points positions rows}
  let z : α := ⟨0, by
    intro i
    funext r
    simp [rs8Column]⟩
  have hzero :
      (Finset.univ.filter (fun x : α => x.val = 0)) = {z} := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_singleton]
    exact ⟨fun hx => Subtype.ext hx, fun hx => by rw [hx]⟩
  have hpart := Finset.card_filter_add_card_filter_not
    (s := (Finset.univ : Finset α)) (p := fun x : α => x.val = 0)
  rw [hzero] at hpart
  simp only [Finset.card_singleton, Finset.card_univ] at hpart
  have hcard := rs8ThreeZero_card points hpoints positions hpositions
  change 1 + (Finset.univ.filter (fun x : α => x.val ≠ 0)).card =
    Fintype.card α at hpart
  change Fintype.card α = Fintype.card F ^ 8 at hcard
  change (Finset.univ.filter (fun x : α => x.val ≠ 0)).card =
    Fintype.card F ^ 8 - 1
  omega

omit [Fintype F] in
/-- The nonzero messages of this shortened code have exactly five active
columns, so the previous count is an exact five-column support count. -/
theorem rs8ThreeZero_nonzero_support_five
    [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions)
    (rows : RS8Group F) (hzero : rs8ZeroOnThree points positions rows)
    (hrows : rows ≠ 0) :
    (rs8Support points rows).card = 5 := by
  classical
  let zpos : Finset (Fin 8) :=
    Finset.univ.image (fun i : Fin 3 => positions ⟨i.val, by omega⟩)
  have hzpos : zpos.card = 3 := by
    rw [Finset.card_image_of_injective]
    · simp
    · intro i j hij
      exact Fin.ext (Fin.mk.inj (hpositions hij))
  let zeros : Finset (Fin 8) :=
    Finset.univ.filter (fun j => rs8Column points rows j = 0)
  have hsubset : zpos ⊆ zeros := by
    intro j hj
    obtain ⟨i, _, rfl⟩ := Finset.mem_image.mp hj
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hzero i⟩
  have hz : 3 ≤ zeros.card := by
    simpa [hzpos] using Finset.card_le_card hsubset
  have hpart : zeros.card + (rs8Support points rows).card = 8 := by
    simpa [zeros, rs8Support] using
      (Finset.card_filter_add_card_filter_not
        (s := (Finset.univ : Finset (Fin 8)))
        (p := fun j => rs8Column points rows j = 0))
  have hfloor := rs8Support_ge_five points hpoints rows hrows
  omega

#print axioms rs8ThreeZeroEquiv
#print axioms rs8ThreeZero_card
#print axioms rs8ThreeZero_nonzero_card
#print axioms rs8ThreeZero_nonzero_support_five

end Spin.Turbo
