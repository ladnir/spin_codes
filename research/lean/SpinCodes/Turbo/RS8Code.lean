import SpinCodes.Turbo.RSRowCode

/-!
# Eight parallel RS[8,4] rows

One group has eight degree-below-four rows. At a fixed evaluation
position their eight symbols form one column. A nonzero group has a
nonzero row, hence at least five nonzero columns. This is the exact
outer support floor used by the finite Turbo certificate.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F]

abbrev RS8Group (F : Type*) [Field F] := Fin 8 → RSRow F

def rs8Column (points : Fin 8 → F) (rows : RS8Group F)
    (j : Fin 8) : Fin 8 → F :=
  fun r => (rows r).1.eval (points j)

def rs8Support [DecidableEq F]
    (points : Fin 8 → F) (rows : RS8Group F) : Finset (Fin 8) :=
  Finset.univ.filter (fun j => rs8Column points rows j ≠ 0)

/-- Eight row messages together contain 32 symbols over the base field. -/
theorem rs8Group_card [Fintype F] :
    Fintype.card (RS8Group F) = (Fintype.card F) ^ 32 := by
  simp [RS8Group, rsRow_card, ← pow_mul]

/-- Every nonzero eight-row group has at least five active columns. -/
theorem rs8Support_ge_five [DecidableEq F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (rows : RS8Group F) (hrows : rows ≠ 0) :
    5 ≤ (rs8Support points rows).card := by
  classical
  have hex : ∃ r : Fin 8, rows r ≠ 0 := by
    by_contra h
    apply hrows
    funext r
    by_contra hr
    exact h ⟨r, hr⟩
  obtain ⟨r, hr⟩ := hex
  let rowSupport : Finset (Fin 8) :=
    Finset.univ.filter (fun j => (rows r).1.eval (points j) ≠ 0)
  have hsubset : rowSupport ⊆ rs8Support points rows := by
    intro j hj
    have hj' : (rows r).1.eval (points j) ≠ 0 :=
      (Finset.mem_filter.mp hj).2
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    intro hcol
    apply hj'
    have h := congrFun hcol r
    simpa [rs8Column] using h
  have hrow : 5 ≤ rowSupport.card :=
    rsRow_support_ge_five points hpoints (rows r) hr
  exact hrow.trans (Finset.card_le_card hsubset)

/-- Four distinct columns determine the complete eight-row message. -/
theorem rs8Group_eq_of_four_columns
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions)
    (rows rows' : RS8Group F)
    (heval : ∀ i : Fin 4,
      rs8Column points rows (positions i) =
      rs8Column points rows' (positions i)) : rows = rows' := by
  funext r
  apply rsRow_eq_of_four_evals_code (fun i => points (positions i))
    (hpoints.comp hpositions)
  intro i
  have hi := congrFun (heval i) r
  simpa [rs8Column] using hi

/-- In particular, the outer group encoder is injective. -/
theorem rs8Columns_injective
    (points : Fin 8 → F) (hpoints : Function.Injective points) :
    Function.Injective (fun rows : RS8Group F => rs8Column points rows) := by
  intro rows rows' heval
  apply rs8Group_eq_of_four_columns points hpoints
    (fun i : Fin 4 => ⟨i.val, by omega⟩)
  · intro i j hij
    exact Fin.ext (Fin.mk.inj hij)
  · intro i
    exact congrFun heval _

/-- Any four columns form a complete information set for all eight rows. -/
noncomputable def rs8_four_eval_equiv [Fintype F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions) :
    RS8Group F ≃ (Fin 4 → Fin 8 → F) :=
  (Equiv.piCongrRight (fun _ : Fin 8 =>
    rsRow_four_eval_equiv (fun i => points (positions i))
      (hpoints.comp hpositions))).trans
    (Equiv.piComm (fun (_ : Fin 8) (_ : Fin 4) => F))

theorem rs8_four_eval_equiv_apply [Fintype F]
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (positions : Fin 4 → Fin 8) (hpositions : Function.Injective positions)
    (rows : RS8Group F) (i : Fin 4) (r : Fin 8) :
    rs8_four_eval_equiv points hpoints positions hpositions rows i r =
      rs8Column points rows (positions i) r := by
  rfl

#print axioms rs8Group_card
#print axioms rs8Support_ge_five
#print axioms rs8Group_eq_of_four_columns
#print axioms rs8Columns_injective
#print axioms rs8_four_eval_equiv
#print axioms rs8_four_eval_equiv_apply

end Spin.Turbo
