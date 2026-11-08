import SpinCodes.LinearCircuit

/-!
# Chronological circuits with a carried state

The scan is a circuit built from composition, parallel execution, and wires.
It introduces no new cost primitive. Its semantic recurrence retains each
output and passes only the state to the next local step.
-/
noncomputable section
namespace Spin.LinearCircuit

variable {S X Y : Type} [Fintype S] [Fintype X] [Fintype Y]

theorem apply_sum_elim {A B C D : Type} (f : C → D) (g : A → C) (h : B → C) :
    (fun i => f (Sum.elim g h i)) =
      Sum.elim (fun a => f (g a)) (fun b => f (h b)) := by
  funext i
  cases i <;> rfl

def scanIn (R : ℕ) :
    Sum S (Fin (R + 1) × X) ≃ Sum (Sum S X) (Fin R × X) where
  toFun := Sum.elim (fun s => .inl (.inl s))
    (fun ix => Fin.cases (.inl (.inr ix.2)) (fun i => .inr (i, ix.2)) ix.1)
  invFun := Sum.elim (Sum.elim Sum.inl (fun x => .inr (0, x)))
    (fun ix => .inr (ix.1.succ, ix.2))
  left_inv := by
    rintro (s | ⟨i, x⟩)
    · rfl
    · refine Fin.cases ?_ (fun i => ?_) i <;> rfl
  right_inv := by
    rintro ((s | x) | ⟨i, x⟩) <;> rfl

def scanOut (R : ℕ) :
    Sum Y (Sum (Fin R × Y) S) ≃ Sum (Fin (R + 1) × Y) S where
  toFun := Sum.elim (fun y => .inl (0, y))
    (Sum.elim (fun iy => .inl (iy.1.succ, iy.2)) Sum.inr)
  invFun := Sum.elim
    (fun iy => Fin.cases (.inl iy.2) (fun i => .inr (.inl (i, iy.2))) iy.1)
    (fun s => .inr (.inr s))
  left_inv := by rintro (y | (⟨i, y⟩ | s)) <;> rfl
  right_inv := by
    rintro (⟨i, y⟩ | s)
    · refine Fin.cases ?_ (fun i => ?_) i <;> rfl
    · rfl

def scanNil : Sum S (Fin 0 × X) ≃ Sum (Fin 0 × Y) S where
  toFun := Sum.elim Sum.inr (fun ix => Fin.elim0 ix.1)
  invFun := Sum.elim (fun iy => Fin.elim0 iy.1) Sum.inl
  left_inv := by rintro (s | ⟨i, x⟩); rfl; exact Fin.elim0 i
  right_inv := by rintro (⟨i, y⟩ | s); exact Fin.elim0 i; rfl

def fullScan : {R : ℕ} → (Fin R → Circuit (Sum S X) (Sum Y S)) →
    Circuit (Sum S (Fin R × X)) (Sum (Fin R × Y) S)
  | 0, _ => .reindex scanNil
  | R + 1, cs => .comp (.reindex (scanOut R))
      (.comp (.parallel (.reindex (Equiv.refl Y)) (fullScan (Fin.tail cs)))
        (.comp (.reindex (Equiv.sumAssoc Y S (Fin R × X)))
          (.comp (.parallel (cs 0) (.reindex (Equiv.refl (Fin R × X))))
            (.reindex (scanIn R)))))

def runScan : {R : ℕ} → (Fin R → Circuit (Sum S X) (Sum Y S)) →
    (Fin R → X → Bit) → (S → Bit) → (Fin R → Y → Bit) × (S → Bit)
  | 0, _, _, q => (Fin.elim0, q)
  | _R + 1, cs, xs, q =>
      let v := (cs 0).eval (Sum.elim q (xs 0))
      let tail := runScan (Fin.tail cs) (Fin.tail xs) (fun s => v (.inr s))
      (Fin.cons (fun y => v (.inl y)) tail.1, tail.2)

theorem fullScan_size {R : ℕ} (cs : Fin R → Circuit (Sum S X) (Sum Y S)) :
    (fullScan cs).size = ∑ i, (cs i).size := by
  induction R with
  | zero => simp [fullScan, Circuit.size]
  | succ R ih =>
    simp only [fullScan, Circuit.size, Nat.zero_add, Nat.add_zero, ih,
      Fin.sum_univ_succ, Fin.tail_def]
    omega

theorem fullScan_eval {R : ℕ} (cs : Fin R → Circuit (Sum S X) (Sum Y S))
    (xs : Fin R → X → Bit) (q : S → Bit) :
    (fullScan cs).eval (Sum.elim q (fun ix => xs ix.1 ix.2)) =
      Sum.elim (fun iy => (runScan cs xs q).1 iy.1 iy.2) (runScan cs xs q).2 := by
  induction R generalizing q with
  | zero =>
    funext z
    rcases z with ⟨i, y⟩ | s
    · exact Fin.elim0 i
    · rfl
  | succ R ih =>
    funext z
    rcases z with ⟨i, y⟩ | s
    · refine Fin.cases ?_ (fun i => ?_) i
      · simp only [fullScan, Circuit.eval, scanIn, scanOut, Equiv.sumAssoc, runScan,
          Equiv.symm, Equiv.refl, Equiv.coe_fn_mk, Sum.elim_inl, Sum.elim_inr,
          Fin.cases_zero, Fin.cases_succ, Fin.cons_zero, Fin.cons_succ, Function.comp_apply, id_eq]
        congr 1
        funext a
        cases a <;> rfl
      · simpa only [fullScan, Circuit.eval, scanIn, scanOut, Equiv.sumAssoc, runScan,
          Equiv.symm, Equiv.refl, Equiv.coe_fn_mk, Sum.elim_inl, Sum.elim_inr,
          Fin.cases_zero, Fin.cases_succ, Fin.cons_zero, Fin.cons_succ, Function.comp_apply, id_eq,
          apply_sum_elim, Fin.tail_def] using
          congrFun (ih (Fin.tail cs) (Fin.tail xs)
            (fun s => (cs 0).eval (Sum.elim q (xs 0)) (.inr s))) (.inl (i, y))
    · simpa only [fullScan, Circuit.eval, scanIn, scanOut, Equiv.sumAssoc, runScan,
        Equiv.symm, Equiv.refl, Equiv.coe_fn_mk, Sum.elim_inl, Sum.elim_inr,
        Fin.cases_zero, Fin.cases_succ, Fin.cons_zero, Fin.cons_succ, Function.comp_apply, id_eq,
        apply_sum_elim, Fin.tail_def] using
        congrFun (ih (Fin.tail cs) (Fin.tail xs)
          (fun s => (cs 0).eval (Sum.elim q (xs 0)) (.inr s))) (.inr s)

def emptyLeft (A : Type) : A ≃ Sum Empty A where
  toFun := Sum.inr
  invFun := Sum.elim Empty.elim id
  left_inv := fun _ => rfl
  right_inv := by rintro (e | a); exact e.elim; rfl

def emptyRight (A : Type) : Sum A Empty ≃ A where
  toFun := Sum.elim id Empty.elim
  invFun := Sum.inl
  left_inv := by rintro (a | e); rfl; exact e.elim
  right_inv := fun _ => rfl

def zeroScan {R : ℕ} (cs : Fin R → Circuit (Sum S X) (Sum Y S)) :
    Circuit (Fin R × X) (Fin R × Y) :=
  .comp (.reindex (emptyRight _))
    (.comp (.parallel (.reindex (Equiv.refl _))
      (.dense (0 : Matrix Empty S Bit)))
      (.comp (fullScan cs)
        (.comp (.parallel (.dense (0 : Matrix S Empty Bit))
          (.reindex (Equiv.refl _))) (.reindex (emptyLeft _)))))

theorem zeroScan_size {R : ℕ} (cs : Fin R → Circuit (Sum S X) (Sum Y S)) :
    (zeroScan cs).size = ∑ i, (cs i).size := by
  simp [zeroScan, Circuit.size, fullScan_size]

theorem zeroScan_eval {R : ℕ} (cs : Fin R → Circuit (Sum S X) (Sum Y S))
    (xs : Fin R → X → Bit) :
    (zeroScan cs).eval (fun ix => xs ix.1 ix.2) =
      fun iy => (runScan cs xs 0).1 iy.1 iy.2 := by
  funext iy
  have h := congrFun (fullScan_eval cs xs 0) (.inl iy)
  simpa only [zeroScan, Circuit.eval, emptyLeft, emptyRight, Equiv.symm,
    Matrix.zero_mulVec, Equiv.refl, Equiv.coe_fn_mk, Sum.elim_inl, Sum.elim_inr, id_eq] using h

#print axioms fullScan_eval
#print axioms zeroScan_eval
#print axioms zeroScan_size

end Spin.LinearCircuit
