import SpinCodes.Structured.BAGolay

namespace Spin.Structured.ConcreteOuter
open Finset

/-- A tuple is zero exactly when its Hamming weight is zero. -/
theorem wtF_eq_zero {b : ℕ} (x : Fin b → Bool) : wtF x = 0 ↔ x = 0 := by
  rw [wtF_eq_card, Finset.card_eq_zero]
  constructor
  · intro h
    funext i
    have hi : x i ≠ true := by
      intro hx
      have : i ∈ (∅ : Finset (Fin b)) := h ▸ (by simp [hx] : i ∈ univ.filter (fun j => x j = true))
      exact Finset.notMem_empty i this
    exact Bool.eq_false_of_not_eq_true hi
  · intro h
    simp [h]

theorem accLAux_injective (p : Bool) : Function.Injective (accLAux p) := by
  intro xs
  induction xs generalizing p with
  | nil =>
      intro ys h
      cases ys <;> simp_all [accLAux]
  | cons x xs ih =>
      intro ys h
      cases ys with
      | nil => simp [accLAux] at h
      | cons y ys =>
          have hh := List.cons.inj h
          have hxy : x = y := by
            cases p <;> cases x <;> cases y <;> simp_all [accLAux]
          subst y
          exact congrArg (List.cons x) (ih (xor p x) hh.2)

theorem accF_injective {b : ℕ} : Function.Injective (@accF b) := by
  intro x y h
  have hh := congrArg List.ofFn h
  simp only [ofFn_accF] at hh
  have he := accLAux_injective false hh
  exact List.ofFn_injective he

theorem accLAux_false_replicate (n : ℕ) :
    accLAux false (List.replicate n false) = List.replicate n false := by
  induction n with
  | zero => rfl
  | succ n ih => simp [List.replicate_succ, accLAux, ih]

@[simp] theorem accF_zero (b : ℕ) : accF (0 : Fin b → Bool) = 0 := by
  apply List.ofFn_injective
  rw [ofFn_accF]
  change accL (List.ofFn (fun _ : Fin b => false)) = List.ofFn (fun _ : Fin b => false)
  simpa only [List.ofFn_const, accL] using accLAux_false_replicate b

@[simp] theorem accF_eq_zero {b : ℕ} (x : Fin b → Bool) : accF x = 0 ↔ x = 0 :=
  ⟨fun h => accF_injective (h.trans (accF_zero b).symm), fun h => h ▸ accF_zero b⟩

theorem ofFn_mem_allWords {n : ℕ} (x : Fin n → Bool) : List.ofFn x ∈ allWords n := by
  induction n with
  | zero => simp [allWords]
  | succ n ih =>
      rw [← Fin.snoc_init_self x, ofFn_snoc]
      simp only [allWords, List.mem_flatMap]
      refine ⟨List.ofFn (Fin.init x), ih _, ?_⟩
      cases x (Fin.last n) <;> simp

end Spin.Structured.ConcreteOuter

