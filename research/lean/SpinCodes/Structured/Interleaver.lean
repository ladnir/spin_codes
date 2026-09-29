/-
The uniform interleaver.

The paper's `P_b(a,c) = T_b(a,c)/C(b,a)` is "the exact weight-transition
probability of a uniformly interleaved accumulator", and the spectrum proof
turns on one asserted sentence: "the second permutation makes the word uniform
on its Hamming slice".

That sentence is the content of this file.  For a fixed `u` of weight `a` and
a uniform permutation `τ`, the word `u ∘ τ` is uniform on the weight-`a`
slice.  The argument is orbit-stabilizer in concrete form: the fibers of
`τ ↦ u ∘ τ` all have the same size, because right-translating by a permutation
carrying one slice element to another is a bijection between fibers.
-/
import SpinCodes.Structured.AccTuple

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-! ## Weight as a support cardinality -/

lemma wtL_eq_sum_map (l : List Bool) :
    wtL l = (l.map (fun x => if x then 1 else 0)).sum := by
  induction l with
  | nil => simp [wtL]
  | cons x t ih => simp [wtL, ih]

lemma wtF_eq_sum {b : ℕ} (v : Fin b → Bool) :
    wtF v = ∑ i, (if v i then 1 else 0) := by
  unfold wtF
  rw [wtL_eq_sum_map, List.map_ofFn]
  simp [List.sum_ofFn]

lemma wtF_eq_card {b : ℕ} (v : Fin b → Bool) :
    wtF v = (univ.filter (fun i => v i = true)).card := by
  classical
  rw [wtF_eq_sum, Finset.card_filter]

/-- Weight is invariant under permuting positions. -/
lemma wtF_comp_perm {b : ℕ} (v : Fin b → Bool) (ρ : Equiv.Perm (Fin b)) :
    wtF (v ∘ ρ) = wtF v := by
  classical
  rw [wtF_eq_sum, wtF_eq_sum]
  exact Fintype.sum_equiv ρ _ _ (fun i => rfl)

/-! ## The slice, and its size -/

/-- The number of tuples of weight `a`. -/
def sliceCount (b a : ℕ) : ℕ := (univ.filter (fun v : Fin b → Bool => wtF v = a)).card

lemma sliceCount_succ_zero (b : ℕ) : sliceCount (b + 1) 0 = sliceCount b 0 := by
  unfold sliceCount
  rw [card_filter_succ_split]
  have h0 : ∀ v : Fin b → Bool, wtF (Fin.snoc v false) = 0 ↔ wtF v = 0 := by
    intro v; rw [wtF_snoc]; simp
  have h1 : ∀ v : Fin b → Bool, ¬ (wtF (Fin.snoc v true) = 0) := by
    intro v; rw [wtF_snoc]; simp
  rw [Finset.filter_congr (fun v _ => by simpa using h0 v),
    Finset.filter_false_of_mem (fun v _ => h1 v)]
  simp

lemma sliceCount_succ (b a : ℕ) :
    sliceCount (b + 1) (a + 1) = sliceCount b (a + 1) + sliceCount b a := by
  unfold sliceCount
  rw [card_filter_succ_split]
  congr 1
  · refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
    rw [wtF_snoc]; simp
  · refine congrArg Finset.card (Finset.filter_congr fun v _ => ?_)
    rw [wtF_snoc]; simp

/-- The weight-`a` slice has `C(b,a)` elements. -/
theorem sliceCount_eq_choose (b a : ℕ) : sliceCount b a = b.choose a := by
  induction b generalizing a with
  | zero =>
      have hzero : ∀ v : Fin 0 → Bool, wtF v = 0 := by
        intro v
        unfold wtF
        rw [show List.ofFn v = [] by simp]
        rfl
      unfold sliceCount
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [Finset.filter_true_of_mem (fun v _ => hzero v)]
        simp
      · rw [Finset.filter_false_of_mem (fun v _ => by rw [hzero v]; omega),
          Nat.choose_eq_zero_of_lt (by omega)]
        simp
  | succ n ih =>
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [sliceCount_succ_zero, ih]
        simp
      · obtain ⟨a', rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : a ≠ 0)
        rw [sliceCount_succ, ih, ih, Nat.choose_succ_succ]
        simp only [Nat.succ_eq_add_one]
        omega

/-! ## Transitivity on a slice -/

/-- Two tuples of the same weight are related by a permutation of positions. -/
lemma exists_perm_comp {b : ℕ} {v v' : Fin b → Bool} (h : wtF v = wtF v') :
    ∃ ρ : Equiv.Perm (Fin b), v ∘ ρ = v' := by
  classical
  have hcard : Fintype.card {i : Fin b // v' i = true}
      = Fintype.card {i : Fin b // v i = true} := by
    rw [Fintype.card_subtype, Fintype.card_subtype, ← wtF_eq_card, ← wtF_eq_card, h]
  have hccard : Fintype.card {i : Fin b // ¬ (v' i = true)}
      = Fintype.card {i : Fin b // ¬ (v i = true)} := by
    have e1 : Fintype.card {i : Fin b // v' i = true}
        + Fintype.card {i : Fin b // ¬ (v' i = true)} = b := by
      rw [← Fintype.card_sum]
      simpa using Fintype.card_congr (Equiv.sumCompl (fun i : Fin b => v' i = true))
    have e2 : Fintype.card {i : Fin b // v i = true}
        + Fintype.card {i : Fin b // ¬ (v i = true)} = b := by
      rw [← Fintype.card_sum]
      simpa using Fintype.card_congr (Equiv.sumCompl (fun i : Fin b => v i = true))
    omega
  let e1 : {i : Fin b // v' i = true} ≃ {i : Fin b // v i = true} :=
    Fintype.equivOfCardEq hcard
  let e2 : {i : Fin b // ¬ (v' i = true)} ≃ {i : Fin b // ¬ (v i = true)} :=
    Fintype.equivOfCardEq hccard
  refine ⟨((Equiv.sumCompl (fun i : Fin b => v' i = true)).symm.trans
    ((Equiv.sumCongr e1 e2).trans (Equiv.sumCompl (fun i : Fin b => v i = true)))), ?_⟩
  funext i
  by_cases hi : v' i = true
  · have : ((Equiv.sumCompl (fun i : Fin b => v' i = true)).symm i)
        = Sum.inl ⟨i, hi⟩ := by
      simp [Equiv.sumCompl, hi]
    simp only [Function.comp_apply, Equiv.trans_apply, this, Equiv.sumCongr_apply,
      Sum.map_inl, Equiv.sumCompl_apply_inl]
    exact ((e1 ⟨i, hi⟩).2).trans hi.symm
  · have : ((Equiv.sumCompl (fun i : Fin b => v' i = true)).symm i)
        = Sum.inr ⟨i, hi⟩ := by
      simp [Equiv.sumCompl, hi]
    simp only [Function.comp_apply, Equiv.trans_apply, this, Equiv.sumCongr_apply,
      Sum.map_inr, Equiv.sumCompl_apply_inr]
    have h2 := (e2 ⟨i, hi⟩).2
    simp only [Bool.not_eq_true] at h2 hi
    rw [h2, hi]


/-! ## Fibers of the interleaving map -/

lemma comp_perm_mul {b : ℕ} (u : Fin b → Bool) (τ ρ : Equiv.Perm (Fin b)) :
    u ∘ (τ * ρ) = (u ∘ τ) ∘ ρ := rfl

/-- All fibers of `τ ↦ u ∘ τ` over a single weight slice have the same size.
Right-translating by a permutation that carries one slice element to another
is a bijection between the corresponding fibers. -/
lemma card_fiber_eq {b : ℕ} (u : Fin b → Bool) {v v' : Fin b → Bool}
    (h : wtF v = wtF v') :
    (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v)).card
      = (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v')).card := by
  classical
  obtain ⟨ρ, hρ⟩ := exists_perm_comp h
  have hmap : (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v'))
      = (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v)).map
          (Equiv.mulRight ρ).toEmbedding := by
    ext σ
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_map,
      Equiv.coe_toEmbedding, Equiv.coe_mulRight]
    constructor
    · intro hσ
      refine ⟨σ * ρ⁻¹, ?_, by group⟩
      rw [comp_perm_mul, hσ, ← hρ]
      funext i
      simp
    · rintro ⟨τ, hτ, rfl⟩
      rw [comp_perm_mul, hτ, hρ]
  rw [hmap, Finset.card_map]

/-! ## The transition law for a uniform interleaver -/

/-- **The uniform interleaver.**  For a fixed `u` of weight `a` and a uniformly
random permutation of positions, the accumulated weight equals `c` with
probability `T_b(a,c)/C(b,a)`.  Stated multiplicatively to stay in `ℕ`. -/
theorem card_perm_accWt_eq {b a c : ℕ} (u : Fin b → Bool) (hu : wtF u = a) :
    (univ.filter (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ) = c)).card * b.choose a
      = accT b a c * Fintype.card (Equiv.Perm (Fin b)) := by
  classical
  set K := (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = u)).card with hK
  -- every fiber over the weight-`a` slice has size `K`
  have hfib : ∀ v : Fin b → Bool, wtF v = a →
      (univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v)).card = K := by
    intro v hv
    rw [hK]
    exact card_fiber_eq u (by rw [hv, hu])
  -- the whole permutation group fibers over the slice
  have htot : Fintype.card (Equiv.Perm (Fin b)) = K * b.choose a := by
    have hmaps : ∀ τ ∈ (univ : Finset (Equiv.Perm (Fin b))),
        (u ∘ τ) ∈ univ.filter (fun v : Fin b → Bool => wtF v = a) := by
      intro τ _
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      rw [wtF_comp_perm, hu]
    have := Finset.card_eq_sum_card_fiberwise hmaps
    rw [Finset.card_univ] at this
    rw [this, Finset.sum_congr rfl (fun v hv => by
      simp only [Finset.mem_filter] at hv
      exact hfib v hv.2)]
    rw [Finset.sum_const, smul_eq_mul, mul_comm]
    congr 1
    exact sliceCount_eq_choose b a
  -- the event fibers over the slice intersected with the accumulated-weight level
  have hev : (univ.filter (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ) = c)).card
      = K * accT b a c := by
    have hmaps : ∀ τ ∈ univ.filter (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ) = c),
        (u ∘ τ) ∈ univ.filter (fun v : Fin b → Bool => wtF v = a ∧ accWtF v = c) := by
      intro τ hτ
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hτ ⊢
      exact ⟨by rw [wtF_comp_perm, hu], hτ⟩
    rw [Finset.card_eq_sum_card_fiberwise hmaps]
    have hinner : ∀ v ∈ univ.filter (fun v : Fin b → Bool => wtF v = a ∧ accWtF v = c),
        ((univ.filter (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ) = c)).filter
            (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v)).card = K := by
      intro v hv
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hv
      have hset : ((univ.filter (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ) = c)).filter
          (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v))
          = univ.filter (fun τ : Equiv.Perm (Fin b) => u ∘ τ = v) := by
        ext τ
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        constructor
        · rintro ⟨-, h2⟩; exact h2
        · intro h2; exact ⟨by rw [h2]; exact hv.2, h2⟩
      rw [hset]
      exact hfib v hv.1
    rw [Finset.sum_congr rfl hinner, Finset.sum_const, smul_eq_mul, mul_comm]
    congr 1
    exact accCountF_eq_accT b a c
  rw [hev, htot]
  ring

end Spin.Structured
