import SpinCodes.Structured.LowCancellationMoments

noncomputable section
namespace Spin.Structured.LowCancellation
open PackedMap ConcreteMaps

theorem sum_list_finset {α β : Type*} (xs : List α) (s : Finset β) (f : α → β → ℝ) :
    (∑ b ∈ s, (xs.map fun a => f a b).sum) =
      (xs.map fun a => ∑ b ∈ s, f a b).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih => simp [Finset.sum_add_distrib, ih]

theorem shell_exponents_sum (gs : List Group) (w : Nat) (z : ℝ) :
    ((shellExponents gs w).map (fun e => z ^ e)).sum =
      (gs.map fun g => if g.shell = w then groupMoment z g else 0).sum := by
  induction gs with
  | nil => rfl
  | cons g gs ih =>
    simp only [shellExponents, List.flatMap_cons, List.map_append, List.sum_append,
      List.map_cons, List.sum_cons] at *
    rw [ih]
    by_cases h : g.shell = w <;> simp [h, groupMoment, List.map_map, Function.comp_def]

theorem shell_cancellation_eq {j : Nat} {gs : List Group}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j) (z : ℝ) (i : Fin 5) :
    shellCancellationMoment j z i =
      ((shellExponents gs (shellWeight i)).map (fun e => z ^ e)).sum /
        Nat.choose 128 j / shellCount i := by
  unfold shellCancellationMoment
  simp_rw [cancellation_eq hv hn hl]
  rw [← Finset.sum_div, sum_list_finset, shell_exponents_sum]
  congr 2
  apply congrArg List.sum
  apply List.map_congr_left
  intro g hg
  rw [Finset.sum_ite_eq]
  have hw : support 19 g.syndrome ∈ weightShell i ↔ g.shell = shellWeight i := by
    simp only [weightShell, Finset.mem_filter, Finset.mem_univ, true_and,
      Aset_support (hv g hg).1, ← weight_eq_card_support, (hv g hg).2.2.1]
  simp only [hw]

theorem shell_cancellation_pattern {j : Nat} {gs : List Group} {p : List (Nat × Nat)}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j) (i : Fin 5)
    (hs : certSort (shellExponents gs (shellWeight i)) = expand p) (z : ℝ) :
    shellCancellationMoment j z i =
      Spin.Imt.Occupation.Sparse.pattern j p z / shellCount i := by
  rw [shell_cancellation_eq hv hn hl]
  have hp : (shellExponents gs (shellWeight i)).Perm (expand p) := by
    rw [← hs]
    exact (certSort_perm _).symm
  rw [(hp.map (fun e => z ^ e)).sum_eq, expand_sum]
  rfl

end Spin.Structured.LowCancellation
