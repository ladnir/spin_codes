import SpinCodes.Turbo.RS8SpectrumCount

/-!
# The RS8 support spectrum at sixteen field elements

The preceding algebraic proof reduces each exact-support count to a
finite sum over subsets of eight columns. This module checks that this
small sum equals the paper's MDS `A_h` coefficient for every possible
support. The computation involves only eight-element finite sets and
integers; it does not enumerate the `2^128` group messages.
-/

namespace Spin.Turbo

set_option maxRecDepth 4096

private def supportCountSum (s : Finset (Fin 8)) : ℤ :=
  ∑ t ∈ s.powerset, (-1 : ℤ) ^ t.card *
    ((if (sᶜ ∪ t).card ≤ 4 then
        16 ^ (8 * (4-(sᶜ ∪ t).card))
      else 1 : ℕ) : ℤ)

private def supportCountClosed (h : ℕ) : ℤ :=
  if h = 0 then 1 else
    ∑ j ∈ Finset.range (h-4), (-1 : ℤ) ^ j * (h.choose j : ℤ) *
      (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ)

/-- The exact-support expression depends only on support size. The empty
support contains the zero message, and sizes five through eight equal the
closed RS MDS count. -/
theorem supportCountSum_eq_closed (s : Finset (Fin 8)) :
    supportCountSum s = supportCountClosed s.card := by
  have hall : ∀ u : Finset (Fin 8),
      supportCountSum u = supportCountClosed u.card := by decide
  exact hall s

/-- The count of group messages with a specified support is the paper's
MDS count per support, provided the base field has sixteen elements. -/
theorem rs8_exact_support_card_closed
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF : Fintype.card F = 16)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (s : Finset (Fin 8)) :
    ((Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card : ℤ) =
      supportCountClosed s.card := by
  rw [rs8_exact_support_formula points hpoints]
  simp only [hF]
  exact supportCountSum_eq_closed s

/-- For support size five or more, the zero-message exception disappears
and the paper's alternating binomial expression is exact. -/
theorem rs8_exact_support_card_MDS
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF : Fintype.card F = 16)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (s : Finset (Fin 8)) (hs : 5 ≤ s.card) :
    ((Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card : ℤ) =
      ∑ j ∈ Finset.range (s.card-4),
        (-1 : ℤ) ^ j * (s.card.choose j : ℤ) *
          (((2^32 : ℕ) ^ (s.card-4-j) - 1 : ℕ) : ℤ) := by
  rw [rs8_exact_support_card_closed hF points hpoints s]
  have hne : s.card ≠ 0 := by omega
  simp [supportCountClosed, hne]

/-- Sum the exact-support counts over all supports of size `h`. -/
noncomputable def rs8ColumnWeightCount
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (points : Fin 8 → F) (h : ℕ) : ℤ :=
  ∑ s ∈ (Finset.univ : Finset (Fin 8)).powersetCard h,
    ((Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card : ℤ)

private theorem rs8_local_exact_eq_full
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (points : Fin 8 → F) (s : Finset (Fin 8)) :
    (Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card =
    (Finset.univ.filter (fun rows : RS8Group F =>
      rs8Support points rows = s)).card := by
  classical
  let α := {rows : RS8Group F // rs8SupportedIn points s rows}
  let localSet : Finset α :=
    Finset.univ.filter (fun x => rs8Support points x.val = s)
  let full : Finset (RS8Group F) :=
    Finset.univ.filter (fun rows => rs8Support points rows = s)
  have himage : localSet.image Subtype.val = full := by
    ext rows
    constructor
    · intro hr
      obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp hr
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _,
        (Finset.mem_filter.mp hx).2⟩
    · intro hr
      have hs := (Finset.mem_filter.mp hr).2
      have hsub : rs8SupportedIn points s rows := by
        unfold rs8SupportedIn
        rw [hs]
      exact Finset.mem_image.mpr
        ⟨⟨rows, hsub⟩, Finset.mem_filter.mpr ⟨Finset.mem_univ _, hs⟩, rfl⟩
  calc
    localSet.card = (localSet.image Subtype.val).card :=
      (Finset.card_image_of_injective localSet Subtype.val_injective).symm
    _ = full.card := by rw [himage]

/-- The support-class sum is the actual number of RS8 group messages of
column weight `h`, including the zero message when `h=0`. -/
theorem rs8ColumnWeightCount_eq_card
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (points : Fin 8 → F) (h : ℕ) :
    rs8ColumnWeightCount points h =
      ((Finset.univ.filter (fun rows : RS8Group F =>
        (rs8Support points rows).card = h)).card : ℤ) := by
  classical
  let supports := (Finset.univ : Finset (Fin 8)).powersetCard h
  calc
    rs8ColumnWeightCount points h =
        ∑ s ∈ supports,
          ((Finset.univ.filter (fun rows : RS8Group F =>
            rs8Support points rows = s)).card : ℤ) := by
      unfold rs8ColumnWeightCount
      apply Finset.sum_congr rfl
      intro s _
      exact congrArg (fun n : ℕ => (n : ℤ)) (rs8_local_exact_eq_full points s)
    _ = ((∑ s ∈ supports,
        (Finset.univ.filter (fun rows : RS8Group F =>
          rs8Support points rows = s)).card : ℕ) : ℤ) := by simp
    _ = ((Finset.univ.filter (fun rows : RS8Group F =>
        rs8Support points rows ∈ supports)).card : ℤ) := by
      congr 1
      exact Finset.sum_card_fiberwise_eq_card_filter
        (Finset.univ : Finset (RS8Group F)) supports
        (rs8Support points)
    _ = _ := by
      congr 1
      congr 1
      apply Finset.filter_congr
      intro rows _
      simp [supports, Finset.mem_powersetCard]

/-- The paper's `A_h` formula for the RS8 outer, as a count of messages
over all supports of a specified weight. -/
theorem rs8ColumnWeightCount_MDS
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF : Fintype.card F = 16)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (h : ℕ) (hh : 5 ≤ h) :
    rs8ColumnWeightCount points h =
      (Nat.choose 8 h : ℤ) *
        ∑ j ∈ Finset.range (h-4),
          (-1 : ℤ) ^ j * (h.choose j : ℤ) *
            (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ) := by
  unfold rs8ColumnWeightCount
  calc
    _ = ∑ s ∈ (Finset.univ : Finset (Fin 8)).powersetCard h,
        ∑ j ∈ Finset.range (h-4),
          (-1 : ℤ) ^ j * (h.choose j : ℤ) *
            (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ) := by
      apply Finset.sum_congr rfl
      intro s hs
      have hcard : s.card = h := (Finset.mem_powersetCard.mp hs).2
      simpa [hcard] using
        rs8_exact_support_card_MDS hF points hpoints s (by omega)
    _ = _ := by
      simp [Finset.card_powersetCard]

/-- The paper's `A_h` is the actual number of eight-row messages whose
outer word has exactly `h` active columns. -/
theorem rs8_actual_column_weight_count_MDS
    {F : Type*} [Field F] [Fintype F] [DecidableEq F]
    (hF : Fintype.card F = 16)
    (points : Fin 8 → F) (hpoints : Function.Injective points)
    (h : ℕ) (hh : 5 ≤ h) :
    ((Finset.univ.filter (fun rows : RS8Group F =>
      (rs8Support points rows).card = h)).card : ℤ) =
      (Nat.choose 8 h : ℤ) *
        ∑ j ∈ Finset.range (h-4),
          (-1 : ℤ) ^ j * (h.choose j : ℤ) *
            (((2^32 : ℕ) ^ (h-4-j) - 1 : ℕ) : ℤ) := by
  rw [← rs8ColumnWeightCount_eq_card]
  exact rs8ColumnWeightCount_MDS hF points hpoints h hh

#print axioms supportCountSum_eq_closed
#print axioms rs8_exact_support_card_MDS
#print axioms rs8ColumnWeightCount_MDS
#print axioms rs8ColumnWeightCount_eq_card
#print axioms rs8_actual_column_weight_count_MDS

end Spin.Turbo
