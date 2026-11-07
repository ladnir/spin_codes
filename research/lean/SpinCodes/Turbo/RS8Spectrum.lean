import SpinCodes.Turbo.RS8ShorteningGeneral
import Mathlib.Combinatorics.Enumerative.InclusionExclusion

/-!
# Exact support counting interface for Turbo RS8

The completed Structured SPIN proof counted outer words by weight before
applying its first-moment kernel. Turbo needs the analogous column-support
counts. This file first derives the inclusion-exclusion identity for an
arbitrary set of active columns. The shortened-code theorem supplies the
intersection counts in that identity.
-/

noncomputable section

namespace Spin.Turbo

variable {F : Type*} [Field F] [Fintype F] [DecidableEq F]

def rs8SupportedIn (points : Fin 8 → F) (s : Finset (Fin 8))
    (rows : RS8Group F) : Prop := rs8Support points rows ⊆ s

noncomputable instance (points : Fin 8 → F) (s : Finset (Fin 8)) :
    Fintype {rows : RS8Group F // rs8SupportedIn points s rows} :=
  Fintype.ofFinite _

def rs8MissingAt (points : Fin 8 → F) (s : Finset (Fin 8))
    (j : Fin 8) :
    Finset {rows : RS8Group F // rs8SupportedIn points s rows} :=
  Finset.univ.filter (fun x => rs8Column points x.val j = 0)

/-- Inclusion-exclusion for the messages with exactly a specified column
support. The right-hand intersections still need the shortening count. -/
theorem rs8_exact_support_inclusion_exclusion
    (points : Fin 8 → F) (s : Finset (Fin 8)) :
    ((Finset.univ.filter (fun x :
      {rows : RS8Group F // rs8SupportedIn points s rows} =>
      rs8Support points x.val = s)).card : ℤ) =
    ∑ t ∈ s.powerset, (-1 : ℤ) ^ t.card *
      (t.inf (rs8MissingAt points s)).card := by
  classical
  let α := {rows : RS8Group F // rs8SupportedIn points s rows}
  let bad : Fin 8 → Finset α := rs8MissingAt points s
  have heq :
      s.inf (fun j => (bad j)ᶜ) =
        Finset.univ.filter (fun x : α => rs8Support points x.val = s) := by
    ext x
    simp only [Finset.mem_inf, Finset.mem_compl, Finset.mem_filter,
      Finset.mem_univ, true_and]
    constructor
    · intro hx
      apply Finset.Subset.antisymm x.property
      intro j hj
      have hjnot := hx j hj
      have hjnonzero : rs8Column points x.val j ≠ 0 := by
        simpa [bad, rs8MissingAt] using hjnot
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, hjnonzero⟩
    · intro hx j hj
      have hmem : j ∈ rs8Support points x.val := hx.symm ▸ hj
      have hjnonzero : rs8Column points x.val j ≠ 0 :=
        (Finset.mem_filter.mp hmem).2
      simpa [bad, rs8MissingAt] using hjnonzero
  have hie := Finset.inclusion_exclusion_card_inf_compl s bad
  rw [heq] at hie
  simpa [bad] using hie

#print axioms rs8_exact_support_inclusion_exclusion

end Spin.Turbo
