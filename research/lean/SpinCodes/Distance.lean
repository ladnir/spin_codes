/-
The bad event of the first-moment framework is exactly `Z_d ≥ 1`
(paper: proof of `thm:class-first-moment`).

"`E` is not injective, or `d_min(E(𝔽₂^K)) ≤ d`" is not literally a statement
about nonzero *messages*, so the reduction is worth recording: a nonzero
kernel message contributes output weight zero, and every nonzero image word
of weight at most `d` has a nonzero preimage.
-/
import SpinCodes.Framework

set_option linter.unusedSectionVars false

namespace Spin

open Finset

variable {M W : Type*} [Fintype M] [DecidableEq M] [AddGroup M] [AddGroup W]

/-- An additive encoder sends `0` to `0`. -/
lemma map_zero_of_sub (E : M → W) (hsub : ∀ x y, E (x - y) = E x - E y) : E 0 = 0 := by
  have := hsub 0 0
  simpa using this

/-- Non-injectivity or small minimum distance both produce a nonzero message
of small encoded weight. -/
theorem exists_nonzero_light (wt : W → ℕ) (E : M → W) (d : ℕ)
    (hsub : ∀ x y, E (x - y) = E x - E y) (hwt0 : wt 0 = 0)
    (h : ¬ Function.Injective E ∨ ∃ y ∈ Set.range E, y ≠ 0 ∧ wt y ≤ d) :
    ∃ x : M, x ≠ 0 ∧ wt (E x) ≤ d := by
  rcases h with hinj | ⟨y, ⟨x, hx⟩, hy0, hyd⟩
  · -- a nonzero kernel message has encoded weight zero
    rw [Function.not_injective_iff] at hinj
    obtain ⟨a, b, hab, hne⟩ := hinj
    refine ⟨a - b, sub_ne_zero.mpr hne, ?_⟩
    rw [hsub, hab, sub_self, hwt0]
    exact Nat.zero_le d
  · -- a light nonzero codeword has a nonzero preimage
    refine ⟨x, ?_, ?_⟩
    · rintro rfl
      exact hy0 (by rw [← hx, map_zero_of_sub E hsub])
    · rw [hx]; exact hyd

/-- `Z_d ≥ 1` is exactly the existence of a light nonzero message. -/
lemma one_le_Z_iff {Ωout Ωin : Type*} [Fintype Ωout] [Fintype Ωin]
    (S : Setup M W Ωout Ωin) (d : ℕ) (ω : Ωout × Ωin) :
    1 ≤ S.Z d ω ↔ ∃ x : M, x ≠ 0 ∧ S.innerWt ω.2 (S.outer ω.1 x) ≤ d := by
  classical
  unfold Setup.Z
  constructor
  · intro h
    obtain ⟨x, hx⟩ := Finset.card_pos.mp h
    simp only [Finset.mem_filter, nonzeroMsgs, Finset.mem_univ, true_and] at hx
    exact ⟨x, hx.1, hx.2⟩
  · rintro ⟨x, hx0, hxd⟩
    refine Finset.card_pos.mpr ⟨x, ?_⟩
    simp only [Finset.mem_filter, nonzeroMsgs, Finset.mem_univ, true_and]
    exact ⟨hx0, hxd⟩

end Spin
