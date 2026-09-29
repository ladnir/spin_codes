/-
Composing two interleaved accumulator stages.

`lem:structured-exact-ba-spectrum` reasons: "Fix a Golay direct-sum word of
weight `a`.  The first uniform permutation and accumulator produce weight `c`
with probability `P_b(a,c)`.  Conditioned on that weight, the second
permutation makes the word uniform on its Hamming slice.  The second
accumulator then produces weight `w` with probability `P_b(c,w)`."

The middle step is what makes the composition work, and it is *not* that the
first stage's output is uniform on its slice — it is not.  What is true is
that the second permutation is independent of the first, so conditioned on the
first stage's output weight, the second stage sees a uniform slice element.
Formally: `prob_accWt` applies to *any* fixed word of weight `c`, so the inner
probability depends on the first stage only through `accWtF (u ∘ τ₁)`.

That is exactly what `prob_two_stage` below exploits.
-/
import SpinCodes.Prob
import SpinCodes.Structured.Interleaver

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-- `P_b(a,c) = T_b(a,c)/C(b,a)`, the weight-transition probability of a
uniformly interleaved accumulator. -/
noncomputable def Pt (b a c : ℕ) : ℝ := (accT b a c : ℝ) / (b.choose a : ℝ)

lemma accWtF_le {b : ℕ} (v : Fin b → Bool) : accWtF v ≤ b := by
  rw [← wtF_accF]
  exact wtF_le _

/-- **The interleaver law, in probability form.** -/
theorem prob_accWt {b a c : ℕ} (u : Fin b → Bool) (hu : wtF u = a) :
    (FinPMF.uniform (Equiv.Perm (Fin b))).prob (fun τ => accWtF (u ∘ τ) = c)
      = Pt b a c := by
  classical
  have hab : a ≤ b := hu ▸ wtF_le u
  have hC : (0:ℝ) < (b.choose a : ℝ) := by exact_mod_cast Nat.choose_pos hab
  have hN : (0:ℝ) < (Fintype.card (Equiv.Perm (Fin b)) : ℝ) := by
    exact_mod_cast Fintype.card_pos
  rw [FinPMF.uniform_prob, Pt, div_eq_div_iff (ne_of_gt hN) (ne_of_gt hC)]
  exact_mod_cast card_perm_accWt_eq u hu

/-- **Two interleaved accumulator stages compose by the Chapman–Kolmogorov
sum**: for independent uniform permutations `τ₁, τ₂`,

  `Pr[ wt(Acc(τ₂ · Acc(τ₁ · u))) = w ] = Σ_c P_b(a,c) · P_b(c,w)`,

where `a = wt(u)`.  This is the inner identity of
`eq:structured-exact-ba-spectrum`. -/
theorem prob_two_stage {b a w : ℕ} (u : Fin b → Bool) (hu : wtF u = a) :
    ((FinPMF.uniform (Equiv.Perm (Fin b))).prod
        (FinPMF.uniform (Equiv.Perm (Fin b)))).prob
        (fun τ => accWtF ((accF (u ∘ τ.1)) ∘ τ.2) = w)
      = ∑ c ∈ range (b + 1), Pt b a c * Pt b c w := by
  classical
  rw [FinPMF.prob_prod_eq_expect _ _
    (fun τ₁ τ₂ : Equiv.Perm (Fin b) => accWtF ((accF (u ∘ τ₁)) ∘ τ₂) = w)]
  have hinner : ∀ τ₁ : Equiv.Perm (Fin b),
      (FinPMF.uniform (Equiv.Perm (Fin b))).prob
          (fun τ₂ => accWtF ((accF (u ∘ τ₁)) ∘ τ₂) = w)
        = Pt b (accWtF (u ∘ τ₁)) w := by
    intro τ₁
    have := prob_accWt (c := w) (accF (u ∘ τ₁)) (wtF_accF (u ∘ τ₁))
    exact this
  simp only [hinner]
  rw [FinPMF.expect_comp _ (fun τ₁ : Equiv.Perm (Fin b) => accWtF (u ∘ τ₁))
    (range (b + 1)) (fun τ₁ => by
      simp only [Finset.mem_range, Nat.lt_succ_iff]
      exact accWtF_le _) (fun c => Pt b c w)]
  exact Finset.sum_congr rfl fun c _ => by rw [prob_accWt u hu]

/-- The transition probabilities out of a fixed weight sum to one. -/
theorem sum_Pt {b a : ℕ} (hab : a ≤ b) : ∑ c ∈ range (b + 1), Pt b a c = 1 := by
  classical
  obtain ⟨u, hu⟩ : ∃ u : Fin b → Bool, wtF u = a := by
    have hpos : 0 < sliceCount b a := by
      rw [sliceCount_eq_choose]
      exact Nat.choose_pos hab
    unfold sliceCount at hpos
    obtain ⟨u, hu⟩ := Finset.card_pos.mp hpos
    exact ⟨u, by simpa using hu⟩
  have hcover : ∀ τ : Equiv.Perm (Fin b), accWtF (u ∘ τ) ∈ range (b + 1) := by
    intro τ
    simp only [Finset.mem_range, Nat.lt_succ_iff]
    exact accWtF_le _
  have h1 : ∑ c ∈ range (b + 1),
      (FinPMF.uniform (Equiv.Perm (Fin b))).prob (fun τ => accWtF (u ∘ τ) = c) = 1 := by
    have := FinPMF.expect_comp (FinPMF.uniform (Equiv.Perm (Fin b)))
      (fun τ : Equiv.Perm (Fin b) => accWtF (u ∘ τ)) (range (b + 1)) hcover
      (fun _ => (1:ℝ))
    rw [FinPMF.expect_const] at this
    simpa using this.symm
  rw [← h1]
  exact Finset.sum_congr rfl fun c _ => (prob_accWt u hu).symm

end Spin.Structured
