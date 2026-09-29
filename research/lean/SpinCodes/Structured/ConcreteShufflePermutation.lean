import SpinCodes.Structured.ConcreteShuffleMixture

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.Routing

open Finset

/-- Relabel positions by a permutation, using the same pullback convention as the encoder. -/
def shuffleSupport {n : ℕ} (σ : Equiv.Perm (Fin n)) (S : Finset (Fin n)) :
    Finset (Fin n) := S.map σ.symm.toEmbedding

@[simp] lemma mem_shuffleSupport {n : ℕ} (σ : Equiv.Perm (Fin n))
    (S : Finset (Fin n)) (i : Fin n) : i ∈ shuffleSupport σ S ↔ σ i ∈ S := by
  simp [shuffleSupport]

@[simp] lemma shuffleSupport_card {n : ℕ} (σ : Equiv.Perm (Fin n))
    (S : Finset (Fin n)) : (shuffleSupport σ S).card = S.card := Finset.card_map _

def supportBool {n : ℕ} (S : Finset (Fin n)) (i : Fin n) : Bool := decide (i ∈ S)

lemma wtF_supportBool {n : ℕ} (S : Finset (Fin n)) : wtF (supportBool S) = S.card := by
  rw [wtF_eq_card]
  congr 1
  ext i
  simp [supportBool]

lemma supportBool_shuffle {n : ℕ} (S T : Finset (Fin n)) (σ : Equiv.Perm (Fin n)) :
    supportBool S ∘ σ = supportBool T ↔ shuffleSupport σ S = T := by
  classical
  constructor
  · intro h
    ext i
    have hi := congrFun h i
    simpa [supportBool, Function.comp_def] using hi
  · intro h
    funext i
    have hi := Finset.ext_iff.mp h i
    simp only [mem_shuffleSupport] at hi
    simpa only [supportBool, Function.comp_apply, decide_eq_decide] using hi

lemma shuffle_fiber_card_eq {n : ℕ} (S T U : Finset (Fin n)) (h : T.card = U.card) :
    (univ.filter (fun σ : Equiv.Perm (Fin n) => shuffleSupport σ S = T)).card =
      (univ.filter (fun σ : Equiv.Perm (Fin n) => shuffleSupport σ S = U)).card := by
  have hh := card_fiber_eq (supportBool S)
    (v := supportBool T) (v' := supportBool U) (by simpa only [wtF_supportBool] using h)
  simpa only [supportBool_shuffle] using hh

lemma shuffle_fiber_total {n : ℕ} (S T : Finset (Fin n)) (h : T.card = S.card) :
    Fintype.card (Equiv.Perm (Fin n)) =
      (univ.filter (fun σ : Equiv.Perm (Fin n) => shuffleSupport σ S = T)).card *
        n.choose S.card := by
  have maps : ∀ σ ∈ (univ : Finset (Equiv.Perm (Fin n))),
      shuffleSupport σ S ∈ univ.filter (fun U : Finset (Fin n) => U.card = S.card) := by
    intro σ _
    simp
  have hh := Finset.card_eq_sum_card_fiberwise maps
  rw [Finset.card_univ] at hh
  rw [hh]
  rw [Finset.sum_congr rfl (fun U hU =>
    shuffle_fiber_card_eq S U T ((Finset.mem_filter.mp hU).2.trans h.symm))]
  rw [Finset.sum_const, smul_eq_mul, card_layer_filter, mul_comm]

/-- The closed-form layer law is the pushforward of a genuinely uniform permutation. -/
theorem uniform_permutation_shuffleLaw {n : ℕ} (S : Finset (Fin n)) :
    (FinPMF.uniform (Equiv.Perm (Fin n))).map (fun σ => shuffleSupport σ S) =
      shuffleLaw S := by
  classical
  ext T
  have hp : ((FinPMF.uniform (Equiv.Perm (Fin n))).map
      (fun σ => shuffleSupport σ S)).p T =
      ((FinPMF.uniform (Equiv.Perm (Fin n))).map
      (fun σ => shuffleSupport σ S)).prob (fun U => U = T) := by
    simp [FinPMF.prob, Finset.sum_filter]
  rw [hp, FinPMF.map_prob, FinPMF.uniform_prob, shuffleLaw_apply]
  by_cases ht : T.card = S.card
  · rw [if_pos ht]
    have htotal := shuffle_fiber_total S T ht
    have htotalR : (Fintype.card (Equiv.Perm (Fin n)) : ℝ) =
        ((univ.filter (fun σ : Equiv.Perm (Fin n) => shuffleSupport σ S = T)).card : ℝ) *
          (n.choose S.card : ℝ) := by exact_mod_cast htotal
    have hpos : (0 : ℝ) < (Fintype.card (Equiv.Perm (Fin n)) : ℝ) := by
      exact_mod_cast Fintype.card_pos
    have hc := choose_card_pos S
    apply (div_eq_div_iff hpos.ne' hc.ne').mpr
    simpa using htotalR.symm
  · rw [if_neg ht]
    have he : univ.filter (fun σ : Equiv.Perm (Fin n) => shuffleSupport σ S = T) = ∅ := by
      apply Finset.filter_eq_empty_iff.mpr
      intro σ _ hs
      exact ht (by rw [← hs, shuffleSupport_card])
    rw [he]
    simp

/-- Shuffling a random support uses a fresh permutation independent of that support. -/
noncomputable def permutedLaw {n : ℕ} (P : FinPMF (Finset (Fin n))) :
    FinPMF (Finset (Fin n)) :=
  (P.prod (FinPMF.uniform (Equiv.Perm (Fin n)))).map
    (fun x => shuffleSupport x.2 x.1)

theorem permutedLaw_eq_shuffledLaw {n : ℕ} (P : FinPMF (Finset (Fin n))) :
    permutedLaw P = shuffledLaw P := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [permutedLaw, FinPMF.expect_map, FinPMF.expect_prod, shuffledLaw, FinPMF.expect_bind]
  unfold FinPMF.expect
  apply Finset.sum_congr rfl
  intro S _
  congr 1
  have hh := congrArg (fun Q : FinPMF (Finset (Fin n)) => Q.expect f)
    (uniform_permutation_shuffleLaw S)
  rw [FinPMF.expect_map] at hh
  exact hh

end Spin.Structured.Routing




