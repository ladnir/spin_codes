import SpinCodes.Structured.DenseOccupationProfileCount
import SpinCodes.Structured.DenseOccupationOuterStrip
import SpinCodes.Structured.ConcreteOuterCountingProbability

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter ConcreteRoute

/-- Independent row shuffling makes the actual failure probability depend only on row weights. -/
theorem failureProbability_eq_of_card_eq {L b R : ℕ}
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ)
    (rows rows' : Fin L → Finset (Fin b))
    (h : ∀ i, (rows i).card = (rows' i).card) :
    failureProbability e rows d = failureProbability e rows' d := by
  have hr : rowLaw rows = rowLaw rows' := by
    apply FinPMF.ext
    funext out
    simp only [rowLaw, piPMF_apply]
    apply Finset.prod_congr rfl
    intro i _
    simp only [Routing.shuffleLaw_apply, h i]
  have hl : law rows = law rows' := by unfold law; rw [hr]
  simp only [failureProbability_eq_route, hl]

/-- Counting an actual fixed weight profile multiplies one routed probability by its exact size. -/
theorem profile_failure_eq {L k R : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (rows : Fin L → Finset (Fin (k*24)))
    (hz : ∀ i, i ∉ S → rows i = ∅)
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) :
    (∑ x ∈ profileMessages seed S (fun i => (rows i).card),
      failureProbability e (rowSupports seed x) d) =
    (∏ i ∈ S, (spectrum seed (rows i).card:ℝ))*failureProbability e rows d := by
  have he : ∀ x ∈ profileMessages seed S (fun i => (rows i).card),
      failureProbability e (rowSupports seed x) d = failureProbability e rows d := by
    intro x hx
    apply failureProbability_eq_of_card_eq
    intro i
    have hi := (Finset.mem_filter.mp hx).2 i
    by_cases hs : i ∈ S
    · simp only [hs, ite_true] at hi
      simpa only [rowSupports, support_card] using hi.2
    · simp only [hs, ite_false] at hi
      simp only [rowSupports, hi, encode_zero, hz i hs, support_card, Finset.card_empty]
      exact (wtF_eq_zero _).mpr rfl
  rw [Finset.sum_congr rfl he, Finset.sum_const, nsmul_eq_mul, profileMessages_card]

#print axioms failureProbability_eq_of_card_eq
#print axioms profile_failure_eq
theorem selected_profile_failure_rate {L k R d : ℕ} (hL : 0 < L) (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (S : Finset (Fin L)) (rows : Fin L → Finset (Fin (k*24)))
    (hz : ∀ i, i ∉ S → rows i = ∅)
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128))
    (hw : ∀ i ∈ S, (rows i).card ∈ W)
    (hr : ∀ i ∈ S, ((rows i).card:ℝ)/(k*24:ℕ) ∈ Set.Icc (13/125) (112/125))
    {α x : ℝ} (hmean : ∑ i ∈ S, ((rows i).card:ℝ)/(k*24:ℕ) = (S.card:ℝ)*x)
    (hαQ : (L:ℝ)*α = S.card)
    (hα : α ∈ Set.Icc (1/10000) (10031/320000))
    (hx : x ∈ Set.Icc (13/125) (18296026121/120000000000))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*(k*24))
    (hden : density rows = α*x) (hround : 128*R = L*(k*24))
    (hd : (d:ℝ) ≤ (11/100)*((L:ℝ)*(k*24:ℕ))) :
    (∑ x ∈ profileMessages seed S (fun i => (rows i).card),
      failureProbability e (rowSupports seed x) d) ≤
    outerCost (k*24)^S.card *
      (((((k*24:ℕ):ℝ)+1)^activeRows rows*((L:ℝ)+1)^(k*24))*568*
        Real.exp (-(4/10000000)*((L:ℝ)*(k*24:ℕ)))) := by
  rw [profile_failure_eq seed S rows hz e d]
  exact selected_profile_routed_rate hL hk seed W hg S rows e hw hr hmean hαQ hα hx hw0 hw1 hden hround hd

#print axioms selected_profile_failure_rate
end Spin.Structured.DenseOccupationFixed
