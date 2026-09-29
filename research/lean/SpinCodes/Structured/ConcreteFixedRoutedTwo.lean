import SpinCodes.Structured.ConcreteFixedRoutedOne
import SpinCodes.Structured.ConcreteFixedTwoMoment

noncomputable section
namespace Spin.Structured.Placement
open Finset Filter ConcreteRoute ConcreteEncoder ConcreteMarked

/-- The checked two-mark continuum bound now applies to the actual routed experiment. -/
theorem actual_two_fair_moment :
    ∀ᶠ R : ℕ in atTop, ∀ b : ℕ, ∀ S : Finset (Fin (128*R)), S.card=2 →
      (fairExperimentLaw S b (b*R)).expect
        (fun ω => (Real.exp (-(4/(128*R))))^ConcreteRoutedEncoder.weight regionWiring ω.1 ω.2) ≤
          1000*(1011/4000:ℝ)^b := by
  filter_upwards [actual_two_mark_moment] with R hR
  intro b S hS
  let e : Fin 2 ↪ Fin (128*R) := (S.orderEmbOfFin hS).toEmbedding
  have he : univ.map e = S := by
    rw [Finset.map_eq_image]
    exact Finset.image_orderEmbOfFin_univ S hS
  rw [fair_regionWiring_moment]
  have hh := hR b (fun _ => e)
  simpa only [he,regionsLaw] using hh

theorem actual_two_fair_probability :
    ∀ᶠ R : ℕ in atTop, ∀ b : ℕ, ∀ S : Finset (Fin (128*R)), S.card=2 → ∀ d : ℕ,
      (fairExperimentLaw S b (b*R)).prob
        (fun ω => ConcreteRoutedEncoder.weight regionWiring ω.1 ω.2≤d) ≤
          (1000*(1011/4000:ℝ)^b)/(Real.exp (-(4/(128*R))))^d := by
  filter_upwards [actual_two_fair_moment] with R hR
  intro b S hS d
  apply FinPMF.prob_weight_le_of_moment _ _ d (Real.exp_pos _)
    (Real.exp_le_one_iff.mpr (by apply neg_nonpos.mpr; positivity))
  exact hR b S hS

theorem actual_two_fair_probability_regionMajor :
    ∀ᶠ R : ℕ in atTop, ∀ b T : ℕ, ∀ h : b*(128*R)=T*128,
      ∀ S : Finset (Fin (128*R)), S.card=2 → ∀ d : ℕ,
      (fairExperimentLaw S b T).prob
        (fun ω => ConcreteRoutedEncoder.weight (regionMajor h) ω.1 ω.2≤d) ≤
          (1000*(1011/4000:ℝ)^b)/(Real.exp (-(4/(128*R))))^d := by
  filter_upwards [actual_two_fair_moment] with R hR
  intro b T h S hS d
  apply FinPMF.prob_weight_le_of_moment _ _ d (Real.exp_pos _)
    (Real.exp_le_one_iff.mpr (by apply neg_nonpos.mpr; positivity))
  rw [fair_regionMajor_moment,←fair_regionWiring_moment]
  exact hR b S hS

#print axioms actual_two_fair_moment
#print axioms actual_two_fair_probability
#print axioms actual_two_fair_probability_regionMajor
end Spin.Structured.Placement
