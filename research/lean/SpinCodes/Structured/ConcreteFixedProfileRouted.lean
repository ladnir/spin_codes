import SpinCodes.Structured.ConcreteFixedRowEmbedding

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteRoute ConcreteEncoder Filter

theorem route_embed_expect {Q L b : ℕ} (e : Fin Q ↪ Fin L)
    (rows : Fin Q → Finset (Fin b)) (F : (Fin b → Finset (Fin L))→ℝ) :
    (law (embedRows e rows)).expect F = (rowLaw rows).expect (fun shuffled =>
      (Spin.piPMF (fun j => shuffleLaw ((transpose shuffled j).map e))).expect F) := by
  rw [law,FinPMF.expect_bind,FinPMF.expect_map,rowLaw_embed,FinPMF.expect_map]
  simp only [embedRows,transpose_transpose,regionKernel]

theorem actual_moment_route {L b T : ℕ}
    (e : (Fin b × Fin L) ≃ (Fin T × Fin 128)) (rows : Fin L → Finset (Fin b)) (z : ℝ) :
    (ConcreteRoutedEncoder.experimentLaw L b T).expect
      (fun ω => z^ConcreteRoutedEncoder.weight e rows ω) =
    (law rows).expect (fun regions => inputMoment z (reshape e regions) ∅) := by
  rw [ConcreteRoutedEncoder.moment_eq]
  rw [←FinPMF.expect_map (seedLaw L b) (routeEval rows)
    (fun regions => inputMoment z (reshape e regions) ∅),permutation_law_eq]

/-- Prescribed active rows, full actual two-stage route, and exact region-major IMT serialization. -/
theorem actual_profile_routed {Q : ℕ} {θ v c δ : ℝ} (hθ : 0≤θ) (hv : 0<v)
    (hc : 0≤c) (hδ : 0<δ) (u : Fin Q → ℝ) (hu : ∀ i,0<u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : ℕ in atTop, ∀ b T : ℕ, ∀ h : b*(128*R)=T*128,
      ∀ e : Fin Q ↪ Fin (128*R), ∀ rows : Fin Q → Finset (Fin b),
      (ConcreteRoutedEncoder.experimentLaw (128*R) b T).expect
        (fun ω => Real.exp (-(θ/(128*R)))^ConcreteRoutedEncoder.weight (regionMajor h) (embedRows e rows) ω) ≤
      ((c+δ)^b/min 1 v)/(profileChoices rows*profileCoefficient u (fun i => (rows i).card)) := by
  filter_upwards [actual_rows_regionStream_margin hθ hv hc hδ u hu hK] with R hR
  intro b T h e rows
  rw [actual_moment_route,route_embed_expect]
  simp_rw [inputMoment_regionMajor]
  exact hR b (fun _ => e) rows

theorem actual_profile_routed_probability {Q : ℕ} {θ v c δ : ℝ} (hθ : 0≤θ) (hv : 0<v)
    (hc : 0≤c) (hδ : 0<δ) (u : Fin Q → ℝ) (hu : ∀ i,0<u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : ℕ in atTop, ∀ b T : ℕ, ∀ h : b*(128*R)=T*128,
      ∀ e : Fin Q ↪ Fin (128*R), ∀ rows : Fin Q → Finset (Fin b), ∀ d : ℕ,
      (ConcreteRoutedEncoder.experimentLaw (128*R) b T).prob
        (fun ω => ConcreteRoutedEncoder.weight (regionMajor h) (embedRows e rows) ω≤d) ≤
      (((c+δ)^b/min 1 v)/(profileChoices rows*profileCoefficient u (fun i => (rows i).card)))/
        Real.exp (-(θ/(128*R)))^d := by
  filter_upwards [actual_profile_routed hθ hv hc hδ u hu hK] with R hR
  intro b T h e rows d
  apply FinPMF.prob_weight_le_of_moment _ _ d (Real.exp_pos _)
    (Real.exp_le_one_iff.mpr (by apply neg_nonpos.mpr; positivity))
  exact hR b T h e rows

#print axioms route_embed_expect
#print axioms actual_profile_routed
#print axioms actual_profile_routed_probability
end Spin.Structured.Placement
