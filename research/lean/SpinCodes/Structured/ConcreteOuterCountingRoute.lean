import SpinCodes.Structured.ConcreteOuterCountingProduct

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset ConcreteMarked

/-- The finite counting comparison survives the entire actual structured route. -/
theorem activeMessages_route_expect_le {L k : ℕ} (seed : Seed k) (S : Finset (Fin L)) {B : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (F : (Fin (k * 24) → Finset (Fin L)) → ℝ) (hF : ∀ regions, 0 ≤ F regions) :
    (∑ x ∈ activeMessages S, (ConcreteRoute.law (rowSupports seed x)).expect F) ≤
      ((2 : ℝ)^(k * 24) * B)^S.card * (regionsLaw S (k * 24)).expect F := by
  let G : (Fin L → Finset (Fin (k * 24))) → ℝ :=
    fun rows => (ConcreteRoute.regionKernel (ConcreteRoute.transpose rows)).expect F
  have hG : ∀ rows, 0 ≤ G rows := fun rows => FinPMF.expect_nonneg _ hF
  have hr (rows : Fin L → Finset (Fin (k * 24))) :
      (ConcreteRoute.law rows).expect F = (ConcreteRoute.rowLaw rows).expect G := by
    rw [ConcreteRoute.law, FinPMF.expect_bind, FinPMF.expect_map]
  have hfair : (fairRows S (k * 24)).expect G = (regionsLaw S (k * 24)).expect F := by
    calc (fairRows S (k * 24)).expect G =
        ((fairRows S (k * 24)).bind ConcreteRoute.rowLaw).expect G := by rw [fairRows_rowShuffle]
      _ = (fairRows S (k * 24)).expect (fun rows => (ConcreteRoute.rowLaw rows).expect G) :=
        FinPMF.expect_bind _ _ _
      _ = (fairRows S (k * 24)).expect (fun rows => (ConcreteRoute.law rows).expect F) := by
        simp only [hr]
      _ = ((fairRows S (k * 24)).bind ConcreteRoute.law).expect F :=
        (FinPMF.expect_bind _ _ _).symm
      _ = (regionsLaw S (k * 24)).expect F := by rw [fairRows_route]
  simp_rw [hr]
  rw [← hfair]
  exact activeMessages_expect_le seed S hB G hG

/-- No realized spectrum bound is assumed when the exact finite maximum is retained. -/
theorem activeMessages_route_expect_le_max {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (F : (Fin (k * 24) → Finset (Fin L)) → ℝ) (hF : ∀ regions, 0 ≤ F regions) :
    (∑ x ∈ activeMessages S, (ConcreteRoute.law (rowSupports seed x)).expect F) ≤
      ((2 : ℝ)^(k * 24) * maxShellRatio seed)^S.card * (regionsLaw S (k * 24)).expect F :=
  activeMessages_route_expect_le seed S (spectrum_le_max seed) F hF

/-- Selection of the shared outer gives the comparison with its expected spectrum. -/
theorem good_activeMessages_route_expect_le {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (W : Finset ℕ) {B : ℝ} (hB0 : 0 ≤ B)
    (hgood : Spin.Good (seedLaw k) spectrum (k * 24) W seed)
    (hB : ∀ w ∈ W, Spin.Abar (seedLaw k) spectrum w ≤ B * ((k * 24).choose w : ℝ))
    (F : (Fin (k * 24) → Finset (Fin L)) → ℝ) (hF : ∀ regions, 0 ≤ F regions) :
    (∑ x ∈ activeMessages S, (ConcreteRoute.law (rowSupports seed x)).expect F) ≤
      ((2 : ℝ)^(k * 24) * (((k * 24 : ℕ) : ℝ)^2 * B))^S.card *
        (regionsLaw S (k * 24)).expect F :=
  activeMessages_route_expect_le seed S (good_spectrum_bound seed W hB0 hgood hB) F hF

end Spin.Structured.ConcreteOuter
