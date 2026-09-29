import SpinCodes.Structured.ConcreteOuterCountingRoute
import SpinCodes.Structured.ConcreteMarkedReferenceMoment

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset ConcreteMarked

/-- The actual routed encoder's low-weight probability with its input rows fixed. -/
def failureProbability {L b R : ℕ}
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (d : ℕ) : ℝ :=
  (ConcreteRoutedEncoder.experimentLaw L b R).prob
    (fun ω => ConcreteRoutedEncoder.weight e rows ω ≤ d)

theorem failureProbability_eq_route {L b R : ℕ}
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (d : ℕ) :
    failureProbability e rows d = (ConcreteRoute.law rows).expect
      (fun regions => (piPMF (fun _ : Fin R => ConcreteEncoder.transvectionLaw)).prob
        (fun ts => ConcreteEncoder.outputWeight (ConcreteRoute.reshape e regions) ts ∅ ≤ d)) := by
  unfold failureProbability ConcreteRoutedEncoder.experimentLaw ConcreteRoutedEncoder.weight
  rw [FinPMF.prob_prod_eq_expect _ _ (fun seed ts =>
    ConcreteEncoder.outputWeight (ConcreteRoute.reshape e (ConcreteRoute.routeEval rows seed)) ts ∅ ≤ d)]
  rw [← FinPMF.expect_map (ConcreteRoute.seedLaw L b) (ConcreteRoute.routeEval rows)
    (fun regions => (piPMF (fun _ : Fin R => ConcreteEncoder.transvectionLaw)).prob
      (fun ts => ConcreteEncoder.outputWeight (ConcreteRoute.reshape e regions) ts ∅ ≤ d)),
    ConcreteRoute.permutation_law_eq]

theorem failureProbability_fair_eq {L b R : ℕ} (S : Finset (Fin L))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) :
    (fairExperimentLaw S b R).prob (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2 ≤ d) =
      (regionsLaw S b).expect
      (fun regions => (piPMF (fun _ : Fin R => ConcreteEncoder.transvectionLaw)).prob
        (fun ts => ConcreteEncoder.outputWeight (ConcreteRoute.reshape e regions) ts ∅ ≤ d)) := by
  rw [fairExperimentLaw, FinPMF.prob_prod_eq_expect _ _
    (fun rows ω => ConcreteRoutedEncoder.weight e rows ω ≤ d)]
  change (fairRows S b).expect (fun rows => failureProbability e rows d) = _
  simp only [failureProbability_eq_route]
  rw [← FinPMF.expect_bind, fairRows_route]

/-- Finite fair-row comparison of the actual low-output probabilities. -/
theorem activeMessages_failure_le_fair {L k R : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) {B : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ)) :
    (∑ x ∈ activeMessages S, failureProbability e (rowSupports seed x) d) ≤
      ((2 : ℝ)^(k * 24) * B)^S.card *
        (fairExperimentLaw S (k * 24) R).prob
          (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2 ≤ d) := by
  simp only [failureProbability_eq_route, failureProbability_fair_eq]
  exact activeMessages_route_expect_le seed S hB _
    (fun _ => FinPMF.prob_nonneg _ _)

/-- The complete finite sparse bound with row-counting and conditioning costs retained. -/
theorem activeMessages_sparse_failure {L k R : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) {B α : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ activeMessages S, failureProbability e (rowSupports seed x) d) ≤
      ((2 : ℝ)^(k * 24) * B)^S.card *
        ((((1 / markMass L S.card ((8/5)*α)) ^ (k * 24)) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d)) := by
  have hB0 : 0 ≤ B := by
    have hh := hB 0 (by omega)
    simp only [Nat.choose_zero_right, Nat.cast_one, mul_one] at hh
    exact (Nat.cast_nonneg _).trans hh
  exact (activeMessages_failure_le_fair seed S e d hB).trans
    (mul_le_mul_of_nonneg_left (fairExperiment_sparse_probability S e hα0 hα1 d)
      (by positivity))

/-- Active positions partition all message tuples exactly. -/
def activePositions {L k : ℕ} (x : Fin L → LocalMessage k) : Finset (Fin L) :=
  univ.filter (fun i => x i ≠ 0)

theorem activeMessages_eq_filter {L k : ℕ} (S : Finset (Fin L)) :
    activeMessages (k := k) S = univ.filter (fun x => activePositions x = S) := by
  ext x
  simp only [activeMessages, Finset.mem_filter, Finset.mem_univ, true_and]
  constructor
  · intro h
    ext i
    simpa only [activePositions, mem_filter, mem_univ, true_and] using h i
  · intro h i
    have hi := Finset.ext_iff.mp h i
    simpa only [activePositions, mem_filter, mem_univ, true_and] using hi

theorem sum_occupation_eq_activeSets {L k : ℕ} (Q : ℕ)
    (f : (Fin L → LocalMessage k) → ℝ) :
    (∑ x ∈ univ.filter (fun x => occupation x = Q), f x) =
      ∑ S ∈ (univ : Finset (Fin L)).powersetCard Q, ∑ x ∈ activeMessages S, f x := by
  simp only [activeMessages_eq_filter, Finset.sum_filter]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro x _
  simp only [Finset.sum_ite_eq, Finset.sum_ite_eq', Finset.mem_powersetCard, Finset.subset_univ, true_and]
  rfl

/-- Summing the finite sparse bound over Q active positions costs choose(L,Q). -/
theorem occupation_sparse_failure {L k R : ℕ} (seed : Seed k) (Q : ℕ)
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) {B α : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ univ.filter (fun x => occupation x = Q), failureProbability e (rowSupports seed x) d) ≤
      (L.choose Q : ℝ) * (((2 : ℝ)^(k * 24) * B)^Q *
        ((((1 / markMass L Q ((8/5)*α)) ^ (k * 24)) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d))) := by
  let C : ℝ := ((2 : ℝ)^(k * 24) * B)^Q *
    ((((1 / markMass L Q ((8/5)*α)) ^ (k * 24)) *
      (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d))
  rw [sum_occupation_eq_activeSets]
  change (∑ S ∈ (univ : Finset (Fin L)).powersetCard Q,
    ∑ x ∈ activeMessages S, failureProbability e (rowSupports seed x) d) ≤ (L.choose Q : ℝ) * C
  calc
    (∑ S ∈ (univ : Finset (Fin L)).powersetCard Q,
      ∑ x ∈ activeMessages S, failureProbability e (rowSupports seed x) d) ≤
      (∑ _S ∈ (univ : Finset (Fin L)).powersetCard Q, C) := by
      apply Finset.sum_le_sum
      intro S hS
      have hc := (Finset.mem_powersetCard.mp hS).2
      simpa only [hc] using activeMessages_sparse_failure seed S e d hB hα0 hα1
    _ = (L.choose Q : ℝ) * C := by simp [Finset.card_powersetCard]

end Spin.Structured.ConcreteOuter

