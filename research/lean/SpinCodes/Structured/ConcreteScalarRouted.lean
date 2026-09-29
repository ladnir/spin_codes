import SpinCodes.Structured.ConcreteScalarMoment
import SpinCodes.Structured.PositiveFiniteTilted

noncomputable section
namespace Spin.Structured.ConcreteScalar
open Finset ConcreteRoute ConcreteRoutedEncoder ConcreteEncoder

/-- The scalar alternative, with the same actual route and exact likelihood cost. -/
theorem routed_scalar_moment {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 ≤ z) :
    (ConcreteRoutedEncoder.experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) ≤
      ((((b : ℝ)+1)^activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (density rows) p)) * (scalarBound p z)^R := by
  rw [moment_eq]
  rw [← FinPMF.expect_map (seedLaw L b) (routeEval rows)
    (fun regions => inputMoment z (reshape e regions) ∅)]
  refine ((PositiveFinite.routed_dominates_tilted hL hb rows hw0 hw1 hp0 hp1).expect_le
    (fun regions => inputMoment_nonneg hz (reshape e regions) ∅)).trans ?_
  rw [iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact bernoulli_encoder_scalar hp0.le hp1.le hz R ∅

theorem routed_scalar_probability {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 < z) (hz1 : z ≤ 1) (d : ℕ) :
    (ConcreteRoutedEncoder.experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      ((((b : ℝ)+1)^activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (density rows) p)) * (scalarBound p z)^R / z^d :=
  FinPMF.prob_weight_le_of_moment _ _ d hz hz1 (routed_scalar_moment hL hb e rows hw0 hw1 hp0 hp1 hz.le)

end Spin.Structured.ConcreteScalar
