import SpinCodes.Structured.ConcreteFourierTransfer
import SpinCodes.Structured.PositiveFiniteTilted

noncomputable section
namespace Spin.Structured.ConcreteFourier
open Spin.Imt ConcreteRoute ConcreteRoutedEncoder ConcreteEncoder ConcreteMaps

theorem iterate_le {p z lam wmin : ℝ} {w : Coords 5}
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hz : 0 ≤ z)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((matrix p z).applyCol w).le (Coords.smul lam w)) (R : ℕ) :
    (((matrix p z).apply)^[R] (Coords.eZ 5)).total ≤ lam^R*w.Z/wmin := by
  apply (le_div_iff₀ hwmin).mpr
  exact collatz_total (matrix_nonneg hp0 hp1 hz) hw hwmin hwZ hwD hwS R

/-- The Fourier transfer bounds the actual independent inputs and transvections. -/
theorem encoder_moment {p z lam wmin : ℝ} {w : Coords 5}
    (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 < z)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((matrix p z).applyCol w).le (Coords.smul lam w)) (R : ℕ) :
    averagedMoment (poissonBinom (fun _ : Fin 128 => p) (fun _ => hp0.le) (fun _ => hp1.le))
      R z ∅ ≤ lam^R*w.Z/wmin := by
  unfold averagedMoment
  rw [averagedMoment_eq_iterate _ p z (fun x => poissonBinom_const_apply _ _ x)]
  exact (moment_bound hp0 hp1 hz R).trans (iterate_le hp0.le hp1.le hz.le hwmin hwZ hwD hwS hw R)

/-- Actual routed moment, with the exact fixed-weight likelihood factor. -/
theorem routed_moment {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z lam wmin : ℝ} {w : Coords 5} (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 < z)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((matrix p z).applyCol w).le (Coords.smul lam w)) :
    (ConcreteRoutedEncoder.experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) ≤
      ((((b : ℝ)+1)^activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (density rows) p)) * (lam^R*w.Z/wmin) := by
  rw [moment_eq]
  rw [← FinPMF.expect_map (seedLaw L b) (routeEval rows)
    (fun regions => inputMoment z (reshape e regions) ∅)]
  refine ((PositiveFinite.routed_dominates_tilted hL hb rows hw0 hw1 hp0 hp1).expect_le
    (fun regions => inputMoment_nonneg hz.le (reshape e regions) ∅)).trans ?_
  rw [iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact encoder_moment hp0 hp1 hz hwmin hwZ hwD hwS hw R

theorem routed_probability {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z lam wmin : ℝ} {w : Coords 5} (hp0 : 0 < p) (hp1 : p < 1)
    (hz : 0 < z) (hz1 : z ≤ 1)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((matrix p z).applyCol w).le (Coords.smul lam w)) (d : ℕ) :
    (ConcreteRoutedEncoder.experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      ((((b : ℝ)+1)^activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (density rows) p)) * (lam^R*w.Z/wmin) / z^d :=
  FinPMF.prob_weight_le_of_moment _ _ d hz hz1
    (routed_moment hL hb e rows hw0 hw1 hp0 hp1 hz hwmin hwZ hwD hwS hw)

end Spin.Structured.ConcreteFourier
