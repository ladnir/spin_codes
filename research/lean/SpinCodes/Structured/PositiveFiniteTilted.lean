import SpinCodes.Structured.PositiveFiniteDensity
import SpinCodes.Structured.PositiveFiniteCollatz

noncomputable section
namespace Spin.Structured.PositiveFinite
open Finset Spin.Imt ConcreteRoute ConcreteRoutedEncoder ConcreteEncoder

/-- Changing the iid comparison density costs the exact fixed-weight KL factor. -/
theorem routed_dominates_tilted {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < L*b) {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    Spin.Dominates ((seedLaw L b).map (routeEval rows)) (iidLaw L b p hp0.le hp1.le)
      ((((b:ℝ)+1)^activeRows rows * ((L:ℝ)+1)^b) *
        Real.exp (((L:ℝ)*b)*binKL (density rows) p)) := by
  have hq0 : 0 < density rows := div_pos (by exact_mod_cast hw0) (by positivity)
  have hq1 : density rows < 1 := (div_lt_one (by positivity)).mpr (by exact_mod_cast hw1)
  intro x
  by_cases hx : totalWeight x = totalWeight rows
  · have hr := iidLaw_weight_ratio hw0 hw1 hp0 hp1 x hx
    change (iidLaw L b (density rows) _ _).p x = _ at hr
    have hd := permutation_dominates hL hb rows hq0 hq1 x
    rw [hr] at hd
    simpa only [density, mul_assoc] using hd
  · have hz : ((seedLaw L b).map (routeEval rows)).p x = 0 := by
      rw [FinPMF.map_p]
      apply Finset.sum_eq_zero
      intro seed _
      have hne : routeEval rows seed ≠ x := by
        intro he
        exact hx (he ▸ routeEval_totalWeight rows seed)
      simp [hne]
    rw [hz]
    exact mul_nonneg (mul_nonneg (by positivity) (Real.exp_pos _).le)
      ((iidLaw L b p hp0.le hp1.le).nonneg x)

/-- Finite moment at an arbitrary reference density, with exact KL and route costs. -/
theorem routed_tilted_moment_le {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z lam wmin : ℝ} {w : Coords 5} (hp0 : 0 < p) (hp1 : p < 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((Occupation.Sparse.numericalMatrix p z).applyCol w).le (Coords.smul lam w)) :
    (experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) ≤
      ((((b:ℝ)+1)^activeRows rows * ((L:ℝ)+1)^b) *
        Real.exp (((L:ℝ)*b)*binKL (density rows) p)) * (lam^R*w.Z/wmin) := by
  rw [moment_eq]
  rw [← FinPMF.expect_map (seedLaw L b) (routeEval rows)
    (fun regions => inputMoment z (reshape e regions) ∅)]
  refine ((routed_dominates_tilted hL hb rows hw0 hw1 hp0 hp1).expect_le
    (fun regions => inputMoment_nonneg hz0 (reshape e regions) ∅)).trans ?_
  rw [iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  refine (encoder_moment_bound
    (poissonBinom (fun _ : Fin 128 => p) (fun _ => hp0.le) (fun _ => hp1.le))
    (fun x => poissonBinom_const_apply _ _ x) hp0.le hp1.le hz0 hz1 R).trans ?_
  exact occupation_iterate_le hp0.le hp1.le hz0 hwmin hwZ hwD hwS hw R

/-- Actual lower-weight probability with all finite prefactors and a tilted Collatz witness. -/
theorem routed_tilted_probability_le {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z lam wmin : ℝ} {w : Coords 5} (hp0 : 0 < p) (hp1 : p < 1)
    (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((Occupation.Sparse.numericalMatrix p z).applyCol w).le (Coords.smul lam w)) (d : ℕ) :
    (experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      ((((b:ℝ)+1)^activeRows rows * ((L:ℝ)+1)^b) *
        Real.exp (((L:ℝ)*b)*binKL (density rows) p)) * (lam^R*w.Z/wmin) / z^d :=
  FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
    (routed_tilted_moment_le hL hb e rows hw0 hw1 hp0 hp1 hz0.le hz1 hwmin hwZ hwD hwS hw)

end Spin.Structured.PositiveFinite
