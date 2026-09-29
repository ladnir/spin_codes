import SpinCodes.Structured.ConcreteImpulseCoarseProductApprox

/-! A uniform deterministic-gap error estimate, ready for good-placement events. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder

theorem emptyError_gap_bound {c : ℝ} (hc : 0 ≤ c) {H G g : Nat}
    (hH : H ≤ g) (hG : g ≤ G) :
    emptyError c g ≤ 524288 * (c * (128 * Real.sqrt (3 * G)) + (1 / 2 : ℝ) ^ H) := by
  have hs : Real.sqrt (3 * (g : ℝ)) ≤ Real.sqrt (3 * (G : ℝ)) := by
    apply Real.sqrt_le_sqrt
    have hcast : (g : ℝ) ≤ G := by exact_mod_cast hG
    linarith
  have hp : (1 / 2 : ℝ) ^ g ≤ (1 / 2 : ℝ) ^ H :=
    pow_le_pow_of_le_one (by norm_num) (by norm_num) hH
  unfold emptyError
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  exact add_le_add (mul_le_mul_of_nonneg_left
    (mul_le_mul_of_nonneg_left hs (by norm_num)) hc) hp

theorem impulseErrorBudget_gap_bound {a : Nat} {c : ℝ} (hc : 0 ≤ c)
    (gaps : Fin a → Nat) (last H G : Nat)
    (hH : ∀ i, H ≤ gaps i) (hG : ∀ i, gaps i ≤ G) (hlH : H ≤ last) (hlG : last ≤ G) :
    impulseErrorBudget c gaps last ≤
      ((a : ℝ) + 1) * (524288 * (c * (128 * Real.sqrt (3 * G)) + (1 / 2 : ℝ) ^ H)) +
        a * (128 * c) := by
  rw [impulseErrorBudget_sum]
  let E := 524288 * (c * (128 * Real.sqrt (3 * G)) + (1 / 2 : ℝ) ^ H)
  have he (i : Fin a) : emptyError c (gaps i) ≤ E := emptyError_gap_bound hc (hH i) (hG i)
  have hl : emptyError c last ≤ E := emptyError_gap_bound hc hlH hlG
  calc emptyError c last + ∑ i, (emptyError c (gaps i) + 128 * c) ≤
      E + ∑ _i : Fin a, (E + 128 * c) :=
        add_le_add hl (Finset.sum_le_sum (fun i _ => add_le_add (he i) le_rfl))
    _ = _ := by
      simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
      dsimp only [E]
      ring

theorem endpointKernel_uniform_impulses_good_gaps {a : Nat} {θ L : ℝ}
    (hθ : 0 ≤ θ) (hL : 0 < L) (gaps : Fin a → Nat) (last H G : Nat)
    (hH : ∀ i, H ≤ gaps i) (hG : ∀ i, gaps i ≤ G) (hlH : H ≤ last) (hlG : last ≤ G)
    (q r : State) :
    |(matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
        (fun coords => endpointKernel (Real.exp (-(θ / L))) (impulseInputs gaps coords last).get)) q r -
      (liveLift * coarseImpulseProduct (θ / L) gaps last * liveProjection) q r| ≤
      ((a : ℝ) + 1) * (524288 * ((θ / L) * (128 * Real.sqrt (3 * G)) + (1 / 2 : ℝ) ^ H)) +
        a * (128 * (θ / L)) :=
  (endpointKernel_uniform_impulses_coarse_error hθ hL gaps last q r).trans
    (impulseErrorBudget_gap_bound (div_nonneg hθ hL.le) gaps last H G hH hG hlH hlG)

end Spin.Structured.ConcreteEncoder
