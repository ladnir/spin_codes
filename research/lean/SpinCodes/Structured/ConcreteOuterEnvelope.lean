import SpinCodes.Structured.ConcreteOuter
import SpinCodes.Structured.ConcreteOuterEnvelopeCoefficient
import SpinCodes.Structured.ConcreteOuterEntropyTransition

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

/-- Exactly the nonzero finite summands in the Golay-BAA spectrum formula. -/
def supportedPaths (k w : ℕ) : Finset (ℕ × ℕ) :=
  ((Icc 1 (k * 24)) ×ˢ (range (k * 24 + 1))).filter (fun p =>
    golayCoefficient k p.1 ≠ 0 ∧ accT (k * 24) p.1 p.2 ≠ 0 ∧ accT (k * 24) p.2 w ≠ 0)

def pathEntropy (k w : ℕ) (p : ℕ × ℕ) : ℝ :=
  ((k * 24 : ℕ) : ℝ) * gBA ((p.1 : ℝ) / (k * 24 : ℕ)) +
    transitionEntropy (k * 24) p.1 p.2 + transitionEntropy (k * 24) p.2 w

/-- The maximum uses only actual nonzero finite paths; zero-support cases get value zero. -/
def finiteEnvelope (k w : ℕ) : ℝ :=
  if h : (supportedPaths k w).Nonempty then (supportedPaths k w).sup' h (pathEntropy k w) else 0

theorem pathEntropy_le_envelope {k w : ℕ} {p : ℕ × ℕ} (hp : p ∈ supportedPaths k w) :
    pathEntropy k w p ≤ finiteEnvelope k w := by
  have hn : (supportedPaths k w).Nonempty := ⟨p, hp⟩
  rw [finiteEnvelope, dif_pos hn]
  exact Finset.le_sup' _ hp

theorem expected_spectrum_coefficient (k w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) =
      ∑ a ∈ Icc 1 (k * 24), (golayCoefficient k a : ℝ) *
        ∑ c ∈ range (k * 24 + 1), Pt (k * 24) a c * Pt (k * 24) c w := by
  simpa only [golayCoefficient, golayPolynomial, Golay.enumerator_eq] using expected_spectrum k w

theorem spectrum_term_le_envelope {k : ℕ} (hk : 0 < k) (w a c : ℕ)
    (ha : a ∈ Icc 1 (k * 24)) (hc : c ∈ range (k * 24 + 1)) :
    (golayCoefficient k a : ℝ) * Pt (k * 24) a c * Pt (k * 24) c w ≤
      (((k * 24 : ℕ) : ℝ) + 1)^2 * Real.exp (finiteEnvelope k w) := by
  have haB := (Finset.mem_Icc.mp ha).2
  have hcB : c ≤ k * 24 := by simpa only [mem_range, Nat.lt_succ_iff] using hc
  have hpnon (v t : ℕ) : 0 ≤ Pt (k * 24) v t := by unfold Pt; positivity
  by_cases hp : (a,c) ∈ supportedPaths k w
  · have h₁ := golayCoefficient_le_exp_gBA hk a
    have h₂ := transition_le_exp_entropy haB c
    have h₃ := transition_le_exp_entropy hcB w
    calc (golayCoefficient k a : ℝ) * Pt (k * 24) a c * Pt (k * 24) c w ≤
        Real.exp (((k * 24 : ℕ) : ℝ) * gBA ((a : ℝ) / (k * 24 : ℕ))) *
          ((((k * 24 : ℕ) : ℝ) + 1) * Real.exp (transitionEntropy (k * 24) a c)) *
          ((((k * 24 : ℕ) : ℝ) + 1) * Real.exp (transitionEntropy (k * 24) c w)) := by
          exact mul_le_mul (mul_le_mul h₁ h₂ (hpnon _ _) (by positivity)) h₃
            (hpnon _ _) (by positivity)
      _ = (((k * 24 : ℕ) : ℝ) + 1)^2 * Real.exp (pathEntropy k w (a,c)) := by
          rw [pathEntropy, Real.exp_add, Real.exp_add]
          ring
      _ ≤ _ := mul_le_mul_of_nonneg_left
        (Real.exp_le_exp.mpr (pathEntropy_le_envelope hp)) (sq_nonneg _)
  · have hz : golayCoefficient k a = 0 ∨ accT (k * 24) a c = 0 ∨ accT (k * 24) c w = 0 := by
      simpa only [supportedPaths, mem_filter, mem_product, ha, hc, true_and, not_and_or, not_not] using hp
    rcases hz with hz | hz | hz <;> simp only [Pt, hz, Nat.cast_zero, zero_mul, mul_zero, zero_div] <;> positivity

/-- A fully finite spectral envelope with the explicit polynomial factor b(b+1)^3. -/
theorem expected_spectrum_finite_bound {k : ℕ} (hk : 0 < k) (w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) ≤
      ((k * 24 : ℕ) : ℝ) * (((k * 24 : ℕ) : ℝ) + 1)^3 * Real.exp (finiteEnvelope k w) := by
  rw [expected_spectrum_coefficient]
  simp only [Finset.mul_sum]
  have hh : (∑ a ∈ Icc 1 (k * 24), ∑ c ∈ range (k * 24 + 1),
      (golayCoefficient k a : ℝ) * (Pt (k * 24) a c * Pt (k * 24) c w)) ≤
      ∑ _a ∈ Icc 1 (k * 24), ∑ _c ∈ range (k * 24 + 1),
        (((k * 24 : ℕ) : ℝ) + 1)^2 * Real.exp (finiteEnvelope k w) := by
    apply Finset.sum_le_sum
    intro a ha
    apply Finset.sum_le_sum
    intro c hc
    simpa only [mul_assoc] using spectrum_term_le_envelope hk w a c ha hc
  refine hh.trans_eq ?_
  simp only [Finset.sum_const, Finset.card_range, Nat.card_Icc, Nat.add_sub_cancel,
    nsmul_eq_mul, Nat.cast_add, Nat.cast_one]
  ring

/-- Normalizing the finite envelope hides no constant: the remaining loss is at most (b+1)^4. -/
theorem expected_spectrum_finite_bound_fourth {k : ℕ} (hk : 0 < k) (w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) ≤
      ((((k * 24 : ℕ) : ℝ) + 1)^4) * Real.exp (finiteEnvelope k w) := by
  refine (expected_spectrum_finite_bound hk w).trans ?_
  have hb : (0 : ℝ) ≤ (k * 24 : ℕ) := by positivity
  have hh := mul_le_mul_of_nonneg_right
    (show ((k * 24 : ℕ) : ℝ) ≤ ((k * 24 : ℕ) : ℝ) + 1 by linarith)
    (show (0 : ℝ) ≤ (((k * 24 : ℕ) : ℝ) + 1)^3 * Real.exp (finiteEnvelope k w) by positivity)
  convert hh using 1 <;> ring

/-- The finite remainder is explicitly at most 4 log(b+1) in the exponent. -/
theorem expected_spectrum_explicit_log_bound {k : ℕ} (hk : 0 < k) (w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) ≤
      Real.exp (finiteEnvelope k w + 4 * Real.log (((k * 24 : ℕ) : ℝ) + 1)) := by
  refine (expected_spectrum_finite_bound_fourth hk w).trans_eq ?_
  have he := Real.exp_nat_mul (Real.log (((k * 24 : ℕ) : ℝ) + 1)) 4
  norm_num only [Nat.cast_ofNat] at he
  rw [Real.exp_add, he, Real.exp_log (by positivity)]
  ring

end Spin.Structured.ConcreteOuter



