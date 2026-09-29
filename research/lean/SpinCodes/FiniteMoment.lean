import SpinCodes.Framework

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin
open Finset

namespace FinPMF
variable {Ω M : Type*} [Fintype Ω] [Fintype M]

/-- A lower output-weight tail is bounded by its generating function. -/
theorem prob_weight_le (P : FinPMF Ω) (W : Ω → ℕ) (d : ℕ)
    {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    P.prob (fun ω => W ω ≤ d) ≤ P.expect (fun ω => z ^ W ω) / z ^ d := by
  refine (P.prob_mono (fun ω h => pow_le_pow_of_le_one hz0.le hz1 h)).trans ?_
  exact P.markov (fun ω => pow_nonneg hz0.le (W ω)) (pow_pos hz0 d)

theorem prob_weight_le_of_moment (P : FinPMF Ω) (W : Ω → ℕ) (d : ℕ)
    {z B : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hB : P.expect (fun ω => z ^ W ω) ≤ B) :
    P.prob (fun ω => W ω ≤ d) ≤ B / z ^ d :=
  (P.prob_weight_le W d hz0 hz1).trans (div_le_div_of_nonneg_right hB (pow_nonneg hz0.le d))

theorem expect_card_filter (P : FinPMF Ω) (s : Finset M) (E : M → Ω → Prop) :
    P.expect (fun ω => ((s.filter (fun x => E x ω)).card : ℝ)) =
      ∑ x ∈ s, P.prob (E x) := by
  have he : (fun ω => ((s.filter (fun x => E x ω)).card : ℝ)) =
      fun ω => ∑ x ∈ s, if E x ω then (1 : ℝ) else 0 := by
    funext ω
    simp only [Finset.card_filter, Nat.cast_sum]
    apply sum_congr rfl
    intro x _
    split_ifs <;> simp
  rw [he, expect_sum]
  apply sum_congr rfl
  intro x _
  exact (prob_eq_expect_indicator _ _).symm

end FinPMF

namespace Setup
variable {M W Ωout Ωin : Type*}
variable [Fintype M] [DecidableEq M] [Zero M] [Fintype Ωout] [Fintype Ωin]

/-- The occupation count first moment retains the shared outer realization. -/
theorem expect_ZQ_eq (S : Setup M W Ωout Ωin) (d : ℕ)
    (occ : Ωout → M → ℕ) (Q : ℕ) :
    S.joint.expect (fun ω => (S.ZQ d occ Q ω : ℝ)) =
      S.Pout.expect (fun ω => ∑ x ∈ (nonzeroMsgs M).filter (fun x => occ ω x = Q),
        S.qd d (S.outer ω x)) := by
  rw [joint, FinPMF.expect_prod]
  apply congrArg S.Pout.expect
  funext ω
  have hc (ωin : Ωin) : S.ZQ d occ Q (ω, ωin) =
      (((nonzeroMsgs M).filter (fun x => occ ω x = Q)).filter
        (fun x => S.innerWt ωin (S.outer ω x) ≤ d)).card := by
    simp only [ZQ, Finset.filter_filter]
  simp only [hc]
  convert FinPMF.expect_card_filter S.Pin
    ((nonzeroMsgs M).filter (fun x => occ ω x = Q))
    (fun x ωin => S.innerWt ωin (S.outer ω x) ≤ d) using 1
  · apply congrArg S.Pin.expect
    funext ωin
    apply congrArg (fun n : ℕ => (n : ℝ))
    apply congrArg Finset.card
    ext x
    simp
  · apply sum_congr rfl
    intro x _
    unfold qd FinPMF.prob
    apply sum_congr
    · ext ωin
      simp
    · intro ωin _
      rfl

theorem expect_ZQ_le (S : Setup M W Ωout Ωin) (d : ℕ)
    (occ : Ωout → M → ℕ) (Q : ℕ) (B : Ωout → M → ℝ)
    (hB : ∀ ω x, x ∈ nonzeroMsgs M → occ ω x = Q → S.qd d (S.outer ω x) ≤ B ω x) :
    S.joint.expect (fun ω => (S.ZQ d occ Q ω : ℝ)) ≤
      S.Pout.expect (fun ω => ∑ x ∈ (nonzeroMsgs M).filter (fun x => occ ω x = Q), B ω x) := by
  rw [expect_ZQ_eq]
  apply FinPMF.expect_mono
  intro ω
  apply sum_le_sum
  intro x hx
  exact hB ω x (mem_filter.mp hx).1 (mem_filter.mp hx).2

end Setup
end Spin
