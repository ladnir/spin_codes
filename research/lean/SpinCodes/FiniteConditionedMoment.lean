import SpinCodes.FiniteMoment

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin
open Finset
namespace FinPMF
variable {Ω : Type*} [Fintype Ω]

theorem expect_condition_mono (P : FinPMF Ω) (G : Ω → Prop) (hG : 0 < P.prob G)
    {f g : Ω → ℝ} (h : ∀ ω, G ω → f ω ≤ g ω) :
    (P.condition G hG).expect f ≤ (P.condition G hG).expect g := by
  unfold expect
  apply sum_le_sum
  intro ω _
  by_cases hg : G ω
  · exact mul_le_mul_of_nonneg_left (h ω hg) ((P.condition G hG).nonneg ω)
  · simp [condition_p_apply, hg]

end FinPMF
namespace Setup
variable {M W Ωout Ωin : Type*}
variable [Fintype M] [DecidableEq M] [Zero M] [Fintype Ωout] [Fintype Ωin]

theorem expect_ZQ_eq_outerLaw (S : Setup M W Ωout Ωin) (P : FinPMF Ωout) (d : ℕ)
    (occ : Ωout → M → ℕ) (Q : ℕ) :
    (P.prod S.Pin).expect (fun ω => (S.ZQ d occ Q ω : ℝ)) =
      P.expect (fun ω => ∑ x ∈ (nonzeroMsgs M).filter (fun x => occ ω x = Q),
        S.qd d (S.outer ω x)) := by
  exact expect_ZQ_eq { S with Pout := P } d occ Q

/-- Conditioning restricts which outer realizations need the pointwise inner envelope. -/
theorem expect_ZQ_condition_le (S : Setup M W Ωout Ωin) (d : ℕ)
    (occ : Ωout → M → ℕ) (Q : ℕ) (G : Ωout → Prop) (hG : 0 < S.Pout.prob G)
    (B : Ωout → M → ℝ)
    (hB : ∀ ω, G ω → ∀ x, x ∈ nonzeroMsgs M → occ ω x = Q →
      S.qd d (S.outer ω x) ≤ B ω x) :
    ((S.Pout.condition G hG).prod S.Pin).expect (fun ω => (S.ZQ d occ Q ω : ℝ)) ≤
      (S.Pout.condition G hG).expect
        (fun ω => ∑ x ∈ (nonzeroMsgs M).filter (fun x => occ ω x = Q), B ω x) := by
  rw [expect_ZQ_eq_outerLaw]
  apply FinPMF.expect_condition_mono
  intro ω hg
  apply sum_le_sum
  intro x hx
  exact hB ω hg x (mem_filter.mp hx).1 (mem_filter.mp hx).2

end Setup
end Spin
