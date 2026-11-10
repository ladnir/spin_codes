import SpinCodes.Turbo.PhysicalGenericLive22
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# Feedback-event Cauchy bounds for the actual physical live rows

The old centered-square argument bounds an unrestricted first moment by its
second moment. Here finite weighted Cauchy additionally retains the probability
of the feedback-nonzero event. No feedback/output independence is assumed.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer
attribute [local instance] Classical.propDecidable
attribute [local irreducible] stepWeight stepFeedback packetStepInput22
local instance physicalState22DecidableEq : DecidableEq State22 :=
  fun a b => Classical.propDecidable (a = b)

/-- Weighted finite Cauchy with an event factor. -/
theorem physical_expect_event_sq_le {α : Type*} [Fintype α]
    (P : FinPMF α) (f : α → ℝ) (E : α → Prop) [DecidablePred E] :
    (P.expect (fun x => if E x then f x else 0)) ^ 2 ≤
      P.expect (fun x => f x ^ 2) * P.prob E := by
  rw [FinPMF.prob_eq_expect_indicator]
  unfold FinPMF.expect
  apply Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul
  · intro x _; exact mul_nonneg (P.nonneg x) (sq_nonneg _)
  · intro x _
    change 0 ≤ P.p x * (if E x then 1 else 0)
    split_ifs
    · simpa only [mul_one] using P.nonneg x
    · simp only [mul_zero, le_refl]
  · intro x _
    change (P.p x * (if E x then f x else 0)) ^ 2 ≤
      (P.p x * f x ^ 2) * (P.p x * (if E x then 1 else 0))
    apply le_of_eq
    by_cases hx : E x <;> simp only [hx, ite_true, ite_false] <;> ring
theorem physical_expect_event_le_sqrt {α : Type*} [Fintype α]
    (P : FinPMF α) (f : α → ℝ) (E : α → Prop) [DecidablePred E] :
    P.expect (fun x => if E x then f x else 0) ≤
      Real.sqrt (P.expect (fun x => f x ^ 2) * P.prob E) :=
  Real.le_sqrt_of_sq_le (physical_expect_event_sq_le P f E)

def physicalFeedbackProbability22 (P : FinPMF PacketStep16) : ℝ :=
  P.prob (fun x => stepFeedback (packetStepInput22 x) ≠ 0)

def physicalInputSecondMoment22 (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) : ℝ := P.expect (fun x => (z ^ stepWeight (packetStepInput22 x) q) ^ 2)

theorem physicalInputSecondMoment22_nonneg (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) : 0 ≤ physicalInputSecondMoment22 P z q :=
  P.expect_nonneg (fun _ => sq_nonneg _)

theorem physicalInputSecondMoment22_eq (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) : physicalInputSecondMoment22 P z q = physicalInputMoment22 P (z*z) q := by
  unfold physicalInputSecondMoment22 physicalInputMoment22
  apply congrArg
  funext x
  simp only [mul_pow, pow_two]

theorem physicalInputReturn22_le_cauchy (P : FinPMF PacketStep16) (z : ℝ)
    (q : State22) :
    physicalInputReturn22 P z q ≤ Real.sqrt
      (physicalInputSecondMoment22 P z q * physicalFeedbackProbability22 P) := by
  have he : physicalInputReturn22 P z q =
      P.expect (fun x => if stepFeedback (packetStepInput22 x) ≠ 0 then
        z ^ stepWeight (packetStepInput22 x) q else 0) := by
    unfold physicalInputReturn22
    apply congrArg
    funext x
    split_ifs <;> first | rfl | contradiction
  rw [he]
  exact physical_expect_event_le_sqrt P
    (fun x => z ^ stepWeight (packetStepInput22 x) q)
    (fun x => stepFeedback (packetStepInput22 x) ≠ 0)
theorem physicalInputReturn22_le_first (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (q : State22) : physicalInputReturn22 P z q ≤ physicalInputMoment22 P z q := by
  apply FinPMF.expect_mono
  intro x
  split_ifs
  · exact pow_nonneg hz _
  · exact le_rfl

theorem physicalInputReturn22_nonneg (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (q : State22) : 0 ≤ physicalInputReturn22 P z q := by
  apply FinPMF.expect_nonneg
  intro x
  split_ifs
  · exact le_rfl
  · exact pow_nonneg hz _

/-- The retained numeric row uses first/second-moment and feedback caps on
one physical input law; these inequalities are still explicit obligations. -/
theorem physicalInputReturn22_le_numeric_cauchy (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (q : State22) (M S F : ℝ)
    (hS : 0 ≤ S) (hfirst : physicalInputMoment22 P z q ≤ M)
    (hsecond : physicalInputSecondMoment22 P z q ≤ S)
    (hfeedback : physicalFeedbackProbability22 P ≤ F) :
    physicalInputReturn22 P z q ≤ min M (Real.sqrt (S*F)) := by
  apply le_min
  · exact (physicalInputReturn22_le_first P z hz q).trans hfirst
  · apply (physicalInputReturn22_le_cauchy P z q).trans
    apply Real.sqrt_le_sqrt
    exact mul_le_mul hsecond hfeedback
      (FinPMF.prob_nonneg P _) hS

/-- Complete exact-nonzero live row from classwise moment/census caps. -/
theorem turboPhysicalKernel22_live_row_of_max_cauchy (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ)
    (M S F : ℝ) (hM : 0 ≤ M) (hS : 0 ≤ S)
    (hfirst : ∀ q, μ q ≠ 0 → physicalInputMoment22 (packetStepLaw22 j) z q ≤ M)
    (hsecond : ∀ q, μ q ≠ 0 → physicalInputSecondMoment22 (packetStepLaw22 j) z q ≤ S)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤ F) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k
        (b * physicalLiveCapRow22 M (min M (Real.sqrt (S*F))) k) (parts k) := by
  apply turboPhysicalKernel22_live_row_of_caps j z hz i hi b μ hμ
  · exact physicalLivePiece22_max_bound hi hμ _ M hM hfirst
  · apply physicalLivePiece22_max_bound hi hμ _ _ (le_min hM (Real.sqrt_nonneg _))
    intro q hq
    exact physicalInputReturn22_le_numeric_cauchy _ z hz q M S F hS
      (hfirst q hq) (hsecond q hq) hfeedback

private theorem weighted_event_sq_le {α : Type*} [Fintype α]
    (w f : α → ℝ) (E : α → Prop) [DecidablePred E] (hw : ∀ x, 0 ≤ w x) :
    (∑ x, w x * (if E x then f x else 0)) ^ 2 ≤
      (∑ x, w x * f x ^ 2) * ∑ x, w x * (if E x then 1 else 0) := by
  apply Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul
  · intro x _; exact mul_nonneg (hw x) (sq_nonneg _)
  · intro x _
    by_cases hx : E x <;> simp only [hx, ite_true, ite_false, mul_one, mul_zero]
    · exact hw x
    · exact le_rfl
  · intro x _
    apply le_of_eq
    by_cases hx : E x <;> simp only [hx, ite_true, ite_false] <;> ring

/-- Cauchy is applied jointly to the entering measure and complete input word.
In particular this covers a uniformly capped source without replacing its
mean return moment by a pointwise maximum. -/
theorem physicalInputReturn22_weighted_sq_le (P : FinPMF PacketStep16) (z : ℝ)
    (μ : State22 → ℝ) (hμ : ∀ q, 0 ≤ μ q) :
    (∑ q, μ q * physicalInputReturn22 P z q) ^ 2 ≤
      (∑ q, μ q * physicalInputSecondMoment22 P z q) *
        ((∑ q, μ q) * physicalFeedbackProbability22 P) := by
  have h := weighted_event_sq_le
    (fun qx : State22 × PacketStep16 => μ qx.1 * P.p qx.2)
    (fun qx => z ^ stepWeight (packetStepInput22 qx.2) qx.1)
    (fun qx => stepFeedback (packetStepInput22 qx.2) ≠ 0)
    (fun qx => mul_nonneg (hμ qx.1) (P.nonneg qx.2))
  have he (q : State22) : physicalInputReturn22 P z q =
      P.expect (fun x => if stepFeedback (packetStepInput22 x) ≠ 0 then
        z ^ stepWeight (packetStepInput22 x) q else 0) := by
    unfold physicalInputReturn22
    apply congrArg
    funext x
    split_ifs <;> first | rfl | contradiction
  have ha : (∑ qx : State22 × PacketStep16,
      (μ qx.1 * P.p qx.2) * (if stepFeedback (packetStepInput22 qx.2) ≠ 0 then
        z ^ stepWeight (packetStepInput22 qx.2) qx.1 else 0)) =
      ∑ q, μ q * physicalInputReturn22 P z q := by
    rw [Fintype.sum_prod_type]
    simp only [mul_assoc, ← Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro q _
    change μ q * P.expect (fun x => if stepFeedback (packetStepInput22 x) ≠ 0 then
      z ^ stepWeight (packetStepInput22 x) q else 0) = μ q * physicalInputReturn22 P z q
    rw [he q]
  have hb : (∑ qx : State22 × PacketStep16,
      (μ qx.1 * P.p qx.2) * (z ^ stepWeight (packetStepInput22 qx.2) qx.1)^2) =
      ∑ q, μ q * physicalInputSecondMoment22 P z q := by
    rw [Fintype.sum_prod_type]
    simp only [mul_assoc, ← Finset.mul_sum]
    rfl
  have hc : (∑ qx : State22 × PacketStep16,
      (μ qx.1 * P.p qx.2) * (if stepFeedback (packetStepInput22 qx.2) ≠ 0 then (1 : ℝ) else 0)) =
      (∑ q, μ q) * physicalFeedbackProbability22 P := by
    rw [Fintype.sum_prod_type]
    simp only [mul_assoc, ← Finset.mul_sum]
    unfold physicalFeedbackProbability22
    rw [FinPMF.prob_eq_expect_indicator, Finset.sum_mul]
    rfl
  rw [ha, hb, hc] at h
  exact h
/-- The retained Cauchy cap from weighted first/second moments of one source. -/
theorem physicalInputReturn22_weighted_cauchy (P : FinPMF PacketStep16) (z : ℝ)
    (hz : 0 ≤ z) (μ : State22 → ℝ) (hμ : ∀ q, 0 ≤ μ q)
    (b S F : ℝ) (hb : 0 ≤ b) (hS : 0 ≤ S)
    (hmass : (∑ q, μ q) ≤ b)
    (hsecond : (∑ q, μ q * physicalInputSecondMoment22 P z q) ≤ b*S)
    (hfeedback : physicalFeedbackProbability22 P ≤ F) :
    (∑ q, μ q * physicalInputReturn22 P z q) ≤ b * Real.sqrt (S*F) := by
  have hp : 0 ≤ physicalFeedbackProbability22 P := FinPMF.prob_nonneg P _
  have hF : 0 ≤ F := hp.trans hfeedback
  have hevent : (∑ q, μ q) * physicalFeedbackProbability22 P ≤ b*F :=
    mul_le_mul hmass hfeedback hp hb
  have hevent0 : 0 ≤ (∑ q, μ q) * physicalFeedbackProbability22 P :=
    mul_nonneg (Finset.sum_nonneg (fun q _ => hμ q)) hp
  have hsquare : (∑ q, μ q * physicalInputReturn22 P z q)^2 ≤ (b*S)*(b*F) :=
    (physicalInputReturn22_weighted_sq_le P z μ hμ).trans
      (mul_le_mul hsecond hevent hevent0 (mul_nonneg hb hS))
  calc
    _ ≤ Real.sqrt ((b*S)*(b*F)) := Real.le_sqrt_of_sq_le hsquare
    _ = b * Real.sqrt (S*F) := by
      rw [show (b*S)*(b*F) = b^2*(S*F) by ring,
        Real.sqrt_mul (sq_nonneg b), Real.sqrt_sq_eq_abs, abs_of_nonneg hb]

/-- A complete actual exact-nonzero row from weighted moment caps. This is the
joint Cauchy interface used by both class-max and uniformly capped mean rows. -/
theorem turboPhysicalKernel22_live_row_of_weighted_cauchy (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ)
    (M S F : ℝ) (hS : 0 ≤ S)
    (hfirst : (∑ q, μ q * physicalInputMoment22 (packetStepLaw22 j) z q) ≤ b*M)
    (hsecond : (∑ q, μ q * physicalInputSecondMoment22 (packetStepLaw22 j) z q) ≤ b*S)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤ F) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k
        (b * physicalLiveCapRow22 M (min M (Real.sqrt (S*F))) k) (parts k) := by
  have hb : 0 ≤ b := (Finset.sum_nonneg (fun q _ => hμ.1 q)).trans hμ.2.2.1
  apply turboPhysicalKernel22_live_row_of_caps j z hz i hi b μ hμ
  · exact hfirst
  · by_cases hMS : M ≤ Real.sqrt (S*F)
    · rw [min_eq_left hMS]
      apply le_trans (Finset.sum_le_sum (fun q _ =>
        mul_le_mul_of_nonneg_left (physicalInputReturn22_le_first _ z hz q) (hμ.1 q))) hfirst
    · rw [min_eq_right (le_of_not_ge hMS)]
      exact physicalInputReturn22_weighted_cauchy _ z hz μ hμ.1 b S F hb hS
        hμ.2.2.1 hsecond hfeedback

/-- The uniform-live source uses means of both moments, with Cauchy applied
jointly before the complete row is emitted. -/
theorem turboPhysicalKernel22_uniform_row_of_mean_cauchy (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 2 b μ) (M S F : ℝ) (hS : 0 ≤ S)
    (hfirst : turboMeanLive22 (physicalInputMoment22 (packetStepLaw22 j) z) ≤ M)
    (hsecond : turboMeanLive22 (physicalInputSecondMoment22 (packetStepLaw22 j) z) ≤ S)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤ F) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22 (Matrix.of (turboPhysicalKernel22 j z)) μ r = ∑ k, parts k r) ∧
      ∀ k, physicalPiece22 k
        (b * physicalLiveCapRow22 M (min M (Real.sqrt (S*F))) k) (parts k) := by
  apply turboPhysicalKernel22_live_row_of_weighted_cauchy j z hz 2 (by decide) b μ hμ M S F hS
  · apply physicalUniformPiece22_mean_bound hμ _ _ M hfirst
    intro q
    exact (packetStepLaw22 j).expect_nonneg (fun x => pow_nonneg hz _)
  · exact physicalUniformPiece22_mean_bound hμ _
      (physicalInputSecondMoment22_nonneg _ z) S hsecond
  · exact hfeedback

#print axioms physicalInputReturn22_weighted_sq_le
#print axioms turboPhysicalKernel22_live_row_of_weighted_cauchy
#print axioms turboPhysicalKernel22_uniform_row_of_mean_cauchy#print axioms physical_expect_event_sq_le
#print axioms physicalInputReturn22_le_cauchy
#print axioms turboPhysicalKernel22_live_row_of_max_cauchy
end Spin.Turbo








