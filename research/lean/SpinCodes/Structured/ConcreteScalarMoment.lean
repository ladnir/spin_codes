import SpinCodes.Structured.ConcreteScalarBits
import SpinCodes.Structured.ConcreteEncoderMoment

noncomputable section
namespace Spin.Structured.ConcreteScalar
open Finset ConcreteEncoder ConcreteMaps
open scoped symmDiff

def scalarF (β z : ℝ) (d : ℕ) : ℝ := (1-β+β*z)^(128-d) * (β+(1-β)*z)^d

def scalarBound (β z : ℝ) : ℝ := max (scalarF β z 0)
  (univ.sup' univ_nonempty (fun i : Fin 5 => scalarF β z (shellWeight i)))

theorem scalarF_le_bound (β z : ℝ) (q : State) : scalarF β z (Aset q).card ≤ scalarBound β z := by
  by_cases hq : q = ∅
  · subst q
    simp only [Aset_empty, card_empty]
    exact le_max_left _ _
  · have hc := weightShell_covers actual_shell_counts
    have hmem : q ∈ (univ : Finset (Fin 5)).biUnion weightShell := by
      rw [hc]
      exact Spin.mem_nonzeroStates hq
    obtain ⟨i, _, hi⟩ := mem_biUnion.mp hmem
    rw [weightShell_weight hi]
    exact (le_sup' (fun i : Fin 5 => scalarF β z (shellWeight i)) (mem_univ i)).trans (le_max_right _ _)

theorem emission_expect (β z : ℝ) (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q : State) :
    (poissonBinom (fun _ : Fin 128 => β) (fun _ => hβ0) (fun _ => hβ1)).expect
      (emitted z q) = scalarF β z (Aset q).card :=
  bernoulli_xor_moment β z hβ0 hβ1 (Aset q)

theorem scalarBound_nonneg {β z : ℝ} (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) :
    0 ≤ scalarBound β z := by
  apply le_trans _ (le_max_left _ _)
  unfold scalarF
  have : 0 ≤ 1-β := by linarith
  positivity

theorem averagedMoment_succ_expect (P : FinPMF Input) (R : ℕ) (z : ℝ) (q : State) :
    averagedMoment P (R+1) z q = P.expect (fun x => emitted z q x *
      transvectionLaw.expect (fun t => averagedMoment P R z (nextState q x t))) := by
  rw [averagedMoment_succ]
  simp_rw [nextState_expect]
  simp only [FinPMF.expect, mul_assoc]

/-- A scalar one-round envelope propagates through the actual dependent-state recurrence. -/
theorem averagedMoment_le_scalar (P : FinPMF Input) {z lam : ℝ} (hz : 0 ≤ z)
    (hlam : 0 ≤ lam) (hrow : ∀ q, P.expect (emitted z q) ≤ lam) (R : ℕ) (q : State) :
    averagedMoment P R z q ≤ lam^R := by
  induction R generalizing q with
  | zero => rw [averagedMoment_zero, pow_zero]
  | succ R ih =>
    rw [averagedMoment_succ_expect]
    calc
      _ ≤ P.expect (fun x => emitted z q x * lam^R) := by
        apply FinPMF.expect_mono
        intro x
        apply mul_le_mul_of_nonneg_left _ (pow_nonneg hz _)
        exact (FinPMF.expect_mono _ (fun t => ih (nextState q x t))).trans_eq (FinPMF.expect_const _ _)
      _ = P.expect (emitted z q) * lam^R := by
        simp only [FinPMF.expect, Finset.sum_mul, mul_assoc]
      _ ≤ lam * lam^R := mul_le_mul_of_nonneg_right (hrow q) (pow_nonneg hlam _)
      _ = lam^(R+1) := by ring

theorem bernoulli_encoder_scalar {β z : ℝ} (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z)
    (R : ℕ) (q : State) :
    averagedMoment (poissonBinom (fun _ : Fin 128 => β) (fun _ => hβ0) (fun _ => hβ1)) R z q ≤
      (scalarBound β z)^R := by
  apply averagedMoment_le_scalar _ hz (scalarBound_nonneg hβ0 hβ1 hz)
  intro q
  rw [emission_expect]
  exact scalarF_le_bound β z q

/-- At fair input every entering state has the same exact emission moment. -/
theorem scalarF_half (z : ℝ) (d : ℕ) (hd : d ≤ 128) :
    scalarF (1/2) z d = ((1+z)/2)^128 := by
  unfold scalarF
  have h0 : (1 : ℝ)-1/2+1/2*z = (1+z)/2 := by ring
  have h1 : (1/2 : ℝ)+(1-1/2)*z = (1+z)/2 := by ring
  rw [h0, h1, ← pow_add, Nat.sub_add_cancel hd]

end Spin.Structured.ConcreteScalar

