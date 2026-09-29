/-
The common first-moment framework (paper: `sec:framework`).

This is `thm:class-first-moment` together with the two remarks the structured
proof depends on:

* the bad event (non-injectivity *or* small minimum distance) is exactly
  `Z_d ≥ 1`;
* the theorem applies verbatim under a law conditioned on any event that is
  determined by the outer setup alone.

Nothing here is specific to the structured construction.  The structured
instantiation chooses the class map `c` to be the ordered row-weight profile.
-/
import SpinCodes.Prob

set_option linter.unusedSectionVars false

namespace Spin

open Finset

/-- The nonzero messages of a finite message space. -/
def nonzeroMsgs (M : Type*) [Fintype M] [DecidableEq M] [Zero M] : Finset M :=
  univ.filter (fun x : M => x ≠ 0)

variable {M W Ωout Ωin : Type*}
variable [Fintype M] [DecidableEq M] [Zero M]
variable [Fintype Ωout] [Fintype Ωin]

/-- A sampled SPIN encoder `E = I ∘ Π ∘ O`.

`outer ω₁` is the realized outer encoder `O`; `innerWt ω₂ u` is the weight
`wt(I(Π(u)))` of the routed and inner-encoded word.  Bundling the route and
the inner stage into a single weight function is faithful to the framework:
it uses `Π` and `I` only through that composite.

The two laws are separate fields, and all statements below are about their
independent product.  That is how `O ⫫ (Π, I)` is expressed. -/
structure Setup (M W Ωout Ωin : Type*) [Fintype M] [Fintype Ωout] [Fintype Ωin] where
  /-- Law of the outer setup. -/
  Pout : FinPMF Ωout
  /-- Law of the route and inner setup, independent of `Pout`. -/
  Pin : FinPMF Ωin
  outer : Ωout → M → W
  innerWt : Ωin → W → ℕ

namespace Setup

variable (S : Setup M W Ωout Ωin) (d : ℕ)

/-- `Z_d`: the number of nonzero messages whose encoding has weight at most `d`. -/
def Z (ω : Ωout × Ωin) : ℕ :=
  ((nonzeroMsgs M).filter (fun x => S.innerWt ω.2 (S.outer ω.1 x) ≤ d)).card

/-- `q_d(u)`: the conditional failure probability of a *fixed* outer word `u`,
over the route and inner setup only (paper: `eq:conditional-inner-failure`). -/
def qd (u : W) : ℝ := S.Pin.prob (fun ωin => S.innerWt ωin u ≤ d)

lemma qd_nonneg (u : W) : 0 ≤ S.qd d u := FinPMF.prob_nonneg _ _

lemma qd_le_one (u : W) : S.qd d u ≤ 1 := FinPMF.prob_le_one _ _

/-- The joint law of the whole construction. -/
def joint : FinPMF (Ωout × Ωin) := S.Pout.prod S.Pin

/-- **First-moment identity.**  Condition on the outer realization and use
independence of `O` from `(Π, I)`:

  `𝔼[Z_d] = 𝔼_O[ Σ_{x ≠ 0} q_d(O(x)) ]`. -/
theorem expect_Z_eq :
    S.joint.expect (fun ω => (S.Z d ω : ℝ))
      = S.Pout.expect (fun ω₁ => ∑ x ∈ nonzeroMsgs M, S.qd d (S.outer ω₁ x)) := by
  classical
  unfold joint
  rw [FinPMF.expect_prod]
  refine Finset.sum_congr rfl fun ω₁ _ => ?_
  refine congrArg _ ?_
  have hcard : (fun ω₂ : Ωin => ((S.Z d (ω₁, ω₂) : ℝ)))
      = fun ω₂ : Ωin =>
          ∑ x ∈ nonzeroMsgs M, (if S.innerWt ω₂ (S.outer ω₁ x) ≤ d then (1:ℝ) else 0) := by
    funext ω₂
    simp only [Z, Finset.card_filter, Nat.cast_sum]
    exact Finset.sum_congr rfl fun x _ => by split <;> simp
  show S.Pin.expect (fun ω₂ => ((S.Z d (ω₁, ω₂) : ℝ)))
      = ∑ x ∈ nonzeroMsgs M, S.qd d (S.outer ω₁ x)
  rw [hcard, FinPMF.expect_sum]
  exact Finset.sum_congr rfl fun x _ => (FinPMF.prob_eq_expect_indicator _ _).symm

/-- **Markov.**  The bad event is exactly `Z_d ≥ 1`. -/
theorem prob_bad_le_expect_Z :
    S.joint.prob (fun ω => 1 ≤ S.Z d ω) ≤ S.joint.expect (fun ω => (S.Z d ω : ℝ)) :=
  FinPMF.markov_one _ _

section Classes

variable {T : Type*} [Fintype T] [DecidableEq T] (c : W → T)

/-- `A_τ^out`: the expected number of nonzero messages whose outer word lies in
class `τ` (paper: `eq:outer-class-multiplicity`). -/
def Aout (τ : T) : ℝ :=
  S.Pout.expect
    (fun ω₁ => (((nonzeroMsgs M).filter (fun x => c (S.outer ω₁ x) = τ)).card : ℝ))

lemma Aout_nonneg (τ : T) : 0 ≤ S.Aout c τ :=
  FinPMF.expect_nonneg _ fun _ => by positivity

/-- **Theorem (class-indexed first-moment envelope).**

If `Q` is a certified envelope for the conditional failure probability that is
uniform over each class (`eq:class-failure-envelope`), then

  `𝔼[Z_d] ≤ Σ_τ A_τ^out · Q_τ`.

The class map need not determine `q_d` exactly; it need only expose enough
information for the uniform bound. -/
theorem expect_Z_le_class_sum (Q : T → ℝ) (hQ : ∀ u : W, S.qd d u ≤ Q (c u)) :
    S.joint.expect (fun ω => (S.Z d ω : ℝ)) ≤ ∑ τ, S.Aout c τ * Q τ := by
  classical
  rw [S.expect_Z_eq d]
  have hstep : ∀ ω₁ : Ωout,
      (∑ x ∈ nonzeroMsgs M, S.qd d (S.outer ω₁ x))
        ≤ ∑ τ, (((nonzeroMsgs M).filter (fun x => c (S.outer ω₁ x) = τ)).card : ℝ) * Q τ := by
    intro ω₁
    calc ∑ x ∈ nonzeroMsgs M, S.qd d (S.outer ω₁ x)
        ≤ ∑ x ∈ nonzeroMsgs M, Q (c (S.outer ω₁ x)) :=
          Finset.sum_le_sum fun x _ => hQ _
      _ = ∑ τ, ∑ x ∈ (nonzeroMsgs M).filter (fun x => c (S.outer ω₁ x) = τ),
              Q (c (S.outer ω₁ x)) :=
          (Finset.sum_fiberwise_of_maps_to (fun x _ => Finset.mem_univ _) _).symm
      _ = ∑ τ, (((nonzeroMsgs M).filter (fun x => c (S.outer ω₁ x) = τ)).card : ℝ) * Q τ := by
          refine Finset.sum_congr rfl fun τ _ => ?_
          rw [Finset.sum_congr rfl (fun x hx => by
                simp only [Finset.mem_filter] at hx
                rw [hx.2]), Finset.sum_const, nsmul_eq_mul]
  calc S.Pout.expect (fun ω₁ => ∑ x ∈ nonzeroMsgs M, S.qd d (S.outer ω₁ x))
      ≤ S.Pout.expect (fun ω₁ =>
          ∑ τ, (((nonzeroMsgs M).filter (fun x => c (S.outer ω₁ x) = τ)).card : ℝ) * Q τ) :=
        FinPMF.expect_mono _ hstep
    _ = ∑ τ, S.Aout c τ * Q τ := by
        rw [FinPMF.expect_sum]
        refine Finset.sum_congr rfl fun τ _ => ?_
        unfold Aout FinPMF.expect
        rw [Finset.sum_mul]
        exact Finset.sum_congr rfl fun ω₁ _ => by ring

/-- The form consumed downstream: a bound on the probability of the bad event. -/
theorem prob_bad_le_class_sum (Q : T → ℝ) (hQ : ∀ u : W, S.qd d u ≤ Q (c u)) :
    S.joint.prob (fun ω => 1 ≤ S.Z d ω) ≤ ∑ τ, S.Aout c τ * Q τ :=
  le_trans (S.prob_bad_le_expect_Z d) (S.expect_Z_le_class_sum d c Q hQ)

end Classes

end Setup

end Spin

namespace Spin

open Finset

namespace Setup

variable {M W Ωout Ωin : Type*}
variable [Fintype M] [DecidableEq M] [Zero M]
variable [Fintype Ωout] [Fintype Ωin]
variable (S : Setup M W Ωout Ωin) (d : ℕ)

section Occupation

variable (occ : Ωout → M → ℕ)

/-- `Z_{d,Q}`: nonzero messages with exactly `Q` active outer blocks whose
encoding has weight at most `d`.  In the structured proof `occ ω₁ x` is the
number of active outer positions of the outer word `O(x)`. -/
def ZQ (Q : ℕ) (ω : Ωout × Ωin) : ℕ :=
  ((nonzeroMsgs M).filter
    (fun x => occ ω.1 x = Q ∧ S.innerWt ω.2 (S.outer ω.1 x) ≤ d)).card

/-- The occupation refinement is a partition of `Z_d`. -/
lemma Z_eq_sum_ZQ (Lmax : ℕ) (hocc : ∀ ω₁ x, x ∈ nonzeroMsgs M → occ ω₁ x ∈ Ico 1 (Lmax + 1))
    (ω : Ωout × Ωin) :
    S.Z d ω = ∑ Q ∈ Ico 1 (Lmax + 1), S.ZQ d occ Q ω := by
  classical
  unfold Z ZQ
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun x => occ ω.1 x) (t := Ico 1 (Lmax + 1))
    (fun x hx => hocc ω.1 x (Finset.mem_of_mem_filter x hx))]
  refine Finset.sum_congr rfl fun Q _ => ?_
  congr 1
  ext x
  simp only [Finset.mem_filter]
  tauto

/-- `𝔼[Z_d] = Σ_Q 𝔼[Z_{d,Q}]` under any law. -/
lemma expect_Z_eq_sum_expect_ZQ (P : FinPMF (Ωout × Ωin)) (Lmax : ℕ)
    (hocc : ∀ ω₁ x, x ∈ nonzeroMsgs M → occ ω₁ x ∈ Ico 1 (Lmax + 1)) :
    P.expect (fun ω => (S.Z d ω : ℝ))
      = ∑ Q ∈ Ico 1 (Lmax + 1), P.expect (fun ω => (S.ZQ d occ Q ω : ℝ)) := by
  classical
  rw [← FinPMF.expect_sum]
  refine Finset.sum_congr rfl fun ω _ => ?_
  refine congrArg _ ?_
  show ((S.Z d ω : ℝ)) = ∑ i ∈ Ico 1 (Lmax + 1), ((S.ZQ d occ i ω : ℝ))
  rw [S.Z_eq_sum_ZQ d occ Lmax hocc ω]
  push_cast
  rfl

/-- **Linkage.**  The whole structured argument in one inequality:

  `Pr[bad] ≤ Pr[¬𝒢] + Σ_{Q=1}^{L} 𝔼[Z_{d,Q} | 𝒢]`.

The conditioning is on an event `𝒢` of the *outer setup only*, and the proof
goes through `condition_prod_left`: the route and inner setup remain
independent of the conditioned outer encoder.  This is the step the paper
records as "the pair `(Π, I)` remains independent of the conditioned outer
encoder", and it is what licenses using `𝒢_b`-conditional regime bounds in
the final union. -/
theorem prob_bad_le_cond_sum (G : Ωout → Prop) [DecidablePred G]
    (hG : 0 < S.Pout.prob G) (Lmax : ℕ)
    (hocc : ∀ ω₁ x, x ∈ nonzeroMsgs M → occ ω₁ x ∈ Ico 1 (Lmax + 1)) :
    S.joint.prob (fun ω => 1 ≤ S.Z d ω)
      ≤ S.Pout.prob (fun ω₁ => ¬ G ω₁)
        + ∑ Q ∈ Ico 1 (Lmax + 1),
            ((S.Pout.condition G hG).prod S.Pin).expect
              (fun ω => (S.ZQ d occ Q ω : ℝ)) := by
  classical
  have hGjoint : 0 < S.joint.prob (fun ω : Ωout × Ωin => G ω.1) := by
    unfold joint
    rwa [FinPMF.prob_prod_left]
  have hsplit := FinPMF.prob_le_compl_add_condition S.joint
    (fun ω => 1 ≤ S.Z d ω) (fun ω : Ωout × Ωin => G ω.1) hGjoint
  have hcompl : S.joint.prob (fun ω : Ωout × Ωin => ¬ G ω.1)
      = S.Pout.prob (fun ω₁ => ¬ G ω₁) := by
    unfold joint
    exact FinPMF.prob_prod_left S.Pout S.Pin (fun ω₁ => ¬ G ω₁)
  have hcond : S.joint.condition (fun ω : Ωout × Ωin => G ω.1) hGjoint
      = (S.Pout.condition G hG).prod S.Pin := by
    unfold joint at hGjoint ⊢
    exact FinPMF.condition_prod_left _ _ _ hG hGjoint
  rw [hcompl, hcond] at hsplit
  have hkey : ((S.Pout.condition G hG).prod S.Pin).prob (fun ω => 1 ≤ S.Z d ω)
      ≤ ∑ Q ∈ Ico 1 (Lmax + 1), ((S.Pout.condition G hG).prod S.Pin).expect
          (fun ω => (S.ZQ d occ Q ω : ℝ)) :=
    le_trans (FinPMF.markov_one _ _)
      (le_of_eq (S.expect_Z_eq_sum_expect_ZQ d occ _ Lmax hocc))
  linarith

end Occupation

end Setup

end Spin
