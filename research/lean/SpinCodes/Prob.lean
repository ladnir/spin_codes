/-
Finite probability layer for the SPIN formalization.

Every source of randomness in the structured SPIN construction is finite:
two BA permutations in `S_b`, one block shuffle per outer position, one
region shuffle per region, and one transvection per IMT step.  We therefore
model probability with explicit mass functions on `Fintype`s rather than
with `MeasureTheory`.  This keeps expectations as `Finset` sums, which is
what the first-moment argument manipulates.
-/
import Mathlib

namespace Spin

open Finset

/-- A probability mass function on a finite type. -/
structure FinPMF (Ω : Type*) [Fintype Ω] where
  p : Ω → ℝ
  nonneg : ∀ ω, 0 ≤ p ω
  total : ∑ ω, p ω = 1

namespace FinPMF

variable {Ω Ω₁ Ω₂ : Type*} [Fintype Ω] [Fintype Ω₁] [Fintype Ω₂]

/-- Expectation of a real-valued random variable. -/
def expect (P : FinPMF Ω) (f : Ω → ℝ) : ℝ := ∑ ω, P.p ω * f ω

/-- Probability of a decidable event. -/
def prob (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] : ℝ :=
  ∑ ω ∈ univ.filter s, P.p ω

lemma prob_eq_expect_indicator (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] :
    P.prob s = P.expect (fun ω => if s ω then 1 else 0) := by
  classical
  unfold prob expect
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun ω _ => ?_
  by_cases h : s ω <;> simp [h]

lemma expect_nonneg (P : FinPMF Ω) {f : Ω → ℝ} (hf : ∀ ω, 0 ≤ f ω) :
    0 ≤ P.expect f :=
  Finset.sum_nonneg fun ω _ => mul_nonneg (P.nonneg ω) (hf ω)

lemma expect_mono (P : FinPMF Ω) {f g : Ω → ℝ} (h : ∀ ω, f ω ≤ g ω) :
    P.expect f ≤ P.expect g :=
  Finset.sum_le_sum fun ω _ => mul_le_mul_of_nonneg_left (h ω) (P.nonneg ω)

lemma expect_const (P : FinPMF Ω) (c : ℝ) : P.expect (fun _ => c) = c := by
  unfold expect
  rw [← Finset.sum_mul, P.total, one_mul]

lemma expect_add (P : FinPMF Ω) (f g : Ω → ℝ) :
    P.expect (fun ω => f ω + g ω) = P.expect f + P.expect g := by
  unfold expect
  rw [← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun ω _ => by ring

lemma expect_sum {ι : Type*} (P : FinPMF Ω) (s : Finset ι) (f : ι → Ω → ℝ) :
    P.expect (fun ω => ∑ i ∈ s, f i ω) = ∑ i ∈ s, P.expect (f i) := by
  unfold expect
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun ω _ => by rw [Finset.mul_sum]

lemma prob_nonneg (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] : 0 ≤ P.prob s :=
  Finset.sum_nonneg fun ω _ => P.nonneg ω

lemma prob_le_one (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] : P.prob s ≤ 1 := by
  unfold prob
  rw [← P.total]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
    (fun ω _ _ => P.nonneg ω)

lemma prob_mono (P : FinPMF Ω) {s t : Ω → Prop} [DecidablePred s] [DecidablePred t]
    (h : ∀ ω, s ω → t ω) : P.prob s ≤ P.prob t := by
  unfold prob
  refine Finset.sum_le_sum_of_subset_of_nonneg ?_ (fun ω _ _ => P.nonneg ω)
  intro ω hω
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hω ⊢
  exact h ω hω

/-- Complement probability. -/
lemma prob_compl (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] :
    P.prob (fun ω => ¬ s ω) = 1 - P.prob s := by
  classical
  have h := Finset.sum_filter_add_sum_filter_not (univ : Finset Ω) s P.p
  rw [P.total] at h
  unfold prob
  linarith

/-- Equal events have equal probability (used to move across coercions). -/
lemma prob_congr (P : FinPMF Ω) {s t : Ω → Prop} [DecidablePred s] [DecidablePred t]
    (h : ∀ ω, s ω ↔ t ω) : P.prob s = P.prob t :=
  le_antisymm (prob_mono P fun ω => (h ω).mp) (prob_mono P fun ω => (h ω).mpr)

/-- The union bound over a finite index set. -/
lemma prob_biUnion_le {ι : Type*} (P : FinPMF Ω) (s : Finset ι) (E : ι → Ω → Prop)
    [∀ i, DecidablePred (E i)] [DecidablePred fun ω => ∃ i ∈ s, E i ω] :
    P.prob (fun ω => ∃ i ∈ s, E i ω) ≤ ∑ i ∈ s, P.prob (E i) := by
  classical
  rw [prob_eq_expect_indicator,
      Finset.sum_congr rfl (fun i (_ : i ∈ s) => prob_eq_expect_indicator P (E i)),
      ← expect_sum]
  refine expect_mono P fun ω => ?_
  by_cases h : ∃ i ∈ s, E i ω
  · obtain ⟨i, hi, hEi⟩ := h
    have hone : (if (∃ i ∈ s, E i ω) then (1:ℝ) else 0) = 1 := if_pos ⟨i, hi, hEi⟩
    have hEi' : (if E i ω then (1:ℝ) else 0) = 1 := if_pos hEi
    rw [hone]
    calc (1:ℝ) = (if E i ω then (1:ℝ) else 0) := hEi'.symm
      _ ≤ ∑ j ∈ s, (if E j ω then (1:ℝ) else 0) :=
          Finset.single_le_sum (f := fun j => if E j ω then (1:ℝ) else 0)
            (fun j _ => by positivity) hi
  · rw [if_neg h]
    exact Finset.sum_nonneg fun j _ => by positivity

/-- Subadditivity for two events. -/
lemma prob_or_le (P : FinPMF Ω) (s t : Ω → Prop) [DecidablePred s] [DecidablePred t]
    [DecidablePred fun ω => s ω ∨ t ω] :
    P.prob (fun ω => s ω ∨ t ω) ≤ P.prob s + P.prob t := by
  classical
  rw [prob_eq_expect_indicator, prob_eq_expect_indicator, prob_eq_expect_indicator,
      ← expect_add]
  refine expect_mono P fun ω => ?_
  by_cases hs : s ω <;> by_cases ht : t ω <;> simp [hs, ht]

/-- Markov's inequality.  The form used throughout: a nonnegative random
variable exceeds `t > 0` with probability at most `𝔼[f] / t`. -/
theorem markov (P : FinPMF Ω) {f : Ω → ℝ} (hf : ∀ ω, 0 ≤ f ω) {t : ℝ} (ht : 0 < t) :
    P.prob (fun ω => t ≤ f ω) ≤ P.expect f / t := by
  classical
  rw [le_div_iff₀ ht]
  unfold prob expect
  rw [Finset.sum_mul]
  calc ∑ ω ∈ univ.filter (fun ω => t ≤ f ω), P.p ω * t
      ≤ ∑ ω ∈ univ.filter (fun ω => t ≤ f ω), P.p ω * f ω := by
        refine Finset.sum_le_sum fun ω hω => ?_
        simp only [Finset.mem_filter] at hω
        exact mul_le_mul_of_nonneg_left hω.2 (P.nonneg ω)
    _ ≤ ∑ ω, P.p ω * f ω :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
          (fun ω _ _ => mul_nonneg (P.nonneg ω) (hf ω))

/-- Markov for a counting variable: the probability that a nonnegative
integer-valued variable is at least one is bounded by its mean.  This is the
step `Pr[Z_d ≥ 1] ≤ 𝔼[Z_d]` of the first-moment framework. -/
theorem markov_one (P : FinPMF Ω) (Z : Ω → ℕ) :
    P.prob (fun ω => 1 ≤ Z ω) ≤ P.expect (fun ω => (Z ω : ℝ)) := by
  classical
  have h := markov P (f := fun ω => (Z ω : ℝ)) (fun ω => by positivity) (t := 1) one_pos
  rw [div_one] at h
  refine le_trans (le_of_eq (prob_congr P (fun ω => ?_))) h
  exact (Nat.one_le_cast (α := ℝ)).symm

/-- The independent product of two finite mass functions.  This is how the
paper's independence of the outer encoder from `(Π, I)` is expressed. -/
def prod (P : FinPMF Ω₁) (Q : FinPMF Ω₂) : FinPMF (Ω₁ × Ω₂) where
  p := fun ω => P.p ω.1 * Q.p ω.2
  nonneg := fun ω => mul_nonneg (P.nonneg ω.1) (Q.nonneg ω.2)
  total := by
    rw [Fintype.sum_prod_type]
    simp_rw [← Finset.mul_sum, Q.total, mul_one, P.total]

/-- Fubini for the product mass function. -/
lemma expect_prod (P : FinPMF Ω₁) (Q : FinPMF Ω₂) (f : Ω₁ × Ω₂ → ℝ) :
    (P.prod Q).expect f = P.expect (fun ω₁ => Q.expect (fun ω₂ => f (ω₁, ω₂))) := by
  unfold expect prod
  rw [Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun ω₁ _ => ?_
  rw [Finset.mul_sum]
  exact Finset.sum_congr rfl fun ω₂ _ => by ring

/-- Conditioning on an event of positive probability.  The framework theorem
is stated for an arbitrary outer law, so conditioning on a setup-only event
`𝒢` and re-applying it is exactly the remark following
`thm:class-first-moment` in the paper. -/
noncomputable def condition (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] (h : 0 < P.prob s) :
    FinPMF Ω where
  p := fun ω => if s ω then P.p ω / P.prob s else 0
  nonneg := fun ω => by
    by_cases hs : s ω
    · simp only [if_pos hs]; exact div_nonneg (P.nonneg ω) h.le
    · simp [hs]
  total := by
    classical
    rw [← Finset.sum_filter, ← Finset.sum_div]
    exact div_self h.ne'

lemma condition_p_apply (P : FinPMF Ω) (s : Ω → Prop) [DecidablePred s] (h : 0 < P.prob s)
    (ω : Ω) : (P.condition s h).p ω = if s ω then P.p ω / P.prob s else 0 := rfl


/-- Two mass functions with the same weights are equal. -/
@[ext] lemma ext {P Q : FinPMF Ω} (h : P.p = Q.p) : P = Q := by
  cases P; cases Q; simp_all

lemma prob_and_add_prob_and_not (P : FinPMF Ω) (s G : Ω → Prop)
    [DecidablePred s] [DecidablePred G] :
    P.prob (fun ω => s ω ∧ G ω) + P.prob (fun ω => s ω ∧ ¬ G ω) = P.prob s := by
  classical
  unfold prob
  rw [← Finset.sum_filter_add_sum_filter_not (univ.filter s) G P.p]
  congr 1 <;> · congr 1; ext ω; simp [Finset.mem_filter, and_comm]

lemma condition_prob (P : FinPMF Ω) (s G : Ω → Prop) [DecidablePred s] [DecidablePred G]
    (h : 0 < P.prob G) :
    (P.condition G h).prob s = P.prob (fun ω => s ω ∧ G ω) / P.prob G := by
  classical
  have h1 : (P.condition G h).prob s
      = ∑ x ∈ univ.filter s, (if G x then P.p x / P.prob G else 0) := rfl
  rw [h1, ← Finset.sum_filter, Finset.filter_filter, ← Finset.sum_div]
  rfl

/-- Splitting on a conditioning event: the paper's "condition on any outer
satisfying `𝒢_b`" step.  The conditional probability is *not* multiplied back
by `Pr[𝒢]`; we only use `Pr[𝒢] ≤ 1`, which is the direction the union needs. -/
lemma prob_le_compl_add_condition (P : FinPMF Ω) (s G : Ω → Prop)
    [DecidablePred s] [DecidablePred G] (h : 0 < P.prob G) :
    P.prob s ≤ P.prob (fun ω => ¬ G ω) + (P.condition G h).prob s := by
  classical
  rw [condition_prob P s G h]
  have hsplit := prob_and_add_prob_and_not P s G
  have h1 : P.prob (fun ω => s ω ∧ ¬ G ω) ≤ P.prob (fun ω => ¬ G ω) :=
    prob_mono P fun ω hω => hω.2
  have h2 : P.prob (fun ω => s ω ∧ G ω)
      ≤ P.prob (fun ω => s ω ∧ G ω) / P.prob G := by
    rw [le_div_iff₀ h]
    nlinarith [prob_le_one P G, prob_nonneg P (fun ω => s ω ∧ G ω)]
  linarith

lemma prob_prod_left (P : FinPMF Ω₁) (Q : FinPMF Ω₂) (G : Ω₁ → Prop) [DecidablePred G] :
    (P.prod Q).prob (fun ω => G ω.1) = P.prob G := by
  classical
  rw [prob_eq_expect_indicator, expect_prod]
  have hc : ∀ ω₁ : Ω₁,
      Q.expect (fun _ : Ω₂ => (if G ω₁ then (1:ℝ) else 0)) = if G ω₁ then 1 else 0 :=
    fun ω₁ => expect_const Q _
  simp_rw [hc]
  exact (prob_eq_expect_indicator P G).symm

/-- Conditioning a product law on an event of the first coordinate conditions
only the first factor: `(O, (Π, I)) | 𝒢 = (O | 𝒢, (Π, I))`.  This is the
formal content of "the pair `(Π, I)` remains independent of the conditioned
outer encoder". -/
lemma condition_prod_left (P : FinPMF Ω₁) (Q : FinPMF Ω₂) (G : Ω₁ → Prop)
    [DecidablePred G] (h : 0 < P.prob G)
    (h' : 0 < (P.prod Q).prob (fun ω => G ω.1)) :
    (P.prod Q).condition (fun ω => G ω.1) h' = (P.condition G h).prod Q := by
  classical
  have hpp : (P.prod Q).prob (fun ω => G ω.1) = P.prob G := prob_prod_left P Q G
  ext ω
  show (if G ω.1 then (P.p ω.1 * Q.p ω.2) / (P.prod Q).prob (fun ω => G ω.1) else 0)
      = (if G ω.1 then P.p ω.1 / P.prob G else 0) * Q.p ω.2
  rw [hpp]
  by_cases hG : G ω.1
  · simp only [if_pos hG]
    ring
  · simp [hG]


/-- The uniform law on a nonempty finite type. -/
noncomputable def uniform (Ω : Type*) [Fintype Ω] [Nonempty Ω] : FinPMF Ω where
  p := fun _ => 1 / (Fintype.card Ω : ℝ)
  nonneg := fun _ => by positivity
  total := by
    have hpos : (0 : ℝ) < (Fintype.card Ω : ℝ) := by
      exact_mod_cast Fintype.card_pos
    rw [Finset.sum_const, nsmul_eq_mul, Finset.card_univ]
    field_simp

lemma uniform_prob (Ω : Type*) [Fintype Ω] [Nonempty Ω] (s : Ω → Prop) [DecidablePred s] :
    (uniform Ω).prob s = ((univ.filter s).card : ℝ) / (Fintype.card Ω : ℝ) := by
  unfold prob uniform
  simp only
  rw [Finset.sum_const, nsmul_eq_mul]
  ring

/-- Expectation of a function of a finite-valued random variable. -/
lemma expect_comp {ι : Type*} [DecidableEq ι] (P : FinPMF Ω) (X : Ω → ι) (t : Finset ι)
    (hX : ∀ ω, X ω ∈ t) (g : ι → ℝ) :
    P.expect (fun ω => g (X ω)) = ∑ i ∈ t, P.prob (fun ω => X ω = i) * g i := by
  classical
  unfold expect prob
  rw [← Finset.sum_fiberwise_of_maps_to (g := X) (t := t) (fun ω _ => hX ω)]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl fun ω hω => ?_
  simp only [Finset.mem_filter] at hω
  show P.p ω * g (X ω) = P.p ω * g i
  rw [hω.2]

/-- Conditioning the product law on the first coordinate: the probability of a
joint event is the expectation over the first factor of the conditional
probability under the second. -/
lemma prob_prod_eq_expect (P : FinPMF Ω₁) (Q : FinPMF Ω₂) (R : Ω₁ → Ω₂ → Prop)
    [∀ a b, Decidable (R a b)] :
    (P.prod Q).prob (fun ω => R ω.1 ω.2) = P.expect (fun ω₁ => Q.prob (R ω₁)) := by
  rw [prob_eq_expect_indicator, expect_prod]
  unfold expect
  refine Finset.sum_congr rfl fun ω₁ _ => ?_
  have h : (∑ ω₂ : Ω₂, Q.p ω₂ * (if R ω₁ ω₂ then (1:ℝ) else 0)) = Q.prob (R ω₁) :=
    (prob_eq_expect_indicator Q (R ω₁)).symm
  simpa using congrArg (fun z : ℝ => P.p ω₁ * z) h

end FinPMF

end Spin
