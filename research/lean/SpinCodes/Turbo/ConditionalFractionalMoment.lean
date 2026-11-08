import SpinCodes.FiniteLaw
import SpinCodes.Turbo.Occupation22

/-!
# Route-conditioned probability capping before fractional powers

The retained K18 replay conditions on potential-slot step counts. Its
conditional first-moment estimate must hold in every fiber of an explicit
disintegration of the original setup law. An unconditional mean does not
supply that hypothesis. The active-set union stays outside the power.

Old counterpart: `Setup.prob_bad_le_cond_sum` conditions an ordinary first
moment on an outer event. There is no fractional counterpart in the old
Structured proof. This module isolates the additional probability argument;
it does not authenticate physical matrices or their regional placement.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset

/-- Capping uses the probability's upper bound one *before* the fractional
power; no independence among the counted messages is assumed. -/
theorem prob_count_le_fractional_mean {Ω : Type*} [Fintype Ω]
    (P : FinPMF Ω) (Z : Ω → ℕ) {A α : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : P.expect (fun ω => (Z ω : ℝ)) ≤ A) :
    P.prob (fun ω => 1 ≤ Z ω) ≤ A ^ α := by
  have hp0 := P.prob_nonneg (fun ω => 1 ≤ Z ω)
  have hp1 := P.prob_le_one (fun ω => 1 ≤ Z ω)
  have hpA := (P.markov_one Z).trans hmean
  calc
    P.prob (fun ω => 1 ≤ Z ω) =
        (P.prob (fun ω => 1 ≤ Z ω)) ^ (1 : ℝ) := (Real.rpow_one _).symm
    _ ≤ (P.prob (fun ω => 1 ≤ Z ω)) ^ α :=
      Real.rpow_le_rpow_of_exponent_ge' hp0 hp1 hα0.le hα1
    _ ≤ A ^ α := Real.rpow_le_rpow hp0 hpA hα0.le

/-- A supplied conditional law retains the routing average outside the power. -/
theorem bind_prob_count_le_fractional_mean
    {J Ω : Type*} [Fintype J] [Fintype Ω]
    (P : FinPMF J) (Q : J → FinPMF Ω) (Z : Ω → ℕ)
    (A : J → ℝ) {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ j, (Q j).expect (fun ω => (Z ω : ℝ)) ≤ A j) :
    (P.bind Q).prob (fun ω => 1 ≤ Z ω) ≤
      P.expect (fun j => (A j) ^ α) := by
  rw [FinPMF.bind_prob]
  exact P.expect_mono (fun j => prob_count_le_fractional_mean
    (Q j) Z hα0 hα1 (hmean j))

/-- Each active-set event is capped separately; the finite union is outside
all fractional powers and all routing expectations. -/
theorem bind_prob_union_count_le_fractional_mean
    {J Ω I : Type*} [Fintype J] [Fintype Ω]
    (P : FinPMF J) (Q : J → FinPMF Ω) (sets : Finset I)
    (Z : I → Ω → ℕ) (A : I → J → ℝ)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ i ∈ sets, ∀ j,
      (Q j).expect (fun ω => (Z i ω : ℝ)) ≤ A i j) :
    (P.bind Q).prob (fun ω => ∃ i ∈ sets, 1 ≤ Z i ω) ≤
      ∑ i ∈ sets, P.expect (fun j => (A i j) ^ α) := by
  apply (FinPMF.prob_biUnion_le _ sets (fun i ω => 1 ≤ Z i ω)).trans
  exact sum_le_sum (fun i hi => bind_prob_count_le_fractional_mean
    P Q (Z i) (A i) hα0 hα1 (hmean i hi))

/-- Bad-message count for one fixed active group set in the actual encoder. -/
def turboActiveBadCount22 (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (S : Finset (Fin (R * 16)))
    (ω : TurboOuterSeed22 R × RouteInnerSeed22 R) : ℕ :=
  ((turboActiveMessages22 S).filter (fun msg =>
    (turboSetup22 R points).innerWt ω.2
      ((turboSetup22 R points).outer ω.1 msg) ≤ d)).card

/-- Conditional whole-encoder failure bound under the original setup law.
The disintegration and every fiberwise active-set mean are explicit premises.
In K18, `J` should record potential-slot counts and `A` is the chronological
matrix bound after outer domination. This theorem does not derive either. -/
theorem turboSetup22_bad_prob_le_conditional_fractional
    {J : Type*} [Fintype J]
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (P : FinPMF J)
    (Q : J → FinPMF (TurboOuterSeed22 R × RouteInnerSeed22 R))
    (hlaw : (turboSetup22 R points).joint = P.bind Q)
    (sets : Finset (Finset (Fin (R * 16))))
    (hcover : ∀ msg : TurboMessage22 R, msg ≠ 0 →
      turboActivePositions22 msg ∈ sets)
    (A : Finset (Fin (R * 16)) → J → ℝ)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ S ∈ sets, ∀ j,
      (Q j).expect (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) ≤ A S j) :
    (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      ∑ S ∈ sets, P.expect (fun j => (A S j) ^ α) := by
  rw [hlaw]
  apply (FinPMF.prob_mono _ (fun ω hω => ?_)).trans
    (bind_prob_union_count_le_fractional_mean P Q sets
      (turboActiveBadCount22 R points d) A hα0 hα1 hmean)
  obtain ⟨msg, hm⟩ := card_pos.mp hω
  have hn : msg ≠ 0 := (mem_filter.mp (mem_filter.mp hm).1).2
  refine ⟨turboActivePositions22 msg, hcover msg hn, ?_⟩
  apply card_pos.mpr
  exact ⟨msg, mem_filter.mpr ⟨mem_filter.mpr ⟨mem_univ _, rfl⟩,
    (mem_filter.mp hm).2⟩⟩

/-- The retained K18 application conditions separately for each active set.
Its potential-slot counts depend on `S`; this interface does not require a
single joint conditioning that reveals the placements of other active sets. -/
theorem turboSetup22_bad_prob_le_activeSet_conditional_fractional
    {J : Type*} [Fintype J]
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (sets : Finset (Finset (Fin (R * 16))))
    (hcover : ∀ msg : TurboMessage22 R, msg ≠ 0 →
      turboActivePositions22 msg ∈ sets)
    (P : Finset (Fin (R * 16)) → FinPMF J)
    (Q : Finset (Fin (R * 16)) → J →
      FinPMF (TurboOuterSeed22 R × RouteInnerSeed22 R))
    (hlaw : ∀ S ∈ sets, (turboSetup22 R points).joint = (P S).bind (Q S))
    (A : Finset (Fin (R * 16)) → J → ℝ)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ S ∈ sets, ∀ j,
      (Q S j).expect (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) ≤ A S j) :
    (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      ∑ S ∈ sets, (P S).expect (fun j => (A S j) ^ α) := by
  have hsub : (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      (turboSetup22 R points).joint.prob
        (fun ω => ∃ S ∈ sets, 1 ≤ turboActiveBadCount22 R points d S ω) := by
    apply FinPMF.prob_mono
    intro ω hω
    obtain ⟨msg, hm⟩ := card_pos.mp hω
    have hn : msg ≠ 0 := (mem_filter.mp (mem_filter.mp hm).1).2
    refine ⟨turboActivePositions22 msg, hcover msg hn, ?_⟩
    apply card_pos.mpr
    exact ⟨msg, mem_filter.mpr ⟨mem_filter.mpr ⟨mem_univ _, rfl⟩,
      (mem_filter.mp hm).2⟩⟩
  apply hsub.trans
  apply (FinPMF.prob_biUnion_le _ sets
    (fun S ω => 1 ≤ turboActiveBadCount22 R points d S ω)).trans
  apply sum_le_sum
  intro S hS
  rw [hlaw S hS]
  exact bind_prob_count_le_fractional_mean (P S) (Q S)
    (turboActiveBadCount22 R points d S) (A S) hα0 hα1 (hmean S hS)

#print axioms turboSetup22_bad_prob_le_activeSet_conditional_fractional
#print axioms bind_prob_count_le_fractional_mean
#print axioms bind_prob_union_count_le_fractional_mean
#print axioms turboSetup22_bad_prob_le_conditional_fractional
end Spin.Turbo
