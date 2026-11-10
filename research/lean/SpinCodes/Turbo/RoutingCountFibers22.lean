import SpinCodes.Turbo.RoutingCountDisintegration22

/-!
# Fibers of the actual S-specific routing-count statistic

Finite conditioning is performed on the original shared setup, not a
resampled setup per message. Zero-probability fibers use a harmless fallback;
their mixing weight is zero. The count marginal is the already checked
regional profile law, with all 64 independent regional count profiles.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset

private theorem finite_map_mass_eq_fiber {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) (j : J) :
    (P.map f).p j = P.prob (fun ω => f ω = j) := by
  simp only [FinPMF.map_p, FinPMF.prob, sum_filter]

/-- Conditional law on one finite statistic; unused fibers retain P. -/
def finiteStatisticFiber {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) (j : J) : FinPMF Ω :=
  if hp : 0 < P.prob (fun ω => f ω = j) then
    P.condition (fun ω => f ω = j) hp else P

theorem finiteStatisticFiber_disintegration {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) :
    P = (P.map f).bind (finiteStatisticFiber P f) := by
  classical
  apply FinPMF.ext
  funext ω
  have hterm (j : J) : (P.map f).p j * (finiteStatisticFiber P f j).p ω =
      if f ω = j then P.p ω else 0 := by
    rw [finite_map_mass_eq_fiber]
    by_cases hp : 0 < P.prob (fun x => f x = j)
    · rw [finiteStatisticFiber, dite_eq_left hp, FinPMF.condition_p_apply]
      split_ifs with he
      · field_simp [hp.ne']
      · simp
    · have hz : P.prob (fun x => f x = j) = 0 :=
        le_antisymm (le_of_not_gt hp) (P.prob_nonneg _)
      rw [hz, zero_mul]
      by_cases he : f ω = j
      · rw [ite_eq_left he]
        have hle : P.p ω ≤ (P.map f).p j := by
          rw [FinPMF.map_p]
          have h := single_le_sum
            (f := fun x => if f x = j then P.p x else 0)
            (fun x (_ : x ∈ (univ : Finset Ω)) => by
              split_ifs <;> first | exact P.nonneg x | exact le_rfl)
            (mem_univ ω)
          simpa only [ite_eq_left he] using h
        rw [finite_map_mass_eq_fiber, hz] at hle
        exact (le_antisymm hle (P.nonneg ω)).symm
      · rw [ite_eq_right he]
  rw [FinPMF.bind_p]
  simp_rw [hterm]
  simp

def turboSetupCountFiber22 (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R) :
    FinPMF (TurboOuterSeed22 R × RouteInnerSeed22 R) :=
  finiteStatisticFiber (turboSetup22 R points).joint
    (turboSetupPotentialCounts22 S) counts

/-- Actual original-law disintegration, separately for each active set S.
No matrix estimate or conditional independence is assumed in this identity. -/
theorem turboSetup22_count_disintegration (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) :
    (turboSetup22 R points).joint =
      (turboPotentialCountLaw22 R S).bind (turboSetupCountFiber22 R points S) := by
  rw [← turboSetupPotentialCounts22_law R points S]
  exact finiteStatisticFiber_disintegration _ _

/-- Zero-weight profiles need no conditional estimate. -/
theorem bind_prob_count_le_fractional_mean_on_support
    {J Ω : Type*} [Fintype J] [Fintype Ω]
    (P : FinPMF J) (Q : J → FinPMF Ω) (Z : Ω → ℕ)
    (A : J → ℝ) {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ j, 0 < P.p j →
      (Q j).expect (fun ω => (Z ω : ℝ)) ≤ A j) :
    (P.bind Q).prob (fun ω => 1 ≤ Z ω) ≤
      P.expect (fun j => (A j) ^ α) := by
  rw [FinPMF.bind_prob]
  unfold FinPMF.expect
  apply sum_le_sum
  intro j _
  by_cases hp : 0 < P.p j
  · exact mul_le_mul_of_nonneg_left
      (prob_count_le_fractional_mean (Q j) Z hα0 hα1 (hmean j hp)) (P.nonneg j)
  · have hz : P.p j = 0 := le_antisymm (le_of_not_gt hp) (P.nonneg j)
    simp only [hz, zero_mul, le_refl]

/-- Original setup failure bounded by S-specific routing-conditioned moments.
Only the conditional moment estimate on feasible profiles remains a premise;
the exact original-law disintegration and profile weights are proved above. -/
theorem turboSetup22_bad_prob_le_count_fractional
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (sets : Finset (Finset (Fin (R * 16))))
    (hcover : ∀ msg : TurboMessage22 R, msg ≠ 0 →
      turboActivePositions22 msg ∈ sets)
    (A : Finset (Fin (R * 16)) → TurboPotentialCounts22 R → ℝ)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hmean : ∀ S ∈ sets, ∀ counts,
      0 < (turboPotentialCountLaw22 R S).p counts →
      (turboSetupCountFiber22 R points S counts).expect
        (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) ≤ A S counts) :
    (turboSetup22 R points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 R points).Z d ω) ≤
      ∑ S ∈ sets, (turboPotentialCountLaw22 R S).expect
        (fun counts => (A S counts) ^ α) := by
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
  rw [turboSetup22_count_disintegration R points S]
  exact bind_prob_count_le_fractional_mean_on_support
    (turboPotentialCountLaw22 R S) (turboSetupCountFiber22 R points S)
    (turboActiveBadCount22 R points d S) (A S) hα0 hα1 (hmean S hS)

#print axioms turboSetup22_bad_prob_le_count_fractional
#print axioms finiteStatisticFiber_disintegration
#print axioms turboSetup22_count_disintegration
end Spin.Turbo


