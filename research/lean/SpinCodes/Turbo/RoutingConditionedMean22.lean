import SpinCodes.Turbo.RoutingCountFibers22
import SpinCodes.Turbo.RS8ReferenceLaw22

/-!
# Independence and outer domination inside actual routing-count fibers

The count statistic uses only regional permutations. Conditioning retains
the original independent shared outer, row shuffles, and GL update law.
This is the actual-law counterpart of the old condition_prod_left argument.
-/
noncomputable section
set_option maxHeartbeats 1000000
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset

private theorem map_mass_eq_fiber {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) (j : J) :
    (P.map f).p j = P.prob (fun ω => f ω = j) := by
  simp only [FinPMF.map_p, FinPMF.prob, sum_filter]

def turboRegionCountFiber22 (R : ℕ) (S : Finset (Fin (R * 16)))
    (counts : TurboPotentialCounts22 R) : FinPMF (TurboRegionPerms22 R) :=
  finiteStatisticFiber (turboRegionPermLaw22 R) (turboPotentialCounts22 S) counts

/-- Count conditioning changes only the regional-permutation law. -/
theorem turboSetupCountFiber22_eq_product (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (counts : TurboPotentialCounts22 R) :
    turboSetupCountFiber22 R points S counts =
      (turboOuterLaw22 R).prod
        (((turboRowPermLaw22 R).prod (turboRegionCountFiber22 R S counts)).prod
          (turboStateLaw22 R)) := by
  letI : DecidableEq (TurboPotentialCounts22 R) :=
    fun a b => Classical.propDecidable (a = b)
  have hactual : (turboSetup22 R points).joint.prob
      (fun ω => turboSetupPotentialCounts22 S ω = counts) =
      (turboPotentialCountLaw22 R S).p counts := by
    exact (map_mass_eq_fiber (turboSetup22 R points).joint
      (turboSetupPotentialCounts22 S) counts).symm.trans
      (congrArg (fun P : FinPMF (TurboPotentialCounts22 R) => P.p counts)
        (turboSetupPotentialCounts22_law R points S))
  have hregion : (turboRegionPermLaw22 R).prob
      (fun σ => turboPotentialCounts22 S σ = counts) =
      (turboPotentialCountLaw22 R S).p counts := by
    exact (map_mass_eq_fiber (turboRegionPermLaw22 R)
      (turboPotentialCounts22 S) counts).symm.trans
      (congrArg (fun P : FinPMF (TurboPotentialCounts22 R) => P.p counts)
        (turboPotentialCounts22_law R S))
  unfold turboSetupCountFiber22 turboRegionCountFiber22 finiteStatisticFiber
  by_cases hp : 0 < (turboPotentialCountLaw22 R S).p counts
  · rw [dite_eq_left (by exact hactual.symm ▸ hp), dite_eq_left (by exact hregion.symm ▸ hp)]
    apply FinPMF.ext
    funext ω
    rw [FinPMF.condition_p_apply]
    change (if turboSetupPotentialCounts22 S ω = counts then
      (turboOuterLaw22 R).p ω.1 *
        (((turboRowPermLaw22 R).p ω.2.1.1 * (turboRegionPermLaw22 R).p ω.2.1.2) *
          (turboStateLaw22 R).p ω.2.2) /
          (turboSetup22 R points).joint.prob
            (fun ω => turboSetupPotentialCounts22 S ω = counts) else 0) =
      (turboOuterLaw22 R).p ω.1 *
        (((turboRowPermLaw22 R).p ω.2.1.1 *
          (if turboPotentialCounts22 S ω.2.1.2 = counts then
            (turboRegionPermLaw22 R).p ω.2.1.2 /
              (turboRegionPermLaw22 R).prob
                (fun σ => turboPotentialCounts22 S σ = counts) else 0)) *
          (turboStateLaw22 R).p ω.2.2)
    erw [hactual, hregion]
    change (if turboPotentialCounts22 S ω.2.1.2 = counts then _ else _) = _
    split_ifs <;> ring
  · rw [dite_eq_right (by exact hactual.symm ▸ hp), dite_eq_right (by exact hregion.symm ▸ hp)]
    rfl

/-- The arbitrary fiber observable can be a continuous inner lower-tail
probability or a tilted chronological full-state moment. Outer domination
is applied before averaging that observable over the conditional route. -/
theorem turboShuffledCountMoment22_le_reference_observable
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (f : TurboOuterWord22 R → ℝ)
    (hf : ∀ words, 0 ≤ f words) :
    (∑ words : TurboOuterWord22 R,
      turboShuffledTupleCount22 R points S words * f words) ≤
      turboOldUniformFactor22 ^ S.card *
        ∑ words : TurboOuterWord22 R,
          (∏ i : Fin (R * 16),
            turboReferenceRowDensity22 points (i ∈ S) (words i)) * f words := by
  simp only [turboShuffledTupleCount22_eq_prod]
  exact finite_weighted_sum_le_factor _ _ f _
    (fun words => turboPositionProduct22_le_reference R points hpoints S words) hf

/-- Regional routing conditional on potential counts, with unchanged GL maps. -/
def turboFiberRemainingFailure22 (R : ℕ)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R)
    (words : TurboOuterWord22 R) (d : ℕ) : ℝ :=
  ((turboRegionCountFiber22 R S counts).prod (turboStateLaw22 R)).prob
    (fun seed => Spin.UniformMixer.sampledOutputWeight
      (flatInput (turboRouteRegions22 words seed.1)) seed.2 0 ≤ d)

private theorem expect_card_filter_decidable
    {Ω M : Type*} [Fintype Ω] [Fintype M]
    (P : FinPMF Ω) (s : Finset M) (E : M → Ω → Prop)
    [∀ x ω, Decidable (E x ω)] :
    P.expect (fun ω => ((s.filter (fun x => E x ω)).card : ℝ)) =
      ∑ x ∈ s, P.prob (E x) := by
  have he (ω : Ω) : ((s.filter (fun x => E x ω)).card : ℝ) =
      ∑ x ∈ s, if E x ω then (1 : ℝ) else 0 := by
    simp only [Finset.card_filter, Nat.cast_sum]
    exact sum_congr rfl (fun x _ => by split_ifs <;> simp)
  simp_rw [he]
  rw [FinPMF.expect_sum]
  exact sum_congr rfl (fun x _ => (FinPMF.prob_eq_expect_indicator _ _).symm)
/-- Exact active-set count under the actual conditional setup. Row shuffles
remain outside the shared outer message sum, as in the old ordinary proof. -/
theorem turboActiveBadCount22_fiber_eq_shuffled_count
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R) (d : ℕ) :
    (turboSetupCountFiber22 R points S counts).expect
      (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) =
    ∑ words : TurboOuterWord22 R,
      turboShuffledTupleCount22 R points S words *
        turboFiberRemainingFailure22 R S counts words d := by
  rw [turboSetupCountFiber22_eq_product]
  unfold turboActiveBadCount22
  rw [expect_card_filter_decidable]
  calc
    _ = ∑ msg ∈ turboActiveMessages22 S,
          (turboRowShuffledLaw22 R points msg).expect
            (fun words => turboFiberRemainingFailure22 R S counts words d) := by
      apply sum_congr rfl
      intro msg _
      simp only [turboRowShuffledLaw22, FinPMF.expect_map,
        FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod,
        turboFiberRemainingFailure22]
      rfl
    _ = _ := by
      simp only [FinPMF.expect, turboShuffledTupleCount22]
      rw [sum_comm]
      simp only [← sum_mul]

/-- Pointwise outer domination remains valid inside each actual routing
fiber. This is a conditional mean bound, not the unconditional first moment. -/
theorem turboActiveBadCount22_fiber_le_reference
    (R : ℕ) (points : Fin 8 → Field16)
    (hpoints : Function.Injective points)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R) (d : ℕ) :
    (turboSetupCountFiber22 R points S counts).expect
      (fun ω => (turboActiveBadCount22 R points d S ω : ℝ)) ≤
    turboOldUniformFactor22 ^ S.card *
      (turboReferenceArrayLaw22 R points S).expect
        (fun words => turboFiberRemainingFailure22 R S counts words d) := by
  rw [turboActiveBadCount22_fiber_eq_shuffled_count]
  have h := turboShuffledCountMoment22_le_reference_observable R points hpoints S
    (fun words => turboFiberRemainingFailure22 R S counts words d)
    (fun words => FinPMF.prob_nonneg _ _)
  simpa only [FinPMF.expect, turboReferenceArrayLaw22_apply] using h

#print axioms turboActiveBadCount22_fiber_eq_shuffled_count
#print axioms turboActiveBadCount22_fiber_le_reference
#print axioms turboSetupCountFiber22_eq_product
#print axioms turboShuffledCountMoment22_le_reference_observable
end Spin.Turbo






