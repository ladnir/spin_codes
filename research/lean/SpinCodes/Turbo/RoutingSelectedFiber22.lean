import SpinCodes.Turbo.RoutingConditionedMean22

/-!
# Conditional latent selected subsets under actual potential-slot counts

The selected sets are latent and include zero-valued packets. The proof
reuses the checked profile mixture and the actual uniform permutation
pushforward, rather than recounting placements.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.Structured.Routing

private theorem map_mass_eq_fiber {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) (j : J) :
    (P.map f).p j = P.prob (fun ω => f ω = j) := by
  simp only [FinPMF.map_p, FinPMF.prob, sum_filter]

private theorem fiber_expect {Ω J : Type*} [Fintype Ω] [Fintype J]
    (P : FinPMF Ω) (f : Ω → J) (j : J)
    (hp : 0 < P.prob (fun ω => f ω = j)) (g : Ω → ℝ) :
    (finiteStatisticFiber P f j).expect g =
      P.expect (fun ω => if f ω = j then g ω else 0) /
        P.prob (fun ω => f ω = j) := by
  simp only [finiteStatisticFiber, dite_eq_left hp, FinPMF.expect,
    FinPMF.condition_p_apply]
  rw [sum_div]
  apply sum_congr rfl
  intro ω _
  split_ifs <;> ring

/-- Conditioning commutes with a pushforward when the statistic factors
through that pushforward. Only positive fibers are used. -/
theorem finiteStatisticFiber_map {Ω X J : Type*}
    [Fintype Ω] [Fintype X] [Fintype J]
    (P : FinPMF Ω) (g : Ω → X) (f : X → J) (j : J)
    (hp : 0 < P.prob (fun ω => f (g ω) = j)) :
    (finiteStatisticFiber P (fun ω => f (g ω)) j).map g =
      finiteStatisticFiber (P.map g) f j := by
  have hp' : 0 < (P.map g).prob (fun x => f x = j) := by
    rw [FinPMF.map_prob]
    exact hp
  apply FinPMF.eq_of_expect_eq
  intro k
  rw [FinPMF.expect_map, fiber_expect P _ j hp, fiber_expect (P.map g) f j hp',
    FinPMF.expect_map, FinPMF.map_prob]

/-- A supported mixture's positive fibers are exactly its component laws. -/
theorem finiteStatisticFiber_supported_mixture
    {Ω J : Type*} [Fintype Ω] [Fintype J]
    (L : FinPMF Ω) (P : FinPMF J) (Q : J → FinPMF Ω) (f : Ω → J)
    (hmix : L = P.bind Q) (hmarg : L.map f = P)
    (hsupport : ∀ j ω, f ω ≠ j → (Q j).p ω = 0)
    (j : J) (hp : 0 < P.p j) : finiteStatisticFiber L f j = Q j := by
  classical
  have hprob : L.prob (fun ω => f ω = j) = P.p j :=
    (map_mass_eq_fiber L f j).symm.trans (congrArg (fun K => K.p j) hmarg)
  have hpos : 0 < L.prob (fun ω => f ω = j) := hprob.symm ▸ hp
  rw [finiteStatisticFiber, dite_eq_left hpos]
  apply FinPMF.ext
  funext ω
  rw [FinPMF.condition_p_apply]
  by_cases he : f ω = j
  · rw [ite_eq_left he, hprob, hmix, FinPMF.bind_p]
    have hs : (∑ k, P.p k * (Q k).p ω) = P.p j * (Q j).p ω := by
      apply sum_eq_single j
      · intro k _ hk
        have hn : f ω ≠ k := by rw [he]; exact Ne.symm hk
        rw [hsupport k ω hn, mul_zero]
      · intro hn
        exact (hn (mem_univ j)).elim
    rw [hs]
    field_simp [hp.ne']
  · rw [ite_eq_right he, hsupport j ω he]

abbrev TurboSelectedStepSets22 (R : ℕ) := Fin 64 → Fin R → Finset (Fin 16)

def turboSelectedStepSets22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (perms : TurboRegionPerms22 R) : TurboSelectedStepSets22 R :=
  fun region => regionalStepSelectedSet22 (shuffleSupport (perms region) S)

def turboSelectedStepSetCounts22 {R : ℕ} (sets : TurboSelectedStepSets22 R) :
    TurboPotentialCounts22 R := fun region => regionalStepCounts22 (sets region)

def turboSelectedStepSetsLaw22 (R : ℕ) (S : Finset (Fin (R * 16))) :
    FinPMF (TurboSelectedStepSets22 R) :=
  Spin.piPMF (fun _ : Fin 64 => regionalSelectedStepSetsLaw22 R S)

def turboUniformSelectedStepSets22 {R : ℕ} (counts : TurboPotentialCounts22 R) :
    FinPMF (TurboSelectedStepSets22 R) :=
  Spin.piPMF (fun region : Fin 64 => uniformStepSelectedSets22 (counts region))

theorem turboSelectedStepSets22_law (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (turboRegionPermLaw22 R).map (turboSelectedStepSets22 S) =
      turboSelectedStepSetsLaw22 R S := by
  change (Spin.piPMF (fun _ : Fin 64 =>
    FinPMF.uniform (Equiv.Perm (Fin (R * 16))))).map
      (fun perms region => regionalStepSelectedSet22 (shuffleSupport (perms region) S)) = _
  rw [Spin.piPMF_map (fun _ : Fin 64 =>
    FinPMF.uniform (Equiv.Perm (Fin (R * 16))))
    (fun (_ : Fin 64) (σ : Equiv.Perm (Fin (R * 16))) =>
      regionalStepSelectedSet22 (R := R) (shuffleSupport σ S))]
  apply congrArg Spin.piPMF
  funext region
  rw [regionalSelectedStepSetsLaw22]
  rw [← uniform_permutation_shuffleLaw, FinPMF.map_comp]
  rfl

theorem turboSelectedStepSets22_count_law (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (turboSelectedStepSetsLaw22 R S).map turboSelectedStepSetCounts22 =
      turboPotentialCountLaw22 R S := by
  change (Spin.piPMF (fun _ : Fin 64 => regionalSelectedStepSetsLaw22 R S)).map
    (fun sets region => regionalStepCounts22 (sets region)) = _
  rw [Spin.piPMF_map]
  simp only [regionalSelectedStepSets22_count_law, turboPotentialCountLaw22]

theorem turboSelectedStepSets22_profile_mixture (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    turboSelectedStepSetsLaw22 R S =
      (turboPotentialCountLaw22 R S).bind turboUniformSelectedStepSets22 := by
  unfold turboSelectedStepSetsLaw22 turboPotentialCountLaw22 turboUniformSelectedStepSets22
  simp only [regionalSelectedStepSetsLaw22_eq_profile_mixture]
  exact (Spin.piPMF_bind _ _).symm

theorem turboUniformSelectedStepSets22_support {R : ℕ}
    (counts : TurboPotentialCounts22 R) (sets : TurboSelectedStepSets22 R)
    (hne : turboSelectedStepSetCounts22 sets ≠ counts) :
    (turboUniformSelectedStepSets22 counts).p sets = 0 := by
  have hn : ¬ ∀ region : Fin 64, regionalStepCounts22 (sets region) = counts region := by
    intro h
    exact hne (funext h)
  obtain ⟨region, hregion⟩ := not_forall.mp hn
  rw [turboUniformSelectedStepSets22, Spin.piPMF_apply]
  apply prod_eq_zero (mem_univ region)
  have hnot : ¬ ∀ step : Fin R,
      (sets region step).card = (counts region step).val := by
    intro h
    apply hregion
    funext step
    exact Fin.ext (h step)
  rw [uniformStepSelectedSets22_apply_profile, ite_eq_right hnot]

/-- The actual regional permutation fiber has independent uniform selected
subsets in every physical step, at the observed potential-slot counts. -/
theorem turboRegionCountFiber22_selected_sets (R : ℕ)
    (S : Finset (Fin (R * 16))) (counts : TurboPotentialCounts22 R)
    (hp : 0 < (turboPotentialCountLaw22 R S).p counts) :
    (turboRegionCountFiber22 R S counts).map (turboSelectedStepSets22 S) =
      turboUniformSelectedStepSets22 counts := by
  letI : DecidableEq (TurboPotentialCounts22 R) :=
    fun a b => Classical.propDecidable (a = b)
  have hprob : (turboRegionPermLaw22 R).prob
      (fun σ => turboPotentialCounts22 S σ = counts) =
      (turboPotentialCountLaw22 R S).p counts := by
    exact (map_mass_eq_fiber (turboRegionPermLaw22 R)
      (turboPotentialCounts22 S) counts).symm.trans
      (congrArg (fun K : FinPMF (TurboPotentialCounts22 R) => K.p counts)
        (turboPotentialCounts22_law R S))
  change (finiteStatisticFiber (turboRegionPermLaw22 R)
    (fun σ => turboSelectedStepSetCounts22 (turboSelectedStepSets22 S σ)) counts).map _ = _
  rw [finiteStatisticFiber_map _ _ _ _ (by exact hprob.symm ▸ hp),
    turboSelectedStepSets22_law]
  exact finiteStatisticFiber_supported_mixture _ _ _ _
    (turboSelectedStepSets22_profile_mixture R S)
    (turboSelectedStepSets22_count_law R S)
    turboUniformSelectedStepSets22_support counts hp

#print axioms turboRegionCountFiber22_selected_sets
end Spin.Turbo


