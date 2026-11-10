import SpinCodes.Turbo.RS8ReferenceStepArray22
import SpinCodes.Turbo.ConditionalFractionalMoment

/-!
# Actual S-specific potential-slot count law

The old uniform-permutation support theorem and the checked RS8 profile
mixture determine these marginals. Selected slots include zero-valued
packets; the statistic depends only on regional routing and the fixed S.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.Structured.Routing

/-- Conditional step subsets have their displayed potential-slot counts. -/
theorem uniformStepSelectedSets22_count_expect {R : ℕ}
    (counts : Fin R → Fin 17) (f : (Fin R → Fin 17) → ℝ) :
    (uniformStepSelectedSets22 counts).expect
      (fun sets => f (regionalStepCounts22 sets)) = f counts := by
  unfold FinPMF.expect
  calc
    _ = ∑ sets, (uniformStepSelectedSets22 counts).p sets * f counts := by
      apply sum_congr rfl
      intro sets _
      change (uniformStepSelectedSets22 counts).p sets * f (regionalStepCounts22 sets) =
        (uniformStepSelectedSets22 counts).p sets * f counts
      by_cases he : regionalStepCounts22 sets = counts
      · rw [he]
      · have hn : ¬ ∀ step : Fin R, (sets step).card = (counts step).val := by
          intro h
          apply he
          funext step
          exact Fin.ext (h step)
        rw [uniformStepSelectedSets22_apply_profile, ite_eq_right hn]
        simp
    _ = f counts := by rw [← sum_mul, FinPMF.total, one_mul]

/-- The pre-existing profile mixture is the actual pushforward count law. -/
theorem regionalSelectedStepSets22_count_law (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (regionalSelectedStepSetsLaw22 R S).map regionalStepCounts22 =
      regionalStepProfileLaw22 R S := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [FinPMF.expect_map, regionalSelectedStepSetsLaw22_eq_profile_mixture,
    FinPMF.expect_bind]
  simp only [uniformStepSelectedSets22_count_expect]

/-- Counts of potential S-slots after one actual regional permutation. -/
def turboRegionPotentialCounts22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (σ : Equiv.Perm (Fin (R * 16))) : Fin R → Fin 17 :=
  regionalStepCounts22 (regionalStepSelectedSet22 (shuffleSupport σ S))

theorem turboRegionPotentialCounts22_law (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (FinPMF.uniform (Equiv.Perm (Fin (R * 16)))).map
      (turboRegionPotentialCounts22 S) = regionalStepProfileLaw22 R S := by
  calc
    _ = ((FinPMF.uniform (Equiv.Perm (Fin (R * 16)))).map
          (fun σ => shuffleSupport σ S)).map
          (fun T => regionalStepCounts22 (regionalStepSelectedSet22 T)) := by
        rw [FinPMF.map_comp]
        rfl
    _ = (regionalSelectedStepSetsLaw22 R S).map regionalStepCounts22 := by
        rw [uniform_permutation_shuffleLaw, regionalSelectedStepSetsLaw22,
          FinPMF.map_comp]
        rfl
    _ = _ := regionalSelectedStepSets22_count_law R S

abbrev TurboPotentialCounts22 (R : ℕ) := Fin 64 → Fin R → Fin 17

def turboPotentialCounts22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (perms : TurboRegionPerms22 R) : TurboPotentialCounts22 R :=
  fun region => turboRegionPotentialCounts22 S (perms region)

def turboPotentialCountLaw22 (R : ℕ) (S : Finset (Fin (R * 16))) :
    FinPMF (TurboPotentialCounts22 R) :=
  Spin.piPMF (fun _ : Fin 64 => regionalStepProfileLaw22 R S)

theorem turboPotentialCounts22_law (R : ℕ)
    (S : Finset (Fin (R * 16))) :
    (turboRegionPermLaw22 R).map (turboPotentialCounts22 S) =
      turboPotentialCountLaw22 R S := by
  change (Spin.piPMF (fun _ : Fin 64 =>
    FinPMF.uniform (Equiv.Perm (Fin (R * 16))))).map
      (fun perms region => turboRegionPotentialCounts22 S (perms region)) = _
  rw [Spin.piPMF_map]
  simp only [turboRegionPotentialCounts22_law, turboPotentialCountLaw22]

/-- Extract only the regional routing coordinate from the actual shared setup. -/
def turboSetupPotentialCounts22 {R : ℕ} (S : Finset (Fin (R * 16)))
    (ω : TurboOuterSeed22 R × RouteInnerSeed22 R) : TurboPotentialCounts22 R :=
  turboPotentialCounts22 S ω.2.1.2

theorem turboSetupPotentialCounts22_law (R : ℕ)
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16))) :
    ((turboSetup22 R points).joint).map (turboSetupPotentialCounts22 S) =
      turboPotentialCountLaw22 R S := by
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [FinPMF.expect_map]
  change ((turboOuterLaw22 R).prod
    ((Spin.Structured.ConcreteRoute.seedLaw (R * 16) 64).prod
      (turboStateLaw22 R))).expect
        (fun ω => f (turboPotentialCounts22 S ω.2.1.2)) = _
  simp only [Spin.Structured.ConcreteRoute.seedLaw, FinPMF.expect_prod,
    FinPMF.expect_const]
  rw [← FinPMF.expect_map]
  exact congrArg (fun P => P.expect f) (turboPotentialCounts22_law R S)

#print axioms turboRegionPotentialCounts22_law
#print axioms turboSetupPotentialCounts22_law
end Spin.Turbo

