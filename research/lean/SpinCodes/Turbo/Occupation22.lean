import SpinCodes.Turbo.ActualSetup22

/-!
# Occupation of the actual Turbo 22-state encoder

The finite first moment is partitioned by the number of nonzero RS8
group messages. This is the Turbo instantiation of the old generic
`Framework.Setup.ZQ` partition; it refers to the same sampled outer,
route, and state-map setup as `ActualSetup22`.
-/

noncomputable section

namespace Spin.Turbo

open Finset

/-- Number of groups carrying a nonzero input message. -/
def turboOccupation22 {R : ℕ} (msg : TurboMessage22 R) : ℕ :=
  (Finset.univ.filter (fun i => msg i ≠ 0)).card

/-- The precise active group set, retained when the one- and two-group
folds account for placements and their common regional permutation. -/
def turboActivePositions22 {R : ℕ}
    (msg : TurboMessage22 R) : Finset (Fin (R * 16)) :=
  Finset.univ.filter (fun i => msg i ≠ 0)

def turboActiveMessages22 {R : ℕ}
    (S : Finset (Fin (R * 16))) : Finset (TurboMessage22 R) :=
  Finset.univ.filter (fun msg => turboActivePositions22 msg = S)

theorem turboOccupation22_eq_active_card {R : ℕ}
    (msg : TurboMessage22 R) :
    turboOccupation22 msg = (turboActivePositions22 msg).card := rfl

theorem turboOccupation22_pos {R : ℕ} (msg : TurboMessage22 R)
    (hne : msg ≠ 0) : 0 < turboOccupation22 msg := by
  apply Finset.card_pos.mpr
  obtain ⟨i, hi⟩ : ∃ i, msg i ≠ 0 := by
    by_contra h
    simp only [not_exists, not_not] at h
    exact hne (funext h)
  exact ⟨i, by simp [hi]⟩

theorem turboOccupation22_le {R : ℕ} (msg : TurboMessage22 R) :
    turboOccupation22 msg ≤ R * 16 := by
  simpa only [turboOccupation22, card_univ, Fintype.card_fin] using
    Finset.card_le_card
      (filter_subset (fun i : Fin (R * 16) => msg i ≠ 0) univ)

theorem turboOccupation22_mem_Ico {R : ℕ}
    (msg : TurboMessage22 R) (hne : msg ≠ 0) :
    turboOccupation22 msg ∈ Finset.Ico 1 (R * 16 + 1) := by
  simp only [Finset.mem_Ico]
  exact ⟨turboOccupation22_pos msg hne,
    Nat.lt_succ_of_le (turboOccupation22_le msg)⟩

/-- Every occupation sum splits exactly over active group sets. This is
the Turbo counterpart of `Structured/ConcreteOuterCountingProbability`,
with RS8 group inputs in place of its BA/Golay local messages. -/
theorem turboSumOccupation_eq_activeSets {R : ℕ} (q : ℕ)
    (f : TurboMessage22 R → ℝ) :
    (∑ msg ∈ (Finset.univ : Finset (TurboMessage22 R)).filter
        (fun msg => turboOccupation22 msg = q), f msg) =
      ∑ S ∈ (Finset.univ : Finset (Fin (R * 16))).powersetCard q,
        ∑ msg ∈ turboActiveMessages22 S, f msg := by
  classical
  calc
    _ = ∑ msg : TurboMessage22 R,
          if turboOccupation22 msg = q then f msg else 0 := by
        rw [Finset.sum_filter]
    _ = ∑ msg : TurboMessage22 R,
          ∑ S ∈ (Finset.univ : Finset (Fin (R * 16))).powersetCard q,
            if turboActivePositions22 msg = S then f msg else 0 := by
        apply Finset.sum_congr rfl
        intro msg _
        simp only [turboOccupation22_eq_active_card,
          Finset.sum_ite_eq,
          Finset.mem_powersetCard, Finset.subset_univ, true_and]
        rfl
    _ = _ := by
        rw [Finset.sum_comm]
        simp only [turboActiveMessages22, Finset.sum_filter]

/-- Exact partition of the realized bad-message count. -/
theorem turboSetup22_bad_count_by_occupation
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ)
    (ω : TurboOuterSeed22 R × RouteInnerSeed22 R) :
    (turboSetup22 R points).Z d ω =
      ∑ q ∈ Finset.Ico 1 (R * 16 + 1),
        (turboSetup22 R points).ZQ d
          (fun _ msg => turboOccupation22 msg) q ω := by
  apply Spin.Setup.Z_eq_sum_ZQ
  intro seed msg hmsg
  exact turboOccupation22_mem_Ico msg (Finset.mem_filter.mp hmsg).2

/-- The exact occupation decomposition is an identity under the one
actual sampled setup law. Each occupation term remains a shared-setup
first moment, not a product of probabilities for different messages. -/
theorem turboSetup22_expected_bad_count_by_occupation
    (R : ℕ) (points : Fin 8 → Field16) (d : ℕ) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).Z d ω : ℝ)) =
      ∑ q ∈ Finset.Ico 1 (R * 16 + 1),
        (turboSetup22 R points).joint.expect
          (fun ω => ((turboSetup22 R points).ZQ d
            (fun _ msg => turboOccupation22 msg) q ω : ℝ)) := by
  apply Spin.Setup.expect_Z_eq_sum_expect_ZQ
  intro seed msg hmsg
  exact turboOccupation22_mem_Ico msg (Finset.mem_filter.mp hmsg).2

/-- Each occupation contribution reduces to an average over the actual
outer setup of the conditional route-and-inner failure probability. -/
theorem turboSetup22_occupation_first_moment
    (R : ℕ) (points : Fin 8 → Field16) (d q : ℕ) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) q ω : ℝ)) =
      (turboOuterLaw22 R).expect
        (fun seed => ∑ msg ∈
          (Spin.nonzeroMsgs (TurboMessage22 R)).filter
            (fun msg => turboOccupation22 msg = q),
          (turboSetup22 R points).qd d
            (turboOuterWord22 points seed msg)) :=
  Spin.Setup.expect_ZQ_eq _ _ _ _

/-- The occupation first moment keeps each selected active set and the
same sampled encoder. This is the starting point for the exact `U₁` and
`U₂` packet folds and for the larger-occupation transfer. -/
theorem turboSetup22_occupation_first_moment_activeSets
    (R : ℕ) (points : Fin 8 → Field16) (d q : ℕ) (hq : 0 < q) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) q ω : ℝ)) =
      (turboOuterLaw22 R).expect
        (fun seed =>
          ∑ S ∈ (Finset.univ : Finset (Fin (R * 16))).powersetCard q,
            ∑ msg ∈ turboActiveMessages22 S,
              (turboSetup22 R points).qd d
                (turboOuterWord22 points seed msg)) := by
  rw [turboSetup22_occupation_first_moment]
  apply congrArg
  funext seed
  have hnonzero (msg : TurboMessage22 R)
      (hqmsg : turboOccupation22 msg = q) : msg ≠ 0 := by
    intro hm
    subst msg
    simp [turboOccupation22] at hqmsg
    omega
  have hfilter :
      (Spin.nonzeroMsgs (TurboMessage22 R)).filter
        (fun msg => turboOccupation22 msg = q) =
      (Finset.univ : Finset (TurboMessage22 R)).filter
        (fun msg => turboOccupation22 msg = q) := by
    ext msg
    simp only [Finset.mem_filter, Spin.nonzeroMsgs,
      Finset.mem_univ, true_and]
    exact ⟨fun h => h.2, fun h => ⟨hnonzero msg h, h⟩⟩
  rw [hfilter]
  exact turboSumOccupation_eq_activeSets q _

#print axioms turboSetup22_bad_count_by_occupation
#print axioms turboSetup22_expected_bad_count_by_occupation
#print axioms turboSetup22_occupation_first_moment
#print axioms turboSumOccupation_eq_activeSets
#print axioms turboSetup22_occupation_first_moment_activeSets

end Spin.Turbo
