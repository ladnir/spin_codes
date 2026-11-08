import SpinCodes.Turbo.RowShuffledArray22

/-!
# Regional routing and the continuous Turbo inner after row shuffling

The first route stage shuffles complete packet rows independently.
The remaining regional permutations and the 22-state inner are kept as
one joint experiment. In particular, the state is never reset between
regions, and no product of per-region failure probabilities is used.

The old counterpart is `Structured/ConcreteFixedRegionComposition`,
whose support-valued region kernels do not model Turbo packet labels.
-/

noncomputable section

namespace Spin.Turbo

open Spin.Structured.ConcreteRoute Spin.UniformMixer

abbrev TurboRegionPerms22 (R : ℕ) :=
  Fin 64 → Equiv.Perm (Fin (R * 16))

def turboRegionPermLaw22 (R : ℕ) : FinPMF (TurboRegionPerms22 R) :=
  piPMF (fun _ : Fin 64 =>
    FinPMF.uniform (Equiv.Perm (Fin (R * 16))))

def turboStateLaw22 (R : ℕ) :
    FinPMF (Fin (64 * R) → StateMap) :=
  piPMF (fun _ : Fin (64 * R) => FinPMF.uniform StateMap)

def turboRemainingLaw22 (R : ℕ) :
    FinPMF (TurboRegionPerms22 R ×
      (Fin (64 * R) → StateMap)) :=
  (turboRegionPermLaw22 R).prod (turboStateLaw22 R)

def turboShuffleRows22 {R : ℕ}
    (rows : TurboOuterWord22 R) (perm : TurboRowPerms22 R) :
    TurboOuterWord22 R :=
  fun i => shufflePacketWord (perm i) (rows i)

def turboRouteRegions22 {R : ℕ}
    (rows : TurboOuterWord22 R)
    (perm : TurboRegionPerms22 R) :
    Fin 64 → Fin (R * 16) → Packet4 :=
  fun j i => rows (perm j i) j

theorem turboRoute_two_stages22 {R : ℕ}
    (rows : TurboOuterWord22 R)
    (rowPerm : TurboRowPerms22 R)
    (regionPerm : TurboRegionPerms22 R) :
    routePackets rows (rowPerm, regionPerm) =
      turboRouteRegions22 (turboShuffleRows22 rows rowPerm)
        regionPerm := by
  rfl

/-- The conditional failure probability after row shuffling. The
regional route and all sampled state updates are one joint experiment. -/
def turboRemainingFailure22 {R : ℕ}
    (rows : TurboOuterWord22 R) (d : ℕ) : ℝ :=
  (turboRemainingLaw22 R).prob (fun seed =>
    sampledOutputWeight
      (flatInput (turboRouteRegions22 rows seed.1))
      seed.2 0 ≤ d)

/-- Fubini for the *actual* route-and-inner setup: first average the
independent row permutations, then test the one continuous regional
encoder. -/
theorem turboSetup22_qd_row_shuffled
    (R : ℕ) (points : Fin 8 → Field16)
    (rows : TurboOuterWord22 R) (d : ℕ) :
    (turboSetup22 R points).qd d rows =
      (turboRowPermLaw22 R).expect (fun rowPerm =>
        turboRemainingFailure22
          (turboShuffleRows22 rows rowPerm) d) := by
  change (((turboRowPermLaw22 R).prod
      (turboRegionPermLaw22 R)).prod
      (turboStateLaw22 R)).prob
      (fun seed => sampledOutputWeight
        (flatInput (routePackets rows seed.1)) seed.2 0 ≤ d) = _
  rw [FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod,
    FinPMF.expect_prod]
  apply congrArg
  funext rowPerm
  rw [turboRemainingFailure22, turboRemainingLaw22,
    FinPMF.prob_eq_expect_indicator, FinPMF.expect_prod]
  apply congrArg
  funext regionPerm
  apply congrArg
  funext maps
  rw [turboRoute_two_stages22]

/-- The exact active-set first moment is a nonnegative test of the
complete shuffled packet-array count. The common regional permutations
and continuous inner state remain inside `turboRemainingFailure22`.
This is the semantic bridge required before numerical `W_v` or
`V_{u,v}` endpoints can bound the one- and two-group terms. -/
theorem turboSetup22_activeSet_expect_eq_shuffled_count
    (R : ℕ) (points : Fin 8 → Field16)
    (S : Finset (Fin (R * 16))) (d : ℕ) :
    (turboOuterLaw22 R).expect (fun outerSeed =>
      ∑ msg ∈ turboActiveMessages22 S,
        (turboSetup22 R points).qd d
          (turboOuterWord22 points outerSeed msg)) =
    ∑ words : TurboOuterWord22 R,
      turboShuffledTupleCount22 R points S words *
        turboRemainingFailure22 words d := by
  calc
    _ = ∑ msg ∈ turboActiveMessages22 S,
          (turboOuterLaw22 R).expect (fun outerSeed =>
            (turboSetup22 R points).qd d
              (turboOuterWord22 points outerSeed msg)) := by
        rw [FinPMF.expect_sum]
    _ = ∑ msg ∈ turboActiveMessages22 S,
          (turboRowShuffledLaw22 R points msg).expect
            (fun words => turboRemainingFailure22 words d) := by
        apply Finset.sum_congr rfl
        intro msg _
        simp only [turboSetup22_qd_row_shuffled]
        rw [turboRowShuffledLaw22, FinPMF.expect_map,
          FinPMF.expect_prod]
        rfl
    _ = _ := by
        simp only [FinPMF.expect, turboShuffledTupleCount22]
        rw [Finset.sum_comm]
        simp only [← Finset.sum_mul]

/-- Exact finite occupation formula for the actual 22-state Turbo setup.
The local factors are complete packet-word counts after row shuffling;
the final factor still tests the entire regional route and one
continuous state evolution. This formula is valid before any finite
numerical envelope is supplied. -/
theorem turboSetup22_occupation_complete_formula
    (R : ℕ) (points : Fin 8 → Field16)
    (d q : ℕ) (hq : 0 < q) :
    (turboSetup22 R points).joint.expect
      (fun ω => ((turboSetup22 R points).ZQ d
        (fun _ msg => turboOccupation22 msg) q ω : ℝ)) =
    ∑ S ∈ (Finset.univ : Finset (Fin (R * 16))).powersetCard q,
      ∑ words : TurboOuterWord22 R,
        (∏ i : Fin (R * 16),
          turboShuffledPositionCount22 points (i ∈ S) (words i)) *
          turboRemainingFailure22 words d := by
  rw [turboSetup22_occupation_first_moment_activeSets R points d q hq]
  rw [FinPMF.expect_sum]
  apply Finset.sum_congr rfl
  intro S _
  rw [turboSetup22_activeSet_expect_eq_shuffled_count]
  simp only [turboShuffledTupleCount22_eq_prod]

#print axioms turboRoute_two_stages22
#print axioms turboSetup22_qd_row_shuffled
#print axioms turboSetup22_activeSet_expect_eq_shuffled_count
#print axioms turboSetup22_occupation_complete_formula

end Spin.Turbo
