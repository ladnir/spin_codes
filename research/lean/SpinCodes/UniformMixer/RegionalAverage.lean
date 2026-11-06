import SpinCodes.UniformMixer.RegionalTrace

/-!
# Averaging a uniformly shuffled regional kernel

The chronological kernel has nonnegative entries and row sums at most one
for a tilt in `[0,1]`. This lets us charge discarded placements by their
probability before comparing the good placements with the two-state kernel.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Routing
attribute [local instance] Classical.propDecidable

def Substochastic (K : State22 → State22 → ℝ) : Prop :=
  (∀ q r, 0 ≤ K q r) ∧ (∀ q, (∑ r, K q r) ≤ 1)

theorem identity_substochastic :
    Substochastic (identity (S := State22)) := by
  constructor
  · intro q r
    by_cases h : q = r <;> simp [identity, h]
  · intro q
    simp [identity]

theorem stepKernel_substochastic (z : ℝ)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) (x : Fin 64 → ZMod 2) :
    Substochastic (stepKernel z x) := by
  constructor
  · intro q r
    exact mul_nonneg (pow_nonneg hz0 _) (update_nonneg _ _ _)
  · intro q
    calc
      (∑ r, stepKernel z x q r) =
          z ^ stepWeight x q * ∑ r : State22,
            update (stepFeedback x) q r := by
          simp [stepKernel, Finset.mul_sum]
      _ = z ^ stepWeight x q := by rw [update_sum]; ring
      _ ≤ 1 := pow_le_one₀ hz0 hz1

theorem compose_substochastic (K L : State22 → State22 → ℝ)
    (hK : Substochastic K) (hL : Substochastic L) :
    Substochastic (compose K L) := by
  constructor
  · intro q r
    exact compose_nonneg K L hK.1 hL.1 q r
  · intro q
    calc
      (∑ r, compose K L q r) =
          ∑ v : State22, K q v * ∑ r : State22, L v r := by
        simp only [compose, Finset.mul_sum]
        rw [Finset.sum_comm]
      _ ≤ ∑ v : State22, K q v := by
        apply Finset.sum_le_sum
        intro v _
        have h := mul_le_mul_of_nonneg_left (hL.2 v) (hK.1 q v)
        simpa using h
      _ ≤ 1 := hK.2 q

theorem inputListKernel_substochastic (z : ℝ)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (xs : List (Fin 64 → ZMod 2)) :
    Substochastic (inputListKernel z xs) := by
  induction xs with
  | nil => exact identity_substochastic
  | cons x tail ih =>
      exact compose_substochastic _ _
        (stepKernel_substochastic z hz0 hz1 x) ih

theorem substochastic_entry_le_one (K : State22 → State22 → ℝ)
    (hK : Substochastic K) (q r : State22) : K q r ≤ 1 := by
  have hsingle : K q r ≤ ∑ v : State22, K q v := by
    exact Finset.single_le_sum
      (fun v _ => hK.1 q v) (Finset.mem_univ r)
  exact hsingle.trans (hK.2 q)

theorem chronologicalRegionKernel_unit_interval {R : ℕ}
    (z : ℝ) (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (T : Finset (Fin (64 * R))) (q r : State22) :
    0 ≤ chronologicalRegionKernel z T q r ∧
      chronologicalRegionKernel z T q r ≤ 1 := by
  have h := inputListKernel_substochastic z hz0 hz1
    (chronologicalInputList T)
  exact ⟨h.1 q r, substochastic_entry_le_one _ h q r⟩

/-- The actual regional kernel averaged over a uniformly shuffled support
of the given cardinality. -/
def averagedRegionKernel {R : ℕ} (S : Finset (Fin (64 * R)))
    (z : ℝ) (q r : State22) : ℝ :=
  (shuffleLaw S).expect (fun T => chronologicalRegionKernel z T q r)

theorem averagedRegionKernel_bad_error {R : ℕ} (hR : 0 < R)
    (S : Finset (Fin (64 * R)))
    (z : ℝ) (hz0 : 0 ≤ z) (hz1 : z ≤ 1) (q r : State22) :
    averagedRegionKernel S z q r ≤
      (shuffleLaw S).expect (fun T =>
        if ¬ badRegionalSites T then chronologicalRegionKernel z T q r
        else 0) +
      128 * (S.card : ℝ) / (64 * R) +
      257 * (S.card : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1) := by
  let P := shuffleLaw S
  let f : Finset (Fin (64 * R)) → ℝ :=
    fun T => chronologicalRegionKernel z T q r
  have hpoint (T : Finset (Fin (64 * R))) :
      f T ≤ (if ¬ badRegionalSites T then f T else 0) +
        (if badRegionalSites T then 1 else 0) := by
    by_cases hb : badRegionalSites T
    · simpa [hb] using
        (chronologicalRegionKernel_unit_interval z hz0 hz1 T q r).2
    · simp [hb]
  have hbase : P.expect f ≤
      P.expect (fun T => if ¬ badRegionalSites T then f T else 0) +
        P.prob badRegionalSites := by
    calc
      P.expect f ≤ P.expect (fun T =>
          (if ¬ badRegionalSites T then f T else 0) +
            (if badRegionalSites T then 1 else 0)) :=
        P.expect_mono hpoint
      _ = _ := by rw [Spin.FinPMF.expect_add,
        Spin.FinPMF.prob_eq_expect_indicator]
  have hbad := shuffle_prob_badRegionalSites
    (by omega : 2 ≤ 64 * R) S
  unfold averagedRegionKernel
  dsimp [P, f] at hbase
  have hfinal := hbase.trans (add_le_add_right hbad
    ((shuffleLaw S).expect (fun T =>
      if ¬ badRegionalSites T then chronologicalRegionKernel z T q r else 0)))
  simpa only [Nat.cast_mul, Nat.cast_ofNat, add_assoc] using hfinal

/-- The two-state product envelope, with its full-state boundary maps. -/
def coarseRegionalEnvelope {R : ℕ} (θ : ℝ)
    (T : Finset (Fin (64 * R))) (q r : State22) : ℝ :=
  let L : ℝ := 64 * R
  let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
  (1 + C) ^ (T.card + 1) *
    (liftM * coarseRegionProduct
      (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
      (supportMarks R (T.sort (· ≤ ·)) 0)
      (supportLastGap R (T.sort (· ≤ ·)) 0) * projectM) q r

theorem averagedRegionKernel_coarse_upper {R : ℕ} (hR : 0 < R)
    (S : Finset (Fin (64 * R))) (θ : ℝ) (hθ : 0 ≤ θ)
    (q r : State22) :
    let L : ℝ := 64 * R
    averagedRegionKernel S (Real.exp (-(θ / L))) q r ≤
      (shuffleLaw S).expect (fun T =>
        if ¬ badRegionalSites T then coarseRegionalEnvelope θ T q r
        else 0) +
      128 * (S.card : ℝ) / L +
      257 * (S.card : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1) := by
  have hL : (0 : ℝ) < 64 * R := by
    exact_mod_cast Nat.mul_pos (by omega : 0 < 64) hR
  have hz0 : 0 ≤ Real.exp (-(θ / (64 * R : ℝ))) :=
    (Real.exp_pos _).le
  have hz1 : Real.exp (-(θ / (64 * R : ℝ))) ≤ 1 :=
    Real.exp_le_one_iff.mpr (neg_nonpos.mpr (div_nonneg hθ hL.le))
  have hbad := averagedRegionKernel_bad_error hR S
    (Real.exp (-(θ / (64 * R : ℝ)))) hz0 hz1 q r
  have hgood :
      (shuffleLaw S).expect (fun T =>
        if ¬ badRegionalSites T then
          chronologicalRegionKernel (Real.exp (-(θ / (64 * R : ℝ)))) T q r
        else 0) ≤
      (shuffleLaw S).expect (fun T =>
        if ¬ badRegionalSites T then coarseRegionalEnvelope θ T q r
        else 0) := by
    apply Spin.FinPMF.expect_mono
    intro T
    by_cases hg : ¬ badRegionalSites T
    · simpa [hg, coarseRegionalEnvelope] using
        good_chronologicalRegion_coarse_upper hR T hg θ hθ q r
    · simp [hg]
  have hbad' : averagedRegionKernel S (Real.exp (-(θ / (64 * R : ℝ)))) q r ≤
      (shuffleLaw S).expect (fun T =>
        if ¬ badRegionalSites T then
          chronologicalRegionKernel (Real.exp (-(θ / (64 * R : ℝ)))) T q r
        else 0) +
      (128 * (S.card : ℝ) / (64 * R : ℝ) +
        257 * (S.card : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1)) := by
    simpa only [add_assoc] using hbad
  have hfinal := hbad'.trans (add_le_add_left hgood
    (128 * (S.card : ℝ) / (64 * R : ℝ) +
      257 * (S.card : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1)))
  simpa only [add_assoc] using hfinal

end Spin.UniformMixer
