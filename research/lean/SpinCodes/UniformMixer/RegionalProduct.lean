import SpinCodes.UniformMixer.Regional22

/-!
# A separated sequence of active regional steps

Each `(g,j)` records `g` empty steps followed by one active bit at local
coordinate `j`. A final empty gap follows the last active bit. The upper
comparison retains the full entering and exiting states, so it can be
composed across regions without resampling a boundary state.
-/

noncomputable section

namespace Spin.UniformMixer

theorem production_emptyGap_nonneg (z : ℝ) (hz : 0 ≤ z)
    (g : ℕ) (q r : State22) :
    0 ≤ emptyGap z (coordinateWeight coordMap) g q r := by
  induction g generalizing q r with
  | zero =>
      simp only [emptyGap]
      split_ifs <;> norm_num
  | succ g ih =>
      rw [emptyGap]
      apply Finset.sum_nonneg
      intro v _
      have he : 0 ≤ emptyStep z (coordinateWeight coordMap) v r :=
        mul_nonneg (pow_nonneg hz _) (update_nonneg 0 v r)
      exact mul_nonneg (ih q v) he

theorem productionImpulse_nonneg (z : ℝ) (hz : 0 ≤ z)
    (j : Fin 64) (q r : State22) :
    0 ≤ productionImpulse z j q r :=
  mul_nonneg (pow_nonneg hz _) (update_nonneg (feedbackColumn j) q r)

/-- Exact weighted kernel of a region with separated active steps. -/
def markedRegionKernel (z : ℝ) :
    List (ℕ × Fin 64) → ℕ → State22 → State22 → ℝ
  | [], g => emptyGap z (coordinateWeight coordMap) g
  | (g, j) :: tail, last =>
      compose (emptyGap z (coordinateWeight coordMap) g)
        (compose (productionImpulse z j) (markedRegionKernel z tail last))

/-- Corresponding full-state envelope after replacing every nonempty gap
by a zero/live lift and every active step by its unweighted transition. -/
def markedRegionEnvelope (α : ℕ → ℝ) :
    List (ℕ × Fin 64) → ℕ → State22 → State22 → ℝ
  | [], g => liftedGap (α g)
  | (g, j) :: tail, last =>
      compose (liftedGap (α g))
        (compose (update (feedbackColumn j))
          (markedRegionEnvelope α tail last))

theorem markedRegionKernel_nonneg (z : ℝ) (hz : 0 ≤ z)
    (marks : List (ℕ × Fin 64)) (last : ℕ) (q r : State22) :
    0 ≤ markedRegionKernel z marks last q r := by
  induction marks generalizing q r with
  | nil => exact production_emptyGap_nonneg z hz last q r
  | cons p tail ih =>
      rcases p with ⟨g, j⟩
      change 0 ≤ compose (emptyGap z (coordinateWeight coordMap) g)
        (compose (productionImpulse z j) (markedRegionKernel z tail last)) q r
      apply compose_nonneg
      · exact production_emptyGap_nonneg z hz g
      · intro u v
        apply compose_nonneg
        · exact productionImpulse_nonneg z hz j
        · exact ih

theorem markedRegionEnvelope_nonneg (α : ℕ → ℝ)
    (hα : ∀ g, 0 ≤ α g)
    (marks : List (ℕ × Fin 64)) (last : ℕ) (q r : State22) :
    0 ≤ markedRegionEnvelope α marks last q r := by
  induction marks generalizing q r with
  | nil => exact liftedGap_nonneg (α last) (hα last) q r
  | cons p tail ih =>
      rcases p with ⟨g, j⟩
      change 0 ≤ compose (liftedGap (α g))
        (compose (update (feedbackColumn j))
          (markedRegionEnvelope α tail last)) q r
      apply compose_nonneg
      · exact liftedGap_nonneg (α g) (hα g)
      · intro u v
        apply compose_nonneg
        · exact update_nonneg (feedbackColumn j)
        · exact ih

/-- An entrywise bound for each empty gap and active step composes through
any fixed separated placement. No boundary-state averaging is assumed. -/
theorem markedRegionKernel_upper (z a : ℝ) (α : ℕ → ℝ)
    (gapOk : ℕ → Prop) (hz : 0 ≤ z) (ha : 0 ≤ a)
    (hα : ∀ g, 0 ≤ α g)
    (hgap : ∀ g, gapOk g → ∀ q r : State22,
      emptyGap z (coordinateWeight coordMap) g q r ≤
        a * liftedGap (α g) q r)
    (hstep : ∀ j : Fin 64, ∀ q r : State22,
      productionImpulse z j q r ≤ update (feedbackColumn j) q r)
    (marks : List (ℕ × Fin 64)) (last : ℕ)
    (hmarks : ∀ p ∈ marks, gapOk p.1) (hlast : gapOk last)
    (q r : State22) :
    markedRegionKernel z marks last q r ≤
      a ^ (marks.length + 1) * markedRegionEnvelope α marks last q r := by
  induction marks generalizing q r with
  | nil =>
      simpa [markedRegionKernel, markedRegionEnvelope] using hgap last hlast q r
  | cons p tail ih =>
      rcases p with ⟨g, j⟩
      have hg : gapOk g := hmarks (g, j) (by simp)
      have htail : ∀ p ∈ tail, gapOk p.1 := by
        intro p hp
        exact hmarks p (by simp [hp])
      have htailbound (u v : State22) :
          markedRegionKernel z tail last u v ≤
            a ^ (tail.length + 1) * markedRegionEnvelope α tail last u v :=
        ih htail u v
      have hinner (u v : State22) :
          compose (productionImpulse z j) (markedRegionKernel z tail last) u v ≤
            a ^ (tail.length + 1) *
              compose (update (feedbackColumn j))
                (markedRegionEnvelope α tail last) u v := by
        have h := compose_scaled_bound
          (productionImpulse z j) (markedRegionKernel z tail last)
          (update (feedbackColumn j)) (markedRegionEnvelope α tail last)
          1 (a ^ (tail.length + 1)) (by norm_num)
          (markedRegionKernel_nonneg z hz tail last)
          (update_nonneg (feedbackColumn j))
          (by intro x y; simpa using hstep j x y) htailbound u v
        simpa using h
      have houter := compose_scaled_bound
        (emptyGap z (coordinateWeight coordMap) g)
        (compose (productionImpulse z j) (markedRegionKernel z tail last))
        (liftedGap (α g))
        (compose (update (feedbackColumn j))
          (markedRegionEnvelope α tail last))
        a (a ^ (tail.length + 1)) ha
        (compose_nonneg (productionImpulse z j)
          (markedRegionKernel z tail last)
          (productionImpulse_nonneg z hz j)
          (markedRegionKernel_nonneg z hz tail last))
        (liftedGap_nonneg (α g) (hα g))
        (hgap g hg) hinner q r
      change compose (emptyGap z (coordinateWeight coordMap) g)
        (compose (productionImpulse z j) (markedRegionKernel z tail last)) q r ≤
        a ^ (((g, j) :: tail).length + 1) *
          compose (liftedGap (α g))
            (compose (update (feedbackColumn j))
              (markedRegionEnvelope α tail last)) q r
      convert houter using 1
      simp only [List.length_cons, pow_succ]
      ring

/-- The paper's finite-length upper comparison for one fixed good
placement. Every gap is nonempty and occupies at most one region. -/
theorem production_markedRegion_upper (θ L : ℝ)
    (hθ : 0 ≤ θ) (hL : 0 < L)
    (marks : List (ℕ × Fin 64)) (last : ℕ)
    (hmarks : ∀ p ∈ marks, 0 < p.1 ∧ (p.1 : ℝ) ≤ L / 64)
    (hlast : 0 < last ∧ (last : ℝ) ≤ L / 64)
    (q r : State22) :
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    markedRegionKernel (Real.exp (-(θ / L))) marks last q r ≤
      (1 + C) ^ (marks.length + 1) *
        markedRegionEnvelope
          (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
          marks last q r := by
  dsimp
  let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
  have hC : 0 ≤ C := by
    dsimp [C]
    positivity
  have hz0 : 0 ≤ Real.exp (-(θ / L)) := (Real.exp_pos _).le
  have hz1 : Real.exp (-(θ / L)) ≤ 1 :=
    Real.exp_le_one_iff.mpr (neg_nonpos.mpr (div_nonneg hθ hL.le))
  apply markedRegionKernel_upper
    (Real.exp (-(θ / L))) (1 + C)
    (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
    (fun g => 0 < g ∧ (g : ℝ) ≤ L / 64)
    hz0 (by linarith) (by intro g; positivity)
  · intro g hg q r
    obtain ⟨hgpos, hgcap⟩ := hg
    obtain ⟨k, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hgpos)
    simpa [C] using
      (production_emptyGap_scaled_sandwich θ L k hθ hL hgcap q r).2
  · intro j q r
    exact (productionImpulse_bounds (Real.exp (-(θ / L))) j q r hz0 hz1).2
  · exact hmarks
  · exact hlast

end Spin.UniformMixer
