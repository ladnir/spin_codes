import SpinCodes.UniformMixer.ScalableFixedTwoState

/-! The full-state continuum weight action is exactly two scalar rows.
The resulting fixed-profile bound needs only a two-state numerical gate. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.Placement Filter

def scalableStateWeight (v : ℝ) (q : State22) : ℝ :=
  if q = 0 then 1 else v

theorem scalableStateWeight_pos {v : ℝ} (hv : 0 < v)
    (q : State22) : 0 < scalableStateWeight v q := by
  unfold scalableStateWeight
  split_ifs <;> positivity

theorem scalableStateWeight_min {v : ℝ} (q : State22) :
    min 1 v ≤ scalableStateWeight v q := by
  unfold scalableStateWeight
  split_ifs
  · exact min_le_left _ _
  · exact min_le_right _ _

theorem projectM_stateWeight_sum (v : ℝ) (b : Bool) :
    (∑ r : State22, projectM b r * scalableStateWeight v r) =
      if b then v else 1 := by
  cases b
  · simpa [projectM, scalableStateWeight] using
      (project_zero_const (S := State22) (scalableStateWeight v))
  · have h : ∀ r : State22, r ≠ 0 → scalableStateWeight v r = v := by
      intro r hr
      simp [scalableStateWeight, hr]
    simpa [projectM] using
      (project_live_const (S := State22) (scalableStateWeight v) v h)

theorem lifted_stateWeight_action (K : Matrix Bool Bool ℝ)
    (v : ℝ) (q : State22) :
    (∑ r : State22,
      (liftM * K * projectM) q r * scalableStateWeight v r) =
    K (if q = 0 then false else true) false +
      K (if q = 0 then false else true) true * v := by
  calc
    _ = ∑ r : State22,
        (∑ b : Bool, (liftM * K) q b * projectM b r) *
          scalableStateWeight v r := by
      simp only [Matrix.mul_apply]
    _ = ∑ b : Bool,
        (liftM * K) q b *
          (∑ r : State22, projectM b r * scalableStateWeight v r) := by
      simp only [Finset.sum_mul]
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro b _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro r _
      ring
    _ = ∑ b : Bool,
        K (if q = 0 then false else true) b *
          (if b then v else 1) := by
      simp only [liftM_mul_apply, projectM_stateWeight_sum]
    _ = _ := by simp; ring

def scalableTwoRowNormLe (v c : ℝ) (K : Matrix Bool Bool ℝ) : Prop :=
  K false false + K false true * v ≤ c ∧
    K true false + K true true * v ≤ c * v

theorem scalableFugacityContinuum_weighted_bound {Q : ℕ}
    (θ : ℝ) (u : Fin Q → ℝ) {v c : ℝ}
    (hK : scalableTwoRowNormLe v c (scalableFugacityTwo θ u))
    (q : State22) :
    (∑ r : State22,
      scalableFugacityContinuum θ u q r * scalableStateWeight v r) ≤
      c * scalableStateWeight v q := by
  rw [scalableFugacityContinuum_eq_lift_two]
  rw [lifted_stateWeight_action]
  by_cases hq : q = 0
  · simpa [scalableStateWeight, hq, scalableTwoRowNormLe] using hK.1
  · simpa [scalableStateWeight, hq, scalableTwoRowNormLe] using hK.2

theorem scalable_fixed_profile_two_state_margin {Q : ℕ}
    {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 < u i)
    (hK : scalableTwoRowNormLe v c (scalableFugacityTwo θ u)) :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin (bsched m) → (Fin Q ↪ Fin (Lsched m)),
        ∀ rowWeights : Fin Q → ℕ,
          (∑ A : Fin (bsched m) → Finset (Fin Q),
            if rowCount A = rowWeights then
              (ConcreteRoute.regionKernel (fun j => (A j).map (marks j))).expect
                (fun regions =>
                  ∑ r : State22,
                    inputListKernel (Real.exp (-(θ / (Lsched m : ℝ))))
                      (List.ofFn (scalableFlatInput m regions)) 0 r)
            else 0) ≤
          ((c + δ) ^ bsched m / min 1 v) /
            profileCoefficient u rowWeights := by
  have hm : 0 < min 1 v := lt_min (by norm_num) hv
  filter_upwards [scalable_fixed_profile_margin hθ hc hδ u hu
    (scalableStateWeight v) (scalableStateWeight_pos hv)
    hm scalableStateWeight_min
    (scalableFugacityContinuum_weighted_bound θ u hK)] with m hmargin
  intro marks rowWeights
  simpa [scalableStateWeight] using
    hmargin marks rowWeights

end Spin.UniformMixer
