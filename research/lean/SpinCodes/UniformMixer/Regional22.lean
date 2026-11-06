import SpinCodes.UniformMixer.Tilt

/-!
# One-bit regional step for the production 22-state map

An active input coordinate `j` changes the emitted 64-bit vector at
coordinate `j` and feeds back column `j` of `Aᵀ`. The kernel below uses
the same update law as an empty step, with that nonzero feedback value.
-/

noncomputable section

namespace Spin.UniformMixer

/-- Column `j` of the paper's feedback map `C=Aᵀ`. -/
def feedbackColumn (j : Fin 64) : State22 :=
  fun i => coordCoeff j i

theorem feedbackColumn_nonzero (j : Fin 64) : feedbackColumn j ≠ 0 := by
  obtain ⟨i, hi⟩ := every_output_covered j
  intro h
  have h' := congrFun h i
  simp [feedbackColumn, coordCoeff, hi] at h'

/-- Output weight of a step with a single active bit at coordinate `j`. -/
def activeWeight (j : Fin 64) (q : State22) : ℕ :=
  ∑ k : Fin 64,
    if coordMap k q + (if k = j then (1 : ZMod 2) else 0) = 1
      then 1 else 0

theorem activeWeight_zero (j : Fin 64) : activeWeight j 0 = 1 := by
  simp [activeWeight]

theorem activeWeight_le (j : Fin 64) (q : State22) :
    activeWeight j q ≤ 64 := by
  unfold activeWeight
  calc
    (∑ k : Fin 64,
      if coordMap k q + (if k = j then (1 : ZMod 2) else 0) = 1
        then 1 else 0) ≤ ∑ _k : Fin 64, 1 := by
          apply Finset.sum_le_sum
          intro k _
          split_ifs <;> omega
    _ = 64 := by simp

/-- Weighted full-state kernel for one active bit. -/
def productionImpulse (z : ℝ) (j : Fin 64) : State22 → State22 → ℝ :=
  weightedImpulse z (activeWeight j) (feedbackColumn j)

theorem productionImpulse_bounds (z : ℝ) (j : Fin 64)
    (q r : State22) (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    z ^ 64 * update (feedbackColumn j) q r ≤
      productionImpulse z j q r ∧
    productionImpulse z j q r ≤ update (feedbackColumn j) q r := by
  exact weightedImpulse_bounds z (activeWeight j) (feedbackColumn j)
    q r 64 hz0 hz1 (activeWeight_le j q)

/-- Every input coordinate induces the same zero/live transition matrix. -/
theorem productionImpulse_class_matrix (j : Fin 64) :
    ((classImpulse (feedbackColumn j) false false,
      classImpulse (feedbackColumn j) false true),
     (classImpulse (feedbackColumn j) true false,
      classImpulse (feedbackColumn j) true true)) =
      (((0 : ℝ), 1),
       (1 / (liveCount State22 : ℝ),
        1 - 1 / (liveCount State22 : ℝ))) :=
  classImpulse_matrix (feedbackColumn j) (feedbackColumn_nonzero j)

end Spin.UniformMixer
