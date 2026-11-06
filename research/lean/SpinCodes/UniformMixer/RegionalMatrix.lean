import SpinCodes.UniformMixer.RegionalProduct
import Mathlib.Data.Matrix.Mul

/-!
# Full-state and two-state matrix interfaces

The lift and projection retain the state at region boundaries. Between
two nonempty gaps, a one-bit step reduces exactly to the same two-state
matrix for every local input coordinate.
-/

noncomputable section

namespace Spin.UniformMixer

def liftM : Matrix State22 Bool ℝ := classLift
def projectM : Matrix Bool State22 ℝ := classProject
def gapM (α : ℝ) : Matrix State22 State22 ℝ := liftedGap α
def impulseM (j : Fin 64) : Matrix State22 State22 ℝ :=
  update (feedbackColumn j)

def coarseD (α : ℝ) : Matrix Bool Bool ℝ :=
  fun b d => if b = d then (if b then α else 1) else 0

def coarseP : Matrix Bool Bool ℝ :=
  fun b d => if b then (if d then 1 - 1 / (liveCount State22 : ℝ)
      else 1 / (liveCount State22 : ℝ))
    else if d then 1 else 0

theorem project_liftM : projectM * liftM = 1 := by
  ext b d
  change (∑ q : State22, classProject b q * classLift q d) =
    if b = d then 1 else 0
  exact project_lift (S := State22) b d

theorem coarseP_eq_classImpulse (j : Fin 64) (b d : Bool) :
    coarseP b d = classImpulse (feedbackColumn j) b d := by
  cases b <;> cases d
  · exact (classImpulse_zero_zero (feedbackColumn j) (feedbackColumn_nonzero j)).symm
  · exact (classImpulse_zero_live (feedbackColumn j) (feedbackColumn_nonzero j)).symm
  · exact (classImpulse_live_zero (feedbackColumn j) (feedbackColumn_nonzero j)).symm
  · exact (classImpulse_live_live (feedbackColumn j) (feedbackColumn_nonzero j)).symm

theorem project_impulse_lift (j : Fin 64) :
    projectM * impulseM j * liftM = coarseP := by
  rw [Matrix.mul_assoc]
  ext b d
  change classImpulse (feedbackColumn j) b d = coarseP b d
  exact (coarseP_eq_classImpulse j b d).symm

theorem gapM_factor (α : ℝ) :
    gapM α = liftM * coarseD α * projectM := by
  ext q r
  change liftedGap α q r = (liftM * coarseD α * projectM) q r
  rw [liftedGap_as_JDR]
  change (∑ b : Bool, classLift q b * (if b then α else 1) *
      classProject b r) =
    ∑ d : Bool, (∑ b : Bool, classLift q b * coarseD α b d) *
      classProject d r
  simp [coarseD]

def coarseRegionProduct (α : ℕ → ℝ) :
    List (ℕ × Fin 64) → ℕ → Matrix Bool Bool ℝ
  | [], last => coarseD (α last)
  | (g, _) :: tail, last =>
      coarseD (α g) * coarseP * coarseRegionProduct α tail last

theorem compose_eq_matrixMul (K L : Matrix State22 State22 ℝ) :
    compose K L = K * L := by
  funext q r
  rfl

theorem matrix_sandwich
    (J : Matrix State22 Bool ℝ) (R : Matrix Bool State22 ℝ)
    (D A : Matrix Bool Bool ℝ) (V : Matrix State22 State22 ℝ) :
    (J * D * R) * V * (J * A * R) =
      J * (D * (R * V * J) * A) * R := by
  simp only [Matrix.mul_assoc]

theorem markedRegionEnvelope_factor (α : ℕ → ℝ)
    (marks : List (ℕ × Fin 64)) (last : ℕ) :
    (markedRegionEnvelope α marks last : Matrix State22 State22 ℝ) =
      liftM * coarseRegionProduct α marks last * projectM := by
  induction marks with
  | nil => exact gapM_factor (α last)
  | cons p tail ih =>
      rcases p with ⟨g, j⟩
      change compose (liftedGap (α g))
        (compose (update (feedbackColumn j))
          (markedRegionEnvelope α tail last)) =
        liftM * (coarseD (α g) * coarseP *
          coarseRegionProduct α tail last) * projectM
      let E : Matrix State22 State22 ℝ := markedRegionEnvelope α tail last
      change gapM (α g) * (impulseM j * E) = _
      have hE : E = liftM * coarseRegionProduct α tail last * projectM := ih
      rw [hE, ← Matrix.mul_assoc]
      rw [gapM_factor, matrix_sandwich, project_impulse_lift]

theorem production_markedRegion_coarse_upper (θ L : ℝ)
    (hθ : 0 ≤ θ) (hL : 0 < L)
    (marks : List (ℕ × Fin 64)) (last : ℕ)
    (hmarks : ∀ p ∈ marks, 0 < p.1 ∧ (p.1 : ℝ) ≤ L / 64)
    (hlast : 0 < last ∧ (last : ℝ) ≤ L / 64)
    (q r : State22) :
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    markedRegionKernel (Real.exp (-(θ / L))) marks last q r ≤
      (1 + C) ^ (marks.length + 1) *
        (liftM * coarseRegionProduct
          (fun g => Real.exp (-((θ / L) * liveWeightAverage coordMap)) ^ g)
          marks last * projectM) q r := by
  have h := production_markedRegion_upper θ L hθ hL
    marks last hmarks hlast q r
  simpa only [markedRegionEnvelope_factor] using h

end Spin.UniformMixer
