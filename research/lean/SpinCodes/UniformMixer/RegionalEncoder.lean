import SpinCodes.UniformMixer.RegionalSupport

/-!
# Stepwise encoder law for a regional input

For an arbitrary 64-bit input group, the output is `Aq+x` and the
feedback is `Aᵀx`. The uniform mixer then acts on the current state.
This definition specializes exactly to the empty-step and one-bit
kernels used by the regional product proof.
-/

noncomputable section

namespace Spin.UniformMixer

def stepWeight (x : Fin 64 → ZMod 2) (q : State22) : ℕ :=
  ∑ j : Fin 64, if coordMap j q + x j = 1 then 1 else 0

def stepFeedback (x : Fin 64 → ZMod 2) : State22 :=
  fun i => ∑ j : Fin 64, x j * coordCoeff j i

def stepKernel (z : ℝ) (x : Fin 64 → ZMod 2) :
    State22 → State22 → ℝ :=
  fun q r => z ^ stepWeight x q * update (stepFeedback x) q r

def singletonInput (j : Fin 64) : Fin 64 → ZMod 2 :=
  fun k => if k = j then 1 else 0

theorem stepWeight_empty (q : State22) :
    stepWeight 0 q = coordinateWeight coordMap q := by
  simp only [stepWeight, coordinateWeight, Pi.zero_apply, add_zero]

theorem stepFeedback_empty : stepFeedback 0 = 0 := by
  funext i
  simp only [stepFeedback, Pi.zero_apply, zero_mul, Finset.sum_const_zero]

theorem stepKernel_empty (z : ℝ) :
    stepKernel z 0 = emptyStep z (coordinateWeight coordMap) := by
  funext q r
  simp [stepKernel, emptyStep, stepWeight_empty, stepFeedback_empty]

theorem stepWeight_singleton (j : Fin 64) (q : State22) :
    stepWeight (singletonInput j) q = activeWeight j q := rfl

theorem stepFeedback_singleton (j : Fin 64) :
    stepFeedback (singletonInput j) = feedbackColumn j := by
  funext i
  simp [stepFeedback, singletonInput, feedbackColumn,
    Finset.sum_ite_eq']

theorem stepKernel_singleton (z : ℝ) (j : Fin 64) :
    stepKernel z (singletonInput j) = productionImpulse z j := by
  funext q r
  simp [stepKernel, productionImpulse, weightedImpulse,
    stepWeight_singleton, stepFeedback_singleton]

/-- Bit `j` of step `i` is bit `64i+j` of the routed region support. -/
def regionBit {R : ℕ} (i : Fin R) (j : Fin 64) : Fin (64 * R) :=
  ⟨64 * i.val + j.val, by have hi := i.isLt; have hj := j.isLt; omega⟩

def regionInput {R : ℕ} (T : Finset (Fin (64 * R)))
    (i : Fin R) : Fin 64 → ZMod 2 :=
  fun j => if regionBit i j ∈ T then 1 else 0

def regionStep {R : ℕ} (x : Fin (64 * R)) : Fin R :=
  ⟨x.val / 64, by have hx := x.isLt; omega⟩

def regionLocal {R : ℕ} (x : Fin (64 * R)) : Fin 64 :=
  ⟨x.val % 64, Nat.mod_lt _ (by omega)⟩

theorem regionBit_step {R : ℕ} (i : Fin R) (j : Fin 64) :
    regionStep (regionBit i j) = i := by
  ext
  simp only [regionStep, regionBit]
  have hj := j.isLt
  omega

theorem regionBit_local {R : ℕ} (i : Fin R) (j : Fin 64) :
    regionLocal (regionBit i j) = j := by
  ext
  simp only [regionLocal, regionBit]
  have hj := j.isLt
  omega

theorem regionBit_reconstruct {R : ℕ} (x : Fin (64 * R)) :
    regionBit (regionStep x) (regionLocal x) = x := by
  ext
  simp only [regionBit, regionStep, regionLocal]
  exact Nat.div_add_mod x.val 64

theorem regionInput_empty_step {R : ℕ}
    (T : Finset (Fin (64 * R))) (i : Fin R)
    (hi : ∀ x ∈ T, regionStep x ≠ i) :
    regionInput T i = 0 := by
  funext j
  have hnot : regionBit i j ∉ T := by
    intro h
    exact (hi (regionBit i j) h) (regionBit_step i j)
  simp [regionInput, hnot]

theorem regionInput_singleton_step {R : ℕ}
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (x : Fin (64 * R)) (hx : x ∈ T) :
    regionInput T (regionStep x) = singletonInput (regionLocal x) := by
  funext j
  by_cases hj : j = regionLocal x
  · subst j
    have hbit : regionBit (regionStep x) (regionLocal x) ∈ T := by
      rw [regionBit_reconstruct]
      exact hx
    simp [regionInput, singletonInput, hbit]
  · have hnot : regionBit (regionStep x) j ∉ T := by
      intro hy
      have hneq : regionBit (regionStep x) j ≠ x := by
        intro heq
        have heq' := congrArg regionLocal heq
        rw [regionBit_local] at heq'
        exact hj heq'
      have hsep := goodRegionalSites_step_separation T hgood
        (regionBit (regionStep x) j) x hy hx hneq
      have hs : (regionBit (regionStep x) j).val / 64 = x.val / 64 := by
        have h := congrArg Fin.val (regionBit_step (regionStep x) j)
        exact h
      simp [hs] at hsep
    simp [regionInput, singletonInput, hj, hnot]

theorem stepKernel_empty_region {R : ℕ} (z : ℝ)
    (T : Finset (Fin (64 * R))) (i : Fin R)
    (hi : ∀ x ∈ T, regionStep x ≠ i) :
    stepKernel z (regionInput T i) =
      emptyStep z (coordinateWeight coordMap) := by
  rw [regionInput_empty_step T i hi]
  exact stepKernel_empty z

theorem stepKernel_singleton_region {R : ℕ} (z : ℝ)
    (T : Finset (Fin (64 * R))) (hgood : ¬ badRegionalSites T)
    (x : Fin (64 * R)) (hx : x ∈ T) :
    stepKernel z (regionInput T (regionStep x)) =
      productionImpulse z (regionLocal x) := by
  rw [regionInput_singleton_step T hgood x hx]
  exact stepKernel_singleton z (regionLocal x)

end Spin.UniformMixer
