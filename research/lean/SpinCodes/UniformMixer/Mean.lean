import SpinCodes.UniformMixer.Kernel
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum

/-!
# Live-state mean of binary linear output coordinates

A nonzero binary linear functional takes value one on exactly half of all
states. The zero state contributes zero, so the same numerator applies to
the uniform distribution on nonzero states.
-/

noncomputable section

namespace Spin.UniformMixer

variable {V : Type*} [AddCommGroup V] [Module (ZMod 2) V] [Fintype V]

theorem zmod2_zero_or_one (a : ZMod 2) : a = 0 ∨ a = 1 := by
  fin_cases a
  · exact Or.inl rfl
  · exact Or.inr rfl

omit [Fintype V] in
theorem binary_functional_toggle (ℓ : V →ₗ[ZMod 2] ZMod 2)
    (u v : V) (hu : ℓ u = 1) :
    (ℓ (v + u) = 1 ↔ ℓ v = 0) := by
  rw [map_add, hu]
  constructor
  · intro h
    have h' : ℓ v + (1 : ZMod 2) = 0 + (1 : ZMod 2) := by simpa using h
    exact add_right_cancel h'
  · intro h
    simp [h]

theorem binary_functional_half (ℓ : V →ₗ[ZMod 2] ZMod 2)
    (u : V) (hu : ℓ u = 1) :
    (∑ v : V, if ℓ v = 1 then (1 : ℝ) else 0) =
      (Fintype.card V : ℝ) / 2 := by
  classical
  have hpair :
      (∑ v : V, if ℓ v = 0 then (1 : ℝ) else 0) =
        (∑ v : V, if ℓ v = 1 then (1 : ℝ) else 0) := by
    apply Fintype.sum_equiv (Equiv.addRight u)
    intro v
    change (if ℓ v = 0 then (1 : ℝ) else 0) =
      if ℓ (v + u) = 1 then 1 else 0
    simp only [binary_functional_toggle ℓ u v hu]
  have hpoint (v : V) :
      (if ℓ v = 0 then (1 : ℝ) else 0) +
        (if ℓ v = 1 then 1 else 0) = 1 := by
    rcases zmod2_zero_or_one (ℓ v) with h0 | h1
    · simp [h0]
    · simp [h1]
  have htotal :
      (∑ v : V, if ℓ v = 0 then (1 : ℝ) else 0) +
        (∑ v : V, if ℓ v = 1 then 1 else 0) =
        (Fintype.card V : ℝ) := by
    rw [← Finset.sum_add_distrib]
    simp_rw [hpoint]
    simp
  rw [hpair] at htotal
  linarith

/-- Weight emitted by a family of binary linear output coordinates. -/
def coordinateWeight {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) (v : V) : ℕ :=
  ∑ j : Fin t, if A j v = 1 then 1 else 0

omit [Fintype V] in
theorem coordinateWeight_zero {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) :
    coordinateWeight A 0 = 0 := by
  simp [coordinateWeight]

omit [Fintype V] in
theorem coordinateWeight_cast {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) (v : V) :
    (coordinateWeight A v : ℝ) =
      ∑ j : Fin t, if A j v = 1 then (1 : ℝ) else 0 := by
  simp [coordinateWeight]

theorem coordinateWeight_sum {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2)
    (hA : ∀ j : Fin t, ∃ u : V, A j u = 1) :
    (∑ v : V, (coordinateWeight A v : ℝ)) =
      (t : ℝ) * (Fintype.card V : ℝ) / 2 := by
  classical
  have heach (j : Fin t) :
      (∑ v : V, if A j v = 1 then (1 : ℝ) else 0) =
        (Fintype.card V : ℝ) / 2 := by
    obtain ⟨u, hu⟩ := hA j
    exact binary_functional_half (A j) u hu
  calc
    (∑ v : V, (coordinateWeight A v : ℝ)) =
        ∑ j : Fin t, ∑ v : V,
          if A j v = 1 then (1 : ℝ) else 0 := by
            simp_rw [coordinateWeight_cast]
            rw [Finset.sum_comm]
    _ = ∑ _j : Fin t, (Fintype.card V : ℝ) / 2 := by
          apply Finset.sum_congr rfl
          intro j _
          exact heach j
    _ = (t : ℝ) * (Fintype.card V : ℝ) / 2 := by
          simp
          ring

omit [Fintype V] in
theorem coordinateWeight_le {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) (v : V) :
    coordinateWeight A v ≤ t := by
  unfold coordinateWeight
  calc
    (∑ j : Fin t, if A j v = 1 then 1 else 0) ≤
        ∑ _j : Fin t, 1 := by
          apply Finset.sum_le_sum
          intro j _
          by_cases hj : A j v = 1 <;> simp [hj]
    _ = t := by simp

/-- The output-weight mean under the uniform distribution on nonzero states.
The zero state has zero output weight, so its omission affects only the
denominator. -/
def liveWeightAverage {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) : ℝ :=
  (∑ v : V, (coordinateWeight A v : ℝ)) / (liveCount V : ℝ)

theorem liveWeightAverage_erase [DecidableEq V] {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2) :
    liveWeightAverage A =
      ((Finset.univ : Finset V).erase 0).sum
        (fun v => (coordinateWeight A v : ℝ)) / (liveCount V : ℝ) := by
  classical
  have hsum := Finset.add_sum_erase (Finset.univ : Finset V)
    (fun v => (coordinateWeight A v : ℝ)) (Finset.mem_univ (0 : V))
  simp only [coordinateWeight_zero, Nat.cast_zero, zero_add] at hsum
  rw [liveWeightAverage, hsum]

theorem liveWeightAverage_formula {t : ℕ}
    (A : Fin t → V →ₗ[ZMod 2] ZMod 2)
    (hA : ∀ j : Fin t, ∃ u : V, A j u = 1) :
    liveWeightAverage A =
      (t : ℝ) * ((Fintype.card V : ℝ) / (2 * (liveCount V : ℝ))) := by
  rw [liveWeightAverage, coordinateWeight_sum A hA]
  ring

theorem binary_vector_liveWeightAverage_formula (m t : ℕ)
    (A : Fin t → (Fin m → ZMod 2) →ₗ[ZMod 2] ZMod 2)
    (hA : ∀ j : Fin t, ∃ u, A j u = 1) :
    liveWeightAverage A =
      (t : ℝ) * (((2 ^ m : ℕ) : ℝ) /
        (2 * (((2 ^ m - 1 : ℕ)) : ℝ))) := by
  rw [liveWeightAverage_formula A hA, card_binary_vector,
    liveCount_binary_vector]

theorem state22_liveWeightAverage_formula
    (A : Fin 64 → State22 →ₗ[ZMod 2] ZMod 2)
    (hA : ∀ j : Fin 64, ∃ u, A j u = 1) :
    liveWeightAverage A = 64 * ((2 ^ 21 : ℝ) / (2 ^ 22 - 1)) := by
  rw [binary_vector_liveWeightAverage_formula 22 64 A hA]
  norm_num

end Spin.UniformMixer
