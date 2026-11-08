import SpinCodes.Turbo.PhysicalCone22

/-!
# Exact birth measure for the Turbo physical step

An entering zero state is not refreshed by the sampled invertible map.
Its next state is exactly the feedback of the complete four-bit input,
including selected slots whose value is zero. This is the source row
that the retained 11-coordinate producer represents by zero-state,
uniform-density, and expansion-weight-class masses.
-/

noncomputable section

namespace Spin.Turbo

open Spin.UniformMixer

def turboBirthMass22 (j : Fin 17) (z : ℝ)
    (r : State22) : ℝ :=
  (selectedStepLaw22 j).expect (fun x =>
    if stepFeedback (packetStepInput22 x) = r then
      z ^ stepWeight (packetStepInput22 x) 0 else 0)

theorem turboSelectedKernel22_zero_entry (j : Fin 17)
    (z : ℝ) (r : State22) :
    turboSelectedKernel22 j z 0 r = turboBirthMass22 j z r := by
  unfold turboSelectedKernel22 turboBirthMass22
  apply congrArg
  funext x
  simp [stepKernel, update_zero, eq_comm]

theorem turboBirthMass22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (r : State22) :
    0 ≤ turboBirthMass22 j z r := by
  unfold turboBirthMass22
  apply FinPMF.expect_nonneg
  intro x
  by_cases h : stepFeedback (packetStepInput22 x) = r
  · simp only [h, ite_true]
    exact pow_nonneg hz _
  · simp only [h, ite_false]
    exact le_refl 0

theorem turboBirthMass22_total (j : Fin 17) (z : ℝ) :
    (∑ r : State22, turboBirthMass22 j z r) =
      (selectedStepLaw22 j).expect (fun x =>
        z ^ stepWeight (packetStepInput22 x) 0) := by
  simp only [turboBirthMass22]
  rw [← FinPMF.expect_sum]
  apply congrArg
  funext x
  simp [Finset.mem_univ]

/-- Deterministic source-row class assignment. The fallback coordinate
is the arbitrary nonzero class; the checked production spectrum shows
that it is never reached by an actual nonzero state. -/
def physicalBirthIndex22 (r : State22) : Fin 11 :=
  if r = 0 then 0 else
  if coordinateWeight coordMap r = 16 then 3 else
  if coordinateWeight coordMap r = 24 then 4 else
  if coordinateWeight coordMap r = 28 then 5 else
  if coordinateWeight coordMap r = 32 then 6 else
  if coordinateWeight coordMap r = 36 then 7 else
  if coordinateWeight coordMap r = 40 then 8 else
  if coordinateWeight coordMap r = 48 then 9 else
  if coordinateWeight coordMap r = 64 then 10 else 1

theorem physicalBirthIndex22_allows (r : State22) :
    physicalClassAllows22 (physicalBirthIndex22 r) r := by
  by_cases hr : r = 0
  · subst r
    simp [physicalBirthIndex22, physicalClassAllows22]
  · have hw := productionWeight_mem_spectrumWeights r
    simp only [spectrumWeights, Finset.mem_insert,
      Finset.mem_singleton] at hw
    rcases hw with h | h | h | h | h | h | h | h | h
    · exact (hr (productionWeight_zero_only r h)).elim
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]
    · simp [physicalBirthIndex22, physicalClassAllows22,
        physicalClassWeight22, hr, h]

theorem physicalBirthIndex22_ne_two (r : State22) :
    physicalBirthIndex22 r ≠ 2 := by
  unfold physicalBirthIndex22
  split_ifs <;> decide

theorem physicalBirthIndex22_ne_one (r : State22) :
    physicalBirthIndex22 r ≠ 1 := by
  by_cases hr : r = 0
  · simp [physicalBirthIndex22, hr]
  · have hw := productionWeight_mem_spectrumWeights r
    simp only [spectrumWeights, Finset.mem_insert,
      Finset.mem_singleton] at hw
    rcases hw with h | h | h | h | h | h | h | h | h
    · exact (hr (productionWeight_zero_only r h)).elim
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]
    · simp [physicalBirthIndex22, hr, h]

def turboBirthClassMass22 (j : Fin 17) (z : ℝ)
    (i : Fin 11) : ℝ :=
  ∑ r : State22,
    if physicalBirthIndex22 r = i then turboBirthMass22 j z r else 0

theorem turboBirthClassMass22_nonneg (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (i : Fin 11) :
    0 ≤ turboBirthClassMass22 j z i := by
  unfold turboBirthClassMass22
  apply Finset.sum_nonneg
  intro r _
  split_ifs
  · exact turboBirthMass22_nonneg j z hz r
  · exact le_refl 0

theorem turboBirthClassMass22_total (j : Fin 17) (z : ℝ) :
    (∑ i : Fin 11, turboBirthClassMass22 j z i) =
      (∑ r : State22, turboBirthMass22 j z r) := by
  unfold turboBirthClassMass22
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro r _
  simp [Finset.mem_univ]

theorem turboBirthClassMass22_one_zero (j : Fin 17) (z : ℝ) :
    turboBirthClassMass22 j z 1 = 0 := by
  unfold turboBirthClassMass22
  apply Finset.sum_eq_zero
  intro r _
  exact if_neg (physicalBirthIndex22_ne_one r)

theorem turboBirthClassMass22_two_zero (j : Fin 17) (z : ℝ) :
    turboBirthClassMass22 j z 2 = 0 := by
  unfold turboBirthClassMass22
  apply Finset.sum_eq_zero
  intro r _
  exact if_neg (physicalBirthIndex22_ne_two r)

theorem physicalZeroPiece22_supported {b : ℝ}
    {μ : State22 → ℝ}
    (h : physicalPiece22 0 b μ)
    (q : State22) (hq : q ≠ 0) : μ q = 0 := by
  by_contra hn
  have ha := h.2.1 q hn
  exact hq (by simpa [physicalClassAllows22] using ha)

theorem physicalZeroPiece22_mass {b : ℝ}
    {μ : State22 → ℝ}
    (h : physicalPiece22 0 b μ) :
    μ 0 ≤ b := by
  have hsum : (∑ q : State22, μ q) = μ 0 := by
    apply Finset.sum_eq_single 0
    · intro q _ hq
      exact physicalZeroPiece22_supported h q hq
    · intro hzero
      exact (hzero (Finset.mem_univ 0)).elim
  rw [← hsum]
  exact h.2.2.1

def turboBirthParts22 (j : Fin 17) (z b : ℝ)
    (i : Fin 11) (r : State22) : ℝ :=
  if physicalBirthIndex22 r = i then b * turboBirthMass22 j z r else 0

theorem turboBirthParts22_sum (j : Fin 17) (z b : ℝ)
    (r : State22) :
    (∑ i : Fin 11, turboBirthParts22 j z b i r) =
      b * turboBirthMass22 j z r := by
  simp [turboBirthParts22, Finset.mem_univ]

theorem turboBirthParts22_piece (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) (b B : ℝ)
    (hb0 : 0 ≤ b) (hbB : b ≤ B)
    (i : Fin 11) :
    physicalPiece22 i (B * turboBirthClassMass22 j z i)
      (turboBirthParts22 j z b i) := by
  constructor
  · intro r
    unfold turboBirthParts22
    split_ifs
    · exact mul_nonneg hb0 (turboBirthMass22_nonneg j z hz r)
    · exact le_refl 0
  constructor
  · intro r hnonzero
    have hi : physicalBirthIndex22 r = i := by
      by_contra hne
      exact hnonzero (by simp [turboBirthParts22, hne])
    rw [← hi]
    exact physicalBirthIndex22_allows r
  constructor
  · have hsum :
        (∑ r : State22, turboBirthParts22 j z b i r) =
          b * turboBirthClassMass22 j z i := by
        unfold turboBirthParts22 turboBirthClassMass22
        calc
          (∑ r : State22,
            if physicalBirthIndex22 r = i then
              b * turboBirthMass22 j z r else 0) =
            ∑ r : State22,
              b * (if physicalBirthIndex22 r = i then
                turboBirthMass22 j z r else 0) := by
              apply Finset.sum_congr rfl
              intro r _
              split_ifs <;> ring
          _ = _ := (Finset.mul_sum ..).symm
    rw [hsum]
    exact mul_le_mul_of_nonneg_right hbB
      (turboBirthClassMass22_nonneg j z hz i)
  · intro hi r
    have hindex : physicalBirthIndex22 r ≠ i := by
      have hi2 : i = 2 := Fin.ext hi
      rw [hi2]
      exact physicalBirthIndex22_ne_two r
    have hpart : turboBirthParts22 j z b i r = 0 := by
      simp [turboBirthParts22, hindex]
    rw [hpart]
    exact div_nonneg
      (mul_nonneg (hb0.trans hbB)
        (turboBirthClassMass22_nonneg j z hz i))
      (by norm_num)

/-- Exact symbolic birth-row domination for every nonnegative tilt.
The selected row consists of the complete class masses; no endpoint
from the external producer is trusted here. -/
theorem turboSelectedKernel22_zero_row_sound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z)
    (b : ℝ) (μ : State22 → ℝ)
    (hμ : physicalPiece22 0 b μ) :
    ∃ parts : Fin 11 → State22 → ℝ,
      (∀ r, physicalKernelAction22
        (Matrix.of (turboSelectedKernel22 j z)) μ r =
          ∑ i : Fin 11, parts i r) ∧
      ∀ i, physicalPiece22 i
        (b * turboBirthClassMass22 j z i) (parts i) := by
  let c := μ 0
  have hc0 : 0 ≤ c := hμ.1 0
  have hcb : c ≤ b := physicalZeroPiece22_mass hμ
  refine ⟨turboBirthParts22 j z c, ?_, ?_⟩
  · intro r
    have hk : physicalKernelAction22
        (Matrix.of (turboSelectedKernel22 j z)) μ r =
        c * turboBirthMass22 j z r := by
      unfold physicalKernelAction22
      rw [Finset.sum_eq_single 0]
      · exact congrArg (c * ·)
          (turboSelectedKernel22_zero_entry j z r)
      · intro q _ hq
        rw [physicalZeroPiece22_supported hμ q hq]
        ring
      · intro hzero
        exact (hzero (Finset.mem_univ 0)).elim
    rw [hk, turboBirthParts22_sum]
  · intro i
    exact turboBirthParts22_piece j z hz c b hc0 hcb i

#print axioms turboSelectedKernel22_zero_entry
#print axioms turboBirthMass22_nonneg
#print axioms turboBirthMass22_total
#print axioms physicalBirthIndex22_allows
#print axioms physicalBirthIndex22_ne_one
#print axioms turboBirthClassMass22_total
#print axioms turboSelectedKernel22_zero_row_sound

end Spin.Turbo
