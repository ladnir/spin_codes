import SpinCodes.Turbo.PhysicalBirth22

/-!
# One joint capped-density birth decomposition

The retained producer's classwise minimum uses two bounds on the SAME
residual `max (μ-d) 0`, after the SAME uniform cap `min μ d` is removed.
These law-independent lemmas apply to exact-nonzero as well as selected-slot
birth laws. They authenticate no external census or numerical endpoint.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer
attribute [local instance] Classical.propDecidable

def cappedBirthUniform22 (μ : State22 → ℝ) (d : ℝ) (r : State22) : ℝ :=
  if r = 0 then 0 else min (μ r) d

def cappedBirthResidual22 (μ : State22 → ℝ) (d : ℝ) (r : State22) : ℝ :=
  if r = 0 then 0 else max (μ r - d) 0

def cappedBirthClassPart22 (μ : State22 → ℝ) (d : ℝ)
    (i : Fin 11) (r : State22) : ℝ :=
  if physicalBirthIndex22 r = i then cappedBirthResidual22 μ d r else 0

def cappedBirthClassBound22 (μ U : State22 → ℝ) (d : ℝ) (i : Fin 11) : ℝ :=
  min (∑ r : State22, if r ≠ 0 ∧ physicalBirthIndex22 r = i then μ r else 0)
    (∑ r : State22,
      if r ≠ 0 ∧ physicalBirthIndex22 r = i then max (U r - d) 0 else 0)

theorem cappedBirth22_decomposition (μ : State22 → ℝ) (d : ℝ) (r : State22) :
    μ r = (if r = 0 then μ 0 else 0) + cappedBirthUniform22 μ d r +
      ∑ i : Fin 11, cappedBirthClassPart22 μ d i r := by
  have hsum : (∑ i : Fin 11, cappedBirthClassPart22 μ d i r) =
      cappedBirthResidual22 μ d r := by
    simp [cappedBirthClassPart22]
  rw [hsum]
  by_cases hr : r = 0
  · subst r; simp [cappedBirthUniform22, cappedBirthResidual22]
  · simp only [hr, ite_false, cappedBirthUniform22, cappedBirthResidual22]
    by_cases hd : μ r ≤ d
    · rw [min_eq_left hd, max_eq_right (sub_nonpos.mpr hd)]
      ring
    · have hd' : d ≤ μ r := le_of_not_ge hd
      rw [min_eq_right hd', max_eq_left (sub_nonneg.mpr hd')]
      ring

private theorem uniformCap_mass_le {S : Type*} [Fintype S] [DecidableEq S] [Zero S]
    (μ : S → ℝ) (d : ℝ) :
    (∑ r : S, if r = 0 then 0 else min (μ r) d) ≤ (liveCount S : ℝ) * d := by
  have h0 : (if (0 : S) = 0 then (0 : ℝ) else min (μ 0) d) = 0 := if_pos rfl
  rw [← Finset.sum_erase_add _ _ (Finset.mem_univ (0 : S)), h0, add_zero]
  calc
    (∑ r ∈ (Finset.univ : Finset S).erase 0, if r = 0 then 0 else min (μ r) d) ≤
        ∑ _r ∈ (Finset.univ : Finset S).erase 0, d := by
      apply Finset.sum_le_sum
      intro r hr
      rw [if_neg (Finset.mem_erase.mp hr).1]
      exact min_le_right _ _
    _ = (liveCount S : ℝ) * d := by
      rw [Finset.sum_const, nsmul_eq_mul,
        Finset.card_erase_of_mem (Finset.mem_univ (0 : S)), Finset.card_univ]
      rfl

theorem cappedBirthUniform22_piece (μ : State22 → ℝ) (d : ℝ)
    (hμ : ∀ r, 0 ≤ μ r) (hd : 0 ≤ d) :
    physicalPiece22 2 ((liveCount State22 : ℝ) * d) (cappedBirthUniform22 μ d) := by
  constructor
  · intro r
    unfold cappedBirthUniform22
    split_ifs
    · exact le_refl 0
    · exact le_min (hμ r) hd
  constructor
  · intro r hr
    have hne : r ≠ 0 := by
      intro h; subst r; exact hr (by simp [cappedBirthUniform22])
    simpa [physicalClassAllows22] using hne
  constructor
  · exact uniformCap_mass_le μ d
  · intro _ r
    have hN : (liveCount State22 : ℝ) = (2 : ℝ)^22 - 1 := by
      rw [state22_liveCount]
      norm_num
    rw [hN, mul_div_cancel_left₀ d (by norm_num : (2 : ℝ)^22 - 1 ≠ 0)]
    unfold cappedBirthUniform22
    split_ifs
    · exact hd
    · exact min_le_right _ _
/-- The baseline class mass and the upper-census excess both bound one residual. -/
theorem cappedBirthClassPart22_mass_le (μ U : State22 → ℝ) (d : ℝ)
    (hμ : ∀ r, 0 ≤ μ r) (hU : ∀ r, μ r ≤ U r) (hd : 0 ≤ d)
    (i : Fin 11) :
    (∑ r : State22, cappedBirthClassPart22 μ d i r) ≤
      cappedBirthClassBound22 μ U d i := by
  apply le_min
  · apply Finset.sum_le_sum
    intro r _
    by_cases hr : r = 0
    · subst r; simp [cappedBirthClassPart22, cappedBirthResidual22]
    · by_cases hi : physicalBirthIndex22 r = i
      · simp only [cappedBirthClassPart22, cappedBirthResidual22, hr, hi,
          ite_true, ite_false, ne_eq, not_false_eq_true, true_and]
        exact max_le (sub_le_self _ hd) (hμ r)
      · simp [cappedBirthClassPart22, cappedBirthResidual22, hi]
  · apply Finset.sum_le_sum
    intro r _
    by_cases hr : r = 0
    · subst r; simp [cappedBirthClassPart22, cappedBirthResidual22]
    · by_cases hi : physicalBirthIndex22 r = i
      · simp only [cappedBirthClassPart22, cappedBirthResidual22, hr, hi,
          ite_true, ite_false, ne_eq, not_false_eq_true, true_and]
        exact max_le_max (sub_le_sub_right (hU r) d) le_rfl
      · simp [cappedBirthClassPart22, cappedBirthResidual22, hi]

/-- Every class bound belongs to the complete, jointly selected capped row. -/
theorem cappedBirthClassPart22_piece (μ U : State22 → ℝ) (d : ℝ)
    (hμ : ∀ r, 0 ≤ μ r) (hU : ∀ r, μ r ≤ U r) (hd : 0 ≤ d)
    (i : Fin 11) :
    physicalPiece22 i (cappedBirthClassBound22 μ U d i)
      (cappedBirthClassPart22 μ d i) := by
  have hnonneg (r : State22) : 0 ≤ cappedBirthClassPart22 μ d i r := by
    simp only [cappedBirthClassPart22, cappedBirthResidual22]
    split_ifs <;> positivity
  have hmass := cappedBirthClassPart22_mass_le μ U d hμ hU hd i
  refine ⟨hnonneg, ?_, hmass, ?_⟩
  · intro r hr
    have hi : physicalBirthIndex22 r = i := by
      by_contra h
      exact hr (by simp [cappedBirthClassPart22, h])
    rw [← hi]
    exact physicalBirthIndex22_allows r
  · intro hi r
    have hi' : i = 2 := Fin.ext hi
    have hz : cappedBirthClassPart22 μ d i r = 0 := by
      simp [cappedBirthClassPart22, hi', physicalBirthIndex22_ne_two]
    rw [hz]
    exact div_nonneg ((Finset.sum_nonneg (fun q _ => hnonneg q)).trans hmass)
      (by norm_num)

/-- The birth measure under the actual exact-nonzero physical input law. -/
def turboExactBirthMass22 (j : Fin 17) (z : ℝ) (r : State22) : ℝ :=
  (packetStepLaw22 j).expect (fun x =>
    if stepFeedback (packetStepInput22 x) = r then
      z ^ stepWeight (packetStepInput22 x) 0 else 0)

theorem turboPhysicalKernel22_zero_entry (j : Fin 17) (z : ℝ) (r : State22) :
    turboPhysicalKernel22 j z 0 r = turboExactBirthMass22 j z r := by
  unfold turboPhysicalKernel22 turboExactBirthMass22
  apply congrArg
  funext x
  simp [stepKernel, update_zero, eq_comm]

theorem turboExactBirthMass22_nonneg (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (r : State22) : 0 ≤ turboExactBirthMass22 j z r := by
  apply FinPMF.expect_nonneg
  intro x
  split_ifs
  · exact pow_nonneg hz _
  · exact le_refl 0

/-- One cap simultaneously gives the uniform and every residual-class piece
for the actual exact-occupancy zero-source row. The upper function `U` still
requires authentication from the retained census. -/
theorem turboExactBirth22_capped_pieces (j : Fin 17) (z : ℝ) (hz : 0 ≤ z)
    (d : ℝ) (hd : 0 ≤ d) (U : State22 → ℝ)
    (hU : ∀ r, turboExactBirthMass22 j z r ≤ U r) :
    physicalPiece22 2 ((liveCount State22 : ℝ) * d)
      (cappedBirthUniform22 (turboExactBirthMass22 j z) d) ∧
    (∀ i, physicalPiece22 i
      (cappedBirthClassBound22 (turboExactBirthMass22 j z) U d i)
      (cappedBirthClassPart22 (turboExactBirthMass22 j z) d i)) := by
  exact ⟨cappedBirthUniform22_piece _ d (turboExactBirthMass22_nonneg j z hz) hd,
    fun i => cappedBirthClassPart22_piece _ U d
      (turboExactBirthMass22_nonneg j z hz) hU hd i⟩
#print axioms cappedBirth22_decomposition
#print axioms cappedBirthUniform22_piece
#print axioms cappedBirthClassPart22_piece
end Spin.Turbo
