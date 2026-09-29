import SpinCodes.Structured.ConcreteFixedProfile
import SpinCodes.Structured.ConcreteFixedRoutedSerialization

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteRoute ConcreteEncoder Filter
attribute [local instance] Classical.propDecidable

def profileChoices {Q b : ℕ} (rows : Fin Q → Finset (Fin b)) : ℝ :=
  ∏ i, (b.choose (rows i).card:ℝ)

theorem profileChoices_pos {Q b : ℕ} (rows : Fin Q → Finset (Fin b)) :
    0<profileChoices rows := Finset.prod_pos (fun i _ => choose_card_pos (rows i))

/-- After independent row shuffles, every subset pattern of the prescribed profile has equal mass. -/
theorem transposed_rowLaw_mass {Q b : ℕ} (rows : Fin Q → Finset (Fin b))
    (A : Fin b → Finset (Fin Q)) :
    ((rowLaw rows).map transpose).p A =
      if rowCount A = (fun i => (rows i).card) then 1/profileChoices rows else 0 := by
  change ((rowLaw rows).map (transposeEquiv Q b)).p A = _
  rw [FinPMF.map_equiv_apply]
  change (∏ i, if (rowCount A i)=(rows i).card then 1/(b.choose (rows i).card:ℝ) else 0) = _
  by_cases h : rowCount A = (fun i => (rows i).card)
  · rw [if_pos h]
    simp only [h,ite_true,profileChoices,prod_div_distrib,prod_const_one]
  · rw [if_neg h]
    obtain ⟨i,hi⟩ := not_forall.mp (fun he => h (funext he))
    exact prod_eq_zero (mem_univ i) (if_neg hi)

theorem rowLaw_pattern_expect {Q b : ℕ} (rows : Fin Q → Finset (Fin b))
    (F : (Fin b → Finset (Fin Q)) → ℝ) :
    (rowLaw rows).expect (fun shuffled => F (transpose shuffled)) =
      (∑ A, if rowCount A = (fun i => (rows i).card) then F A else 0)/profileChoices rows := by
  rw [←FinPMF.expect_map (rowLaw rows) transpose F]
  simp only [FinPMF.expect,transposed_rowLaw_mass,ite_mul,zero_mul,sum_div]
  apply sum_congr rfl
  intro A _
  split_ifs <;> ring

/-- The actual row-shuffle and region-shuffle experiment inherits the fixed-profile margin. -/
theorem actual_rows_regionStream_margin {Q : ℕ} {θ v c δ : ℝ} (hθ : 0≤θ) (hv : 0<v)
    (hc : 0≤c) (hδ : 0<δ) (u : Fin Q → ℝ) (hu : ∀ i, 0<u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : ℕ in atTop, ∀ b : ℕ, ∀ marks : Fin b → (Fin Q ↪ Fin (128*R)),
      ∀ rows : Fin Q → Finset (Fin b),
      (rowLaw rows).expect (fun shuffled =>
        (Spin.piPMF (fun j => shuffleLaw ((transpose shuffled j).map (marks j)))).expect
          (fun regions => inputMoment (Real.exp (-(θ/(128*R)))) (regionStream regions).get ∅)) ≤
      ((c+δ)^b/min 1 v)/(profileChoices rows*profileCoefficient u (fun i => (rows i).card)) := by
  filter_upwards [actual_profile_regionStream_margin hθ hv hc hδ u hu hK] with R hR
  intro b marks rows
  rw [rowLaw_pattern_expect rows (fun A =>
    (Spin.piPMF (fun j => shuffleLaw ((A j).map (marks j)))).expect
      (fun regions => inputMoment (Real.exp (-(θ/(128*R)))) (regionStream regions).get ∅))]
  have hh := div_le_div_of_nonneg_right (hR b marks (fun i => (rows i).card))
    (profileChoices_pos rows).le
  convert hh using 1
  field_simp

#print axioms transposed_rowLaw_mass
#print axioms rowLaw_pattern_expect
#print axioms actual_rows_regionStream_margin
end Spin.Structured.Placement
