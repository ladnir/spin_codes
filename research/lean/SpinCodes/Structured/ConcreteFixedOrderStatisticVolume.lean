import SpinCodes.Structured.ConcreteFixedLargeSpacings
import Mathlib.MeasureTheory.Measure.Haar.NormedSpace
import Mathlib.MeasureTheory.Integral.Prod

noncomputable section
namespace Spin.Structured.Placement
open MeasureTheory Set

def orderedBetween (n : ℕ) (lo hi : ℝ) : Set (Fin n→ℝ) :=
  {x | (∀ i, lo≤x i ∧ x i<hi) ∧ StrictMono x}

theorem orderedBetween_measurable (n : ℕ) (lo hi : ℝ) : MeasurableSet (orderedBetween n lo hi) := by
  unfold orderedBetween StrictMono
  simp only [Set.ofPred_and, Set.ofPred_forall]
  exact MeasurableSet.inter
    (MeasurableSet.iInter (fun i => (measurableSet_le measurable_const (measurable_pi_apply i)).inter
      (measurableSet_lt (measurable_pi_apply i) measurable_const)))
    (MeasurableSet.iInter (fun i => MeasurableSet.iInter (fun j => MeasurableSet.iInter (fun _ =>
      measurableSet_lt (measurable_pi_apply i) (measurable_pi_apply j)))))

theorem orderedBetween_normalized {n : ℕ} {lo hi : ℝ} (h : lo<hi) (x : Fin n→ℝ) :
    (hi-lo)⁻¹ • (x-(fun _=>lo)) ∈ orderedSiteDomain n ↔ x∈orderedBetween n lo hi := by
  have hw : 0<hi-lo := sub_pos.mpr h
  change ((∀i,0≤(hi-lo)⁻¹*(x i-lo) ∧ (hi-lo)⁻¹*(x i-lo)<1) ∧
    StrictMono (fun i => (hi-lo)⁻¹*(x i-lo))) ↔ _
  simp only [←div_eq_inv_mul]
  constructor
  · rintro ⟨hx,hm⟩
    constructor
    · intro i
      have hl := (le_div_iff₀ hw).mp (hx i).1
      have hu := (div_lt_one hw).mp (hx i).2
      constructor <;> linarith
    · intro i j hij
      have hh := (div_lt_div_iff_of_pos_right hw).mp (hm hij)
      linarith
  · rintro ⟨hx,hm⟩
    constructor
    · intro i
      exact ⟨(le_div_iff₀ hw).mpr (by have := (hx i).1; linarith),
        (div_lt_one hw).mpr (by have := (hx i).2; linarith)⟩
    · intro i j hij
      apply (div_lt_div_iff_of_pos_right hw).mpr
      exact sub_lt_sub_right (hm hij) lo

/-- Exact volume of a translated and scaled ordered-site simplex. -/
theorem orderedBetween_integral_one (n : ℕ) {lo hi : ℝ} (h : lo<hi) :
    (∫ _x in orderedBetween n lo hi, (1:ℝ))=(hi-lo)^n/(n.factorial:ℝ) := by
  let f : (Fin n→ℝ)→ℝ := (orderedSiteDomain n).indicator (fun _ => (1:ℝ))
  have he (x : Fin n→ℝ) : f ((hi-lo)⁻¹ • (x-(fun _=>lo)))=
      (orderedBetween n lo hi).indicator (fun _ => (1:ℝ)) x := by
    by_cases hx : x∈orderedBetween n lo hi
    · simp [f,hx,(orderedBetween_normalized h x).mpr hx]
    · simp [f,hx,mt (orderedBetween_normalized h x).mp hx]
  rw [←integral_indicator (orderedBetween_measurable n lo hi)]
  simp_rw [←he]
  rw [integral_sub_right_eq_self (fun x : Fin n→ℝ => f ((hi-lo)⁻¹ • x)) (fun _=>lo)]
  have hs := Measure.integral_comp_inv_smul_of_nonneg (volume:Measure (Fin n→ℝ)) f (sub_pos.mpr h).le
  simp only [Module.finrank_fintype_fun_eq_card,Fintype.card_fin,smul_eq_mul] at hs
  rw [hs]
  rw [show (∫ x, f x)=(∫ x in orderedSiteDomain n,(1:ℝ)) from integral_indicator (orderedSiteDomain_measurable n)]
  rw [orderedSite_integral_const_one]
  ring

#print axioms orderedBetween_integral_one
end Spin.Structured.Placement
