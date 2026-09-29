import SpinCodes.Structured.ConcreteFixedOneMoment
import Mathlib.MeasureTheory.Integral.Prod

/-! Converting the ordered two-site integral to the triangular gap integral. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory
attribute [local instance] Classical.propDecidable

theorem orderedSiteDomain_two (x : Fin 2 → ℝ) :
    x ∈ orderedSiteDomain 2 ↔ 0 ≤ x 0 ∧ x 0 < x 1 ∧ x 1 < 1 := by
  constructor
  · intro hx; exact ⟨(hx.1 0).1, hx.2 (by decide : (0:Fin 2) < 1), (hx.1 1).2⟩
  · rintro ⟨h0,h01,h1⟩
    constructor
    · intro i; fin_cases i <;> simp only [Fin.zero_eta, Fin.mk_one]
      · exact ⟨h0, h01.trans h1⟩
      · exact ⟨h0.trans h01.le,h1⟩
    · apply Fin.strictMono_iff_lt_succ.mpr
      intro i; fin_cases i; exact h01

theorem orderedSite_indicator_two (F : (Fin 2 → ℝ) → ℝ) (u v : ℝ) :
    (orderedSiteDomain 2).indicator F ![u,v] =
      (Set.Ico (0:ℝ) 1).indicator (fun u => (Set.Ioo u 1).indicator (fun v => F ![u,v]) v) u := by
  have he : (![u,v] : Fin 2 → ℝ) ∈ orderedSiteDomain 2 ↔
      u ∈ Set.Ico (0:ℝ) 1 ∧ v ∈ Set.Ioo u 1 := by
    rw [orderedSiteDomain_two]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Set.mem_Ico, Set.mem_Ioo]
    constructor
    · rintro ⟨h0,huv,hv⟩; exact ⟨⟨h0,huv.trans hv⟩,huv,hv⟩
    · rintro ⟨⟨h0,_⟩,huv,hv⟩; exact ⟨h0,huv,hv⟩
  simp only [Set.indicator, he]
  split_ifs <;> simp_all

theorem orderedSite_integral_two (F : (Fin 2 → ℝ) → ℝ) (hF : Continuous F) :
    (∫ x in orderedSiteDomain 2, F x) =
      ∫ u in (0:ℝ)..1, ∫ v in u..1, F ![u,v] := by
  let e := MeasurableEquiv.piFinTwo (fun _ : Fin 2 => ℝ)
  let G : ℝ × ℝ → ℝ := fun p => (orderedSiteDomain 2).indicator F (e.symm p)
  have hFD : IntegrableOn F (orderedSiteDomain 2) := by
    apply hF.integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun i => (hx.1 i).1, fun i => (hx.1 i).2.le⟩
  have hG : Integrable G (volume.prod volume) := by
    have h := (volume_preserving_piFinTwo (fun _ : Fin 2 => ℝ)).symm.integrable_comp_of_integrable
      (hFD.integrable_indicator (orderedSiteDomain_measurable 2))
    exact h
  have hcomp := (volume_preserving_piFinTwo (fun _ : Fin 2 => ℝ)).symm.integral_comp'
    ((orderedSiteDomain 2).indicator F)
  rw [← integral_indicator (orderedSiteDomain_measurable 2)]
  have heq (u v : ℝ) : e.symm (u,v) = ![u,v] := by
    ext i; fin_cases i <;> rfl
  have hg (u v : ℝ) : G (u,v) = (Set.Ico (0:ℝ) 1).indicator
      (fun u => (Set.Ioo u 1).indicator (fun v => F ![u,v]) v) u := by
    dsimp only [G]
    rw [heq, orderedSite_indicator_two]
  calc
    _ = ∫ p, G p ∂volume.prod volume := by convert! hcomp.symm using 1
    _ = ∫ u, ∫ v, G (u,v) := integral_prod G hG
    _ = ∫ u in Set.Ico (0:ℝ) 1, ∫ v in Set.Ioo u 1, F ![u,v] := by
      rw [← integral_indicator measurableSet_Ico]
      apply integral_congr_ae
      filter_upwards [] with u
      by_cases hu : u ∈ Set.Ico (0:ℝ) 1
      · simp only [hg, Set.indicator_of_mem hu]
        exact integral_indicator measurableSet_Ioo
      · simp only [hg, Set.indicator_of_notMem hu, integral_zero]
    _ = _ := by
      rw [intervalIntegral.integral_of_le (by norm_num : (0:ℝ) ≤ 1), ← integral_Ico_eq_integral_Ioc]
      apply setIntegral_congr_fun measurableSet_Ico
      intro u hu
      dsimp only
      rw [← integral_Ioc_eq_integral_Ioo, intervalIntegral.integral_of_le hu.2.le]

theorem orderedSite_integral_two_gaps (F : (Fin 2 → ℝ) → ℝ) (hF : Continuous F) :
    (∫ x in orderedSiteDomain 2, F x) =
      ∫ u in (0:ℝ)..1, ∫ v in (0:ℝ)..(1-u), F ![u,u+v] := by
  rw [orderedSite_integral_two F hF]
  apply intervalIntegral.integral_congr
  intro u hu
  have h := intervalIntegral.integral_comp_add_left (fun v => F ![u,v]) (a := (0:ℝ)) (b := 1-u) u
  simpa only [add_zero, add_sub_cancel] using h.symm

end Spin.Structured.Placement


