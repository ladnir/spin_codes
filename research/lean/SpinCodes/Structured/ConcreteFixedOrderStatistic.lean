import SpinCodes.Structured.ConcreteFixedOrderStatisticSplit

noncomputable section
namespace Spin.Structured.Placement
open MeasureTheory Set

theorem orderedSite_orderStatistic_fiber (a b : ℕ) (F : ℝ→ℝ) (y : ℝ)
    (hy : y∈Ioo (0:ℝ) 1) :
    (∫ z : Fin (a+b)→ℝ,
      (orderedSiteDomain (a+b+1)).indicator
        (fun x => F (x ⟨a,by omega⟩))
        ((⟨a,by omega⟩ : Fin (a+b+1)).insertNth y z)) =
    (y^a/(a.factorial:ℝ))*((1-y)^b/(b.factorial:ℝ))*F y := by
  let f : (Fin (a+b)→ℝ)→ℝ := fun z =>
    (orderedSiteDomain (a+b+1)).indicator (fun x => F (x ⟨a,by omega⟩))
      ((⟨a,by omega⟩ : Fin (a+b+1)).insertNth y z)
  have hc := (appendMeasurableEquiv_preserving a b).integral_comp' f
  change (∫ z, f z)=_
  rw [←hc]
  have hne : ∀ᵐ v : Fin b→ℝ, ∀j,y≠v j := by
    rw [ae_all_iff]
    intro j
    filter_upwards [Measure.ae_eval_ne (fun _ : Fin b => (volume:Measure ℝ)) j y] with v hv
    exact hv.symm
  have hp : ∀ᵐ p : (Fin a→ℝ)×(Fin b→ℝ), ∀j,y≠p.2 j :=
    (Measure.quasiMeasurePreserving_snd (μ := (volume : Measure (Fin a→ℝ)))
      (ν := (volume : Measure (Fin b→ℝ)))).ae hne
  have he : (fun p => f (appendMeasurableEquiv a b p)) =ᵐ[volume]
      (fun p : (Fin a→ℝ)×(Fin b→ℝ) =>
        (orderedBetween a 0 y).indicator (fun _ => (1:ℝ)) p.1 *
        (orderedBetween b y 1).indicator (fun _ => (1:ℝ)) p.2 * F y) := by
    filter_upwards [hp] with p hp
    rcases p with ⟨u,v⟩
    change ∀j,y≠v j at hp
    rw [appendMeasurableEquiv_apply]
    have hm := orderStatisticJoin_mem y u v
    simp only [hy.1.le,hy.2,mem_Ico,true_and,and_true,hp,orderStatisticJoin] at hm
    change (orderedSiteDomain (a+b+1)).indicator (fun x => F (x ⟨a,by omega⟩))
      (orderStatisticJoin y u v)=_
    by_cases hu : u∈orderedBetween a 0 y <;> by_cases hv : v∈orderedBetween b y 1 <;>
      simp [Set.indicator,hm,hu,hv,orderStatisticJoin,hp]
  rw [integral_congr_ae he,integral_mul_const]
  rw [show (volume : Measure ((Fin a→ℝ)×(Fin b→ℝ)))=volume.prod volume from rfl,
    integral_prod_mul,integral_indicator (orderedBetween_measurable a 0 y),
    integral_indicator (orderedBetween_measurable b y 1),
    orderedBetween_integral_one a hy.1,orderedBetween_integral_one b hy.2]
  simp

theorem orderedSite_orderStatistic_integral (a b : ℕ) (F : ℝ→ℝ) (hF : Continuous F) :
    (∫ x in orderedSiteDomain (a+b+1), F (x ⟨a,by omega⟩)) =
      (1/((a.factorial:ℝ)*b.factorial)) *
        (∫ y in (0:ℝ)..1, y^a*(1-y)^b*F y) := by
  let k : Fin (a+b+1) := ⟨a,by omega⟩
  let f : (Fin (a+b+1)→ℝ)→ℝ :=
    (orderedSiteDomain (a+b+1)).indicator (fun x => F (x k))
  let e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (a+b+1) => ℝ) k
  let G : ℝ×(Fin (a+b)→ℝ)→ℝ := fun p => f (e.symm p)
  have hFD : IntegrableOn (fun x : Fin (a+b+1)→ℝ => F (x k)) (orderedSiteDomain (a+b+1)) := by
    apply (hF.comp (continuous_apply k)).integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun i => (hx.1 i).1,fun i => (hx.1 i).2.le⟩
  have hG : Integrable G (volume.prod volume) :=
    (volume_preserving_piFinSuccAbove (fun _ : Fin (a+b+1) => ℝ) k).symm.integrable_comp_of_integrable
      (hFD.integrable_indicator (orderedSiteDomain_measurable _))
  have hcomp := (volume_preserving_piFinSuccAbove (fun _ : Fin (a+b+1) => ℝ) k).symm.integral_comp' f
  rw [←integral_indicator (orderedSiteDomain_measurable _)]
  change (∫x,f x)=_
  rw [←hcomp]
  change (∫p,G p ∂volume.prod volume)=_
  rw [integral_prod G hG]
  have he : (fun y => ∫z,G (y,z)) =ᵐ[volume]
      (Ioo (0:ℝ) 1).indicator
        (fun y => (y^a/(a.factorial:ℝ))*((1-y)^b/(b.factorial:ℝ))*F y) := by
    filter_upwards [Measure.ae_ne volume (0:ℝ)] with y hy0
    by_cases hy : y∈Ioo (0:ℝ) 1
    · rw [indicator_of_mem hy]
      exact orderedSite_orderStatistic_fiber a b F y hy
    · rw [indicator_of_notMem hy]
      apply integral_eq_zero_of_ae
      filter_upwards [] with z
      change (orderedSiteDomain (a+b+1)).indicator (fun x => F (x k))
        (k.insertNth y z)=0
      apply indicator_of_notMem
      intro hx
      have hk := hx.1 k
      simp only [Fin.insertNth_apply_same] at hk
      apply hy
      exact ⟨lt_of_le_of_ne hk.1 (Ne.symm hy0),hk.2⟩
  rw [integral_congr_ae he,integral_indicator measurableSet_Ioo,
    ←integral_Ioc_eq_integral_Ioo,
    ←intervalIntegral.integral_of_le (by norm_num : (0:ℝ)≤1),
    ←intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_congr
  intro y hy
  ring

#print axioms orderedSite_orderStatistic_integral
end Spin.Structured.Placement
