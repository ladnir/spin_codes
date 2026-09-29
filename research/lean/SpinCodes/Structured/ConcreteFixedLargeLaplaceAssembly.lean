import SpinCodes.Structured.ConcreteFixedLargeLaplaceEndpoints

/-! Assembly of the spacing Laplace bound from the canonical order-statistic integral. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Set MeasureTheory
attribute [local instance] Classical.propDecidable

def OrderStatisticIntegral : Prop := ∀ (a b : ℕ) (F : ℝ → ℝ), Continuous F →
  (∫ x in orderedSiteDomain (a+b+1), F (x ⟨a,by omega⟩)) =
    (1/((a.factorial:ℝ)*(b.factorial:ℝ))) *
      (∫ y in (0:ℝ)..1, y^a*(1-y)^b*F y)

theorem liveDuration_laplace_le_of_orderStatistic (hMarginal : OrderStatisticIntegral)
    {Q : ℕ} {sig tau G : ℝ} (hQ : 0 < Q) (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (hG : ∀ y ∈ Icc (0:ℝ) 1, -tau*y+(1+1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G)
    (S : Finset (Fin (Q+1))) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      Real.exp (-(Q:ℝ)*tau*liveDuration S x)) ≤ Real.exp ((Q:ℝ)*G)*sig^S.card := by
  by_cases hS : S = ∅
  · subst S; exact liveDuration_laplace_empty_le hG
  by_cases hSu : S = Finset.univ
  · subst S; exact liveDuration_laplace_univ_le hQ hsig0 hsig1 hG
  have hc0 : 0 < S.card := card_pos.mpr (nonempty_iff_ne_empty.mpr hS)
  have hc1 : S.card < Q+1 := by simpa only [Fintype.card_fin] using (card_lt_iff_ne_univ S).mpr hSu
  obtain ⟨a,ha⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt hc0)
  obtain ⟨b,hb⟩ := Nat.exists_eq_add_of_le (show a+1 ≤ Q by omega)
  have hdim : Q = a+b+1 := by omega
  clear hb
  subst Q
  let i : Fin (a+b+1) := ⟨a,by omega⟩
  have hcard : S.card = (Finset.Iic i.castSucc).card := by
    rw [Fin.card_Iic]
    simpa only [Fin.val_castSucc,i] using ha
  let F : ℝ → ℝ := fun y => Real.exp (-((a+b+1:ℕ):ℝ)*tau*y)
  have hFc : Continuous F := by dsimp [F]; fun_prop
  have hFdiff : ContDiffOn ℝ 1 F (Icc 0 1) := by dsimp [F]; fun_prop
  have hI : (∫ x in orderedSiteDomain (a+b+1), F (liveDuration S x)) =
      (1/((a.factorial:ℝ)*(b.factorial:ℝ)))*
        (∫ y in (0:ℝ)..1, y^a*(1-y)^b*F y) := by
    rw [liveDuration_integral_card_eq S (Finset.Iic i.castSucc) hcard hFdiff]
    simp_rw [liveDuration_prefix_site]
    exact hMarginal a b F hFc
  have hN := hMarginal a b (fun _ => (1:ℝ)) continuous_const
  simp only [mul_one] at hN
  let A : ℝ := ((a+b+1).factorial:ℝ)*(1/((a.factorial:ℝ)*(b.factorial:ℝ)))
  have hNorm : A*(∫ y in (0:ℝ)..1,y^a*(1-y)^b) = 1 := by
    have h := orderedSite_normalized_one (a+b+1)
    rw [hN] at h
    simpa only [A,mul_assoc] using h
  have h := beta_laplace_normalized hsig0 hsig1 a b (show 0 ≤ A by dsimp [A]; positivity) hNorm hG
  change ((a+b+1).factorial:ℝ)*(∫ x in orderedSiteDomain (a+b+1), F (liveDuration S x)) ≤ _
  rw [hI,← mul_assoc]
  simpa only [A,F,ha,Nat.succ_eq_add_one] using h

#print axioms liveDuration_laplace_le_of_orderStatistic
end Spin.Structured.Placement
