import SpinCodes.Structured.ConcreteNativeFixedLargeProduct

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteFixedNumeric Placement

def fixedProfileCountBound (Q b L d : ℕ) : ℝ := (1600/3)*fixedShellCost b^Q*
  Real.exp (-(b:ℝ)*(Q:ℝ)*fixedRowCharge+(b:ℝ)*(Q:ℝ)*(gThree+1/10000)+
    (((Q:ℝ)*(133/125)/(262144/524287))/(L:ℝ))*(d:ℝ))

theorem fixed_counted_profile_bound {Q k L d : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ) (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (rows : Fin Q → Finset (Fin (k*24))) (u : FixedFugacityChoice Q)
    (hw : ∀ i, (rows i).card∈W)
    (hr : ∀ i, ((rows i).card:ℝ)/(k*24:ℕ)∈Set.Icc (13/125) (112/125))
    (hu : ∀ i, ConcreteFixedNumeric.rate (((rows i).card:ℝ)/(k*24:ℕ)) (fixedFugacity u i)≤-(8679/10000000))
    {P : ℝ} (hP : P ≤
      (((fixedNormCost u+fixedNormSlack u)^(k*24)/min 1 (3/1600:ℝ))/
        (profileChoices rows*profileCoefficient (fixedFugacity u) (fun i => (rows i).card)))/
        Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(L:ℝ)))^d) :
    (∏ i, (spectrum seed (rows i).card:ℝ))*P ≤ fixedProfileCountBound Q (k*24) L d := by
  have hs := normalizedFixedProfile_bound hk seed W hg rows u hw hr hu
  calc
    _ ≤ (∏ i, (spectrum seed (rows i).card:ℝ))*
      ((((fixedNormCost u+fixedNormSlack u)^(k*24)/min 1 (3/1600:ℝ))/
        (profileChoices rows*profileCoefficient (fixedFugacity u) (fun i => (rows i).card)))/
        Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(L:ℝ)))^d) :=
      mul_le_mul_of_nonneg_left hP (prod_nonneg (fun _ _ => Nat.cast_nonneg _))
    _ = ((1600/3)*normalizedFixedProfile seed rows u*
      Real.exp (((k*24:ℕ):ℝ)*(Q:ℝ)*(gThree+1/10000)))/
        Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(L:ℝ)))^d := by
      rw [←mul_div_assoc,fixed_profile_cancellation]
    _ ≤ ((1600/3)*(fixedShellCost (k*24)^Q*Real.exp (-((k*24:ℕ):ℝ)*(Q:ℝ)*fixedRowCharge))*
      Real.exp (((k*24:ℕ):ℝ)*(Q:ℝ)*(gThree+1/10000)))/
        Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(L:ℝ)))^d := by
      apply div_le_div_of_nonneg_right _ (by positivity)
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hs (by norm_num)) (Real.exp_pos _).le
    _ = _ := by
      unfold fixedProfileCountBound
      rw [←Real.exp_nat_mul,div_eq_mul_inv,←Real.exp_neg]
      simp only [Real.exp_add]
      ring_nf

#print axioms fixed_counted_profile_bound
end Spin.Structured.ConcreteNativeFamily
