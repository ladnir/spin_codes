import SpinCodes.Structured.ConcreteNativeFixedLargeSpectrum
import SpinCodes.Structured.ConcreteNativeFixedLargeRoute

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteFixedNumeric Placement

def normalizedFixedProfile {Q k : ℕ} (seed : Seed k) (rows : Fin Q → Finset (Fin (k*24)))
    (u : FixedFugacityChoice Q) : ℝ :=
  ∏ i, ((spectrum seed (rows i).card:ℝ)/((k*24).choose (rows i).card:ℝ))*
    ((Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i))^(k*24)/(fixedFugacity u i)^(rows i).card)

theorem normalizedFixedProfile_bound {Q k : ℕ} (hk : 0 < k) (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k*24) W seed)
    (rows : Fin Q → Finset (Fin (k*24))) (u : FixedFugacityChoice Q)
    (hw : ∀ i, (rows i).card∈W)
    (hr : ∀ i, ((rows i).card:ℝ)/(k*24:ℕ)∈Set.Icc (13/125) (112/125))
    (hu : ∀ i, ConcreteFixedNumeric.rate (((rows i).card:ℝ)/(k*24:ℕ)) (fixedFugacity u i)≤-(8679/10000000)) :
    normalizedFixedProfile seed rows u ≤
      fixedShellCost (k*24)^Q*Real.exp (-((k*24:ℕ):ℝ)*(Q:ℝ)*fixedRowCharge) := by
  unfold normalizedFixedProfile
  have h := prod_le_prod₀ (fun i (_ : i∈(univ:Finset (Fin Q))) => by
    have := (fixedFugacity_pos u i).le
    have := (fixedMv_pos u i).le
    positivity) (fun i (_ : i∈(univ:Finset (Fin Q))) => weighted_shell_bound hk seed W hg (hw i) (hr i) u i (hu i))
  convert h using 1
  simp only [prod_const,card_univ,Fintype.card_fin,mul_pow,←Real.exp_nat_mul]
  congr 2
  ring

theorem fixed_profile_cancellation {Q k : ℕ} (seed : Seed k)
    (rows : Fin Q → Finset (Fin (k*24))) (u : FixedFugacityChoice Q) :
    (∏ i, (spectrum seed (rows i).card:ℝ))*
      (((fixedNormCost u+fixedNormSlack u)^(k*24)/min 1 (3/1600:ℝ))/
        (profileChoices rows*profileCoefficient (fixedFugacity u) (fun i => (rows i).card))) =
    (1600/3)*normalizedFixedProfile seed rows u*
      Real.exp (((k*24:ℕ):ℝ)*(Q:ℝ)*(gThree+1/10000)) := by
  have hC := (profileChoices_pos rows).ne'
  have hU := (profileCoefficient_pos (fixedFugacity_pos u) (fun i => (rows i).card)).ne'
  have he : Real.exp (((k*24:ℕ):ℝ)*(Q:ℝ)*(gThree+1/10000))=
      Real.exp ((Q:ℝ)*gThree)^(k*24)*Real.exp ((Q:ℝ)/10000)^(k*24) := by
    rw [←Real.exp_nat_mul,←Real.exp_nat_mul,←Real.exp_add]
    congr 1
    ring
  rw [he,fixedNormCost_add_slack]
  norm_num only [min_eq_right (by norm_num : (3/1600:ℝ)≤1)]
  unfold fixedNormCost normalizedFixedProfile
  simp only [prod_mul_distrib,prod_div_distrib,mul_pow,prod_pow]
  unfold profileChoices profileCoefficient at *
  field_simp
  <;> ring

#print axioms normalizedFixedProfile_bound
#print axioms fixed_profile_cancellation
end Spin.Structured.ConcreteNativeFamily
