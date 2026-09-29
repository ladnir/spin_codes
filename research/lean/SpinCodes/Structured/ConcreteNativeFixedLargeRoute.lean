import SpinCodes.Structured.ConcreteNativeFixedLargeDefs
import SpinCodes.Structured.ConcreteFixedProfileRouted

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteRoute ConcreteOuter Placement

def FixedProfileBound (Q R : ℕ) (u : FixedFugacityChoice Q) : Prop :=
  ∀ b T : ℕ, ∀ h : b*(128*R)=T*128,
    ∀ e : Fin Q ↪ Fin (128*R), ∀ rows : Fin Q → Finset (Fin b), ∀ d : ℕ,
    (ConcreteRoutedEncoder.experimentLaw (128*R) b T).prob
      (fun ω => ConcreteRoutedEncoder.weight (regionMajor h) (embedRows e rows) ω ≤ d) ≤
      (((fixedNormCost u+fixedNormSlack u)^b/min 1 (3/1600:ℝ))/
        (profileChoices rows*profileCoefficient (fixedFugacity u) (fun i => (rows i).card)))/
        Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(128*R)))^d

theorem fixed_profile_uniform {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    ∀ᶠ R : ℕ in atTop, ∀ u : FixedFugacityChoice Q, FixedProfileBound Q R u := by
  classical
  apply Filter.eventually_all.mpr
  intro u
  exact actual_profile_routed_probability (by positivity) (by norm_num)
    (fixedNormCost_pos u).le (fixedNormSlack_pos hQ u) (fixedFugacity u)
    (fixedFugacity_pos u) (hK u)

#print axioms fixed_profile_uniform
end Spin.Structured.ConcreteNativeFamily
