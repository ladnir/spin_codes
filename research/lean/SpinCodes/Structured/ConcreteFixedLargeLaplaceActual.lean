import SpinCodes.Structured.ConcreteFixedLargeLaplaceAssembly
import SpinCodes.Structured.ConcreteFixedOrderStatistic

/-! The spacing Laplace bound for the actual ordered placement law. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Set MeasureTheory

theorem liveDuration_laplace_le {Q : ℕ} {sig tau G : ℝ}
    (hQ : 0 < Q) (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (hG : ∀ y ∈ Icc (0:ℝ) 1, -tau*y+(1+1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G)
    (S : Finset (Fin (Q+1))) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      Real.exp (-(Q:ℝ)*tau*liveDuration S x)) ≤ Real.exp ((Q:ℝ)*G)*sig^S.card := by
  apply liveDuration_laplace_le_of_orderStatistic ?_ hQ hsig0 hsig1 hG S
  intro a b F hF
  exact orderedSite_orderStatistic_integral a b F hF

theorem paper_liveDuration_laplace_le {Q : ℕ} (hQ : 3 ≤ Q) (S : Finset (Fin (Q+1))) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      Real.exp (-(Q:ℝ)*(133/125)*liveDuration S x)) ≤
      Real.exp ((Q:ℝ)*paperLargeExponent)*(127/250:ℝ)^S.card := by
  exact liveDuration_laplace_le (by omega) (by norm_num) (by norm_num)
    (paperLargeExponent_bound hQ) S

#print axioms liveDuration_laplace_le
#print axioms paper_liveDuration_laplace_le
end Spin.Structured.Placement
