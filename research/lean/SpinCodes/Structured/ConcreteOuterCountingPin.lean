import SpinCodes.Structured.ConcreteOuterCountingRoute

namespace Spin.Structured.ConcreteOuter
open Finset

example {k : ℕ} (seed : Seed k) (S : Finset (Fin (k * 24))) :
    (∑ x ∈ univ.filter (fun x : LocalMessage k => x ≠ 0),
      (Routing.shuffleLaw (support (encode seed x))).p S) =
        (spectrum seed S.card : ℝ) / ((k * 24).choose S.card : ℝ) := shuffledCount_eq seed S

example {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (F : (Fin (k * 24) → Finset (Fin L)) → ℝ) (hF : ∀ regions, 0 ≤ F regions) :
    (∑ x ∈ univ.filter (fun x : Fin L → LocalMessage k => ∀ i, x i ≠ 0 ↔ i ∈ S),
      (ConcreteRoute.law (rowSupports seed x)).expect F) ≤
        ((2 : ℝ)^(k * 24) * maxShellRatio seed)^S.card *
          (ConcreteMarked.regionsLaw S (k * 24)).expect F :=
  activeMessages_route_expect_le_max seed S F hF

#print axioms Spin.Structured.ConcreteOuter.shuffledCount_eq
#print axioms Spin.Structured.ConcreteOuter.shuffledCount_le_fair
#print axioms Spin.Structured.ConcreteOuter.shuffledCount_le_max_fair
#print axioms Spin.Structured.ConcreteOuter.good_spectrum_bound
#print axioms Spin.Structured.ConcreteOuter.tupleCount_eq_prod
#print axioms Spin.Structured.ConcreteOuter.tupleCount_le_fairRows
#print axioms Spin.Structured.ConcreteOuter.activeMessages_expect_le
#print axioms Spin.Structured.ConcreteOuter.activeMessages_route_expect_le
#print axioms Spin.Structured.ConcreteOuter.activeMessages_route_expect_le_max
#print axioms Spin.Structured.ConcreteOuter.good_activeMessages_route_expect_le

end Spin.Structured.ConcreteOuter
