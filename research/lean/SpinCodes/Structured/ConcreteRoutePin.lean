import SpinCodes.Structured.ConcreteRoutePermutation

namespace Spin.Structured.ConcreteRoute

example {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    (seedLaw L b).map (routeEval rows) =
      (((piPMF (fun i => Routing.shuffleLaw (rows i))).map transpose).bind
        (fun regions => piPMF (fun j => Routing.shuffleLaw (regions j)))) :=
  permutation_law_eq rows

example {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < L * b)
    (F : (Fin b → Finset (Fin L)) → ℝ) (hF : ∀ x, 0 ≤ F x) :
    (seedLaw L b).expect (fun seed => F (routeEval rows seed)) ≤
      (((b : ℝ) + 1) ^ (Finset.univ.filter (fun i => rows i ≠ ∅)).card *
        ((L : ℝ) + 1) ^ b) *
      (piPMF (fun _ : Fin b => poissonBinom
        (fun _ : Fin L => ((∑ i, (rows i).card : ℕ) : ℝ) / ((L : ℝ) * b))
        (fun _ => density_nonneg rows)
        (fun _ => density_le_one hL hb rows))).expect F :=
  permutation_expect_le_of_weight hL hb rows hw0 hw1 F hF

#print axioms Spin.Structured.ConcreteRoute.transpose_bernoulli
#print axioms Spin.Structured.ConcreteRoute.law_dominates
#print axioms Spin.Structured.ConcreteRoute.permutation_law_eq
#print axioms Spin.Structured.ConcreteRoute.permutation_dominates
#print axioms Spin.Structured.ConcreteRoute.permutation_expect_le
#print axioms Spin.Structured.ConcreteRoute.permutation_expect_le_of_weight
#print axioms Spin.Structured.ConcreteRoute.routeEval_totalWeight

end Spin.Structured.ConcreteRoute
