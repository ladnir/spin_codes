import SpinCodes.Structured.ConcretePlacementFibers

open Spin.Structured.Placement
attribute [local instance] Classical.propDecidable

example {R a : Nat} (S : Finset (Fin (128 * R))) (hS : S.card = a)
    (blocks : Fin a ↪ Fin R) (H : Nat) (hg : ¬ badBlockIndices blocks H)
    (f : Finset (Fin (128 * R)) → ℝ) :
    (Spin.Structured.Routing.shuffleLaw S).expect
      (fun T => if T ∈ singletonFiber blocks ∧ ¬ badBlockPlacement H T then f T else 0) =
    (Spin.Structured.Routing.shuffleLaw S).prob (fun T => T ∈ singletonFiber blocks) *
      (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
        (fun coords => f (singletonSupport blocks coords)) :=
  shuffle_good_singletonFiber_coordinates S hS blocks H hg f

#print axioms shuffle_prob_mem
#print axioms shuffle_prob_pair
#print axioms permutation_prob_pair
#print axioms shuffle_prob_close
#print axioms shuffle_prob_boundary
#print axioms shuffle_prob_badBlockPlacement
#print axioms singletonSupport_injective
#print axioms shuffle_prob_singletonFiber
#print axioms coordinate_tuple_law
#print axioms shuffle_singletonFiber_independent_coordinates
#print axioms mem_singletonFiber_iff
#print axioms badBlockPlacement_singletonSupport
#print axioms shuffle_good_singletonFiber_coordinates

