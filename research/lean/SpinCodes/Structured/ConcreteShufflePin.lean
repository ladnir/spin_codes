import SpinCodes.Structured.ConcreteShufflePermutation
import SpinCodes.Structured.ConcreteRowLaw

open Spin Spin.Structured.Routing

example {n : ℕ} (S : Finset (Fin n)) :
    (FinPMF.uniform (Equiv.Perm (Fin n))).map (fun σ => shuffleSupport σ S) =
      shuffleLaw S := uniform_permutation_shuffleLaw S

example {n : ℕ} (hn : 0 < n) (S : Finset (Fin n)) :
    Dominates (shuffleLaw S) (rowBernoulli S)
      (if S = ∅ then 1 else (n : ℝ) + 1) := shuffleLaw_dominates_rowBernoulli hn S

example {n : ℕ} (P : FinPMF (Finset (Fin n))) (k : ℕ) :
    (permutedLaw P).prob (fun T => T.card = k) = P.prob (fun T => T.card = k) := by
  rw [permutedLaw_eq_shuffledLaw, shuffledLaw_prob_card]

#print axioms uniform_permutation_shuffleLaw
#print axioms permutedLaw_eq_shuffledLaw
#print axioms shuffleLaw_dominates_rowBernoulli
#print axioms shuffleLaw_prob_card
#print axioms shuffledLaw_fiber_uniform
#print axioms shuffledLaw_prob_card
