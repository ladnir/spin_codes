import SpinCodes.Structured.SparseConditioningMoment

open Spin Spin.Structured.ConcreteMarked

example {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L) :
    1/(8*Real.sqrt (Q:ℝ)) ≤ markMass L Q ((Q:ℝ)/L) := markMass_mode_lower hQ hQL

example {L Q : ℕ} (hQ : 0 < Q) (hQL : Q ≤ L)
    {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (b : ℕ) :
    (1/markMass L Q p)^b ≤ (8*Real.sqrt (Q:ℝ))^b *
      Real.exp (((L:ℝ)*b)*Spin.binKL ((Q:ℝ)/L) p) :=
  markMass_inverse_pow_le hQ hQL hp0 hp1 b

#print axioms binPMF_moments
#print axioms binPMF_variance
#print axioms mode_mass_of_variance
#print axioms markMass_mode_lower
#print axioms markMass_likelihood_ratio
#print axioms markMass_lower
#print axioms markMass_inverse_pow_le
#print axioms fairExperiment_sparse_probability_exp
