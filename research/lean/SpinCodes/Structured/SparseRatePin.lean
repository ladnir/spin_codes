import SpinCodes.Structured.SparseRateOuter

open Spin.Structured.SparseRate

example : leadingConstant < -(1340384:ℝ)/100000000 := leadingConstant_lt

example {L Q b R d : ℕ} (hQ : 4096 ≤ Q) (hQL : Q ≤ L)
    (hb : 1000 ≤ b) {α ε : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1/10000)
    (hα : (L:ℝ)*α = Q) (hrounds : 128*R = L*b)
    (hd : (d:ℝ) ≤ (11/100)*(L:ℝ)*b)
    (hs : Real.log (L:ℝ) ≤ (4*Real.log 2/39)*(b:ℝ)) (hε : ε ≤ 1/500) :
    finiteBound L Q b R d α ε ≤ Real.exp (-(3/500)*(Q:ℝ)*b) :=
  finiteBound_le hQ hQL hb hα0 hα1 hα hrounds hd hs hε

#print axioms leadingConstant_lt
#print axioms sparse_kl_le
#print axioms sparse_neg_log_le
#print axioms log_conditioning_le
#print axioms remainder_le
#print axioms exponent_le
#print axioms finiteBound_le
#print axioms countedMassBound_le
#print axioms native_countedMassBound_eventually
#print axioms outer_occupation_sparse_rate
