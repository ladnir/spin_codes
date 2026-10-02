# Retaining the Random Refresh in the Inner Bound

The initial two-state bound charges every nonzero state as if it were the
worst expansion word. The new bound retains a component whose density is
bounded by a multiple of the uniform nonzero distribution. This is an
upper envelope; it does not assert that the actual state is uniform after
conditioning on previous output weights.

## One Step of the Comparison Experiment

Fix the binary expansion A and feedback C, with state space F_2^s and
step length t. Write L=2^s-1. The next input X consists of independent
four-bit packets. Each is zero with probability 1-p and otherwise uniform
over the 15 nonzero values. Inputs at different steps are independent.
These are the iid comparison inputs used by the dense proof, not the
actual routed BCH inputs.

For entering state a, the output is Y=Aa+X. The next state is T(a)+CX,
where T is the composition of r independent sampled transvections.
On a fixed nonzero state its distribution is the mixture

    alpha * point_mass(a) + beta * Uniform(F_2^s excluding zero),
    alpha=2^(-r), beta=1-alpha.

T fixes zero. This identity follows from the already verified one-update
kernel; it does not add a new random operation to the encoder. T and X are
independent of the entering state.

Fix 0<z<1. Define h_a=E[z^wt(Aa+X)], H=max_{a!=0} h_a, and
Hbar=L^(-1) sum_{a!=0} h_a. Define H2 and Hbar2 in the same way with z^2.
For zero entering state, put

    m0=E[z^wt(X)],
    q0=E[z^wt(X) 1_{CX=0}].

The existing Fourier calculation computes q0 and bounds a nonzero
feedback atom under the original input distribution. Call the latter
bound b. It also gives a common bound kappa on a feedback atom under
each distribution proportional to Pr[X=x] z^wt(Aa+x).

Consequently, for every a!=0,

    E[z^wt(Aa+X) 1_{CX=a}] <= c,
    c=min(H, sqrt(b H2), kappa H).

The average cancellation contribution has the stronger bound

    L^(-1) sum_{a!=0} E[z^wt(Aa+X) 1_{CX=a}] <= cbar,
    cbar=min(Hbar, 1/L, sqrt(Hbar2/L), kappa Hbar).

For the 1/L bound, each X satisfies CX=a for at most one nonzero a, and
z^wt(Aa+X)<=1. Cauchy--Schwarz on the joint sum over a and X gives
sqrt(Hbar2/L). The remaining two bounds follow from the total moment
and the common tilted atom bound.

## The Three Coordinates

Represent an unnormalized incoming state measure by three nonnegative
coordinates (Z,M,U). Z bounds its mass at zero. Its nonzero part is
dominated by a sum of two measures: one of total mass at most M, and one
whose density at every nonzero state is at most U/L. Thus total mass is
at most Z+M+U.

Multiplying by z^wt(Y) and averaging the next input and update preserves
such a representation, with row-vector coordinates bounded by

```text
             zero                    arbitrary       uniform-density
zero         q0                      m0-q0           0
arbitrary    alpha*c + beta*H/L       alpha*H         beta*H
uniform      alpha*cbar + beta*Hbar/L alpha*Hbar      beta*Hbar
```

To justify the arbitrary row, the lazy branch has nonzero output mass at
most alpha H. On the refresh branch, for each fixed input and entering
state, every possible next state has probability at most 1/L. Averaging
the output weight therefore gives density at most beta H/L. The zero
entry separately bounds lazy and refreshed cancellation. Counting zero
again in a nonzero upper bound is conservative.

For the uniform-density row, dominate the incoming measure by U times
the uniform nonzero distribution before applying the nonnegative
transition. Its total output moment is bounded using Hbar, and its lazy
zero contribution using cbar. The same refresh argument supplies the
pointwise density bound. It is important to use domination here, not to
assume that the actual conditioned state is uniform.

Starting from (1,0,0), the matrix power for all 16384 steps, followed by
the terminal vector (1,1,1), bounds the full iid comparison moment.
The existing outer and shuffle factors then give bounds on message
counts over sampled setup.

## Computation and Checks

`refresh_kernel.py` evaluates the matrix with outward Arb arithmetic.
It retains the multiplicity of each expansion-image histogram, so Hbar
is an actual uniform average over all nonzero states. Integer histogram
encoding accelerates preparation without changing the counted objects.
The original two-state kernel remains available as a control.

The tests enumerate every eight-bit input, every three-bit entering
state, four activity parameters, and one through three updates. They
check lazy-mass and refreshed-density bounds separately and compare
eight successive bounded steps against the exact transition matrix.
Additional tests compare histogram generation against the previous
implementation and floating proposals against outward evaluation.

This bound improves the first proof attempt but is not by itself a
distance certificate. A claim for the code additionally requires a
complete dense cover, matching sparse coverage, and their summed bound.

## Optional Nonzero-Feedback Condition

The `feedback-refresh` variant also computes nu=Pr[CX!=0] by the exact
unweighted Fourier sum. A refresh can return a nonzero state to zero only
when CX is nonzero. Its zero-entry multipliers can therefore use
min(H,sqrt(H2 nu)) and min(Hbar,sqrt(Hbar2 nu)), respectively, in place
of H and Hbar. The uniform lazy cancellation bound uses nu/L and
sqrt(Hbar2 nu/L) in place of 1/L and sqrt(Hbar2/L).

These follow from Cauchy--Schwarz on the weighted feedback event. They do
not assume independence between that event and output weight. In
particular, an empty input cannot cause a return to zero. The original
`refresh` choice remains available for reproducing the 4% certificate.

## GF-Specific Cancellation in the IID Comparison

The `gf-refresh` variant sharpens the lazy return by conditioning on the
number J of active packets in one step. Under the iid comparison input,
J has the binomial distribution with W=t/4 trials and probability p.
Conditional on J=j, positions are a uniform j-subset and values are
independent uniform nonzero packets. This is exactly the distribution in
[the fixed-occupancy calculation](OCCUPANCY_BOUND.md).

Write b_j for its nonzero feedback-atom bound and nu_j for its probability
of nonzero feedback. Let d_A be the minimum nonzero expansion weight and
put e_j=z^max(0,d_A-4j). For each fixed nonzero entering state, the
weighted lazy return is bounded by

    c_GF = sum_{j=0}^W binom(W,j) p^j (1-p)^(W-j) e_j b_j.

For the uniform-density component the corresponding bound is

    cbar_GF = L^-1 sum_{j=0}^W binom(W,j) p^j (1-p)^(W-j) e_j nu_j.

These are averages of valid conditional bounds, not an independence
approximation between output and feedback. The kernel uses exact feedback
counts through eight packets and rational atom bounds at larger j.

In the three-coordinate matrix, the arbitrary-to-zero entry can therefore
also use alpha c_GF + beta H/L. The uniform-to-zero entry can use
alpha cbar_GF + beta Hbar/L. The implementation takes each minimum with
the existing `feedback-refresh` entry and leaves all other entries unchanged.
Thus the state-domination invariant is unchanged.

The outward evaluator checks the input probabilities exactly: conditional
on activity, packet weights must have probabilities (4,6,4,1)/15. The
floating evaluator only proposes witnesses. Small-model tests enumerate
every input and compare eight consecutive transitions for one, two, and
three updates, including activity probabilities zero and one.

## Retain Input Weight in the Feedback Bound

The occupancy floor charges all j packets as weight four. The
`weighted-return` kernel instead uses the deterministic inequality

    z^wt(Aa+X) <= z^wt(Aa) z^-wt(X)       (a!=0).

The right side may exceed one; it is still a valid upper bound. Set
h=1/z. For feedback character u, let r_w(u) be the weight of its four
feedback-column pairings at packet position w. Define

    f_r = (1-p) + p ((1+h)^(4-r) (1-h)^r - 1)/15,
    F(u) = product_{w=1}^W f_{r_w(u)}.

Fourier inversion gives the unnormalized tilted feedback mass

    m(b) = E[h^wt(X) 1_{CX=b}]
         = 2^-s sum_u (-1)^(u dot b) F(u).

For every nonzero b and every real c, character orthogonality implies
m(b)<=2^-s sum_u |F(u)-c|. The implementation takes the smaller bound
from c=0 and c=(1-p)^W; call it B_h. It also computes

    V_h = sum_{b!=0} m(b) = f_0^W - 2^-s sum_u F(u).

Thus the lazy-return multipliers can use

    c_weighted = z^d_A B_h,
    cbar_weighted = min(z^d_A V_h/L, B_h E_{a!=0}[z^wt(Aa)]).

For the first uniform bound, fix X and use its unique cancelling nonzero
state a=CX. For the second, bound each state's tilted feedback mass by
B_h and then average its expansion-weight factor. Neither argument
assumes that feedback and output weight are independent.

As above, multiply these lazy bounds by alpha, add the existing refresh
upper bound, and take the minimum with the earlier matrix entry. The
expansion average is computed from the exact histogram multiplicities.
This refinement changes neither the sampled setup nor the encoder.
