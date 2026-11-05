# A limit of separable bounds and a delayed pairing gain

The [causal flow](RANDOM_PAIRING_FLOW.md) reduces a selected structural
diagnostic to bounding overflow loss. It does not yet bound code distance.
This continuation asks which tail arguments can exploit random re-pairing.
The encoder, route, pair arithmetic, and prior certificates are unchanged.

## A separable one-step bound loses too much

Fix the target-row count s=1974 and set p=s/8192. First use IID target
labels with probability p. The actual regional model is obtained by
conditioning each region to contain exactly s labels. All auxiliary states
below belong to the five-state flow, not the field-valued encoder state.

Consider a positive, time-independent product potential

    V(h) = product_(a=-2)^2 v(a)^h_a.

A one-step exponential envelope seeks B>0 such that, for every histogram h,

    E[exp(theta*loss + eta*(J-512p)) V(h_next) | h]
       <= B^256 V(h),

where theta>=0, eta is real, and J is the number of target labels used
in that update. This condition averages the actual uniform matching.
It does not replace the matching by an adversarial pairing.

Nevertheless, homogeneous histograms impose a strong restriction. If all
coordinates have auxiliary state a, every matched pair has equal states.
For such a pair, j target labels give

    a_next = clamp(a+2j-1,-2,2),
    ell(a,j) = 2*1[a=-2 and j=0].

Let P_j=binom(2,j)p^j(1-p)^(2-j). Setting f(a)=v(a)^2, the envelope implies

    sum_j P_j exp(theta*ell(a,j)+eta*(j-2p)) f(a_next) <= B f(a)

for every a. Thus an equal-partner five-state process restricts even a
bound that averages uniform random pairings correctly.

### An exact witness against this proof template

Choose any positive transition probabilities Q_a(j) for the equal-partner
process, and let mu be a stationary distribution. Define

    ell_Q = sum_(a,j) mu_a Q_a(j) ell(a,j),
    j_Q   = sum_(a,j) mu_a Q_a(j) j,
    D_Q   = sum_(a,j) mu_a Q_a(j) ln(Q_a(j)/P_j).

Apply Jensen's inequality to the preceding envelope, then average over mu.
Stationarity cancels the log f terms, giving

    ln B >= theta*ell_Q + eta*(j_Q-2p) - D_Q.

The helper `random_pairing_flow_tail_screen.py --verify` supplies an exact
rational stationary flow with

    j_Q = 2p,
    ell_Q >= 371901/(4095*256),
    (4095*256) D_Q < 5500 ln 2.

The witness is a rational mixture of two rational stationary flows. The
mixture weight enforces the target-label mean exactly. The verifier checks
stationarity and the loss inequality without floating-point arithmetic.
For the entropy inequality it uses

    D_Q <= sum_(a,j) mu_a (Q_a(j)-P_j)^2/P_j,
    ln 2 > 56/81 > 69/100.

The two transition tables are fixed integer constants in the helper; the
verification path does not invoke the floating-point search that found them.

For T=4095 updates, the resulting Chernoff exponent satisfies

    -theta*371901 + 256*T ln B > -5500 ln 2.

Meanwhile binom(8192,1974)>2^6519, verified by an exact integer comparison.
Consequently this envelope followed by a separate union over target-row
sets yields an expression greater than 2^1019. The initial/final product-
potential factor is at least one under the usual minimum-potential bound.
The ignored initial input batch has a centered IID label moment at least
one. Conditioning back to exact regional counts adds cost, not a saving.

This is a proof-template obstruction, not a code counterexample. It does
not exclude histogram-dependent potentials, multi-update envelopes,
time-dependent envelopes, direct treatment of remaining regional counts,
or an argument that avoids this separate union over target sets.

For context only, a floating spectral screen finds an optimistic exponent
about -2226.37 bits before the union. The row-set union costs about6519.51
bits, and conditioning256 regions costs about1689.70 bits. These diagnostic
values are not used in the exact obstruction above.

## Where re-pairing first changes the mean loss

The separable obstruction suggests retaining more than a one-step
coordinate-wise estimate. A short calculation identifies the first
re-pairing benefit in the causal flow.

Fix an incoming histogram h on n coordinates, with n even. Write
a=h_-2, b=h_-1, and c=h_0. Suppose three consecutive updates lie within
one region. Before those updates, condition on the exposed history and
on m remaining regional positions containing r target labels. Assume
m>=3n. The remaining target labels are a uniform r-subset; the upcoming
state matchings are fresh uniform matchings, independent of that subset.

Let (x)_k=x(x-1)...(x-k+1), with (x)_0=1 and (x)_k=0 for integer 0<=x<k.
Write L_0,L_1,L_2 for the losses in these three updates. Then

    E[L_0 | h,m,r] = ((m-r)_2/(m)_2) a,

    E[L_1 | h,m,r] = ((m-r)_4/(m)_4) (a+b),

    E[L_2 | h,m,r] = ((m-r)_6/(m)_6) (a+b+c)
       + (2r (m-r)_5/(m)_6) a(a+b-1)/(n-1).

These formulas use sampling without replacement, not an independent-label
approximation. They do not cover a three-update block crossing a region
boundary without a corresponding change in the label factors.

To derive the last identity, observe that a zero-target pair decrements
each auxiliary state by one, clipping at -2. A pair containing at least
one target cannot produce state -2. Therefore loss at the third update
requires zero-target pairs in the second and third updates, and state at
most -1 after the first update.

The first update can produce such a state in two ways. A zero-target pair
does so at coordinates initially in states -2,-1,0. Those trajectories
require six specified labels to be zero. A one-target pair produces a
total of two states -1 from input (-2,-2), or one from input (-2,-1).
Under a uniform matching, the expected number is

    a(a-1)/(n-1) + ab/(n-1) = a(a+b-1)/(n-1).

Those trajectories require one target among the first two specified
positions and zero targets among the next four. Their probability is
2r(m-r)_5/(m)_6. Adding the two cases proves the formula.

For comparison, if every histogram count is even, pair equal states and
keep those partners fixed. The last term then contains a instead of
a(a+b-1)/(n-1). Random re-pairing reduces the third update's mean loss by

    (2r (m-r)_5/(m)_6) a(n-a-b)/(n-1).

Thus pairing a low auxiliary state with a state outside {-2,-1} has an
explicit delayed benefit. The first two loss means do not see that benefit.
For IID labels, replace the six-position factors by (1-p)^6 and
2p(1-p)^5, respectively. Neither mean identity implies the needed tail.

## Validation and next discriminating step

The exact barrier verifier passes in0.005s. The finite-population kernel
suite passes its previous36 profile,8 full-transition, and local-flow
checks, plus5 IID delayed-gain checks and9 regional three-step mean
checks. The last checks compose the exact histogram/remaining-label
transition and compare its means with the displayed formulas.
The full suite took0.895s. Runs were serial, local, Python-B, and file-free;
local free space was22.56GiB before the numerical pass. No native encoder
benchmark, construction change, paper edit, commit, or promotion was made.

Next test a multi-update exponential bound retaining histogram information
and regional counts. The quadratic term above identifies information that
must not be discarded. Begin with the1974-row diagnostic; reject a proposed
envelope if it cannot pay that class's union cost. A successful structural
bound must still compose with field-cancellation constraints and cover
other support histories before it becomes a whole-code distance result.
