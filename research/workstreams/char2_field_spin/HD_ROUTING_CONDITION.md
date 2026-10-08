# Cap rare routing failures before counting their codewords

The [dense-class obstruction](HD_COMPONENT_DENSE.md) shows why the raw
first moment cannot certify the current high-distance candidate. A rare
route creates a large low-weight subspace. Counting every line in that
subspace overwhelms the route's small probability. Exact MDS source counts
do not remove this issue.

This note bounds a family of these routing events directly. It does not
change the construction. At K20, the general gate below has probability
less than 2^-155.23606481. It covers the exhibited obstruction, but it does
not bound all causes of low distance.

## Components and the event to exclude

Use R=32768 source rows, b=64 regions, and N=bR. Each outer row is an
MDS[64,32] code over GF(65536). Put rho=32. The prefix and both post states
carry across region boundaries. The route independently permutes the R
positions within each region, uniformly.

Fix a set J of a source rows. Requiring prefix differences outside J to
vanish leaves d=ba free component values. Each component has a nonempty
interval as its support in the global prefix stream. These supports are
disjoint, and a component may cross a region boundary. There are exactly
a component-start positions in each region, for every fixed row placement.
Nonzero prefix weights do not change this component space.

For an output-position set W, let X_J(W) count components whose entire
support is routed into W. Their values span an X_J(W)-dimensional input
subspace supported on W after routing.

Suppose W is a union of g nonempty intervals. Imposing two zero post-state
coordinates at the end of each interval costs at most 2g dimensions.
The actual outer checks on J cost at most rho*a more dimensions. Thus

    X_J(W) >= rho*a + 2g + 1

is sufficient to construct a nonzero actual codeword supported on W, for
every coefficient setup. This is a sufficient dimension argument, not a
characterization of all bad codewords. Dependent outer checks, cancellation
outside W, and other subspaces require separate analysis.

## A bound for each fixed J and W

Let w_c be the size of W in region c, and assume |W|<=floor(0.4N).
A fully contained component must have its start position in W. Consequently
X_J(W)<=Y, where Y counts the routed component starts in W.

There are a starts in each region. Their image is a uniform a-subset, so
Y is a sum of independent regional hypergeometric counts Y_c. Set p_c=w_c/R.
For an integer k>=0, sampling without replacement gives

    E[C(Y_c,k)] = C(a,k) * (w_c)_k/(R)_k <= C(a,k)*p_c^k,

where (n)_k denotes the falling factorial. Expanding t^Y_c=(1+(t-1))^Y_c
for t>=1 therefore yields

    E[t^Y] <= product_c (1+p_c*(t-1))^a
             <= (1+0.4*(t-1))^(ba).

The second inequality is AM-GM and uses the total weight of W. No bound on
any individual w_c is needed. This is a containment probability bound,
not a count of field-valued vectors.

At t=3/2 and b=2*rho=64, Markov's inequality gives

    Pr[X_J(W) >= rho*a+2g+1]
      <= (24/25)^(32a) * (2/3)^(2g+1).

The route is the only randomness in this event. The inequality holds for
every fixed regionwise placement, including the measured group-of-four
placement and the individual-row affine control.

## Union over row sets and interval sets

The number of subsets of [N] with exactly g runs is C(N+1,2g), by choosing
their interval endpoints. Dropping the output-weight restriction increases
this count. Define c=(24/25)^32. For any threshold A with

    r_A = (R-A)/(A+1) * c < 1,

the successive ratios of C(R,a)c^a decrease with a. Hence the probability
that the sufficient obstruction occurs for some a>=A and some W with at
most G runs and weight at most floor(0.4N) is bounded above by

    C(R,A)*c^A/(1-r_A)
      * sum_(g=1)^G C(N+1,2g)*(2/3)^(2g+1).

With A=18300 and G=64, directed Arb arithmetic gives

    log2(bad-route upper) < -155.23606481189150.

This is a uniform bound over all row sets and interval locations, not a
test of selected cuts. The exhibited high-dimensional obstruction has
a=19660 and at most 64 output runs, so it lies within this event family.

The same formula also gives directed bounds below 2^-1700.48 for
(A,G)=(19660,128), below 2^-2141.30 for (21000,256), and below 2^-9096.78
for (25000,512). These extend the stated event family, not the whole-code
distance claim. The helper accepts both A and G explicitly.

## A second gate with restricted windows

For comparison, give each region one interval of length at most 13107.
There are at most R^64 choices of their starting positions. Disjoint
component supports imply a second direct factorial-moment bound:

    E[C(X_J(W),k)] <= C(64a,k)*(0.4)^k.

For any selected components, joint containment is at most 0.4 raised to
their total support size, which is at least the number selected. This
proof does not assume independent component-containment events.

At the conservative threshold X_J(W)>=32a+128, union over all a>=17820
and these windows gives probability below 2^-119.67218549779845. The
general interval gate above has broader location/length scope but starts
at a larger row-class threshold.

## Scope and next proof step

Both gates bound a structural event once, rather than multiplying its
probability by the dimension of the bad subspace it creates. They show
that the obstruction defeating the raw first moment can be rare enough
for a 40-bit target.

They do not supply a whole-code certificate. Remaining work includes
smaller row classes, output sets with more runs, and low-weight vectors
created by cancellations rather than entirely contained components. A
routing-conditioned rank argument must cover those cases before these
exceptional-event bounds can be used in a complete proof.

Six exact/directed tests cover disjoint and cross-region components,
containment factorial moments, unequal regional occupancies, interval
counts, geometric class tails, and the K20 bounds at two precisions.
Both gates were independently audited, including a 256-bit replay of the
general interval gate.

```text
python -B research/workstreams/char2_field_spin/theory/hd_routing_condition_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_routing_condition.py
```

No implementation, coefficient distribution, or benchmark changed.
