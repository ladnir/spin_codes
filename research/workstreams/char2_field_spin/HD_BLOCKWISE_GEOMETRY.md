# Direct-2K blockwise SPIN: finite window exclusions

The direct blockwise SPIN restart has a useful necessary geometry test.
At K=2^20 and 40% distance, outer GRS[64,32] requires at least 46 state
elements to escape the test. GRS[128,64] is not excluded with four state
elements. At 45%, the corresponding thresholds are 246--247 and 47.
GRS[256,128] requires at least six state elements at 45%.

These are thresholds for a failure-probability lower bound, not sufficient
conditions for distance. In particular, a zero lower bound does not certify
a geometry. The test contains no accumulator or other mixing pre-stage.

This note extends the exact finite argument in
[HIGH_DISTANCE_WINDOW_OBSTRUCTION.md](HIGH_DISTANCE_WINDOW_OBSTRUCTION.md).
Two refinements remove avoidable losses: a final suffix window saves one
state-exit constraint, and unequal regional windows use the packet budget
without losing up to one packet in every region.

## Construction and the dimensional obstruction

Let b be one of 64,128,256,512, and set k=b/2. Encode R=2^20/k independent
message rows with a [b,k] MDS code over a field F. Collect each coordinate
into a region of R symbols. Sample one independent uniform permutation within
each region. The outer code and all other maps are fixed for this argument.

A causal linear inner consumes w symbols and emits w symbols per packet.
It has s field elements of state, initially zero, carried across regions.
Zero input and zero state produce zero output and zero next state. Assume
the complete unpunctured encoder is injective. An invertible input-to-output
block in every inner packet suffices. The terminal state is unrestricted,
and no extra termination symbols are appended.

The output length is N=bR=2^21. All menu geometries have R divisible by
w in {16,32,64,128}; there is no message padding or output puncturing.
For b=512, an ordinary GRS outer requires at least 512 field elements.
GF(65536) accommodates every outer in this menu. The window argument does
not otherwise depend on field size or the inner's coefficient distribution.
For 128-bit records encoded as parallel field lanes with shared coefficients,
activate just one lane. The same witness then bounds nonzero record weight.

Choose packet-aligned window lengths W_c in the b regions. In the first
b-1 regions, use a prefix window. In the final region, use a suffix window.
Let T=sum_c W_c. For source row i, let t_i count how many of its coordinates
route into these windows. MDS shortening gives the exact available dimension

    S = sum_{i=1}^R max(0, t_i-(b-k)).

Impose zero inner state at the exit of each of the first b-1 windows. There
are at most C=(b-1)s homogeneous equations. If S>C, their common kernel
contains a nonzero message. Zero state and zero input then keep the output
zero through each gap. The final suffix ends at the encoder boundary and
needs no exit constraint. Its final state may be nonzero.

Thus S>C implies a nonzero output of weight at most T, for every choice
of the inner maps. This is a statement about actual codewords, not a failed
upper-bound proof. The inherited choice of b prefix windows remains valid,
but costs bs equations instead of (b-1)s.

## Exact failure lower bound

For a fixed source row, the region-membership indicators are independent
Bernoulli variables with probabilities p_c=W_c/R. Their sum t has the
Poisson-binomial law

    Pr[t=j] = [x^j] product_c ((1-p_c)+p_c*x).

Distinct source rows need not be independent. Linearity of expectation gives

    mu = E[S] = R * sum_{j=b-k+1}^b (j-b+k)*Pr[t=j].

Exactly T coordinates enter the windows. The inequality
(t-(b-k))_+ <= kt/b therefore gives S<=kT/b=T/2. Write M=T/2, which is
an integer for this packet menu. For M>C,

    Pr_route[S>C] >= max(0, (mu-C)/(M-C)).

For M<=C the helper returns zero. If T<delta*N, the displayed expression
is a lower bound on the probability that minimum distance is below delta*N.
No independence between source rows, random inner coefficients, or asymptotic
approximation enters the argument.

Let epsilon=2^-40. A target is excluded only if the exact rational lower
bound is strictly greater than epsilon. For fixed windows, the least integer
state dimension not excluded by this test is

    s_min = max(0, ceil((mu-epsilon*M)/((b-1)*(1-epsilon)))).

Every s<s_min is ruled out at the requested distance and setup margin.
Every s>=s_min is only unexcluded by this particular test.

## Packet-balanced windows

Using the same W in every region can discard almost b packets from the
strict target budget. That produces a misleading improvement at large w.
Instead, define

    P = floor((ceil(delta*N)-1)/w),
    P = b*l+h,                  0<=h<b.

Use W_c=w*(l+1) in h regions and W_c=w*l in the others. Then T=wP<delta*N,
and adding one more packet would violate the strict budget. The final
region still uses its assigned length as a suffix.

Put D=R/w. The helper computes the integer coefficients of

    (D-l+l*x)^(b-h) * (D-l-1+(l+1)*x)^h

and divides by D^b. Thus mu and every rejection comparison are exact
rationals. Displayed decimal means and probabilities are for readability.
For the same total T, balancing the regional lengths maximizes this mean:
conditioning on the other regions, the two-coordinate smoothing term is
p_i*p_j times a nonnegative second difference of (t-(b-k))_+.

Packet width affects this obstruction only through admissible window lengths.
It does not change the number of carried state elements or add expansion.
After balancing, the integer state thresholds are almost width-independent.

## The bounded K20 menu

The table gives s_min for a 40-bit setup margin. Each entry is
`40% threshold / 45% threshold`. The helper also reproduces the inherited
equal-window, bs-cost gate as a regression.

| Outer | w=16 | w=32 | w=64 | w=128 |
| --- | ---: | ---: | ---: | ---: |
| GRS[64,32] | 46 / 247 | 46 / 247 | 46 / 247 | 46 / 246 |
| GRS[128,64] | 3 / 47 | 3 / 47 | 3 / 47 | 3 / 47 |
| GRS[256,128] | 1 / 6 | 1 / 6 | 1 / 6 | 1 / 6 |
| GRS[512,256] | 1 / 1 | 1 / 1 | 1 / 1 | 1 / 1 |

Within the requested state menu {4,8,16,32,64}, the consequences hold for
all four packet widths:

| Outer | States excluded at 40% | States excluded at 45% |
| --- | --- | --- |
| GRS[64,32] | 4,8,16,32 | 4,8,16,32,64 |
| GRS[128,64] | none | 4,8,16,32 |
| GRS[256,128] | none | 4 |
| GRS[512,256] | none | none |

The exclusions are not numerically marginal. Conservative rounded lower
bounds from several examples are:

| Outer, width, state, target | Failure lower bound | Setup margin cannot exceed |
| --- | ---: | ---: |
| [64,32], 16, 32, 40% | 0.0019873 | 8.975 bits |
| [64,32], 128, 64, 45% | 0.0244971 | 5.352 bits |
| [128,64], 16, 32, 45% | 0.0038984 | 8.003 bits |
| [256,128], 16, 4, 45% | 0.0010143 | 9.946 bits |

In particular, [256,128], width128, state4 is also excluded at 45%.
The equal-window test alone missed that case because it discarded too much
of the weight budget.

The cause is outer shortening inside routed windows. Increasing b at fixed
outer rate concentrates each row's regional occupancy below the b/2
shortening threshold. Increasing s instead adds exit equations. A wider
packet by itself does neither. Changing the diagonal alphabet cannot repair
any excluded case, because the dimensional witness holds for every local
coefficient realization.

## What remains plausible and what remains unproved

For the 40% restart, GRS[128,64] with width64 and state16 or32 is a reasonable
next proof target: it escapes this necessary gate without requiring the
state64 geometry forced by GRS[64,32]. The test gives no performance estimate
and does not choose between state16 and state32.

For a later 45% target, GRS[256,128] with state8 or larger escapes this test.
GRS[128,64] instead requires state64 within the requested menu. GRS[512,256]
also escapes, but these bounds alone do not justify its larger outer code.
The next step should analyze the actual chosen inner's rank and support
behavior for all message supports, rather than treating an unexcluded
window test as a proof pass.

## Reproduction

`theory/hd_blockwise_geometry.py` evaluates all 160 requested combinations
with exact rational comparisons. Ten bounded tests cover strict packet
budgets, the old formula, the terminal refinement, exhaustive small routing
laws, Poisson-binomial coefficients, every frozen threshold, and API guards.
A small causal GF4 example has nonzero final state and zero output outside
the prefix/suffix windows, checking why the final exit constraint is absent.
All ten tests pass.

```text
python -B research/workstreams/char2_field_spin/theory/hd_blockwise_geometry.py
python -B research/workstreams/char2_field_spin/theory/hd_blockwise_geometry_tests.py
```

Only this note and the two new helpers were written. No existing research
file was changed, and no native build, remote call, benchmark, or commit
was performed.
