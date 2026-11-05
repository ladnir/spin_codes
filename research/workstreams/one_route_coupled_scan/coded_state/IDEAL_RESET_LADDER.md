# Separate state size, reset algebra, and outer concentration

Started 2026-10-22. The target remains K=2^20, rate one half, symbol
distance greater than 49.1% over F=GF(2^128), and setup failure below
2^-40. This is proof-model work, not a native implementation or a timing
claim. Preserve the [window certificate](../NEAR_HALF_WINDOW.md) and every
previous partial proof.

Follow-up: [cheaper reset coordinates](CHEAP_RESET_COORDINATES.md) replace
the dense H_t application by structured transforms, then normalize the
state to remove products on common transitions. Both changes preserve the
two certificates below with uniform feedback. Exact algebra tests and
operation counts are available; native speedups have not been measured.

## What survives the 256-state structural screen

The proposed bridge has batches and state of 256 symbols:

    Y_t = X_t + T_t Q_(t-1),
    Q_t = Q_(t-1) + Y_t,       Q_0=0.

The T_t are independent uniform 256-by-256 matrices over F, sampled once
and fixed. The outer is still RS[256,128], with 8192 blocks. Region j
independently permutes coordinate j of those blocks. The terminal state
is free.

The earlier deterministic obstruction uses a union U of c batch-aligned
intervals: if its shortened outer dimension D(U) exceeds the number of
reset equations, a nonzero codeword is supported inside U. This argument
excluded states 64 and 128 regardless of feedback randomness.

For state 256, the first expectation-based test no longer excludes the
target. Every candidate support of at most W=1029701 symbols has at most
4022 complete batches. There are 32 batches per region. At most 2c regions
are partial. Convex equalization bounds the expected shortened dimension
by maximizing, over the possible full-region count f,

    8192 E[(f+Bin(r,p)-128)_+],
    r=min(2c,256), p=(4022/32-f)/r.

Using the largest allowed cardinality is conservative. Pad fewer partial
regions with deterministic zero/one indicators before equalization. An
exact rational calculation covers c=1,...,128; the c=128 envelope and an
increasing reset budget cover larger c. Even allowing the last reset to
be free, the maximum expected dimension per required reset is less than
184.224, attained by this envelope at c=26. It is below 256.

This result only says the expectation lower-bound method does not reject
state 256. It does NOT bound the probability of a large dimension, cover
setup-adaptively chosen interval unions, or certify minimum distance.
No structural counterexample for this 256-state ensemble was established.

## The running sum still has weak reset constraints

Two equal nonzero vectors emitted in different batches cancel the state.
They can be arbitrarily far apart. A single lane can therefore create a
weight-two closed output segment. This is an inner trajectory, not a claim
that its inverse belongs to the outer code.

The distinction matters for counting. For a length-L segment and one lane,
let P_L(z) count zero-sum outputs by symbol weight, with q=|F|. For 0<=z<=1,

    P_L(z) = ((1+(q-1)z)^L + (q-1)(1-z)^L)/q,
    P_L(z) <= (1+(q-1)z^2)(1+(q-1)z)^(L-2), L>=2.

The equality follows by counting sums with the additive characters of F.
For the inequality, put S=1+(q-1)z and use (1-z)^(L-2)<=S^(L-2).
For c episodes occupying v batches, leave the last episode unconstrained.
The number of outputs of weight at most W is at most

    inf_(0<z<=1) z^-W
      (1+(q-1)z^2)^[256(c-1)]
      (1+(q-1)z)^[256(v-2(c-1))].

Each closed episode has at least two batches. Its lane checks act on
disjoint variables from the other episodes. Dropping the nonzero-state
conditions only enlarges this bound. This improves the earlier blanket
support factor, but selected complete-class bounds remain positive.

For example, the FLOAT proposal at v4096,c64 remains about +1090 log2
units after optimizing the output marker and testing rational route
markers. This is a failed upper bound, not a codeword. The outer-shortening
moment and interval-placement count remain substantial, even after giving
the inner full feedback randomness.

## An ideal reset control, with an explicit construction

Keep N=2^21 and the outer/route unless a profile explicitly changes them.
Choose batch size b and state size s=b. Fix N distinct points xi_0,...,
xi_(N-1) in the subfield GF(2^32) of F. Define the b-by-b matrix

    H_t[r,j] = xi_(bt+j)^r,       0<=r,j<b.

Every H_t is invertible. Sample independent uniform b-by-b T_t over F,
independently of the route, once during setup. Encode with

    Y_t = X_t + T_t Q_t,
    Q_(t+1) = Q_t + H_t Y_t,     Q_0=0, 0<=t<N/b.

Emit all N symbols; the terminal state is free.
The state update is deterministic and linear. For a fixed output, states
are determined without knowing any T_t, so the inner is invertible.

For a segment starting at zero, its ending state is the sum of globally
distinct Vandermonde columns times its outputs. Any w selected columns
have rank min(b,w). Therefore:

- Every nonzero closed segment has at least b+1 nonzero output symbols.
- For each valid support, a reset removes exactly b amplitude dimensions.
- For each fixed input and nonzero entering state, a prescribed reset
  has probability q^-b, because feedback is uniform and H_t is invertible.

These are precisely the inverse and forward interfaces of the retained
[persistent-MDS proof](../PERSISTENT_MDS.md), Sections 3--6. That proof's
earlier multiplicative schedule used b<=255. We do not extend that schedule
past its domain. Globally distinct Vandermonde columns give a new explicit
schedule for b=256 or 512; the counting theorem needs only the interfaces
above. No efficient implementation of these H_t matrices is assumed.

For c episodes and v total batches, put u=bv and a=b(v-c). The inverse
input is uniform on an affine coordinate space of dimension a. Its outer
membership probability is at most q^[D(A)-a]. Valid closed segments remove
b(c-1) output dimensions. The complementary forward bound counts the
projective outer subspace on U and charges q^[-b(c-1)]. The regional
moment and representation count are unchanged. Use the better complete
bound for each class, then sum classes; never multiply independently
averaged forward and inverse bounds.

## Closed control: state 512, original outer and route

The profile b=s=512 with RS[256,128] closes the complete >49.1%/40-bit
target. The probability is over the route and inner setup, fixed for all
messages; messages may be selected after setup. This is not a certificate
for the 256-state running sum or the old structured inner.

Separate directed 80- and 120-digit replays passed. There are 535 cover
rectangles and 4194712 feasible (v,c) pairs. Every class contribution,
including its interval placements, is at most 2^-63. The conservative
rectangular count n*cmax is below 2^23, yielding failure below 2^-40.
The weakest directed rectangle is below -63.33389521022477 log2 units.

An independent integer adjacency check covers every c slice. Removing,
duplicating, or moving a cover boundary rejects. The cover SHA-256 is
`0e1d475d76e4a686c70153e51e01ec391c24b83174d5d8d39cc0635c8fe604ef`.
The two complete runs took about 19.2 and 20.6 seconds; these are verifier
times, not encoding times. A 34-marker proposal exhausted its 45-second
cap and made no claim; the retained seven-marker cover suffices.

## Why the earlier uniform-inner proof does not settle state 256

There are three different meanings of uniform here:

1. A global uniform invertible inner sends each fixed nonzero input to a
   uniform nonzero N-vector. Ordinary random-code counting suffices, but
   that is not a low-work recursive encoder.
2. The certified window has uniform feedback on each nonzero state and
   stores the last 1024 outputs. Returning to zero forces a 1024-symbol
   zero tail. Its >49.1%/40-bit proof already works with the current outer
   and route.
3. The running-sum bridge also has uniform feedback on every nonzero
   state, but stores only a 256-symbol linear summary. Returning to zero
   can instead happen by cancellation. Uniform feedback does not erase
   this distinction.

For the first control, a direct union bound is at most
`2 * 2^N * q^(K+W-N) = 2^-318847`. The factor two bounds
`q^N/(q^N-1)`. Thus the alphabet, rate, and target distance are not the
problem; retaining enough of that behavior with local, low-work operations
is the problem.

The original outer and route are therefore not intrinsically incompatible
with the target: both the window and the new 512-state control certify it.
They do create shortened outer subspaces concentrated on fragmented
interval unions. The inner must impose enough independent constraints
to eliminate those subspaces. The 64/128 failures were real joint
outer/route/memory obstructions. At 256, the current failures are still
failures of the bound, and both reset algebra and route counting matter.

The complementary test keeps state 256 but doubles the
outer to RS[512,256], over GF(2^16) extended to F. It has 4096 blocks and
512 independently permuted regions. That tests whether greater outer
redundancy per block can substitute for inner memory. It changes the
regional geometry, not the number of global routing stages. Any successful
proof or timing must be labeled separately from the original outer.

Both 80- and 120-digit complete directed replays pass for this MDS-reset variant:
517 rectangles cover 16773302 feasible pairs, with a worst log2 upper bound
below -75.60640837817. The per-class acceptance threshold is 2^-65, and
the conservative n*cmax count is below 2^25. Thus this profile also gives
>49.1% distance with failure below 2^-40. Its cover SHA-256 is
`5171c28f0b54ef5bbd1cf4307956c46cb3bfdf1aeac26b3b86e2198cda2663c2`.
Both replays have the same cover hash and verify every rectangle. The
80-digit run took65.45 seconds and the120-digit run113.15 seconds, including
proposal and replay. These are verifier times under varying machine load,
not encoding times. Interrupted earlier attempts are not passing receipts.

The original regional moment theorem permits arbitrary even B and any
fixed [B,B/2] MDS outer. Its proof does not depend on GF256 coefficients.
The larger outer uses enough distinct points in GF(2^16); extending to F
preserves each Vandermonde minor. The ideal reset columns remain in GF(2^32).

For comparison, keeping the running-sum update with this larger outer gives
a selected FLOAT bound near -5.82 at v4096,c64. That is neither the required
per-class margin nor a full cover. Stronger reset algebra still makes a
material difference in the present proof. Do not label the larger-outer
MDS certificate as a certificate for the running sum.

## Next reduction, not another blind state-size sweep

Keep both ideal MDS controls. The original state256/running-sum/RS256
combination remains open, not rejected. Its two explicit proof losses are
the weight-two reset paths and maximizing a route moment separately from
counting all interval placements. Closing those losses might avoid either
construction change; no such improvement is claimed here.

First seek a cheap realization or weakening of the MDS reset interface,
with full feedback randomness retained as a control. Compare the same
weakened interface on both certified outer/state profiles. This isolates
whether the extra outer work can buy cheaper inner work. Only then reduce
feedback randomness, and require a proof of that actual distribution.
Do not infer speed from state size, matrix names, or these certificates.

## Reproduce, serially

From the repository root:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/uniform256_probe.py
python -B research/workstreams/one_route_coupled_scan/coded_state/test_uniform256.py
python -B research/workstreams/one_route_coupled_scan/coded_state/ideal_reset_certificate.py --precision 80 --check-controls --seconds 50
python -B research/workstreams/one_route_coupled_scan/coded_state/ideal_reset_certificate.py --precision 120 --check-controls --seconds 50
python -B research/workstreams/one_route_coupled_scan/coded_state/ideal_reset_certificate.py --outer 512 --batch 256 --precision 80 --check-controls --seconds 50
python -B research/workstreams/one_route_coupled_scan/coded_state/ideal_reset_certificate.py --outer 512 --batch 256 --precision 120 --check-controls --seconds 120
```

The first command separates an exact expectation screen from FLOAT-only
proposal scores. Only the directed commands accept whole-code bounds.
The wrapper imports frozen helpers without editing them. Runs retain only
bounded in-memory coefficient/cover arrays and stdout, never generated
reports or raw data. No library or paper changes, native benchmarks,
commits, or promotions are part of this checkpoint.
Exact checks pass40 parity-enumerator inequalities,1883 Vandermonde rank
fixtures,65536 small-field output assignments,15 exhaustive full-feedback
laws,64 stream inverses,16 partition fixtures, and48 corrupt-cover controls.
The earlier random-feedback/compact-state regression also passes. An
immediate-timeout control exits without a certificate.
Time caps apply separately to proposal and replay. The initial combined
50-second cap stopped the larger-outer run without a claim; separating
the two bounded phases allowed its complete directed replay.
Higher-precision attempts under 50/55-second replay caps also stopped
without claims. The 120-digit reproduction allows up to 120 seconds per
phase. Tool calls yield during execution; all numerical jobs remain serial.
