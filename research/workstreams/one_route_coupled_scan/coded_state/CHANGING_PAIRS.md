# Nonzero updates between the paired observations

The [actual-support route bound](FEEDBACK_SUPPORT.md) extends to a
full-length live trajectory with **up to 32 nonzero intervening batches**.
At K=2^20 and rate 1/2 over GF(2^128), the probability that such a trajectory
has weight at most floor(.491 N) is below 2^-61. This is a restricted-class
setup-failure bound, not a whole-code distance certificate.

The [construction](ALTERNATING_MAPS.md), including the outer, route, maps,
diagonal distributions, and free terminal state, is unchanged. No encoder
timing or implementation change accompanies this result.

## Class and probability space

Retain n=16384 batches, state size s=64, batch size 128, and Q_0=0.
Require Q_1,...,Q_(n-1) to be nonzero; Q_n is free. Pair feedback
observations (2,3),(4,5),...,(n-2,n-1), giving P=8191 pairs.
An intervening batch is the first output in one such pair: y_2,y_4,...,
y_(n-2). Previously all of these outputs were zero. Here at most d=32
may be nonzero, at arbitrary positions, with arbitrary supports and values.
All odd outputs and y_n remain unrestricted apart from the total weight.

Choose d allowed exceptional positions. This covers any trajectory with
fewer exceptions by padding its allowed set. We union over binom(P,d)
sets; we do not separately sum every smaller d. Fix the nonzero compression
scalings R throughout the count. All bounds below hold for each such R.

## Fresh variables for two different states

For one exceptional pair, fix holes H_0,H_1 in opposite evaluation phases,
with sizes h_0,h_1<=63. Let u and v be the weights of the preceding odd
output and the intervening even output. Their variables are distinct.
Put r_u=min(s,u), r_v=min(s,v), and h=h_0+h_1.

With the earlier state fixed, the coefficient matrix on these fresh
variables has the block form

    [ A_i[H_0] B_U R_U                  0          ]
    [ A_(1-i)[H_1] B_U R_U   A_(1-i)[H_1] B_V R_V ].

Two valid rank lower bounds are

    r_sep = (h_0+r_u-s)_+ + (h_1+r_v-s)_+,
    r_joint = (min(s,h)+r_u-s)_+ + (h+r_v-2s)_+.
    r_pair = max(r_sep,r_joint).                         (1)

The first uses diagonal blocks in temporal order. For the second, initially
treat the incoming state Q and the even output as variables. The rank on Q
is min(s,h), since the two evaluation sets are disjoint. The full system
has at least (h+r_v-2s)_+ independent rows after eliminating Q, by the
two-state rank argument in ALTERNATING_MAPS.md. These rows survive when Q
is restricted to the affine image of B_U. The remaining Q-dependent rows
retain rank at least (min(s,h)+r_u-s)_+. This proves (1) without charging Q
as an additional random source.

Apply (1) to consecutive, disjoint pairs of fresh output blocks. Earlier
constraints do not contain these variables, so triangular minors add their
rank lower bounds. Shared states are retained; they are not resampled.

## A complete count using the simpler part of (1)

For a nonexceptional pair the same state is observed twice, so its total
hole count is at most 63. For an exceptional pair each state may have 63
holes separately. The last observation has at most 63 holes. Therefore

    H_max = 63(8192+d),       |S^c| = 128+H.

For fixed exceptional positions the number of masks of total size H is
bounded by

    L_(d,H) = [x^H] P_256(x)^(8191-d) P_128(x)^(2d+1),
    P_l(x) = sum_(j=0)^63 binom(l,j) x^j.               (2)

There are U_d=128(8193+d) potentially nonzero output coordinates. Use the
unchanged-state combined constraint for nonexceptional pairs, and r_sep
for exceptional pairs. The scalar tradeoff from FEEDBACK_SUPPORT.md applies
to each fresh block and gives total remaining dimension at most

    v_(d,H) = min(W, U_d-H-ceil(64(U_d-W)/127)),
    W = floor(.491 N) = 1029701.

For each fixed mask the number of output candidates is thus at most

    C_(d,H) = (U_d+1) binom(U_d,max(v_(d,H),U_d/2)) q^v_(d,H).

Using the actual feedback support in the regional route moment gives

    failure <= binom(8191,d) sum_(H=0)^H_max
       L_(d,H) C_(d,H)^(log_q tau) tau^(-K)
       G_bin(128+H,tau)^8192.                          (3)

Here G_bin is the exact positive binomial moment defined in
FEEDBACK_SUPPORT.md; tau may be fixed separately for each H. Both the
exception-position count and mask count remain OUTSIDE the fractional
power. Output supports and amplitudes are inside C. Thus the bound covers
witnesses selected after setup, not just a fixed input word.

The directed verifier at d=32 covers all 518113 integer H values in 74
disjoint intervals at both 80 and 120 decimal digits. Every summand,
including binom(8191,32), is at most 2^-80. Their sum is below 2^-61.
The weakest accepted interval is H=515584..515710, tau=10, with log2
upper bound approximately -81.6787. This is a conservative checkpoint,
not an assertion that 32 is the largest provable number of exceptions.

The d=34 replay fails the chosen -80-per-term acceptance rule at H=515390
(upper approximately -79.5893). That failure alone does not rule out a
40-bit bound, a different cover, or the construction's actual distance.

## Additional structure of the actual maps

Each of the two 64-point cosets in Y_1 gives an MDS restriction of A_1 B_W
for ANY support W of size w<=64, including supports crossing both input
halves. This does not assert that arbitrary minors across both output
cosets have full rank.

To prove it, let V be the degree-64 subspace polynomial taking values 0,1
on the two input halves and alpha on Z. Write g_W for the vanishing
polynomial of W. The interpolating input polynomial is
V(V+1)u/g_W with deg(u)<w. If d is its degree-<64 representative on Z, then

    d g_W = alpha(alpha+1)u + (V-alpha)v,   deg(v)<w.

On a coset V=eta in Y_1 the numerator on the right has degree<w and the
denominator g_W is nonzero. The map from supported inputs to that numerator
is injective: a zero numerator would make d vanish on all 64 coset points,
then d=0, then the input is zero because B_W has rank w. Both spaces have
dimension w, proving the claimed GRS image. Nonzero column scalings R do
not change it.

Consequently rank(A_1[H] B_W)>=min(ceil(|H|/2),w,64); for w>=64 one can
also use the exact full-state image. In our canonical pairs the first
observation uses A_1. The separate bound in (1) can therefore be improved
using this floor for its first block. It passes all 988 square-minor checks
within the two cosets of the state4 GF16 analogue, and the corresponding
multitime rank tests. The proof above, not the toy enumeration, establishes
the statement over GF(2^128). The directed result (3) does not need it.

## What the stronger pair count does and does not fix

The separate FLOAT diagnostic retains r_pair rather than replacing it by
one scalar slope. For a positive output marker z=q^(-beta), bound the local
support sum by 2^128 q^max(w-r-beta*w), or 2^256 times the analogous
two-output maximum. Raise this bound to log_q(tau), then count hole masks
with a positive polynomial marker. This respects the same ordering of
conditioning, output union, fractionalization, and mask union as (3).

At selected near-maximum H values, the stronger count still gives positive
upper scores for d=64 and larger. At d=64,H=518080 its best sampled score
is approximately +9189 log2 units. Including the within-coset rank floor
does not improve these selected optimized scores. These are failed FLOAT
upper-bound diagnostics, not counterexamples, certificates, or a proof
that all possible markers fail.

A concrete loss is visible in the coarse rank envelope. At h_0=h_1=63
and beta=64/127, its exact maximizers are (u,v)=(1,128) and (128,128), both
with scaled objective 126/127. Repeating the sparse maximizer would place
a weight-one odd update between two adjacent states with 63 holes each.
The retained one-half-update theorem forbids that: such an update permits
at most 106 holes across the two observations, not 126. The coarse pair
count missed this because those observations straddle its pair boundary.

This is an explicit relaxation artifact, not an exhibited bad codeword.
The improved within-coset floor removes some local extrema, but maximizing
independent pair envelopes still loses boundary compatibility. The next
discriminating test is an overlapping/transfer count retaining the last
hole profile and update support information. Apply the sparse-update
restriction on BOTH parities before any pairwise maximization. Do not
change the construction to repair a profile already known impossible.
Arbitrary changing trajectories, later activations, and multiple episodes
remain outside the new certificate.

## Reproduce serially

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_changing_pairs.py
python -B research/workstreams/one_route_coupled_scan/coded_state/changing_pair_certificate.py --precision 80
python -B research/workstreams/one_route_coupled_scan/coded_state/changing_pair_certificate.py --precision 120
python -B research/workstreams/one_route_coupled_scan/coded_state/joint_pair_probe.py
python -B research/workstreams/one_route_coupled_scan/coded_state/joint_pair_probe.py --phase1
```

The tests check 2048 three-pair GF16 trajectories, the within-coset minors,
local dimension bounds, small mask-polynomial domains, exact envelope
extremizers, the d=0 regression, and timeout rejection. The joint rank
floor improves 1125 local pairs over the separate floor in those fixtures.
The FLOAT probe uses NumPy and cannot accept a certificate. Jobs have
45-second limits and small bounded arrays; none writes generated data.
