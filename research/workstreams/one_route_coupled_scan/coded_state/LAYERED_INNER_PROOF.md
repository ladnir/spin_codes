# Layered inner: structural checks and the first two closed classes

Keep the pair-shear construction in [the requirements review](INNER_REQUIREMENTS_REVIEW.md)
unchanged. The GF32 sandwich remains the certified 46%/40-bit control.
The layered candidate does not yet have a whole-code certificate.
This note establishes an affine counting interface and two restricted results.
It also identifies a real weakness of the singleton-diffusion argument.

## Objects and probability space

The outer is RS[256,128], with 8192 blocks. Each of the 256 regions independently
permutes those blocks uniformly. Thus K=2^20 and N=2^21. There are 4096 batches
of 512 symbols. Keep the coefficient field J=GF(2^32), and optionally extend
the sampled generator matrix to F=GF(2^128).

Let delta be a fixed element of GF256 other than 0 and 1. In batch t, pair
indices differing in bit t mod 9, orienting the bit-zero index first.
For each pair independently, sample a from J* and c from

    S = J \ {0, 1, 1+delta, delta+delta^(-1)},

then put k=c/a. All keys are fixed at setup and independent of the route.
For a state pair (x,y), the two shears compute

    p=x+a*y,             r=y+k*p.

Writing A_t for the pair maps, the encoder emits and updates

    Y_t=X_t+(A_t+delta I)Q_t,
    Q_(t+1)=A_t Q_t+X_t=delta Q_t+Y_t,       Q_0=0.

The terminal state is free. The inverse reconstructs Q from Y without keys.
Both A_t and A_t+delta I are invertible and pairwise superregular.
This note uses symbol Hamming weight, not bit weight.

## Structural checks: what diffusion does not establish

On zero input, a singleton spreads to all 512 coordinates in nine batches.
However, an input can prevent that expansion at very low local weight.
Start with Q=u*e_0, where u is nonzero. In batch t, let j=2^(t mod 9).
The pair map gives

    A_t Q=u*e_0+k_t*u*e_j.

Choose X_t=k_t*u*e_j. Then Q_(t+1)=Q and Y_t=(1+delta)u*e_0.
Only one input and one output symbol are nonzero per batch, for every setup.
Activating the state initially costs one additional singleton input.

This is not a codeword counterexample: these prescribed inputs need not be
in the routed RS outer. It proves that arbitrary-input diffusion, or a
deterministic lower bound from active-state duration alone, is unavailable.
The outer and route must exclude or charge such trajectories.

There is also a cheap local reset. Activate Q=u*e_0 in one batch, then input
u*e_0+v*e_j in the next matching, where u,v are fixed nonzero values. The reset
occurs precisely when k=v/u, with probability 1/(q-1), q=|J|. Its two output
batches have total weight two. This disproves a uniform q^(-512) reset bound
for the layered inner. An RS outer word cannot itself have this three-symbol
input support, so the reset is not a whole-code distance obstruction.

Finally, a complete nine-layer empty-input cycle is not a uniform state
refresh. Starting from u*e_0, coordinate zero remains u along its unique
all-bit-zero path. General-state cancellation also invalidates extrapolation
from singleton expansion. An exact GF4 width-four census has maximum point
mass 9/81; the singleton (1,0,0,0) instead has 27 equiprobable outputs.
The tiny-field result tests the shortcut, not a GF32 asymptotic density bound.

## Keep the fixed coordinates in the inverse count

Fix the route and a candidate J-valued output Y. Its states Q are now fixed.
For each incoming pair, the feedback f=(A+delta I)(x,y) has four cases:

| Incoming pair and branch | Free feedback coordinates | Fixed feedback coordinates | Point probability bound |
|---|---|---|---|
|x=y=0|none|both zero|1|
|x!=0, y=0|second|first=(1+delta)x|1/(q-1)|
|y!=0, p!=0|both|none|1/((q-1)|S|)|
|x,y!=0, p=0|none|(delta*x,(1+delta)*y)|1/(q-1)|

For the regular branch, f determines a and k uniquely. For the exceptional
branch, a=x/y is fixed and all |S| choices of c give the same f. When y=0,
k is uniform on J*: each k has exactly |S| allowed a values. These facts
prove the table without an independence assumption about the two keys.

Choose a regular or exceptional branch at each pair with x,y!=0. Let F_0
be the set of free scalar coordinates in the table, D=|F_0|, and e the number
of exceptional pairs. Outside F_0 the inverse input X=Y+f is fixed.
Let d_pi(F_0) be the dimension of the routed outer shortened to F_0. The
intersection with the prescribed affine coordinate space is either empty
or has q^d_pi(F_0) points. Independence of distinct pair keys therefore gives

    Pr_keys[X is in the outer and this branch occurs]
       <= q^d_pi(F_0) (q-4)^(-D-e).

Keep an empty intersection as zero when it can be detected. Replacing its
indicator by one is an upper bound, but can lose important constraints.
This is an inverse-output bound for a fixed Y, not a distance theorem by itself.

Put eta=(q/(q-4))^N. Here eta<=512/511<2, by Bernoulli's inequality.
Thus the displayed bound is at most eta*q^(d_pi(F_0)-D-e).
Exceptional pairs have less randomness but also *two fixed coordinates*.
Treating them as ordinary free coordinates with a density penalty discards
the feature that can compensate for the weaker distribution.

## First restricted certificate: nearly full state histories

For any route and any coordinate set F_0, RS shortening gives

    d_pi(F_0)=sum_i (f_i-128)_+ <= sum_i f_i/2=D/2,

where f_i is the number of selected coordinates in block i.
Count incoming state pairs of types 10, 01, and 11 by n10,n01,n11 over all
batches. Let H=n10+n01+2*n11 be the number of nonzero incoming state symbols.
For any exceptional-branch choice,

    D=n10+2*n01+2*(n11-e),
    D/2+e=(H+n01)/2.

Summing at most 2^n11 branches proves the fixed-output bound

    Pr_keys[X(Y) is in the outer]
       <= eta * 2^n11 * q^(-(H+n01)/2).

This holds uniformly over every route. There are at most 2^N*q^W outputs
of weight at most W. At W=floor(.45N)=943718, sum over all such J-valued
outputs whose incoming histories have at most 8192 zero symbols in total.
The resulting probability is strictly below

    2^(N + N/2 + 32W - 16(N-8192) + 1) = 2^(-78655).

This is a complete union over that restricted class, including messages
chosen after setup. It is not a certificate for all messages. In particular,
an extended initial zero prefix falls outside the class.

The class is defined over J. A dense F-valued state can have sparse J-components;
do not relabel this restricted statement as a dense-history theorem over F.
Once all J-valued classes are covered, equality of minimum distance under
scalar extension transfers the complete certificate to F with no lane penalty.

## Second restricted certificate: fixed small spatial support

The sustained singleton example is excluded much more strongly by the route.
In fact, fix any set S_0 of g state-coordinate indices. Suppose every output
batch is supported in S_0. Then every Q_t is supported there, and every input
batch is supported in S_0 union its partners in the current matching.
This set has at most 2g positions per batch, independent of the pair keys.

Each region has 16 batches, so at most 32g of its 8192 positions are allowed.
For any fixed outer block, its positions in distinct regions are independent.
Its chance of landing in the allowed set in at least 129 regions is at most

    binom(256,129) * (g/256)^129.

A nonzero RS word needs at least 129 coordinates. Union over the 8192 blocks
and binom(512,g) sets S_0. Padding a smaller set to size g covers all spatial
supports of size at most g. With g=23, the exact integer comparison proves

    binom(512,23)*8192*binom(256,129)*23^129 / 256^129 < 2^(-51).

This excludes every nonzero output confined to at most 23 fixed coordinate
indices, regardless of its weight, amplitudes, or pair keys. It holds directly
over F as well as J. The same coarse estimate for g=24 does not reach 40 bits.
It is a failed bound, not an exhibited bad codeword.

## Forward companion and the remaining middle

The exceptional branch can also be charged forward. Fix a route, a set U of
allowed input batches, and a fixed nonzero outer word supported in U. Before
exposing each pair's keys, its incoming state is determined. The event
x,y!=0 and p=0 has probability at most 1/(q-1). Consequently a specified
set of e exceptional pair locations costs at most (q-1)^(-e), even though
the states depend on earlier keys. If P pair locations are available, a union
bound over their choices and outer messages gives

    binom(P,e) * q^d_pi(U) * (q-1)^(-e).

This bounds at least e exceptions for inputs supported in U. It does not
charge regular updates or assume independent states. For a suffix U of at
least N/2 scalar coordinates, d_pi(U)=|U|-N/2 for every route: the suffix
contains at least 128 full regions. Shorter suffixes have dimension zero.

Neither this bound nor the coarse inverse dimension bound closes the middle.
For illustration, ignore lower-order activation and entropy terms and take
|U|=3N/4, H approximately 3N/4, and e=3N/20. Then D approximately 9N/20.
The route-uniform bound d_pi(F_0)<=D/2 permits an inverse exponent about
0.075N at 45% output weight; the forward exponent is about 0.10N.
Both are positive. These are failures of the bounds, not codeword witnesses.

The inverse bound would instead be strongly negative if the actual route
forced d_pi(F_0) to be small. This identifies the next task: retain the
regional occupancies of the free-coordinate set, together with the history
constraints and exceptional-branch count. The fixed-coordinate consistency
equations may also supply rank that the affine-fiber upper bound discards.
Do not replace one-batch feedback by a uniform refresh, or multiply separately
averaged routing and state-event probabilities.

## A rank-based alternative keeps both shears honest

The affine method can discard consistency equations on fixed coordinates.
The following interface retains them. Fix a message space V and one common
span W of all zero constraints already imposed. Keep a list of forms required
nonzero, and discard a branch if any such form lies in W. For tau>=1, use
the potential 1[branch retained]*tau^(-dim W).

Expose a first. Its multiplier is the incoming form y. After this exposure,
the intermediate form p=x+a*y is known. Conditional on a, k is uniform on
S/a, a set of size |S|; its multiplier is p. Thus each shear admits the usual
rank-or-root argument, but the second key must use its conditional set.

Specifically, impose u+alpha*z=0 with z required nonzero modulo W. If u is
outside W+J*z, the equation increases rank by one for every alpha. Otherwise,
at most one alpha avoids forcing z into W and killing the branch. For a key
uniform on L distinct values, expected potential contracts by at most
1/min(tau,L). This works for L=q-1 in the first shear and L=|S| in the second.
If the multiplier is zero, impose the resulting deterministic equation instead.

Several zero requests using the same key do not give several root charges.
Subtract their equations first and charge the actual deterministic rank.
For example, zero output and zero next state in the first coordinate imply
delta*x=0. In the second coordinate they imply delta*y=0. Keep all earlier
nonzero requirements when imposing these differences. Branch requirements
are not conditioning events that renormalize the setup distribution.

This gives a valid sequential interface for the actual correlated pair keys.
It does not yet supply a tractable continuation count for all state histories.
Use it when fixed-coordinate outer consistency is the loss; use the affine
method when the regional free-set dimension is already sufficiently small.

## Validation and next test

test_layered_inner.py checks the affine bound on 9728 GF16 affine fibers,
exhausts 255 GF4 inputs and 81 width-four cycle setups per input, verifies the
sparse reset and sustained path, and checks both displayed integer receipts.
The g=24/40-bit negative control is rejected. No native timing is assigned.

The subsequent [matching-graph argument](LAYERED_ROUTE_GRAPH.md) carries out
the next routing step. It bounds all paired free sets simultaneously, closes
the illustrative intermediate class at45%, and covers every exception count
across all regional suffixes at a lower threshold. Its warm-up extension
reaches about40.18%/52 bits for a restricted GF32 family. Long sparse histories
remain open; that result is not a whole-code certificate either.
