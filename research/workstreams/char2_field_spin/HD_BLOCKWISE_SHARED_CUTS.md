# Source and route cuts for the shared-transform candidate

The shared eight-point feedthrough removes the support-preservation step in the
old width-64/state-4 quiet-subcode rejection. That change is not a 40% distance
proof. A new, source-aware bound does settle one complete class: the probability
that any nonzero quiet word uses at most 256 outer rows is below
2^-2236.70950574. This bound covers every output weight in that class.

The result needs only the outer distance and the feedback-kernel distance. It
does not depend on the feedthrough translation or the distribution of mixers.
It therefore applies to the saved beta8 construction and to a beta64 variant
that preserves the stated feedback map. Larger source sets and nonquiet
histories remain open.

## Model and the quiet event

Use R=8192 independent GRS[256,128] source rows, one independent uniform row
permutation in each of the 256 coordinate regions, and packets of 64 positions.
Every region contains complete packets. The encoder has four state coordinates,
initial state zero, and free terminal state. Its packet equations are

    Y_i = F X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The fixed maps are defined in [HD_BLOCKWISE_SHARED_MAPS.md](HD_BLOCKWISE_SHARED_MAPS.md).
In particular, ker C is an MDS[64,60,5] code. All claims below concern field
symbols; the shared-lane implementation inherits any one-lane obstruction.

Call a word quiet when every emitted state Q_i is zero. This is equivalent to
CX_i=0 for the first P-1 packets, where P=32768. The final feedback syndrome is
unrestricted. A source row is active when its message, and hence its outer
codeword, is nonzero. Each active row contributes at least 129 nonzero input
symbols before the inner map.

The event bounded here is the existence of a nonzero quiet word with at most
256 active source rows. We do not condition on its output being low weight.
An upper bound on this larger event also covers all low-weight quiet words
in the same source classes.

## A necessary routed-support condition

Fix a source-row set J of cardinality a. In one region, mark the a routed
positions belonging to J, regardless of whether a particular source word is
zero at some of those coordinates. Call a packet crowded if it contains at
least five marked positions. Let Z_c count marked positions in crowded packets
of region c, and let Z=sum_c Z_c.

Every nonzero nonterminal quiet packet contains at least five actual nonzero
inputs, because ker C has minimum distance five. Such a packet is crowded
with respect to the active row set. All its actual nonzero positions are among
the positions counted by Z.

There is only one free terminal packet. It lies inside the final coordinate
region, so it contains at most min(a,64) positions from J. Consequently every
nonzero quiet word with active set J satisfies

    Z >= 129a-min(a,64).

This is a necessary condition only. It does not replace row-code membership,
field values, or the four exact feedback equations by independent events.

For a<=4, no nonterminal packet can be crowded, so all nonterminal inputs must
be zero. An active source row would then have at most one nonzero coordinate,
contradicting its distance 129. Thus these four source classes are impossible
for every route, independently of the mixers and feedthrough.

## Exact regional moment and the source-set union

For a fixed J, its marked positions in each region form a uniform a-subset of
{0,...,R-1}, independently across regions. There are B=R/64=128 packets in a
region. For a mark u>=1, define

    G_u(x) = sum_(j=0)^64 C(64,j) x^j u^(j * 1_[j>=5]).

Selecting j positions from a packet contributes j crowded positions exactly
when j>=5. Therefore the regional moment is the exact coefficient identity

    E[u^Z_c] = [x^a] G_u(x)^B / C(R,a).

All polynomial coefficients are nonnegative. For every x>0,

    E[u^Z_c] <= G_u(x)^B / (x^a C(R,a)).

Let d_a=129a-min(a,64). Markov's inequality, regional independence, and a union
over all C(R,a) possible active row sets give

    Pr[there exists a quiet word with exactly a active rows]
      <= U_a(x,u)
      := G_u(x)^(256B) /
         (x^(256a) C(R,a)^255 u^d_a).

There is no codeword-count multiplier. A route satisfying the support condition
is counted once per candidate row set, even if it supports many field-valued
words. In particular, this argument does not encounter the uncapped projective
first-moment obstruction merely because a bad subspace has high dimension.

## Directed finite receipt

For every integer a from 5 through 256, use x=a/16384 and the following fixed
integer value of u. No optimizer is part of the proof calculation.

| Active row count a | u |
| --- | ---: |
| 5 through 8 | 44 |
| 9 through 16 | 22 |
| 17 through 32 | 12 |
| 33 through 64 | 7 |
| 65 through 128 | 4 |
| 129 through 256 | 2 |

The helper evaluates G_u(x) as an exact rational number. It then sums the
252 positive class bounds using 160-bit Arb interval arithmetic. The result is

    log2(sum_(a=5)^256 U_a)
      = [-2236.709505746191320086292440957051669696193465
         +/- 5.92e-43].

Thus the probability of any nonzero quiet word with at most 256 active rows
is below 2^-2236.70950574. The largest class bound occurs at a=8, with
x=1/2048, u=44, d_a=1024, and log2 U_a below -2236.7095059642604.

This is an unconditional route bound, uniform over the mixer setup. The field
size does not enter it. It only requires the two stated minimum distances and
the independent uniform regional permutations.

## Joint tile-support ranks

The shared polynomial representation gives another exact interface for the
remaining quiet supports. Fix one tile, an allowed input support S of size w,
and an allowed mixed-output support T of size u. The eight input nodes and
eight output nodes must be disjoint, as they are for beta8 and beta64.

A degree-below-eight tile polynomial must vanish at the 16-w-u excluded nodes.
The space L_(S,T) of allowed tile inputs therefore has dimension

    d = max(0,w+u-8).

The tile's contribution to C extracts its four high global-basis coefficients.
Its kernel consists of degree-below-four polynomials. Restricting this map to
L_(S,T) gives

    rank(C_tile | L_(S,T)) = min(4,d),
    dim(L_(S,T) intersect ker C_tile) = max(0,w+u-12).

These formulas depend on both input and output supports, not only on the
number of output zeros. They hold for every choice of the supports, because
all excluded evaluation nodes are distinct.

In particular, a nonzero quiet packet supported in one tile satisfies
wt(X)+wt(FX)>=13. If it uses at least two tiles, the tile branch number nine
gives a stronger lower bound of 18. The minimum combined weight of the local
quiet-packet code is exactly 13: in one tile, the polynomial X(X+1)(X+2) has
input weight five, output weight eight, and zero high coefficients. This is
a local witness, not a routed outer-code word.

Feedback is summed across tiles. The rank of that sum is not the sum of these
individual ranks; feedback-image subspaces can intersect. A subsequent global
bound must retain those intersections or use a proved lower bound on their
combined rank. In particular, changing beta requires checking cross-tile
relations even though the displayed one-tile formulas remain valid.

## What remains at 40%

The quiet source space still has dimension at least

    K-4(P-1) = 917508.

Its emitted word is now FX rather than X. Consequently the old row-support
upper bound does not follow. Conversely, the large dimension alone does not
prove distance at least 40%. The sparse receipt above only shows that its
nonzero words almost surely require more than 256 active source rows.

Two coarse cuts remain inconclusive at the new geometry. A single global
prefix and one zero-state exit give the upper cap 1048640/2097152, slightly
above 50%. The paired-window mean calculation at outer[256,128], width64,
state4 gives expected shortened dimension 10.491849386572055 against 512
exit constraints. Its failure lower bound is zero, not a certificate.

For quiet words with more active source rows, the relevant event combines
outer-row membership, the routed input/output tile-support spaces L_(S,T),
and the rank of each packet's summed feedback. Dropping any of these can
lose the benefit of the shared feedthrough. Nonquiet output histories require
the state-dependent emission constraints as well. Neither class is covered
by this note.

The next bounded proof step is to extend the source-aware support or rank
bound beyond 256 active rows, retaining the cross-tile feedback ranks. The
crowded-position condition alone becomes too coarse before the guaranteed
large quiet subspace is reached; no all-density conclusion is claimed.

## Quiet distance descends to GF256

For every fixed route, the quiet subcode is defined over the embedded GF256:
the outer generator, route, preterminal feedback constraints, and feedthrough
all have coefficients in that field. No mixer appears in those constraints.
Fix a GF256 basis {1,eta} of GF65536. Write any extension-field quiet input as
X=X_0+eta X_1. Both components satisfy the base-field outer and quiet
constraints, and

    FX = FX_0 + eta FX_1,
    support(FX_j) is contained in support(FX),     j=0,1.

A nonzero input has a nonzero component, whose output is nonzero because F
is invertible. Thus its base-field quiet distance is at most its weight.
The reverse inequality follows by embedding every base-field quiet word in
the extension. Therefore the two quiet subcodes have exactly the same minimum
symbol distance. Existence with at most a active source rows also descends;
the exact number of active rows need not be preserved by a component.

Quiet-word counting can consequently use GF256 rather than GF65536. This
statement does not apply to nonquiet histories with general extension-field
mixers, whose equations need not preserve either base-field component.

The smaller counting field does not rescue the uncapped quiet first moment.
Fix 640 source rows and complete-packet windows of lengths 3264 in 205 regions
and 3328 in 51 regions. Use prefixes in the first 255 regions and a terminal
suffix in the last. There are 13107 selected packets and 838848 output positions.
On the event E that all coordinates of these rows lie in the windows, impose
zero feedback in every selected packet except the global final packet. This
costs at most 4*(13107-1)=52424 base-field equations and leaves a quiet subspace
of dimension at least 128*640-52424=29496. All other packet inputs are zero.
Tile-local feedthrough keeps its outputs inside the windows.

The exact route probability is

    Pr[E] = C(3264,640)^205 C(3328,640)^51 / C(8192,640)^256,
    log2 Pr[E] = -231841.400845614123084037... .

If Z counts low-weight GF256 projective quiet-code lines, then

    E[Z] >= Pr[E] * (256^29496-1)/255
         > Pr[E] * 256^29494
         > 2^4110.5991543858769.

This is an all-quiet witness, not an imported exit-only window argument.
It holds for either tile translation and every mixer setup. It does not
establish substantial distance-failure probability: the route event is tiny.
Base-field descent remains useful inside a capped or route-conditioned proof,
but the unconditional raw first moment still cannot certify the whole quiet
subcode.

## Reproduction

[theory/hd_blockwise_shared_cuts.py](theory/hd_blockwise_shared_cuts.py) contains
the directed receipt, exact small regional moments, and local rank checks.
Its [nine tests](theory/hd_blockwise_shared_cuts_tests.py) check:

- coefficient moments against exhaustive regional subset enumeration;
- the nonnegative-coefficient upper bound;
- all 252 fixed rational mark choices and the directed class union;
- a GF4 toy with actual outer and feedback constraints;
- the joint tile-support ranks and an attaining local quiet witness.
- exhaustive GF4-to-GF16 quiet descent and the GF256 quiet first-moment trap.

The GF4 toy has 63 nonzero quiet words: 12 use three source rows and 51 use
four. All satisfy the source/route condition, with the terminal packet free.
The separate descent toy checks all 255 nonzero extension-field quiet words
and their base components. Both fields give minimum output weight two.

From the theory directory:

    python -B -m unittest hd_blockwise_shared_cuts_tests
    python -B hd_blockwise_shared_cuts.py

All nine tests pass. No native files, benchmarks, or root indexes were changed.
