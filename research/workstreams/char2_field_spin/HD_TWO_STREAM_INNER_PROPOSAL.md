# Two interleaved inner-state streams: proof-first proposal

The current G16 proof comparison is dominated by repeated activation and
return of one inner state. This note proposes a different inner schedule
to test that mechanism. It does not propose promoting a construction or
transferring the existing certificates.

## Fixed schedule and cost hypothesis

The current packet map uses a sixteen-symbol state $Q$, a 64-symbol input
$X$, and fixed linear maps $A$ and $C$ with $CA=0$:

$$
Y=X+AQ,\qquad Q'=M(Q+CX).
$$

Here $M=D_2H_{16}D_1$ is the current invertible MDS refresh, with fresh
independent nonzero diagonal keys for each packet. This order matches
`finishStriped` in `HdGroupedInner.h`.

The proposed schedule starts two independent states $Q^{(0)}=Q^{(1)}=0$.
Packet $j$ uses only $Q^{(j\bmod2)}$ and updates that state by the same
packet map. The other state remains unchanged. Both states persist across
region boundaries. The routed 64-symbol packets and their internal pair
geometry remain intact; a routed pair must not be split between streams.

The number of packet projections, emissions, refreshes, and keys is
unchanged. Persistent state doubles from sixteen to thirty-two records.
This is an operation-count hypothesis, not a timing result. Keeping both
states in one loop could increase register pressure. Two separate passes
could keep only one state live, with a fixed output permutation restoring
the required order. Their routing, key, and memory costs still need review.
There is no native implementation or benchmark of this proposal.

## Required joint proof interface

Fix the source set, every constituent row set, and one complete regional
pair permutation and swap assignment. The two state trajectories then use
disjoint refresh keys. For this fixed placement, their chronological
transfer products act on a four-state product cone. Its initial state is
$(0,0)$ and its terminal vector has four ones.

For a packet assigned to the first stream, the product-cone transfer is
$T_j\otimes I$; the second stream uses $I\otimes T_j$. These are fixed-input
packet caps, not already averaged regional matrices. The route average
must retain their joint occupancy. Multiplying two independently averaged
half-region operators would generally resample the same pair permutation
and is not justified.

The joint regional provider need not enumerate every packet assignment.
Let B(x,y) be the existing two-state packet polynomial, whose variables
count singleton and full-pair slots. Alternating packet parity gives the
four-state polynomial product

$$
[(B\otimes I)(I\otimes B)]^{64}=B^{64}\otimes B^{64}.
$$

This identity holds over the shared polynomial ring: the two tensor
factors commute because they act on different states. Extract the joint
coefficient of x^s y^f only after this product, and use the original
4096-slot normalizer and orientation factors once. In particular, the
tensor entries require convolution of occupancy coefficients; multiplying
two separately normalized half-region averages would be a different law.

Each region has 128 packets and 4096 routed pair slots. Packet parity
divides these into two sets of 2048 slots. A source histogram with $s$
singletons and $f$ full pairs occupies $s+f$ distinct pair slots. For that
nonempty fixed potential source, the exact probability that all occupied slots
land in one stream is

$$
2\frac{(2048)_{s+f}}{(4096)_{s+f}},
$$

where $(n)_k=n(n-1)\cdots(n-k+1)$. This is not $2^{1-(s+2f)}$:
both active coordinates of a full pair enter the same packet and stream.
Even when the entire potential source reaches both streams, an actual
output subset can occupy only one. A common-source comparison must retain
that distinction across every mask in its column sum.

If both streams activate and subsequently return, a product-cone bound
can contain two separate flat-state return factors of size
$(256^{16}-1)^{-1}$. This is a prospective advantage, not a universal
$256^{-32}$ return assertion. One-stream paths and their exact placement
masses must remain in the same positive comparison.

## Falsifiable next proof gate

The next test should use one frozen weight marker and one fixed source
class, with no construction or witness grid. Build a joint two-stream
regional operator from the paired-slot generating function. First recover
the one-stream control when the second stream is disabled. Then compare
the two-stream operator on the same homogeneous arity-eight/F1 class,
retaining fixed row sets, all four boundary states, and the whole-word
input-weight floor.

Only a verified mathematical improvement would justify a native schedule
prototype. A future implementation needs a distinct construction identity
and fresh distance receipts. The current complete certified checkpoint
and the faster partially proved research candidate remain unchanged.
