# Reusing the K18 Tail Bound with a Larger Outer

The closed RS16/S22 construction has 2048 outer groups and 64 routing regions.
A larger constituent keeps those 2048 groups while increasing the number of regions.
The inner remains the same t64/s22 map, with independent uniform GL(22,2) updates.
Its state continues across every step and region; there is no added reset or flush.

This note transfers only the uniform-input outer bound, used for occupancies q >= 3.
It does not transfer the exact outer-shell computations for q = 1 or q = 2.
The result is over the ideal setup distribution, not a guarantee for every seed.

## Why Zero State Bounds Every Entering State

Fix the active group identities and consider the uniform-input outer majorant.
Each marked group's entire output is replaced by independent uniform bits.
Its active label is retained even when some or all substituted bits are zero.
Different regions consequently have independent uniform packet values and independently sampled regional routing.
The within-group packet shuffle does not change this product distribution.

The expansion A and feedback C are fixed maps shared by every step.
Fix the routing and the per-step linear state updates in a future block of 64 regions.
Fix an arbitrary entering state a.
The marked input bits in that block are uniform, so the block output is uniform on a coset V + h(a).
Here V is the image of those input bits under the fixed block encoder, and h is the linear contribution of entering state.
Neither V nor h depends on the random input bits.

For 0 < z <= 1, let f(x) = z^wt(x) on n-bit outputs.
Its normalized Fourier coefficients are

\[
\widehat f(u)=((1+z)/2)^{n-\operatorname{wt}(u)}
              ((1-z)/2)^{\operatorname{wt}(u)}\geq0.
\]

For every binary subspace V and shift b,

\[
\mathbb E_{v\gets V}f(v+b)
=\sum_{u\in V^\perp}\widehat f(u)(-1)^{u\cdot b}
\leq\sum_{u\in V^\perp}\widehat f(u)
=\mathbb E_{v\gets V}f(v).
\]

Thus zero entering state maximizes the block's actual weight moment.
Average this pointwise inequality over the independent future routing and inner maps.
Let M(q,z) denote the resulting zero-start moment of one 64-region block.
For every entering state a, the corresponding moment is at most M(q,z).

Now concatenate m such blocks without resetting state.
Condition on the complete past before each block, including its entering state.
Future input bits, routing, and inner maps retain their original distributions.
The conditional moment for that block is at most M(q,z).
Successive conditioning therefore bounds the full moment by M(q,z)^m.
The multiplication applies to actual moments; it does not assert a row-sum property of the proof's matrix envelope.

Uniformity is essential here.
One must first apply the uniform-input majorant, retaining artificial zero values.
Conditioning all marked packets to be nonzero would destroy the subspace argument.

## The Two Larger Outers

The source outer contains four parallel RS[16,8] rows over GF(16).
At each symbol position, an independent uniform GL(16,2) map mixes the four symbols.
Its group dimension is 128 bits, its output length is 256 bits, and it has 64 four-bit packets.
Its pointwise majorant is beta times the uniform binary distribution, where

\[
\beta=\frac{2^{256}}{(2^{16}-1)^8}.
\]

For m = 2, use four parallel RS[16,8] rows over GF(256), followed by independent uniform GL(32,2) symbol maps.
Each group maps 256 bits to 512 bits and supplies 128 four-bit packets.
For m = 4, use four parallel RS[32,16] rows over GF(256), with the same independent GL(32,2) symbol maps.
Each group maps 512 bits to 1024 bits and supplies 256 packets.
Both are rate one-half, with 2048 groups at K = 2^19 and K = 2^20 respectively.
Each group is independently shuffled into its packet regions, and each region is independently shuffled across groups.

All four RS rows use the same evaluation points.
For 0 < h < n-k+1, there are no nonzero aligned words with support of size h.
For h >= n-k+1, each row shortened to a fixed h-position support has dimension r = h-(n-k) over GF(256).
A common r-position information set injects the four aligned rows into r tuples from GF(256)^4.
Exact support requires every such aligned tuple to be nonzero.
The number of fully nonzero aligned words is therefore at most (2^32-1)^r.
Independent uniform GL(32,2) maps then give pointwise expected multiplicity at most (2^32-1)^(-(n-k)).
The zero word is omitted from this active-group measure.
This proves the claimed uniform majorant, including its deliberately added artificial zero mass.

Write beta_m for the new majorant constant.
In these two cases,

\[
\frac{\beta_m}{\beta^m}
=\left(\frac{65535}{65537}\right)^{4m}<1.
\]

## Transfer of a Numerical Endpoint

Let L = 2048, D = floor(524288/10) = 52428, and lambda > 0.
Suppose a source endpoint U_q(lambda) bounds

\[
{L\choose q}\beta^q e^{\lambda D}M(q,e^{-\lambda}).
\]

The new bad-weight cutoff is D_m = floor(m*524288/10), so D_m-mD = m-1 for m = 2,4.
The preceding block argument gives the new bound

\[
U_{m,q}(\lambda)
\leq e^{\lambda(m-1)}
\left(\frac{\beta_m}{\beta^m}\right)^q
\frac{U_q(\lambda)^m}{{L\choose q}^{m-1}}.
\]

Each occupancy can use its own source-selected positive rational lambda.
Any upper bound for the actual source moment is sufficient.
In particular, the conditioned-iid method bounds that same moment by dividing an unconditional nonnegative moment by the marker-count event probability.
It need not itself satisfy any matrix row-sum or factorization identity.
Its already valid endpoint can therefore be substituted in the displayed inequality.

`packet_outer_block_transfer.py` pins the reviewed S22/K18 receipt, verifies its current mathematical source hashes and component hashes, and checks selected tail endpoints and witnesses.
It evaluates this formula with outward Arb arithmetic, then sums exact dyadic endpoints.
Its output is a derived tail receipt, not a fresh census or a complete certificate.
Separate q1/q2 bounds and their union with this tail are still required.

The retained source is `tmp/rs16-s22-k262144-whole-fresh-p256.json`, with raw SHA256
`20cd124704ff83d7ac0e098f5b80ce97e59aabf48f547bd118b1b3b957d5b60b`.
Its complete fresh-replay recipe is stored in that receipt and reproduced in the
S22/K18 section of this workstream's README.
The transfer intentionally rejects a replacement receipt until its numerical evidence has been audited.
The source receipt and its four components must therefore be preserved locally.
Numerical receipts remain outside the committed research sources.
