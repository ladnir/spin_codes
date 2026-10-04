# RS[8,4] over GF16 with eight rows

This candidate replaces the parity circuit while retaining the certified
15-bit inner map and routing geometry. Its complete setup-failure bound is
**61.8840972765688264 bits** at K=65,536 and N=131,072. Runtime requires a
separate measurement; the smaller parity circuit also needs larger mixers.

Each of 512 independent groups contains eight rows of four GF16 message
symbols. Use polynomial-basis GF16 with modulus 0x13. Interpolate each row
at 0, 1, 2, 3 and evaluate it at 0,...,7, giving eight RS[8,4] rows.

At each of the eight symbol positions, collect the eight row symbols into
32 binary coordinates. Coordinate 4*r+b is bit b of row r, with b=0,...,3.
Apply an independent invertible 32-bit map that sends every fixed nonzero
input uniformly to the nonzero vectors. Divide each result into eight
contiguous four-bit packets. Thus every group supplies 64 packets.

Independently shuffle each group's 64 packets into 64 routing regions. Each
region independently shuffles the 512 group slots. Retain the exact paired
GL15 inner defined in
[`PAIRED_S15_CONSTRUCTION.md`](../k16_codesign_100us/proof/PAIRED_S15_CONSTRUCTION.md).
For input x_i and state a_i, it emits x_i+A15*a_i and updates to
M_i*a_i+A15^T*x_i. The M_i are independent uniform GL(15,2) maps. State starts
at zero, continues across all 2,048 physical 64-bit steps, and is discarded
without a final flush. All setup choices are mutually independent.

The outer is injective, and the inner is lower triangular with identity
diagonal. This defines a binary [131072,65536] linear code. Over the stated
ideal setup distribution,

\[
\Pr[d_{\min}\le13107]<2^{-61.8840972765688264}<2^{-40}.
\]

Outside that event, minimum distance is at least 13,108. This is an ensemble
bound, not a distance guarantee for any particular deterministic seed.

## Field32 implementation and one omitted mixer

Fix any binary basis of GF(2^32), and let L_c be multiplication by c in that
basis. Independently sample c uniformly from the nonzero field elements at
each symbol. The forward family L_c^T is transitive on nonzero vectors, so
the transposed encoder may use ordinary field multiplication L_c.

For fixed x!=0, equality L_c^T*x=L_d^T*x would imply L_(c+d)^T*x=0.
For c!=d, multiplication by c+d is invertible. Therefore c maps bijectively
to the nonzero outputs. Fixed invertible maps before or after the family
preserve this law. No special alignment of this field basis with GF16 is
needed for this argument.

The native 32-bit layout above uses no inherited sixteen-coordinate
permutation. Its scalar transpose applies L_c directly to the contiguous
32 coordinates. The existing byte-tower circuit uses nine byte-field
products per 64-bit payload half; this operation count predicts no timing.

Exactly one symbol mixer per group may instead be any fixed invertible
binary map, including the identity. The whole first-moment bound stays
unchanged. [`OMIT_ONE_RANDOMIZER.md`](OMIT_ONE_RANDOMIZER.md) proves equality
of the complete expected outer measures; it does not assume independence
from equal one-symbol marginals.

## Complete proof composition

Let Q=2^32. Eight parallel GF16 RS[8,4] rows form an MDS[8,4] code over a
degree-eight extension of GF16. If A_h is its symbol-weight enumerator,
the exact expected packet-support polynomial of one nonzero group is

\[
\sum_{h=5}^{8} A_h
\left(\frac{(1+15z)^8-1}{Q-1}\right)^h.
\]

For each fixed message, the active randomized symbols are independent
uniform nonzero vectors. Conditional on the induced packet-support pattern,
all active packet labels are independent uniform nonzero four-bit vectors.
Before that conditioning, packets within a symbol need not be independent.
The group shuffle makes the support uniform conditional on its size.
Fresh 256-bit outward calculations give q=1 margin 62.1003225016 bits and
q=2 margin 115.0512266186 bits. The q=2 calculation
includes all ordered support pairs and retains state across every region.

For q>=3, the retained proof uses a pointwise outer-measure majorant.
The old RS[16,8] construction has
beta_old=2^256/(2^16-1)^8. The new construction has
beta_new=2^256/(2^32-1)^4. Their exact ratio is

\[
\rho=\frac{\beta_{\rm new}}{\beta_{\rm old}}
=\left(\frac{65535}{65537}\right)^4<1.
\]

The comparison inputs, regional geometry, and inner operators are identical.
For every q>=3, multiplying any authenticated old endpoint by rho^q gives
a valid new endpoint. This uses independent complete group measures, not
only symbol marginals. The old endpoint's selected tilt remains a valid
witness for the new construction.

`certify_rs8x8_s15.py` verifies the old receipt hash, all 106 old source pins,
its exact final sum, and the complete geometry. It freshly regenerates the
s15 map and checks exact equality with the old map record. Each transferred
tail endpoint is scaled as an exact rational and rounded upward to 256 bits.
The weakest transferred individual tail contribution is q=41 at
66.2433763960 bits.

The final union sums all 512 positive dyadic endpoints exactly, rounds once
upward, and compares with 2^-40 using integer arithmetic. The result is a
complete bound composition. Its `fresh_replay=false` flag records that the
tail endpoints were authenticated and transferred, not recomputed.

The receipt is `rs8x8-s15-complete-p256.json`, SHA-256
`34f38f58b4e29c5c50245f80d8a654db216671637ecdd254fdb95608ec707b2a`.
It references the frozen s15 receipt and pins the new coordinator separately.
The old 106-source receipt and every old source remain unchanged.

Run from the repository root with an unused output filename:

```text
python -B research/workstreams/k16_outer_routing/certify_rs8x8_s15.py --output research/workstreams/k16_outer_routing/new-rs8x8-s15-p256.json
python -B -m unittest discover -s research/workstreams/k16_outer_routing -p test_*.py -v
```

Independent tests check every rational tail conversion, exact union and
rounding, source pins, and q=1 through 32 physical steps per region instead
of the coordinator's 16 two-step macros. The direct physical calculation
agrees within 2^-200 relative error.
