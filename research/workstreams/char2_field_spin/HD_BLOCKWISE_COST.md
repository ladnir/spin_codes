# Cost interface for a blockwise SPIN restart

The reusable performance structure is a packet recurrence with one fused
route, fixed structured maps, and a small carried state. Do not carry the
scalar pre-prefix into the restart by default. Keep rectangular packet
geometries in scope: input width t and state size s serve different purposes.

This note gives implementation costs and design choices, not new timing or
distance results. No native build, remote command, performance benchmark,
or existing-file edit was performed.

## What the previous kernels actually establish

The frozen GF256 winner uses four input and four state records, keeps three
output records, and uses a normalized full-Cauchy map. Its approximately
3.2 ms receipt is for the K20 **transpose**, at about 83.33% nominal rate
and 10% certified distance. It is not a half-rate, high-distance forward
baseline. `RateSpeed4Inner.h` and `RateSpeed4InnerNext.h` are useful for their
explicit four-record packing and bounded packet pipeline, not as a prediction
for larger states.

The old wide-field prototype in `HighDistanceWide.h` also uses width/state
four. Its keep-two and keep-one controls compensate for puncturing with
4K and 8K raw outer expansion at final rate one half. That expansion, and
the associated route and outer work, should not be inherited merely because
the arithmetic is reusable.

The original IMT implementation has a different recurrence:

    y = x + A*q,       q' = F*q + B*x.

`spin/src/kernels/Inner.h` retains raw packet inputs for the separate feedback
map B. `ImtRounds.h` implements binary rank-one transvections with XOR
reductions and XOR updates. This is not the same randomness as fresh
coordinatewise field diagonals. Both the fixed maps and refresh law are
valid choices to revisit; neither should be silently replaced by a
full-Cauchy packet.

## Field arithmetic that can be reused directly

Use 128-bit records containing eight parallel GF65536 coordinates. The
existing tower is AES GF256[u]/(u^2+u+0x20). A genuine GF65536 coefficient
a+b*u acts as

    (x0,x1) -> (a*x0 + 0x20*b*x1, b*x0 + (a+b)*x1).

`WideDiagonal16.h::multiply512` processes four records with two GFNI products,
two coefficient expansions, a byte-pair swap, and an XOR. Setup stores four
bytes per coefficient. The raw scalar is two bytes, but it is not the
representation consumed by this kernel.

A fixed coefficient in the embedded GF256 subfield multiplies both tower
bytes by the same byte. It therefore remains one bytewise GFNI operation,
or the existing constant-affine equivalent. Widening the sampled coefficient
alphabet does **not** require two-GFNI multiplication on every fixed-map edge.
Nonzero minors of fixed GF256 matrices stay nonzero over GF65536. This
algebraic fact does not transfer a whole-code distance certificate.

For a candidate with one fresh input diagonal and one fresh state diagonal
per packet, the variable-product count is

    2*ceil(t/4) + 2*ceil(s/4) GFNI512 instructions per packet.

For the divisible-by-four geometries below, the coefficient traffic is
4(t+s) bytes per packet. The table excludes fixed maps, XORs, routing,
payload access, spills, and any additional refresh coefficients.

| t / s | Variable GFNI per packet | Per input record | Coefficient bytes per input record |
|---|---:|---:|---:|
| 16 / 16 | 16 | 1 | 8 |
| 16 / 32 | 24 | 1.5 | 12 |
| 32 / 16 | 24 | 0.75 | 6 |
| 32 / 32 | 32 | 1 | 8 |
| 64 / 16 | 40 | 0.625 | 5 |
| 64 / 32 | 48 | 0.75 | 6 |

Larger t amortizes state refresh and packet boundaries. It does not make
the state carry stronger, and it must pass the separate window/rank gate.
The largest state in this table occupies only eight ZMM registers, but
maps, coefficients, and packet temporaries also compete for registers.

## Keep three distinct fixed-map choices

**Sparse IMT-style maps.** Binary sparse A and B use XORs and fixed wiring.
They are the cheapest baseline if their required q-ary local ranks can be
proved. They need not be square, symmetric, or Cauchy. A field extension
does not repair a deficient binary minor.

**Two rectangular structured maps.** A useful comparison is

    x = Dx*v,    r = Ds*q,
    y = x + A*r,    q' = r + C*x,

with A of size t-by-s and C of size s-by-t. The existing residue-coordinate
maps supply scaled-Cauchy examples for t>=2s; the symmetric choice has
CA=0 in those particular affine-subspace constructions. Other rectangular
pairs are possible, but their birth, quiet, return, and joint-zero ranks
must be established for the chosen pair. Identity feedthrough gives cheap
causal input/output invertibility, but weak birth diffusion.

**Normalized full-Cauchy packet.** Start from a Cauchy matrix of order t+s
and normalize its state-to-state block to identity:

    [ y  ]   [ F  A ] [ Dx*v ],
    [ q' ] = [ C  I ] [ Ds*q ].

Here F is generally not identity. The input mix F*x is a third fixed map
unless a checked algebraic factorization combines it with the other two.
The small winner exploits such special structure; do not assume its
constant counts extend to t=16,32,64. The selected-output/all-state rank
properties of this normalization are not equivalent to full superregularity
of its identity-containing block matrix.

A full-Cauchy matrix of order t+s can be partitioned for t!=s. For the
geometries above, at most 96 row labels and 96 disjoint column labels are
needed, so GF256 has enough labels. Equal packet and state dimensions are
not a field-size or implementation requirement.

## Fixed-map kernels and a concrete schedule

The useful exact-map starting points are:

- `Cauchy256InnerFft.h`: compile-time forward, inverse, and transposed
  additive-FFT butterflies in the original residue-state coordinates.
- `Cauchy128InnerLaneFft.h`: a representation permutation that makes small
  stages lane-local. Its persistent state and diagonal layouts must agree.
- `HdFastOuter.h` and `RateSpeed4OuterNext.h`: explicitly unrolled
  sixteen-register leaves and four-row interleaving for fixed GF256 maps.
- `WideDiagonal16.h`: the genuine-field input/state diagonal kernel.
- `RateSpeed4InnerNext.h`: bounded independent work around a chronological
  state recurrence. Reuse its scheduling idea, not its packet-specific map.

For interpolation-based rectangular A, factor it as one s-point inverse
transform followed by t/s s-point coset evaluations. Its transpose reverses
those operations and sums the coset contributions. This requires t>=s and
t divisible by s; it is one implementation family, not a restriction on
all admissible packet maps. Constants and inverse transforms are prepared
or compiled once, outside encoding.

Naive dense multiplication of both rectangular maps costs about ts/2 packed
fixed products per packet when all dimensions are divisible by four. A
dense input mix adds about t^2/4. These counts explain why simply enlarging
the dense four-state map is unattractive. The structured transform cores
instead have scalar butterfly work proportional to (t+s)*log2(s), although
packing and shuffle costs prevent dividing that count by four blindly.

Use the following forward schedule for a candidate with diagonal refresh:

1. Gather a packet from the outer planes and apply Dx immediately.
2. Compute its input-only fixed maps while the previous packet uses the state.
3. Apply Ds to the carried state once per packet.
4. Emit from the required old/refreshed state and update the carried state.
5. Store output once; retain the state across region boundaries.

Use four chronological packets only where register pressure permits.
At t=64, four uncompressed packet inputs alone need 64 ZMM registers.
Instead, process sixteen-coordinate tiles with one or two live packet
workspaces, and keep only the s-coordinate moments needed by the recurrence.
Fixed-size aligned scratch is acceptable when the map genuinely needs it;
a speculative all-register schedule that spills invisibly is not better.

One regional route can be fused with input reads and Dx. There is no
pre-prefix stage, no extra routed-output pass, and no per-packet allocation.
For K20 direct 2K expansion, a 64-plane outer has 512 KiB of payload per
region; a 128-plane outer has 256 KiB. This is a layout fact, not a claim
that either region fits all concurrent working data in cache.

An implementation should start with four-record grouped addressing and
explicit butterfly leaves, then inspect register lifetimes. Do not substitute
generic callback loops or runtime matrix dispatch in the inner kernel.

## Refresh alternatives and the next gate

Fresh state diagonals cost s full-field record products per packet. A
q-ary rank-one update F=I+u*v^T needs s products for the dot product and s
for its broadcast update, plus a reduction. Its invertibility condition is
1+v^T*u!=0. Under this direct schedule, one refresh costs about s packed
GFNI instructions, twice the variable product count of a state diagonal.
Sparse u or v can reduce work but changes the distribution. The binary
IMT implementation avoids these products through binary masks; its local
probability law must not be reused for a different q-ary sampler.

The first bounded comparison should retain at least one rectangular choice,
such as t/s=32/16 or 64/16, alongside any proof-motivated larger-state choice.
Compare sparse IMT-style fixed maps, two-map residue Cauchy, and full-Cauchy
with its feedthrough cost explicitly charged. Do not choose the family
solely by a low multiplication count: the earlier exact residue FFT reduced
GFNI operations yet lost to the dense kernel because shuffles and dependency
chains dominated.

Before native work, freeze the packet equation, field, coefficient sampler,
state basis, output count, and required coordinate ranks. Then validate
dense scalar basis cases and nonzero-boundary forward/adjoint identities.
Only a geometry/proof survivor should receive a new serial native cost
measurement. The target below 5 ms remains unmeasured for every restart
candidate described here.
