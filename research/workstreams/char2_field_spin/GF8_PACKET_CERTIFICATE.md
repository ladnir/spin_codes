# A cheaper eight-record GF256 inner: complete K20 certificate

The repaired cubic inner closes a new higher-rate checkpoint at K20:
**10% symbol distance, requested rate 64/99 (64.646%), and at least
65 bits of whole-code setup-failure margin.** Its proposed packed inner
uses 18 GFNI512 instructions per sixteen outputs instead of the current
70. This is a complete distance certificate and an algebraic instruction
schedule, **not a measured encoding-speed improvement**.

The new point trades a little rate for substantially less proposed inner
work. It does not replace the existing GRS[128,87] or [128,88] checkpoints,
and it does not meet the aspirational 3/4-rate target. No production code,
native encoder, benchmark, or existing certificate was changed.

## Construction and probability space

Use the AES representation of GF256. Repeat the GRS[128,83] outer across
L=12672 rows, padding a requested 1048576-symbol message to 1051776 source
symbols. Collect each outer coordinate into its own region and sample
one independent uniform permutation of the L positions in each region.
This is one regional routing stage. There is no second global route.

Process the resulting N=1622016 symbols in packets of eight. The state
has four field coordinates. Cycle through the eight phases in the
[exact map reference](theory/gf8_packet_maps.py), starting at phase zero.
For each packet phase, its input labels u and residue nodes v define

\[
A_{u,v}=(u+v)^3+1,\qquad C_{v,u}=(u+v)^3+3.
\]

Here addition of byte labels is XOR, and the constants are field elements
with those byte representations. In particular, C is not the transpose
of A. Let X be the packet after its input-coordinate diagonal and let R
be the incoming state after its state-coordinate diagonal. The update is

\[
Y=X+2AR,\qquad Q'=R+CX.
\]

Every diagonal entry is independently uniform in GF256 excluding zero.
All diagonals and regional permutations are independent, sampled once
during setup, and fixed across encodings. Initialize the state to zero,
carry it across region boundaries, and discard the final state.

Over this setup distribution, the probability that some nonzero scalar
message produces fewer than 162202 nonzero output symbols is at most
2^-65. Since 162202 is greater than .1N, this proves the stated distance.
The same scalar code may operate on all sixteen byte lanes of each
128-bit record. A nonzero record message has a nonzero scalar lane, so
the scalar distance guarantee also bounds record support. This uses one
shared good-code event, not a sixteen-fold failure union.

## Why the new maps fit the proof

The [map audit](GF8_PACKET_MAPS.md) verifies all 7904 square minors of A
and C across the eight phases. They are nonzero, and CA=0 in every phase.
These properties give the emission, return-rank, and packet-invertibility
conditions used by the generic transfer. They do not require C=A^T.

When a nonzero state returns to zero, its output is `(I+2AC)X`.
All eight phases have the same returned matrix in their local coordinate
ordering. Exact support/rank enumeration gives minimum output weights
`[8,6,5,4,3,2,2,1]` for exact input weights one through eight. These are
finite algebraic statements, not sampled output statistics.

This repairs a genuine interface mismatch in an earlier cheap cubic
candidate: its old injection map had singular coordinate minors. Merely
reusing that circuit with the newer proof would not have been justified.
Changing the injection constant to three repairs the minors while
retaining the five-product packed circuit.

The [transfer note](GF8_PACKET_WIDTH_SCREEN.md) explains the joint return
and output-zero rank bound. An independent audit checked that the source
rank count does not assume random outer images, and that state and input
diagonals are charged separately without reusing their randomness.

## Complete coverage

Let a count active outer rows and let W count nonzero outer coordinates.
Every nonzero GRS row has weight at least 46.

| Component | Complete covered classes | Directed margin |
|---|---|---:|
| Sparse | a=1 through 1057 | 65 bits |
| Middle | W=48661 through 389283 | 3303 bits |
| Dense | W=389284 through 1622016 | 212483 bits |

If a is outside the sparse interval, then W is at least
46*1058=48668, which is already inside the middle interval. The middle
and dense ranges are adjacent and reach N. Thus every nonzero message is
covered. The [whole-code replay](theory/gf8_packet_certificate.py) checks
the construction metadata, both joins, and the probability sum; it does
not report the minimum component margin as a union bound.

The sparse replay uses exact rational packet matrices and outward Arb
arithmetic at 4096-bit precision. It covers all active counts with 147
positive intervals. The middle replay covers 210 cells with outward
80-digit decimal arithmetic. The dense replay covers 152 cells and uses
the inner-independent causal-inversion theorem. See the
[component note](GF8_PACKET_MIDDLE_CERTIFICATE.md).

Partition identifiers, in sparse/middle/dense order:

```text
c98422c5020d7f2615af386c9605d20aaa207d1200d12315ceca1b794842f235
a6bb9147e34372abfd664a6c75d5ea2bcbdec8c6fb3dbcf7bd99d37828216c9f
efcdb9e1ded74c57ee496378e5665d96820755cf5b3aa2079aab01e5a2628b5b
```

## Performance hypothesis, not timing

The proposed paired circuit packs two injection blocks and two emission
blocks together. The scalar reference also checks the actual record-lane
transpose schedule and the asymmetric adjoint.

| Adjoint inner arithmetic per sixteen outputs | Current dense width16/state8 | New packed width8/state4 |
|---|---:|---:|
| GFNI512 instructions, including diagonals | 70 | 18 |
| XOR/XOR3 instructions | 38 | 32 |
| Record-lane shuffles | 24 | at most 32 |

The new forward map uses thirty XOR/XOR3 instructions per sixteen outputs.
Its adjoint needs two more: each packet combines its injections before eta scaling, then adds the incoming state.
Both directions retain the stated GFNI512 and record-shuffle counts.

The new profile has 4.762% more output records than the certified K20
GRS[128,87] point: 1622016 versus 1548288. After that expansion, the proposed
inner GFNI count is about 73% lower. Outer work, routing, coefficient
access, register pressure, and memory traffic are not included in that
reduction. The extra shuffles are a particular native-performance risk;
earlier FFT schedules showed why fewer field products alone are not enough.

The next implementation should preserve the explicit five-product
packing, use the asymmetric transpose formulas, and benchmark the whole
encoder against both the current GF256 and binary baselines. Start with
this certified profile before optimizing more outer dimensions. K16/K18
can then be checked and tuned separately; no new certificates or timings
for those sizes are claimed here.

## What the other screens showed

The repaired maps substantially improve the sampled middle-weight
frontier, but sparse classes still limit the present proof. At K20,
dimension 84 has only 22.70 floating bits at a=128; dimension 85 gives
-359.82 there. These are upper-bound failures, not low-distance witnesses.
Dimension 83 was selected to close a useful implementation candidate,
not to exhaust the rate frontier.

Width 32 with state 8 failed the tested middle classes. Doubling its state
helped the proof but raised arithmetic costs. The separate
[local-mixing study](GF8_LOCAL_MIXING.md) found a precise pair-mixing law,
but output-only mixing adds work and does not repair deficient return
equations. Neither alternative is selected for a new kernel.

## Replay

From the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/gf8_packet_certificate.py
```

The full replay rechecks the finite maps and all three coverage components.
From the theory directory, the new test suite is:

```text
python -B -m unittest gf8_packet_maps_tests gf8_packet_width_tests gf8_packet_portability_audit_tests gf8_packet_middle_certificate_tests gf8_local_mixing_tests test_gf8_packet_sparse_screen test_gf8_packet_certificate
```

All 46 tests pass. They include exact map and circuit checks, independent
rank/multiplicity tests, directed arithmetic, and rejection of mismatched
geometry, incomplete receipts, and coverage gaps. Native timings remain
the next gate; all benchmarks must run serially.
