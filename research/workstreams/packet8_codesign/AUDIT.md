# Independent audit of the byte-packet prototype

This audit covers the initial implementation and the additive variants listed below
in `spin/experiments/packet8_codesign`.
It checks algebra and layout, not a distance certificate or performance target.
No earlier four-bit certificate is asserted for this construction.

## The exact local maps

Let F be GF(256), represented in the polynomial basis modulo `0x11b`.
A physical input consists of eight field symbols x_0,...,x_7, or 64 binary coordinates.
The persistent state is (a,b) in F². The evaluation points are alpha_i=i.
Each physical step emits

```
y_i = x_i + a + alpha_i*b
(a_next,b_next) = M(a,b) + (sum_i x_i, sum_i alpha_i*x_i).
```

The initial state is zero, state persists across regions, and the final state is discarded.
Every physical step has its own sampled matrix M in GL₂(F).
All additions in this display are field additions, hence binary XORs.

Write the fixed expansion as A and feedback as C.
Both have binary rank 16, and CA=0 for points 0 through 7.
Indeed, the entries of CA are sums of 1, alpha_i, and alpha_i²; all these sums vanish.
Each one-byte restriction of C has rank 8; each two-byte restriction has rank 16.
Distinct evaluation points establish the second statement through the determinant alpha_i+alpha_j.

For any three distinct packet positions, the feedback kernel has binary dimension 8.
Every nonzero vector in that kernel has all three byte labels nonzero.
Thus independent uniform nonzero labels cancel with probability 255/255³=1/65025.
This is a local cancellation probability, not a whole-code failure probability.

The exact expansion spectrum begins with 8 words of weight 8 and 48 words of weight 12.
There are 124 words of weight 16. The full 65,536-state census is available from
`independent_checks.py`; its coefficients sum to 65,536.
Changing evaluation points alone cannot eliminate the weight-eight constant words:
for b=0 and a a binary unit vector, all eight output symbols equal a.

## Binary adjoints and chronology

The inner's transpose is taken under the binary coordinate dot product.
Field multiplication is not self-adjoint in the AES polynomial basis.
For multiplication by c, its binary adjoint maps a byte z to

```
adjoint_c(z)_j = dot_binary(c * 2^j, z),  j=0,...,7.
```

The conventional GFNI matrix operand places `c*2^j` in byte 7-j.
The reverse state update swaps the off-diagonal field-matrix blocks and applies
this binary adjoint to every scalar multiplication.
For a forward step numbered t, the reverse step uses the matrix M_t, not M_(t-1).
Its incoming reverse state is the coefficient of forward state t+1.
The update after reverse step zero is discarded.

The native prototype's forward feedback is the field transpose of its expansion,
not the binary transpose. Its scalar and SIMD reverse paths account for this distinction.
The inspected pair-sum feedback computes the same weighted raw-input sums as the literal oracle.

## State-update distribution

Under ideal independent uniform input words, sampling four independent field bytes
and rejecting zero determinant samples exactly uniformly from GL₂(F).
This group has `(256²-1)(256²-256)` elements.
For every fixed nonzero state and every nonzero target, exactly `256²-256`
group elements map that state to that target.
One proof completes each nonzero state to a field basis and maps one basis to the other.
Thus a fresh update sends any fixed nonzero state uniformly to the 65,535 nonzero states.

Consequently, for fixed A, C, and input stream, fresh GL₂(F) and GL(16,2)
updates induce the same fixed-message output distribution.
The statement follows step by step, since both maps fix zero and have the same nonzero transition law.
It does not identify the distributions over entire binary codes or over multiple messages jointly.
Seeded implementation streams follow the existing ideal-randomness convention;
the audit does not claim independent ideal samples from every deterministic seed.

## Routing and native outer layout

Every 128-to-256-bit outer group has 32 eight-bit packets.
The setup independently shuffles these packets and independently shuffles group slots
in each of 32 regions. It does not join entries from an existing 64-region route.
Each group contributes exactly one packet per region.
The padded group stride is 260 binary coordinates. Only coordinates 0 through 255 are written.

The retained native outer combines four GF(16) RS[16,8] rows with independent
nonzero GF(2¹⁶) symbol multipliers. The fixed four-by-four coordinate transpose
is already absorbed into its literal symbol matrices.
Independent tests confirm its coordinate cancellation and the resulting binary-adjoint symbol map.
The fast outer receives two adjacent packed byte vectors per symbol and does not repack the entire scratch.

The audit found one packing-oracle mismatch. The fast GFNI transpose writes payload
bits 7 through 0 into successive lanes; the initial scalar pack/unpack used bits 0 through 7.
The implementation author corrected only those scalar helpers. The fast layout was retained.
An independent GFNI emulation now checks all 1,024 raw basis bits and their inverse packing.

## Optional scaled expansion: not implemented

A bounded exact census examined fixed scales

```
c_i = 3^i = [1,3,5,15,17,51,85,255] in F.
A'_i(a,b) = c_i * (a + i*b).
```

This expansion has binary minimum weight 17, with three words attaining that weight.
It does not constitute a whole-code distance result.
Two feedback choices deserve separate treatment:

1. Keep the original C. Its one- and two-packet ranks stay 8 and 16; CA' has rank 16.
   The feedthrough encoder remains invertible because its diagonal input block is the identity.
   Its reverse emission is unchanged; its raw-input feedback acquires the fixed scales.
2. Set C'_i(x)=(c_i^-1*x, i*c_i^-1*x). Then C'A'=0 and the same packet ranks hold.
   The reverse emission and feedback both acquire fixed scales.

The first choice needs up to 16 extra 512-bit GFNI operations per physical step
in a direct implementation: eight packets, two payload halves, one scaled feedback input.
The second needs up to 32 additional GFNI operations. Scheduling may reduce these counts.
No such variant has been implemented or timed in this audit.

A cheaper A-only choice uses scales `[1,1,1,1,15,15,15,15]`.
Its exact expansion minimum is 16, attained by eight words.
Writing L_c for binary-adjoint multiplication by c, its reverse feedback can use
five pair-affine operations, instead of the baseline's two. Define

```
H = x4+x5+x6+x7,  P = x6+x7,
O_lo = x1+x3,     O_hi = x5+x7,
h = L_15(H),
D0 = x0+x1+x2+x3+h,
D1 = O_lo+L_15(O_hi)+L_2(x2+x3+L_15(P))+L_4(h).
```

The three additional pair-affine operations represent six extra 512-bit GFNI
instructions per physical step. This is an instruction-count estimate, not a timing.

Postcomposing a uniform nonzero 16-bit symbol image with fixed invertible byte maps
preserves its fixed-message uniform law. This fact alone does not remove instruction cost.
Arbitrary independent byte scalings need not remain in the sampled GF(2¹⁶)-multiplier family.
Any proposed absorption into that outer must implement the composed binary map exactly,
account for the actual route-dependent byte positions, and measure the resulting kernel.

## Checks and source scope

Run the independent tests from the repository root:

```text
python -B -m unittest discover -s research/workstreams/packet8_codesign -p test_independent.py -v
```

The tests use an independent Python implementation and import no native kernel or proof evaluator.
They cover every AES multiplier on binary basis inputs, the binary matrix adjoint,
exact small-field GL₂ transitivity, multi-step chronology, local ranks, a complete
expansion census, native outer layout, true byte routing, packing, and a whole
RS-outer/route/inner forward-transpose dot-product check.
The small-field enumeration illustrates the separate group-count argument above.
All 14 independent tests pass, including the added half-payload store and prefetch-boundary checks.

Historical first-snapshot hashes, after the scalar packing correction:

| File under `spin/experiments/packet8_codesign` | SHA-256 |
|---|---|
| `Packet8.h` | `7b74e223658bb4fdbff3e6a37ea303b8a72b42164b507333db71b0640d3a9bb4` |
| `Setup.cpp` | `af1cb030bca3300208977c9ec13458611aeaa3604b5406408aebf8bce335a556` |
| `Scalar.cpp` | `c2fd1995069592b596ff6b0f315345889922633f3f2d9a30848263a0f705ac55` |
| `Fast.cpp` | `9b0487054f8bf98ad4d58a5ca6d5e9110da56a564dbc7dd10d8c769c127446bf` |

These first-snapshot hashes do not identify the later header and prefetch additions.
The final reviewed additive snapshot has the following hashes:

| File under `spin/experiments/packet8_codesign` | SHA-256 |
|---|---|
| `Packet8.h` | `1614e9b5daf0fe89ec196b8e6c1872fafc858eb4cc9b10b1169d2018f3825423` |
| `Setup.cpp` | `af1cb030bca3300208977c9ec13458611aeaa3604b5406408aebf8bce335a556` |
| `Scalar.cpp` | `c2fd1995069592b596ff6b0f315345889922633f3f2d9a30848263a0f705ac55` |
| `Fast.cpp` | `dc4ad57ba1d5b8bcc637dc93fb5e2ac8c63c4ddad55cd4ce92907d00c3d12c7b` |
| `OuterShared.cpp` | `3a6f93b5107f661c9e4885b0d18ece09e5e7356557b9d0f2ad6e57fbb64f76bc` |
| `OuterHalf.cpp` | `a362fbd752b7edd03bd688f26f07d7625f38bd60450cb0b6a23fd636f8c1333d` |

`OuterShared.cpp` retains the baseline field randomizer and four-plane layout.
It then calls the frozen shared-parity circuit on planes 0/1 and 2/3.
The reused `k16_codesign_100us/outer_variants/OuterVariants.cpp` has SHA-256
`f31c1ca780d5cab7aecb6d83e41186990167ebaf0b7fa78112d0b8e9a74e9006`.
The shared circuit replaces repeated products by their exact GF(16)-linear sums.

`OuterHalf.cpp` applies the same circuit separately to the two 64-bit payload halves.
Its compact plane index is `2*symbol+plane`, corresponding to original plane `2*plane+half`.
All input, parity, and systematic offsets follow this correspondence.
The two masked stores write disjoint low and high qwords of each output block.
Single-input VBMI ignores selector bit 6; masks `0x55` and `0xaa` select the intended source half.
The independent selector test compares these stores with the original two-source store
and verifies that all 2,048 output bytes are written exactly once.
The retained `k16_codesign_100us/kernel/NativeOuterFinish.h` has SHA-256
`1bd655e93756c1fac0bc6e281b512ec98308fa4967665796c4c01be02ca52678`.

The prefetch variants preserve the baseline arithmetic, update indices, and reverse traversal.
For lookahead P packets, the prefetched epoch starts at `8*epoch-P` after decrementing the epoch.
The loop guard ensures this index is nonnegative; the final loop visits every remaining epoch once.
The independent schedule test covers all three lookaheads and lengths on both sides of each boundary.
Scratch remains disjoint from input and output, as required by the native interface.

No remaining algebra or layout issue was found in the reviewed snapshots.
Native scalar/SIMD tests and timings are owned by the parent task and are not reported here.
The next implementation decision should use those correctness results and matched timings.

## Why exact return exclusion cannot close the first middle-occupancy gap

The first finite-family screen fails near occupancy 119. The following diagnostic
isolates which local relaxation could explain that result. All reported numbers
are floating evaluations, not outward certificate endpoints.

At tilt 0.4, the full comparison gives -4150.93464 bits. Its all-zero-state
paths alone give +1051.75184 bits. Keeping zero states and one-step birth families,
but removing the uniform-state coordinate, still gives -2083.72832 bits.
Keeping only returns at input occupancies one and two gives +1012.87802 bits.
Suppressing all returns gives +880.39300 bits, but this suppression is not a valid bound.
These observations concern the same beta-uniform outer comparison; they are not codeword evidence.

The current nonzero-source operator bounds return mass by T/(2^16-1), where T
is its total weighted emission. An exact return also requires CX nonzero.
The omitted condition has a uniformly small relative effect at this tilt.

To see this, fix an entering state and an active packet support. Write z=exp(-tilt).
Conditional on weighting by z to the emitted binary weight, the active bytes remain independent.
For byte offset v, their normalizer is

```
Z_v = sum_{x != 0} z^wt(x+v) = (1+z)^8-z^wt(v) >= Z_min = (1+z)^8-1.
```

Every tilted byte atom is at most eta=1/Z_min. Fixing every byte except two
determines at most one pair that gives CX=0, because each two-byte feedback
restriction is invertible. Consequently the omitted weighted CX=0 mass is at
most eta² times T. Occupancies one and two actually have zero omitted mass.
Mixtures over supports and source-state families preserve this inequality.

Replacing the return entry by its exact value therefore changes each local
matrix entry by a factor no smaller than 1-eta². There are 2048 physical steps.
At tilt 0.4, the total possible improvement from this subtraction is at most
`-2048*log2(1-eta²) = 0.83220` bits. At tilt 0.6 it is at most 2.86665 bits.

There is a broader lower comparison for the entire local closure. One injective
active byte bounds every tilted syndrome atom by eta. For a nonzero source, the
true weighted density at target u equals `(T-weighted_mass(CX=u))/(2^16-1)`.
Thus every nonzero target receives at least 1-eta times the uniform density
used by the closure. The zero target satisfies the same inequality when a return
is permitted. At zero occupancy, the fresh uniform nonzero state law is exact.
Zero-source births are represented exactly by the retained birth families.

These pointwise inequalities hold for each represented source family. Positivity
then propagates the lower comparison through any physical-step sequence and its
placement average. The complete local closure can lose at most
`-2048*log2(1-eta) = 50.00396` bits at tilt 0.4, or 93.47278 bits at tilt 0.6.
The missing thousands of bits cannot come from this local relaxation at those tilts.
They must be addressed through the outer comparison or construction rather than
by an expensive exact count of the excluded return events alone.

The argument does not upper-bound slack from replacing the true outer distribution
by the beta-uniform envelope. It does not prove any actual code has low distance.
`zero_sector.py` reproduces the independent scalar all-zero path calculation,
selected operator diagnostics, and the analytic slack estimates.
The independent tests exhaustively verify the local lower sandwich in a GF(4) fixture.
