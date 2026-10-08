# Cost of the shared eight-point packet circuit

The proposed circuit adds dense mixing within each eight-record input tile
while sharing the transforms used for emission and feedback. It is a new
code, not an instruction-selection change to the measured identity-feedthrough
code. The later [native comparison](HD_BLOCKWISE_SHARED_NATIVE.md) records
the beta64 FFT and cube controls; this note retains the initial cost model.

This note retains the first FFT-based cost control. The subsequent
[cube circuit](HD_BLOCKWISE_SHARED_FACTOR.md) implements the shift64 map
with59 packed GFNI arithmetic calls, versus82 for its FFT control. It shares
cube sums with feedback and uses a packet-wide state-emission polynomial.
Both counts exclude packing, shuffles, and memory costs; neither is a timing.

The measured reference remains the GRS[256,128], t16/s4 control at
6.352 ms for K20 forward encoding. Its full call includes the outer,
regional routing, staging, and output stores; see the
[native checkpoint](HD_BLOCKWISE_NATIVE.md). Neither that code nor this
proposal has a complete 40%-distance certificate.

## Exact circuit and representation

Let T be16,32, or64, and let the state contain four GF65536 symbols.
Each128-bit record carries eight independent lanes of this scalar code.
Let phi_0,...,phi_7 denote the global normalized additive-polynomial basis
already used by the FFT helpers. For alpha in{0,8,...,T-8}, define

    E_alpha[i,j] = phi_j(alpha XOR i),       0 <= i,j < 8.

The four state coordinates Q are **global polynomial coefficients**, not
the previous values at nodes128..131. For the routed tile X_alpha, compute

    z_alpha = inverse(E_alpha) X_alpha;
    h       = XOR_alpha z_alpha[4:8];
    Y_alpha = E_(alpha XOR 8) (z_alpha + (Q,0));
    Q_next  = M Q + h.

Use the same incoming Q for every tile. Capture each high coefficient vector
before adding Q; only the low four coefficients change. Initialize Q=0,
carry it across regions, and discard the terminal state.

The mixer is M=I+u*v^T over GF65536. Independently for each packet, sample
u uniformly from the nonzero vectors and v uniformly subject to v^T*u!=1.
Sampling and coefficient preparation remain outside encoding.

Equivalently, Y=F X+A Q and Q_next=M Q+C X. The blocks of F are
E_(alpha XOR8)*inverse(E_alpha); A stacks the first four columns of the
shifted evaluation matrices; C concatenates the last four inverse rows.
All blocks of F equal the same invertible eight-point Lagrange translation.
It has branch number9 and satisfies F8^2=I. This is not a binary Hadamard
transform or a separate permutation of the records.

The [shared-map analysis](HD_BLOCKWISE_SHARED_MAPS.md) gives C=H*A^T with
an invertible fixed H, CA=0, emission parameters[T,4,T-3], and check-kernel
parameters[T,T-4,5]. An independent GF256 matrix replay checked the common
F8, its involution, the high-coefficient formula, H invertibility, and CA=0
for all eight tiles of T64. These identities persist over GF65536.

## Reused primitives and counts

Use the packed `transform<8,Alpha,Kind>` schedule in
`HdBlockwiseSmallInner.h`: two ZMMs per tile, with the existing compile-time
butterflies. Global offset0 costs three GFNI products per transform;
the other multiples of8 below64 cost five. Both the input offsets and
shifted output offsets contain zero exactly once per packet.

Thus the fixed maps cost 10*(T/8)-4 GFNI products. One four-state rank-one
refresh adds four, using two calls to `WideDiagonal16.h::multiply512`.
The counts below exclude gathers, XORs, shuffles, loads, stores, and setup.

| Circuit | Outputs per packet | Fixed GFNI | Refresh GFNI | Total GFNI |
| --- | ---: | ---: | ---: | ---: |
| Measured identity-feedthrough t16/s4 | 16 | 18 | 4 | 22 |
| Shared eight-point circuit, T16 | 16 | 16 | 4 | 20 |
| Shared eight-point circuit, T32 | 32 | 36 | 4 | 40 |
| Shared eight-point circuit, T64 | 64 | 76 | 4 | 80 |

For64 outputs, the measured reference executes four packets:72 fixed-map
and16 refresh products, totaling88. The T64 proposal uses76+4=80,
an8-product reduction, about9.1%. The structural gain is dense input
mixing at this count, plus one refresh boundary instead of four.
It is not a large arithmetic reduction by itself.

Prepared refresh coefficients occupy32 bytes per packet: four u entries
and four v entries, each using the existing four-byte tower representation.
For64 outputs, T64 therefore reads32 coefficient bytes instead of128.
The old identity-feedthrough t64/s4 count of70 is not a viable-distance
reference: that geometry is deterministically capped at12.5% distance.

## Schedule and limitations

Gather each eight-record tile once, compute its inverse transform, accumulate
its high coefficients, add the unchanged Q, and run the shifted forward
transform. Store the eight outputs once. Keep M*Q and the accumulated feedback
until the packet ends. This needs only bounded register temporaries; it
does not require a full-packet heap allocation, an independent C transpose
transform, or a separate state-emission transform.

The existing lane shuffles and GFNI dependency chains remain. Their compiled
counts, constant reuse, register pressure, and overlap with routed loads are
unmeasured. Wider packets reduce refresh boundaries, but do not change the
1.25 GFNI-products-per-output count in this family. Four old refreshes with
intervening inputs are not equivalent to one new refresh.

The identity-feedthrough quiet-subcode argument no longer preserves output
support: a nonzero input coordinate now spreads through its eight-record
tile. This removes that argument's premise, not every possible obstruction.
The many-row proof must retain the new feedthrough, state basis, A/C relation,
and refresh frequency. Routing and outer work are unchanged, so the9.1%
product reduction does not imply9.1% less full-encoder time or a result below5 ms.

Next test the exact candidate against structural cuts and derive its joint
output/feedback law. Only a surviving map should receive native basis tests
and a complete forward timing. No header, executable, benchmark, or frozen
checkpoint was changed for this note.
