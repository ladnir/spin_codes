# Native timing of the certified two-diagonal inner

The running-sum construction from [the 45% certificate](GROUPED_45_CERTIFICATE.md)
now has a checked full encoder. At K=2^20 it takes about34.5ms on Peach,
not the desired3--5ms. The certificate remains >45% symbol distance with
at least40 bits of setup-failure margin. No construction shortcut was used.

## Measured construction

The full timing uses RS[256,128],8192 blocks, one uniform permutation per
output region, and state/batch size512. Each of the N=2^21 records is one
GF(2^128) symbol. Setup is precomputed and excluded from encoding time.
The recurrence is

    Y_t = X_t + diag(u_t) D diag(v_t) Q_t,
    Q_(t+1) = Q_t + Y_t.

D is the fixed Cauchy map between disjoint additive cosets. Its application
uses a shifted transpose FFT followed by the inverse unshifted transpose FFT.
All u_t and v_t entries are nonzero full-field elements. Raw diagonal keys
occupy64MiB; route labels occupy8MiB. The benchmark uses a deterministic
MT19937 fixture, not a production cryptographic sampler. The proof concerns
the ideal independent uniform setup distribution.

The implementation uses a quadratic tower to retain the existing GFNI outer
and small-coefficient FFT arithmetic without hot-path basis conversions.
Let J be the existing GF(2^32) tower in MediumField32.h and mu=0x20000000.
Then

    F = J[w,z],  w^2+w=mu,  z^2+z=mu*w.

Trace_J/GF2(mu)=1. Each quadratic is irreducible by trace transitivity.
The record stores coefficients in the basis[1,w,z,wz], with four J
coordinates. Full multiplication uses packed Karatsuba and GFNI, not four
independent GF32 products. The reused Element type supplies storage only;
its original sparse-polynomial multiplication is not used by this encoder.

Only the RS256/state512 profile has a full native timing here. State256
Cauchy maps pass algebra checks, but this file does not implement and time
the other certified RS512/state256 full encoder.

## Correctness checks

`native/sandwich_native.cpp --check` checks1024 full-field products against
an independent scalar oracle, including basis inputs and random keys.
It checks the trace, selected Cauchy columns, and D^2=I at state256/512.
Random vectors agree between the specialized and retained runtime FFT maps.

At K16, the inner inverse recovers arbitrary input. Direct gathers,
outlined transforms, batch staging, and materialized routing agree exactly.
A separate scalar dense-Cauchy computation checks two full512-symbol batches,
including the first nonzero persistent state. Optimized and UBSan builds pass.
At K20, every timed full-encoder path agrees with the materialized reference.
These checks validate implementations; they are not empirical distance tests.

## Serial measurements

Ryzen7950X on Peach, GCC15.2, `-O3 -std=c++20 -march=znver4`, CPU2.
The final column below is the median of31 samples after3 warmups.
All phases run serially. Shared-host noise remains possible. Phase timings
are standalone measurements and must not be added to predict a fused run.

| Operation, processing N=2^21 records | Milliseconds |
|---|---:|
| One diagonal, direct flat tower kernel |10.9464|
| One diagonal through batch workspace |12.2036|
| Fixed MDS map, retained runtime control |9.6954|
| Fixed MDS map, compile-time specialization |6.9037|
| Both diagonals plus MDS map |27.2439|
| Complete inner, contiguous input |27.2597|
| Outer only |3.2063|
| Outer plus standalone materialized route |6.6805|
| Inner with direct fused gathers |37.4602|
| Inner with batch-staged gathers |31.4369|
| Full encoder, direct fused gathers |38.2086|
| Full encoder, materialized route |34.4928|
| Full encoder, batch-staged route |34.7299|

The first checked control took43.0179ms end to end in a7-sample run.
Specializing fixed FFT constants gives a reproducible map improvement;
retaining general full-field multiplication for those constants was costly.
Specializing the two fixed tower-reduction coefficients gave no measurable
standalone-product improvement. Outlining the map did not fix routing cost.

Staging16 gathers at a time into a512-record batch decouples random reads
from the long arithmetic loop. It substantially reduces the fused-gather
penalty. It is approximately tied with a materialized route end to end,
while avoiding that route's separate32MiB buffer. The small7-/31-sample
ordering difference is not evidence that staging is faster than materializing.

A separate representation diagnostic applies the existing sparse-byte
CLMUL/GFNI full-field product to a flat stream:6.6481ms per diagonal,
versus10.9464ms for the tower product in the same run. These outputs use
different field coordinates. The diagnostic is not a compatible replacement
inside the current encoder, nor a projected full-encoder timing. Compatible
small-subfield maps or explicit conversion costs must be included first.

## Interpretation and next gate

The proof reduction saves reset transforms and avoids dense random matrices,
but those operation-count savings do not establish a fast implementation.
Two full-field products per output dominate this tower realization. The
fixed MDS map alone also exceeds the entire3--5ms target.

Do not promote this kernel. Keep it as a certified, checked timing control.
The next discriminating experiment should jointly measure a faster field
representation and the actual fixed MDS map, including every conversion.
If that combined cost remains too high, revisit the feedback distribution
to reduce full-field products, with a new proof for that distribution.
Replacing either random diagonal by subfield keys does not preserve the
point-density lemma.

## Reproduction and preserved sources

Compile from a scratch directory containing the existing flat native header
set, including the unchanged reset and optimized outer sources:

    g++ -O3 -std=c++20 -march=znver4 sandwich_native.cpp -o sandwich_native
    taskset -c 2 ./sandwich_native --check
    taskset -c 2 ./sandwich_native --bench 20 31
    g++ -O1 -g0 -std=c++20 -march=znver4 -fsanitize=undefined -fno-sanitize-recover=all sandwich_native.cpp -o sandwich_ubsan
    taskset -c 2 ./sandwich_ubsan --check

Builds ran under a2GiB virtual-memory limit,60/90-second build timeouts,
and45/90-second checking timeouts. Benchmark timeout was60 seconds.
Remote scratch had63GiB free; local disk retained over20GiB free.
No raw timing logs, setup tapes, or binaries were added to the repository.

SHA256 of the checked source:

    sandwich_native.cpp: 684f6e222ad8f9ea5987288f5e13699a1c4e35fac5de36d8950db455a9a6cb40
    reset_native.cpp: 95f77f80d4d38050069acc30ed7eda5fd7e08624c8705c656a0c363702df8dc6
    RsOuterOptimized.h: a1c15859700fbec68681baedbe0ea1caf076a99b5dcf927edb6708d68bde0829

The reset source is included with its unused main renamed, preserving that
control unchanged. GCC warns that the renamed function lacks a return;
neither checking nor benchmarking calls it.
