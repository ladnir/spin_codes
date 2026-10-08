# Characteristic-two field encoders

Research prototypes only. The production binary SPIN and packet profiles are
unchanged. [Research status](../../../research/workstreams/char2_field_spin/RESULTS.md)
separates complete finite certificates from unproved batched candidates.

## What is implemented

`spin_char2_hd_sumrank_forward` contains the current high-distance forward
encoders: K20, rate 1/2, 40% GF65536-symbol distance, and eight field lanes
per 128-bit record. Setup is precomputed and excluded from encoding time.
The [optimization ledger](../../../research/workstreams/char2_field_spin/HD_ENCODING_OPTIMIZATION.md)
records exact constructions, certificates, and frozen performance controls.
`bucket8-rank2-nt` and `bucket48-fused-nt` have independently replayed
certificates above 50.99 and 51.04 bits, respectively, at about 5.4--5.5 ms.
`bucket4448-fused-nt` has an independently replayed certificate above
56.65 bits. Its timing is similar, not sub-5 ms. The
[final checkpoint](../../../research/workstreams/char2_field_spin/HD_OPTIMIZED_CHECKPOINT.md)
records the selected designs and the limits of the exploration.
Run `--check` for the dense scalar, basis, alignment, and alias tests.
Run `--large-check` for full K20 comparisons against the independent dense
scalar inner for all three leading designs, including alignment and aliasing.
These correctness checks are outside benchmark clocks.
A full encoding benchmark is `--bench 501 full 20 bucket8-rank2-nt`.
Run benchmarks strictly one at a time. Optional `HD_SUMRANK_*` compile
controls preserve negative experiments; leave them disabled for the
documented default timings. No production code or API is changed.

`spin_char2_high_distance_wide` adds full precomputed transposed cost
controls for the [lower-expansion search](../../../research/workstreams/char2_field_spin/HIGH_DISTANCE_LOWER_EXPANSION.md).
The length-64 variants have independent scalar/adjoint checks and pass
Release plus ASAN/UBSAN. Their first full timings confirm that 8K raw
expansion is expensive; these are not new certified production profiles.
The new GF(2^16), GRS[256,64]/keep2 variant uses a pruned additive FFT;
`--check256` validates its outer basis and scalar/adjoint behavior.
`--bench 16 2 101 alias 256` measures repeated in-place calls, with no
restoration. The corresponding `full` mode is out of place; `phases` and
`alias-phases` add diagnostic clocks. Its first full time is about 15.6 ms.

`spin_char2_wide_diagonal16` is a research arithmetic comparison for the
[high-distance exploration](../../../research/workstreams/char2_field_spin/HIGH_DISTANCE_INITIAL_RESULTS.md),
not a complete encoder. It compares byte, genuine 16-bit, old tower 32-bit,
and new precomputed-matrix 32-bit multiplication on 128-bit records.
`--check` runs independent arithmetic, alias, alignment, and tail checks.
`--bench all 101 all auto 113` measures all controls and sizes; run only
one benchmark at a time. Release, CTest, and ASAN/UBSAN pass. Larger-field
distance screens remain separate from complete certificates.

`spin_char2_rate_speed4_profile` isolates per-region compute and routing
costs without changing their execution order. `phases` adds diagnostic
clocks; `full` and `fresh` measure complete calls without those clocks.
Its sixteen alternatives remain experimental:
[the phase audit](../../../research/workstreams/char2_field_spin/GF8_SPEED_PROFILE.md)
finds no material win over `rs4TransposeRouted16`. Release and ASAN/UBSAN
checks pass. Run all benchmark commands serially.

`spin_char2_rate_speed4_next` is the latest exact-map optimization of the
certified83.3333% profile. The
[new native report](../../../research/workstreams/char2_field_spin/GF8_SPEED_ROUTED.md)
selects `RateSpeed4RoutedPlan` and `rs4TransposeRouted16` in
`RateSpeed4RouteNext.h`: about3.2ms at K20, about20% less time than the
matched earlier entry point. `1048576 151 selected` benchmarks those two
entry points with the selected plan only; `selected-fresh` restores input
before every timed call. The other modes preserve29 A/B alternatives.
Independent small/large checks and ASAN/UBSAN pass. The construction,
certificate, and installed library are unchanged.

`spin_char2_rate_speed4` preserves the first speed-first checkpoint: GRS[64,40],
width4/state4, keep3, nominal83.3333% rate. The
[complete certificate](../../../research/workstreams/char2_field_spin/GF8_SPEED_DESIGN.md)
gives10% distance and57-bit K20 setup-failure margin. The
[native report](../../../research/workstreams/char2_field_spin/GF8_SPEED_NATIVE.md)
recommends `rs4TransposeCached`, about4.1ms at K20. It includes independent
scalar checks, aliases, natural sizes, and four/eight-packet A/Bs.
`1048576 101 fresh` additionally restores input outside timing before every
call. `spin_char2_rate_speed8` preserves the older-map optimization controls.
Neither changes the installed library.

`spin_char2_rate_climb8` implements the certified fixed-six, 84%, and 84.09%
K20 profiles alongside the unchanged 75% control. The
[native comparison](../../../research/workstreams/char2_field_spin/GF8_RATE_NATIVE.md)
recommends the 84% profile at 4.807--4.849 ms per full precomputed transposed
call on 128-bit records. Its certificate gives 10% symbol distance and
115-bit setup-failure margin. `--check` uses an independent scalar oracle;
`--large-check` checks K16/K18/K20 schedules and aliases. These smaller-size
checks do not certify distance. The selected K20 benchmark command is:

```text
spin_char2_rate_climb8 1048576 101 full auto page folded final-store
```

The optional selectors retain routing, column-stride, layout-folding, and
outer-store A/B controls. Run benchmarks strictly one at a time. The new
headers and driver are research-only; no installed profile changes.

`spin_char2_punctured8` is the preserved three-quarter research checkpoint:
nominal rate3/4,
GRS[128,84], full-Cauchy width8/state8, and seven retained packet outputs.
The [K20 certificate](../../../research/workstreams/char2_field_spin/GF8_PUNCTURED_WHOLE_CERTIFICATE.md)
gives10% symbol distance with116-bit setup-failure margin. The
[native report](../../../research/workstreams/char2_field_spin/GF8_PUNCTURED_NATIVE.md)
records independent scalar checks and serial timings. It is not an installed
profile. `--check` validates small natural sizes; `--large-check` validates
K16/K18/K20 schedules; `1048576 101 full` times precomputed full calls.
Only K20 is newly certified. All older controls below are preserved.

`spin_char2_full` encodes over native GF(2^128), with modulus
X^128 + X^7 + X^2 + X + 1. It composes a GRS[8,4] or GRS[16,8] outer
with two independent symbol permutations and two input-weighted accumulators.
Each scan coefficient is a polynomial bit pattern from 1 through 128.
The corresponding ideal ensemble has the certificates in
[SMALL_COEFFICIENTS.md](../../../research/workstreams/char2_field_spin/SMALL_COEFFICIENTS.md).

`spin_char2_batched` uses the GRS[8,4] outer, eight fixed-coordinate regions,
independent group permutations inside regions, and a two-symbol state.
It supports t=8 or t=16 symbols per step. Three independently sampled
small-coefficient shears mix the state. The available modes are:

| Mode | State refresh | Independent input diagonal |
|---|---|---|
| `post` | After output | No |
| `pre` | Before output | No |
| `post-d` | After output | Yes |
| `pre-d` | Before output | Yes |

These are distinct constructions, not interchangeable implementation settings.
None has a whole-code distance certificate. In particular, the exact GF256
uniform-SL2 proof model in the research notes is **not** the small-polynomial
three-shear distribution implemented here.

`spin_char2_subfield` implements that distinct GF256 model: a GRS[16,8]
outer (GRS[8,4] is optional), one region per outer coordinate, independent
nonzero input scalars, and an exact uniform-SL2 refresh before each t=8 or
t=16 emission. There are two state symbols. Its hot kernel uses GFNI byte
products, not general full-width multiplication. Its K16 proof currently
covers sparse and dense classes separately, not every message.

The prototypes accept complete field symbols and compute the field transpose,
not the transpose of the expanded binary matrix. They replace
the first K records in an N-record buffer and preserve the remaining records.
Scratch is supplied by the caller; encoding makes no heap allocations.
Setup, allocation, and page preparation are separate from encoding.

`spin_char2_one_route` preserves the GF256 construction but adds exact-map
execution schedules. The best K20 schedule reorders one region at a time and
uses non-temporal stores; scatter remains faster at K16/K18. The `global`
route option is a negative control, not a certified replacement for regional
routing: a small outer followed by a uniformly global route has an actual
late-activation distance obstruction.

`spin_char2_subfield16` is the new proved one-route candidate: GRS[16,8]
over GF16, independent regional permutations, nonzero GF16 input scalars,
and a width-16/four-symbol-state inner with uniform GL4 refreshes.
Its maps are exactly A_j=(1,j,j^2,j^3), C=A^T. The
[K16 proof](../../../research/workstreams/char2_field_spin/ONE_ROUTE_PROGRESS.md)
certifies 10% symbol distance and 42 bits of ideal-setup failure margin.
This is a research prototype, not a new installed binary profile.
`spin_char2_subfield_wide` changes only the inner to six monomial state
coordinates and uniform GL6 refreshes. Its separate K18 certificate gives
10% symbol distance and 43 bits of ideal-setup failure margin.
GF16 uses X^4+X+1 and embeds its generator as byte 0x5c in AES GF256;
the basis bytes are 01,5c,e0,50. All field arithmetic is bytewise GFNI,
not general multiplication in GF(2^32), GF(2^64), or GF(2^128).

`spin_char2_field_frontend` also checks and times small-coefficient arithmetic
for m=8,32,64,128. Optimized complete encoders currently use 128-bit records.
The subfield prototype additionally has full scalar forward and transpose
references for all four field widths.

### Subfield representation

GF256 uses the AES polynomial. An m-bit record represents a degree-below-m/8
polynomial in Z with GF256 coefficients; byte i is its Z^i coefficient.
The extension moduli have the following hexadecimal GF256 coefficients:

| Field | Modulus over GF256 |
|---|---|
| GF(2^32) | Z^4 + d4 Z^3 + d1 Z + 3b |
| GF(2^64) | Z^8 + 0f Z^7 + 48 Z + 57 |
| GF(2^128) | Z^16 + 41 Z^15 + a5 Z + 63 |

Exact Rabin tests verify irreducibility. Multiplication by a GF256 scalar
therefore acts bytewise while remaining linear over the full extension field.
The code has the same minimum symbol distance in all four fields. This basis
differs from the native polynomial basis used by `spin_char2_full` and
`spin_char2_batched`. No implicit conversion is performed or timed.

## Build and validate

On an x86 host with AVX512, GFNI, and VPCLMULQDQ:

```sh
cmake -S spin/experiments/char2_field_spin -B out/char2-field \
  -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_CHAR2_SCALAR_LOOKUP=ON -DSPIN_CHAR2_CLMUL_REDUCTION=ON
cmake --build out/char2-field -j4
ctest --test-dir out/char2-field --output-on-failure
```

Both options default to OFF. Scalar lookup replaces affine-table gathers
with scalar loads. CLMUL reduction replaces two GFNI operations with a tiny
overflow-polynomial product and a fixed reduction constant. It is not a
general 128-by-128-bit field multiplication.

The arithmetic frontend dispatches among available ISA kernels. Complete
prototype executables require AVX512/GFNI and exit with code 77 when unsupported.
The CMake tests treat that code as skipped. These experimental executables
are not installed or linked into the production library.

Checks compare direct polynomial evaluation with additive-FFT kernels,
independent scalar forward/transpose implementations, full-field linearity,
field adjoints, in-place operation, suffix preservation, and regional routes.
The batched check covers all four constructions, both outer kernels, all
three inner kernels, and sizes whose steps cross region boundaries.
The subfield check covers both outers and step widths, its SL2 sampler,
general extension-field adjoints and linearity, baseline/prefetch agreement,
and natural sizes including steps that cross regions.

An additional Linux ASAN/UBSAN build passed all seven test targets:

```sh
cmake -S spin/experiments/char2_field_spin -B out/char2-san \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS=-fsanitize=address,undefined \
  -DSPIN_CHAR2_SCALAR_LOOKUP=ON -DSPIN_CHAR2_CLMUL_REDUCTION=ON
cmake --build out/char2-san -j4
ctest --test-dir out/char2-san --output-on-failure
```

## Serial measurements

Never run two benchmarks concurrently. Example commands on Linux:

```sh
taskset -c 2 out/char2-field/spin_char2_full 8 65536 lookup 113 51
taskset -c 2 out/char2-field/spin_char2_full 16 1048576 lookup 113 51
taskset -c 2 out/char2-field/spin_char2_batched 65536 16 113 51 huge post wide prefetch
taskset -c 2 out/char2-field/spin_char2_batched 65536 16 113 51 huge pre-d wide base
taskset -c 2 out/char2-field/spin_char2_subfield 65536 16 16 113 51 huge full base
taskset -c 2 out/char2-field/spin_char2_subfield 1048576 16 16 7919 51 huge phases base
taskset -c 2 out/char2-field/spin_char2_one_route 1048576 16 16 7919 51 huge full regional region-nt-all
taskset -c 2 out/char2-field/spin_char2_subfield16 65536 113 51 huge full scatter
taskset -c 2 out/char2-field/spin_char2_subfield16 1048576 7919 51 huge full region-nt-all
taskset -c 2 out/char2-field/spin_char2_subfield16 65536 113 101 huge full scatter packed
taskset -c 2 out/char2-field/spin_char2_subfield_wide 262144 113 101 huge full scatter direct
taskset -c 2 out/char2-field/spin_char2_subfield_wide 262144 113 101 huge full scatter fft
```

The full encoder syntax is `OUTER K MODE [SEED] [CALLS]`.
Its modes are `lookup`, `prepared`, `gather`, and `phases`.
The first three preserve the same seeded map. `prepared` expands coefficient
tables in setup; `gather` uses inverse routing and additional sequential passes.
`phases` measures the two scans and the outer separately using `lookup`.

The batched syntax is
`K T [SEED] [CALLS] [huge|normal] [MODE] [wide|scalar|phases] [base|unroll|prefetch]`.
`wide` processes four outer groups together; `scalar` retains the single-group
SIMD reference. These two kernels preserve the same map. `phases` separates
the inner plus routing from the wide outer finish. The final selector chooses
the original inner loop, explicitly unrolled batches, or unrolled batches
with destination prefetch. Prefetch is not a universal win across constructions.

The subfield syntax is
`K T [OUTER=16] [SEED=113] [CALLS=51] [huge|normal] [full|phases] [base|prefetch]`.
The baseline is the default. `prefetch` requests routed destinations 128
records ahead without changing the code. It regressed the small sizes and
helped K20 by less than 1% in the retained matched run.

The GF16 syntax is
`K [SEED=113] [CALLS=51] [huge|normal] [full|phases] [scatter|region-nt|region-nt-all] [direct|fft|packed]`.
It requires positive K divisible by eight within the 32-bit routing limit.
Region-based schedules fall back to exact full-code scatter, without an
additional regional temporary, when a region does not align to a complete
inner step. These accepted
natural sizes are not automatically covered by the K16 certificate.
The optional FFT factorization preserves the same monomial maps and sampled
code. It saves GFNI multiplications but adds lane shuffles; the measured
full encode regresses about 10--14%, so `direct` remains the default.
The `packed` option reuses compact coefficient loads and byte shuffles.
It is the fastest four-state variant in the latest matched campaign,
without changing setup storage or the code's distribution.

The separate `spin_char2_subfield_wide` executable fixes the state at six
symbols and uses the same syntax except its last selector is `direct|fft`.
Both maps already use packed coefficients. This is the profile covered by
the [K18 certificate](../../../research/workstreams/char2_field_spin/GF16_SIZES.md),
not an automatic claim for every accepted natural size. The six-state FFT
measures about 1.014 ms at K18, improving 3--5.5% across matched campaigns.
Both pass independent basis-vector tests and scalar/SIMD adjoints.

The separate `spin_char2_toeplitz` executable uses independent systematic
Toeplitz[64,32] outers over GF16, one 64-region route, and a three-symbol
state. Its default `nod` profile omits input scalars, as justified by the
[complete K20 proof](../../../research/workstreams/char2_field_spin/K20_RANK_CONTINUATION.md).
The syntax is
`K [SEED=113] [CALLS=51] [huge|normal] [full|phases] [scatter|region-nt|region-nt-all] [kara|dense] [nod|diag]`.
It accepts positive K divisible by 32 within the same uint32 route limit;
unaligned regions fall back to scatter. `--check` covers scalar references,
recursive/dense agreement, samplers, all inner basis vectors, field adjoints,
alignments, guards, and aliasing. `--check-large` compares complete outputs
at K16/K18/K20. Both pass in Release and ASAN/UBSAN; the research CMake now
has eight test targets. The initial timings and tradeoffs are in
[RESULTS.md](../../../research/workstreams/char2_field_spin/RESULTS.md).

Each measurement performs five warmups and then reports median, p10, and p90
over the requested calls. The input is not reset between calls. Checksums
permit comparisons between exact implementation variants at equal call counts.
There is no allocation, setup, data fill, or copy in the timed region.
Payloads at least 16 MiB request transparent huge pages under `huge`;
the request is a preference, not a guarantee from the OS.

Natural sizes: the full prototype requires positive K divisible by half the
outer length. The batched prototype requires positive K divisible by four
and 2K divisible by t. The subfield prototype requires positive K divisible
by half its outer length and 2K divisible by t. All reject sizes beyond
their 32-bit route offsets.
These input checks are implementation constraints, not certified size ranges.

## Higher-rate GRS[256,176] checkpoint

`spin_char2_cauchy256` implements the one-route, width-sixteen/state-eight
map with a complete K20 certificate for rate 11/16 and 10% distance.
This is separate from the half-rate controls above and from production.
It accepts positive requested K within its uint32 routing limit, pads to
a multiple of 176*64, and returns exactly the requested transposed prefix.
The K20 certificate does not cover every accepted size.

Its syntax is
`[--compare] K [SEED=113] [CALLS=51] [scatter|region-nt|region-nt-all] [full|phases] [basis|dense|fft|dense-fft|dense-fft-fused|fft-inner-fused]`.
The default `basis` retains the original factorization. Select
`dense-fft-fused` for the best measured implementation. All six variants
compute the same map, with compile-time kernels and one dispatch per call.
`--compare` requires `SPIN_CHAR2_PUBLIC_ROOT` at configure time and `full`
mode; it alternates batches with the public PacketRs8T64 encoder using
the same payload address and initial bytes. Setup and restores are untimed.

`--check` and `--check-large` pass in Release and ASAN/UBSAN for all variants.
The former includes independent scalar interpolation, packet state/output
checks, parity-only checks, field adjoints, padding, alignment and aliases.
The latter compares full native outputs at K16/K18/K20.
[Native results and reproduction](../../../research/workstreams/char2_field_spin/FIELD256_NATIVE.md)
record a higher-rate proof but no performance win over binary SPIN.

## Smaller GRS[128,87] checkpoint

`spin_char2_cauchy128` keeps the dense width16/state8 inner and one regional
route, with a smaller outer and an explicit four-row-interleaved source
order. K20 has a complete 10% certificate with 67 bits. Natural positive sizes
pad to a multiple of 87*64. Separate [K16/K18 certificates](../../../research/workstreams/char2_field_spin/FIELD128_SIZE_CERTIFICATES.md)
now give 63/79 bits at the same 10% target and requested rate 2/3.
Other accepted sizes still need their own proof.
The encoder uses caller-owned scratch and preserves the in-place tail.

Syntax:
`[--compare|--variants] K [SEED=113] [CALLS=51] [scatter|region-nt|region-nt-all] [full|phases] [baseline|compact|batch|compact-batch|lane-fft|compact-lane-fft]`.
`--compare` rotates five-call batches among binary SPIN, the optimized
GRS256/176 reference, and GRS128/87 on the same payload address. It requires
the optional public library and full mode. `--variants` compares all six
exact-map implementations on one payload with matching checksums. It uses
the diagnostic allocator rather than the public comparison's page policy.
All benchmarks must run serially. The default remains `baseline`.

Release and ASAN/UBSAN pass `--check` and `--check-large`; the registered
CTest also passes for all six variants. Both FFT layouts regress; keep the
dense inner. `compact-batch` is the latest exact-map candidate, with uint16
regional indices and a uint32 fallback. Both route tables remain in setup
for experimental controls. Latest K20 matched launches vary from 5.73 to
7.46 ms; the variation is unresolved. See [validation and timings](../../../research/workstreams/char2_field_spin/FIELD128_NATIVE.md)
for all retained results, including the earlier 6.01--6.10 ms checkpoint.
No production profile changes or performance promotion.

The proof scripts and reproduction commands live in the linked research
directory. Raw measurement logs and generated witness receipts remain ignored
under `tmp/char2-field-spin/`.

## Cubic width-eight research control

`spin_char2_cubic8` implements the separately certified GRS[128,83] /
width8-state4 map in coefficient source coordinates. It uses the full
length128 additive FFT and repaired asymmetric cubic maps. The complete
10%/65-bit certificate covers K20 only. Other natural positive requested
sizes pad to a multiple of83*64 and are not implicitly certified.

Syntax: `K [seed=113] [calls=101] [full|phases] [all|gather|fused|paired]`.
The optional public root adds the binary library control. The benchmark
rotates five-call batches on shared input and field scratch buffers;
setup, allocation, and input restoration remain outside timing. Use
`--check` and `--large-check` for correctness, not timing.

Release and ASAN/UBSAN pass independent scalar/direct-outer, shortened
adjoint, alias/tail/guard, stride, and large-size tests. Both outer schedules
compute the same map. The [native receipt](../../../research/workstreams/char2_field_spin/GF8_CUBIC_NATIVE.md)
records the first packing regression and its 14--16% paired-schedule fix.
`paired` compares the two exact-map inner schedules using the same outer
winner and existing controls. It remains a research control; do not
replace production or old benchmarks with it.
