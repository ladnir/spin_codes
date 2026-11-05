# Native persistent-state gate

This checkpoint implements the scalar-route GF(2^32) profiles from
[Persistent MDS](PERSISTENT_MDS.md) and the three separately certified
[packet-route profiles](PACKET_ROUTE.md). Four field lanes share one code in
each 128-bit record. These are forward encoders; no transposed timing is
claimed. All setup and allocation are outside encoding measurements.

The current implementation is substantially slower than the preserved
3.798 ms GF256 reference. Its purpose is to expose the remaining cost of
a construction with a complete proof for its actual distribution.

## Multiplicative-coset control

`native/Persistent32.h` implements the exact spectral feedback, fixed
Vandermonde moment map, and diagonal persistent update. The field uses
the existing GF256-based GF32 tower. A bounded exact-order check selects
alpha=65538, a primitive element. This is not the incompatible binary
polynomial representation in the same field header.

The two spectral cosets require four b-point transforms and 2b products
by independent uniform GF32 scalars. The moment map uses interpolation,
an exact novel-to-monomial conversion, reversal, and sparse corrections.
The persistent update additionally multiplies by the fixed powers of alpha.

Setup uses reproducible MT19937 output for experiments. The certificates
concern independent ideal uniform coefficients and permutations, not a
proof that a particular benchmark seed produces a good code.

### Verification

Both profiles pass independent root-product/Lagrange feedback oracles,
direct Vandermonde moment maps, direct unrolled global-state identities,
fixed-field multiplication checks through those oracles, regional
bijections, the RS outer oracle, and a complete small free-terminal stream.
Fused routing agrees with a materialized route. Both also pass UBSan with
recovery disabled. The large-field kernels are distinct from the earlier
GF16 toy proof tests.

### Serial measurements

Host: Peach, AMD Ryzen 7950X, GCC15.2, CPU2 pinned. Build uses
`-O3 -std=c++20 -march=znver4`; three warmups, median of31 samples.
The memory cap was1 GiB per run, with90-second process timeouts. No
benchmarks ran concurrently. Each profile processes K=2^20 input records
and N=2^21 output records. Component measurements are not additive.

| Phase | B64, state64 | B128, state32 |
|---|---:|---:|
| Outer |3.851 ms|6.912 ms|
| Materialized route only |6.450 ms|6.338 ms|
| Inner, contiguous input |23.033 ms|15.928 ms|
| Inner, fused route |30.487 ms|24.615 ms|
| Full, first run |34.218 ms|29.630 ms|
| Full, second run |34.526 ms|31.338 ms|

The same-session spectral192 GF128 window control took42.436 ms. The
preserved3.798 ms reference was not rerun; it remains a historical target,
not a matched comparison from this session. Report the new full timings as
roughly34--35 ms and30--31 ms rather than selecting only the best sample.

Each new profile stores64 MiB of prepared feedback keys. Setup took about
47--48 ms in these runs. This is deliberately a precomputed encoding test.

Full output checksums are10680270870258502386 for B64/state64 and
9179641585574512085 for B128/state32. They identify reproducibility under
the benchmark's fixed seeds; they are not cryptographic hashes or certificates.

Reproduction with the source and its included headers on the include path:

```text
g++ -O3 -std=c++20 -march=znver4 -I. persistent_native.cpp -o persistent_native
./persistent_native --check
taskset -c 2 ./persistent_native --bench 64 64 20 31 all
taskset -c 2 ./persistent_native --bench 128 32 20 31 all
```

Run the commands serially. This default build accepts only the two
certified scalar-route profiles and bounded K/repetition ranges.
It writes no timing dumps.

## Additive-coset coordinates

The [additive-coset variant](ADDITIVE_MDS.md) removes the separate moment
map and diagonal state multiplication. The two interpolation maps used by
feedback also supply the state update. It preserves the exact-uniform
feedback and global MDS reset interfaces, so the scalar-route certificates
transfer: >35%/46 bits for B64/state64 and >35%/44 bits for B128/state32.

Both profiles pass direct Cauchy-matrix oracles for every scheduled shift
through field bit21, involution checks, RS and route checks, and independent
full-stream recurrences. After every checked batch, the stored state also
matches the conceptual global Vandermonde syndrome. UBSan passes with
recovery disabled. The checks use GF32 arithmetic, not only toy fields.

Matched serial measurements use the same host, compiler, affinity, warmups,
and sample count as above:

| Phase | B64, state64 | B128, state32 |
|---|---:|---:|
| Outer |3.824 ms|6.853 ms|
| Materialized route only |6.473 ms|6.353 ms|
| Inner, contiguous input |10.401 ms|7.810 ms|
| Inner, fused route |19.353 ms|17.865 ms|
| Full, first run |22.973 ms|25.507 ms|
| Full, second run |24.246 ms|25.369 ms|

The best scalar-route profile is therefore approximately23--24 ms, not
near the historical3.798 ms reference. Prepared keys still occupy64 MiB;
setup took roughly47--54 ms. Full output checksums are5265152156448912871
and6785076828369530904, respectively.

Compile with `-DSPIN_OR_ADDITIVE=1` and run the same check and benchmark
commands. Compiling without that flag retains the multiplicative control.
Neither distribution is silently replaced by small random coefficients.

## Whole-packet routing

The [packet-route construction](PACKET_ROUTE.md) moves fixed groups of two
or four adjacent outer rows together in each region. Its distance bounds
use the actual group distribution, not the scalar-route certificate.
The native implementation stores one label per packet and reads complete
aligned packets in the fused encoder. Width two uses two 256-bit loads
per four records; width four uses one 512-bit load. The default width-one
gather loops are preserved.

### Verification

All three packet profiles pass `--check`. The check uses an independent
scalar source-index formula for group bijections and lane order. Fused
encoding agrees with a materialized route. Direct checks cover the RS
outer, Cauchy shifts and involution, the complete K256 stream, and the
global Vandermonde state after every batch.

Additional K1024 fixtures ensure nontrivial packet permutations. This
matters for B128/p4, whose K256 regions contain only one packet.
K=B*b fixtures exercise two batches per region and the fast fused path.
The larger fixtures check routes and fused/materialized agreement without
repeating the quadratic direct field oracles. All three packet profiles
also pass their own UBSan checks with recovery disabled.

### Serial measurements

The host, compiler, CPU2 affinity, three warmups, and median-of-31 protocol
are the same as above. All runs are serial, at K=2^20 and N=2^21.
Setup and allocation remain outside encoding measurements.

| Phase | p2, B128, state64 | p4, B128, state128 | p4, B64, state128 |
|---|---:|---:|---:|
| Outer |6.847151 ms|6.798981 ms|3.728449 ms|
| Materialized route only |4.062343 ms|3.155800 ms|3.154939 ms|
| Inner, contiguous input |10.474301 ms|14.321332 ms|13.710633 ms|
| Inner, fused route |16.962221 ms|18.382182 ms|18.484083 ms|
| Full, first run |23.489125 ms|25.127033 ms|23.184977 ms|
| Full, repeat run |23.447340 ms|25.826732 ms|24.051388 ms|
| Setup, approximately |45.9 ms|54.2 ms|54.2 ms|
| Prepared keys |64 MiB|64 MiB|64 MiB|
| Route labels |4 MiB|2 MiB|2 MiB|

Component timings are not additive. The smaller route cost does not by
itself establish a faster encoder: the certified packet profiles also
change state width or outer size. A same-session scalar B64/state64 control
takes 23.325273 ms. The packet profiles remain in that performance range
or slower; the repeats do not establish an end-to-end packet-routing win.
Do not select a winner from the lowest median. The historical 3.798 ms
GF256 reference is still substantially faster.

Full output checksums, in table order, are 11346476762432146158,
8398293696577212920, and 4303464318288633739. These identify the fixed-seed
benchmark outputs; they are neither cryptographic hashes nor distance
certificates.

Source and binary SHA-256 receipts for these runs:

| Artifact | SHA-256 |
|---|---|
| `persistent_native.cpp` | `f6e0f3882c9a4fe0790824e5e9c4f7b7b5873373280db3584b394c5773870533` |
| `AdditivePersistent32.h` | `dd6412a144e1695bb0d53a4d420622c6c451271344938d4f772befc2f51cf331` |
| Packet-two binary | `15cf08ea4730674b41d43c78a5b373eaf58300be1c8576006b542b157a40f35d` |
| Packet-four binary | `8ab1aed7f1dece9f5ab8108eba2f95bde801a5bad2959373911695fe08f55359` |

With the source and included headers on the include path, reproduce serially:

```text
g++ -O3 -std=c++20 -march=znver4 -I. -DSPIN_OR_ADDITIVE=1 -DSPIN_OR_PACKET=2 persistent_native.cpp -o persistent_packet2
./persistent_packet2 --check
taskset -c 2 ./persistent_packet2 --bench 128 64 20 31 all
g++ -O3 -std=c++20 -march=znver4 -I. -DSPIN_OR_ADDITIVE=1 -DSPIN_OR_PACKET=4 persistent_native.cpp -o persistent_packet4
./persistent_packet4 --check
taskset -c 2 ./persistent_packet4 --bench 128 128 20 31 all
taskset -c 2 ./persistent_packet4 --bench 64 128 20 31 all
```

The packet-two build accepts only B128/state64. The packet-four build
accepts only B128/state128 and B64/state128. Both compile-time variants
retain the same bounded benchmark interface. To reproduce the sanitizer
gate, add `-fsanitize=undefined -fno-sanitize-recover=undefined` to each
build and rerun `--check`; do not use those builds for timing.

Next target the remaining inner arithmetic and fused-route overhead with
a structural change that preserves or replaces the actual proof interface.
These measurements do not justify promoting packet routing as a speedup.

## Full Cantor schedule

Compile-time `SPIN_OR_FULL_CANTOR=1` extends the existing GF256 Cantor
chain to32 directions in GF32. A bounded binary linear solve prepares
the chain; no solve occurs in encoding. The proof transfer and the exact
byte-field shortcut are specified in ADDITIVE_MDS.md, Section8.

All four exposed scalar/packet-four profiles pass their native oracles.
An additional independent tower-squaring oracle verifies all32 chain
relations, and binary elimination checks rank32. The complete stream
and global syndrome checks use the new deterministic directions. Changing
those directions changes benchmark checksums, not the uniformity or MDS
properties used by the certificates.

All four also pass UBSan with recovery disabled. The packet-four build
with debug information exceeded the2 GiB compiler cap; rebuilding with
`-g0` completed within that cap and passed. This did not change the timed
optimized binaries or disable sanitizer checks.

The first component runs, under the same serial31-sample protocol:

| Phase | p1/B64/b64 | p4/B64/b128 | p1/B128/b32 |
|---|---:|---:|---:|
| Outer |4.881 ms|4.798 ms|6.770 ms|
| Route only |6.543 ms|3.193 ms|6.413 ms|
| Contiguous inner |10.144 ms|14.537 ms|7.353 ms|
| Routed inner |19.223 ms|18.006 ms|15.504 ms|
| Full |24.685 ms|22.722 ms|21.648 ms|

Repeated serial raw/full/full/raw comparisons give:

| Profile | Raw control medians | Full-chain medians |
|---|---|---|
| p1/B128/b32 |22.727,25.470 ms|21.827,21.513 ms|
| p4/B64/b128 |23.387,23.534 ms|21.972,22.265 ms|

Rebuilding both controls from the same final source, changing only the
full-chain flag, then repeating the same serial comparison gives:

| Profile | Raw control medians | Full-chain medians |
|---|---|---|
| p1/B128/b32 |22.214,22.264 ms|21.441,21.687 ms|
| p4/B64/b128 |23.072,23.781 ms|22.426,22.074 ms|

These rebuilt controls support a modest improvement, not the larger
gain one might infer from comparing only the earlier slow raw runs.

The full-chain B128/b32 profile is about21.5--21.8 ms in these runs and
retains >35%/44 bits. Full-chain packet4 B64/b128 is about22.0--22.7 ms
and retains >35%/59 bits. The scalar B64/b64 run did not improve; the
flag remains optional, and no universal speedup is claimed. Variability
in outer/control times also argues against attributing every difference
to arithmetic counts. All results remain far above3.798 ms.

Full checksums in the first table's order are8385321041039357848,
3585798834595335514, and10103151014887362579. Setup remains precomputed;
keys still occupy64 MiB. Add `-DSPIN_OR_FULL_CANTOR=1` to the corresponding
additive build command and rerun `--check` before timing it.

Final source/binary SHA-256 receipts:

| Artifact | SHA-256 |
|---|---|
| Runner | `2ac22266d22722ad9537f5ad9064c7fbecabe507732d415f58b93824a138ec4b` |
| Additive header | `71e9cf3cb9530d1546d38b0d154ca6abc459b4cd8178f3043ce8630246d3ac6c` |
| Full-chain scalar binary | `ae62ca16b086653315bc0129287317640c5172f6836a39014f122fe5b70f7c83` |
| Full-chain packet-four binary | `964086ac9f5b921f55c050f25303adfd5efdee30ea20b7843f6079cca3a1062e` |

The next exact-distribution A/B should test compact random keys. The
existing tower multiplier can consume raw32-bit scalars instead of16-byte
prepared matrices. This would reduce random-key storage from64 to16 MiB
but add arithmetic; no speedup is assumed. A separate one-diagonal
construction could remove more work but still lacks a whole-code proof.
The earlier [field-kernel comparison](../char2_field_spin/HIGH_DISTANCE_INITIAL_RESULTS.md)
already found raw tower multiplication slower in isolation: about1.44 ms
versus1.07 ms per2^20-record pass on the same host. Therefore permit one
bounded whole-encoder A/B to test the different cache pressure, not an
open-ended retuning of a kernel already known to lose in isolation.

## Compact-key gate: no reliable encoding win

`SPIN_OR_RAW_KEYS=1` stores each random GF32 coefficient in four bytes,
instead of its sixteen-byte prepared matrix. Fixed state coefficients
remain prepared. Setup consumes the same random draws in the same order;
the route and encoded codewords are unchanged. The flag defaults to zero
and is restricted to additive profiles.

All four scalar/packet-four profiles pass the original native checks in
both representations. An additional test compares raw, prepared, and
selected kernels with the independent scalar oracle on1024 binary-basis
products,128 dense records, and zero fixtures.

The following matched runs use the full Cantor schedule, the same compiler
and CPU2 protocol above, three warmups, and31 samples per median. The
first pass measures components and full encoding; the reverse pass repeats
full encoding. Every run is serial and excludes setup/allocation.

| Profile | Prepared full, first/repeat | Compact full, first/repeat |
|---|---|---|
| p1/B128/b32 |21.386,22.278 ms|24.712,24.968 ms|
| p4/B64/b128 |22.935,23.153 ms|24.314,22.802 ms|

Compact keys reduce random-key storage from64 to16 MiB. They slow the
scalar profile and give mixed packet results, not a reproducible encoding
improvement. Setup takes about25 versus47 ms for the scalar profile,
and33 versus55 ms for the packet profile; setup is not this goal's metric.
The first-pass contiguous inner costs are7.419 versus10.986 ms (scalar)
and14.201 versus16.535 ms (packet), prepared versus compact.
Checksums match exactly:10103151014887362579 and3585798834595335514.

Stop this implementation branch after the bounded comparison. Keep the
optional representation for reproducibility, but retain prepared keys as
the control. The next structural gate removes a feedback diagonal and a
state transform; it does not merely trade coefficient traffic for work.

## Experimental one-diagonal gate

`SPIN_OR_SINGLE_DIAGONAL=1` implements the recurrence in
[One-diagonal feedback](SINGLE_DIAGONAL.md). It retains the fixed RS outer,
full Cantor schedule, and unshuffled packet route. Each batch computes
V=DZ, Y=X+r.*V, and Z'=V+Y. The second transform and random diagonal are
absent at compile time. Prepared random keys occupy32 MiB rather than64.
The default remains the certified two-diagonal construction.

This changes the code distribution. None of the two-diagonal distance
certificates transfers. The benchmark and check output explicitly labels
this mode experimental and uncertified.

All four exposed scalar/packet-four profiles pass direct Cauchy and
scalar field oracles, fused/materialized routing checks, complete small
streams, and the global Vandermonde syndrome identity. They also pass
UBSan with recovery disabled, built using `-O1 -g0` within a2 GiB cap.
Only nonsanitized optimized binaries were timed.

| Phase | p1/B128/b32 | p4/B64/b128 | p1/B64/b64 |
|---|---:|---:|---:|
| Outer |7.019 ms|4.140 ms|4.385 ms|
| Route only |6.409 ms|3.223 ms|6.477 ms|
| Contiguous inner |4.356 ms|9.311 ms|5.938 ms|
| Routed inner |12.086 ms|12.075 ms|12.549 ms|
| Full |18.529 ms|16.132 ms|17.149 ms|
| Repeated full |18.162 ms|16.015 ms|not repeated|

These use K20, four GF32 lanes per128-bit record, three warmups,
31 samples, GCC15.2 `-O3 -std=c++20 -march=znver4`, and Peach CPU2.
Setup/allocation remain excluded; runs are serial and component times
are not additive. Full checksums in table order are2384260574386945032,
15732821993527327231, and13434656586125930495.

Rebuilding the packet control from the same final source and running
control/single/single/control gives 22.250/17.034/16.934/22.536 ms.
Thus the current one-diagonal range is about 16--17 ms; the matched final
comparison supports roughly a quarter less time, not a claim based only
on the fastest sample. The packet result remains far above the preserved
3.798 ms reference. It is an
implementation measurement, not a newly certified operating point.
The partial-support and delayed-rank proofs are indexed in
SINGLE_DIAGONAL.md; completing their whole-code count is the next gate.

Reproduce by adding `-DSPIN_OR_SINGLE_DIAGONAL=1` to the full-Cantor
additive scalar or packet-four build commands, running `--check`, then
the corresponding bounded `--bench` command. Keep the flag zero for
the certified control. SHA-256 receipts for this gate:

| Artifact | SHA-256 |
|---|---|
| Runner | `96bbecdccf4b56e9ab537cd3e70c9b0316c53d5b994e9197c67c90afe4da9545` |
| One-diagonal scalar binary | `6cfdafb533c3c1be1bcc5ba7d36daa2ece1fd6973b96b303797a01e14f59a228` |
| One-diagonal packet-four binary | `834b831f322c64f65bb09be0986e98b8ef86213bf444429a5d3f85782f6beb5d` |

## Outer stores: isolated gain, limited full-encoder benefit

The default-off `SPIN_OR_OUTER_TILE` and `SPIN_OR_OUTER_STREAM` options
change only the outer's memory schedule. The fixed RS map, input basis,
random draws, route, and encoded outputs remain identical. The default
four-row cached-store implementation is retained verbatim as the control.

A tile evaluates several four-row SIMD groups in bounded aligned scratch,
then finishes their adjacent cache lines in one output plane before moving
to the next plane. Streaming stores use an explicit fence before the route
reads the workspace. Tile sizes are compile-time constants, not a dynamic
dispatch in the hot loop. No extra global permutation is introduced.

Serial K20, B64/state128, packet-four runs use the same compiler, CPU2,
three warmups, and 31 samples as above. All timings exclude setup.

| Outer schedule | Outer alone | Certified two-diagonal full |
|---|---:|---:|
| Four rows, cached control | 4.140 ms | 22.696 ms |
| Sixteen rows, cached | 4.336 ms | 22.472 ms |
| Sixteen rows, streaming | 2.583 ms | 22.367 ms |
| Four rows, streaming | 4.266 ms | 23.796 ms |

The isolated outer improvement does not translate into a comparable full
gain. Streaming stores also change the cache state seen by the following
route, so component savings cannot simply be added to the full timing.
The small two-diagonal full differences do not establish a reliable win.

For the one-diagonal candidate, a control/tile/tile/control replay gives
17.247/16.473/16.544/17.070 ms. That initially suggests a modest local
improvement, but earlier untiled runs already reached 16 ms. A final
tile/control/tile replay gives 16.797/16.489/16.732 ms: the advantage does
not persist. Do not claim a reliable full-encoder speedup from this change.
Retain the option for later comparisons, with no default or production
change. The candidate still has no complete distance certificate.

The initial measurements used SFENCE after the stores. Final source uses
one MFENCE at outer completion to make store-before-load ordering explicit
for immediate read-back; the final replay above uses that source. This is
one barrier per complete outer, not per tile. Both fences are documented
in Intel's [optimization reference](https://cdrdv2-public.intel.com/671488/248966-Software-Optimization-Manual-V1-048.pdf).

All four tested schedule/build combinations pass the direct checks for
both B64 and B128 packet-four profiles. The sixteen-row streaming,
one-diagonal build also passes UBSan with recovery disabled before the
single fence-intrinsic change; final source passes the direct native checks.
The outer checksum is 1605201245610991117; full two-/one-diagonal checksums are
3585798834595335514 and 15732821993527327231. Thus the memory change does
not alter either sampled code. Builds and runs were bounded and serial;
no timing dumps were retained.

Reproduce with `-DSPIN_OR_OUTER_TILE=16 -DSPIN_OR_OUTER_STREAM=1`, adding
`-DSPIN_OR_SINGLE_DIAGONAL=1` only for the uncertified candidate. Run
`--check` before `--bench 64 128 20 31 full`. Do not select the four-row
streaming variant on the strength of a presumed non-temporal-store win.
