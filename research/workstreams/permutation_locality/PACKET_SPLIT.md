# Independent BCH-row shuffles with cache-line packets

## Result

Independently shuffling each BCH row, then forming packets from four adjacent
rows, preserves the fast global route. The best tested transposed encoder
takes **6.50 ms**, versus **5.42 ms** with a coordinate permutation shared by
the four rows: an additional **1.08 ms (19.8%)**. The difference is local BCH
input rearrangement, not long-range routing. This is a viable performance
candidate for the next proof investigation, not a completed distance proof.

All measurements below use K = 2^20, N = 2^21, BCH[256,128], IMT(128,19)
with two transvection updates, and 128-bit elements. They measure precomputed,
single-threaded, in-place transposed encoding, excluding setup and allocation.
Production code and paper numbers are unchanged.

## What changes in the permutation

There are 8192 BCH rows, partitioned into 2048 groups of four adjacent rows.
For each row, independently sample a uniform permutation of its 256
coordinates. At each shuffled position, collect one coordinate from each of
the four rows into a packet. The existing region shuffle moves these packets
as units; a lane shuffle randomizes the four positions within each packet.
Each row still contributes exactly one coordinate to each of the 256 regions.

Previously, the four rows shared their coordinate permutation. Equal BCH
words in a group therefore had equal active-region sets. Independent row
permutations remove this forced alignment: for fixed row weights a and b,
the expected number of regions occupied by both is ab/256. Under a shared
permutation, identical words of weight w overlap in all w active regions.
This observation motivates the change but does not prove the required
distance. Dense words, including all-one words, can still align. Nor does
the change separate the four coordinates after their packet has been formed.

For a fixed message, the four active-region sets are now independent uniform
subsets of their respective sizes. Positions within each row remain sampled
without replacement. A proof must retain these dependencies and the later
packet routing; replacing everything by independent Bernoulli bits would
analyze a different ensemble.

## Where to split the packet

In the forward direction, the proposed order is BCH encoding, independent
row-coordinate shuffles, packet formation, then the packet route and inner.
Transposed encoding reverses this: inner transpose, packet routing, local
inverse row shuffles, then BCH transpose.

The implementation keeps each long-range transfer at 64 bytes (four 128-bit
elements). It stores all 256 packets of a four-row group in a padded tile,
using row zero's original coordinate as the physical packet index. Row zero
is consequently already in canonical order. A precomputed table restores
the other three rows' coordinate order before the unchanged four-row GFNI
BCH circuit.

The winning variant first copies the 16 KiB tile sequentially into aligned
local storage, then gathers within that copy to form canonical BCH columns.
The gather pattern is local instead of a sequence of cold, scattered loads.
The source tile stride remains 1028 elements, including cache-set padding.
The extra row-offset table occupies 4 MiB at this K (one uint16 per outer
coordinate). Existing packet bases and lane controls occupy 6 MiB. These
are metadata sizes, not the total memory footprint of the benchmark harness.

## Measurements

Confirmation used two route seeds (1 and 17), three runs per seed, and 101
timed encodings per run after three warmups. Run order was reversed in the
middle repetition. The table reports the mean and range of the six run
medians, not a confidence interval. Inner masks were fixed across variants.

| Layout | Mean of run medians | Range | Change from baseline |
|---|---:|---:|---:|
| Shared row permutation, existing packet route | 5.420 ms | 5.413–5.429 ms | — |
| Independent rows, copy tile then repack | **6.496 ms** | **6.471–6.521 ms** | **+19.8%** |
| Independent rows, copy tile then fuse gathers into GFNI | 6.601 ms | 6.578–6.609 ms | +21.8% |

Phase means corroborate the location of the overhead: the baseline spends
about 2.72 ms in inner-plus-route and 2.71 ms in BCH. The winning independent
variant spends about 2.74 ms in inner-plus-route and 3.76 ms in local
rearrangement-plus-BCH. Phase means include warmups and need not sum exactly
to the full-encoder median.

Earlier two-seed screens used 31 timed encodings per run:

| Alternative | Full encoder | Interpretation |
|---|---:|---|
| Direct 128-bit global scatter into canonical BCH coordinates | 14.76–14.79 ms | Breaking cache-line stores is expensive |
| Packet route, then gather/repack directly from routed tiles | 8.28–8.30 ms | Local indexing alone does not eliminate cold-load cost |
| Packet route, direct gathers fused into GFNI | 8.13–8.14 ms | Fusion alone is insufficient |
| Prefetch routed tile before repacking | 8.19–8.20 ms | Tested prefetch loop does not recover the copy benefit |

The copy's benefit is measured; attributing it to better cache residency and
more predictable loading is an implementation-level explanation, not a
hardware-counter diagnosis. No setup-time claim follows from these tests:
the new ensemble samples four coordinate permutations per group instead of
one, and builds the additional offset table before timing.

## Checks and reproduction

`experiment.cpp` checks route bijectivity, region occupancy, four-row packet
membership, and offset completeness. Each tested mode compares the complete
output against an explicit dense two-update inner, a materialized route, and
the original BCH transpose; it also checks the inner adjoint identity and
unchanged output suffix. Both route seeds passed these checks at K = 2^20.
The GFNI generator independently checks 65,536 submatrix actions and 512
transpose basis vectors. These are implementation checks, not distance tests.

All six independent-row modes also passed AddressSanitizer and
UndefinedBehaviorSanitizer at K = 2^14 for both route seeds (12 cases), with
the same full-output checks. Sanitizer timings are not performance results.
The build used `-O1 -fsanitize=address,undefined -fno-omit-frame-pointer`.
The four existing SPIN library tests passed in both Release and sanitizer
builds as well.

On the existing Linux experiment build, run sequentially:

```sh
bash research/workstreams/permutation_locality/run_packet_split.sh ROOT confirm
cmake --build ROOT/build-sanitize -j 4
bash research/workstreams/permutation_locality/run_packet_split.sh ROOT sanitize
```

Here ROOT contains the experiment's `build` and `build-sanitize` directories.
The script pins timing runs to CPU 15 and acquires the three shared benchmark
locks. Do not run competing benchmarks or rebuild during measurements.
The measurement machine was Peach (Ryzen 7950X), GCC 15.2, Release build tuned
for znver4. Input is reinitialized before warmups but not between timed
in-place calls; setup, correctness checks, and allocation are outside timing.

Remote root: `/tmp/spin-locality-JjY7yR`.
Confirmation logs: `measurements/packet-split-confirm-a2l3GT`.
Exploratory logs: `packet-split-screen-it6Y6L`, `packet-split-warm-ZPP7Ws`,
and `packet-split-final-uobnPP` under `measurements`.
Raw samples stay outside the tracked report. Measured executable SHA-256:
`c1a215d6f527c731133adaacdf719bb3584a1e64b5f4debffb79cf2fc23e626a`.
Linked SPIN archive SHA-256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Next step

Analyze independent row support sets before further tuning. In particular,
test whether averaging their intersections removes the current four-row
pattern-count obstruction. Keep all-one and dense-row cases explicit. The
shared-permutation q = 1 through 14 certificates do not automatically apply
to this changed ensemble, and neither ensemble currently has a full-code
certificate. If the new proof becomes tractable, 6.5 ms is the measured
starting point for closing the remaining performance gap toward 5 ms.
