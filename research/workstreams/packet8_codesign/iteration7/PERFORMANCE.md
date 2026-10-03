# K16--K20 byte-packet performance

Measured 2026-10-03 UTC on Peach, Ryzen 9 7950X, CPU15, GCC15.2 Release,
AVX-512/VBMI/GFNI, tuned for znver4. Calls encode 128-bit XOR payloads in
place, with precomputed setup and preallocated storage. Five warmups precede
each timed run. Every benchmark was serialized under the shared locks.
Correctness, setup, and allocation are outside the clock; required store
fences are inside it. No production or frozen source was changed.

## Matched holdout

Each result is the median of eight run medians, using four fresh seeds
(271,419,557,863), two runs per seed, with alternating ABBA/BAAB order.
The K16/K18 campaigns compare retained and copied flat9 cached circuits,
using 1001 calls per run. K20 compares byte packets against the archived
four-bit s20 kernel, using 301 calls per run. These codes have different
inner/route distributions; equal seeds do not imply equal sampled maps.

| K | Kernel | Page request | Median, us | Run-median range, us |
|---|---|---|---:|---:|
| 2^16 | Retained flat9, cached | Normal | **94.2260** | 93.795--98.594 |
| 2^16 | Copied flat9, cached | Normal | 96.1495 | 94.035--98.303 |
| 2^16 | Retained flat9, cached | Huge | 92.5975 | 91.992--96.490 |
| 2^16 | Copied flat9, cached | Huge | 92.3770 | 92.172--96.720 |
| 2^18 | Retained flat9, cached | Normal | 390.8935 | 386.420--409.894 |
| 2^18 | Copied flat9, cached | Normal | 387.0365 | 383.946--406.868 |
| 2^18 | Retained flat9, cached | Huge | **373.2610** | 372.675--411.097 |
| 2^18 | Copied flat9, cached | Huge | 373.6770 | 372.444--378.866 |
| 2^20 | Byte packet, streaming route/output | Huge | **3388.7655** | 3345.334--3563.291 |
| 2^20 | Prior four-bit s20, streaming route/output | Huge | **3267.5895** | 3231.823--3303.136 |

The byte K20 kernel is **3.71% slower**, not a speedup over the prior code.
K16's small page-policy difference does not justify changing that frozen
checkpoint. The copied cached circuit is not consistently faster than the
retained circuit; keep the retained one as the small-size control.

All sorted run medians (us), including high observations:

```text
K16 normal retained: 93.795 94.036 94.126 94.156 94.296 94.406 98.103 98.594
K16 normal copied:   94.035 94.156 94.266 94.607 97.692 97.802 97.812 98.303
K16 huge retained:   91.992 92.122 92.242 92.412 92.783 92.883 96.310 96.490
K16 huge copied:     92.172 92.182 92.312 92.372 92.382 92.433 92.583 96.720
K18 normal retained: 386.420 387.853 388.995 389.496 392.291 392.491 405.445 409.894
K18 normal copied:   383.946 384.908 385.679 386.731 387.342 389.566 391.089 406.868
K18 huge retained:   372.675 372.814 372.995 373.176 373.346 374.338 377.493 411.097
K18 huge copied:     372.444 373.185 373.185 373.226 374.128 378.205 378.605 378.866
K20 byte: 3345.334 3357.888 3367.095 3387.914 3389.617 3393.975 3408.573 3563.291
K20 old:  3231.823 3234.288 3235.069 3267.229 3267.950 3272.209 3294.751 3303.136
```

## Store-policy screen

These are exploratory run medians at seed41, not fresh holdout aggregates.
K16/K18 use normal pages and 1001/501 calls; K20 uses 201 calls.

| Policy | K16, us | K18, us | K20 normal, us | K20 huge, us |
|---|---:|---:|---:|---:|
| Cached route/output | 98.384 | 411.538 | 5613.872 | 5101.405 |
| Streaming route only | 193.992 | 747.937 | 4178.161 | 3499.133 |
| Streaming output only | 154.138 | 611.892 | 4772.191 | 4621.799 |
| Streaming route/output | 233.677 | 907.875 | 3836.413 | 3406.901 |

Reverse-order seed113 confirms K20's ordering: both-streaming 3374.477us,
route-only3491.806us, output-only4611.052us, cached5077.371us. Streaming
stores are harmful while the smaller working sets fit in cache; at K20
the measured store/page policies reduce the initial full-call latency by
about 39%. Byte packets alone do not halve payload traffic.

Same-call K20 phase diagnostics at seed113, huge-page request:

| Policy | Full call, us | Inner plus route, us | Outer, us |
|---|---:|---:|---:|
| Cached | 5084.674 | 3311.220 | 1769.748 |
| Both streaming | 3337.018 | 1764.668 | 1570.367 |

Phase medians need not sum to the full-call median. Headline comparisons
disable phase clocks. Both remaining stages are substantial; the data do
not support attributing all cost to field multiplication.

Normal-policy processes report zero anonymous huge-page KiB. Huge-policy
processes report 6144KiB at K16, 18432KiB at K18, and 67584KiB at K20.
These are process totals, not allocation-by-allocation residency proofs.

## Verification and identities

All ten native correctness configurations pass, including natural K=6144,
K18, K20, literal scalar/adjoint checks, all modes, misalignment fallback,
in-place encoding, guards, and padding. The old four-bit K20 binary also
passes its literal/adjoint/in-place checks at seed11 before the comparison.
An independent source review checked geometry, exact-map store changes,
completion fences, and test coverage.

```text
new binary SHA256:
a7b2e783b9e0907f12eb08aab2c36a8f047d88ec14758791a21c3e71df1a98e1
archived four-bit K20 binary SHA256:
37c553278e15873d990aeb3009d9828b79880ee99b66ba8a5b3d33da58451dec
new experiment source archive SHA256:
e525c123dabf5b8a84b252a74247babb3c8dff8df2345a3b74fd6886e15b80af
raw result archive SHA256:
124ba615aec9db669a86cd4727e930c8c360b07d0934e468302734e8403c3354
```

The ignored archives are `tmp/packet8-scaling-v1-source.tar.gz` and
`tmp/packet8-scaling-results-v1.tar.gz`. The source archive contains the
new experiment before documentation/comparator additions; frozen source
dependencies are the unchanged iteration5/6 inputs. Native sources and
reproduction commands live in the
[experiment](../../../../spin/experiments/packet8_wider24_scaling/README.md).

This branch is parked by user decision. If resumed, the practical next step
would be completing K18's proof coverage, not promoting the byte K20 encoder.
The larger-size proof screens are recorded separately and must not be
described as complete certificates.
