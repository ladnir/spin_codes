# Wider-outer, 24-bit encoder: K16 performance

The isolated precomputed transpose takes **98.1625 microseconds** at
K=65536 in the fresh-seed holdout, measured as the median of eight run
medians. This meets the 100-microsecond median target on this machine.
Seven run medians were below 100 microseconds; one was 102.011 microseconds.
No sample or run was discarded.

The matching [outward certificate](README.md) gives 68.893748662558 bits
of whole-code margin at the 10% distance target. Its exact union has also
passed an independent global replay. Timing and proof remain separate checks
of this same construction.
The construction has not been promoted into the production library.

## Workload and measurement

The measured operation is the complete precomputed transpose on 128-bit XOR
elements: N=131072 input elements produce K=65536 output elements.
It includes inner evaluation, routing, outer randomization, RS transpose,
and final output unpacking. Setup and allocation are excluded.

The host is an AMD Ryzen 9 7950X, with the process pinned to logical CPU 15.
The build uses GCC 15 and the experiment's Release CMake configuration.
The userspace governor permits frequencies from 3.0 to 4.5 GHz;
the experiment did **not** establish a fixed 4.5 GHz clock.
All three implementations use normal pages.

The retained [campaign script](../../../../spin/experiments/packet8_wider24/run_ab.sh)
runs seeds 41, 113, 257, and 997. For each seed, it executes the three
implementations in A, B, C, C, B, A order. A is the certified nibble control,
B is the earlier byte-native control, and C is the wider-outer candidate.
Each run performs five warmups and 2001 timed calls.
The script holds the shared encoder benchmark locks throughout the campaign.
All runs are sequential. Headline timings have intermediate phase clocks
disabled.

The timed candidate executable has SHA256
`2c8192bc3ba85195bcd4eddd52103af2e531cb62fe79d73776e144abc364fb7e`.
The retained [raw holdout log](measurements/holdout-abc.txt) contains all
24 runs, executable hashes, per-run quantiles, and output checksums.
The [implementation audit](IMPLEMENTATION_AUDIT.md) records the source hashes
and the construction-to-kernel mapping.
All 19 local C++ source/header dependencies were also compared with the
remote timing build. Fourteen were byte-identical; five differed only in
LF versus CRLF line endings. The proof receipt pins the local file bytes.

## Complete holdout results

Each cell is one run's median latency in microseconds. Pass 1 is the forward
half of the A, B, C, C, B, A sequence; pass 2 is its reverse half.

| Seed | Pass | Certified nibble control | Earlier byte control | Wider outer, 24-bit state |
|---:|---:|---:|---:|---:|
| 41 | 1 | 99.867 | 92.433 | 97.933 |
| 41 | 2 | 106.409 | 92.553 | 97.442 |
| 113 | 1 | 99.316 | 90.920 | 99.846 |
| 113 | 2 | 99.045 | 91.851 | 98.383 |
| 257 | 1 | 99.666 | 90.829 | 102.011 |
| 257 | 2 | 99.516 | 91.200 | 97.733 |
| 997 | 1 | 106.038 | 91.300 | 98.453 |
| 997 | 2 | 99.446 | 91.501 | 97.942 |
| **Median of run medians** | | **99.5910** | **91.4005** | **98.1625** |

The candidate's summary latency is 1.43% lower than the certified nibble
control and 7.40% higher than the earlier byte-native control.
These percentages compare medians of run medians, not paired confidence
intervals or a pooled median over all calls.

The earlier byte control has a smaller outer and a 16-bit state. It does
not have the candidate's positive all-occupancy proof proposal. The new
implementation pays for the larger construction without crossing the
100-microsecond median target in this campaign.

## Diagnostic split

A separate seed-41 run enabled phase clocks. These measurements diagnose
the stages; they do not replace the headline campaign.

| Implementation | Total, microseconds | Inner and routing | Outer |
|---|---:|---:|---:|
| Earlier byte control | 91.301 | 44.864 | 46.416 |
| Wider outer, 24-bit state | 97.932 | 45.876 | 52.048 |

Stage medians need not sum to the total median. In this diagnostic, the
outer accounts for most of the added latency: approximately 5.63 microseconds,
compared with 1.01 microseconds for the inner and routing.

The generated `reverseRouteFast` object contains no stack-pointer references
and therefore no stack spills in that function. There is no evidence here
for changing the inner merely to shorten its register lifetimes. The outer
uses an explicit 8 KiB packed workspace, which is distinct from accidental
register spills.

## Correctness and next step

Native CTest passed all six combinations of K in {2048,6144,65536} and seed
in {1,17}. Tests compare the optimized stages against literal scalar maps,
check routing and memory guards, exercise in-place output, and verify the
complete forward/transpose identity independently in all 128 payload lanes.
The fresh-seed timing runs execute those correctness checks before timing.

The measured baseline and completed certificate should remain frozen.
A later optimization should be a separately measured variant,
motivated by stage or assembly evidence. No production promotion is implied
by the latency result alone.
