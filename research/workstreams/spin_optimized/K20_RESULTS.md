# K20 routing exploration

2026-09-19. Retain the existing tiled implementation at K=2^20.
This bounded screen found no improvement from direct routing, streaming
stores, buffered bucket writes, or the tested alternative tile sizes.
The successful K18 direct route should not be extended to K20 by default.

The workload remains the rate-1/2 BCH [256,128] encoder with the
weight-five IMT inner, (t,s,r)=(128,19,1), at 10% relative distance.
Its full certified margin remains **50.0620882264 bits**; see
[the certificate ledger](../inner_design/finite_migration/PAPER_RESULTS.md).
These experiments change storage and scheduling, not the binary linear map,
setup distribution, or code parameters. No new distance bound is claimed.

## Fresh confirmation of the retained default

Each entry is the median of three process medians, with three warmups and
101 measured calls per process. The machine, workload, and exclusions match
the screen described below. Library order reverses in the middle repetition.

| Route seed | Transpose-only | Bidirectional transpose |
|---|---:|---:|
| 1 | 9.119352 ms | 9.072505 ms |
| 17 | 9.093224 ms | 9.118912 ms |

The current implementation therefore remains approximately 9.1 ms.
The paper reports 10.110 ms on the same host in a separate campaign:
approximately 1.11 times the current latency. That gain comes from the
earlier BCH optimization, not a new improvement in this routing screen.

The fresh default build passed all six expanded tests. The object-code
audit confines EVEX instructions to the separately dispatched AVX-512
objects. Final matched output hashes agree across libraries and repetitions.

## Candidate screen

Peach: Ryzen 7950X, CPU 15, GCC 15.2, Release, znver4 tuning, AVX-512 BCH.
Times measure the complete in-place transpose on 128 parallel binary instances.
Setup, workspace allocation, input initialization, and hashing are excluded.
Each process performs three warmups and 51 timed calls without resetting input.
Every paired candidate uses route seed 1 and coefficient seed 2.
Benchmarks run serially under both project locks.

Each table entry is the median of two process medians, with candidate order
reversed in the second repetition. These are selection measurements, not
distribution-wide or cross-machine performance guarantees.

| Routing candidate | Transpose-only, ms | Bidirectional transpose, ms |
|---|---:|---:|
| Tiled control, direct-route screen | 9.2557 | 9.0584 |
| Direct 32-bit | 16.3117 | 16.2299 |
| Direct 32-bit, write-prefetch 16 | 15.1798 | 15.0245 |
| Direct 32-bit, write-prefetch 64 | 13.7321 | 13.6806 |
| Tiled control, streaming-store screen | 9.2444 | 9.0283 |
| Tiled, 16-byte streaming stores | 11.9402 | 11.6730 |
| Tiled control, buffered-store screen | 9.1804 | 9.0151 |
| Buffer four blocks per bucket, then streaming stores | 16.3434 | 16.3064 |
| Buffer four blocks per bucket, then cached stores | 13.5591 | 13.4755 |

All matched output hashes agree across candidates and libraries.
Each of the seven builds passed all six tests. The buffered-store builds
also passed expanded K20 tests: alternative tiles, both layouts, in-place
suffix preservation, dense forward/transpose references, adjoint equality,
packed-bit forward encoding, and unchanged wide-view schedules.
The rejected candidates have not received a new sanitizer campaign.

The tiled control also tried 128, 256, 512, 1024, and 4096 rows with both
index formats. Single-process medians ranged from 9.2809 to 10.4306 ms.
The default 2048-row geometry remains selected. The streaming-store variant
also tried 512, 1024, 2048, and 4096 rows; none beat the tiled control.

## What the experiments establish

At K18, direct routing scatters into 8 MiB. At K20, that destination grows
to 32 MiB, in addition to the input, routing metadata, and output.
Direct routing avoids an intermediate pass but loses badly in this screen.
The larger scattered working set is a plausible explanation, not a measured
cache-miss attribution: this campaign did not collect hardware counters.

The tiled path first writes four bucket streams, each 8 MiB, then scatters
one bucket into an 8 MiB tile before applying BCH. Its default workspace
remains 40 MiB. After Packed24 compaction, retained setup is 12,713,992 bytes
for transpose-only and 27,525,132 bytes for the bidirectional library.

Streaming stores preserve the same bucket addresses and fence before bucket
reads. The buffered variants use 256 bytes of aligned stack scratch, with
four blocks per bucket. They apply only to the default four-bucket geometry.
Each bucket's slots arrive in descending order, so all four pending blocks
are initialized before a group is written. Other geometries retain the
ordinary kernel. No variant allocates memory in the encoding loop.

Reducing the number of full-buffer passes is therefore not sufficient for
this workload. Neither software buffering nor bypassing the cache solved
the problem. This does not rule out a different exact routing schedule.

## Replay and scope

`SPIN_K20_EXPERIMENT=0` remains the default. Modes 1, 3, and 4 select direct
routing variants; mode 2 exposes the packed direct-route implementation but
was not part of this K20 screen. Modes 5, 6, and 7 select streaming stores,
buffered streaming stores, and buffered cached stores. They are research
overrides, not recommended deployment settings.

Run these scripts serially from the repository root, supplying its path:

```sh
bash workstreams/spin_optimized/k20_screen.sh "$PWD"
bash workstreams/spin_optimized/k20_stream.sh "$PWD"
bash workstreams/spin_optimized/k20_batch.sh "$PWD"
bash workstreams/spin_optimized/k20_check.sh "$PWD"
```

The scripts use the upstream snapshot at `ROOT/hypercat/native/spin` without
editing it. The last script rebuilds and checks the retained default, audits
ISA isolation, then measures three processes of 101 calls for each of two
route seeds and both libraries. Logs, samples, and source/binary receipts
remain ignored under `measurements/k20-*`.

K16 and K18 retain their existing optimized paths. The submitted paper,
artifact, and certificate producers are unchanged by this exploration.
No data files were committed.

## Next step

Use the improved K16/K18 paths and keep tiled routing at K20. Before another
K20 schedule search, profile the current fused kernel to distinguish inner
arithmetic, routing stalls, and BCH spill traffic. Earlier scheduling and
prefetch screens are recorded in
[the scheduling report](../inner_design/scheduling_20260919/README.md).
Repeating those variants without a new bottleneck measurement has low value.

The follow-up [phase and hardware-counter profile](K20_PROFILE.md) is now
complete. It supports investigating BCH register pressure next, while
retaining the existing tiled route and inner parameters.
