# Four-row BCH locality screen

2026-09-19. Keep the existing globally shared four-row BCH circuit.
Five alternatives reduced sharing lifetimes or forced input reloading;
none improved complete K20 encoding. The default remains
`SPIN_BCH_SCHEDULE=global`. No new speedup is claimed.

The [K20 profile](K20_PROFILE.md) found substantial BCH stack traffic.
This screen tests whether avoiding that traffic is profitable after accounting
for recomputation, repeated input loads, and extra calls. The result narrows
the optimization hypothesis: fewer spills alone do not imply a faster circuit.

## Candidates and measurements

All candidates process four independent BCH rows with AVX-512 vectors.
The input layout, output order, inner, routing, and code parameters are unchanged.
The grouped candidates synthesize shared XORs independently within consecutive
groups of 8, 16, or 32 outputs. Separate non-inlined functions bound sharing
across groups. The single-output candidate expands the existing circuit one
output at a time, discarding all sharing between outputs.

The reload candidate keeps the original shared XOR circuit. It replaces
retained input-leaf vectors with explicit reads of immutable input data.
A volatile unaligned SIMD read prevents the compiler from merging those
reads across outputs. This research variant uses the GNU SIMD type's
unaligned, may-alias semantics; it does not insert hardware fences.

Measurements use Peach's Ryzen 7950X, CPU 15, GCC 15.2, Release, znver4
tuning, and the same AVX-512 backend as the control. The workload is the
complete in-place transpose at K=2^20 on 128 parallel binary instances.
Setup, workspace allocation, initialization, and output hashing are excluded.
Each process performs three warmups and 51 timed calls without input reset.
Route seed is 1 and coefficient seed is 2.

Each entry is the median of two process medians. Candidate order reverses
in the second repetition, and all benchmarks run serially under both locks.
These are bounded selection measurements, not confirmed performance gains.

| BCH schedule | Transpose-only | Bidirectional transpose |
|---|---:|---:|
| Global control, grouped screen | 9.1790 ms | 9.0345 ms |
| Local groups of 32 outputs | 9.4243 ms | 9.3796 ms |
| Local groups of 16 outputs | 9.6545 ms | 9.5060 ms |
| Local groups of 8 outputs | 9.7866 ms | 9.5829 ms |
| One output per function | 13.2369 ms | 13.1136 ms |
| Global control, reload screen | 9.1963 ms | 9.0358 ms |
| Shared circuit with input reloading | 10.9608 ms | 10.7315 ms |

Even the closest candidate, groups of 32, is slower in both libraries and
both process repetitions. No candidate advanced to a larger confirmation
campaign or became the default.

## What changed in the generated code?

The assembly counts below cover the complete standalone BCH object,
including its wrappers. They count static instruction lines, not dynamic
executions. Stack references include memory-operand arithmetic as well as
explicit loads and stores; they are not all independently removable spills.

| Schedule | Static instructions | Instructions referencing the stack |
|---|---:|---:|
| Global | 6,317 | 4,567 |
| Local 32 | 6,969 | 4,493 |
| Local 16 | 7,623 | 4,299 |
| Local 8 | 7,507 | 3,075 |
| Single output | 9,207 | 0 |
| Input reloading | 7,384 | 3,938 |

The single-output version removes stack references but increases complete
bidirectional latency by about 45%. It repeats input loads and loses shared
work between outputs. Conversely, local-32 uses 3,388 symbolic XORs versus
3,758 in the existing circuit, but emits more instructions and runs slower.
Its source loads 640 input vectors across groups instead of retaining 256
input vectors for the full circuit. Source load counts are not dynamic
hardware load counts.

The existing sharing policy therefore represents a useful measured tradeoff,
despite its unattractive spill profile. The profile identified work worth
investigating, not an automatically recoverable performance loss. A better
schedule might still exist, but this screen supplies no evidence that more
group-size tuning would deliver a substantial gain.

## Correctness and proof scope

The generator checks each grouped output as a formal binary linear form
against the original BCH matrix. The reload transformation preserves the
verified XOR dependency graph and substitutes reads of the same input leaf.
All five screen builds passed the existing six full-encoder tests. The later
reload build passed seven tests, including the new direct BCH basis test.

The new test also ran against all six actual release libraries: control and
five candidates. It exhausts 256 input coordinates separately in each of
four SIMD lanes, checks 12 dense inputs, and checks input preservation and
output canaries. Every test passed. Complete benchmark output hashes agree
across candidates and libraries for matching iteration counts and seeds.

The full-encoder tests include K16, K18, and K20, both routing layouts,
alternative tiles, compaction, in-place suffix preservation, dense forward
and transpose references, and the adjoint identity. The rejected candidates
did not receive a new sanitizer campaign.

No construction parameter changes. K20 retains rate 1/2, 10% relative
distance, and the existing full 50.062-bit margin. The paper, certificate
producers, upstream native sources, and submitted artifact remain unchanged.

## Reproduction

`bch_local.py` generates experimental four-row kernels without changing
the pinned source. `SPIN_BCH_SCHEDULE` accepts `global`, `local8`, `local16`,
`local32`, `single`, or `reload`; only `global` is recommended.
The bidirectional generator copies the same selected BCH source.
Normal builds still default to AVX-512 disabled unless explicitly enabled.

Run from the repository root, serially:

```sh
bash workstreams/spin_optimized/bch_local_screen.sh "$PWD"
bash workstreams/spin_optimized/bch_reload_screen.sh "$PWD"
bash workstreams/spin_optimized/bch_basis_check.sh "$PWD"
python3 workstreams/spin_optimized/bch_local_summary.py \
  workstreams/spin_optimized/measurements/bch-local-screen
python3 workstreams/spin_optimized/bch_local_summary.py \
  workstreams/spin_optimized/measurements/bch-reload-screen
```

The scripts use the existing snapshot at `ROOT/hypercat/native/spin`.
`bch_basis_check.sh` links its small test against the already measured
libraries, without rebuilding their kernels. Future AVX-512 CMake test
builds include this test automatically and skip it when the backend is
unavailable. The ordinary caller remains compiled at the baseline ISA.

Samples, logs, assembly, and basis-test results remain ignored in
`measurements/bch-local-screen`, `measurements/bch-reload-screen`, and
`measurements/bch-basis`. No data files were committed.

## Recommendation

Stop this K20 micro-tuning branch for now. Retain its approximately 9.1 ms
kernel and the successful K16/K18 specializations. The next useful step is
to consolidate those implementations and measure their downstream benefit,
rather than repeat tile, prefetch, or output-group sweeps without a new idea.
