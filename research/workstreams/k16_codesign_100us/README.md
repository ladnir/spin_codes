# K16 construction / implementation / proof co-design

The [frozen research checkpoint](../../../spin/experiments/k16_codesign_100us/checkpoint/README.md)
pins the certified winner and its source dependencies. It also records
larger-size timing probes without extending the K16 certificate.

Target: precomputed transposed encoding at K=65,536 and rate 1/2, carrying
128-bit XOR elements, in at most **0.100 ms** on Peach. The same construction
must have a complete **10% distance / 40-bit setup-failure** certificate.
**Target met:** the final s15 kernel measures **0.098910 ms**, with a separate
holdout measurement of **0.098950 ms**, and a complete **62.04117-bit** bound.
This is a research branch; production dispatch and prior proofs are unchanged.

## Current best certified construction

The [15-bit paired-map construction](proof/PAIRED_S15_CONSTRUCTION.md) has a fresh,
complete **62.0411732686-bit** bound at output-weight cutoff 13,107. Its
minimum distance is therefore at least 13,108 outside the stated setup-failure
event. The bound covers every occupancy from 1 through 512.

Each 128-to-256-bit outer group uses four parallel GF(16) RS[16,8] rows.
Independent 16-bit symbol randomizers precede four-bit packet routing into
64 regions. The inner processes 64 bits per step with a 15-bit persistent
state. Its fixed expansion uses seven affine functions, three single quadratic
monomials, and five disjoint pairs. The feedback is the transpose of this
expansion. Independent uniform GL15 state randomizers are sampled at setup,
once per physical step.

Replacing uniform GL16 matrices by independent nonzero GF(2^16) scalar maps
preserves the entire fixed-message output distribution, including fixed
invertible changes of coordinates. The [transfer argument](proof/TRANSITIVE_FIELD16.md)
justifies the cheaper native-layout outer. The final 15-bit inner uses GL15,
not the field16 replacement. The [independent review](review/PAIRED15_REVIEW.md)
checks its map interface, all occupancy endpoints, source pins, and exact
final dyadic sum. The prior [s16 certificate](proof/DISJOINT_PAIR_CONSTRUCTION.md)
remains intact at 62.8401523533 bits.

## Measured progress

Each row reports the median of eight process medians: four setup seeds in
forward and reverse order, 2,001 calls each, after five warmups. Timings
exclude setup, allocation, scalar checks, and diagnostic stage clocks. Runs
are serial on Peach's Ryzen 7950X, CPU15, with normal pages and GCC15.2.
The final session used the userspace governor with a 4.5 GHz configured
maximum; instantaneous frequency readings are not a fixed-clock guarantee.

| Construction / implementation | Whole encoding | Complete margin |
|---|---:|---:|
| Retained smaller RS/s16, final matched control (mode 0) | 0.125089 ms | 49.7217 bits |
| Initial paired t64/s16 implementation (mode 12) | 0.110516 ms | 62.8402 bits |
| Paired s15, peeled loop (mode 50) | 0.099807 ms | 62.0412 bits |
| Paired s16, reduced feedback shuffles (mode 51) | 0.100012 ms | 62.8402 bits |
| **Paired s15, reduced feedback shuffles (mode 52)** | **0.098910 ms** | **62.0412 bits** |

The final comparison uses seeds 1, 17, 41, and 113. Mode 52's process medians
range from 98.454 to 101.419 microseconds. A separate holdout uses seeds
913, 1123, 65,537, and 7,771 in both orders. It gives **98.950 microseconds**,
with process medians ranging from **98.754 to 99.466 microseconds**. The
matched holdout baseline is 124.908 microseconds, a **20.8% latency reduction**.
These are latency measurements, not a worst-case runtime guarantee.

The retained production control measured 0.152574 ms in the earlier matched
session, and the first native-field s16 candidate measured 0.115471 ms.
Those earlier measurements are useful context, not the final matched controls.

The final kernel changes the construction and the data representation together.
It retains the cheaper small outer, uses native-layout field16 symbol maps,
and removes one paired quadratic state coordinate. Exact implementation
changes keep feedback in wide registers, change the state basis, use raw-input
feedback via CA=0, peel the boundary steps, and reduce feedback shuffles.
The approximate stage split is 47 microseconds for inner plus routing and
51 microseconds for the outer. Cached stores remain faster than non-temporal
stores for this repeated K16 workload.

The [isolated C++ harness](../../../spin/experiments/k16_codesign_100us/)
contains both retained controls and candidate kernels. Its 319 compiled
checks cover scalar equality, forward/transpose adjoints, in-place use,
four alignments, and output guards at three lengths and two seeds. Every
confirmation process also runs the independent scalar and adjoint checks.
All **28 independent algebra, implementation-interface, and receipt-review
tests pass**, including the final state permutation and folded feedback.
The measured executable SHA256 is
`5678ca1235eade2b86843ff422d450638f39d2b889ff529d465635fa89024d33`.
Raw timing logs and numeric receipts remain ignored.

## Continuing experiments

- Alternating identity and randomized inner steps could remove half the
  refresh cost. This changes the construction and needs a new full bound.
- The 15-bit candidate is now fully certified and measured. A separate
  14-bit pairing search remains an optional source of further headroom.
- A 10-bit bilinear inner improves the expansion minimum weight to 28,
  but the current middle-occupancy bounds fail. See [its screen](proof/bilinear10/).
- Pure monomial maps reached roughly 0.105 ms but did not close the middle
  occupancies. They are not a certified alternative.
- Larger steps, smaller outer blocks, and eight-bit packets have partial
  screens or unsuccessful bounds. No such bound failure proves the code bad.

The goal is complete. The recommended next step is controlled promotion of
mode 52 as a K16 profile, retaining its exact hot loop and adding regression
benchmarks before generalizing the interface. Production promotion has not
been performed in this workstream.

## Reproduce the complete s15 certificate

From the repository root, with Python and python-flint installed:

```text
python -B research/workstreams/k16_codesign_100us/proof/reproduce_s15_drop10.py --output tmp/k16-paired15-fresh.json
```

The output path must be unused. The full replay regenerates local profiles
and computes every occupancy contribution with 256-bit outward arithmetic;
it does not read a saved numerical receipt. The review tests accept
`SPIN_PAIRED15_RECEIPT` as a path to this fresh file; otherwise they use the
author-side default receipt. Run `test_paired15_interface.py` and
`test_paired15_receipt.py` in the review directory through unittest discovery.
