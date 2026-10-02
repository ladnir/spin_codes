# Precomputed SPIN performance update

The paper uses the optimized implementation's retained confirmation campaign,
not the newer standalone library's slower consolidation build. All numbers below
are complete encoding times for 128 parallel binary instances, with setup and
workspace allocation excluded. Input is preserved and output is separate.
No heuristic permutation bank or fresh-code refresh is used.

| K | Inner (t,s,rounds) | Forward ms | Transpose ms |
|---|---|---:|---:|
| 2^16 | (64,12,2) | 0.348080 | 0.339835 |
| 2^18 | (128,19,1) | 1.601299 | 1.459816 |
| 2^20 | (128,19,1) | 9.030879 | 9.289382 |

The host is a Ryzen 9 7950X, CPU 15, with boost disabled and GCC 15.2.
Each cell is the median of three process medians, each with three warmups and
101 timed calls. Buffers are initialized once, without input resets between
calls. Directions alternate order across repeats. Benchmarks run serially.
These are warm steady-state measurements, not cold-cache or protocol timings.

At K=2^20, the transpose process medians are 9.158006, 9.289382, and
9.320550 ms. Forward medians are 9.010111, 9.030879, and 9.035999 ms.
The separately measured in-place transpose is about 9.135 ms; it is not the
number used in the paper's separate-output comparison.

At K=2^18, the same campaign's in-place transpose median is 1.418608 ms
(rounded to 1.419 ms in text). This gives 262144 / 0.001418608 = about
184.8 million 128-bit output blocks/s. It is the rate of the linear-encoding
stage, not full OT throughput: tree expansion and output hashing remain excluded.
The table keeps the separate-output figure, 1.459816 ms. We did not find a
retained validated 1.3 ms complete-encoder result at exactly K=2^18.

## Implementation and proof scope

The measured implementation and generator settings are documented in
[FORWARD_SCHEDULE.md](../workstreams/spin_optimized/FORWARD_SCHEDULE.md).
The selected options include AVX-512 BCH, four-row forward evaluation,
direct routing at small lengths, and the `dfs` forward schedule. At K=2^20,
forward uses 256-row tiles and transpose uses 2048-row tiles. The baseline
compiler target is x86-64-v3 with Zen 4 tuning and carry-less multiplication;
AVX-512 BCH is separately compiled and dispatched.

The one-round finite theorem and its five-length certificate ladder are
unchanged. K=2^16 uses a separately certified two-round inner, described
informally in the implementation appendix rather than folded into that theorem.
Its full bound is 49.3275773868 bits at relative distance 0.10, covering all
occupancies 1 through 512. See the
[full result](../workstreams/inner_design/k16_parameters_20260919/FULL49_RESULT.md)
and [integration checks](../workstreams/spin_optimized/R2_RESULTS.md).
The checker authenticates that certificate and compares its exact inner-map
columns with the current implementation. Authentication is not an independent
proof of the certificate producer's inequalities.

## Retained evidence and table generation

`precomputed_results.py` authenticates the 18 separate-output confirmation files
under `research/workstreams/spin_optimized/measurements/forward_schedule`.
In sorted filename order it hashes each filename, a newline, and the raw bytes:

```
fe4afcc33fec10d1a802500394991ad6091b03fef646bbd0f4f62f68c540b087
```

It recomputes the medians from all 1818 samples and checks configurations,
sample counts, and consistent setup/workspace metadata. The separate R2 receipt
is pinned by SHA-256:

```
dd171397b8ba21860aea9566f2c155b6410947d195d76eea45747d1cfb355c15
```

From the repository root, with the retained local evidence available:

```sh
python -B research/paper/precomputed_results.py
python -B research/paper/build_imt_comparison.py --check
python -B research/paper/build_application_tables.py --check
python -B research/paper/check_imt_integration.py
python -B research/paper/check_finite_integration.py
```

The finite integration check retains the seven original distance targets and
authenticates the additional K=2^16 two-round certificate. It checks six
optimized forward/transpose timing cells and the unchanged quarter-rate cell,
including the displayed timing ranges and memory sizes. Its 746-file count
describes the original evidence set; the 18 optimized timing files and separate
two-round certificate are authenticated additionally.

Omit `--check` on either table generator to regenerate its tables. These commands
do not benchmark. The historical build and serial benchmark scripts are linked
from FORWARD_SCHEDULE.md; they require the frozen research source snapshot and
its original directory layout. Raw measurements are intentionally not committed.
A clean code checkout therefore cannot authenticate these author-side receipts.

## Standalone-library discrepancy

A September 22 rerun of `spin/bench/paper.cpp` measured 9.661047 ms at K=2^20
with the correct 2048-row tile. Alternating the retained optimized executable
and the public library reproduced 9.210--9.261 ms versus 9.600--9.763 ms.
This is an unresolved roughly 5% consolidation regression, not evidence that
the historical 9.289 ms measurement was wrong. The paper reports the optimized
research implementation; it does not promise that timing from the current
standalone API. Fixing that regression is separate work.

The diagnostic raw-run digest is
`ad1cdff13eb4bd6927a0ae64f817171a750db902017f5a9741383cb1e6e1d2c9`;
the uploaded source archive digest is
`f77dc3a6f93bb56b6ed66e8e8e33dd6b05236c3962a41db19f394b4f5cb4a3c9`.
Those local files are under `measurements/paper-precomputed-20260922` in the
optimized workstream and are not inputs to the accepted paper tables.

## Unchanged measurements

Quarter-rate SPIN remains 15.254 ms. External encoder baselines retain their
comparison campaign. OT, PCS, and Flock retain their original implementation
campaigns; none of their application timings is rescaled from the new kernel
numbers. The abstract's approximate original-BAA comparison becomes 3.4x;
chosen Golay and RM BAA take 2.58x and 2.98x as long, respectively.
