# Measured SPIN Silent OT performance

Paper labels: `sec:spin-pcg`, `app:pcg-details`, `tab:spin-ot`.
These measurements replace the earlier encoder-plus-estimated-tree projection.
They measure one party's computation, not a two-party network session.

## Results and aggregation

At K=2^20, each entry is the median of nine calls after one excluded warmup.

| Timed work | Sender (ms) | Receiver (ms) |
| --- | ---: | ---: |
| Expansion + SPIN | 19.530826 | 20.625179 |
| Expansion + SPIN + hashing | 22.098610 | 21.841110 |

For each row of `rot-final.log`, add that party's `expand_map_ms` and
`hash_ms` before taking the median. Do not add independent phase medians.
Discard trial -1. All nine retained trials in each mode have `valid=1` and
check all 1,048,576 outputs; the warmups also pass.

The slower hash-enabled endpoint corresponds to 1,048,576 / 0.022098610,
approximately 47.45 million outputs per second of single-core computation.
Do not sum sender and receiver latencies and call the result one-party cost,
or describe this number as measured network throughput.

## Timing and parameter contract

- Ryzen 9 7950X, CPU 15, requested 4.5 GHz, boost disabled; GCC 15.2.0,
  `-O3 -march=znver4`.
- Base correlations are already installed. The harness creates matching
  synthetic bases before timing; it does not measure base OT generation.
- Both endpoints use `SilentSecType::SemiHonest`; no malicious-security
  consistency check is included in these timings.
- Public SPIN setup, independent endpoint workspaces, buffer allocation and
  touching occur before timing.
- Sender `silentSendInplace` runs first. Its actual buffered transcript feeds
  receiver `silentReceiveInplace`. Endpoints run serially with separate timers.
- Receiver choice packing is timed. Network latency, transcript handoff,
  output validation and object destruction are excluded.
- No-hash output packs choices into low bits and checks the correlated-OT
  relation in the other 127 bits. The hash-enabled mode checks every selected
  receiver message against the corresponding sender message.
- BCH [256,128], IMT t=128, s=19, weight-five feedback; regular binary noise
  with 400 partitions of 5,242 coordinates. The final 352 of N=2^21 coordinates
  are zero and charged against the distance input used in noise configuration.
- Online communication: 179,624 bytes, excluding bases.

This is a performance integration at K=2^20, not a production-interface release
or a new structured-LPN hardness theorem. Encoder dense-oracle tests also pass
at K=2^16, 2^18 and 2^20. No integrated baseline speedup is claimed.

## Source and evidence pins

libOTe branch `codex/spin-ot-perf`, commit
`910f88337fda607616bf2cd99ce2a0a974226472`, based on `2b686ad`.
The source repository is https://github.com/osu-crypto/libOTe; the benchmark
commit is currently retained in the local branch, not asserted to be upstream
or publicly released. Source checkout on the author machine:
`C:/Users/peter/.codex/worktrees/spin-ot-perf/libOTe`.

Dependencies: cryptoTools `e4ba77e1009f94af6099aa3000e465b1d01dd132`,
coproto `ded64cbc51c041ee534e5c34a5a81135af455ce7`,
macoro `91fc6b42ff719c681713c6d9f7476b94fd822983`.
The branch's `perf/README.md` documents the encoder import and build contract.

Measured executable SHA-256:
`63d96c98350ef8a936da321ee5ecd60ad43b4fe87605479a078d94607d3cddd6`.
Correctness executable SHA-256:
`3533059bc7fd6285e01555d16111c19b2610e82bbe4bc932970e59c8e3cef26b`.

Logs are retained in that checkout's ignored `perf/results/` directory:

| File | SHA-256 |
| --- | --- |
| `cot-final.log` | `9ff9740cb322adee7cfce767361d12c4478ddbc46cc30482620940b90955eaca` |
| `rot-final.log` | `3c5cff0358fa56f8919c2a16048192b138a09bced9faaac781ef36092332e383` |
| `cot-phases.log` | `ae221cced84f055f65d6898def899b008871a7c8cf7378257c0b4636147450a9` |
| `rot-phases.log` | `dff0d669391881c7f0f292484dda1bf1fdea56844ec02d193efaa6d46e6a854c` |
| `correctness.log` | `f71f06a1259e089e3cfddb641a0f5b4df3149859a70a5ab3404814d9a54baa7d` |

On 2026-09-17 the paper integration recomputed all four medians from these
logs and verified the trial validity flags. No benchmark was rerun and no
raw data was copied into the paper source tree. External/anonymous release
must include access to the pinned source and logs; this document alone is not
a complete reproduction archive. Strip author paths from any anonymous package.

## Reproduce on the benchmark host

From the pinned libOTe checkout, with no other benchmark running:

```sh
git submodule update --init --recursive
cmake -S . -B out/build/linux -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-march=znver4 \
  -DENABLE_SILENTOT=ON -DENABLE_MOCK_OT=ON -DENABLE_LOGVOLE=OFF \
  -DENABLE_BOOST=OFF -DFETCH_AUTO=ON
cmake --build out/build/linux --target spin_ot_perf spin_perf_correctness -j 8
bash perf/run-spin.sh
```

The runner holds the shared benchmark lock and runs modes sequentially.
`ENABLE_MOCK_OT` is for this measurement harness; timed calls assert that
bases are already installed.
