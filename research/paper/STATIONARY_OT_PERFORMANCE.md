# Stationary Silent OT timing update

The paper reports **3.100658 ms** for a sender batch of 262144 hashed OTs:
**84.5446 million OTs/s**, rounded to 3.101 ms and 84.5 million in prose.
This is measured production-path computation, not a projection from encoder speed.

## Scope

- Semi-honest stationary Silent OT, MultType::Spin, packed choices, K=2^18,
  N=2^19, T128S19, weight-five feedback, BankedHeuristic refresh.
- 312 trees of 1680 leaves and 128 zero-padded coordinates. Each party caches
  8,386,560 bytes of leaf seeds. Each batch consumes 312 fresh base-VOLE correlations.
- Timed: fresh leaves and sums, padding, code refresh, transposed encoding,
  AES hashing of both sender messages, and routine allocation/message enqueueing
  within those production calls.
- Untimed: initial tree and code preparation, workspace/output allocation,
  generation of supplied base correlations, network transport and transcript
  draining. Output consumption is timed only in the explicitly labeled scan control.
- The permutation bank is a heuristic code family, not the uniform setup
  covered by the paper's distance certificates. Implementation agreement is
  checked; no new distance or LPN-hardness claim is made for this family.
- The isolated sender is timed after an untimed two-party correctness check and
  receiver cleanup. Separate paired runs check all OTs across refreshes.

## Results and aggregation

Ryzen 9 7950X, CPU 15, requested 4.5 GHz, boost disabled. libOTe/cryptoTools
use GCC 13.4 Release with VAES; SPIN uses the existing GCC 15.2 Zen 4 archive
with runtime-dispatched AVX-512 BCH. Every process is run serially under the
shared benchmark locks. Each cell is the median of three process medians,
each with five warmups and 101 timed batches. Stage medians need not add to totals.

| Sender stage | Cached output ms | Streaming output ms |
|---|---:|---:|
| Fresh leaves | 0.941248 | 0.917282 |
| Refresh and compression | 2.213212 | 1.752481 |
| Hash both messages | 0.423460 | 0.428821 |
| Total | 3.576958 | 3.100658 |
| Million OTs/s | 73.2869 | 84.5446 |

Total process medians: cached 3.579172, 3.576958, 3.557601 ms;
streaming 3.132598, 3.100658, 3.075361 ms. Reversing order in a second campaign
gave 3.571167 and 3.068929 ms; the paper uses the primary campaign above.

An immediate sequential scan of all output bytes gives totals 3.884151 ms
cached and 3.424693 ms streaming (67.49 and 76.55 million OTs/s). This includes
both scan time and next-batch cache effects, but is not a real protocol consumer.
Streaming stores remain opt-in. They barely help at K=2^16 in this experiment.

The 1.418608 ms / 184.8 million blocks/s fixed-code encoder measurement and
9.289382 ms K=2^20 table entry remain unchanged; see
[PRECOMPUTED_PERFORMANCE.md](PRECOMPUTED_PERFORMANCE.md).
The historical regular-noise 47-million-OT/s measurement uses a different K,
noise configuration, and implementation. It remains in the appendix as historical
context, not a matched baseline or claimed speedup.

## Retained receipt

The benchmark is `libOTe_Tests/SpinStationaryBench.cpp` in the libOTe checkout.
`run_spin_store_bench.sh` runs store/scan modes serially. Raw samples are retained
outside source control at `/tmp/libote-spin-ssd-pq75QW/results/nt/` on the measurement
host. The paper and core-code artifact do not require these raw files to build.

SHA-256 identities of the measured snapshot and primary campaign:

| Item | SHA-256 |
|---|---|
| Benchmark source | 511d31b3a8b8963af537d14229b45a76d281618216f9f03bab2d06271feb02d3 |
| Sender source | 5043efef985c657fb59fe337ed14e1b51148da1bb938fd14b5319494b31f26c3 |
| AESVaes.h | 4a697635583b2449ccc750b972925dbe077b61f917966d7ab936f38c9d8d2f2f |
| Benchmark executable | f62027eb97a673e1f9fa1e3d1c085cf72313b232e249fe25498fc5b27d3c7afe |
| SPIN archive | ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7 |
| cached-1.jsonl | 6abe4c858181d33373420bf829bf166b1ebe1f612fcd77bcc010c6cf0d8eec38 |
| cached-2.jsonl | 29c56ae4bd37c2db64f8ee7f01200ebfe05041120057a4962b382fcf5e5466f0 |
| cached-3.jsonl | 8401b177ccf4d3389aed58e3ccc2f16b316c70e079c4b09b3c8776a4727fc9ec |
| nt-1.jsonl | 7833ade976062c407ac1d9c272d0d7b06c65808dbead0c4ece1d4222ce5229e5 |
| nt-2.jsonl | 2648f58c245e5df6b25be155ae2419e7336375d696ec5d04fa484b9f26ef3d81 |
| nt-3.jsonl | 0160e4adbfa27ff41f64760e69df8e00feb7d13cfb5446b6e0c9d460c8ae9437 |

Code correctness checks passed with GCC VAES on and off and MSVC Release.
The streaming hash tests cover lengths 0 through 33, tails, output guards,
16-byte offset alignment, and equality to scalar AES. No raw measurement data
is added to the supplemental core-code artifact by this paper update.
