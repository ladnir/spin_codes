# Regular-noise Silent OT paper update

The abstract, introduction, and PCG section lead with **89.0 million hashed
OTs/s** at K=2^18: 2.943876 ms per sender batch, rounded to 2.944 ms.
This uses the full precomputed BCH--IMT (128,19) code with weight-five feedback,
VAES, and opt-in streaming output stores. No heuristic permutation refresh.
The appendix also records K=2^20: 17.185599 ms, or 61.0 million OTs/s.

These are local sender timings with fresh base correlations supplied beforehand,
not network throughput. They include full tree expansion, compression, and
hashing both sender messages. Setup, base-correlation generation, workspace and
output allocation, transport, and output consumption are excluded. The noise
configurator's effective-distance value 0.25 is a tuning input, not a certificate.

The matched stationary run is 3.103495 ms (84.47 million OTs/s); the preceding
stationary table is retained as a separately identified campaign. Different base
correlation requirements prevent inferring complete-protocol cost from this comparison.

Measurements: Ryzen 9 7950X, CPU 15, requested 4.5 GHz, boost disabled;
libOTe/cryptoTools GCC 13.4 with VAES, SPIN GCC 15.2 Zen 4 archive. Each figure
is the median of three process medians, each from 101 batches after five warmups.
All runs were serial. Separate paired tests check the OT outputs across batches.

The source report and reproduction commands are in the adjacent libOTe checkout:
`libOTe/Tools/Spin/REGULAR_PERFORMANCE.md`. Its summarizer is
`libOTe_Tests/summarize_spin_noise_bench.py`. Retained raw results are on the
measurement host under `/tmp/libote-spin-ssd-pq75QW/results/regular-20260923`;
no raw samples are committed.

| Measured snapshot | SHA-256 |
|---|---|
| Benchmark source | 47823a52f1c491e542ebbce525c611d3c6dee71afd9bdfd7d2d99e76c6e4ee61 |
| Summarizer | 7a883140ecb7b77db4c2818224d4cabbce0508b90dfeaaa5261cc7e4eb696199 |
| Benchmark executable | 589c2f646aaa3c64847c907563ebcc2530ab30340a97733695ae0f96aa693c09 |
| SPIN archive | ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7 |

The precomputed encoder-only results and historical OT baseline are unchanged.
In particular, the integrated K=2^20 compression time is not replaced by the
9.29 ms optimized-research encoder timing.
