# Inner-design investigation

This is a code-and-documentation checkpoint. New data files, including map
JSON, numerical certificates, replay receipts, and benchmark records, remain
local and are not included in this commit. Reported results describe those
local runs. Receipt-dependent verification needs the data or a regenerated
artifact; this checkout alone is not the complete verification bundle.
The latest paused investigation is in [the BCH-256 overlap note](asymmetric/bch256/overlap/README.md).

The certified quarter-rate implementation is unchanged. The two-part goal is
recorded in [GOAL.md](GOAL.md).

1. **Equivalent implementations:** practical search completed, no candidate
   adopted. See [BASIS_PERFORMANCE.md](BASIS_PERFORMANCE.md).
2. **Different inner mixer:** the balanced-image one-round rank-one update
   has a full 256-bit outward certificate at K=2^20: 41.048 bits at 16.5%
   distance and 30.033 bits at 19%. Its exact implemented map is about 4%
   faster in serial comparisons. The 512-bit replay passes every saved bound.
   See [BALANCED_RESULT.md](BALANCED_RESULT.md) for the consolidated result.

A same-size map without an all-one expansion state is now exactly audited.
Its image weights lie in [48,80], and its kernel minimum weight is five.
[SPECTRUM_REDESIGN.md](SPECTRUM_REDESIGN.md) explains the construction and tradeoff.
Both isolated implementations pass the dense reference. The original-map
one-round mixer remains a separate, Q=1-only certified experiment; its timings
are retained in [MIXER_PERFORMANCE.md](MIXER_PERFORMANCE.md).

## Layout

| Files | Purpose |
|---|---|
| `search_basis.py`, `BASIS_SEARCH*.json` | Exact state-basis search and compact candidates |
| `generate_candidates.py`, `generated/`, `CMakeLists.txt` | Isolated equivalent C++ encoders; original dense oracle retained |
| `measurements/`, `summarize_measurements.py`, `BASIS_PERFORMANCE.json` | Serial full-encoder measurements and checked summary |
| `transvection.py`, `TRANSVECTION_Q1.json` | Exact mixer law; deliberately coarse three-state screen |
| `weight_memory.py`, `WEIGHT_MEMORY_Q1.json` | Stronger pointwise weight-class Q1 screen |
| `WEIGHT_NEIGHBORS.json` | Exhaustive 67,108,864 state/input-pair diagnostic; not required by the successful proof |
| `certify_q1.py`, `TRANSVECTION_Q1_CERTIFICATE.json` | Outward arithmetic and independent-precision replay of Q1 |
| `generate_mixer.py`, `MixerInner.h`, `MIXER_PERFORMANCE.*` | Three equivalent implementations of the new mixer and serial timings |
| `general_occupancies.py`, `FIBER_BOUNDS.md` | Exact syndrome-fiber caps and sparse/dense diagnostics |
| `refine_dense.py`, `probe_periodic.py` | Bounded alternative dense-proof probes; no certificate claims |
| `remove_constant.py`, `NO_CONSTANT_MAP.json` | Exactly audited balanced-image map at the same t,s |
| `dense_cover.py`, `NO_CONSTANT_DENSE_COVER_*.json` | Adaptive binary64 witness discovery and exact dense partitions |
| `certify_no_constant.py`, `NO_CONSTANT_MARGIN_CERTIFICATE.json` | Full outward all-occupancy certificate and 512-bit replay command |
| `generate_balanced.py`, `BALANCED_IMPLEMENTATION.json` | Final proof-compatible map, symbolic circuit checks, two isolated kernels |
| `summarize_balanced.py`, `BALANCED_PERFORMANCE.json`, `measurements/balanced/` | Final-map serial timing and correctness receipts |
| `verify_balanced_artifact.py` | Fast binding check after full replay: certificate, exact map, outer span, and performance |
| `test_*.py` | Exact matrix identities, exhaustive toy mixer laws, rational toy transfer checks |

The separate neighbor-envelope experiment is retained because it identifies
a trap: maximizing each destination density independently can lose thirteen
bits. Remembering weights on **empty** epochs and forgetting them on nonempty
lazy epochs is substantially tighter here.

## Reproduce the mathematics

From the repository root, with Python, NumPy and python-flint available:

```text
python -m unittest discover -s workstreams/inner_design -p "test_*.py" -v
python workstreams/inner_design/transvection.py
python workstreams/inner_design/weight_memory.py
python workstreams/inner_design/certify_q1.py --verify
python workstreams/inner_design/certify_no_constant.py --verify
python workstreams/inner_design/verify_balanced_artifact.py
```

The two certificate commands replay at 512 bits against saved 256-bit upward
dyadics. The last command checks artifact bindings; it does not substitute for
the full replay. No replay trusts binary64 bounds or the neighbor table.
To produce a new receipt deliberately, omit `--verify`. Source-byte checks
reject a replay if a bound input or producer has changed.

## Next gate

The two-part investigation is complete. Integrate the selected variant and
its finite-length statement deliberately as the next task.
The production encoder remains unchanged throughout this investigation.
