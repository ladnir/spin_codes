# Bolt-max PCS estimate at 512 MiB

The headline PCS table includes **2,614 ms estimated prover time** and
**11.1 MiB estimated communication** for Bolt-max in the amortized limit.
Both are marked with an asterisk. Neither is a measured complete Bolt proof,
a runtime/communication bound, or a new composed-security claim. Verification
is unmeasured and appears as a dash. The commitment measurement now includes
our x86 expander optimization; the opening and communication models are unchanged.

This row answers the large-proof comparison question: our longer-row SPIN
configuration has a measured 12.4 MiB proof, reasonably close to this model.
The comparison matches 512 MiB input volume, not polynomial dimension or
coefficient field: Bolt has 2^27 GF(2^32) entries; SPIN has 2^32 Boolean entries.
It does not isolate the encoder's contribution from field/layout/hash choices.

## Timing

The table generator reads the existing pinned summaries, without modifying
their measurements:

- Commitment: 1,341.6295525 ms, `bolt-max-xor-d8`, SHA-256,
  from `paper/data/bolt_x86_commitment.json`.
- One amortized-limit opening: 1,272.767206 ms,
  from `paper/data/bolt_opening_projection.json`.
- Sum: 2,614.3967585 ms. This sums independently aggregated phase medians;
  it is not a median of end-to-end runs.

The opening retains the existing row evaluations, random fold, two sumchecks,
and complete Ligerito proofs as proxies for the inner RS work. It does not
time all basis preparation, outer authentication/serialization, or transcript
integration. The non-amortized scenario remains in the appendix (1,728.663130 ms
opening, 3,070.2926825 ms including commitment). We do not give that variant a
proof-size row, since the matrix-proof communication has not been reconstructed.

### x86 commitment refresh (2026-09-17)

We specialize the XOR expander's 128-symbol columns and replace its per-column
heap accumulator with a fixed-size array. LLVM keeps the accumulator in SIMD
registers; no handwritten assembly was added. The graph, code parameters, RS
kernel, and hashing are unchanged. All baseline/optimized commitment roots match.
The appendix describes this implementation change.

At 512 MiB, a fresh comparison reduces expander time from 732.9903245 to
330.9322805 ms and full commitment time from 1757.1455855 to 1341.6295525 ms.
The same optimized binary commits 32 and 128 MiB in 63.0015895 and 307.7880395 ms.
These replace the commitment components of the Flock projections as well:
the current surrounding SPIN work plus unchanged shared-opening estimates gives
228 and 983 ms. The historical input records are retained unchanged.

Each timing is a median of six measured trials across two processes, with one
excluded warmup per process. The order is baseline/optimized/optimized/baseline.
Runs are serial on Peach CPU15, with one worker, requested 4.5 GHz and boost off.
Both binaries use Rust 1.94.1, native CPU features, release LTO, and locked dependencies.
`paper/data/bolt_x86_commitment.json` pins the source and executable hashes and
the three summary receipts. The Bolt worktree's `tools/standalone/X86_EXPANDER.md`
and `run-xor-comparison.py` document and reproduce the optimization campaign.
The starting revision is `9b45089`; the optimization branch is
`codex/bolt-x86-expander`, not yet a published upstream release.

## Communication model

The [Bolt paper](https://eprint.iacr.org/2026/310), Sections 4, 5, 7.3, and 8,
provides the protocol components and query convention. The inspected May 20,
2026 PDF has SHA-256
`e9de9a0a047f1967b84c80171279d30755e8c70d7e60a3facd3455fcdec49f0d`.

Use interleaving width 128, four-byte input symbols, 16-byte challenge-field
elements, 32-byte digests, k=2^20 message columns, and k/4 syndrome columns.
With nominal target 100 and gamma=0.013:

    q_X = ceil(100 / -log2(1 - gamma/3)) = 15962
    q_Y = ceil(100 / -log2(3/4)) = 241

The earlier upstream helper's gamma=0.095 is not Bolt-max's parameter.
Its leaf-plus-path subtotal must not be called a full proof.

The model uses uniform distinct queries in each independent pool, following
the earlier size-accounting convention, not a claim about a finished Fiat--Shamir
implementation. Sampling with replacement and deduplicating would change the
number of transmitted columns slightly.

For a perfect n-leaf tree with q distinct queries, an s-leaf subtree is empty
with probability `binom(n-s,q)/binom(n,q)`. Let O_h be the expected number of
occupied nodes at level h, leaves at zero. The expected sibling count is
`sum(O_h, h=1..log2(n)-1) + 2 - q`. The X and Y trees join with two dummy
sibling digests in the pinned implementation. Exhaustive tests check both
the expectation and the odd-tree join on small trees.

| Component | Bytes |
| --- | ---: |
| Opened outer columns | 8,295,936 |
| Expected outer Merkle payload | 2,726,709.689 |
| Inner message proof proxy | 451,184 |
| Inner syndrome proof proxy | 177,200 |
| Explicit vector of 128 row evaluations | 2,048 |
| Two 20-round sumchecks, three field elements per round | 1,920 |
| Allowance for eight scalar values | 128 |
| **Total scenario** | **11,655,125.689 (11.115 MiB)** |

The scalar allowance accounts for relation/evaluation endpoints rather than
asserting an exact serialized transcript. The proxy byte counts include roots,
opened data, authentication, and their own sumcheck transcripts. They are the
two proxies already used in the runtime model, on extension-field messages of
lengths k and k/8 with 241 queries.

These are not the constrained/virtual/batched inner proofs of Bolt. In particular,
the syndrome proxy materializes a separate commitment; the model also charges
outer Y-column openings. This deliberately retains the standalone proxy instead
of claiming to reconstruct virtual-oracle communication exactly. Constraints,
batching, and different recursion choices can change the actual proof. No
arbitrary multiplicative accuracy guarantee is claimed.

We assign negligible per-opening time and bytes to the matrix proof in the
amortized limit. A finite batch must pay its share and batching work. The model
excludes the original commitment root, public query indices derivable from the
transcript, and serialization/container framing.

## Provenance and reproduction

The communication calculator is `paper/bolt_communication.py`; it writes no
files and runs no benchmarks. The inner sizes are decoded from the `checksum`
field emitted by `proof.size_bytes()` in the existing verified calibration:

- Bolt calibration worktree revision: `9b45089`.
- Upstream: `bcc-research/bolt-rs` revision
  `3832e47b24e7b3e10525c9c5bcfc1cfe66d525f2`.
- Receipt: `tools/standalone/opening-results/27.csv`.
- SHA-256: `e421aa8167a9570a533cea846a5f4d5202985b39fe3ce68bc91bb8f8e3164723`.
- Six identical size receipts per proxy (warmup plus five measured trials).
- Source: `tools/standalone/src/bin/opening-kernels.rs`, `inner_proxy`;
  `src/ligerito_recursive.rs`, `LigeritoProof::size_bytes`.
- The calibration verifies each proxy proof outside its prover timer.

From the manuscript repository:

```sh
python paper/bolt_communication.py
python -m unittest discover -s paper -p 'test_bolt_communication.py'
python paper/build_application_tables.py --check
```

To authenticate the retained external receipt as well:

```sh
python paper/bolt_communication.py --source PATH_TO_BOLT_STANDALONE
```

The outer calculation reuses the method in Hypercat's
`results/bolt-proof-size/account.py`, explained in
`docs/bolt-proof-size-accounting.md`. The earlier paper-anchored estimate
subtracts the calculated change in outer authentication from the reported
12.92 MB at 4 GiB input, retaining the unknown remainder. It yields
10.8587 or 11.4572 MiB here for decimal-MB or binary-MB interpretations.
These are two modeling scenarios, not confidence limits. Agreement with the
11.115 MiB component model is useful context, not independent validation of
a full Bolt proof; the published total has no mode-specific byte receipt.

Only the amortized-limit scenario appears in the headline table. Replacing
proxies with actual constrained-proof byte receipts remains future work.
The older Brakedown-code baseline has not been changed or benchmarked here.
