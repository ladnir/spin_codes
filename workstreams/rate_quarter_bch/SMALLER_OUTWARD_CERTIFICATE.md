# Certified quarter-rate SPIN operating points

The fixed BCH-derived [128,32,32] outer and the existing RM2Sub (t,s)=(128,19)
inner now have outward certificates for both requested points at K=2^20.
The code length is N=2^22. For the setup distribution defined below,

| Target relative distance delta | Certified setup-failure probability | Margin diagnostic |
|---:|---:|---:|
| 0.165 | <2^-40 | 41.083485 bits |
| 0.19 | <2^-30 | 30.052513 bits |

Failure means that some nonzero message produces weight at most floor(delta N).
Thus the respective good-setup conclusions are minimum distance greater than
692,060 and 796,917. The displayed margin diagnostics are truncated decimals;
the inequalities are checked against exact powers of two, not rounded displays.

## Fixed construction and setup

Use the generator in `SMALLER_OUTER_AUDIT.json`, derived and checked in
[SMALLER_OUTER.md](SMALLER_OUTER.md). Encode each of the L=32768 message rows
with this fixed [128,32] outer. Independently permute the 128 coordinates in
each encoded row, transpose into 128 regions, and independently permute the
L coordinates in each region.

The fixed inner is the existing selected t128_s19 map, authenticated against
its retained selection and optimized implementation manifest. Its image
spectrum is enumerated exactly; its kernel spectrum is derived by integer
MacWilliams transform. Rank, BA=0, and map identity are checked again.

Serialize the permuted bits in region order and apply

    q_0=0,   Y_i=X_i+A q_i,   q_(i+1)=alpha_i q_i+A^T X_i.

Each alpha_i is an independent uniform nonzero element of GF(2^19), independent
of the permutations. State continues across region boundaries, without a
final flush. One setup is shared across all messages. Neither outer nor inner
map is sampled in this probability statement.

## Coverage and certified bounds

Let Z_Q count messages with exactly Q nonzero outer rows and bad output weight.
The certificate bounds E[sum_Q Z_Q]; Markov's inequality then bounds setup
failure. The computation covers every integer Q from 1 through L.

| Disjoint contribution | Margin at delta=0.165 | Margin at delta=0.19 |
|---|---:|---:|
| Q=1 | 41.083559 | 30.052988 |
| Q=2 | 55.323238 | 42.487246 |
| Q=3 through sparse endpoint | 66.742759 | 42.771257 |
| Remaining dense range | 307.927332 | 241.780698 |
| Combined | 41.083485 | 30.052513 |

The sparse endpoint is 64 at delta=0.165 and 128 at delta=0.19. Dense boxes
cover Q>=65 and Q>=129, respectively. Coverage of their integer row-count
simplexes is checked exactly, including the all-one outer word.

Q1 uses support averaging and the activation-aware three-state transfer.
Q2 uses 91 compositions of 13 outer-weight bands. This replaces the native
pair-support calculation with a looser bound that is simpler to certify.
For Q>=3 in the sparse range, the kernel-aware adaptive recurrence bounds
all band assignments. Consequently, the native binary64 pair kernel is not
part of the outward arithmetic implementation.

The analytic ingredients are the bridge's
[kernel-aware counting bound](../bch_rm2sub_bridge/GENERAL_OCCUPANCIES.md),
its intersected termination bounds, and the landscape's typed Cauchy bound
implemented in `typed_dense_boxes.py`. The latter holds each witness fixed
over a box and bounds its convex count dependence at every vertex. An exact
integer product of coordinate widths bounds the number of count vectors.

## Arithmetic and replay

`discover_outward_witnesses.py` selects tilts, probabilities, and boxes using
binary64. Its reported margins are not accepted as evidence by the certifier.
`certify_smaller_margins.py` recomputes every transfer, counting factor,
maximum, matrix power, and union with 256-bit Arb interval arithmetic.

Every recorded upper bound is an exact dyadic, rounded upward to a 160-bit
mantissa. Aggregate bounds are also formed outward and stored as dyadics.
The final comparisons with 2^-40 and 2^-30 therefore have no tolerance term.
Imported dense proposals are interpreted as exact binary64 rationals and
renormalized exactly; their density costs are recomputed rather than copied.

The 512-bit replay verifies that every recomputed contribution is at most
its recorded upper bound. Q1 and Q2 regional powers use linear epoch iteration
in replay, versus binary powering in production. Higher occupancies replay
the same recurrence at higher precision. Exact rational toy tests separately
check epoch envelopes, regional convolution, Q1 support averaging, and dyadic
rounding. This is not a claim of two independent derivations of the full proof.

## Artifact and reproduction

- `SMALLER_OUTWARD_WITNESSES.json`: selected witnesses and input hashes.
- `SMALLER_MARGIN_CERTIFICATE.json`: all contribution bounds and final unions.
- `SMALLER_MARGIN_CERTIFICATE_REPLAY.json`: 512-bit replay, bound to the producer receipt.
- `test_outward_margins.py`: exact toy checks and certificate regressions.

Use Python with NumPy, SciPy, and python-flint; production used python-flint
0.9.0. The selected-map implementation manifest must be available; regenerate
it with `python workstreams/bare_bch_rm2sub/generate.py` if needed. From the
repository root:

```sh
python workstreams/rate_quarter_bch/discover_outward_witnesses.py
python workstreams/rate_quarter_bch/certify_smaller_margins.py
python workstreams/rate_quarter_bch/certify_smaller_margins.py --verify
python -m unittest discover -s workstreams/rate_quarter_bch -p 'test_*.py'
```

With retained witnesses and certificates, run only the last two commands
to replay and test. Source hashes bind the producer, inputs, map data, and
numerical witness receipt. Earlier binary64 exploration receipts remain
unchanged; they are not substitutes for this certificate.

Next: implement and benchmark the complete quarter-rate transposed encoder.
No quarter-rate timing is claimed here, and the parked t256 investigation
does not affect these certificates.
