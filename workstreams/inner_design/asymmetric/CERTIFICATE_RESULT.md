# First IMT certificate

The quarter-rate BCH [128,32,32] SPIN instance now passes both requested
operating points with the [IMT inner](../IMT.md). Its parameters
are K=2^20 message bits, N=2^22 output bits, t=128, s=19, balanced expansion A,
feedback map `greedy3_2`, and one sampled transvection per epoch.

| Relative distance | Required margin | Margin from saved upper bound |
|---|---:|---:|
| 33/200 = 16.5% | 40 bits | 41.0481676058 bits |
| 19/100 = 19% | 30 bits | 30.0334910634 bits |

For each row, over the shared random SPIN routing and independent per-epoch
mixer setup described in `TRANSFER_ARGUMENT.md`, the probability that some
nonzero message produces an output of weight at most floor(delta N) is less
than 2^(-required margin). Thus a setup outside that bad event has relative
minimum distance greater than delta. The decimal margins are diagnostic
renderings; the verifier checks the claimed 40- and 30-bit inequalities with
exact dyadic arithmetic.

Every outer occupancy Q=1..32768 is covered. The 16.5% calculation handles
Q=1..64 separately and uses 139 dense cover leaves for the remainder. The
19% calculation handles Q=1..128 separately and uses 318 dense leaves.
The 19% operating point has only about 0.0335 bits of spare certified margin;
it must not be rounded into a stronger claim.

## Checks and artifacts

- `ASYMMETRIC_MARGIN_CERTIFICATE.json`: all union terms evaluated with
  256-bit Arb outward arithmetic and saved as upward dyadic bounds.
- `ASYMMETRIC_MARGIN_CERTIFICATE_REPLAY.json`: every term recomputed at
  512 bits and checked against the saved bound. This uses the same producer.
- `CERTIFICATE_VERIFICATION.json`: source hashes, exact union and coverage,
  reconstructed independent maps, generated map header, outer-generator row
  space, and existing C++ correctness receipts checked by a separate verifier.
- `test_asymmetric_outward.py`: the actual outward fixed-weight and Bernoulli
  transfers checked against exact rational laws in a small state space with
  BA nonzero. Together with `test_asymmetric.py`, all six tests pass.

The producer reconstructs the A and B^T spectra separately. It does not use
the symmetric-map Q=1 or Fourier-intersection identities. Its dense witnesses
choose a complete valid transfer family rather than mixing entries from
different representations. The bound and its assumptions are derived in
`TRANSFER_ARGUMENT.md`.

## Reproduce

From the repository root, with `workstreams/inner_design/requirements.txt`
installed, run the following in order. Do not use Python's `-O` option:
these research verifiers use assertions for failed checks.

```text
python -m unittest discover -s workstreams/inner_design/asymmetric -p "test_*.py" -v
python workstreams/inner_design/asymmetric/certify.py
python workstreams/inner_design/asymmetric/certify.py --verify
python workstreams/inner_design/asymmetric/verify_certificate.py
```

To check the retained artifact without recomputing its numerical terms, run
only the final command. That is an integrity and exact-union check of the
saved replay receipt, not a substitute for rerunning the 512-bit calculation.
Witness discovery and performance benchmarks are not needed for this replay.

## What this enables

This closes the first finite-instance gate toward replacing RM2Sub. It does
not certify other K values, other outers, the mixed feedback candidate, or
the candidate's asymptotic distance. The supported default and paper theorem
remain unchanged. `MIGRATION_PLAN.md` records the theorem dependencies and
the remaining finite-certificate, asymptotic, and forward-performance gates.

Next prioritize the paper's BCH-256, t=128,s=19 cells at K=2^16,2^18,2^20.
Their existing rigorous outer-spectrum envelopes can be reused without an
exact BCH-256 spectrum. A complete optimized forward encoder is a separate
engineering gate; transpose timings alone do not establish its performance.
