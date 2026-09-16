# Quarter-rate SPIN with the existing inner

Status, 2026-09-11: all six tested configurations have full numerical coverage
of the nonzero messages. The existing optimized RM2Sub map, t=128 and s=19,
was held fixed. No inner search or implementation change was performed.

These are nearest-binary64 evaluations of analytic first-moment bounds.
They are not outward-rounded distance certificates or measured failure rates.
The results support proceeding with this inner before revisiting its parameters.

## Results

Let K be the number of message bits and N=4K the code length. For each target
delta, the bad event is an output weight at most floor(delta N) for some
nonzero message. The table gives minus log2 of the computed union bound.

| K | Distance target 0.10 | Distance target 0.19 |
|---|---:|---:|
| 2^16 | 149.91 bits | 66.78 bits |
| 2^18 | 154.66 bits | 79.35 bits |
| 2^20 | 158.92 bits | 80.50 bits |

Displayed values are truncated to two decimal places. Exact floating-point
outputs and input hashes appear in `FIXED_INNER_SUMMARY.json`.

Every cell exceeds the desired 40-bit target in this numerical evaluation.
In particular, the stronger distance target 0.19 clears at K=2^20 without
changing the current inner. Converting the retained bounds to an outward
certificate remains the next proof step.

These values are not a fitted growth curve or optimized margins. The dense
search stops once its bound clears 60 bits, and different cells retain
different amounts of slack. The increasing displayed margin must not be
interpreted as increasing true security with K.

## Construction and setup randomness

The outer is the fixed BCH-derived [256,64,62] code in `CONSTRUCTION_AUDIT.json`.
Its exact spectrum is `BCH256_64.wd`; `README.md` derives it from published
parent spectra. Each message has L=K/64 outer rows.

Independently permute the 256 coordinates of each encoded row. Transpose
into 256 regions of L bits, then independently permute each region. Process
the resulting N bits with the existing RM2Sub recurrence:

    q_0 = 0
    Y_i = X_i + A q_i
    q_(i+1) = alpha_i q_i + A^T X_i

Here X_i and Y_i contain 128 bits, and q_i contains 19 bits. Each alpha_i
is independently uniform in the nonzero elements of GF(2^19). State continues
across region boundaries; there is no final flush. One setup fixes all
permutations and multipliers for every message. Failure probability concerns
this setup, not the choice of the fixed outer or inner maps.

The inner is exactly the selected map in
`../bch_rm2sub_bridge/generated/larger_state_inputs_v1/t128_s19_selection.json`.
`load_inner()` compares its ordered columns and selection hash against the
optimized encoder's generated manifest. It independently enumerates the
image spectrum and computes the kernel spectrum by integer MacWilliams
transform. It also checks rank, transpose identity, BA=0, and distinct
nonzero columns. The image has minimum weight 48; the kernel has minimum
weight 6. These checks authenticate the existing map rather than select a
new one with matching dimensions.

## How the coverage is assembled

An occupation Q is the number of nonzero outer rows. The calculation splits
the range into disjoint parts:

- Q=1 uses the exact outer spectrum and the activation-aware three-class
  transfer. The numerical search uses log-surprisal values from -18 to 0.
- Q=2 through 64 use the kernel-aware regional transfer and deterministic
  outer spectrum bands. Separate composition bounds refine Q=2,3,4.
- Q=65 through L use the existing typed-box inequality. The row types are
  zero, ordinary nonzero, and all-one. The all-one word remains an exact
  atom; it is not absorbed into the ordinary Bernoulli counting measure.

The sparse contribution sums all covered integer occupations. The dense
contribution sums bounds over a disjoint cover of the row-type counts.
`verify_fixed_inner_results.py` checks that cover using integer intervals
for every Q. It then replays every selected dense witness and combines the
sparse and dense contributions.

At K=2^20 and delta=0.19, the sparse and dense contributions separately give
81.887 and 81.199 bits. Their sum gives 80.502 bits. The first coarse dense
attempt did not close; finer numerical witnesses and extending the sparse
range to Q=64 resolved that failure without changing the encoder.

The verifier checks combinatorial coverage exactly, but evaluates moments
in binary64. That distinction is why this work does not yet claim a formal
setup-failure certificate.

## Reproduction

Requirements: Python, NumPy, and SciPy. No benchmark runs are involved.
The existing encoder's generated manifest must be present. If it is absent
in a fresh checkout, reproduce the existing generated code first:

```sh
python workstreams/bare_bch_rm2sub/generate.py
```

Generate the sparse results for all three message sizes and both targets:

```sh
python workstreams/rate_quarter_bch/evaluate_fixed_inner.py
```

Run these commands for exponents 16,18,20, replacing 20 in both the argument
and output filename:

```sh
python workstreams/rate_quarter_bch/evaluate_dense_fixed_inner.py \
  --exponent 20 --distance 1/10 --nodes 127 \
  --output workstreams/rate_quarter_bch/DENSE_E20_D10.json
python workstreams/rate_quarter_bch/evaluate_dense_fixed_inner.py \
  --exponent 20 --distance 19/100 --nodes 511 \
  --output workstreams/rate_quarter_bch/DENSE_E20_D19.json
```

To replay all retained dense witnesses and regenerate the combined summary:

```sh
python workstreams/rate_quarter_bch/verify_fixed_inner_results.py \
  workstreams/rate_quarter_bch/DENSE_E16_D10.json \
  workstreams/rate_quarter_bch/DENSE_E16_D19.json \
  workstreams/rate_quarter_bch/DENSE_E18_D10.json \
  workstreams/rate_quarter_bch/DENSE_E18_D19.json \
  workstreams/rate_quarter_bch/DENSE_E20_D10.json \
  workstreams/rate_quarter_bch/DENSE_E20_D19.json \
  --output workstreams/rate_quarter_bch/FIXED_INNER_SUMMARY.json
python -m unittest discover -s workstreams/rate_quarter_bch -p 'test_*.py'
```

The new workstream has 16 passing tests. The reused Q1, arbitrary-occupation,
and typed-box modules also pass their 12 existing tests, including exact
small-instance comparisons. The full sparse numerical search is reproduced
by its producer; the standalone verifier does not independently recompute it.

## Next steps

First, outward-round the retained witnesses for a 40-bit claim with the
current inner. No exhaustive search over alternative inner maps is needed
for that step. Then implement and measure the quarter-rate transposed
encoder, reusing the existing inner kernel and routing design where applicable.

Revisit inner parameters after establishing this baseline, as requested.
No quarter-rate timing is reported here. The rate-half timing cannot be
reused unchanged: at fixed K, quarter rate processes twice as many output bits.
