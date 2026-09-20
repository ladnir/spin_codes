# Outward arithmetic certificate

## Claim

Fix the explicit shortened-XBCH constituent recorded in
`scripts/xbch64_32_philips_manifest.json`. Use 16 copies per 1024-bit outer
block, followed by two independently interleaved accumulators. Use the
recorded SplitState inner transfer with `t=128` and `s=16`.

For output length `2^21` and target distance

```text
floor(0.09 * 2^21) = 188743,
```

the expected number of nonzero messages whose output has smaller weight is at
most

```text
2^-40.448462678054.
```

The expectation is over the independent interleavers and randomizers specified
by the construction. This statement imports the activation and live-spectrum
bounds named in the receipt.

## Message partition

Set the endpoint-tail width to 97. Each nonzero outer block is a tail block or
a body block according to its output weight. The proof sums five disjoint
classes:

1. exactly one tail block and no body block;
2. exactly two tail blocks and no body block;
3. at least three tail blocks and no body block;
4. at least one body block and no tail block;
5. at least one body block and at least one tail block.

These classes exhaust all nonzero messages.

## Outer spectrum

The constituent spectrum is an exact integer table of mass `2^32`. The direct
sum of 16 constituents is computed as an exact `fmpz` polynomial power.

For input weight `w` and accumulator output weight `h`, the uniform-interleaver
transition is an exact ratio of binomial coefficients. Arb encloses each ratio
at 256-bit precision. The two accumulator stages use the upper binary64
endpoints of these intervals. Every product and sum is rounded toward positive
infinity.

## Inner transfer

For each selected Chernoff parameter, Arb encloses the value of `z` and every
support-averaged transfer entry. The activation table supplies decimal upper
bounds. The live-state spectrum supplies exact integer multiplicities.

The one-tail recurrence averages matrix products over uniform subsets. Its
normalized form avoids a final division by a large binomial coefficient.

The two-tail recurrence first replaces each fixed support by independent
Bernoulli marks. Arb supplies a lower endpoint for each probability of the
conditioning event. Dividing by those lower endpoints gives an upper bound for
the original fixed-support experiment.

For the body classes, let `F[j]` be the transfer matrix conditioned on `j`
active positions in one region. The recurrence

```text
G[b+1,j] = (G[b,j] + G[b,j+1]) / 2
```

averages the `b` fair body bits. The body-only calculation uses `G[b,0]`. The
mixed calculation takes the entrywise maximum over `j`.

Each conditioned matrix carries its own signed power-of-two exponent. Matrix
addition aligns these exponents before rounding upward. Matrix multiplication
adds them exactly. This representation prevents underflow across the 1024-fold
region power.

## Outer change of measure

The body calculation uses the complete outer spectrum through the recorded
Holder change of measure. Each tested `p>1` is a fixed binary64 value and hence
a valid exact choice. Arb encloses its spectrum moment. Taking the minimum of
several upper bounds remains an upper bound.

For `m` available tail blocks and expected tail count `T`, the mixed calculation
uses

```text
(1+T)^m - 1 <= m T (1+T)^(m-1).
```

This positive expression avoids interval subtraction.

## Results

| Message class | Outward lower margin |
|---|---:|
| Exactly one tail block and no body block | 40.464902278356 |
| Exactly two tail blocks and no body block | 65.828399762063 |
| At least three tail blocks and no body block | 46.912124928575 |
| Body blocks and no tail block | 230.086343012485 |
| Body blocks and tail blocks | 237.522280635105 |
| Combined | 40.448462678054 |

The audit recomputes the final class sum and checks that every outward margin is
no larger than its nearest-binary64 diagnostic. It also checks that the retained
receipt used at least 256 Arb bits.

## Artifacts

- Analyzer: `scripts/certify_riffle_xbch_tail_outward.py`.
- Audit: `scripts/audit_riffle_xbch_outward_certificate.py`.
- Receipt: `receipts/outward_full_endpoint97_dba_xbch64_philips.json`.
- Audit receipt:
  `receipts/outward_full_endpoint97_dba_xbch64_philips_audit.json`.

This certificate closes the numerical arithmetic layer under the recorded
inner-transfer model. It does not replace the lemmas that justify the activation
and live-spectrum inputs.
