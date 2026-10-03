# What remains after the H2 route bound

The H2 gate disperses potential slots, but its weighted comparison can
still assign zero values to many of those slots. This distinction explains
why counting occupied potential steps alone did not close the remaining gap.

The diagnostic uses q=119 active outer groups, 32 regions, and 64 physical
steps per region. Each step has eight potential byte slots. Let j be the
number of selected potential slots in a step, and let k be the number
whose byte values are nonzero. These quantities are not interchangeable.

The gate applies the route marker exp(3.5 min(j,2)) and the fixed threshold
H2>=2559. The baseline16 point uses weight tilt theta=0.4; the actual24
point uses theta=0.5 and its authenticated saved local operators.

All expectations below refer to the normalized positive comparison path
sum, including that exponential marker. They describe neither the actual
encoder's failure distribution nor the exact law conditioned on H2>=2559.
The numerical results are floating diagnostics, not a certificate.

## Keeping the two occupancies separate

Write T_k for the existing local comparison operator at k nonzero bytes.
For p=255/256, the marked potential-slot operator is

\[
 Q_j=e^{3.5\min(j,2)}\sum_{k=0}^j\binom jk p^k(1-p)^{j-k}T_k.
\]

To count active occupancy k, the diagnostic multiplies T_k by exp(epsilon)
before this mixture. To count potential occupancy j, it instead multiplies
Q_j afterward. Both operations keep the H2 marker fixed.

The exact potential-slot placement has q selected positions in each
512-slot region. Chronological matrix multiplication preserves state
through all 2048 physical steps. Central differences of the log moment,
with epsilon=1e-4, estimate the corresponding expected counts.

| Quantity | Baseline16 | Actual24 |
|---|---:|---:|
| Good-message margin, bits | -617.493 | -258.391 |
| Potential packets | 3808 | 3808 |
| Nonzero byte packets | 2828.820 | 1834.078 |
| Steps with potential packets, H1 | 1523.104 | 1694.784 |
| Steps with nonzero byte packets | 788.013 | 364.658 |
| H2 | 2522.063 | 2683.692 |
| Zero-source steps | 1390.576 | 1705.430 |
| Nonzero-source steps | 657.424 | 342.570 |

Thus the actual24 comparison has about 1974 zero-valued packets among
3808 potential packets. Potential dispersion does not imply comparable
dispersion of nonzero bytes.

The complete expected histograms, indexed by occupancy 0 through 8, are:

```text
baseline16 potential: 524.90 524.14 403.27 203.99 187.33 128.98 58.34 15.29  1.76
baseline16 active:   1259.99 19.37 178.52 202.08 186.36 127.67 57.37 14.92  1.70
actual24 potential:   353.22 705.88 635.05 44.01  66.65  95.58 88.94 47.55 11.13
actual24 active:     1683.34   2.27  18.74 35.17  67.24  95.72 88.13 46.60 10.79
```

## Which state transitions carry the weight

Source-row markers distinguish zero from nonzero state. More selective
markers distinguish empty zero-to-zero steps, occupied zero-feedback
steps, births, and returns.

| Transition or source class | Baseline16 | Actual24 |
|---|---:|---:|
| Empty input, zero to zero | 1259.969 | 1683.341 |
| Nonzero input, zero feedback, zero to zero | 8.551 | 0.0586 |
| Zero to a birth coordinate | 122.056 | 22.0303 |
| Nonzero coordinate to zero | 121.056 | 21.0304 |
| Uniform-nonzero source coordinate | 535.516 | 320.550 |
| Birth-family source coordinate | 121.909 | 22.0195 |

The return row is the upper operator's filled-in return coefficient. Its
count must not be called an exact count of feasible encoder returns.
The zero-feedback count, by contrast, marks the weighted zero-source
entry T_k(0,0) for k>0. That entry uses the feedback census directly.

Occupied zero-feedback cycles are too rare to explain the actual24 gap.
Its zero-state runs are protected mainly by zero-valued bytes, despite
the presence of potential packets.

The H1 diagnostic also explains the failed conjunction at the tested
H2 multiplier. Adding a multiplier nu1 changes the log bound with
derivative E[H1]-1483. This derivative is positive at nu1=0 for both
points above. Since the log moment is convex in nu1, increasing nu1
cannot improve the bound while the other parameters remain fixed.
This statement does not optimize jointly over the other parameters.

## Numerical and provenance checks

`conditioned_diagnostics.py` batches independent positive perturbations
of the ordered placement recurrence. Its scaled arithmetic can discard
subnormal contributions. Each reported case therefore checks its baseline
and six perturbations against the frozen all-logarithmic evaluator.
The largest observed log-moment difference was 6.4e-12.

Potential and active step counts each sum to 2048. The weighted potential
count equals 32q=3808. Source-state partitions and birth-minus-return
terminal flow also agree. The largest conservation residual was 3.4e-6.

Two portable tests compare batched placement with the log implementation
and validate the marked derivatives with independent block-triangular
matrix transfers. `conditioned_diagnostics_v1.json` pins all sources and
the saved 24-bit receipt. Both cases took about 47 seconds together;
no new large-state census, encoder benchmark, or remote computation ran.

## Next proof direction

H3=sum min(j,3) is a better-supported next gate than another H1 constraint:
H2 stops distinguishing a step once it has two potential packets, whereas
the residual actual24 expression has many one- and two-packet steps with
zero byte values. A third-cap constraint probes the remaining tradeoff
between zero labels and densely packed nonzero bytes.

The separate `capped_gate.py` experiment subsequently obtained a positive
single-occupancy actual24 gate at q119, theta=0.5, and nu3=4. Independent
review verified its strict route threshold, potential-slot marker, and
fixed-route expectation argument. Re-evaluation gave 168.231009 good-message
bits, agreeing with its receipt within 2.8e-12 log units. Its exact route
bound passes 61 bits. This is one floating occupancy gate, not a full-q
outward certificate.

For a cheaper 16-bit construction, the birth count suggests a distinct
target: strengthen what is emitted when zero state first receives input.
A feedforward layer I+AC can address that event without enlarging the
state. It requires new joint local moments; the diagnostic does not
establish that it improves the full bound. See `FEEDFORWARD.md`.
