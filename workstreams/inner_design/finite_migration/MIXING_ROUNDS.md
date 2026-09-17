# IMT mixing strength and message length

The fixed-map study supports a length-dependent choice of mixing rounds.
One transvection per epoch is below the Q1 plateau at K=2^16, but nearly
reaches the same analysis's full-refresh limit at K=2^20. A separate dense
bound remains insensitive to mixing at the tested short-length bottlenecks.

All cells keep the BCH-derived [256,128] outer, t=128, s=19, weight-five
feedback, expansion map, routing distribution, and 10% distance threshold.
Only K and the number r of independent transvections per epoch change.
This is a proof-bound study, not an encoder benchmark or an estimate of
the true failure probability. It does not sweep new t or s values.

## Q1 results

Each cell is minus log2 of an outward upper bound on the one-active-row
contribution. It is not a full-distance margin.

| Transvections per epoch r | K=2^16 | K=2^18 | K=2^20 |
|---|---:|---:|---:|
| 1 | 44.255708 | 50.304186 | 50.076584 |
| 2 | 51.206637 | 51.838187 | 50.347349 |
| 3 | 53.054187 | 52.186953 | 50.411467 |
| 4 | 53.677353 | 52.301965 | 50.435387 |
| 8 | 54.106360 | 52.385324 | 50.454338 |
| Ideal full refresh | 54.128383 | 52.390185 | 50.455500 |

The six selected witnesses per cell emphasize low BCH weights. The earlier
K16 one-round injection-shell calculation retains a stronger 44.385342-bit
bound through a different witness bank; see `Q1_SLACK.md`. No earlier bound
is invalidated. The remaining rows differ from their wider binary64 screens
by less than 0.001 bits. The screen is used only to propose witnesses.

Two rounds recover about 6.82 bits over the best retained one-round K16 bound.
Three rounds recover about 8.67 bits. The benefit of a second round is much
smaller at K18 (1.534 bits) and K20 (0.271 bits).

At K16, three rounds remain 1.074 bits below the full-refresh calculation;
four rounds reduce that gap to 0.451 bits. At K18, two rounds are within
0.552 bits; at K20, one round is within 0.379 bits.
These gaps describe the selected sufficient bounds, not optimality among
all possible proofs or inner maps.

Three rounds restore the ordering in which smaller K gives a larger Q1
margin at these three lengths. In the full-refresh comparison, the decreases
are 1.738 bits from K16 to K18 and 1.935 bits from K18 to K20.
That is approximately 0.87 and 0.97 bits per doubling, respectively.
The full-refresh curve is a reference for these fixed maps and bounds,
not an absolute ceiling imposed by the BCH outer alone.

## What the parameters control

For one fresh transvection M and any fixed nonzero state q, the exact marginal is

    Mq ~ (1/2) delta_q + (1/2) Uniform(nonzero states).

Let P be this transition kernel on the m=2^s-1 nonzero states, and let J
be the uniform kernel. Then P=(I+J)/2 and J^2=J. For r independent rounds,

    P^r = 2^(-r) I + (1-2^(-r)) J.

The zero state remains zero. All r transvections are applied to the old
state before feedback is added:

    Y_i = X_i + A q_i,
    q_(i+1) = (M_(i,r) ... M_(i,1)) q_i + C X_i.

There are no intermediate outputs or feedback additions between these rounds.
The same sampled linear maps act on every message; the marginal argument
does not introduce fresh per-message setup randomness.

The roles are therefore distinct:

| Parameter | Direct effect |
|---|---|
| s | State space size; a fresh nonzero state has point mass 1/(2^s-1). Changing s also requires compatible maps. |
| r | Retention term 2^(-r) in each update, without changing the fixed maps or input partition. |
| t | Bits emitted per epoch; changing t changes update frequency, feedback grouping, and the maps. |
| K | Number of possible active outer rows and number of epochs available within each region. |

Increasing s alone does not reduce the 2^(-r) retention term. At t=128,
each transposed region has K/(128t) epochs: 4, 16, and 64 at the tested lengths.
More epochs provide more opportunities to mix between inputs. The low-weight
tail also depends on input locations and activation; r times the number of
epochs is not a universal sufficient statistic for the margin.

The ideal-refresh row sets the retention term to zero while keeping A and C.
Its marginal can be realized by a fresh uniform invertible linear map per
epoch. It is not a zero-round IMT and is not the previous inner implementation.

## Q1 transfer and checks

`mixing_rounds.py` generalizes the twelve-coordinate representation in
`Q1_SLACK.md`. Put epsilon=2^(-r) and eta=1-epsilon. Every lazy-branch factor
1/2 becomes epsilon; every fresh-branch factor 1/2 becomes eta.
Activation from zero is unchanged. For example, an arbitrary measure of
image weight v obeys

    H_v -> H_v: epsilon z^v                 (zero input),
    H_v -> C_w: eta z^v a_w / m             (zero input),
    H_v -> Z: epsilon c_v + eta f_v / m     (singleton input),
    H_v -> D: epsilon f_v                   (singleton input),
    H_v -> C_w: eta f_v a_w / m             (singleton input).

The definitions of f_v, c_v, and the represented measures are unchanged.
Uniform-shell rows use exact averaged singleton-cancellation weights;
the arbitrary-state row uses their maximum. Positive region and coefficient
products preserve the resulting bounds. Each outer weight may select its own
valid tilt witness before the same deterministic BCH spectrum bound is applied.

The binary64 screen uses a normalized positive coefficient recurrence.
The outward producer independently uses Arb and the existing unnormalized
recurrence. The replay uses 512-bit Arb, linear region products, and native
matrix multiplication for the coefficient layers. It checks every retained
coefficient and the exact rational BCH objective.

Exact small-state tests enumerate the actual transvection distribution and
its products. Additional tests check every represented source-law type for
zero and singleton inputs, including ideal refresh. The native coefficient
recurrence is cross-checked against the scalar recurrence, and the screen is
compared with an outward evaluation on the selected geometry.

## The distinct dense bottleneck

The three largest rectangles in the current K16 dense cover were reevaluated
with their existing witnesses, using only the coupled Bernoulli bound for
each new mixer law. No bound from the old law was used as a fallback.

| Rectangle | One round | Two rounds | Three rounds | Ideal refresh |
|---|---:|---:|---:|---:|
| `0011000000010010000` | 44.004934 | 44.006351 | 44.007060 | 44.007768 |
| `0011000000100000000` | 45.419526 | 45.420942 | 45.421650 | 45.422358 |
| `0010101001110100111` | 45.534943 | 45.536362 | 45.537072 | 45.537781 |

These are individual rectangle contributions, not dense totals or complete
certificates. The tiny differences persist under 512-bit replay.

To extend the existing Bernoulli argument, let F_v be the emitted moment
from image shell v and H_v its bound on a weighted feedback-syndrome point
probability, as defined in `../asymmetric/TRANSFER_ARGUMENT.md`.
Here H_v denotes that density bound, not the Q1 coordinate used above.
For any target state, the new bound is

    F_v (epsilon H_v + eta/m).

This follows by separating the lazy branch and the uniform-nonzero branch.
The zero-source row and the existing activation-density alternative remain
valid because the mixer fixes zero. This argument is implemented separately
in `mixing_dense_probe.py`; no one-round fixed-weight transfer is reused.

At these witnesses, the dense comparison is already nearly insensitive to
the retention term. The full short-length certificate therefore has two
different bottlenecks: limited Q1 mixing, and a dense bound that extra mixing
does not materially improve. Neither result measures the actual distance.

## Reproduction and next experiments

The local data are `MIXING_ROUNDS_SCREEN_v1.json`, `MIXING_ROUNDS_Q1_v1.json`,
and `MIXING_DENSE_PROBE_v1.json`. Keep these data separate from source commits.
The source hashes bind the exact maps and proof inputs. Existing paper
certificates, production encoders, and benchmark numbers are untouched.

```text
python -B workstreams/inner_design/finite_migration/mixing_rounds.py screen --output <fresh-screen.json>
python -B workstreams/inner_design/finite_migration/mixing_rounds.py certify --source <fresh-screen.json> --output <fresh-q1.json>
python -B workstreams/inner_design/finite_migration/replay_mixing_rounds.py --source <fresh-q1.json> --output <fresh-replay.json>
python -B workstreams/inner_design/finite_migration/mixing_dense_probe.py --output <fresh-dense-probe.json>
python -B workstreams/inner_design/finite_migration/mixing_dense_probe.py --output <fresh-dense-probe.json> --verify
python -B -m unittest discover -s workstreams/inner_design/finite_migration -p test_mixing_rounds.py
```

The scalar Q1 replay remains available through `mixing_rounds.py certify
--verify --output <existing-q1.json>`. All replay outputs require fresh names.

The next candidates are r=3 or 4 at K16, r=2 at K18, and r=1 at K20.
These choices come from Q1 proximity to the tested plateau, not full proofs
or an optimized timing comparison. Before selecting them:

1. Improve the dominant dense bounds and reconstruct coverage under each new law.
2. Check the sparse contributions and the complete rational union.
3. Benchmark the additional state updates with unchanged data movement.

Extra rounds do not repeat the BCH encoder, bit transpose, expansion, or
feedback calculation. They do add state parity/XOR operations and setup
storage. Their end-to-end cost must be measured, not inferred from round count.
