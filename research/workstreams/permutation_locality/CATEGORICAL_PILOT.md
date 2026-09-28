# Histogram averaging and a restricted dense certificate

The histogram-averaged transfer did not improve the difficult balanced
64-group point. A separate information-set argument does certify a class
of dense messages. Neither result closes the complete two-update goal.

## Averaging pilot

`categorical_model.py` retains the nine mass, fresh-state, density, and
uniform-class coordinates. It builds a separate local operator for each
input-weight shape. The operator includes the existing window averages,
multi-window moments, fresh collisions, and full feedback distributions
through six windows. It deliberately omits the two mature-tail coordinates.
All numerical transfer arithmetic in this pilot is binary64.

For categorical probabilities theta on weights 1,2,3,4, it averages those
operators with their multinomial shape probabilities. It then performs
the existing without-replacement epoch placement. The region operator
averages its active-group count using an auxiliary Bernoulli probability p.
The final bound divides by the probability of each fixed group histogram,
as justified in ODD_COLUMNS.md. The optimization varies p.

The matched maximum baseline uses these same per-shape operators and
outer count, but takes entrywise shape maxima. That baseline only needs
the Bernoulli support-conditioning factor; it does not pay for conditioning
on the finer histogram. It is not the stronger eleven-coordinate verifier.

The following diagnostics use 64 groups, union support 128 in every group,
lambda = 0.032, and two updates. A profile lists the numbers of columns of
weights 1,2,3,4; the other 128 columns are zero. Counts are integer uppers
intersecting the support, XOR-weight, total-weight, and all-one-column bounds.

| Nonzero-column profile | Log2 outer count | Averaged bound | Matched maximum bound |
|---|---:|---:|---:|
| (32,56,32,8) | 341.2959 | +3821.0744 | +3374.9497 |
| (64,0,64,0) | 274.6728 | -826.6036 | -888.9292 |
| (0,128,0,0) | 278.6920 | -828.1258 | -631.7006 |
| (96,0,32,0) | 243.7910 | -2805.1057 | -2865.3641 |
| (32,80,16,0) | 303.9595 | +1232.4310 | +985.4176 |

Negative entries are selected-profile diagnostics, not certificates. In
particular, their signs do not establish a new complete occupancy result.
For the balanced profile, optimizing theta converges after 78 evaluations
to approximately (0.249456,0.437881,0.250383,0.062280). The bound improves
only to +3821.0576. The histogram-conditioning cost outweighs the averaging
gain in this test. A larger theta sweep is not the next priority.

Checks include exact multinomial normalization through local degree 32,
all four atomic categorical limits, and 270 direct empty/single-window
state inequalities. Every per-shape coefficient is below the existing
nine-coordinate envelope, within numerical tolerance. That comparison
checks consistency, not outward validity of the pilot.

    python -B research/workstreams/permutation_locality/categorical_screen.py --cache tmp/categorical-q64-pilot-checked.npz --profile 32 56 32 8 --optimize-theta
    python -B research/workstreams/permutation_locality/categorical_screen.py --load-cache tmp/categorical-q64-pilot-checked.npz

The optional NumPy cache contains diagnostic arrays and dependency hashes.
It is not a certificate input and stays outside version control.

## Dense messages with bounded histogram-conditioning cost

This argument applies to the specified ideal one-column grouped route.
It requires its independent uniform coordinate and lane permutations.
It does not require a transvection distribution: it holds for every fixed
linear inner of the form y_j = x_j + E q_j, with q_j determined by earlier
inputs and fixed setup. In particular, it applies to the two-update target.
The full length-N inner is invertible because these output equations are
block triangular with identity diagonal. Here N = 2^21 and the initial
state is zero.

Fix a nonzero message with exactly q active four-row outer groups. For
active group i, let n_i,w count its columns of weight w, for w = 0,...,4.
Define

    P_i = 256! product_w binomial(4,w)^n_i,w
          / (16^256 product_w n_i,w!).

This is the probability of its histogram under 256 independent uniform
four-bit columns. Restrict attention to messages satisfying

    product_i P_i >= 2^(-bq).

Equivalently, their total histogram-conditioning cost is at most bq bits.
This is an aggregate condition; individual groups need not each have cost
at most b. The restriction is a property of the message's fixed outer
words, independent of the sampled SPIN setup.

Fix the region group placements and the inner setup. Replace the inputs
of these q groups by independent uniform four-bit columns, leaving all
other groups zero. The resulting N-bit input is uniform on a coordinate
subspace of dimension D = 1024q. The inner sends it to a D-dimensional
binary linear subspace. Some D output coordinates form an information set,
so they are independent unbiased bits. For 0 < z <= 1, dropping the other
output coordinates gives

    E[z^output_weight] <= ((1+z)/2)^D.

Condition the reference input on all the fixed group histograms. Its
distribution then equals the original shared column shuffle followed by
independent lane shuffles, for this fixed message. Nonnegativity therefore
gives

    E[z^output_weight | histograms] <= 2^(bq) ((1+z)/2)^D.

The bound is independent of the fixed group placements and inner setup,
so averaging those choices preserves it. At most binomial(2048,q) 2^(512q)
messages have q active groups. For integer threshold T = 209715, a union
bound over the restricted messages gives

    U_q = binomial(2048,q) 2^((512+b)q)
          z^(-T) ((1+z)/2)^(1024q).

The implementation uses the exact rational witness z = T/(1024q-T), which
lies in (0,1) in the reported range. It evaluates and sums these bounds
with outward Arb arithmetic. No BCH spectrum approximation enters this
argument; it uses the 512 message bits per group and the actual histogram.

## Verified restricted ranges

At both 192 and 384 bits of precision, b = 20 gives

    sum_{q=1978}^{2048} U_q < 8.770366e-86,

or more than 282.553 bits of margin. Occupancy 1977 is not included: its
bound has log2 value +52.7323. The histogram (16,64,96,64,16) has conditioning
cost about 14.6150 bits, so repeating this histogram satisfies the restriction.

At 384-bit precision, the looser aggregate cost b = 30 gives

    sum_{q=2038}^{2048} U_q < 3.773151e-71,

or more than 233.941 bits of margin. The two restricted classes may overlap;
adding their upper bounds still bounds their union. Neither class includes
all messages at its occupancies. Atypical histograms and intermediate
occupancies remain open, as do the required small-occupancy two-update replays.

    python -B research/workstreams/permutation_locality/dense_histogram_certificate.py
    python -B research/workstreams/permutation_locality/dense_histogram_certificate.py --precision 384
    python -B research/workstreams/permutation_locality/dense_histogram_certificate.py --bits 30 --precision 384

The verifier also checks 71 exact information-set and multinomial examples,
and exact histogram-eligibility examples. The proof above, not those small
tests, establishes the universal information-set inequality.

The next substantive task is to bound the complement of this dense class:
messages with unusually concentrated column histograms. They require their
smaller outer multiplicities and an inner bound that uses the concentration.
The intermediate balanced profiles still need a stronger argument; the
negative averaging experiment should not be promoted into a full cover.
