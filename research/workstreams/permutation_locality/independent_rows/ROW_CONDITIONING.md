# Condition on row weights, then average packet shapes

Independent row permutations provide more information than a group's union
support. A row of weight w is uniform over the binomial(256,w) supports.
The four rows are independent. This permits a different coefficient bound:
introduce independent Bernoulli bits for each row and condition only on
its total weight. The reference input then has an explicit packet-weight
distribution, so its inner operators can average shapes rather than
maximize over them.

`row_verify.py` now combines this reduction with outward, shape-resolved
eleven-coordinate operators. At q = 64 it bounds every message whose four
rows in each active group have weights in [64,192], with margin exceeding
79.92 bits. Rows may have different weights in that interval. This is a
restricted message class, not full coverage of occupancy 64.

The optional refinements in [REFINEMENTS.md](REFINEMENTS.md) raise this
restricted-class margin to 209.98 bits. The later density refinements
raise it to 403.96 bits; [DENSE_DIAGNOSIS.md](DENSE_DIAGNOSIS.md) records
their replay and why occupancy 80 still needs a stronger bound.

[REFINEMENTS.md](REFINEMENTS.md) records three later analysis refinements
that raise this same restricted-class margin to 209.98 bits. The baseline
results below are retained for comparison.

The earlier nine-coordinate binary64 pilot remains available for comparison.
Its middle-weight failures no longer describe the strongest current bound.

## The coefficient inequality

Fix q active four-row groups and their locations among 2048 groups. For
each selected group g, fix which of its rows are nonzero and their weights
w_(g,i). Its remaining rows are zero. Let A_w be the BCH weight enumerator;
the computation uses authenticated coefficientwise upper bounds Abar_w.

For each nonzero row, choose a witness 0 < p_(g,i) < 1. In a reference
experiment, replace that row by 256 independent Bernoulli(p_(g,i)) bits.
Leave zero rows identically zero. All reference bits, lane permutations,
region placements, and inner setup choices are independent as specified
by the construction. Define

    beta(w,p) = binomial(256,w) p^w (1-p)^(256-w).

Conditioning every reference row on its specified weight produces exactly
the original independent-row support distribution for the fixed message.
For any nonnegative function F of the encoded output,

    E[F | all row weights]
       <= E[F] / product_(g,i) beta(w_(g,i),p_(g,i)).

This is the elementary inequality E[F 1_event] <= E[F], divided by the
positive event probability. It does not assume that the conditioned bits
remain independent. Independence is used only in the unconditioned
reference experiment.

At output threshold H = 209715, set F = exp(-lambda weight(output)).
Multiplying the conditional moment by exp(lambda H) bounds the failure
probability for this fixed message. Multiplying by the number of messages
in the stated row-weight class gives a first-moment upper for that class.
The witness probabilities need not equal w/256. They select a coefficient
bound, not the actual construction distribution.

## Packet law and local operators

In the reference experiment, the packet-weight distribution in group g is
the Poisson-binomial law

    theta_(g,b) = [x^b] product_i (1-p_(g,i)+p_(g,i)x),   b = 0,...,4,

where fixed zero rows contribute factors 1. The independent lane
permutation makes the mask uniform among the binomial(4,b) masks of weight
b. Before conditioning on row weights, different regions use independent
reference bits and independent lane permutations.

The pilot restricts every active group to the same number r of nonzero
rows, the same nonzero weight w, and a common witness p. Then

    a = 1-(1-p)^r,
    phi_b = binomial(r,b) p^b (1-p)^(r-b) / a   for 1 <= b <= r,

with phi_b = 0 for b > r. Here a is the probability that a selected
group's packet is nonzero in a reference region. Conditional on a packet
being nonzero, phi is its weight distribution.

For j nonzero packets in one inner epoch, the reference input weights are
independent samples from phi. Let T_(n1,n2,n3,n4) be the existing local
operator for the specified weight multiplicities, including the average
over distinct physical windows. The averaged local operator is

    Tbar_j = sum_(n1+...+n4=j)
               j!/product_b(n_b!) product_b(phi_b^n_b) T_(n1,...,n4).

This is a nonnegative average of valid per-shape operators in the same
state-envelope coordinates. It averages the complete operators, including
their zero-return and output terms. It does not replace a conditioned
state distribution by a uniform one.

Let R_k be the region operator obtained by placing k nonzero packets
uniformly without replacement among its 2048 four-bit positions and
propagating the 64 inner epochs. The reference region operator is

    M = sum_(k=0)^q binomial(q,k) a^k (1-a)^(q-k) R_k.

Consequently the reference tilted moment is bounded by
e_zero M^256 terminal. The reference can have zero packets in selected
groups; conditioning the rows on their positive weights restores exactly
q active groups in the original experiment.

## Homogeneous row-weight class

Suppose every selected group has exactly r nonzero rows, each of weight
w > 0. Count all choices of the r row positions. The resulting bound is

    binomial(2048,q) (binomial(4,r) Abar_w^r)^q
      exp(lambda H) e_zero M^256 terminal / beta(w,p)^(rq).

Thus the conditioning denominator has rq factors, not one factor per
group. When all four rows are nonzero, it has 4q factors. A zero row needs
no conditioning factor because it is kept fixed in the reference.

The implementation currently tests 1 <= w <= 255. Weight 256 is not
excluded from the construction: a separate endpoint can set p = 1,
where beta(256,1) = 1. The all-zero four-row tuple belongs to an inactive
group and must not enter an active-group class.

## Selected diagnostic results

The following are optimized binary64 log2 contribution uppers at q = 64
and lambda = .032. They include group locations and row-position choices.
Negative values are not outward certificates.

| Nonzero rows per group r | Row weight w | Optimized witness p | Log2 contribution upper |
|---|---:|---:|---:|
| 1 | 64 | .740721 | -3616.13 |
| 1 | 96 | .759244 | -6851.82 |
| 1 | 128 | .794263 | -10520.04 |
| 2 | 64 | .508897 | -5202.26 |
| 2 | 96 | .561090 | -6006.07 |
| 2 | 128 | .625777 | -8223.62 |
| 4 | 48 | .317283 | -3608.43 |
| 4 | 64 | .351405 | -2179.14 |
| 4 | 80 | .389768 | +932.75 |
| 4 | 96 | .432167 | +3069.23 |
| 4 | 112 | .478340 | +4137.68 |
| 4 | 128 | .527926 | +4069.89 |

Witness optimization matters. At r = 1 and w = 96, using p = w/256 gives
log2 upper +6619.57. Optimizing p changes that diagnostic to -6851.82.
This is improvement within this coefficient bound, not measured code
performance or a distance estimate.

```sh
python -B research/workstreams/permutation_locality/independent_rows/row_weight_pilot.py --active-rows 4 --weights 48 64 80 96 112 128
python -B research/workstreams/permutation_locality/independent_rows/row_weight_pilot.py --active-rows 1 2 --weights 64 96 128
```

The default cache is `tmp/categorical-q64-pilot-checked.npz`. The driver
checks its source hash and tilt, authenticates the BCH caps, and passes
274 exact checks of the packet law and independent row conditioning.
The cache is diagnostic data, not a certificate input. The underlying
categorical model and its reconstruction are described in
[the earlier categorical pilot](../CATEGORICAL_PILOT.md).

These numbers are not directly comparable to the eleven-coordinate
support-shell bounds. The pilot omits mature-tail coordinates and the
later joint-cancellation/window-histogram refinements. It uses all-one
penalty rho = 1, not rho = .95. Its classes fix row weights but average
over their shuffled union supports; a fixed union-support class is a
different event. Its four-row middle-weight failures therefore neither
establish an obstruction nor demonstrate that averaging cannot help.

## Outward shape-averaged bounds

`shape_inner.py` constructs coupled eleven-coordinate operators for each
packet-weight shape through eight occupied windows. Each shape retains the
mature-state mass, density, and low-expansion-weight subset coordinates.
The construction includes the existing six-window feedback census,
eight-window output histograms, and three-window joint cancellation census.
The two-transvection conversion precedes the final shape average.

Above eight occupied windows, the implementation embeds the coarse
nine-coordinate universal bound in eleven coordinates. Each mature subset
receives the whole mature mass as an upper bound. This fallback is valid
but weaker than the strongest universal eleven-coordinate construction.

`row_verify.py` chooses reference probabilities in binary64, then recomputes
the selected bound using exact rational witnesses and outward Arb arithmetic.
The replay reconstructs the exact censuses; it does not consume the pilot
cache. The following margins include every group location and row-position
choice. Each row weight belongs to the indicated interval independently.

| q | Nonzero rows per group | Row-weight interval | Tilt | Reference p | Margin, rounded down |
|---:|---:|---|---:|---|---:|
| 64 | 4 | [80,80] | .048 | 433290547/1000000000 | 1435.309110 |
| 64 | 4 | [96,96] | .048 | 468806721/1000000000 | 267.127089 |
| 64 | 4 | [128,128] | .048 | 275056093/500000000 | 661.450722 |
| 64 | 4 | [112,112] | .052 | 102586993/200000000 | 159.432203 |
| 64 | 4 | [108,116] | .052 | 503217179/1000000000 | 111.808807 |
| 64 | 4 | [64,192] | .052 | 1/2 | 79.928346 |
| 64 | 4 | [80,176] | .052 | 1/2 | 108.608998 |

The last two rows were replayed at 384-bit precision. The [64,192] class
has first-moment upper less than 8.693011e-25. The event is that some
message in this class has output weight at most 209715. Probability is
over the ideal independent-row setup, including the routing and inner
randomness. This does not assert a distance bound for each sampled code.

The classes overlap, so these rows are alternatives, not a disjoint sum.
At tilt .048, the weight-112 bound was still positive (log2 upper 16.825).
Changing the proof tilt to .052 closes it without changing the encoder.
Larger tilts .064 and .08 weakened this bound.

```sh
python -B research/workstreams/permutation_locality/independent_rows/row_verify.py --groups 64 --active-rows 4 --intervals 64:192 80:176 --tilt .052 --p 1/2 --cut 8 --full-feedback 6 --window-histogram 8 --joint-cancellation --precision 384
python -B research/workstreams/permutation_locality/independent_rows/row_verify.py --groups 64 --active-rows 4 --intervals 112 108:116 --tilt .052 --cut 8 --full-feedback 6 --window-histogram 8 --joint-cancellation
```

The verifier passes 156 exact packet-law and thinning/placement checks.
Each cutoff-eight build passes 4356 atomic mixture checks and 726 direct
state inequalities, including selected three-window shapes. Exact census
totals and independent feedback/output marginals are also checked.
`test_row_counts.py` checks scalar domination and label multiplicities in
eight tests. An independent semantic review found no correctness blocker
in the local shape construction or the row-interval replay.

## Toward complete coverage

The remaining task is to cover row intervals and active-row counts that
differ across groups. The fixed-message coefficient inequality remains valid
with separate witnesses p_(g,i); the homogeneous transfer shortcut does
not. Those heterogeneous packet laws need a placement computation that
retains their types or a justified common upper.

One useful way to avoid summing a separate bound for every row weight is
pointwise domination. For a row-weight interval I and fixed p, define

    gamma_I(p) = max_(w in I) Abar_w / beta(w,p).

The averaged row counting measure restricted to I is pointwise at most
gamma_I(p) times the Bernoulli(p) reference distribution. Thus a box of
row-weight intervals can use a product of gamma factors. This bounds all
row-weight vectors in that box without pretending that differences of
cumulative upper bounds are shell counts. It still requires complete box
coverage and the correct heterogeneous inner operator. The outward replay
now implements the homogeneous-reference case; a complete heterogeneous
cover is still open. Occupancies 1 through 40 have separate full
support certificates in [LOW_OCCUPANCIES.md](LOW_OCCUPANCIES.md).

`typed_placement.py` implements exact placement for two fixed group types.
Its seven tests cover exhaustive labeled slot placements, noncommuting
operators, zero packets, and reduction to homogeneous placement. A type
is fixed across all 256 regions; averaging types independently in each
region would be a different reference distribution.

`typed_inner.py` forms the required epoch operators, with six additional
tests against direct packet and slot enumeration. `typed_verify.py` adds
the group-location, type-assignment, and row-position factors. At q = 64,
tilt .052, p = 1/2, and active-row weights in [64,192], its outward replay
gives the following margins at 192-bit precision:

| Four-active-row groups | Three-active-row groups | Margin, rounded down |
|---:|---:|---:|
| 63 | 1 | 137.259540 |
| 62 | 2 | 202.064506 |

Every permitted weight choice, group location, and zero-row position is
included. These classes are disjoint from the all-four-active class above.
The local fallback discards shape information above eight selected groups,
including when fewer than eight packets are actually nonzero. Its bound
is therefore potentially weaker than homogeneous thinning at high local
occupancy. Independent review found no counting or transfer defect.

```sh
python -B research/workstreams/permutation_locality/independent_rows/typed_verify.py --groups 64 --second-counts 1 2 --rows 4 3 --intervals 64:192 64:192 --probabilities 1/2 1/2 --tilt .052 --joint-cancellation
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows -p 'test_*.py'
```

[HOLDER_TYPES.md](HOLDER_TYPES.md) gives a proved product-measure envelope
for collections of persistent group types. It can reduce the number of
types needed by the placement calculation. The scalar diagnostic in
`holder_pilot.py` shows that merging every exception into one family loses
too much. Separating the number and orientation of exceptional rows retains
more information; these scalar costs are not distance margins.

Higher occupancies remain a separate issue. For example, at q = 80 the
same middle-interval bound is positive in a binary64 screen at tilt .052.
Neither selected 64-group classes nor a valid family reduction closes that
gap. The goal still requires all types at every occupancy through 2048.
