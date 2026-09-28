# Covering mixtures of BCH row weights

Sparse support covers retain exact regional placement, but maximize the
inner bound over packet shapes. At larger occupancies that maximum loses
too much information. This method keeps a positive mixture of input
distributions and bounds all mixtures through two coordinates.

The construction remains the two-bit ensemble in [README.md](README.md).
All comparison measures below are proof devices. The actual row shuffles,
regional permutations, lane swaps, and transvections remain independent
setup choices, sampled once and fixed for all messages.

## A positive bound for one shuffled row

Let A_w be the number of BCH words of weight w, and let a_w be its
authenticated upper bound. Averaging a row's coordinate permutation
assigns counting mass A_w/binomial(256,w) to each particular word of
weight w. Choose rational C>0 and 0<theta<1/2. Define

    T = max(0, max_{1<=w<=255}
        (a_w/binomial(256,w) - C/2^256)
        / (theta^w (1-theta)^(256-w)
           + (1-theta)^w theta^(256-w))).

The averaged nonzero-row measure is pointwise bounded by the sum of four
product measures: the deterministic all-one row with coefficient 1,
Bernoulli(1/2) with coefficient C, and Bernoulli(theta) and
Bernoulli(1-theta), each with coefficient T. Keep the zero row separately
with coefficient 1. `row_mixture.py` checks all 257 shells using exact
rational arithmetic. The current choice C=2^130, theta=1/3 gives
log2(T) approximately 82.16812408847.

Multiplying two such row measures gives 15 unordered pair components.
For row bit probabilities p,r, independent lane swaps give categorical
packet probabilities

    pi = ((1-p)(1-r), p(1-r)+(1-p)r, pr)

for zero, single, and double packets. A single packet's orientation is
uniform. The component coefficient is the product of row coefficients,
multiplied by two for distinct row components.

Only the zero/zero component has inactive label 0. Every other component
has active label 1, even when its comparison product distribution produces
an all-zero pair. Consequently an original message with q active pairs
is bounded by components with exactly q active labels. The comparison
does not confuse actual support with this label count.

## Shuffling heterogeneous independent packets

Consider n independent categorical variables with distributions pi_i,
followed by a uniform permutation of their positions. Write bar_pi for
the average distribution. For nonnegative category counts a_j summing
to n, positivity of generating-function coefficients and AM--GM imply

    Pr[counts=a] <= inf_{y_j>0} product_i (pi_i dot y) / product_j y_j^a_j
                 <= inf_{y_j>0} (bar_pi dot y)^n / product_j y_j^a_j.

Taking y_j=a_j/(n bar_pi_j), with limits at zero coordinates, and dividing
by the iid multinomial probability gives

    Pr_heterogeneous[counts=a] / Pr_iid[counts=a]
       <= n^n product_j a_j! / (n! product_j a_j^a_j).

Coordinates absent from bar_pi are absent from every pi_i and may be
discarded. Use 0^0=1. The largest right-hand side occurs at balanced
integer counts. To see this, put g(a)=a!/a^a. The ratio
g(a+1)/g(a)=(a/(a+1))^a decreases with a, including the limiting value
1 at a=0. Transferring one unit from a larger count to a smaller count
therefore increases the product of g values until the counts are balanced.

Let R_n denote this balanced bound for three categories. Conditional on
the counts, both shuffled laws are uniform over their category strings.
Thus the same ratio bounds their probability masses pointwise. Single
packet orientations are independently uniform under both distributions.
For n=4096, log2(R_n) is approximately 12.27428706309. Pay this loss in
each of the 256 independently permuted regions.

## Positive tilts retain only two mixture coordinates

Fix positive category weights z=(1,z1,z2). Component i has coefficient
c_i, categorical distribution pi_i, and active label b_i. Define

    Z_i = pi_i dot z,
    a_ij = pi_ij z_j / Z_i,
    D_i = c_i Z_i^256.

For component counts n_i summing to G=4096, put

    m_j = sum_i (n_i/G) a_ij,
    w_j = m_j/z_j,    S = sum_j w_j.

Applying the preceding density bound to the tilted component laws gives
a comparison with iid packet law w/S. Its multiplier in each region is

    R_G S^G product_i Z_i^n_i.

For a fixed category string, tilting multiplies its probability by
product_j z_j^a_j / product_i Z_i^n_i. Undoing this factor gives the
displayed multiplier; no conditioning on a typical event is used.
Across all regions the normalization is S^(G*256), once per **packet**,
not once per bit. The products of Z_i have already entered D_i.

Since m0+m1+m2=1, only (m1,m2) is needed by the inner bound. This projection
removes the need to enumerate all 15 component counts.

## Summing a complete cell

Fix a rectangle C=[x0,x1] times [y0,y1] for (m1,m2). Each comparison input
in this cell is pointwise bounded by the unnormalized categorical weights

    w_upper = (1-x0-y0, x1/z1, y1/z2).

Let S_upper be their sum. Normalize them to a packet distribution and use
the two-state envelope K(lambda) from [PROFILE_MEASURE.md](PROFILE_MEASURE.md).
For every entering inner state, positivity permits the unnormalized
comparison, including when S_upper is less than one.

It remains to sum the outer weights of all component assignments in C
with at least q_min active labels. For arbitrary real eta1,eta2 and mu>=0,
define

    P = sum_i D_i exp(eta1*a_i1 + eta2*a_i2 + mu*b_i),
    h_C = min_{(x,y) in C} (eta1*x + eta2*y).

The multinomial theorem bounds that entire sum by

    P^G exp(-G*h_C - mu*q_min).

For each retained assignment the exponential marking factor is at least
exp(G*h_C+mu*q_min). Dropping the restrictions after inserting this
factor gives P^G. This includes all assignments and their multinomial
multiplicities; it is not a bound on only one sampled composition.

Put d=209715. A valid first-moment bound for the whole cell is therefore

    exp(lambda*d) R_G^256 S_upper^(G*256)
      * e_zero K(lambda)^16384 (1,1)^T
      * P^G exp(-G*h_C - mu*q_min).

Every cell may use its own rational lambda>0 and dual witness. The cover
adds the resulting bounds; it does not interpolate optimized witnesses.

The driver can also change the input tilt inside a cell without changing
the coordinates that define its partition. Let z' be the alternative tilt
and put Z'_i=pi_i dot z'. For each category j, its unnormalized comparison
weight is the component average of v_ij=pi_ij/Z'_i. To bound that average,
choose real eta1,eta2 and mu>=0, and define

    c0 = max_i (v_ij - eta1*a_i1 - eta2*a_i2 + mu*b_i).

Throughout the cell with at least q_min active labels, a valid upper is

    c0 + max_{(x,y) in C}(eta1*x+eta2*y) - mu*q_min/G.

Apply this independently to each of the three categorical weights.
For the outer sum, replace D_i by c_i*(Z'_i)^256, but retain the original
features a_i that define the cell. The preceding first-moment formula then
holds with these new coefficients and the three upper weights. Numerical
linear programs only propose rational dual witnesses: replay checks the
finite maximum defining c0 exactly and never trusts a solver's optimum.

Near a boundary, the driver also reduces the universal density loss.
Let epsilon_j be the smallest positive feature a_ij. A component capable
of producing category j contributes at least epsilon_j to that feature.
For a cell with coordinate upper m_j^+, at most
floor(G*m_j^+/epsilon_j) slots can produce that category. Thus its realized
count has the same cap, under every positive tilt. Replace R_G by the
maximum of its factorial expression over these capped integer counts.
Log-concavity of a!/a^a makes this maximum a capped balancing problem,
which `capped_density_loss` solves exactly. When double packets are
impossible, it recovers the two-category loss instead of paying for a
third category. The uniform three-category bound remains valid elsewhere.

The feasible projected domain for q>=q_min is the convex hull of
{a_i : b_i=1} and {(q_min/G)*a_i : b_i=1}. Indeed, the inactive component
projects to zero, and the total active fraction lies in [q_min/G,1].
`mixture_cover.py` prunes rectangles only by exact separating inequalities.
Boundary cells remain covered. Dyadic bisection constructs an exact
partition of the root square; shared boundaries may be counted twice.

## Implementation and current scope

`mixture_probe.py` optimizes witnesses for selected component compositions.
At q=1024, the component consisting of a zero row and a Bernoulli(1/3)
row has an outward margin above 8852 bits, including location multiplicity.
This selected result alone says nothing about other mixtures.

`mixture_cover.py` searches the complete projected domain, replays accepted
cells with 192-bit Arb arithmetic, and saves rational witnesses. Its
`--replay` checks the whole partition and recomputes each bound. Stored
upper bounds are not trusted. The verifier rejects incomplete covers.

```sh
python -B research/workstreams/permutation_locality/two_bit/mixture_cover.py --minimum-groups 1024 --max-cells 5000 --output tmp/two-bit-mixture-dense-1024.json
python -B research/workstreams/permutation_locality/two_bit/mixture_cover.py --replay tmp/two-bit-mixture-dense-1024.json
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

Generated witnesses belong under ignored `tmp/`, not in version control.
The first command is a proof attempt; completion must be checked before
using its output as a certificate. Full-code closure also requires every
occupancy below q_min and the final outward sum.

The initial fixed-tilt runs did not close. Changing the relative weights
of single and double packets repairs individual failing cells, but the
adaptive two-coordinate covers at q_min=1024, 2048, and 3072 still leave
cells unresolved. These partial-cell sums must not be reported as bounds
for the entire dense range. A continuous alternative-tilt search repaired
one additional boundary cell, with an outward log2 upper below -90.89478;
that selected repair does not complete any of these covers.

One remaining loss is explicit: alternative tilts bound each categorical
weight separately, then maximize the outer count separately. These maxima
can correspond to different component compositions. For example, many
central/central pairs have large counting mass, whereas deterministic
single-packet pairs can maximize a categorical weight. The present
projection need not distinguish those two explanations of a mean.
The [lifted cover](LIFTED_MIXTURE.md) now retains an additional coordinate
for the central-component fraction. It covers the full composition domain;
a selected high-entropy composition is not a substitute for that cover.
Its first complete vertical slice closes, but the full-domain run must
finish and replay before it establishes a dense-range certificate.

Small exhaustive tests check the heterogeneous density inequality, its
tilted form, the multinomial cell sum, projection pruning, normalization
per packet, and rejection of incomplete or overlapping saved partitions.
The old four-bit proof and production encoder are unchanged.
