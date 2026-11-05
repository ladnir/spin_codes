# Retain the first two age cohorts jointly

This continues [the regional fresh-count analysis](QUARTET_REGIONAL_FRESH_TRANSFER.md).
The previous turn found that the chosen one-cohort upper transfer remains
too loose in the392..424 input-count band. Boundary-potential tuning alone
does not repair that numerical expression. This turn tests a local
constraint that the one-cohort transfer discards. The encoder is unchanged.

## What the omitted cohort does

Keep the five-type convention from the earlier notes. Type0 is an inactive
input coordinate. Types1,2,3,4 have active input and preceding free-run
length0,1,2,>=3. Write a for the active-input count, r=n1, and s=n2.
The former transfer records a,r but maximizes over all divisions of the
remaining a-r coordinates among types2,3,4.

That maximization can favor different ages at adjacent input counts. In
the selected local screen, low input counts benefit from assigning many
unresolved coordinates to type2; high counts can benefit from older types.
The actual transition cannot choose those assignments independently.

Fix the next active-input count l and the number k of positions in fully
active quartets. Let v_i count selected positions of incoming type i. Then

    n1' = l-k+v0,
    n2' = v1,
    n3' = v2,
    n4' = v3+v4.

In particular, a selected fresh coordinate must enter type2. It cannot
skip directly to an older type. Conversely, n2'=0 can be feasible; a proof
must count that event instead of excluding it without justification.

## Exact joint coefficient

Let w be the state width, and use the same positive product potential
phi(n)=product_i h_i^n_i, with h0=1. Set

    d=(h1,h2,h3,h4 b^(3lambda-tau),h4 b^lambda).

The number J_(l,k) is defined by the quartet polynomial in the preceding
notes. For a fixed incoming histogram n, define

    E_n(z)=product_(i=2..4)(1+d_i z)^n_i.

The potential-weighted transition coefficient to (l,r',s') is the sum,
over k,v0 with r'=l-k+v0 and v1=s', of

    J_(l,k) / binom(w,k)
    * binom(w-a,v0) binom(r,v1)
    * d0^v0 d1^v1 h1^(l-k)
    * [z^(k-v0-v1)] E_n(z) / phi(n).

All indices obey their binomial/support ranges. This is an identity for
the actual conditional transition, including its potential weight. It is
not an independent-coordinate approximation or a fixed-message certificate.

For fixed a,r,s, the old coefficient divided by its potential is at most,
for every z>0,

    z^(-(k-v0-v1))
    * ((1+d2 z)/h2)^s
    * max((1+d3 z)/h3,(1+d4 z)/h4)^(a-r-s).

Only the division between types3 and4 is relaxed here. The nonnegative
terms can be bounded separately and then summed. Auxiliary point-mass
factors from the previous note can sharpen the coefficient estimate;
their numerical evaluation still needs directed verification.

## A cheap two-update constraint, and its limitation

If the preceding active/fresh counts are a,r and the new fresh count is
r', then exactly u=l-r' active incoming positions were selected. Therefore

    max(0,l-r'-a+r) <= s' <= min(r,l-r').

These bounds describe a polygon in (r',s'). For a fresh-count bucket and
a second-age bucket, intersect the rectangle with

    l-a+r <= r'+s' <= l,   s'<=r,   r'>=l-a.

After a positive coefficient evaluation is fixed, the logarithm of the
potential factor is affine in these two counts. Its maximum occurs at a
polygon vertex. Take this maximum BEFORE choosing the evaluation point.
An empty polygon may be discarded only if its incoming coefficient is zero.
The implementation checks that condition explicitly.

The range constraint alone helps some two-update paths substantially but
does not settle the problem. For initial a=408,r=40 followed by input
counts424,424, it saves essentially nothing. It permits s'=0 without
accounting for the chance of selecting no fresh coordinates.

## Retain the selected-fresh count's multiplicity

The first joint transition can be factored without a three-dimensional
enumeration. Define

    A_(a,l,u) = sum_(k-v0=u) J_(l,k)/binom(w,k)
                       * binom(w-a,v0) d0^v0 h1^(l-k).

Then its coefficient for r'=l-u,s'=v1 is

    A_(a,l,u) binom(r,v1) d1^v1
    * [z^(u-v1)] E_n(z) / phi(n).

The first factor does not depend on v1. It is computed once per u. The
old-polynomial coefficient depends only on u-v1. This allows the joint
bucket sum to retain the binomial multiplicities without a large dense
two-cohort matrix. All contributions are positive; none are sampled or
dropped because they are small.

The current two-update diagnostic bounds the first transition uniformly
over its unresolved incoming old ages. It then sums the outgoing joint
buckets against the second transition's polygon bound. This preserves
the actual selected-fresh count in the bound; it does not merely substitute
its mean or its feasible range.

## Bounded full-width results

The following FLOAT results use w=512, the unchanged selected potential,
lambda=707/2560, tau=381/512, log2(b)=4, and bucket width8.
They fix the displayed input counts and initial fresh count, then sum all
intermediate joint buckets. The saving compares with the SAME first-step
bound followed by a second-step bound that discards the second-age count.

| Input counts (previous, first, second) | Initial fresh count | Joint two-update log2 upper | Saving |
| --- | ---: | ---: | ---: |
| 408,424,424 | 40 | 878.26548 | 1.11826 |
| 408,424,424 | 160 | 857.30189 | 6.58553 |
| 408,424,392 | 160 | 856.36470 | 5.96890 |

These are local counting factors, before the regional normalization and
whole-stream distance argument. Do not multiply the displayed savings by
the number of updates and call the result a security margin. A full
calculation must retain state constraints at every composition boundary.

The separate range-only test covers twelve selected two-update profiles.
Its savings range from essentially zero to9.29bits. The earlier fixed-s
screen is only a conditional-histogram comparison: it does not prove that
the favorable s value occurs. The joint multiplicity calculation above
addresses that distinction for the stated two-update tests.

## Exact checks and remaining work

`quartet_two_cohort.py` checks every width4 and width8 histogram at exact
toy weights. All13200 joint coefficients equal the independent five-type
kernel. All13200 joint Cauchy upper checks pass. The factored joint formula
equals the unfactored formula;565 unweighted histogram normalizations equal
2^w. Exhaustive polygon tests cover1349 nonempty and2476 empty boxes.

An initial range-screen run stopped on a genuinely unreachable fresh bucket.
The corrected code handles it only after confirming its incoming mass is
zero. The final range and joint-probability screens complete successfully.
Remote runs are serial, single-threaded,90s bounded, below64MiB, and file-free.
The joint screen and exact tests take about3.42s. No encoder timing,
construction change, whole-code certificate, or new failure allowance follows.

The regional extension below performs that composition. The remaining
count tails, covered counts, row counts, and directed arithmetic are still
required for the full45%/40-bit goal. Preserve all previous checkpoints.

## Uniform joint buckets and regional composition

The state label now records a and intervals R,S containing n1,n2. Retain
only boxes intersecting n1+n2<=a. For each exact r in R, intersect S with
0<=s<=a-r. This clipping is part of the incoming state, not a new condition
on the sampled construction.

For count=a-r, let Ebar_q dominate, uniformly over that clipped interval
and all n3+n4=count-s,

    [z^q] (1+d2 z)^s (1+d3 z)^n3 (1+d4 z)^n4
    / (h2^s h3^n3 h4^n4).

At each positive evaluation point, maximize over the endpoints of the
s interval and the two remaining old types BEFORE minimizing the Cauchy
bound. The point-mass correction uses a variance lower bound valid for
every histogram in the box. It does not use the variance of just one
maximizing histogram.

For outgoing intervals R',S', sum the factored coefficient over all
u,v satisfying l-u in R' and v in S':

    sum A_(a,l,u) binom(r,v) d1^v h1^(-r) Ebar_(u-v).

Then maximize this sum over exact r in the incoming R interval. These
maxima give a coefficientwise upper matrix between joint bucket labels.
The maximizing r may differ between outgoing labels; this is a relaxation,
not an independence claim. Induction therefore bounds every composition
of the true histogram kernel by the positive bucket matrix.

The implementation computes the inner sums as positive convolutions.
For each destination S', convolve the selected-fresh binomial coefficients
restricted to S' with Ebar. Reuse this convolution for every next input
count l. A coefficient tilt by zeta^q, zeta=1/4, stabilizes the arithmetic:
the factors zeta^v, zeta^(u-v), zeta^(-u) cancel exactly. Assertions stop
the FLOAT screen if a positive coefficient, convolution, or matrix entry
underflows. No small positive contributions are deliberately discarded.

Apply the sixteen-update regional counter DP from
[the preceding note](QUARTET_REGIONAL_FRESH_TRANSFER.md) to this matrix.
It carries both cohorts across regional boundaries. Only the counter for
the regional input total restarts. Use the actual initial state once and
the free terminal factor once. The payment marker still sums all regional
payment allocations at the fixed global count111572.

The bucket terminal factor maximizes the affine log weight on the vertices
of R x S intersected with r+s<=a. The remaining a-r-s coordinates use the
larger type3/type4 terminal ratio. This is again a uniform upper bound,
not a choice of a new terminal state for each region.

## Regional audit and restricted numerical comparison

`quartet_joint_regional_tests.py` compares the uniform bucket matrix with
the independent exact histogram kernel at widths4 and8. All8692 outgoing
bucket comparisons pass. A separate four-update test keeps two regional
input totals and the inherited state. All81 possible total pairs pass,
including the free terminal factor. Another1155 exact polynomial
coefficients lie below the FLOAT implementation's coefficient bounds at
toy widths; this last comparison is not directed full-width verification.

For the four-count domain {392,407,408,424}, the regional calculation sums
every count path, not just the previously selected words. It also sums
all regional payment allocations. At the same selected weights:

| Transfer | Bucket width | States | FLOAT log2 expression |
| --- | ---: | ---: | ---: |
| One age cohort | 8 | 207 | +2724.82424 |
| Two age cohorts | 32 | 378 | -2697.30787 |
| Two age cohorts | 16 | 1405 | -9577.53846 |

The one-cohort marked upper expression also has a positive FLOAT lower
diagnostic, +1631.90637. The joint-cohort result therefore addresses a
domain that defeats the previous local relaxation, even with much coarser
buckets. It is not a whole-code margin: other input counts, other covered
counts, other row counts, and directed evaluation are still missing.

The four-count joint screen takes about8.52s on Peach, serial and
single-threaded, with no files. Its proposed positive vector does not
fully converge; the bound uses the maximum row ratio, not an assumed
eigenvalue. The exact audit completes in0.19s. A lost tool observation was
replayed only after confirming that no remote computation remained live.

Extending to EVERY count in392..424, with width32 cohort buckets and3129
states, gives+4184.66545. This is better than the previous one-cohort
width4 value+7566.84, but does not close the full band. The marked upper
expression has a FLOAT lower diagnostic+604.29; its balanced-payment
subexpression has lower diagnostic-637.61. These diagnose the chosen
relaxation, not the actual bad-codeword probability.

The four-count comparison above isolates a substantial bucket-width effect:
halving both cohort buckets improves that same restricted expression by
6880.23bits. This cannot be extrapolated to the full band. It motivated
the separate full-band test below before retaining another age.

The full width32 run completes in289.74s. The width16 four-count run takes
20.93s. Caching old coefficients by their actual degree a-r reproduces the
latter expression exactly and reduces its runtime to16.22s. This is a
proof-calculation optimization, not an encoder performance claim.

The full-band driver has an explicit1200s watchdog and an8GiB conservative
array bound for width16. Smaller runs retain90s and768MiB bounds. The
Peach host had over100GiB available RAM before the fine full-band test.
All runs remain serial, single-threaded, and file-free; no local numeric
run is used while the local disk is below its required free-space limit.

## Full392..424 band with finer joint buckets

The width16 full-band run completes successfully in967.50s with11618
states. Its FLOAT log2 expression is **-3907.2296673407545**. This includes
every count path with all4096 input counts in392..424, and every regional
payment allocation at c=111572, for inactive-row count r=2100. Both cohorts
are inherited across all regional boundaries; the terminal state is free.

| Transfer over every count in392..424 | Bucket width | FLOAT log2 expression |
| --- | ---: | ---: |
| One age cohort | 4 | +7566.84 |
| Two age cohorts | 32 | +4184.66545 |
| Two age cohorts | 16 | -3907.22967 |

The fine run uses log2(z)=-0.0693359375 and
log2(y)=-3.9765229800373474. The maximum regional weighted log2 row norm
is-8.367172190867825. The proposed positive vector has residual0.02409;
as before, the calculation uses its maximum row ratio rather than assuming
an eigenvector. The worst boundary label has a=417,n1 in128..143,n2 in288..289.
Observed resident memory during the final computation was about3.1GiB,
within the explicit8GiB array budget. No files were written.

This clears the previously failing band numerically. It is not a directed
restricted certificate or a whole-code certificate. In particular, a path
that leaves392..424 even once is not covered by this result. The other
covered counts and outer-row counts are also open.

The next step is [coarse state cells with exact edge counts](QUARTET_LABELLED_EDGES.md).
That interface keeps the regional quota exact while allowing fewer state
cells outside the middle band. It includes an exact recurrence for routing
coefficients and a uniform coarse-count geometry inequality. Both have
small exact tests; the full-width tail matrix remains to be constructed.

After this full-band run, the routing evaluator was changed to the equivalent
recurrence. The four-count width16 result changes by less than1.4e-8 in
log2, consistent with floating-point evaluation order. The exact recurrence
is checked separately. The full-band value above was obtained before that
evaluation change; no second full-band run or directed validation is claimed.
