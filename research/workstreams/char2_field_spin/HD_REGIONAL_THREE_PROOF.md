# Coordinate-major three-scan codes have distance below 40%

The proposed K20 construction has distance at most 699169, approximately
33.339%. This holds for every choice of nonzero coefficients and every
regional permutation. The stronger proof below also permits an arbitrary
global permutation between the two pre-scans and the final scan.

Thus the proposed construction cannot meet the 40% distance target. A sharper
enumerator or larger field cannot repair this topology. This is a distance
upper bound, not a failed attempt to prove a distance lower bound.

## The exact map and its parity checks

Let the field be F=GF65536. The proofs also hold over any field. Apply a
full-rank linear [b,b/2] outer independently to R source rows. Write
rho=b/2, K=rho*R, and N=b*R. The current parameters are b=64 and R=32768.
The outer is GRS, but these obstructions need only its dimension.

Store the outer coordinates in coordinate-major order. Coordinate c of
source row r occupies position c*R+r. Each of the b coordinate regions
therefore presents the source rows in the same order.

An input-weighted prefix scan A maps x to y by

    y_i = y_(i-1) + w_i*x_i,    y_(-1)=0,

where every w_i is nonzero. All N outputs are retained. Its inverse D is
lower bidiagonal: column j is supported on j and, when present, j+1.
Two consecutive scans have an inverse supported on j,j+1,j+2 in column j.
All scans carry across region boundaries; they do not reset there.

Let B be the ordered outer encoder, and let P act on coordinates. The
proposed encoder and its parity-check matrix are

    E = A_3 P A_2 A_1 B,
    H = H_0 D_1 D_2 P^(-1) D_3.

Here H_0 has R row groups, each containing rho checks for one source row.
Each source-coordinate column of H_0 touches only its corresponding group.
Since all scans and P are invertible, ker(H) is exactly the output code.

More generally, let p be the number of consecutive pre-scans and q the
number of consecutive post-scans. Their inverse bandwidths are at most p
and q. Thus each post-inverse column reaches at most q+1 pre-inverse columns.
Each pre-inverse column reaches at most p+1 consecutive source positions.
Cancellations can shrink these supports but cannot enlarge them.

The current construction has p=2 and q=1. Its setup samples independent
nonzero scan coefficients and independent uniform permutations of the R
positions in each region. The deterministic results below do not require
either sampling assumption. Shared 128-bit records inherit the obstruction:
place a scalar bad codeword in one field lane and set the other lanes to zero.

## A fixed check prefix rejects every regional permutation

Fix a prefix U of a source-row groups, with p<a<=R. It contains rho*a
parity-check rows. Each region contains A=a-p pre-inverse columns whose
indices have row offset less than a-p. Their entire source support lies in U.
Call these columns safe. A regional permutation preserves their count in
each region, although it changes their locations.

First take q=1. In a binary string of length R with A safe positions, the
number of adjacent safe pairs is at least

    A - (R-A+1) = 2*A-R-1.

There are at most R-A+1 runs of safe positions; each run of length l
contributes l-1 pairs. Every within-region safe pair gives one column of H
confined to U. We discard every post-inverse column crossing a region
boundary. We also discard the terminal singleton. Therefore global carry
does not invalidate the count.

Across all regions, at least b*(2*a-2*p-R-1) columns are confined to rho*a
checks. If

    3*a > 2*R+4*p+2,

then at least rho*a+1 confined columns exist. They are linearly dependent,
so a nonzero codeword has weight at most rho*a+1.

For p=2, R=32768, and b=64, choose a=21849. Each region has 21847 safe
positions and at least 10925 adjacent safe pairs. Hence at least 699200
columns lie in only 699168 checks, proving distance at most 699169.

This argument uses neither coefficients nor an average over permutations.
It applies separately to every setup in the proposed ensemble.

For general q, count within-region windows of length q+1. Such a binary
string contains at least

    max(0, (q+1)*A-q*R-q)

wholly safe windows. To see this, start with R-q windows and charge at most
q+1 intersecting windows to each unsafe position. Thus the same dependence
argument works whenever

    (2*q+1)*a > 2*q*R + 2*(q+1)*p + 2*q.

## A cyclic interval also rejects arbitrary global routing

The preceding fixed-prefix proof uses regional routing. We now remove that
restriction. Fix every coefficient and any permutation P of all N positions.
The only randomness in this proof is an auxiliary choice of check groups.

Sample the start of a cyclic interval U of a row groups from {0,...,R-1}.
Row-group indices are taken modulo R; actual scan indices are not wrapped.
A pre-inverse column touches at most p+1 consecutive row groups on this
cycle. This remains true when its support crosses a coordinate-region
boundary. Its support is wholly inside U with probability at least
(a-p)/R. Truncation at the end of the complete scan only improves this bound.

A final column reaches at most q+1 such clusters through P^(-1). A union
bound on their non-containment probabilities gives

    Pr_U[column confined to U] >= (q+1)*(a-p)/R-q.

No independence between these clusters is assumed. The bound may be
negative; only its positive regime will be used. By linearity, the expected
number of confined columns is at least

    b*((q+1)*(a-p)-q*R).

If this expectation exceeds rho*a, some choice of U confines more columns
than checks. Consequently, the following finite condition is sufficient:

    (2*q+1)*a > 2*q*R + 2*(q+1)*p.

Take the smallest integer satisfying this inequality, provided it is at
most R. The resulting distance upper bound is rho*a+1. At the current
parameters it again gives a=21849 and distance at most 699169. Its expected
column lower bound is 699264, slightly stronger than the fixed-prefix count.

For fixed p and large R, the relative ceiling approaches q/(2*q+1).
Increasing only the pre-order leaves the scalar-post ceiling near 1/3.
This rules out replacing the regional permutation with a global permutation
as a repair for the current common coordinate-major ordering.

## Uniform regional routing also limits wider post-scans

The deterministic ceiling approaches 40% for q=2 and 3/7 for q=3. It does
not settle whether these wider maps offer a 40-bit setup guarantee. A second
exact argument rejects that guarantee for q=2 and q=3 at the current size.

Fix the prefix U containing a=26214 row groups. Set p=2, so each region
contains A=26212 designated safe pre-inverse columns. Let Z count the
within-region windows of length q+1 whose routed labels are all safe.
Each counted window again gives a column of H confined to U.

For a uniform regional permutation, a fixed window receives q+1 distinct
uniform labels from its region. Therefore

    mu := E_P[Z] = b*(R-q)*C(A,q+1)/C(R,q+1).

The coefficient assignment is arbitrary and fixed. Independence between
regions is unnecessary for this expectation. All windows together satisfy
Z<=M, where M=b*max(0,A-q): a string with A safe positions has at most
A-q wholly safe windows when A>=q+1.

Write t=rho*a=838848. If Z>t, the code has distance at most t+1=838849,
which is below the required 838861. Splitting the expectation on this event
gives the lower bound

    Pr_P[distance<=838849] >= Pr_P[Z>t]
                           >= max(0,(mu-t)/(M-t)),

when M>t. This is a lower bound on setup failure, not a first-moment upper
bound on success. It is unaffected by the size of the coefficient field.

The helper computes mu and the failure bounds as exact rational numbers.
The following displays round those exact results:

| Post order q | Deterministic cyclic distance ceiling | mu | Setup-failure lower bound |
|---|---:|---:|---:|
| 1 | 699169 | 1341880.335938 | 0.5998077113 |
| 2 | 838945 | 1073356.840874 | 0.2796459314 |
| 3 | 898881 | 858560.991347 | 0.02350904364 |
| 4 | 932161 | 686743.978673 | 0; inconclusive |

For q=1, deterministic failure is already one; the table's weaker
probabilistic value is included only to check the common formula. For q=3,
setup reliability is at most about 5.411 bits under this necessary event,
not 40 bits. These probabilities do not rule out selecting and certifying
a rare good setup. They do rule out a 2^-40 failure claim for the stated
unconditioned random ensemble. Four post-scans merely survive this test;
no distance lower bound is established for them.

## Distinct row orders must break the common blocks

Suppose coordinate c presents source rows in a fixed order pi_c instead of
the common order. A pre-inverse column then touches at most p+1 consecutive
entries of pi_c. Different orders can destroy the common cyclic-interval
property used above. Thus the deterministic 33.34% ceiling does not extend
to arbitrary distinct row orders by the same argument.

There is a performance-relevant restriction, however. Partition the source
rows into common blocks of L rows, where L divides R. Suppose every pi_c
merely permutes these blocks, with arbitrary orders inside each block.
Select m whole source blocks as U, so a=m*L. Each selected block provides
L-p safe pre-inverse columns per region. Therefore each region has at least
A=m*(L-p) safe columns, irrespective of the coordinate's block permutation.
Choose exactly that many designated safe columns if additional ones exist.

With a scalar post-scan and uniform regional routing, the preceding exact
window argument now gives

    mu = b*A*(A-1)/R,
    Pr[distance<=rho*a+1] >= max(0,(mu-rho*a)/(b*(A-1)-rho*a)),

when the denominator is positive. Choose the largest m making rho*m*L+1
strictly less than the required distance. At K20 and p=2, the results are:

| Preserved block size L | a | A | Setup-failure lower bound |
|---|---:|---:|---:|
| 8 | 26208 | 19656 | 0; inconclusive |
| 16 | 26208 | 22932 | 0.2995605096 |
| 32 | 26208 | 24570 | 0.4638684066 |
| 64 | 26176 | 25358 | 0.5326341484 |

This probabilistic gate already rejects intact 16-row SIMD blocks.
A stronger grouped-check argument also rejects two-, four-, and eight-row
blocks deterministically.

To obtain that stronger bound, merge the rho checks of each of the L rows
in a common source block. There are g=R/L groups, each with tau=rho*L checks.
If p<=L, every pre-inverse column intersects at most two such groups.
This includes columns crossing a global region boundary. A scalar post-inverse
column therefore touches at most Delta=4 groups after any intervening
permutation, regional or global.

The grouped-check lemma from `HD_BAA_PROOF.md` applies to these merged groups.
For completeness, select m groups uniformly. Each column is confined with
probability at least C(m,Delta)/C(g,Delta). The expected number of confined
columns exceeds tau*m exactly when

    2*C(m-1,Delta-1) > C(g-1,Delta-1).

The resulting dependence gives distance at most tau*m+1. The smallest
integer m satisfying this strict inequality yields the following exact
ceilings for p=2, q=1, and K20:

| Common block size L | Check groups g | Delta | Distance upper bound |
|---|---:|---:|---:|
| 2 | 16384 | 4 | 832321 |
| 4 | 8192 | 4 | 832385 |
| 8 | 4096 | 4 | 832513 |
| 16 | 2048 | 4 | 832513 |
| 32 | 1024 | 4 | 833537 |
| 64 | 512 | 4 | 833537 |
| 128 | 256 | 4 | 835585 |

Every value is below the required 838861. Consequently, permuting intact
common two-row blocks already fails, regardless of the permutation law.
The earlier eight-row probabilistic test was inconclusive only because it
ignored this stronger grouping.

In general, a pre-inverse column touches at most ceil(p/L)+1 merged groups.
To escape the degree-four obstruction while preserving common L-row blocks,
one must at least take p>=L+1 or change the post-inner. For four-row blocks,
this means at least five pre-scans. Different source partitions for different
coordinates also escape this particular argument. Neither escape establishes
a distance lower bound. Arbitrary distinct Hamiltonian row orders remain a
structurally different candidate, subject to both proof and implementation costs.

## Consequence for the next construction

The proof obstruction comes from the common row order before the first
random route. A cluster of p+1 pre-inverse positions always remains a short
interval of the same outer-row cycle. A single later permutation cannot
erase that correlation from the parity-check matrix.

A viable next investigation must change that property or change the
post-scan inverse. Examples include distinct row orders in different outer
coordinate regions, or a post-inner whose inverse is not a bounded-width
difference operator. Neither modification is certified here. Distinct row
orders must be analyzed jointly with the regional routing law; simply
inserting the global-interleaver enumerator would remain invalid.

No support/rank union bound was pursued for the rejected three-scan map.
Its deterministic low-weight dependence makes such a lower-distance proof
impossible. The next proof step is to formulate and test a changed topology,
not to tighten this construction's first moment.

## Reproduction and boundary tests

The new helper is `theory/hd_regional_three_obstruction.py`. Nine bounded
tests in `theory/hd_regional_three_tests.py` pass. They check all binary
safe/unsafe strings through length nine, exact small uniform-window laws,
cyclic containment including wraps, every permutation of a six-position
global route, and weighted matrix supports with global boundary carry.
Additional tests check sampled regional routes, block-preserving placements,
merged-group degrees across global boundaries, and the exact K20 conclusions.
The matrix tests use GF17 because the support statements are field-independent.

```text
python -B research/workstreams/char2_field_spin/theory/hd_regional_three_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_regional_three_obstruction.py
```

Only this new note and the new `hd_regional_three*` proof helpers were written.
No existing encoder, build, benchmark, or remote command was run.
