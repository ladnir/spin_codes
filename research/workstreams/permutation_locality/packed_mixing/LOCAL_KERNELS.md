# Local Random Maps in the Packed BCH Representation

This note gives exact local distribution laws and an outer-count interface.
It does not give a SPIN distance certificate or a performance claim.

The placement is **BCH encoding, then the local map, then the existing
independent column and regional permutations, then IMT**. The maps below act
on code coordinates. They do not mix the 128 payload bits within an encoded
element. Setup samples every indicated map once and fixes it for all messages.
The analysis fixes one message before sampling setup.

## Which Maps Change the Relevant Support?

One packed block consists of eight consecutive BCH columns in four rows.
Represent its fixed input by a binary matrix X with four rows and eight columns.
Its packet support counts the nonzero columns of X.

An invertible change of the message basis leaves the code unchanged. An
invertible combination of the four rows also leaves each column's zero status
unchanged. Neither operation addresses the current union-support problem.
The maps considered here instead mix columns, or all 32 coordinates together.

| Local distribution | Required input statistic | Expected active packets |
| --- | --- | --- |
| One shared uniform GL8 map on columns | Binary rank r of X | 8(256-2^(8-r))/255 |
| Independent uniform GL8 maps on rows | Number h of nonzero row fragments | 8(1-(127/255)^h) |
| Independent uniform GL16 maps on the two row-pairs | Which row-pairs are nonzero | About 6.0001 for one pair; 7.5000 for both |
| One uniform GL32 map on all coordinates | Whether X is zero | About 7.500000002 if nonzero |

For shared GL8, ranks one through four give means 4.0157, 6.0235, 7.0275,
and 7.5294. The map preserves linear relations among the rows. Independent
row maps can destroy equality among nonzero rows, but cannot activate an
all-zero row fragment.

Uniform nonzero multiplication in GF(256) has exactly the same fixed-input
law as uniform GL8 on one nonzero row fragment. Thus independent GF(256)
multipliers per row and block can replace independent row-wise GL8 maps in
this first-moment analysis. Sharing one multiplier across rows does not have
the same joint law.

Similarly, uniform nonzero multiplication in GF(2^32) has the same fixed-input
law as uniform GL32. This saves setup entropy, not automatically runtime work.
With a degree-four representation over GF(256), its matrix has sixteen 8-bit
linear blocks. A straightforward SIMD schedule has the same four matrix-term
structure as a general GL32 map; a faster factorization requires measurement.

## Exact Local Kernels

For a shared GLm map, choose a basis for the rank-r row space of the input.
The resulting r-by-m output matrix is uniform among all full-row-rank matrices.
There are product(i=0..r-1)(2^m-2^i) such matrices. If S(r,w) counts ordered
nonzero columns spanning binary r-space, then

    Pr[output support = w]
      = binom(m,w) S(r,w) / product(i=0..r-1)(2^m-2^i).

Subspace inclusion-exclusion gives S(r,w) exactly; `shared_gl` implements it.
The distribution depends on rank, not just the number of nonzero input columns.

For h independently randomized nonzero row fragments, the output rows are
independent uniform nonzero m-bit vectors. Inclusion-exclusion over zero rows
and over selected columns gives the exact union-support kernel implemented
by `independent_rows`. `independent_pieces` gives the analogous GL16 row-pair
law. These statements remain valid for the corresponding independent field
scalar maps, because each field action is transitive on nonzero vectors.

For a nonzero four-row block under full GL32, every nonzero 32-bit output is
equally likely. Therefore its eight-packet support W has distribution

    Pr[W=w] = binom(8,w) 15^w / (2^32-1),  1 <= w <= 8.

This is Binomial(8,15/16) conditioned on being nonzero. Unlike the smaller
maps, its law does not depend on the input rank or row-activity pattern.

### Full GL32 Already Supplies the Packet Labels

For each fixed packet support S, a uniform nonzero 32-bit block has exactly
15^|S| equally likely labelings. Its active packet labels are consequently
independent uniform nonzero four-bit vectors, conditional on S. Independent
maps on different blocks give the product law conditional on the full support.
The later independent permutations preserve that law.

Thus a **final** independent full GL32 or GF(2^32) scalar map per block can
replace, rather than precede, the existing packet-wise GF16 randomizers. This
corrects the earlier informal suggestion that those randomizers were needed
for this particular distribution. The smaller shared or row-wise maps still
need the GF16 stage to reuse the current inner interface.

The maps must have the stated distribution. An arbitrary invertible 32-bit
matrix, or a short random shear network, is not automatically uniform GL32.

## Canonical Blocks Versus a Random Partition

Let u be the original union support and H the number of nonempty canonical
eight-column blocks. Knowing u alone gives only H >= ceil(u/8). For example,
112 active columns could fill just fourteen blocks. Full GL32 then has mean
output support 105, so this coarse argument does not establish expansion.

A fresh uniform partition before mixing gives an exact distribution for H:

    Pr[H=h | u]
      = binom(32,h) [z^u] ((1+z)^8-1)^h / binom(256,u).

That distribution is implemented and checked against exhaustive small cases.
It cannot be inferred from the existing permutation, which occurs after the
local map. A random input partition has a separate implementation cost.

### A Stronger Count for the Canonical Partition

The coarse H bound discards the restriction to canonical blocks. Let d(v)
upper-bound the dimension of the BCH code shortened to every fixed set of v
coordinates. A tuple occupying at most h canonical blocks lies in a shortened
code on the union of some h blocks. Consequently, the number of such nonzero
four-row tuples is at most

    binom(32,h) (2^(4 d(8h))-1).

For tuples of binary rank r, replace the parenthesized term by the exact
number of four-by-d(8h) binary matrices of rank r. Clamp these bounds by the
total rank counts. If available, also clamp by the authenticated original
union-support CDF at 8h. Reverse monotonicity repair yields cumulative caps
B(h). `canonical_rank_caps` implements these steps.

The choice factor is binom(32,h), not binom(256,8h). This counts only supports
that are unions of canonical blocks. It can therefore be much stronger than
applying a union-support cap at 8h, without adding a random input partition.

## Transporting Canonical Counts to the Existing Inner Interface

For h nonempty blocks, let F_h(v)=Pr[W_1+...+W_h <= v], using the exact full
GL32 kernel and independent maps. Because every W_i is positive, F_h(v)
decreases with h. Summation by parts gives the expected nonzero-message CDF
bound

    D(v) = sum(h=1..32) (B(h)-B(h-1)) F_h(v).

CDF differences are valid here only inside this decreasing-kernel transform.
They are not bounds on the actual canonical shell counts. Likewise, the
output differences D(v)-D(v-1) are not output-shell bounds.

For the positive-mixture dense interface, D(v) itself bounds the expected
number of words in shell v. After the independent support shuffle, every
particular vector of packet support v has the same probability, including
the product nonzero labels proved above. A checked positive Bernoulli
mixture dominating these shell caps therefore dominates the averaged group
measure pointwise.

This average can be multiplied across groups only because the local maps
and setup randomness are independent across groups. The subsequent regional
routing and inner failure calculation are nonnegative functionals of this
comparison measure. There is no claim that each sampled group has spectrum
D, and no application of a nonlinear spectrum bound before averaging. The
true zero message remains separate from artificial zero draws of positive
comparison components, as in the existing proof interface.

One globally shared random partition or table of sampled local maps would
require a different argument. Raising a mean enumerator to the number of
groups would not justify that construction.

### Tighter Pointwise Shell Bounds

The dense comparison need not use the whole CDF D(v) as its shell-v bound.
Let K_h(v)=Pr[W_1+...+W_h=v] and set M_h(v)=max(j>=h) K_j(v). The suffix
maximum is decreasing in h. The same cumulative-count argument therefore
gives a pointwise expected-shell bound

    A(v) = sum(h=1..32) (B(h)-B(h-1)) M_h(v).

This does not interpret the differences of B as actual shell counts. For any
nonnegative counts with cumulative bounds B, replace K_h by its decreasing
suffix envelope and apply summation by parts. The resulting upper bound is
sharp if only these cumulative constraints are known. Coupling additional
positive W_i shows M_h(v)<=F_h(v), hence A(v)<=D(v) at every support v.
The sum of A(v) need not equal the number of messages.

`monotone.transport_shells` computes these rational bounds. Exhaustive toy
tests compare the transform with every feasible spectrum and check its
linear-program optimum. On the current refined BCH bounds, exact pruning
reduces the positive mixture from 46 components to 42. Its log2 tilted mass
at cost tilt 1/4 decreases from 64.241961 to 63.864705. This is an outer
comparison improvement, not a whole-code margin.

The distinction between pointwise bounds and Laplace-moment bounds matters
here. A(v) dominates each expected shell, so subsequent support/comparison
cells may use different output tilts. A comparison that dominates only a
complete Laplace moment would not justify choosing a different tilt in
each artificial comparison component. The current driver uses the pointwise
route; it does not use Fourier thinning.

### Narrower Full Maps in the Same Packed Layout

The packed 32-coordinate map can be block diagonal: two full maps on four
rows by four columns, four maps on four rows by two columns, or eight maps
on single four-row columns. These are GL16, GL8, and GL4 maps, respectively.
They differ from the earlier GL16 maps on two rows by eight columns.
Each full subblock map retains the conditional product-label property.

Width one preserves packet support exactly. For each fixed input packet,
uniform GL4 has the same law as the old uniform nonzero GF16 multiplier.
Fusing the actual old scalar matrices into the packed BCH operation could
therefore be an implementation-only optimization, with no construction change.

`width_compare.py` regenerated one refined dimension table and reused it for
all four widths. The resulting log2 expected CDF bounds were:

| Packet-support cutoff | Width 1 | Width 2 | Width 4 | Width 8 |
| --- | ---: | ---: | ---: | ---: |
| 32 | zero | 72.8216 | 32.6462 | 13.1183 |
| 64 | 211.2812 | 141.8291 | 116.2760 | 89.7478 |
| 112 | 384.7912 | 336.3230 | 282.2866 | 245.7579 |
| 115 | 393.7722 | 344.1729 | 288.9834 | 252.3099 |
| 128 | 439.6728 | 370.0879 | 311.7615 | 276.5963 |
| 192 | 489.5529 | 442.6277 | 442.0047 | 442.0046 |

This comparison uses canonical shortening bounds alone. In particular, the
width-one column is not the much stronger existing BCH spectrum certificate.
Additional authenticated spectrum bounds can be intersected with these caps.

The current proof method favors width eight at medium supports: width four
loses about 35--37 bits per group at cutoffs 115 and 128. Narrower maps avoid
some new low-support outcomes: the count bounds first become nonzero at
supports 38, 19, 10, and 5 for widths one, two, four, and eight. For width
eight the expected CDF through support five is below 2^-103.93, so that new
tail is small but cannot be omitted from a certificate. The full inner
calculation, rather than these outer bounds alone, must decide the tradeoff.

## Cheap Cross-Block Diffusion

Fresh scalar shears over GF(256) give another exact local interface. For two
row fragments x,y, sample independent nonzero a,b and set

    A = x + a y,   B = y + b A.

This transform is invertible for every a,b. With d=255, input activity 10
always becomes 11; activity 01 becomes 10 with probability 1/d and otherwise
11. Activity 11 becomes 01 with probability 1/d, 10 with probability
(d-1)/d^2, and otherwise 11. `two_shear_masks` tensors these laws over rows
using fresh row-wise scalars.

Thus a nonempty pair remains nonempty, and two nonempty blocks collapse to
one with probability at most (2d-1)/d^2. This does not by itself uniformize
the output labels. A final rank-free block map would supply that property.
Multiple pairings form a different construction, not a uniform global map.

## Verification and Remaining Work

The standard-library tests enumerate all GL3 matrices, all GL4 matrices, small
row-fragment outputs, field scalar orbits, random partitions, small shortened
codes, and exact local convolution results. They also test input rejection
and conditional packet-label independence.

`authenticated_bch_cdf` checks the retained BCH premises and constructs fresh
shortened-dimension bounds before returning the expected CDF. Its printed
logarithms are diagnostic displays of exact rational quantities. The dense
driver now builds and checks a positive mixture from the direct expected-shell
bounds above. A complete occupancy and support cover, followed by outward
replay, is still required for a certificate.

The first fresh refined run authenticated 163 BCH dependency files, accepted
67 exactly checked primal LP witnesses and 155 positive-polynomial witnesses,
and replayed the dual-distance and containment checks. The resulting selected
log2 CDF caps are:

| Canonical active blocks H | Original-tuple CDF bound | Output packet support v | Expected output CDF bound |
| --- | ---: | --- | ---: |
| 5 | 21.5264 | 32 | 13.1183 |
| 8 | 31.3208 | 64 | 89.7478 |
| 12 | 107.7504 | 112 | 245.7579 |
| 14 | 164.8125 | 115 | 252.3099 |
| 16 | 217.1630 | 128 | 276.5963 |
| 24 | 311.3264 | 192 | 442.0046 |

These are counts, not failure probabilities. No distance margin follows from
this table alone. The side-by-side columns are separate selected cutoffs, not
a deterministic mapping from H to v.

**Sparse-verifier interface:** the retained `occupancy_cdf_cover.py` assumes
packet support in [38,256], including a hardcoded probe check and cover-volume
calculation. Full GL32 can reduce output support below 38; the current local
count bound only forces at least five active canonical blocks. That sparse
verifier's 38-based cover cannot be reused unchanged. The new `sparse.py` and
`support_cover.py` authenticate the first positive expected-count entry as
five, keep fractional counts, and verify the entire support-box split tree.
The q=1 fold includes all 257 support entries and clips each conditional
failure bound at one. For q>=2, each support box uses a rational Bernoulli
comparison with directed CDF folding. Independent per-group maps permit
products of these expected-count bounds. Different actual support boxes
may use different output tilts.

The sparse replay reconstructs the BCH bounds and inner operators; saved
floating optimizer results and aggregate bounds are not proof inputs.
The dense replay similarly reconstructs the comparison, verifies each shell
inequality, and recomputes every accepted interval bound. A claim for all
messages additionally requires complete, disjoint coverage of q=1..2048 and
a checked sum of the sparse and dense first moments. Passing selected
occupancies or using a per-cell target is not sufficient.

### Sparse Adapter

`sparse.py` and `support_cover.py` implement the new interface without changing
the earlier BCH-only verifier. They regenerate the refined canonical CDF and
require its first positive entry to be five. Counts remain exact rational
expectations; rounding small counts up to integers would discard useful tail
information. The one-group calculation averages placements separately for all
257 support sizes. The other sparse occupancies cover the full ordered domain
`[5,256]^q`, using symmetric boxes with their exact multiplicities.

The product of expected count measures is valid because mixer setup is
independent across groups. Conditional on the support of each group, the
independent column and region shuffles give the existing fixed-support inner
interface. The proof may choose different output tilts for different support
boxes: these boxes partition actual support vectors, not comparison-mixture
components.

For a support box, the verifier introduces one Bernoulli activity parameter
per group. Nonnegativity bounds each fixed-support moment by the Bernoulli
moment divided by its support probability. A suffix maximum turns this
reciprocal-binomial weight into a nonincreasing function of support. Abel
summation then bounds its integral using the rational expected CDF. Applying
this bound to each independent group gives the box contribution. The factor
`binom(2048,q)` selects the active groups.

Every saved box witness contains exact probability numerators. Its split tree
must reconstruct the complete domain; matching total volume alone is not
accepted. Fresh replay authenticates the count premises, rebuilds inner
operators at the requested precision, and reevaluates each witness. Saved
numeric bounds are never accepted as proof inputs. A sparse receipt covers
only its listed occupancies and is not a whole-code certificate.
