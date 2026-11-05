# A second sum removes the two-row cancellation

The current XOR2 candidate is fast, but two equal active inputs survive its
column mixer with weight two. This local behavior forces the proof to track
which constituent pairs occur across groups. A second sum can remove that
case without returning to a full eight-by-eight matrix product.

This is a new research construction. The old XOR2 and certified Toeplitz
checkpoints remain unchanged. The results below establish exact local algebra,
not whole-code distance. The separate [native comparison](HD_TWO_SUM_NATIVE.md)
measures both default-off variants against unchanged controls.

## Column map and cost

Keep eight RS[32,16] constituents per outer group. For each column, use
independent nonzero input and output diagonals, as in XOR2. Replace only the
fixed matrix between them. Over the AES representation of GF(256), define

    u = (0,1,2,3,4,5,6,7),
    H = I + 11^T + uu^T,
    S = sum_i x_i,   T = sum_i u_i x_i,
    (Hx)_i = x_i + S + u_i T.

The labels here are field bytes, not ordinary real numbers. Their XOR sum
is zero. Thus sum(u_i^2)=0 as well. If U has columns 1 and u, then U^T U=0,
so H^2=(I+UU^T)^2=I. The map is invertible. The same common-scalar gauge
normalizes one D1 entry, without changing the distribution of D2 H D1.

The weighted sum is cheap because u_i uses only three binary basis bits:

    A = x1+x3+x5+x7,
    B = x2+x3+x6+x7,
    C = x4+x5+x6+x7,
    T = A + 2*B + 4*C.

Two fixed field products compute T. Two more compute 2*T and 4*T; XORs
then produce every u_i*T. Shared pair sums also supply S. The explicit
circuit in `theory/hd_two_sum_outer.py` uses four fixed field products and
28 binary XORs per eight outputs. Its normalized variable-product count
remains 15, with no additional sampled coefficients. XOR2 uses the same
15 variable products, no fixed products, and 15 binary XORs.

These are circuit counts, not instruction-throughput or timing predictions.
GFNI can evaluate each fixed byte-linear product lane-wise. The extra work
must be measured in the full encoder before choosing this alternative.

## An affine coefficient set removes the remaining sparse types

A second candidate uses u=(8,9,10,11,12,13,14,15), leaving the definition
H=I+11^T+uu^T unchanged. The same sum identities imply H^2=I. All entries
of this H are nonzero. Its two-column restrictions again have all
two-row minors nonzero: in the mixed-minor argument below, both factors
u_i and u_i+u_j+u_k now have polynomial degree three. Their product has
degree six, with no field reduction, and cannot equal one.

Every singleton therefore has full output support. Every fixed input pair
has exactly the same support distribution: full support with probability
247/255, or any specified one-coordinate hole with probability 1/255.
This removes the distinction involving coordinate zero, not just the old
XOR2 pair-multiplicity statistic. The statement uses independent nonzero
input scalings; it holds for every fixed nonzero pair of input values.

Compute T_low using labels 0 through 7, then set T=T_low+8*S. Use
S'=S+8*T as the shared output offset, followed by i*T for i=0,...,7.
The circuit costs six fixed field products and 30 binary XORs per column.
It has no additional random keys. The affine candidate has its own exact
helper and tests; the earlier kernel artifact describes only labels 0..7.

The uniform sparse law also gives an exact positive column calculation.
Suppose s groups have one active coordinate and d groups have two.
Write b=s+d. Each of the d latter groups independently contributes no
hole with probability 247/255, or one hole in each of eight regions with
probability 1/255. If R_i(m) is a nonnegative regional matrix bound valid
for every m-element active subset, the averaged ordered product is bounded by

    d!/255^d * [v^d] exp(247*v)
       * product_(i=0,...,7) (sum_(h=0,...,d) R_i(b-h)*v^h/h!).

Terms with negative occupancy are zero. The product retains chronological
region order. This follows by summing the exact multinomial distribution
of labelled holes; it does not assume regional active subsets are independent.
Apply the common state potential after this positive mixture. A binomial
recurrence computes the same coefficient without factorial denominators.

This is a proposed sparse proof interface, not an inherited certificate.
The two affine variants must remain separately labelled in timing and
proof artifacts. Dense groups still require analysis.

## Exact sparse behavior

Every one- or two-coordinate nonzero input produces at least seven nonzero
output coordinates. This is stronger than a statement about typical inputs:
there is no two-input equal-value exception.

For two selected input coordinates i,j, every two-row minor of their
eight-by-two image matrix is nonzero. Two rows outside {i,j} give determinant
(u_i+u_j)(u_k+u_l). The two rows inside {i,j} give determinant one. A mixed
minor is 1+u_i(u_i+u_j+u_k), after exchanging i,j if needed.

The latter cannot vanish. All three labels have polynomial degree at most
two, so the product has degree at most four and undergoes no AES reduction.
Two nonzero such polynomials multiply to one only when both are one.
That would require u_i=1 and u_j=u_k, contrary to distinct row indices.
Thus two output zeros force both input values to be zero.

Under independent nonzero D1 values, the exact support distribution for two
active inputs has only two types:

| Fixed input pair | Permitted missing output | Probability of each permitted single missing output | Probability of all eight outputs |
|---|---|---:|---:|
| Contains coordinate 0 | Any of coordinates 1 through 7 | 1/255 | 248/255 |
| Does not contain coordinate 0 | Any of the eight coordinates | 1/255 | 247/255 |

A single input at coordinate 0 misses output 0; every other singleton has
full output support. Conditional on support, the independent output diagonal
again supplies independent nonzero amplitudes. Fixed input row sets remain
fixed across all 32 columns; this argument does not resample them.

Consequently, a sparse-cohort proof can distinguish one versus seven
singleton choices and seven versus 21 pair choices. It need not retain the
28 different constituent-pair labels used by the old XOR2 cancellation
analysis. This reduces proof bookkeeping, but is not an automatic transfer
of any old distance certificate. Regional order and the inner state remain
part of the proof.

## What the extra sum does not solve

Three active inputs can still annihilate both S and T. For distinct labels
i,j,k, set their values to (u_j+u_k,u_i+u_k,u_i+u_j). Every value is nonzero,
both sums vanish, and H fixes this weight-three vector. Therefore the local
branch number min(wt(x)+wt(Hx)) is exactly six, not nine.

For any fixed three-coordinate support, there are 255 such nonzero vectors.
Their probability under independent nonzero input scalings is 1/255^2.
Larger active sets have analogous two-check kernels. They must be counted
as codeword families; they cannot simply be discarded as rare bad setup.

The exact local minimum output weights are:

| Nonzero input coordinates | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Minimum nonzero outputs | 7 | 7 | 3 | 4 | 3 | 3 | 1 | 1 |

Dense cohorts therefore still require a proof. Higher local weight alone
does not establish monotonic improvement of the recursive encoder's bound.

## Exact checks and next experiment

The helper computes every input-support/output-support count. For a uniform
input on coordinate space A, outputs contained in T number
256^(|A|-rank(H_[outside T,A])). Integer subset inversion first enforces exact
output support, then exact nonzero input support. This avoids a numerical
probability approximation. The full kernel has 256 by 256 entries.

Tests independently enumerate every one- and two-input projective vector,
check all 64 binary circuit basis vectors, and verify H^2=I. They also check
full-kernel reversal symmetry, total masses, the branch-six witness, and the
two-type pair law. The exact artifact is `hd_two_sum_outer_exact.json` in
`theory`; it records local counts and explicitly disclaims a whole-code
certificate or benchmark.

Both circuits are now implemented as default-off experiments. Serial matched
full-encoder medians are 3.537 ms for labels 0..7 and 3.574 ms for labels
8..15. Their respective certified controls take 3.663 and 3.668 ms. The
affine circuit costs about 3% over its matched 3.471-ms XOR2 control.
Basis tests and independent full K20 oracles pass; old builds remain
byte-identical. The native note supplies the complete protocol and hashes.

For the affine ensemble, a first FLOAT screen covers all 18 sparse cohort
splits at r=17 and source geometry (8,0). Two frozen witnesses, without
fitting, give a union of -338.941883004469 bits; the first witness wins every
cohort. This is not a directed result, other source geometries were not part
of that initial screen, and dense groups remain outside its scope.

A subsequent inherited-witness-only screen covers all 81 ordered source
geometries and all 18 sparse cohorts at r=17. Its 1,458 cells give a FLOAT
log2 union of -246.921306449020. It takes 20.03 seconds, reusing one raw
regional bank and 81 source-conditioned cap banks. Exact source counts sum
to binomial(8192,17). Every cell passes the 40-bit screen; no fitting or
pair-multiplicity partition is used. The subsequent
[directed replay](HD_TWO_SUM_AFFINE_SPARSE_CERTIFICATE.md) now confirms a
log2 union ceiling of -246.921306449012. All 1,458 cell ceilings agree at
256/384 bits; independent 512- and 1024-bit rounded-cell sums agree.

The isolated wrapper is `theory/hd_two_sum_affine_sparse_all_sources.py`;
its artifact is `hd_two_sum_affine_sparse_r17_all81_inherited_float.json`.
The earlier one-source artifact remains unchanged. Tests check the exact
hole recurrence against direct labelled assignments, noncommuting order,
column coefficients, source coverage, and the bounded execution interface.

Next, assess dense groups and other active-group counts. The performance cost is
small enough to justify this proof-guided alternative. Preserve the existing
XOR2 proof work; the entire code has not yet been certified.

The route-independent terminal argument in
[the coverage plan](HD_XOR_COVERAGE_PLAN.md) also applies to this eight-row
variant. It needs only an invertible column mixer, independent unnormalized
D1 scalars, and the unchanged RS[32,16] membership bound. Those conditions
hold for both two-sum matrices. Thus the same bound covers total active
constituent count at least 61,191, with log2 ceiling -165.775266208552.
This limited transfer does not cover the intervening dense cohorts.

The [sixteen-row proposal](HD_TWO_SUM_GF16_PROPOSAL.md) is a further
construction-level opportunity. An explicit affine GF(16) coset preserves
the uniform singleton/pair law while allowing sixteen RS[16,8] constituents.
Its circuit saves 168 field products per equal-size outer batch. The first
implementation regressed, but eliminating an outlined store helper removes
64 MiB of intermediate stack traffic per K20 encode. A matched three-round
comparison then gives 3.398 ms versus 3.616 ms for the eight-row affine
candidate. The [native record](HD_TWO_SUM_GF16_NATIVE.md) preserves both
results. Its fixed-r17 sparse FLOAT screen clears 511 bits, but its
deterministic low-weight outer family makes a dense inner-aware proof
necessary before choosing that design.

## Closed kernel for the affine variants

The affine label sets permit a stronger simplification than the pair law.
The entire uniform-input support kernel has a closed form. Consequently,
summing fixed row sets of one arity produces a weight depending only on
the output weight, not its labelled positions. This removes an exponential
mask enumeration from the dense proof interface. It does not remove the
chronological regional matrices or certify the complete code.

This argument covers labels 8..15 and the sixteen-label GF16 coset, not
the original labels 0..7. Write their label set as `C=gamma+V`, where V is
a binary vector space, `0` is not in C, and C is disjoint from its inverse.
Fix an input coordinate set A of size a. Let `U=sum_(i in A) u_i` and
`epsilon=a mod 2`. The two-by-two matrix

\[
 D_A=\begin{pmatrix}1+\epsilon&U\\U&1+U^2\end{pmatrix}
\]

has determinant 1 when a is even and U^2 when a is odd. In the odd case,
U lies in C and is nonzero. Thus every principal submatrix H[A,A] is
invertible. Uniform input on the coordinate space F^A therefore gives
independent uniform inside outputs `y_A`, including zeros.

Put `sigma=sum x_i` and `tau=sum u_i*x_i`. The inside outputs determine
these moments through

\[
 D_A\binom{\sigma}{\tau}
 =\binom{\sum_{i\in A}y_i}{\sum_{i\in A}u_i y_i}.
\]

An outside output j is `sigma+u_j*tau`. As a linear form in `y_A`, all
its coefficients are nonzero. For even a, a zero coefficient would require
`u_i=U+1/(U+u_j)`; for odd a it would require
`u_i=U+u_j+1/U`. In each case, closure of V and inverse-disjointness of C
exclude membership in C. Two distinct outside equations force both
moments to zero.

Now fix the nonzero inside-output support B, with b coordinates. Define

\[
 z_b=\frac{255^b+255(-1)^b}{256},\qquad
 m_b=\sum_{j=0}^b(-1)^{b-j}\binom bj256^{\max(j-2,0)}.
\]

There are z_b nonzero assignments satisfying one outside equation. There
are m_b satisfying both moment equations. The second formula follows by
inclusion-exclusion: the two Vandermonde checks have rank `min(2,j)` on
every j-coordinate subset. These counts include `z_0=m_0=1`.

Let O be the complement of A and let c be its size. For this fixed B,
the exact output counts are

| Outside support | Assignment count |
|---|---:|
| Empty | `m_b` |
| O except one specified coordinate | `z_b-m_b` |
| All of O | `255^b-m_b-c*(z_b-m_b)` |

All other outside supports have count zero. A nonzero moment vector has
at most one zero outside output, which proves the partition. When c is
zero or one, coincident support rows in the table must be added before
using the count or raising it to a power.

For any fixed target support T of weight w, sum over
`b=|A intersect T|`. There are `C(w,b)*C(g-w,a-b)` such sets A. Let
`N(a,b,w-b)` be the combined table count for that inside/outside weight.
Then the exact arity-class moment used by the whole-word Holder bound is

\[
 W_a(w)^n=\frac{1}{255^{(n-k)a}}
 \sum_b\binom wb\binom{g-w}{a-b}N(a,b,w-b)^n.
\]

This sums the fixed row sets once; it does not resample A by column.
The weights are comparison measures, not probabilities. For Group8 use
`(g,n,k)=(8,32,16)`; for Group16 use `(16,16,8)`.

The helper `theory/hd_two_sum_closed_kernel.py` implements both the kernel
and these exact integer moments. Six tests compare all 65,536 Group8
kernel entries against independent rank/Mobius enumeration, and all Group8
arity moments at powers 1, 2, and 32. Group16 checks include arbitrary
rectangular containment ranks, exhaustive input-set sums for five target
supports, and all principal determinants and outside coefficients.
An independent algebraic review agrees with the derivation.
The frozen helper has SHA256
`2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979`.

The next proof implementation can track output weight through ordered
regional factors, then apply W_a(w). It need not retain all 2^g labelled
dense masks. Two dense groups require two weight indices; their row-set
assignments remain fixed under the same Holder argument.
