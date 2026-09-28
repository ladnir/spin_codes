# Count grouped BCH supports through shortened subcodes

The current intermediate-occupancy bound needs better counts of four BCH
words with a small union support. Ordinary weight counts alone do not
describe that union. This refinement combines lower-rank support counts
with bounds on the dimension of shortened codes. It improves smaller
supports, but does not close the support-128 diagnostic or the full proof.

The containment identity below is standard: see Proposition 1 in
[Jurrius and Pellikaan, Extended and Generalized Weight Enumerators (WCC 2009)](https://ruudp.win.tue.nl/paper/53.pdf),
with complementary coordinate sets. Their extension-code interpretation
also identifies ordered binary tuples with words over an extension field.
Our calculation uses these identities as upper bounds, not as an exact
enumerator for BCH[256,128].

## An exact containment bound

Fix a binary linear code C of length n. For a coordinate set S, let C_S
be the subcode consisting of words supported inside S. Write d_t for an
upper bound on dim(C_S), valid for every set S of size t.

Let A_r(v) count r-dimensional subcodes of C whose union support has size
exactly v. Let B_r(u) be their cumulative count through u. The number of
r-dimensional subspaces in a k-dimensional binary space is

    G(k,r) = product_{i=0}^{r-1} (2^k-2^i)/(2^r-2^i).

Set G(k,r)=0 when r>k. Double-counting pairs (H,S), with H an r-dimensional
subcode supported inside a t-set S, gives

    sum_{|S|=t} G(dim(C_S),r)
      = sum_{v=0}^t A_r(v) binomial(n-v,t-v).

The weights on the right decrease with v. Therefore a cumulative upper
bound U_r(v)>=B_r(v) gives the upper

    M_r(t) = sum_{v=0}^t (U_r(v)-U_r(v-1)) binomial(n-v,t-v),
    U_r(-1)=0.

This is summation by parts; the increments need not bound individual shells.

For 1<=r<h<=d_t, the ratio G(k,h)/G(k,r) increases with k>=h: after
cancelling factors, its k-dependent factors are 2^k-2^i for r<=i<h.
For k<h, the numerator is zero. Consequently,

    sum_{|S|=t} G(dim(C_S),h) <= G(d_t,h)/G(d_t,r) * M_r(t).

Every h-dimensional subcode supported on at most u coordinates belongs
to at least binomial(n-u,t-u) different t-sets, for t>=u. Thus

    B_h(u) <= min_{t>=u, 1<=r<h}
        G(d_t,h) M_r(t) / (G(d_t,r) binomial(n-u,t-u)).

If d_t<h for any t>=u, then B_h(u)=0. All quantities are integers or exact
rationals. The implementation floors the final upper for an integer count.

To count ordered four-word groups of rank h, multiply B_h by
product_{i=0}^{h-1}(16-2^i). Conversely, divide an existing rank-h tuple
CDF upper by this multiplicity and floor to obtain U_h. The implementation
processes ranks in increasing order, so a refined lower-rank CDF can help
the next rank. It retains every previous upper and enforces monotonicity.

Exhaustive tests on three small codes verify 88 containment identities.
They also compare tuple CDFs with the bound under exact, inflated, and
deliberately loose input CDFs. No random-code assumption enters this bound.

    python -B research/workstreams/permutation_locality/shortening_moments.py

## Better shortening dimensions without a fragile numerical LP

The older shortening routine rejected many numerical proposals above
length 110. shortening_polynomial.py supplies a separate exact witness
family. Numerical eigenvalues only propose coefficients; acceptance uses
rational arithmetic.

Let K_j be the binary Krawtchouk polynomial of length n. Choose
p=sum_j p_j K_j with p_j>=0, and require that g(x)=(d-x)p(x) also has
nonnegative Krawtchouk coefficients. Then f=p*g=(d-x)p(x)^2 has
nonnegative Krawtchouk coefficients: products K_i K_j have nonnegative
integer coefficients in this basis. This follows by expanding each K_i
as a sum of characters indexed by i-subsets and grouping products by the
size of their symmetric difference.

For integer distances x>=d, f(x)<=0. Its constant coefficient is

    f_0 = sum_j p_j g_j binomial(n,j).

When f_0>0, the standard positive-polynomial counting argument gives
|C|<=f(0)/f_0 for every binary code of minimum distance at least d.
Indeed, summing each character over codeword pairs gives a square, so
the double sum of f is at least f_0 |C|^2. The off-diagonal terms are
nonpositive, giving the upper |C| f(0).

The script checks p_j>=0, g_j>=0, and f_0>0 exactly. It also independently
expands f by orthogonality at every distance for five small-code tests and
the nine reported lengths 80 through 192. The dimension scan accepts 155
positive witnesses through length 192, retains the older bounds, and
propagates shortening by one coordinate. It improves 87 dimension caps.
Examples are 50 to 47 at length 128, 66 to 62 at 144, and 81 to 77 at 160.

    python -B research/workstreams/permutation_locality/shortening_polynomial.py
    python -B research/workstreams/permutation_locality/shortening_moments.py --positive-polynomials

## Effect on the current proof

Later refinement: [EVEN_COLUMNS.md](EVEN_COLUMNS.md) adds dual-complement
shortening. With `--dual-shortening --positive-polynomials`, the count
at support 128 improves to 348.9427 bits, and at support 160 to 396.1708.
The historical comparisons below omit that refinement.

These are log2 upper counts of nonzero four-word groups with support at
most u, summed over all ranks. Displayed decimals summarize exact bounds.

| u | Prior count | Containment moments | With positive shortening polynomials |
|---|---:|---:|---:|
| 72 | 180.3190 | 160.6013 | 160.6013 |
| 80 | 214.2098 | 193.2477 | 193.2477 |
| 96 | 262.2825 | 254.7668 | 254.2645 |
| 112 | 307.3985 | 306.6766 | 303.6766 |
| 128 | 350.9023 | 350.9023 | 350.9023 |
| 160 | 431.2912 | 431.2912 | 431.2912 |

The cover accepts --shortening-moments and --positive-shortening as optional
flags. Its default certificate path is unchanged. With containment moments
but without the positive-polynomial extension, the two-update 64-group
screen gives these best log2 contribution uppers for equal support sizes:

| Support per group | Binary64 diagnostic |
|---|---:|
| 72 | -369.9079 |
| 80 | +748.7216 |
| 96 | +2494.8916 |
| 128 | +3809.0603 |

These include group locations and the all-one-column weighting. They are
selected-point diagnostics, not outward certificates or occupancy coverage.
The support-72 result does not imply a complete 64-group certificate.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --joint-witness --penalties .75 1 --tilts .008 .016 .032 --probe-supports 72 80 96 128

This refinement does not determine the joint weights of all fifteen
nonzero combinations in a rank-four subcode. The difficult intermediate
supports still need stronger counting or inner-state bounds. A concrete
next inner test is to recover full feedback distributions from the existing
character polynomials, rather than only their zero coefficients. Their
largest nonzero atoms would tighten several high-occupancy density bounds.
