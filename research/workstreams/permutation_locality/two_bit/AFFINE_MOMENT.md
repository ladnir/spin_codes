# Keeping the inner moment coupled to the outer count

The lifted cover still maximizes an inner bound and an outer count
separately on each cell. Both bounds can be valid yet lose their
correlation. This refinement bounds the logarithm of the inner moment
by an affine function of the same mixture coordinates. The outer sum
then retains that affine function instead of replacing it by a constant.
The SPIN ensemble and comparison components remain unchanged.

Use the definitions from [LIFTED_MIXTURE.md](LIFTED_MIXTURE.md).
In particular, G=4096 pairs occupy each of L=256 regions. A composition
with integer counts n_i has mean m=sum_i(n_i/G)f_i in a rectangular cell C.
It has at least q_min active labels. Set rho=q_min/G.

## The exact moment is a positive polynomial

Fix the output tilt lambda>0. For three positive, unnormalized category
weights w=(w0,w1,w2), define

    H_lambda(w) = sum_x w0^(a0(x)) (w1/2)^(a1(x)) w2^(a2(x))
                       * E_T[exp(-lambda*weight(I_T(x)))].

Here x ranges over all strings of G*L two-bit packets. Category weights
are w0 for 00, w1/2 for each single-bit orientation, and w2 for 11.
The expectation averages the two independent transvections at each step.
The expansion and feedback maps are fixed, and the initial state is zero.

Thus H_lambda is a homogeneous polynomial of degree G*L with nonnegative
coefficients. Its logarithm is convex as a function of log w: after the
substitution w_j=exp(u_j), it is a log-sum-exp of affine functions of u.
This statement concerns the exact moment. It does not require the
two-state envelope, or its minimum of several bounds, to be log-convex.

The existing envelope supplies an upper bound at any positive w:

    H_lambda(w) <= S^(G*L) e_zero K_lambda(w/S)^16384 (1,1)^T,
    S = w0+w1+w2.

The exponent G*L is the number of packets, not the number of bits.

## Turn affine weight bounds into a convex function

Fix an alternative input tilt z'. Let Z'_i=pi_i dot z' and
D'_i=d_i*(Z'_i)^L. The exact unnormalized reference weights for a
composition are averages of v_ij=pi_ij/Z'_i.
For a rational vector e_j and mu_j>=0, define

    h_j = max_i (v_ij - e_j dot f_i + mu_j*b_i),
    g_j(m) = h_j + e_j dot m - mu_j*rho.

Every permitted composition satisfies w_j<=g_j(m), by the same
componentwise argument as the lifted cover. Choose a rational anchor
A_j>0. The tangent inequality log u<=log A_j+u/A_j-1 gives

    w_j <= g_j(m) <= W_j(m),
    W_j(m) = A_j exp((g_j(m)-A_j)/A_j).

The second inequality also holds at g_j=0. At feasible means,
g_j is nonnegative because it bounds a nonnegative weight. No positivity
of g_j is required at infeasible vertices of C; W_j is positive everywhere.

Each log W_j is affine in m. Consequently

    phi(m) = log H_lambda(W(m))

is convex on C. At every vertex v, compute an outward upper B_v for
phi(v), using the existing two-state envelope. Before that evaluation,
round each W_j(v) upward to a rational dyadic weight. Monotonicity of
the positive polynomial justifies this rounding.

For any rational slope vector a, put

    A = max_{vertices v} (B_v - a dot v).

The affine function A+a dot m dominates phi at all vertices. Convexity
therefore makes it an upper bound throughout C. This remains valid when
some cell coordinates are fixed and vertices coincide. The output tilt
lambda must be the same at every vertex.

## Sum the correlated bound

Let R_C be the earlier regional density loss. For any outer dual eta and
mu>=0, define

    P = sum_i D'_i exp((a/G+eta) dot f_i + mu*b_i),
    ell = min_{m in C} eta dot m.

The complete-cell first moment is at most

    exp(lambda*d + A) R_C^L
        * P^G exp(-G*ell - mu*q_min).

Here d is the inclusive output-weight threshold; the original run used
d=209715. Lower thresholds weaken the claimed distance and improve this
bound without changing the code distribution.

Indeed, exp(a dot m) factors across component assignments as
product_i exp((a/G) dot f_i)^(n_i). The remaining exponential dual
bounds the restriction to C and the active-label threshold. The
multinomial theorem then sums every labeled composition in the cell.

All slopes, anchors, tilts, and duals are rational witnesses. Replay
recomputes the componentwise inequalities, every vertex bound, and the
outer sum with Arb. Floating proposals and stored scores are not evidence.
The exact split-tree and integer-contraction checks are inherited from
the lifted cover.

## Current scope and reproduction

At q_min=512 and theta=2/5, the previous method gave log2 upper
+333.69829 on the selected rectangle

    m1 in [3/160,19/1000], m2 in [9/10000,11/10000],
    central fraction in [9/10000,11/10000].

The coupled bound outward-verifies that entire rectangle at log2 upper
below -293.04967. A larger nearby rectangle still fails. Neither result
is a complete dense-range certificate.

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --probe-cell 3/160 19/1000 9/10000 11/10000 9/10000 11/10000
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --resume tmp/two-bit-lifted-dense-512-continued.json --max-cells 3000 --max-depth 48 --output tmp/two-bit-affine-dense-512.json
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --replay tmp/two-bit-affine-dense-512.json
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

The resume command accepts a lifted cover and independently rechecks its
accepted leaves before using the new bound on unresolved cells. Replay
refuses incomplete covers. Full-code closure still needs the complementary
occupancies and the final aggregate. The tests include exhaustive small
composition sums and positive-polynomial checks of the convex vertex bound.

The full-cover pilot exposed an overflow in the floating tangent proposal,
before outward verification. Anchors now also account for the maximum
affine weight at the vertices, and vertex normalization uses log-sum-exp.
Unusable numerical proposals are discarded; the cover then retains its
previous bound or subdivides. A regression test covers a nonpositive
cell-center weight. These changes do not alter the exact inequality.

The optional refinement now retains the
[profile-dependent shuffle loss](DENSITY_TANGENT.md) inside the input
moment instead of maximizing it separately. Lowering the output-weight
threshold closed the full range at 5%, 6%, and 7% using the simpler lifted
bound. The affine refinement subsequently closed 8%, 8.5%, and 9% without changing the
encoder. Independent 256-bit replay confirms the dense certificate and
the aggregate with the prior sparse bounds. See
[COVER_STATUS.md](COVER_STATUS.md).

The optional `--polish-tilt` search optimizes the output tilt against the
final coupled bound. The default search proposes that tilt from an earlier,
uncoupled bound and tests three nearby values. Polishing changes only the
witness search: replay still checks the same rational witness and outward
inequalities. It improved two selected 9.5% cell bounds by about 1.4 and
1.7 bits, respectively; it has not supplied a complete 9.5% certificate.
