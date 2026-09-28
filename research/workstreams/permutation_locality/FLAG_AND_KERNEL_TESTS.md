# Two count refinements that do not close the gap

The dual-moment bound improves the outer support counts, but the difficult
64-group diagnostic remains at +3046.1426 bits. This note records two
subsequent exact counting tests. Neither materially improves that point.
No new occupancy or full-code certificate follows.

## Average extension fibers for all-one columns

Fix an r-dimensional subcode H and write V for its support. The number
of cosets modulo H with a prescribed restriction outside V is at most
2^(dim C_V-r). The old flag bound replaces dim C_V by its worst-case
upper bound for sets of size |V|.

`flag_fibers.py` instead bounds

    F_r(u) = sum_{H: dim H=r, |support(H)|<=u} 2^(dim C_support(H)-r).

For a t-coordinate set T, summing over all r-subcodes of C_T gives

    2^(dim C_T-r) G(dim C_T,r)
      = G(dim C_T,r)+(2^(r+1)-1)G(dim C_T,r+1),

where G(d,r) is the binary Gaussian binomial coefficient. Consequently,
for every t>=u,

    F_r(u) <= (M_r(t)+(2^(r+1)-1)M_(r+1)(t))/binom(n-u,t-u).

Here M_r(t) sums G(dim C_T,r) over all sets T of size t. Each H on the
left belongs to at least the denominator's number of such sets. Its
shortening dimension can only grow when its support is enlarged to T.
These two facts justify the inequality with nonnegative terms.

For a tuple with j>0 all-one columns, its even combinations span H.
If H has support v, there are binom(n-v,j) candidate outside restrictions.
Summation by parts applies because this factor decreases with v.
The resulting bound is intersected with the earlier flag bounds before
maximizing the all-one-weighted count.

The script passes 362 exact small-code inequalities. It tests subcode
fibers, tuple histograms, and inflated upper inputs. However, at penalties
1/2, 3/4, and 9/10, it gives no improvement in the displayed weighted
counts at supports 80, 96, 112, 128, 144, 160, 176, or 192.
For example, the support-128 count at penalty 3/4 stays at log2 347.32235.

The optional `extra_prefix` argument in `occupancy_allones.weighted_cdf`
allows this checked experiment. It is not enabled in the cover by default.

    python -B research/workstreams/permutation_locality/flag_fibers.py

## Weighted kernels for the dual spectrum

Let D=C^perp. It has dimension 128, even weights, minimum distance at least
30, and orthogonal-array strength at least 37. Its zero and all-ones
words are known exactly. For a polynomial p and a weight factor W, define

    L(W p^2) = 2^128 E[W(X)p(X)^2] - W(0)p(0)^2 - W(256)p(256)^2,
    X ~ Bin(256,1/2).

When deg(W p^2)<=37, this is the exact sum over the other codewords of D.
If W is nonnegative on every allowed remaining weight and W(w)p(w)^2>0,
then

    A_D(w) <= L(W p^2)/(W(w)p(w)^2).

`dual_kernel.py` computes rational kernel witnesses for W=1 and
W(x)=(x-a)(256-a-x), with a in {0,16,28,30}. Polynomial degrees are at
most 18 and 17 respectively. All factors are checked on the allowed
weights 30,32,...,226. The Gram matrix contains the exact functional L
on products of Krawtchouk polynomials. An exact inverse proposes p;
the bound evaluates its exact quadratic cost. Singular matrices fall
back to lower degree; a wholly singular case contributes no bound.

Thirty-one small-code checks compare the Gram entries with actual
codeword sums and verify the resulting shell bounds. The production
premise is replayed before the BCH calculation.

The improvement at the relevant low weights is negligible. For example,
the log2 cap at weight 30 changes from 52.44992722827814 to
52.449927228263974; at weight 48 it changes from 64.76081688637483 to
64.76081688636216. These decimal summaries are diagnostics of exact
rational caps. The refinement is not enabled in the default dual-moment
path.

    python -B research/workstreams/permutation_locality/dual_kernel.py

## Next test

The existing BCH sandwich LP contains information beyond these distance-
derived moments. Test whether it gives useful low-weight dual caps,
especially near weights 30--60. The physical objective follows from the
MacWilliams transform of the existing q_w+31h_w spectrum. Any witness must
be checked over the rationals, with empirical orbit-search lower bounds
removed. Numerical solver output alone must not enter the certificate.

These tests change research code only. They do not change the measured
encoder, production defaults, or the paper.
