# Removing all-quiet low-arity groups before the distance bound

The scalar proof loses on sparse quiet words because it charges them the full output-weight mark. A different decomposition pays for certain bad outer setups without that mark. The argument below applies to the unchanged bit-syndrome outer mixer, not the proposed star mixer or a new column-permutation schedule.

Let the field be \(\mathbb F_{256}\). The outer mixer is
\(H_*=I+VBV^T\), where the sixteen rows of \(V\) are the augmented binary labels \((1,i_0,i_1,i_2,i_3)\). The matrix \(B\) is invertible, and \(V\) has rank five. Thus \(H_*x=x\) exactly when \(V^Tx=0\). Call such a column **quiet**, before the outer output diagonal.

There are 8,192 source groups. Each group contains sixteen RS[16,8] rows. In the ideal setup law, every input-diagonal entry is sampled independently from \(\mathbb F_{256}^*\). All entries are sampled once and remain fixed for every encoding. The output diagonals, route, and inner keys do not enter the event defined below.

Fix a nonempty labelled row set \(A\), of size \(a\). Put
\[
 d_A=a-\operatorname{rank}(V_A).
\]
The quiet column space supported on \(A\) has dimension \(d_A\). Define \(E_{g,A}\) as the event that group \(g\) contains a nonzero RS word with active row set exactly \(A\), whose sixteen scaled columns are all quiet.

## Unweighted setup bound

For any prescribed scaled nonzero row \(x\), the expected counting measure of its nonzero RS preimages is at most \(255^{-8}\). Indeed, if \(x\) has weight \(w\), independent nonzero coordinate scalings give a uniform vector on its support. A shortened RS[16,8] code has at most \(255^{w-8}\) words on that exact support. Supports below weight nine contain none.

The measure here is the counting measure on **nonzero** RS rows. Its value at the all-zero prescribed row is zero. Therefore the same upper bound remains valid when an enlarged candidate set includes arrays with some zero rows.

There are \(256^{16d_A}-1\) nonzero arrays whose columns lie in the quiet space supported on \(A\). For each array, row independence gives the upper bound \(255^{-8a}\). A realized nonzero witness contributes at least 255 distinct scalar multiples. Hence
\[
 \Pr[E_{g,A}]
 \le \frac{256^{16d_A}-1}{255^{8a+1}}.
\]
No output-weight mark appears in this inequality.

The exact labelled rank histogram for \(1\le a\le7\) is small:

| Arity | Quiet dimension and number of row sets |
|---|---|
| 1 | 0: 16 |
| 2 | 0: 120 |
| 3 | 0: 560 |
| 4 | 0: 1,680; 1: 140 |
| 5 | 0: 2,688; 1: 1,680 |
| 6 | 1: 7,168; 2: 840 |
| 7 | 2: 11,200; 3: 240 |

The union bound over all 8,192 groups and all these row sets is
\[
 \varepsilon_{\mathrm{quiet}}
 =8192\sum_{1\le |A|\le7}
   \frac{256^{16d_A}-1}{255^{8|A|+1}}.
\]
The isolated exact helper computes this rational number. Independent 256- and 384-bit logarithmic enclosures both give
\[
 \log_2\varepsilon_{\mathrm{quiet}}
 \le -50.771255305346.
\]
The 240 arity-seven sets of quiet dimension three dominate. The arity-four affine-plane contribution alone is
\[
 \frac{8192\cdot140\,(256^{16}-1)}{255^{33}}.
\]

This event is invariant under the existing column-scalar gauge lift: multiplying a column by a common nonzero scalar preserves the quiet subspace. The argument does not treat normalized native diagonal entries as independent.

## How the distance proof may use the event

Outside the global event, every active group with at most seven active rows has at least one nonquiet column. This necessary condition can be retained by a positive quiet/nonquiet flag through the sixteen-column word product.

One must not silently reuse the independent setup law after conditioning on the good event. Instead, split the failure probability into the bad-setup event and its complement. Bound the second term by the **unconditional** first moment restricted to words satisfying the new necessary condition. This avoids conditional renormalization and avoids multiplying \(\varepsilon_{\mathrm{quiet}}\) by \(z^{-734003}\) or by another source-support count.

For a prescribed column-permutation schedule \(P_j\), the same counting argument replaces \(16d_A\) by \(\sum_j d_{P_jA}\). The original RS row set \(A\) remains fixed across columns. Only its coordinate image at the mixer changes. Thus the scheduled event bound is
\[
 8192\sum_A\frac{256^{\sum_{j=1}^{16}d_{P_jA}}-1}{255^{8|A|+1}}.
\]
The sum ranges over whichever active row sets the new event excludes. Stronger deterministic exclusions can use individual RS row support constraints. These are separate construction-bound statements; the saved receipt here covers the unchanged column order and arities at most seven only. In particular, arity-nine quiet dimension is always four. The exceptional five-dimensional zero-outside fibers contain nonquiet vectors and must not replace that quiet dimension.

For the next positive matrix bound, fix \(A\). Write \(q_j\) and \(n_j\) for its quiet and nonquiet raw column kernels at column \(j\), under the common matrix-atom measure. The kernel of words with a nonquiet column is exactly
\[
 F_A=\sum_{j=1}^{16}
       \left(\prod_{i<j}q_i\right)n_j
       \left(\prod_{i>j}(q_i+n_i)\right).
\]
The index specifies the first nonquiet column. For the group exponent \(r\), convexity gives
\[
 F_A^r\le16^{r-1}\sum_{j=1}^{16}
       \left(\prod_{i<j}q_i^r\right)n_j^r
       \left(\prod_{i>j}(q_i+n_i)^r\right).
\]
After integration, retain the chronological product of the quiet, nonquiet, and full two-state column matrices. Do not reset the inner state between factors. The uniform decomposition costs \(4(r-1)\) bits in this per-row-set word bound. Fixed positive weights \(\lambda_j\), summing to one, can replace the uniform factor by \(\lambda_j^{1-r}\) on summand \(j\).

Each branch requires a separately valid positive kernel or upper bound. Subtracting a quiet upper matrix from a full upper matrix is not valid. For arity at most seven, outside zero is equivalent to quiet. A full-fiber comparison can therefore form its nonquiet upper branch by deleting the outside image polynomial's exact constant before any group power. This supplies the required positive split without conditioning the operator on a setup event.

## A cached control that preserves sparse activation

Singletons give a simple full-class test of the same principle. Since every entry of \(H_*\) is nonzero, the exact singleton column kernel is
\[
 K(t_1,\ldots,t_{16};Y)=1+255Y\prod_{i=1}^{16}t_i.
\]
It is independent of the labelled input row. Let \(E(k)\) be the common two-state matrix moment at occupancy \(k\). For \(r\) groups the column coefficient of \(Y^k\) is
\[
 Q_k=\binom rk255^kE(k)^{16}.
\]
The chronological matrix polynomial \(Q(Y)^{16}\) keeps both boundary states. Retaining total degree at least \(9r\) relaxes, but does not remove, the individual RS row constraints. The row-set count \(16^r\), RS factor \(255^{-8r}\), source count, and projective factor are each charged once.

At the saved mark \(9927/10000\), with \(r=18\) and source stratum \(F_{\max}\le1\), the cached FLOAT control gives:

- positive degree tail: \(-4003.457675334371\) bits;
- optimized scalar activity bound: \(-3999.4661718633024\) bits, at \(\log_2Y=67.22765089534391\).

This control uses no uniform information-cover hypothesis and rebuilds no regional provider. Its cached inputs remain FLOAT, including the previously recorded tiny moment-fit rounding deficit. It is not a directed distance certificate.

The files `theory/hd_bit_syndrome_quiet_expurgation.py` and its two tests implement the exact setup-event count. The receipt is `theory/hd_bit_syndrome_quiet_expurgation_le7_exact.json`. The singleton helper and four tests are `theory/hd_bit_syndrome_singleton_cached.py` and its test file; the final cached result is `theory/hd_bit_syndrome_singleton_cached_9927_extended_float.json`.

Next, retain the nonquiet-column constraint in a full sparse-family matrix comparison. The global setup-event bound supplies one probability charge; it does not by itself bound the remaining words or complete whole-code coverage.
