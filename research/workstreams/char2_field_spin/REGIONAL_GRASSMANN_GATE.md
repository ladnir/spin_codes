# Fresh orientations in the sparse rank bound

This note records an exact local improvement for the separate construction that refreshes its state orientation more often. It does not modify the once-per-region construction. Refreshing every 64 symbols retains one regional permutation stage and adds an estimated 4.42% to its scalar-product volume. The root-only improvement is useful, but does not itself close the sparse class.

The later [Cauchy residue-coordinate construction](CAUCHY_DIAGONAL_PACKET_COST.md) uses a cheaper diagonal refresh. Its probability law is different. The exact Grassmann law below applies to a fresh uniform general linear matrix, not to that diagonal matrix.

## Conditional root law

Let the coefficient field be \(\mathbb F_q\), and let the state have dimension \(s\). Fix the outer support, all previously revealed randomness, and the common span of previously imposed zero equations. Let \(S\) be the state-form map on the resulting quotient, with rank \(r\). A positive-state branch requires the witness's incoming state to be nonzero.

Sample \(B\gets\mathrm{GL}(s,q)\), independently of this history and all unrevealed input diagonals. The image \(U=\operatorname{im}(BS)\) is a uniform \(r\)-dimensional subspace of \(\mathbb F_q^s\). State-coordinate changes do not change the row span of the state forms. In particular, they do not randomize whether a source form belongs to that row span.

Suppose \(z\le s\) zero-input, zero-output positions impose independent emission equations. Their common kernel \(K\) has dimension \(s-z\). Put \(u=\dim(U\cap K)\). For every integer \(u\),

\[
 p_{r,z}(u)=
 \frac{{s-z\brack u}_q{z\brack r-u}_q
              q^{(s-z-u)(r-u)}}{{s\brack r}_q},
 \qquad
 \max(0,r-z)\le u\le\min(r,s-z),
\]

and the probability is zero outside that range. The Gaussian binomial \({n\brack k}_q\) counts \(k\)-dimensional subspaces of \(\mathbb F_q^n\).

To count the numerator, first choose \(U\cap K\), then choose the image of \(U\) in the quotient by \(K\). The remaining choices are graphs of linear maps into a space of dimension \(s-z-u\). There are \(q^{(s-z-u)(r-u)}\) such maps. This also proves normalization.

Imposing the root equations adds exactly \(r-u\) to the common equation rank. The case \(u=0\) is discarded on an incoming-nonzero branch. It is not redirected to the zero-state branch: these equations would force the incoming witness state to be zero.

The endpoint conventions matter. For \(z=0\), \(u=r\) with probability one. For \(z\ge s\), every positive-state branch is dead. At incoming state zero, the separate packet condition \(I=J\) applies; the positive-state root action must not replace that condition.

## Composition with the remaining zero equations

Let \(f_0\) designated-fresh inputs and \(d_0\) other required inputs be forced to zero output. Their diagonals are still independent after the root equations. Fresh forms preserve the old-state rank. The other equations can each reduce it by at most one.

Let \(V_h\) be any nonnegative continuation bound, valid uniformly after this zero-output phase. No monotonicity is assumed. The conditional upper action is

\[
 \boxed{
 \sum_{u\ge1}p_{r,z}(u)\,
 \tau^{-(r-u+f_0+d_0)}
 \max_{\max(1,u-d_0)\le h\le u}V_h
 },\qquad 1\le\tau\le q-1.
\]

The factor for the input-only zero equations is the alive diagonal-rank bound. The maximum remains necessary because the root event can correlate the state image with the current source forms. This expression replaces the old maximum over \(u\) by an exact average; it does not replace later conditional geometry by a fresh subspace sample.

The same orientation cannot be averaged again after a zero-output survival test that used it. For a refresh followed by three ordinary packets, backward continuation has the order

\[
 \mathcal T_{\rm GL}
 \bigl(\mathcal T_{\rm plain}
   (\mathcal T_{\rm plain}(\mathcal T_{\rm plain}(V)))\bigr).
\]

These actions can be nonlinear. A positive supersolution bounds their composition; taking a matrix power of a nonlinear action is invalid. The separately owned [packet operator](theory/rank_packet_operator.py) exposes this averaging through `fresh_orientation=True`.

## Exact magnitude and arithmetic cost

At \(q=256\), \(s=4\), and \(\tau=255\), define

\[
 E_{r,z}:=\sum_{u\ge1}p_{r,z}(u)255^{-(r-u)}.
\]

Exact rational evaluation of all twelve cases \(1\le r\le4\), \(1\le z\le3\) gives

\[
 E_{r,z}\le\frac{257}{256}\,255^{-z}.
\]

The largest ratio \(E_{r,z}255^z\) is approximately 1.003875673284. It occurs at \((r,z)=(2,1)\) and \((3,2)\). The multiplicative slack cannot be omitted: these cases exceed \(255^{-z}\).

Charge four additional packed GFNI512 products per orientation refresh. At outer rate \(22/32\), this gives the following conservative product-volume comparison. Boundary-phase savings, instruction scheduling, and coefficient setup are not included.

| Refresh interval | Additional products/source | Total products/source | Increase over \(181/88\) |
|---:|---:|---:|---:|
| 64 symbols | \(1/11\) | \(189/88\) | 4.42% |
| 16 symbols | \(4/11\) | \(213/88\) | 17.68% |

These are processed-byte arithmetic estimates, not measured runtimes. Changing the refresh interval changes the encoder distribution. It is not an exact-map optimization of the archived once-per-region profile.

## A deterministic one-output restriction

There is also a useful restriction that needs no fresh orientation. Let the sole allowed output have label \(j\). If at least one different label \(z\) has both input and output zero, the post-zero state image lies in \(\ker a_z\). The injection column \(c_j\) does not lie there, because

\[
 a_zc_j=\eta(z+j)^3\ne0.
\]

The one-output map \(F_j=I+c_ja_j\) fixes \(c_j\). Hence \(c_j\notin F_jU\). An independent allowed-input form raises rank by exactly one; a dependent form preserves rank. No exceptional diagonal can lower rank. This remains true after earlier zero equations, since they only shrink the state image. The condition strengthens rank continuation as well as excluding a return.

## Replay and scope

```sh
python -B research/workstreams/char2_field_spin/theory/regional_grassmann_gate.py --check
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_regional_grassmann_gate.py
```

The five tests enumerate every subspace of \(\mathbb F_2^4\) and \(\mathbb F_4^4\), compare their intersection counts to the formula, and check the GF256 rational bounds. They also test arbitrary nonmonotone continuations and the one-output geometric restriction. All five passed. No native code or benchmark is included.

The exact local law is established. A whole-class claim still requires a valid support cover and a directed evaluation of its complete bound. The cheaper diagonal-refresh construction is the next candidate to test; it must use its own moment bound rather than the GL distribution above.
