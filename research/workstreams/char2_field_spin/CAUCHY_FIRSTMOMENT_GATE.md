# Fixed-word sparse screen for the Cauchy diagonal-state inner

The requested coupled first-moment envelope does not give a useful bound, including at one active outer row. Exact regional counting and the exact one-input return weight yield a directed one-row certificate. Its margins are 107, 105, and 103 bits at K16, K18, and K20. The multi-row class remains open.

The exploratory margins below are floating screens. The final certificate section uses outward-rounded intervals. A negative screening margin means that the tested upper expression exceeds one; it is not a construction-failure lower bound.

## Fixed-word packet envelope

The packet has sixteen inputs and a four-symbol state over GF256. It refreshes the state with four independent nonzero diagonals. The input diagonal is independently sampled at every input coordinate. Fix the incoming nonzero state and the input values before exposing these packet diagonals. Write \(w\) for the input weight and \(0<y\le1\) for the output-weight mark.

For a chosen set of forced zero outputs, let \(i\) count its active-input positions and \(z\) its zero-input positions. Superregular Cauchy emission gives probability at most \(255^{-(i+z)}\) when \(z\le3\). Four zero-input zero outputs are impossible for a nonzero incoming state.

The identity

\[
 y^{\operatorname{wt}(Y)}
 =\prod_{j=1}^{16}\bigl[y+(1-y)\mathbf1_{Y_j=0}\bigr]
\]

therefore gives the marginal upper bound

\[
 M_w(y):=\sum_{i=0}^{w}\sum_{z=0}^{\min(3,16-w)}
 \binom wi\binom{16-w}z
 \frac{y^{16-i-z}(1-y)^{i+z}}{255^{i+z}}.
\]

For a return to zero state, impose up to four further pivots using the remaining input diagonals. This multiplies each term by \(255^{-\min(w-i,4)}\). Omit the terms \(i=w\): no active input then remains outside the forced zero outputs, and the invertible packet map excludes a return. Call the resulting sum \(R_w(y)\), and set \(R_0=0\).

The two-state upper matrix is

\[
 T_w(y)=\begin{pmatrix}
  \mathbf1_{w=0}+\mathbf1_{w\ge5}y^w/255^4 &
  \mathbf1_{w>0}y^w\\
  R_w(y)&M_w(y)
 \end{pmatrix}.
\]

Its second diagonal entry is a marginal bound, so returns are also counted there. That overcount preserves an entrywise upper bound. All statements hold uniformly over fixed incoming values; no uniform-state law is assumed.

## Outer support counting uses the input diagonal only once

For a fixed exact support \(S\) of a GRS[32,22] row, let \(a_S\) be the number of nonzero words with that support. The deterministic MDS bound is

\[
 a_S\le256^{|S|-10}
 \le256^{-10}(256/255)^{32}255^{|S|}.
\]

Apply this support multiplicity bound to the packet estimate, which is uniform over all fixed nonzero coordinate values. The input diagonal is consumed in the packet estimate. It is not averaged into a counting measure and then resampled.

The relaxed measure permits a selected row to be zero. This only increases the bound. The optional three-state implementation removes the globally all-zero relaxed source without cancellation. It does not enforce that every selected row is nonzero.

## Executed coupled screen

For \(A\) selected rows, use

\[
 B(x,y):=\sum_{w=0}^{16}\binom{16}w
 (255x)^w(1+x)^{16-w}T_w(y).
\]

The factor \(1+x\) distinguishes an unselected zero position from a selected position whose field value is zero. Collapsing all regional selection marks to \(x>0\) gives

\[
 \binom LA^{1-32}
 256^{-10A}(256/255)^{32A}
 \frac{x^{-32A}y^{-e}}{255}
 e_0B(x,y)^{32L/16}\mathbf1.
\]

At \(K=2^{20}\), the padded geometry is \(L=47680\), \(N=1525760\), and \(e=152575\). Executed margins in bits are:

| Active rows | Requested local return bound | One-input bound using weight twelve |
|---:|---:|---:|
| 1 | -142.049 | -131.771 |
| 16 | -1676.441 | -1459.678 |
| 64 | -6277.567 | -5426.170 |
| 256 | -23642.851 | -20684.281 |
| 1024 | -82860.656 | -76962.234 |
| 4096 | -220327.547 | -220315.309 |

Removing the globally all-zero source changes the first requested margin only to -142.046 bits. It does not resolve the failure.

The coupled coefficient bound also permits selected coordinates to be unevenly distributed among regions. Exact regional extraction with the stronger one-input return bound gives -1117.661 bits at \(A=16\) and -5112.818 bits at \(A=64\). These remain unusable.

## A stronger return bound for one input

Every injection column \(c_j\) has four nonzero coordinates. For one active input, a return satisfies

\[
 D_{\rm state}R=c_jd_jX_j.
\]

If \(R\) has a zero coordinate, the equation is impossible. Otherwise, each input-diagonal value fixes all four state-diagonal values. Thus the return probability is at most \(255^{-4}\).

On a return, the output is one input coordinate plus the emission of a nonzero state. MDS emission has weight at least thirteen, so the output has weight at least twelve. Consequently,

\[
 \mathbb E\bigl[y^{\operatorname{wt}(Y)}\mathbf1_{R'=0}\bigr]
 \le255^{-4}y^{12}.
\]

Use the minimum of this expression and the earlier \(R_1(y)\). This is a joint exposure of the input and state diagonals. It is not an extra independent factor multiplied into the earlier return proof.

The same argument cannot be generalized by assuming that the injection matrix is superregular. Some two-by-two and three-by-three injection minors vanish in phases zero and two.

## Exact one-row regional transfer

Let \(H=L/16\), and put \(m=M_0(y)\) and

\[
 g:=\frac{1-m^H}{H(1-m)}.
\]

For a region whose outer coordinate is zero, the matrix is \(E=\operatorname{diag}(1,m^H)\). For a nonzero outer coordinate, average its uniform packet position to obtain

\[
 U=\begin{pmatrix}
 0&yg\\R_1(y)g&M_1(y)m^{H-1}
 \end{pmatrix}.
\]

The coefficient of \(z^w\) in \(e_0(E+zU)^{32}\mathbf1\) sums all exact regional support patterns of size \(w\). The number of projective GRS words on any particular support of size \(w\ge11\) is

\[
 n_w:=\sum_{j=0}^{w-11}(-1)^j\binom{w-1}j256^{w-11-j}.
\]

Multiply that coefficient by \(n_w\), sum over \(w\), and multiply by \(Ly^{-e}\). The implementation evaluates the regional coefficients by positive matrix-polynomial recursion.

With the requested return envelope, this exact-source one-row screen gives -11.592 bits. With the stronger one-input bound it gives +104.110 bits at \(y\approx.999716686\). The simple rational witness \(y=9997/10000\) gives approximately +103.321 bits. These are prospective one-row witnesses, not a directed certificate or a multi-row result.

Using the dominated full-field source rather than the exact GRS support enumerator gives only +37.784 bits. Exact support counting matters even after the regional route and one-input return are handled.

## Replay and next step

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_firstmoment_tests.py
python -B research/workstreams/char2_field_spin/theory/cauchy_firstmoment_screen.py
python -B research/workstreams/char2_field_spin/theory/cauchy_firstmoment_screen.py --nonzero-only
python -B research/workstreams/char2_field_spin/theory/cauchy_firstmoment_screen.py --regional-one --one-input-upgrade
python -B research/workstreams/char2_field_spin/theory/cauchy_firstmoment_screen.py --regional --one-input-upgrade --active 1 16 64
```

Seven tests pass. They include exhaustive GF4 checks of the marginal diagonal bound and the joint one-input return lemma, exact selection counting, global-nonzero tracking, and exact MDS support counts.

## Directed one-row certificate

The actual packet maps give a stronger identity than the weight-twelve bound used in the preceding screens. Let \(A\) denote the unscaled emission matrix and \(C_j\) the injection column at packet label \(u_j\). Their exact pairing is

\[
 (AC)_{i,j}=(u_i+u_j)^3.
\]

Write \(S=D_{\rm state}R\) for the refreshed state, and include the nonzero input diagonal in \(X_j\). On a return from a one-input packet, \(S=C_jX_j\), hence

\[
 Y_i=[\mathbf1_{i=j}+\eta(u_i+u_j)^3]X_j.
\]

At \(i=j\), the coefficient is one. At every other packet coordinate it is nonzero because the labels are distinct and \(\eta\ne0\). Every such return therefore has output weight exactly sixteen. Each \(C_j\) has four nonzero entries. For an incoming state with full support, each input-diagonal value determines exactly one choice of all four state diagonals. An incoming state with a zero coordinate cannot return. The uniform conditional bound is thus

\[
 R_1(y)\le y^{16}/255^4.
\]

The directed replay substitutes this expression into \(U\). It evaluates the geometric average positively as

\[
 g=H^{-1}\sum_{j=0}^{H-1}M_0(y)^j,
\]

without subtracting nearly equal values. Exact integer MDS support counts weight a positive 32-step polynomial recursion. The initial state is zero and the final state is unrestricted. The existing 80-digit decimal interval type rounds every operation outward.

The frozen rational marks and resulting upper bounds are:

| Requested input size | Padded rows | Rational \(y\) | Failure upper, rounded upward here | Certified margin |
|---:|---:|---:|---:|---:|
| \(2^{16}\) | 3008 | 93553/94000 | \(5.197\cdot10^{-33}\) | 107 bits |
| \(2^{18}\) | 11968 | 373553/374000 | \(1.991\cdot10^{-32}\) | 105 bits |
| \(2^{20}\) | 47680 | 9997/10000 | \(7.853\cdot10^{-32}\) | 103 bits |

Each bound covers every nonzero source with exactly one active outer row, including the union over that row's location and all its projective words. The bad-output threshold is \(\lceil N/10\rceil-1\). The source padding only enlarges the counted source space.

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_q1_certificate.py
python -B research/workstreams/char2_field_spin/theory/cauchy_q1_certificate_tests.py
```

The six certificate tests include independent Fraction enclosure of a complete small transfer, independent support-count inclusion--exclusion, positive geometric sums, and all 64 actual packet return columns. They also check 192 scalar returns with weight sixteen and replay all three frozen witnesses.

This is a complete certificate for the one-active-row class, not for the whole code. Multi-row progress still requires joint return bounds for two and three active inputs, together with stronger rowwise source counting. The existing coupled first-moment relaxation does not close those classes.
