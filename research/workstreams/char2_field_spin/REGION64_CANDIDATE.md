# Exact local coordinates for the period-64 regional candidate

**Profile status:** the algebra below is valid, but GRS[32,24] with this
four-state regional inner is rejected at K20. The 37-row window event in
[RATE_DESIGN_REVIEW.md](RATE_DESIGN_REVIEW.md) has failure probability above
2^-38 regardless of the label period or coefficient assignments. Retain
this reference only as a possible building block for a changed outer,
state dimension, or route; its algebra tests are not a distance certificate.

Four local basis choices give 64 distinct injection labels while retaining the fixed width-16, four-state GF16 arithmetic.
The scalar reference verifies the global map, local map, and their adjoints, including the first and last state conversions.
It establishes algebraic equivalence; it does not certify the complete encoder or report native performance.

Use the AES field \(F=\mathrm{GF}(256)\) and its embedded subfield \(H=\mathrm{GF}(16)\).
The reference orders \(H\) with the same polynomial-basis embedding as `Subfield16.h`.
It does not interpret abstract GF16 nibbles as AES-field bytes.
For \(u\in H\), define

\[
v(u)=(1,u,u^2,u^3)^T,\qquad
C_0=(v(u))_{u\in H},\qquad A_0=C_0^T.
\]

The power sums through degree six vanish on \(H\), so \(C_0A_0=0\).
The phase offsets, as AES-field bytes, are

\[
(\beta_0,\beta_1,\beta_2,\beta_3)=(0,2,6,4).
\]

The four cosets \(\beta_j+H\) are disjoint.
A step uses sixteen labels from its coset; consecutive steps cycle through these four offsets.
Thus the label pattern has period 64 individual symbols.
Consecutive offset differences alternate between 2 and 4, including the wrap from the fourth phase to the first.

For \(\beta\in F\), define

\[
T_\beta=
\begin{pmatrix}
1&0&0&0\\
\beta&1&0&0\\
\beta^2&0&1&0\\
\beta^3&\beta^2&\beta&1
\end{pmatrix}.
\]

The binomial identities in characteristic two give
\(T_\beta v(u)=v(u+\beta)\),
\(T_\beta^2=I\), and
\(T_\beta T_\gamma=T_{\beta+\gamma}\).

For each region \(r\), fix an invertible orientation \(B_r\in\mathrm{GL}(4,F)\).
The proposed setup samples these orientations independently and uniformly; the equivalence below holds for every fixed invertible choice.
At phase \(\beta\), the global step is

\[
Y=X+A_0T_\beta B_r^{-1}Q,
\qquad
Q'=Q+B_rT_\beta C_0X.
\]

The global state \(Q\) persists across region boundaries.
The optional input diagonal is applied before this step; its adjoint is applied after the transpose step.
The ordinary encoder starts at \(Q=0\) and discards the final state.

Set \(S=T_\beta B_r^{-1}Q\).
The local step becomes

\[
Y=X+A_0S,\qquad S'=S+C_0X.
\]

Within a region, changing phase from \(\beta\) to \(\gamma\) changes the local state by \(T_{\beta+\gamma}\).
Its adjoint uses \(T_{\beta+\gamma}^T\), not \(T_{\beta+\gamma}\).

## Four-product phase change

The displayed matrix gives a direct five-product formula.
Factoring it into two coordinate shears removes one product.
For old state \((q_0,q_1,q_2,q_3)\) and phase difference \(\delta\), compute

\[
b=q_1+\delta q_0,\quad d=q_3+\delta q_2,
\quad c=q_2+\delta^2q_0,\quad d'=d+\delta^2b.
\]

The new state is \((q_0,b,c,d')\).
The last operation intentionally uses the new \(b\); its extra term is \(\delta^3q_0\).
The adjoint computes

\[
a=q_0+\delta^2q_2,\quad b=q_1+\delta^2q_3,
\quad c=q_2+\delta q_3,\quad a'=a+\delta b,
\]

and returns \((a',b,c,q_3)\).
Each direction uses four fixed GF256 state products and four XORs.
The coefficients are prepared once.
For the two phase differences, only the constants 2, 4, and 16 are needed.
With 128-bit state symbols, this is four GFNI128 instructions and four XOR128 instructions per internal phase boundary.

The direct width-16 emission/moment calculation uses 24 GFNI512 instructions per step; the existing four-state additive-FFT arithmetic uses twelve.
An input diagonal adds four GFNI512 instructions per step.
These counts exclude lane shuffles, sums, routing, and stores.
They are instruction counts, not throughput predictions.

## Region and external boundaries

Every region begins in phase zero.
If the preceding region ends in phase \(\beta\), the new local state is

\[
S_{\mathrm{new}}=B_{r+1}^{-1}B_rT_\beta S_{\mathrm{old}}.
\]

The entire matrix can be prepared once for that boundary.
Its transpose is applied after reversing the first step of the new region.
A four-record packed matrix multiply needs four GFNI512 products per region boundary; no random matrix multiply is required at internal phase boundaries.

For a nonzero external initial state, the first conversion is \(S=B_0^{-1}Q\).
At the end, recover \(Q=B_{\mathrm{last}}T_{\beta_{\mathrm{last}}}S\).
Consequently a terminal global dual state first becomes
\(T_{\beta_{\mathrm{last}}}^TB_{\mathrm{last}}^TH\) in the local adjoint.
After reversing the first region, apply \((B_0^{-1})^T\) to report the initial global dual state.
The zero-initial, discarded-final encoder can omit these external conversions because their states are zero.

The test covers regions of one through seven steps, so a boundary can follow any phase.
It verifies translation identities for every field offset, all 64 distinct labels, \(C_0A_0=0\), and 72 global/local/adjoint cases.
Those cases include nonzero initial and terminal states, nontrivial regional orientations, optional nonzero GF256 input diagonals, and the extended inner-product identity.

```sh
python -B research/workstreams/char2_field_spin/theory/region64_reference.py --check
```

The next gate remains a complete multi-row distance argument for the higher-rate construction.
Random orientation once per region does not supply fresh uniform nonzero state at every inner step.
