# Return caps from factorial moments

The positive return expansion bounds many output-zero events separately.
Their upper bounds need not be jointly attainable by one distribution of
output weight. The refinement below enforces two consequences of that
shared distribution, without changing the construction or its setup event.

## The finite measure

Use the [translated-feedback construction](HD_FIELD128_TRANSLATED_FEEDBACK_RETURN.md):
\(F=\operatorname{GF}(2^{128})\), \(q=|F|\), state dimension 32,
emission \(A=A_{128}\), and feedback
\(C=A_{64}^{T}\). Fix a degree-16 scalar \(\gamma\) satisfying its principal
invertibility event \(\mathcal G\). Fix a potential input support \(S\),
with \(|S|=j\), and sample \(X\) uniformly from \(F^S\). Define

\[
R=CX,\qquad Y=(I+\gamma AC)X,\qquad
K=|\{i:Y_i=0\}|.
\]

All moments below include the indicator \(R\ne0\); they are not
conditional expectations given a nonzero syndrome. Their total mass is

\[
p_0=\Pr[R\ne0]=1-q^{-\min(j,32)}.
\]

Let \(p_{j,c,e}\) be the proved V44 upper probability for \(R\ne0\)
and zeros on specified sets of \(c\) coordinates in \(S\) and \(e\)
coordinates outside \(S\). Keep its original denominator \(q^j\) and
its exact zero-syndrome subtraction. For \(0\le k\le64\), put

\[
H_k=\sum_{c+e=k}\binom jc\binom{64-j}e p_{j,c,e}.
\]

Invalid set sizes contribute zero. The existing exact zeros for \(c=j\)
and \(e\ge32\) can be omitted. Counting the \(k\)-subsets of the actual
zero set gives

\[
m_k:=\mathbb E\!\left[\binom Kk1_{R\ne0}\right]\le H_k,
\qquad H_0=m_0=p_0.
\tag{1}
\]

This argument uses linearity of expectation, not independence of zero
coordinates. In V50, the coefficients grouped by \(c+e\) are exactly
these \(H_k\).

## An exact support cutoff and moment clipping

If every \(H_k\) is zero, the measure has zero mass and its return moment
is zero. Otherwise let \(M=\max\{k:H_k>0\}\). On \(R\ne0\), one has
\(K\le M\) almost surely. Indeed, when \(M<64\), (1) gives
\(\mathbb E[\binom K{M+1}1_{R\ne0}]=0\), which excludes \(K>M\).
For \(M=64\), the claim is the coordinate-count bound. A positive upper
bound does not assert attainability; only exact zeros justify this cutoff.

For integers \(0\le a\le b\le M\) and \(0\le K\le M\),

\[
\frac{\binom Kb}{\binom Mb}
\le \frac{\binom Ka}{\binom Ma}.
\]

For \(K<b\), the left side vanishes. Otherwise the inequality follows
by writing the additional factors as \((K-i)/(M-i)\le1\).
Taking expectations yields the clipped upper moments

\[
\widehat H_b=
\min_{0\le a\le b}
H_a\frac{\binom Mb}{\binom Ma},\qquad
m_b\le\widehat H_b\le H_b.
\tag{2}
\]

Set \(\widehat H_b=0\) for \(b>M\). Each minimum selects among valid
deterministic upper bounds; it introduces no probabilistic conditioning.

## Two positive upper caps

Fix an exact rational \(0<z<1\). The true return moment is
\(J=\mathbb E[z^{64-K}1_{R\ne0}]\).
The binomial expansion and (2) give

\[
J\le J_{\rm clip}:=
\sum_{b=0}^{M}\widehat H_b z^{64-b}(1-z)^b.
\tag{3}
\]

For \(1\le k\le M\), define

\[
T_k=\min_{0\le b\le k}
\frac{\widehat H_b}{\binom kb}.
\]

On \(K\ge k\), one has \(\binom Kb\ge\binom kb\), so (1)--(2)
imply \(\Pr[R\ne0,K\ge k]\le T_k\). The exact telescoping identity

\[
z^{64-K}=z^{64}+(1-z)\sum_{k=1}^{K}z^{64-k}
\]

therefore gives a second cap,

\[
J\le J_{\rm tail}:=
z^{64}p_0+(1-z)\sum_{k=1}^{M}z^{64-k}T_k.
\tag{4}
\]

All coefficients in (3)--(4) are nonnegative. The implemented return cap
is \(\min(J_{\rm V50},J_{\rm clip},J_{\rm tail})\), retaining every
previous basic cap through \(J_{\rm V50}\). It cannot worsen V50.

## Integration and scope

[V57](theory/hd_field128_return_tail_moments_v57.py) implements these
identities with exact fractions and pins the unchanged V50 provider.
[Its small tests](theory/hd_field128_return_tail_moments_v57_tests.py)
include independent finite distributions and a binary linear counting oracle.
No coefficient table or numerical search is part of the provider.

Only the return entry changes: divide the selected cap by \(q^{32}-1\)
once. Quiet, birth, live, refresh density, and the empty-packet live cap
are unchanged. The same principal-failure probability
\(8\binom{64}{32}/(q-2^{64})\) remains a once-global term outside a
future conditional code union. This refinement adds no gamma exclusion.

The remote exact scan reports strict local-cap improvements. They do not
establish a selected-weight bound or a distance certificate. Those claims
still require the chronological coefficient calculation, complete profile
coverage, directed arithmetic, and the appropriate whole-code union.
