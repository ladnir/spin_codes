# Weight-coupled fixed-occupation transfer

## Result

The uniform likelihood comparison in
`SINGLE_SAMPLED_BA_RM2SUB_D11.md` maximizes over the BA row weight before it
applies the RM2Sub transfer. That order charges the endpoint likelihood loss
to every row. This note keeps a row-weight fugacity inside the transfer.

Consider the one-sampled Golay--BA-3/RM2Sub-S19 ensemble from that theorem.
Condition on the proved good-spectrum event \(\mathcal G_B\). For every fixed
integer \(Q\ge3\), the weight-coupled transfer proves

\[
  \Pr[Z_{\lfloor0.11N\rfloor,Q}>0]
  \le
  \exp\!\left(-0.0008689\,QB+o_Q(B)\right)
  \tag{1}
\]

under the schedule

\[
  B=\frac{39}{4}\log_2N+O(1).
  \tag{2}
\]

The probability in (1) is over the route permutations and RM2Sub
multipliers. The sampled BA code is fixed after conditioning on
\(\mathcal G_B\). The existing fixed-occupation certificates cover
\(Q=1,2\) with larger margins under (2).

## Transfer with row fugacities

Fix \(Q\) active outer rows. Let row \(j\) have weight \(w_j\), and set

\[
  x_j:=\frac{w_j}{B}\in[0.104,0.896].
  \tag{3}
\]

The independent coordinate permutation of row \(j\) makes that row uniform
on its weight slice. Introduce a positive fugacity \(u_j\) for its one-bit
choice in each region.

Let

\[
  P:=
  \begin{pmatrix}
    0&1\\
    r&1-r
  \end{pmatrix},
  \qquad
  r:=\frac1{2^{19}-1},
  \qquad
  H_u:=I+uP.
  \tag{4}
\]

Rows and columns index the zero and live RM2Sub state classes. For fixed
\(\theta>0\), let \(K_a(\theta)\) be the continuum region matrix for \(a\)
impulses from `RM2SUB_FIXED_OCCUPATION_CONTINUUM.md`. Define

\[
  T_Q(\boldsymbol u,\theta)
  :=
  \sum_{\boldsymbol b\in\{0,1\}^Q}
  \left(\prod_{j=1}^Q u_j^{b_j}\right)
  K_{\sum_jb_j}(\theta).
  \tag{5}
\]

Equation (5) is a matrix-valued polynomial with nonnegative coefficients.
Coefficient extraction enforces each row weight. For a fixed weight tuple,
the fixed-\(Q\) continuum reduction gives

\[
\begin{split}
 &\mathbb E\!\left[z^{\operatorname{wt}(Y)}
       \mid w_1,\ldots,w_Q\right]\\
 &\quad\le
 \left(\prod_{j=1}^Q\binom B{w_j}^{-1}u_j^{-w_j}\right)
 e_0^{\mathsf T}T_Q(\boldsymbol u,\theta)^B\boldsymbol1
 \exp(o_Q(B)).
 \tag{6}
\end{split}
\]

The expectation in (6) is over the route and inner randomness. The factor
\(\binom B{w_j}^{-1}\) is the probability of the selected row slice.

## A product norm bound

Set

\[
  p:=\frac{2^{18}}{2^{19}-1},
  \qquad
  \theta:=\frac{Q\tau}{p},
  \qquad
  D_s:=\operatorname{diag}(1,s),
  \tag{7}
\]

where \(0<s<1\) and \(\tau>0\). Define

\[
  G_Q(s,\tau)
  :=\sup_{0\le y\le1}
  \left\{-\tau y+left(1+\frac1Q\right)
  \ln(1-y+y/s)\right\}.
  \tag{8}
\]

The beta-integral identity used in the uniform fixed-occupation proof now
gives

\[
  T_Q(\boldsymbol u,Q\tau/p)
  \le
  e^{QG_Q(s,\tau)}
  \frac1{Q!}\sum_{\pi}
  D_sH_{u_{\pi(1)}}D_s\cdots H_{u_{\pi(Q)}}D_s.
  \tag{9}
\]

The sum in (9) is over all permutations of \([Q]\). Let
\(v=(1,v_1)^{\mathsf T}\), where \(v_1>0\), and use the weighted row norm

\[
  \lVert A\rVert_v:=\max_i\frac{(Av)_i}{v_i}.
  \tag{10}
\]

Since \(D_s\le I\) entrywise, submultiplicativity in (9) yields

\[
  \lVert T_Q(\boldsymbol u,Q\tau/p)\rVert_v
  \le
  e^{QG_Q(s,\tau)}\prod_{j=1}^Q m_v(u_j),
  \tag{11}
\]

where

\[
\begin{split}
  m_v(u):=\max\{&1+usv_1,\\
  &ur/v_1+s(1+u(1-r))\}.
  \tag{12}
\end{split}
\]

For \(Q\ge3\), the logarithm in (8) is nonnegative. Therefore

\[
  G_Q(s,\tau)\le G_3(s,\tau).
  \tag{13}
\]

Equations (11)--(13) provide one bound for every fixed \(Q\ge3\).

## Coupling the BA spectrum to the transfer

On \(\mathcal G_B\), the selected BA spectrum and the certified majorant
satisfy

\[
  \frac1B\ln\frac{A_{\mathcal O_B}(w)}{\binom Bw}
  \le
  \widehat a_{\rm BA}(x)-h(x)+o(1),
  \qquad x=\frac wB.
  \tag{14}
\]

For \(u>0\), define the weight-coupled row exponent

\[
  R(x;u):=
  \widehat a_{\rm BA}(x)-h(x)-x\ln u+\ln m_v(u).
  \tag{15}
\]

Unlike the earlier likelihood comparison, equation (15) permits a different
fugacity for each BA weight segment. For a fixed fugacity and one affine
segment of \(\widehat a_{\rm BA}\), the map \(x\mapsto R(x;u)\) is convex.
It is therefore enough to check both segment endpoints.

Combining (6), (11), and (14) gives

\[
\begin{split}
 \ln\mathbb E[Z_{\lfloor0.11N\rfloor,Q}]
 \le{}&Q\ln L\\
 &+B\sum_{j=1}^Q
 \left[R(x_j;u_j)+G_3(s,\tau)+\frac{0.11\tau}{p}\right]
 +o_Q(B).
 \tag{16}
\end{split}
\]

The order of choices in (16) is

\[
  \sup_{x_1,\ldots,x_Q}
  \inf_{u_1,\ldots,u_Q}.
  \tag{17}
\]

The fugacities are coefficient bounds chosen after the weight tuple is
fixed. They are not setup randomness.

## Outward certificate

`certify_golay_ba_rm2sub_weight_coupled_fixed.py` uses the rational common
witness

\[
  s=\frac{127}{250},
  \qquad
  \tau=\frac{133}{125},
  \qquad
  v_1=\frac3{1600}.
  \tag{18}
\]

For each of the 31 affine BA segments, the verifier chooses one rational
fugacity and checks both endpoints with 100-digit outward intervals. The 62
checks prove

\[
  R(x;u)+G_3(s,\tau)+\frac{0.11\tau}{p}
  <-0.07196093457399053
  \tag{19}
\]

for every \(x\in[0.104,0.896]\). Binary64 arithmetic proposes the
fugacities. It does not decide any accepted inequality.

Under (2), choosing the active rows costs

\[
  \frac{\ln2}{39/4}QB+o(QB).
  \tag{20}
\]

The largest outward endpoint after adding (20) is

\[
  -0.0008689160550218037.
  \tag{21}
\]

Equations (16)--(21) prove (1). The factor \(B^2\) in the selected-spectrum
event contributes \(o_Q(B)\) for each fixed \(Q\).

## Scope

This note proves a uniform numerical bound for every fixed integer
\(Q\ge3\). It does not allow \(Q=Q_N\to\infty\). The existing exact
four-state certificate still controls growing sparse occupations. The common
fixed-occupation norm bound is tightest at \(Q=3\). A further reduction of
the block constant therefore starts with a sharper \(Q=3\) transfer.
