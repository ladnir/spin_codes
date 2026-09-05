# Region-shuffled transpose followed by one accumulator

## Conclusion

The region permutations remove the adjacent-row cancellation of a pure
transpose.  For random outer blocks, the resulting first moment has an exact
finite transfer with only two accumulator states after the number of active
blocks is fixed.

The fixed-active-block asymptotics are promising.  At rate one half and target
distance \(0.02\), the sufficient logarithmic block constant rises toward
approximately \(6.9514\) as the number of active blocks grows.  This is the
same constant as the uniform-interleaver accumulator theorem.

The equality is not yet a complete theorem.  The present calculation proves
the finite transfer and every fixed-active-block asymptotic class.  A uniform
bound for a growing sublinear number of active blocks remains open.

## Factored permutation

Let \(N=LB\).  Arrange the outer word as an \(L\)-by-\(B\) binary matrix.
Row \(i\) is the output of outer block \(i\).

For each region \(j\in\{1,\ldots,B\}\), setup independently samples
\(\pi_j\gets S_L\).  Transpose the matrix, and apply \(\pi_j\) within region
\(j\).  Serialize the regions in the order \(1,\ldots,B\), and apply the
length-\(N\) accumulator.

A per-block coordinate permutation is distributionally redundant for uniform
random outer injections and is omitted from this proof model.  Such a
permutation can matter for a structured outer.

## One-region accumulator kernel

Fix a region of length \(L\) containing \(a\) ones.  The region permutation
makes its input uniform on the weight-\(a\) Hamming slice.

Set \(r:=\lceil a/2\rceil\).  Define

\[
  F_{L,a}(z)
  :=
  \frac1{\binom La}
  \sum_{w=0}^L
  \binom{w-1}{r-1}
  \binom{L-w}{a-r}z^w,
  \tag{1}
\]

where the \(a=0\) term is \(F_{L,0}(z):=1\), and infeasible binomial
coefficients are zero.  Equation (1) is the output-weight generating function
when the accumulator enters the region in state zero.

If the entering state is one, every output bit is complemented.  The
generating function becomes

\[
  \widetilde F_{L,a}(z):=z^L F_{L,a}(z^{-1}).
  \tag{2}
\]

The ending state equals the entering state XOR \((a\bmod2)\).  Therefore the
two-state region kernel is

\[
  K_{L,a}(z)
  :=
  \begin{cases}
    \begin{pmatrix}
      F_{L,a}(z)&0\\
      0&\widetilde F_{L,a}(z)
    \end{pmatrix},&a\text{ even},\\[4mm]
    \begin{pmatrix}
      0&F_{L,a}(z)\\
      \widetilde F_{L,a}(z)&0
    \end{pmatrix},&a\text{ odd}.
  \end{cases}
  \tag{3}
\]

Rows index the entering state, and columns index the ending state.

## Active-block reduction

Fix a message whose support contains \(q\) outer blocks.  Under independent
uniform random injections, the corresponding \(q\) active rows are independent
and uniform on \(\mathbb F_2^B\setminus\{0\}\).

First relax each row to a uniform element of \(\mathbb F_2^B\).  In each
region, the number of one bits is then an independent
\(\operatorname{Bin}(q,1/2)\) random variable.  Define the averaged kernel

\[
  M_{L,q}(z)
  :=
  2^{-q}\sum_{a=0}^q\binom qa K_{L,a}(z).
  \tag{4}
\]

The regions are independent under the relaxed law.  Hence the complete
output-weight generating function is

\[
  e_0^{\mathsf T}M_{L,q}(z)^B\boldsymbol1,
  \qquad
  e_0:=(1,0)^{\mathsf T}.
  \tag{5}
\]

The probability that all \(q\) relaxed rows are nonzero is
\((1-2^{-B})^q\).  For every event \(E\),

\[
  \Pr[E\mid\text{all active rows are nonzero}]
  \le
  \frac{\Pr[E]}{(1-2^{-B})^q}.
  \tag{6}
\]

This conditioning relaxation is explicit.  It is not an independence
assumption about the realized outer code.

## Exact finite first-moment bound

Let each random outer block have rate \(R\), and define

\[
  \beta_{B,R}
  :=
  \frac{2^{RB}-1}{1-2^{-B}}.
  \tag{7}
\]

For an integer distance threshold \(D\) and any \(z\in(0,1)\), equations
(5)--(7) give

\[
  \mathbb E[Z_D]
  \le
  z^{-D}
  \sum_{q=1}^L
  \binom Lq
  \beta_{B,R}^{q}
  e_0^{\mathsf T}M_{L,q}(z)^B\boldsymbol1.
  \tag{8}
\]

The expectation is over the random outer injections and region permutations.
Thus equation (8) is a computable finite upper
bound for the declared ensemble.  It uses \(O(L)\) two-by-two matrix powers
after the region polynomials are available.

## Fixed-active-block exponent

Fix \(q\) while \(L\to\infty\), and set \(z=e^{-\theta/L}\).  Conditional on
\(a\) region ones, the normalized occupation converges to a beta random
variable.  For \(a\ge1\), define

\[
  k_a:=\left\lceil\frac a2\right\rceil,
  \qquad
  f_a(\theta):={}_1F_1(k_a;a+1;-\theta),
  \qquad
  g_a(\theta):={}_1F_1(a+1-k_a;a+1;-\theta).
  \tag{9}
\]

For \(a=0\), set \(f_0(\theta):=1\) and
\(g_0(\theta):=e^{-\theta}\).  Replace \(F\) and \(\widetilde F\) in (3)
by \(f_a\) and \(g_a\), and call the resulting binomial average
\(\mathcal M_q(\theta)\).  Let \(\Lambda_q(\theta)\) be its Perron root.

For target relative distance \(\delta\), define

\[
  E_q(R,\delta)
  :=
  Rq+inf_{\theta\ge0}
  \left[
    \log_2\Lambda_q(\theta)
    +\frac{\theta\delta}{\ln2}
  \right],
  \qquad
  \kappa_q:=-E_q.
  \tag{10}
\]

If \(\kappa_q>0\), the contribution from exactly \(q\) active blocks tends
to zero for \(B=c\log_2N\) whenever

\[
  c\kappa_q\ge q.
  \tag{11}
\]

The equality case closes because
\(\binom Lq=O(N^q/B^q)\).

## Numerical fixed-q results

The script `analyze_region_shuffled_transpose_accumulator.py` evaluates
(9)--(11) with binary64 arithmetic.  At \(R=1/2\) and \(\delta=0.02\), it
reports:

| active blocks \(q\) | \(\kappa_q/q\) | required \(c=q/\kappa_q\) |
|---:|---:|---:|
| 1 | 0.3045563 | 3.28347 |
| 2 | 0.2164337 | 4.62035 |
| 4 | 0.1784814 | 5.60282 |
| 8 | 0.1607830 | 6.21956 |
| 16 | 0.1522267 | 6.56915 |
| 32 | 0.1480186 | 6.75591 |
| 64 | 0.1459318 | 6.85252 |
| 128 | 0.1448926 | 6.90166 |
| 256 | 0.1443740 | 6.92645 |

The uniform-interleaver theorem has

\[
  \lambda_A(1/2,0.02)
  =0.14385619\ldots,
  \qquad
  1/\lambda_A=6.95138667\ldots.
\]

The fixed-\(q\) values approach this constant from below.  The same pattern
appears at other distances:

| \(\delta\) | required \(c\), \(q=256\) | uniform-interleaver \(c\) |
|---:|---:|---:|
| 0.030 | 12.9744 | 13.0714 |
| 0.035 | 20.4142 | 20.6627 |
| 0.040 | 42.5026 | 43.6197 |

These values are diagnostics.  They are not outward-rounded certificates.

## The remaining theorem gap

Equation (8) is an exact finite upper bound, apart from its stated conditioning
relaxation.  Equations (10) and (11) close every fixed \(q\).  They do not
justify interchanging the limits \(N\to\infty\) and \(q\to\infty\).

A complete logarithmic-block theorem needs three additional uniform bounds.

1. Prove
   \[
     \kappa_q(R,\delta)
     \ge q\lambda_A(R,\delta)-o(q)
   \]
   uniformly as \(q\to\infty\).
2. Control the finite-\(L\) beta approximation when \(q=o(L)\).
3. Close the linear regime \(q=\Theta(L)\), where the desired margin is
   \(\Theta(N)\).

The numerical results support the same sufficient constants as the uniform
interleaver:

\[
  c\lambda_A(R,\delta)\ge1.
  \tag{12}
\]

Equation (12) is a theorem target, not yet a proved statement for the
region-shuffled transpose.

## Status

- **Proved:** the conditional region kernel (1)--(5).
- **Proved:** the conditioning relaxation (6) and finite first moment (8).
- **Proved:** the fixed-active-block exponent reduction (9)--(11).
- **Diagnostic:** the binary64 constants in the two tables.
- **Conjectured:** the uniform condition (12).
- **Open:** the growing-sparse and linear-active-block ranges.
