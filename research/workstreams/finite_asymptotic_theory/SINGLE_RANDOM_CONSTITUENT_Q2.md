# Two active rows under one random constituent

## Question

The occupation-one calculation does not expose reuse of the outer
constituent. Two active rows are the first case where reuse can correlate
their outer codewords.

This calculation keeps that correlation. It covers one random
\([256,128]\) constituent reused in all rows, followed by the structured
routing and RandomStepConv-M22.

The finite parameters are

\[
 B=256,\quad K=128,\quad L=8192,\quad N=2^{21},
 \quad D=\lceil 0.109N\rceil=228590.
\]

## Probability space

Sample one linear injection

\[
 G\gets\operatorname{Inj}(\mathbb F_2^K,\mathbb F_2^B).
\]

Use the sampled map \(G\) in all \(L\) outer rows. Independently sample a
coordinate permutation for every row and a row permutation for every
transposed region. Sample the RandomStepConv-M22 maps once.

The bad-message count is random over this complete setup. The calculation
below bounds its expectation for messages with exactly two active rows.

## Rank split

Fix two row positions. Write their nonzero local messages as \(u\) and
\(v\). Over \(\mathbb F_2\), exactly two cases occur.

- If \(u=v\), the pair has rank one.
- If \(u\ne v\), the pair has rank two.

Set

\[
 M:=2^K-1,\qquad P:=2^B-1,\qquad C_w:=\binom Bw.
\]

In the rank-one case, \(G(u)\) is uniform over the \(P\) nonzero ambient
words. The expected number of local-message pairs in weight shell \((w,w)\)
is therefore

\[
 \frac{M}{P}C_w.
 \tag{1}
\]

In the rank-two case, \((G(u),G(v))\) is a uniform ordered pair of distinct
nonzero ambient words. The expected number in shell \((a,b)\) is

\[
 \frac{M(M-1)}{P(P-1)}
 \left(C_aC_b-\mathbf 1_{a=b}C_a\right).
 \tag{2}
\]

Equation (2) retains the shared-code dependence. Replacing its factor by a
product of two independent expected spectra would use a different ensemble.

## Conditional routing transfer

Fix ambient words of weights \(a\) and \(b\). Their independent
row-coordinate permutations produce independent uniform supports of those
weights. This statement remains true when the two ambient words are equal.

For \(q\in\{0,1,2\}\), let \(R_q(z)\) be the average transfer matrix for one
length-\(L\) region that contains exactly \(q\) active input bits. The
average includes the region permutation and the RandomStepConv maps.

Let \(P^{(n)}_{a,b}(z)\) be the average ordered product through the first
\(n\) regions. The average ranges over independent uniform supports of sizes
\(a\) and \(b\) in those regions. Initialize

\[
 P^{(0)}_{0,0}(z)=I.
\]

Invalid subscripts denote the zero matrix. The exact recurrence is

\[
\begin{aligned}
P^{(n+1)}_{a,b}
={}&\frac{(n+1-a)(n+1-b)}{(n+1)^2}P^{(n)}_{a,b}R_0\\
&+\frac{a(n+1-b)}{(n+1)^2}P^{(n)}_{a-1,b}R_1\\
&+\frac{(n+1-a)b}{(n+1)^2}P^{(n)}_{a,b-1}R_1\\
&+\frac{ab}{(n+1)^2}P^{(n)}_{a-1,b-1}R_2.
\end{aligned}
\tag{3}
\]

The recurrence preserves the order of the transfer matrices. It does not
assume that \(R_0,R_1,R_2\) commute.

Let \(\Phi_{a,b}(z)\) be the scalar moment obtained from
\(P^{(B)}_{a,b}(z)\), starting in the zero state and summing both final
states. For the tested tilt set

\[
 \mathcal Z:=
 \left\{\exp(-\exp(-8.6)),\exp(-\exp(-7.9))\right\},
\]

define

\[
 \beta_{a,b}:=
 \min_{z\in\mathcal Z}
 \min\{1,z^{-D}\Phi_{a,b}(z)\}.
 \tag{4}
\]

Each value in (4) is a valid Chernoff upper bound. Different weight pairs
may use different tested tilts.

## Occupation-two bound

Equations (1)--(4) give the rank-one upper bound

\[
 U_1=
 \binom L2\frac MP
 \sum_{w=1}^B C_w\beta_{w,w}.
 \tag{5}
\]

They give the rank-two upper bound

\[
 U_2=
 \binom L2\frac{M(M-1)}{P(P-1)}
 \sum_{a=1}^B\sum_{b=1}^B
 \left(C_aC_b-\mathbf 1_{a=b}C_a\right)\beta_{a,b}.
 \tag{6}
\]

The binary64 diagnostic reports

\[
 \log_2 U_1=-64.15590473223627,
\]

\[
 \log_2 U_2=-99.70455835150280,
\]

and

\[
 \log_2(U_1+U_2)=-64.15590473220756.
\]

Thus occupation two has 64.1559 bits of diagnostic margin. Rank one is the
bottleneck. Its largest shell contribution occurs at \(w=6\). The largest
rank-two shell contribution occurs at \((a,b)=(25,25)\).

Combining this result with the occupation-one diagnostic gives

\[
 -\log_2(U_{Q=1}+U_{Q=2})=46.59280766329173.
\]

Occupation one remains the bottleneck.

## Checks and status

The implementation verifies that the rank-two ambient shell multiplicities
sum to

\[
 (2^B-1)(2^B-2).
\]

A separate exhaustive test at \(B=4\) compared recurrence (3) with direct
support-pair enumeration. The largest absolute matrix-entry difference was
\(4.45\times10^{-16}\).

`evaluate_single_random_constituent_q2.py` implements the calculation. The
receipt is `single_random_constituent_B256_q2_s22_d109.json`.

The result has three limitations.

- Nearest binary64 arithmetic is not an outward certificate.
- The result covers only occupation two.
- The expectation averages over the sampled constituent and does not select
  a fixed constituent.

The next shared-code obligation begins at occupation three. It requires
message-rank classes one, two, and three.
