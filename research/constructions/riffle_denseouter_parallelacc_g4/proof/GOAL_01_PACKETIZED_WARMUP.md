# Goal 01 proof: packetized dense-outer warmup

## Result

Riffle DenseOuter-ParallelAcc g=4 has linear minimum distance with high
probability. A coarse one-lane projection proves the claim without analyzing
packet values or prefix-XOR orbits.

The proved relative-distance range is

\[
0<\delta<\delta_4,
\qquad
\delta_4:=\frac{(\sqrt2-1)^8}{16e}
=1.992416147\ldots\times10^{-5}.
\tag{1}
\]

The constant in (1) is a proof floor, not a prediction of the true minimum
distance.

## Fixed-word contraction

Let \(u\in\mathbb F_2^{gN}\) be a fixed nonzero binary word, grouped into
\(N\) packets of width \(g\). Let \(w:=\operatorname{wt}(u)\). A uniform packet
permutation applies the same permutation to all \(g\) bit lanes.

For lane \(j\in[g]\), let \(r_j\) be its input weight. Then

\[
\sum_{j=1}^g r_j=w.
\]

Fix a lane \(j_*\) with maximum input weight and write \(r:=r_{j_*}\). Thus

\[
r\ge \frac wg.
\tag{2}
\]

The uniform packet permutation maps the active positions of lane \(j_*\) to a
uniform \(r\)-subset of \([N]\). Its marginal distribution is therefore the
same Hamming-slice distribution used in the scalar accumulator theorem.

Let \(W\) be the total binary weight produced by all \(g\) accumulator lanes.
Fix an integer threshold \(D\ge1\), and define

\[
\theta:=\frac{4eD}{N}.
\]

Assume \(\theta<1\). If \(r>2D\), the selected scalar lane has output weight at
least \(\lceil r/2\rceil>D\). Hence \(\Pr[W\le D]=0\).

Suppose instead that \(r\le2D\). The inequality \(\theta<1\) gives
\(D/N<1/(4e)<1/4\) and \(r<N/2\). The scalar accumulator contraction bound
then gives

\[
\begin{aligned}
\Pr_\Pi[W\le D]
&\le
\Pr_\Pi[W_{j_*}\le D]\\
&\le
\left(\frac{4eD}{N}\right)^{r/2}\\
&\le
\left(\frac{4eD}{N}\right)^{w/(2g)}.
\end{aligned}
\tag{3}
\]

The last inequality uses (2) and \(\theta<1\). Consequently, with

\[
z_D:=\left(\frac{4eD}{N}\right)^{1/(2g)},
\tag{4}
\]

every fixed word satisfies

\[
\Pr_\Pi[W\le D]\le z_D^w.
\tag{5}
\]

Equation (5) is uniform over the packet layout and packet values. Dependence
between accumulator lanes is irrelevant because the proof uses only one lane.

## Dense outer generating function

The span lemma in `outerDense.tex` fails for messages with a zero gap longer
than the outer memory. The following cluster argument supplies the required
generating-function bound.

Fix a nonzero message with one-positions

\[
1\le p_1<\cdots<p_w\le k.
\]

Its active parity positions form the union of \([p_i,p_i+M]\). Hence their
number is

\[
s(x)=M+1+\sum_{i=2}^w\min\{M+1,p_i-p_{i-1}\}.
\tag{6}
\]

At every active position, the parity bit is either uniform or deterministically
one. The tap vectors are independent across positions. Define

\[
a(z):=\frac{1+z}{2},
\]

For \(0<z<1\), the fixed-message contribution satisfies

\[
\mathbb E\left[z^{\operatorname{wt}(C_{\mathrm{out}}(x))}\right]
\le z^w a(z)^{s(x)}.
\tag{7}
\]

Partition the one-positions into clusters. Two consecutive ones belong to the
same cluster when their gap is at most \(M+1\). For one cluster, summing over
all possible internal gaps gives at most

\[
B_M(z):=
\frac{z a(z)^{M+1}}{1-\kappa_M(z)},
\qquad
\kappa_M(z):=z\sum_{d=1}^{M+1}a(z)^d.
\tag{8}
\]

The first one contributes \(z a(z)^{M+1}\). Every later one at gap
\(d\le M+1\) contributes \(z a(z)^d\).

If \(0<z<\sqrt2-1\), then

\[
\kappa_M(z)
\le \frac{z a(z)}{1-a(z)}
=\frac{z(1+z)}{1-z}
<1.
\tag{9}
\]

A support with \(r\) clusters has at most \(\binom{k}{r}\) choices for its
ordered cluster starts. Ignoring boundary and separation constraints only adds
invalid supports. Therefore

\[
W_{\mathrm{out}}(z)
\le
(1+B_M(z))^k-1
\le
\exp(kB_M(z))-1.
\tag{10}
\]

If \(kB_M(z)\le1\), equation (10) also gives

\[
W_{\mathrm{out}}(z)\le(e-1)kB_M(z).
\tag{11}
\]

## Finite first-moment theorem

Sample the dense outer taps and the packet permutation as specified in
`CONSTRUCTION.md`. If \(z_D<\sqrt2-1\), then

\[
\Pr[d_{\min}\le D]
\le
\exp(kB_M(z_D))-1.
\tag{12}
\]

Indeed, let \(Z_D\) count nonzero messages whose final codeword has weight at
most \(D\). Condition on the dense outer word for each message. Equation (5)
then yields

\[
\mathbb E[Z_D]
\le
\sum_{x\ne0}
\mathbb E_{C_{\mathrm{out}}}
\left[z_D^{\operatorname{wt}(C_{\mathrm{out}}(x))}\right]
=W_{\mathrm{out}}(z_D).
\]

Equation (10) bounds this expectation. Markov's inequality proves (12).
The zero padding used to complete the last packet does not change the outer
weight and therefore does not change the argument.

## Asymptotic theorem

Fix a packet width \(g\ge1\) and a constant

\[
0<\delta<\delta_g,
\qquad
\delta_g:=\frac{(\sqrt2-1)^{2g}}{4eg}.
\tag{13}
\]

Let \(D:=\lfloor\delta gN\rfloor\), and define

\[
z_\delta:=(4eg\delta)^{1/(2g)}.
\]

Then \(z_D\le z_\delta<\sqrt2-1\). Choose

\[
M:=\lceil c\log_2 k\rceil,
\qquad
c>\frac{1}{\log_2(2/(1+z_\delta))}.
\tag{14}
\]

Equations (8)--(10) and \(a(z_\delta)<1\) give

\[
\Pr[d_{\min}\le\delta gN]
\le
\exp(kB_M(z_\delta))-1
=o(1).
\tag{15}
\]

For \(g=4\), (13) gives exactly (1). Therefore the four-bit packet construction
has a positive relative distance with high probability.

## Finite scale at the current length

At binary length \(L=4N=2{,}097{,}408\) and threshold \(D=40\), equation (4)
gives

\[
z_D=0.4119542066\ldots<\sqrt2-1.
\]

Thus the coarse theorem reaches the 40-bit threshold. The explicit bound (12)
depends strongly on the outer memory:

| Outer memory \(M\) | Upper bound from (12) | Failure bits |
|---:|---:|---:|
| 56 | \(1.00\times10^{-1}\) | 3.32 |
| 64 | \(5.92\times10^{-3}\) | 7.40 |
| 80 | \(2.25\times10^{-5}\) | 15.44 |
| 90 | \(6.91\times10^{-7}\) | 20.46 |
| 128 | \(1.24\times10^{-12}\) | 39.55 |
| 129 | \(8.76\times10^{-13}\) | 40.05 |

Here “failure bits” means \(-\log_2\) of the displayed first-moment bound.
For each row, choose \(k=\lfloor(L-M)/2\rfloor\) and append at most one zero
bit before packetization.

These are rigorous ensemble bounds for the corrected terminated outer. They do
not prove that a sampled fixed matrix has the stated distance until that matrix
is separately certified or the ensemble failure probability is accepted.

## What remains open

Goal 01 removes the packet-value obstruction for the asymptotic warmup. It does
not analyze the original BCH outer, a structured riffle permutation, or any
two-lap constituent.

The one-lane projection is deliberately coarse. A joint-lane bound could raise
the relative-distance constant or reduce the outer memory required for a fixed
failure budget.
