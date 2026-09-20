# Linear-time outer certificate

## Candidate and purpose

The bulk transfer in `BULK_OCCUPATION_TRANSFER.md` needs a local outer code
with a pointwise spectrum envelope. A random dense local code has the desired
spectrum, but its ordinary and transposed encoders cost \(\Theta(B^2)\).

This note replaces that reference by a BA-3 ensemble. Let \(B=24m\). Start
with the direct sum of \(m\) extended binary Golay \([24,12,8]\) codes. Apply
an independent uniform interleaver, a length-\(B\) accumulator, a second
independent uniform interleaver, and a second accumulator.

The resulting local code has rate one half. Each constituent block has
constant size. Each interleaver and accumulator costs \(O(B)\) in either
circuit direction. Thus one local encoder and its transpose both cost
\(O(B)\).

The certified envelope is

\[
  \mathcal W=[0.104,0.896],
  \qquad
  d_{\rm out}(1/2)\le 0.36.
\tag{1}
\]

All logarithms in (1) and below are natural. The interval and endpoint checks
below prove these constants; `OUTER_INTERVAL_RECEIPT.json` records the
outward-rounded results and verifier hashes.

## Exact expected spectrum

The Golay weight polynomial is

\[
  G(x)=1+759x^8+2576x^{12}+759x^{16}+x^{24}.
\tag{2}
\]

Let

\[
  G_{B}(u):=[x^u]G(x)^{B/24}.
\tag{3}
\]

For \(1\le u\le B\), set \(r:=\lceil u/2\rceil\). The exact
input--output weight enumerator of the accumulator is

\[
  T_B(u,h)
  :=
  \binom{h-1}{r-1}
  \binom{B-h}{u-r}.
\tag{4}
\]

Set \(T_B(0,0):=1\) and \(T_B(0,h):=0\) for \(h\ne0\). Under a uniform
interleaver, the transition probability from input weight \(u\) to output
weight \(h\) is

\[
  P_B(u,h):=\frac{T_B(u,h)}{\binom Bu}.
\tag{5}
\]

Let \(\overline A_B(w)\) be the expected number of BA-3 codewords of weight
\(w\), where the expectation is over the two interleavers. Equations
(2)--(5) give the exact finite formula

\[
  \boxed{
  \overline A_B(w)
  =
  \sum_{u=0}^B\sum_{h=0}^B
    G_B(u)P_B(u,h)P_B(h,w).
  }
\tag{6}
\]

The term \(u=0\) contributes only to \(w=0\). Thus (6) is also the
nonzero spectrum formula for every \(w>0\).

## Deterministic selection

The complete SPIN construction may reuse one selected local code at every
outer-block position. It need not sample a new BA-3 code at each position.
This avoids multiplying the outer-selection failure probability by \(n/B\).

Fix a tail set \(S_B\subseteq\{1,\ldots,B\}\). Suppose

\[
  \sum_{w\in S_B}\overline A_B(w)=o(1).
\tag{7}
\]

Markov's inequality shows that a sampled BA-3 code has no nonzero word in
\(S_B\) with probability \(1-o(1)\). A second application of Markov and a
union bound show that

\[
  A_B(w)\le B^2\overline A_B(w)
  \quad\text{for every }w\notin S_B
\tag{8}
\]

with probability at least \(1-1/B\). Therefore, for all sufficiently large
\(B\), there is one deterministic BA-3 realization that satisfies (7) and
(8) simultaneously.

Take

\[
  S_B:=
  \{w:1\le w<0.104B\}
  \cup
  \{w:0.896B<w\le B\}.
\tag{9}
\]

It remains to verify the following two limits from (6):

\[
  \sum_{w\in S_B}\overline A_B(w)=o(1),
\tag{10}
\]

\[
  \max_{0.104B\le w\le0.896B}
  \frac1B
  \ln\frac{\overline A_B(w)2^B}{\binom Bw}
  <0.36-o(1).
\tag{11}
\]

For the selected realization, the factor \(B^2\) in (8) contributes
\(o(B)\) to the logarithm. Equations (9)--(11) then prove (1).

## Constant-size asymptotic gate

The exact finite formula has a three-variable asymptotic relaxation. Define

\[
  g(\alpha)
  :=\inf_{r>0}
  \left\{
    \frac1{24}\ln G(r)-\alpha\ln r
  \right\}.
\tag{12}
\]

For feasible \((\alpha,\beta)\), define

\[
  p(\alpha,\beta)
  :=
  \beta h\!\left(\frac{\alpha}{2\beta}\right)
  +(1-\beta)
    h\!\left(\frac{\alpha}{2(1-\beta)}\right)
  -h(\alpha).
\tag{13}
\]

The feasibility range is

\[
  \frac\alpha2\le\beta\le1-\frac\alpha2.
\]

Set \(p(\alpha,\beta):=-\infty\) outside this range. Stirling's
inequality applied to (4)--(6) gives

\[
  \limsup_{B\to\infty}
  \frac1B\ln\overline A_B(\lfloor\omega B\rfloor)
  \le
  a_{\rm BA}(\omega),
\tag{14}
\]

where

\[
  a_{\rm BA}(\omega)
  :=
  \sup_{0\le\alpha,\beta\le1}
  \left[
    g(\alpha)+p(\alpha,\beta)+p(\beta,\omega)
  \right].
\tag{15}
\]

Polynomial factors from the two sums in (6) disappear in (14). The
variational exponent can approach zero at an endpoint when a constant-weight
path has only polynomial probability. Thus a compact exponent check must be
paired with a boundary estimate.

Fix a constant \(0<\eta<0.104\). The required boundary estimates are

\[
  \sum_{1\le w<\eta B}\overline A_B(w)=o(1),
  \qquad
  \sum_{(1-\eta)B<w\le B}\overline A_B(w)=o(1).
\tag{16}
\]

The BA-3 small-output argument supplies the first estimate because the Golay
constituent distance is eight. The second estimate requires the analogous
large-output calculation from the same exact accumulator enumerator. After
choosing \(\eta\), the compact tail checks are

\[
  \sup_{\eta\le\omega\le0.104}a_{\rm BA}(\omega)<0,
  \qquad
  \sup_{0.896\le\omega\le1-\eta}a_{\rm BA}(\omega)<0.
\tag{16a}
\]

Together, (16) and (16a) prove (10). The central spectrum check is

\[
  \sup_{0.104\le\omega\le0.896}
  \left[
    a_{\rm BA}(\omega)-h(\omega)+\ln2
  \right]
  <0.36.
\tag{17}
\]

The large-output endpoint reduces to the small-output accumulator tail. Let
\(e_1=(1,0,\ldots,0)\). For every input \(v\), linearity gives

\[
  \operatorname{Acc}_B(v)+\mathbf1
  =\operatorname{Acc}_B(v+e_1).
\tag{17a}
\]

If \(V\) is uniform on the weight-\(h\) slice, then \(V+e_1\) is a mixture
of a weight-\(h-1\) slice conditioned to have first bit zero and a
weight-\(h+1\) slice conditioned to have first bit one. Therefore

\[
\begin{split}
 &\Pr[\operatorname{wt}(\operatorname{Acc}_B(V))\ge(1-\eta)B]\\
 &\quad\le
 \frac{h}{B-h+1}\,p_{h-1}(\eta)
 +\frac{B-h}{h+1}\,p_{h+1}(\eta),
\end{split}
\tag{17b}
\]

where \(p_j(\eta)\) is the low-output probability for a uniform weight-
\(j\) input. The two prefactors in (17b) are at most \(B\). Thus they do not
change an exponential estimate. Applying the BA-3 small-output summation
with the last input weight shifted by one proves the second estimate in (16)
whenever the first estimate has a strict polynomial margin.

## Endpoint proof at relative weight 0.104

This section proves the first estimate in (16) for

\[
  \delta_{\rm out}:=0.104=\frac{13}{125}.
\]

The large-output estimate then follows from (17a)--(17b).

Let \(b\) be the number of active Golay blocks. At most

\[
  N_b
  :=\binom{B/24}{b}(2^{12}-1)^b
  \le
  \left[
    \frac{eB(2^{12}-1)}{24b}
  \right]^b
\tag{17c}
\]

messages activate exactly \(b\) blocks. Every such message enters the first
accumulator with an even weight \(u=2r\) that satisfies

\[
  4b\le r\le12b.
\tag{17d}
\]

We first bound the low-output tail of one accumulator. Let \(V\) be uniform
on the weight-\(v\) slice, and set
\(H:=\operatorname{wt}(\operatorname{Acc}_B(V))\). The exact enumerator (4)
and the hockey-stick identity give

\[
  \Pr[H\le\delta_{\rm out}B]
  \le c^{\lfloor v/2\rfloor},
  \qquad
  c:=\frac{4\delta_{\rm out}}{1-\delta_{\rm out}}
  =\frac{13}{28}<1.
\tag{17e}
\]

For even \(v=2s\), the bound follows from

\[
\begin{split}
 \Pr[H\le\delta_{\rm out}B]
 &\le
 \frac{\binom Bs\binom{\lfloor\delta_{\rm out}B\rfloor}{s}}
      {\binom B{2s}}\\
 &=
 \binom{2s}{s}
 \frac{\binom{\lfloor\delta_{\rm out}B\rfloor}{s}}
      {\binom{B-s}{s}}
 \le
 \left(\frac{4\delta_{\rm out}}{1-\delta_{\rm out}}\right)^s.
\end{split}
\]

The odd case \(v=2s+1\) has the additional factor
\(2\delta_{\rm out}/(1-\delta_{\rm out})<1\). Thus (17e) holds for both
parities.

Set

\[
  q:=\sqrt c,
  \qquad
  K:=\frac{4q}{1-q}.
\]

Equation (17e) bounds the second accumulator by

\[
  \Pr[\text{final weight}\le\delta_{\rm out}B\mid H_1]
  \le q^{-1}q^{H_1}.
\tag{17f}
\]

For the first accumulator and its even input weight \(u=2r\), equation (4)
gives

\[
\begin{split}
 \mathbb E[q^{H_1}\mid u=2r]
 &\le
 \frac{\binom Br}{\binom B{2r}}
 \sum_{h\ge r}\binom{h-1}{r-1}q^h\\
 &=
 \frac{\binom Br}{\binom B{2r}}
 \left(\frac q{1-q}\right)^r\\
 &\le
 \left[
   \frac{Kr}{B-2r+1}
 \right]^r.
\end{split}
\tag{17g}
\]

Take \(b_0:=10^{-3}\). For \(b\le b_0B\), equation (17d) gives
\(r/B\le0.012\). The last expression in (17g) decreases with \(r\)
throughout this range. Therefore its maximum occurs at \(r=4b\). Equations
(17c), (17f), and (17g) imply

\[
  \mathbb E[Z_{\rm sparse}]
  \le
  q^{-1}
  \sum_{1\le b\le b_0B}
  \left[
    C_*\left(\frac bB\right)^3
  \right]^b,
\tag{17h}
\]

where

\[
  C_*:=
  \frac{e(2^{12}-1)}{24}
  \left[
    \frac{4K}{1-8b_0}
  \right]^4.
\]

Outward-rounded evaluation gives

\[
  C_*b_0^3<0.657<1.
\tag{17i}
\]

For \(b\le\sqrt B\), the sum in (17h) is bounded by a geometric series
whose base is \(O(B^{-3/2})\). For \(b\ge\sqrt B\), equation (17i) bounds
the summand by \(0.657^b\). Hence

\[
  \mathbb E[Z_{\rm sparse}]=o(1).
\tag{17j}
\]

It remains to treat \(b>b_0B\). Equations (17d) and the two accumulator
support constraints imply

\[
  \alpha:=u/B\ge0.008,
  \qquad
  \omega:=h/B\ge0.002.
\]

The outward-rounded branch verifier covers

\[
  0.008\le\alpha\le1,
  \quad
  0\le\beta\le1,
  \quad
  0.002\le\omega\le0.104.
\]

It proves

\[
  g(\alpha)+p(\alpha,\beta)+p(\beta,\omega)<0
\tag{17k}
\]

at every feasible point in this box. The verifier exhausted 6,747 rational
boxes. Its largest accepted upper endpoint was
\(-7.68\times10^{-8}\). Uniform Stirling bounds add only \(O(\log B/B)\)
to (17k). Thus the expected number of dense-block messages with output at
most \(0.104B\) is exponentially small.

Equations (17j) and (17k) prove

\[
  \sum_{1\le w<0.104B}\overline A_B(w)=o(1).
\tag{17l}
\]

Finally, (17a)--(17b) add at most a polynomial factor to the same bounds.
They prove

\[
  \sum_{0.896B<w\le B}\overline A_B(w)=o(1).
\tag{17m}
\]

Therefore the endpoint condition (10) holds for the tail set (9).

Equations (12)--(17) reduce the outer work to two endpoint estimates and a
fixed-dimensional interval-arithmetic task. They are independent of \(n\)
and \(B\).

## Certificate result and diagnostic cross-check

`verify_outer_interval.py` and `verify_outer_sparse_constants.py` certify
(16)--(17). The former exhausts the central and dense-tail rational boxes;
the latter proves the sparse endpoint constants. Their commands, hashes, box
counts, and largest accepted upper endpoints appear in
`OUTER_INTERVAL_RECEIPT.json`.

`analyze_golay_ba_outer.py` separately evaluates (6) in exact combinatorial
form followed by binary64 log arithmetic. The following values are only a
cross-check on the certified result.

At \(B=1920\), the diagnostic reports:

\[
  \frac1B\ln\overline A_B(0.104B)\approx-0.00315,
\]

\[
  \max_{0.104B\le w\le0.896B}
  \frac1B\ln
  \frac{\overline A_B(w)2^B}{\binom Bw}
  \approx0.3577.
\]

The expected spectrum first becomes nonnegative near relative weight
\(0.10625\). The upper crossing is near \(0.89375\). These values support
(16a)--(17). The sublinear-path branch has likelihood exponent at most

\[
  \ln2-h(0.104)\approx0.35936<0.36.
\]

This cutoff is necessary for the stated \(0.36\) envelope: the same boundary
value at relative weight \(0.102\) is about \(0.36369\). This explains why
the certified outer cutoff is \(0.104\), rather than \(0.102\).

## Interaction with the structured inner

Assume (1). Then the one-active spectrum exponent may use

\[
  a(\beta)
  =h(\beta)-\ln2+0.36,
  \qquad \beta\in[0.104,0.896].
\tag{18}
\]

For the frozen \((t,s,d_A)=(128,19,48)\) inner and certified target distance
\(\delta=0.101\), the outward-rounded one-active verifier gives

\[
  \chi<-1.30,
  \qquad
  \sup_{\beta\in[0.104,0.896]}
  [a(\beta)+J(\beta)+\beta\chi]
  <-0.0606.
\tag{19}
\]

Thus the one-active class has substantial slack. A growth law

\[
  B=c\ln n+o(\ln n),
  \qquad c>1/0.0606,
\]

is sufficient after outward rounding. Taking \(B/\ln n\to\infty\) avoids
optimizing this coefficient and remains compatible with \(B=o(n)\).

For the bulk gate, (1) replaces the random-spectrum term
\(\tau(\ln2)/2\) by \(0.36\tau\). The conservative shell transfer is
certified for every \(0<\tau\le0.985\) at \(\delta=0.101\). The direct
high-density bound begins before

\[
  \tau>
  \frac{h(0.101)}{\ln2-0.36}
  <0.982373.
\tag{20}
\]

The certified bands overlap. Together with the fixed- and
vanishing-occupation bounds, they give the complete distance certificate in
`ASYMPTOTIC_DISTANCE_CERTIFICATE.md`.
