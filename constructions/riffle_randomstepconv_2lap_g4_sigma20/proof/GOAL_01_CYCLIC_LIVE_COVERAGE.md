# Goal 01: exact retained-lap gap bound

## Objective

For **Riffle RandomStepConv-2Lap g=4 sigma=20**, derive and implement an
upper bound on

\[
p_h(D):=
\Pr\!\left[\operatorname{wt}(Y)\le D
\;\middle|\;
\operatorname{supp}(X)\text{ is a uniform }h\text{-subset of }[N]
\right].
\]

The probability is over the global packet permutation and both random-map
families. Fix the nonzero values of the \(h\) active packets arbitrarily.
Random-map symmetry makes the probability independent of those values.

This goal concerns only the inner encoder.

## Gap variables

Fix \(h\ge1\) active positions. Let \(p\) and \(t\) be the inactive prefix
and suffix lengths. For \(1\le j<h\), let \(\ell_j\) be the inactive gap
between consecutive active positions. Then

\[
p+t+\sum_{j=1}^{h-1}\ell_j=N-h.
\]

Every \(h\)-subset corresponds to one such weak composition. Hence averaging
over a uniform support divides the composition sum by \(\binom Nh\).

Set

\[
q:=2^{-\sigma},
\qquad
a:=1-q.
\]

Let \(U_\ell\) count off packets in a zero-input run of length \(\ell\) that
follows an active packet. Its probability-generating function is

\[
A_\ell(u):=\mathbb E[u^{U_\ell}]
=a^\ell+\sum_{k=0}^{\ell-1}a^kq\,u^{\ell-k}.
\]

The corresponding length-generating function is

\[
A(r,u):=\sum_{\ell\ge0}A_\ell(u)r^\ell
=\frac{1-aur}{(1-ar)(1-ur)}.
\]

## The gap across the wrap

The retained prefix does not begin from an arbitrary state. The last active
packet of the burn-in lap refreshes the state before the burn-in suffix of
length \(t\).

Suppose a zero run of length \(p\) begins in a nonzero state. Let
\(B_p(u)\) be its off-count generating function. Then

\[
B_0(u)=1
\]

and, for \(p\ge1\),

\[
B_p(u)=a^{p-1}+
\sum_{j=1}^{p-1}a^{j-1}q\,u^{p-j}.
\]

Define

\[
B(r,u):=\sum_{p\ge0}B_p(u)r^p
=1+\frac{r}{1-ar}
+\frac{qur^2}{(1-ar)(1-ur)}.
\]

The probability that the state survives the last burn-in activation and the
following \(t\) burn-in zero steps is \(a^{t+1}\). Therefore the retained
prefix has generating function

\[
V_{t,p}(u):=
a^{t+1}B_p(u)+(1-a^{t+1})u^p.
\]

The retained suffix begins at the corresponding active packet in the second
lap. Its generating function is \(A_t(u)\), independently of the prefix.
Consequently, define

\[
H(r,u):=
\sum_{p,t\ge0}A_t(u)V_{t,p}(u)r^{p+t}.
\]

The two sums give the rational identity

\[
H(r,u)=
aA(ar,u)B(r,u)
+\bigl(A(r,u)-aA(ar,u)\bigr)\frac1{1-ur}.
\]

All coefficients of \(H\) are nonnegative, even though the compact rational
form contains one subtraction.

## Exact retained-weight moment

For \(0<z\le1\), define the live-packet weight moment

\[
b(z):=2^{-g}(1+z)^g,
\qquad
u(z):=b(z)^{-1}.
\]

Conditioned on the live/off pattern, every live output packet is independent
and uniform in \(\mathbb F_2^g\). If \(U\) retained packets are off, then

\[
\mathbb E[z^{\operatorname{wt}(Y)}\mid U]=b(z)^{N-U}.
\]

Combining the gap factors gives the exact identity

\[
\mathbb E[z^{\operatorname{wt}(Y)}\mid |\operatorname{supp}(X)|=h]
=
\frac{b(z)^N}{\binom Nh}
[r^{N-h}]H(r,u(z))A(r,u(z))^{h-1}.
\tag{1}
\]

This identity is the main structural target of Goal 01. It replaces the
one-lap late-suffix decomposition with an off-count calculation on the gaps.

## Proof-facing coefficient bound

All coefficients in (1) are nonnegative. For

\[
0<r<b(z),
\]

coefficient tilting and Markov's inequality give

\[
p_h(D)
\le
\frac{z^{-D}b(z)^N}{\binom Nh}
\frac{H(r,u(z))A(r,u(z))^{h-1}}{r^{N-h}}.
\tag{2}
\]

The evaluator will minimize (2) over \(z\in(0,1)\) and \(r\in(0,b(z))\).

## Acceptance criteria

Goal 01 is complete when all of the following conditions hold:

1. A proof derives the segment laws and identities (1) and (2).
2. An exhaustive checker matches (1) on small instances.
3. A target evaluator remains stable for \(N=524352\).
4. The evaluator covers the one-lap 9% saddle range
   \(16384\le h\le32767\).
5. A receipt separates exact identities, floating diagnostics, and certified
   numerical bounds.

The first implementation may use ordinary floating-point optimization. An
outward-rounded evaluator is a later acceptance requirement.
