# Goal 02: compressed all-support inner envelope

## Objective

Goal 01 bounds the inner failure probability at selected input supports.  The
outer calculation can activate any number of packets.  Goal 02 must therefore
cover every integer support

\[
18\le h\le N,
\qquad N=524352,
\]

without storing or optimizing one certificate for each value of \(h\).

Fix a nonzero value for each active input packet.  A uniform global packet
permutation makes their positions a uniform \(h\)-subset of \([N]\).  The
probability below is over that permutation and the independent setup matrices
of RandomStepConv.

## Shared-radius bound

For \(0<z,r<1\), define

\[
b(z):=2^{-g}(1+z)^g,
\qquad
d(z):=(1-2^{-\sigma})b(z),
\]

and

\[
G_z(r):=\frac{b(z)-d(z)r}{(1-r)(1-d(z)r)}.
\]

Goal 01 proves

\[
[x^h]F_N(x,z)
=
[r^{N-h}]\frac{G_z(r)^h}{1-r}.
\]

Nonnegative coefficients give the bound

\[
B_{z,r,D}(h):=
\frac{z^{-D}}{\binom Nh}
\frac{r^{-(N-h)}G_z(r)^h}{1-r}.
\]

Therefore,

\[
p_h(D)\le B_{z,r,D}(h).
\]

The inequality holds for every admissible pair \((z,r)\).  Numerical
optimization may select the pair, but it does not affect validity.

## Interval compression

For fixed \((z,r,D)\), write

\[
\log B_{z,r,D}(h)
=C+h\log\!\bigl(rG_z(r)\bigr)-\log\binom Nh,
\]

where \(C\) does not depend on \(h\).  The binomial coefficients form a
log-concave sequence.  Hence \(-\log\binom Nh\), and therefore
\(\log B_{z,r,D}(h)\), is convex in the integer variable \(h\).

Fix an integer interval \([a,b]\).  One pair \((z,r)\) certifies

\[
\max_{a\le h\le b}p_h(D)
\le
\max\{B_{z,r,D}(a),B_{z,r,D}(b)\}.
\]

The receipt stores one pair for each interval.  It also stores the two
outward-rounded endpoint bounds.  These data define a functional bound
\(B_{z,r,D}(h)\) for every support in the interval, not only a constant cap.

## Acceptance criteria

Goal 02 is complete when all of the following conditions hold:

1. The proof derives the shared-radius bound from the exact Goal 01 identity.
2. The proof establishes endpoint sufficiency from binomial log-concavity.
3. A checker confirms discrete convexity and compares interval certificates
   with independently optimized pointwise bounds.
4. One receipt covers every integer support from 18 through \(N\).
5. The receipt uses outward-rounded arithmetic for every interval endpoint.

This goal remains an inner-only statement.  It does not sum an outer spectrum
or establish an end-to-end minimum distance.
