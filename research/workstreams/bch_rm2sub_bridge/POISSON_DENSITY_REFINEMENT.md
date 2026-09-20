# A smaller density factor for the next frontier step

The frozen K24 and current K26 producers pay a factor (L+1)^256 when
comparing shuffled regions to independent Bernoulli inputs. The comparison
admits a smaller factor. This lemma is not used by those producers; a
new checker and new receipts are required before it affects a certificate.

Let S be a sum of n independent Bernoulli variables, with possibly
different success probabilities and mean nr. Define

\[
D_n:=\min\{n+1,4\lceil\sqrt n\rceil\}.
\]

For every j=0,...,n,

\[
\Pr[S=j]\le D_n{n\choose j}r^j(1-r)^{n-j}.
\]

For 0<r<1, the existing Chernoff argument bounds Pr[S=j] by
exp(-n*kl(j/n || r)), where kl is binary relative entropy in natural units.
The binomial probability on the right equals that exponential times

\[
b_{n,j}:={n\choose j}(j/n)^j(1-j/n)^{n-j}.
\]

At j=0 or j=n this binomial mode mass is1.

It remains to show b_(n,j)>=1/D_n. Let T have distribution Binomial(n,j/n).
Its mean is j, its variance is at most n/4, and j is a mode. The mode
property follows directly from the ratio of consecutive binomial masses.
Chebyshev's inequality gives

\[
\Pr[|T-j|<\sqrt n]\ge3/4.
\]

This interval contains at most2*sqrt(n)+1 integer values. Each has mass
at most b_(n,j). Therefore

\[
b_{n,j}\ge\frac{3}{4(2\sqrt n+1)}\ge\frac{1}{4\sqrt n}.
\]

The last inequality uses n>=1. The earlier mode bound b_(n,j)>=1/(n+1)
also applies, proving the claimed integer factor. The cases r=0 or r=1
are deterministic and satisfy the same inequality.

After a uniform shuffle, vectors of a fixed weight have equal probability.
Dividing both count probabilities by binomial(n,j) gives the same
pointwise density comparison for whole input vectors. Thus the dense
argument may pay D_L^256 instead of (L+1)^256.

For K26, L=2^19 and D_L=2900. For K28, L=2^21 and D_L=5796. This improves
the dense bound, not the Q1 bound that currently determines the full
margin. It does not add the saved logarithmic density cost to the full
distance margin.

`poisson_density_factor.py` computes D_n with integer square roots.
Tests check binomial mode masses through n=64 and small exact
Poisson-binomial distributions. The general claim rests on the proof
above, not enumeration. Next: implement the factor in a new dense checker,
compare against the frozen baseline, and replay any resulting certificate.
