# Matched whole-class mean gate for two inner streams

This diagnostic asks whether interleaving two inner states improves the
existing homogeneous arity-eight/F1 bound. It first uses one common
potential, not an enumeration of four-state paths. A failure of this
comparison would not establish a low-distance codeword.

Fix a source set $G$ of size $r=18$ whose larger full-pair count in the two
source matchings is exactly one. Fix every active group's eight-row set
through all sixteen outer columns. The outer map and paired-and-swapped
route are unchanged. Only the proposed inner schedule differs, as defined
in [the two-stream proposal](HD_TWO_STREAM_INNER_PROPOSAL.md).

The mark is fixed at $z=691/700$. All local and regional operators below
are unpowered first-moment caps. This gate includes every outer syndrome
and every regional boundary state. It does not cover other source strata,
cohorts, or weights, and it is not a directed distance certificate.

## Joint regions and their exact quotient

For $d\in\{1,2\}$, let $R^{(d)}_{m,f}$ be the regional transfer for a
fixed input support of size $m$ containing $f$ full pairs. It averages
the route, fair swaps, nonzero amplitudes, and refresh keys. The one-stream
states are zero and live. The two-stream states are $00,01,10,11$, with
the first bit referring to even packets.

The provider uses one packet polynomial $B$ and one half-region power
$H=B^{64}$. The old regional control extracts its total coefficient from
$HH$; the proposed region extracts it from $H\otimes H$. Both use the
same total histogram and one 4096-slot normalizer. Neither half is
normalized independently. The two-stream product cone carries the
product of the normalized live-state references, not an assertion that
the two routed occupancies are independent.

Stream exchange commutes with each joint coefficient. Therefore the
partition $\{00\},\{01,10\},\{11\}$ has an exact quotient. For a joint
matrix $R$ indexed by $00,01,10,11$, it is

$$
\overline R=
\begin{pmatrix}
R_{00}&2R_{01}&R_{03}\\
R_{10}&R_{11}+R_{12}&R_{13}\\
R_{30}&2R_{31}&R_{33}
\end{pmatrix}.
$$

The middle reference is the normalized equal mixture of the two one-live
references. The quotient preserves every chronological product starting
in $00$ and ending against the all-ones vector. It does not reset a state
at a region or column boundary.

## One scalar comparison, matched across constructions

Use the same fixed number $p=2^{-64}$ for both constructions:

$$
v^{(1)}=(1,p)^T,\qquad
v^{(2)}=(1,p,p,p^2)^T.
$$

The quotient potential is $(1,p,p^2)^T$. Define

$$
c_d(m)=\max_{0\le f\le\min(1,\lfloor m/2\rfloor)}
\max_i\frac{(R^{(d)}_{m,f}v^{(d)})_i}{v^{(d)}_i}.
$$

The maximum over $f$ is taken after each weighted row sum. This avoids
combining unrelated entrywise maxima. Every actual regional subset of a
source in the declared stratum has an allowed $f$, so the cap is uniform
over its fixed source geometry. Empty input fixes the zero state and
contracts live mass; the implementation checks $c_d(0)\le1$ explicitly.

Let $e_0$ select the all-zero initial state, and let $\mathbf1$ have every
entry equal to one. For any chronological sequence of 256 regional inputs,

$$
e_0^T\left(\prod_iR^{(d)}_{m_i,f_i}\right)\mathbf1
\le p^{-d}\prod_i c_d(m_i).
$$

This follows by repeated nonnegative matrix action on $v^{(d)}$ and
$\mathbf1\le p^{-d}v^{(d)}$. The endpoint factor is charged once: 64 bits
for one stream and 128 bits for two streams. No matrix commutativity,
stationary distribution, or independence between regional subsets is
used in this inequality.

For each construction, take only the secant of $\log_2 c_d(m)$ at $m=2$
and $m=18$. Call its slope $b_d$, and set

$$
\log_2\eta_d=\max_{1\le m\le18}
  (\log_2c_d(m)-mb_d),\quad
B_d=2^{b_d},\quad A_d=\max(1,\eta_d),\quad\theta_d=\eta_d/A_d.
$$

Then every $c_d(m)$ is bounded by the $m$th moment of the finite measure

$$
\mu_d=A_d[(1-\theta_d)\delta_0+\theta_d\delta_{B_d}],
\qquad 0^0:=1.
$$

When $A_d>1$, the zero moment is deliberately relaxed to $A_d$ and the
calculation charges $A_d^{256}$. It does not force an invalid probability
law or try another slope. All occupancies enter the intercept maximum.

## Fixed row sets and the complete word tail

For a labelled regional output mask $T\subseteq[16]$, let $N_A(T)$ count
uniform inputs supported on the fixed row set $A$ that produce precisely
$T$ under the affine two-sum outer map. The closed-kernel count gives

$$
W(w)=\left(\sum_{|A|=8}N_A(T)^{16}\right)^{1/16},\qquad w=|T|.
$$

This is a raw weight: it excludes the RS denominator. Hölder over the
sixteen columns bounds the sum over a row set held fixed across those
columns. Applying it separately to each of the eighteen groups supplies
all row-set choices once. There is no additional
$\binom{16}{8}^{18}$ factor.

For half arity, a quiet column has input and output weight at most eight;
every nonquiet column has output weight at least nine and input weight at
most eight. Thus $\kappa(w)=\min(w,8)$ is an upper surrogate for input
weight, independent of $A$. Every actual word in the class has aggregate
input weight at least $9\cdot8\cdot18=1296$.

Conditional on $j$ nonzero atoms among the sixteen regions, the scalar
group polynomial is

$$
P_{d,j}(X)=\sum_{w=0}^{j}\binom jw W(w)B_d^wX^{\kappa(w)}.
$$

The complete positive column polynomial is

$$
Q_d(X)=\sum_{j=0}^{16}\binom{16}j
\theta_d^j(1-\theta_d)^{16-j}P_{d,j}(X)^{18}.
$$

The binomial coefficient counts labelled atom locations; $\binom jw$
counts labelled output masks. These auxiliary atoms implement a positive
moment expansion. They do not resample the source, row sets, or route.

Let $N_1$ be the exact number of source sets with $F_{\max}=1$. The
requested class first moment is at most

$$
\frac{N_1z^{-734003}}{255\,255^{1152}}\,
p^{-d}A_d^{256}\sum_{k\ge1296}[X^k]Q_d(X)^{16}.
$$

The source count, projective factor, RS normalization, moment mass, and
terminal factor each occur once. This is a direct first moment; no
fractional power is applied after a partial mask sum.

## Test status and limitation

`theory/hd_two_stream_mean_gate.py` implements this matched diagnostic.
Its frozen source SHA256 is
`c77a0b5a7fb5ec50ea49baf7f09ebdb4731a4065d504c70dbcee5f5d6a2178f8`.
Seven small tests pass in 0.014 seconds. They check noncommuting word domination, exact
stream lumping, the order of the full-pair maximum, finite-measure mass,
labelled-mask contraction, and fail-closed guards. Before comparing bounds,
the driver also requires all four one-stream F1 curves to reproduce the
saved control. No regional bank or
quantitative gate has been run by this helper's author.

The scalar potential may erase the proposed improvement by maximizing
over incoming states. In particular, routing all active pairs to one
stream is rare, but that stream has only half as many intervening
packets. Those effects can offset in an averaged return bound. The first
comparison must report both constructions under the same rule before
attempting a more expensive matrix-valued moment envelope.

## Bounded result

Root ran the matched control once in 1.206 seconds. The fresh one-stream
regional entries reproduce the saved four-curve F1 bank exactly in FLOAT.
The resulting log2 bounds are $+4462.310394$ for one stream and
$+8109.759132$ for two streams. Both fail; the proposed interleaving is
3647.45 log2 units worse under this common-potential comparison. Fourteen
combined provider and gate tests pass.

The result is saved in `theory/hd_two_stream_mean_gate_float.json`, SHA256
`fc3dd65eb1288922903d9e8957974f10a8d60ea7644d93166fa37f5acc71781e`.
It does not prove that interleaving reduces actual distance. It also does
not assess a coupled two-state recurrence, which would require new local
reference measures and transition caps rather than the tensor provider.
