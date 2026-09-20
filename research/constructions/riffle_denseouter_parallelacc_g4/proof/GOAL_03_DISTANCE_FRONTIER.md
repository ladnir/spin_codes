# Goal 03 proof: finite distance frontier

## Result

Let (N=524{,}352) be the number of four-bit packets, so the binary length is
(L=4N=2{,}097{,}408). The hybrid argument below proves, for every fixed outer
word of bit weight (w),

\[
\Pr_\Pi[W\le 204]\le \rho^{w},
\qquad
\rho=0.4139892039907315<\sqrt 2-1.
\tag{1}
\]

Composing (1) with the corrected dense-outer generating function gives the
following finite ensemble guarantees.

| Outer memory | Upper bound on failure | Failure bits |
|---:|---:|---:|
| 97 | (2^{-20.4436}) | 20.44 |
| 137 | (2^{-40.4528}) | 40.45 |

Thus memory 97 certifies (d_{\min}>204) with more than 20 failure bits, and
memory 137 certifies the same statement with more than 40 failure bits.

Distance 204 is a proved lower frontier for these parameters. It is not an
estimate of the true minimum distance.

## Baseline moment frontier

Goal 02 applies one optimized gap-moment bound to every support size. A full
sweep at distance 124 gives

\[
\rho_{124}=0.41379133700004356<\sqrt2-1.
\]

The worst profile has input weight (w=8), packet support (H=2), terminal
state weight zero, and state-weight histogram ((1,0,0,0,1)). At distance
125, the same moment calculation gives

\[
\rho_{125}=0.414213803780\ldots>\sqrt2-1.
\]

Therefore, the unmodified Goal 02 proof composes through distance 124 and
stops at 125.

This stopping point comes from the moment relaxation. For two equal packets
of value `0xf`, the exact effective contraction at distance 124 is only
(0.3229197013). The generic moment calculation assigns (0.4137913370) to
the same profile.

## Exact one- and two-packet count

For one packet of weight (a), the accumulator output has weight at most
(D) for exactly (min\{N,\lfloor D/a\rfloor\}) packet positions. This gives
the exact failure probability after division by (N).

For two packet values (x,y\ne0), fix their order and define

\[
a=\operatorname{wt}(x),
\qquad
q=\operatorname{wt}(x+y).
\]

If their ordered positions are (p_1<p_2), set

\[
d=p_2-p_1,
\qquad
t=N-p_2+1.
\]

The four-lane accumulator output weight is exactly

\[
W=ad+qt.
\tag{2}
\]

Consequently, the number of failing supports for this order is the number of
positive integer pairs ((d,t)) satisfying

\[
d+t\le N,
\qquad
ad+qt\le D.
\tag{3}
\]

Dividing this count by \(\binom N2\) gives the exact probability for equal
packet values. For unequal values, averaging (3) for the orders ((x,y)) and
((y,x)) gives the exact probability. Exhausting the 15 nonzero four-bit
values gives the worst fixed-word probability for every (H=2) profile.

The hybrid proof uses these exact probabilities for (H\le2) and retains the
Goal 02 moment relaxation for (H\ge3). At distance 40, this replacement
reduces the uniform contraction from (0.3576340714) to (0.3340364242).
The required memories fall from 72 to 68 for 20 failure bits and from 108 to
102 for 40 failure bits.

## Complete distance-204 sweep

The transition inequality (w\le2W) makes failure impossible when (w>2D).
For (D=204), it is therefore enough to sweep (1\le w\le408). For every
(w), the computation covers

\[
\left\lceil\frac w4\right\rceil\le H\le w
\]

and all terminal weights admitted by the Goal 02 relaxation. The sweep covers
62,832 ((w,H)) profiles.

The maximum in (1) occurs at

\[
w=16,
\qquad
H=4,
\qquad
q=0,
\qquad
(m_0,m_1,m_2,m_3,m_4)=(2,0,0,0,2).
\tag{4}
\]

The moment bound for (4) is

\[
\Pr[W\le204]\le2^{-20.3573591883},
\]

whose sixteenth root is the value in (1). The profile is attainable: four
identical `0xf` packets produce state weights ((4,0,4,0)).

For this particular attainable profile, a direct support count gives exact
log-probability (-24.0992851781) and effective contraction (0.3520359516).
The moment bound loses 3.74 bits on the event. This calculation does not yet
cover every four-packet value multiset, but it shows that the new frontier also
contains substantial analytic slack.

At distance 205, this same valid moment bound has sixteenth root
(0.4142473926>\sqrt2-1). Hence the current hybrid proof does not compose at
205. This fact is not a low-weight codeword and does not refute distance 205.

## Next tightening step

The new bottleneck is the four-packet return in (4). The next proof step should
enumerate four-packet value multisets, average exactly over their distinct
orders, and count the support-gap compositions whose weighted state duration
is at most (D). The generic moment bound should remain in place for larger
supports. This is the direct four-packet analogue of the successful Goal 03
refinement and requires no change to the encoder.
