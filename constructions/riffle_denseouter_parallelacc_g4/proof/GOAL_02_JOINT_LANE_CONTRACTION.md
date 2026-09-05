# Goal 02 proof: joint-lane contraction

## Result

The joint-lane bound passes. At binary length \(2{,}097{,}408\) and distance
threshold \(40\), every fixed outer word of bit weight \(w\) satisfies

\[
\Pr_\Pi[W\le40]\le\rho^{w},
\qquad
\rho=0.3576340715.
\tag{1}
\]

The one-lane proof used \(\rho=0.4119542067\). Inserting (1) into the repaired
dense-outer generating function reduces the required outer memory as follows.

| Failure target | One-lane memory | Joint-lane memory |
|---:|---:|---:|
| 20 bits | 90 | 72 |
| 40 bits | 129 | 108 |

No encoder mechanism or randomness changed.

## Exact gap moment

Fix a packet word with \(H\) nonzero packets. Fix one order
\(v_1,\ldots,v_H\in\mathbb F_2^4\setminus\{0\}\), and define

\[
s_j:=\sum_{i=1}^j v_i,
\qquad
a_j:=\operatorname{wt}(s_j),
\qquad
w:=\sum_{j=1}^H\operatorname{wt}(v_j).
\]

Let \(r_0\) be the number of zero packets before the first active packet. For
\(j\ge1\), let \(r_j\) be the number after active packet \(j\) and before the
next active packet. The final gap is \(r_H\). Then

\[
r_0+\cdots+r_H=N-H,
\qquad
W=\sum_{j=1}^H a_j(1+r_j).
\tag{2}
\]

The active support is a uniform \(H\)-subset of \([N]\). Therefore, for
\(0<z,\tau<1\),

\[
\begin{aligned}
\mathbb E_\Pi[z^W\mid v_1,\ldots,v_H]
&=
\frac{1}{\binom NH}
\sum_{r_1+\cdots+r_H\le N-H}
\prod_{j=1}^H z^{a_j(1+r_j)}\\
&\le
\frac{\tau^{-(N-H)}}{\binom NH}
\prod_{j=1}^H
\frac{z^{a_j}}{1-\tau z^{a_j}}.
\end{aligned}
\tag{3}
\]

The inequality sums the leading gap exactly. For
\(r_1+\cdots+r_H\le N-H\), multiplication by
\(\tau^{r_1+\cdots+r_H-(N-H)}\) is at least one. This avoids the unnecessary
factor \((1-\tau)^{-1}\) in the earlier coefficient bound.

The Chernoff step gives

\[
\Pr_\Pi[W\le D]
\le
z^{-D}
\frac{\tau^{-(N-H)}}{\binom NH}
\prod_{j=1}^H \phi_{z,\tau}(a_j),
\qquad
\phi_{z,\tau}(a):=\frac{z^a}{1-\tau z^a}.
\tag{4}
\]

Equation (4) holds for every active-value order. It therefore also holds after
averaging over the order induced by the packet permutation.

## Uniform relaxation of the state path

Write

\[
A:=\sum_{j=1}^H a_j,
\qquad
q:=a_H.
\]

Every actual state path satisfies four constraints.

1. Since \(v_j=s_{j-1}+s_j\), the triangle inequality gives
   \(w\le2A-q\).
2. Since every \(v_j\ne0\), two consecutive states cannot both be zero.
   Also \(s_1\ne0\). Hence at least \(\lceil H/2\rceil\) states are nonzero.
3. The parity of Hamming weight is linear over XOR, so \(q\equiv w\pmod2\).
4. We have \(0\le q\le\min\{4,w\}\). If \(H=1\), then \(q=w\).

For fixed \(z\) and \(\tau\), the function

\[
f(a):=\log\phi_{z,\tau}(a)
\]

is decreasing and convex on \([0,4]\), because

\[
f'(a)=\frac{\log z}{1-\tau z^a}<0,
\qquad
f''(a)=
\frac{\tau z^a(\log z)^2}{(1-\tau z^a)^2}>0.
\tag{5}
\]

Consequently, the relaxed product in (4) is largest when:

- the number of nonzero states is as small as the constraints permit;
- their total weight is as small as the constraints permit; and
- the remaining weight is concentrated into weights four, one residual
  weight, and weights one.

This extremal histogram is computed directly from \((w,H,q)\). It is an upper
relaxation: the histogram need not correspond to an actual state path.

## Finite profile sweep

If \(w>2D\), equation \(w\le2W\) makes the failure event impossible. For
\(D=40\), it is therefore enough to cover \(1\le w\le80\).

For every \(w\), the script covers

\[
\left\lceil\frac w4\right\rceil\le H\le w
\]

and every terminal weight allowed by the four constraints above. There are
2,480 \((w,H)\) profiles. For each relaxed histogram, the script chooses
\(z,\tau\in(0,1)\) and evaluates (4).

The largest \(w\)-th root is

\[
\rho=0.3576340714375845.
\]

It occurs at

\[
w=8,
\qquad
H=2,
\qquad
q=0,
\qquad
\(m_0,m_1,m_2,m_3,m_4)=(1,0,0,0,1)\).
\tag{6}
\]

The optimized bound for (6) is

\[
\Pr[W\le40]\le2^{-11.8675512617}.
\]

This profile is attainable: two identical weight-four packets produce state
weights \(4,0\). Thus the relaxation identifies the expected packet-clumping
obstruction rather than a fictitious state sequence. For comparison, its
exact failure probability is \(2^{-14.6782603859}\); the moment bound retains
about 2.81 bits of slack.

Every numerical parameter returned by the optimizer gives a valid upper bound
through (4). An optimizer status warning can only leave a weaker bound. It
cannot invalidate the reported uniform maximum.

## Composition with the dense outer

Let \(B_M(\rho)\) be the corrected support-cluster quantity from Goal 01:

\[
B_M(\rho)=
\frac{\rho a(\rho)^{M+1}}{1-\kappa_M(\rho)},
\quad
a(\rho)=\frac{1+\rho}{2},
\quad
\kappa_M(\rho)=\rho\sum_{d=1}^{M+1}a(\rho)^d.
\]

Since (1) holds for \(w\le80\) and failure is impossible for \(w>80\), the
first-moment argument gives

\[
\Pr[d_{\min}\le40]
\le
\exp(kB_M(\rho))-1.
\tag{7}
\]

At total binary length \(2{,}097{,}408\), equation (7) gives

| Outer memory \(M\) | Bound | Failure bits |
|---:|---:|---:|
| 40 | \(2.14\times10^{-1}\) | 2.22 |
| 56 | \(3.95\times10^{-4}\) | 11.31 |
| 72 | \(8.02\times10^{-7}\) | 20.25 |
| 90 | \(7.51\times10^{-10}\) | 30.31 |
| 108 | \(7.04\times10^{-13}\) | 40.37 |

## Scope

Goal 02 proves a finite ensemble bound for the stated length and threshold. It
does not certify a particular sampled matrix. It does not improve the
asymptotic relative-distance constant until the finite profile relaxation is
converted into a uniform asymptotic formula.

The remaining numerical slack is concentrated in the attainable two-packet
return. Improving generic high-support estimates cannot remove that profile.
