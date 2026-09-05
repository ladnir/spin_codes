# Goal 23: close the one-dimensional boundary faces

## Result

Every boundary face with one nonzero packet weight has a closed exact
enumerator. This removes the need for a saddle estimate or a Fourier bound on
these faces.

The exact formula covers 19% of the weight-16 boundary shell at binary length
48. Its coverage rises to 35% at binary length 136.

## Exact formula

Fix a packet weight $k\in\{1,2,3,4\}$. Consider the histogram with $h_0$
zero packets, $h_k$ packets of weight $k$, and no other packets. Let

\[
H=kh_k,
\qquad
D=H/2.
\]

If $h_k$ is odd, no boundary path returns to zero. Hence

\[
A_{N,4}(\mathbf h,D)=0.
\tag{1}
\]

Suppose $h_k=2r$. Then

\[
A_{N,4}(\mathbf h,kr)
=
\binom4k^r\binom{h_0+r}{r}.
\tag{2}
\]

## Proof

On the boundary, consecutive state supports are disjoint and their weights
sum to the input packet weight. Starting from state weight zero, a weight-$k$
packet forces a transition to state weight $k$. The next weight-$k$ packet
forces a return to zero. Thus every nonzero excursion has the form

\[
0\longrightarrow k\longrightarrow0.
\]

The outbound transition chooses one of $\binom4k$ supports. The return is
then determined. The $r$ excursions therefore have $\binom4k^r$ concrete
state realizations.

A zero packet can occur only while the state is zero. The $r$ excursions
create $r+1$ zero-state gaps: before the first excursion, between consecutive
excursions, and after the last excursion. Distributing $h_0$ zero packets
among these gaps gives

\[
\binom{h_0+r}{r}
\]

possibilities. Multiplying the two factors proves (2).

## Exact conditional probability

The number of input sequences with this packet histogram is

\[
T_{N,4}(\mathbf h)
=
\binom{h_0+2r}{2r}\binom4k^{2r}.
\]

Goal 21 shows that output weight at most $kr$ is equivalent to the boundary
event. Therefore

\[
\Pr[W\le kr\mid\mathbf h]
=
\frac{\binom{h_0+r}{r}}
{\binom{h_0+2r}{2r}\binom4k^r}.
\tag{3}
\]

Equation (3) is suitable for direct insertion into the outer histogram sum.

## Coverage in the proof gym

At $H=16$, the supported one-dimensional types use packet weight two or
four. Their combined fraction of the exact weight-16 shell is:

| Binary length | Fraction covered exactly |
|---:|---:|
| 48 | 19.0% |
| 80 | 27.9% |
| 112 | 32.4% |
| 136 | 34.7% |

The weight-two face is the largest individual type at lengths 80, 112, and
136. Closing this face exactly removes the most important sparse-face case
from the local-limit proof.

## Verification

The audit compares (1)--(2) with exact integer transfer recurrences for all
four values of $k$, every length through 64, and every compatible pair
$(h_0,h_k)$. All 8,580 comparisons pass.

Run

```powershell
python scripts/verify_riffle_boundary_one_face_formula.py --maximum-length 64
```

The receipt is `receipts/goal23_boundary_one_face_formula_audit.json`.

## Next proof step

The next simplest face uses packet weights zero, two, and four. Its boundary
state graph has only weights zero, two, and four. Weight-two excursions may
contain a run of weight-four packets before returning to zero. This excursion
grammar should yield another closed coefficient formula and eliminate the
main two-dimensional even face before any Fourier analysis is needed.
