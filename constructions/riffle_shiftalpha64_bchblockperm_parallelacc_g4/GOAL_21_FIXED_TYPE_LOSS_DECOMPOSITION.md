# Goal 21: locate the loss for one fixed packet type

## Result

The scalar output-tail bound is not causing the remaining gap at the
weight-16 shell and output threshold eight. It is exact. All measured loss
comes from replacing one fixed multivariate packet coefficient by evaluation
at positive fugacities.

Direct optimization of the boundary matrix measures losses of 4.715 bits at
binary length 48 and 5.088 bits at binary length 80. Thus further refinement
of the scalar output tilt cannot close the current proof-gym gap. The next
bound must retain a local coefficient factor,
or an equivalent anti-concentration factor, for the packet histogram.

## Boundary identity

Fix an accumulator input with total weight

\[
H:=\sum_{t=1}^N\operatorname{wt}(u_t)
\]

and output weight

\[
W:=\sum_{t=1}^N\operatorname{wt}(q_t).
\]

Let (I_t) be the size of the intersection of the supports of
(q_{t-1}) and (q_t). Since (u_t=q_{t-1}+q_t),

\[
H=2W-\operatorname{wt}(q_N)-2\sum_{t=1}^N I_t.
\tag{1}
\]

Suppose (H=2D) and (W\le D). Equation (1) forces all three conclusions

\[
W=D,
\qquad q_N=0,
\qquad I_t=0\quad\text{for every }t.
\tag{2}
\]

Therefore a bad path on this shell has exactly weight (D), returns to the
zero state, and uses disjoint consecutive state supports. This is stronger
than the earlier cutoff (H\le2W).

For packet width (g), define the boundary matrix

\[
L_g(\mathbf x)_{a,b}:=
\begin{cases}
\binom{g-a}{b}x_{a+b},&a+b\le g,\\
0,&a+b>g.
\end{cases}
\tag{3}
\]

The binomial coefficient chooses the new support disjointly from the old
support. If a packet histogram \(\mathbf h\) has total binary weight
\(\sum_k k h_k=2D\), then (2) gives the exact formula

\[
A_{N,g}(\mathbf h,D)
=[\mathbf x^{\mathbf h}]
e_0^{\mathsf T}L_g(\mathbf x)^N e_0.
\tag{4}
\]

Equation (4) eliminates the output variable (z) and excludes every
non-boundary state path before coefficient extraction.

## Why scalar Chernoff is exact here

For a fixed packet histogram, let

\[
F_{\mathbf h}(z):=\sum_w A_{N,g}(\mathbf h,w)z^w.
\]

When (H=2D), equation (1) implies that (A_{N,g}(\mathbf h,w)=0) for
(w<D). Hence

\[
\inf_{0<z<1}z^{-D}F_{\mathbf h}(z)
=A_{N,g}(\mathbf h,D).
\tag{5}
\]

The equality is a limit as (z\) tends to zero. Thus an exact extraction of
the packet coefficient followed by scalar Chernoff has no loss on this
shell.

## Exact-versus-bound experiment

Both instances use double parity and (D=8). The eligible remainder has
(H=16=2D) and fifteen packet histograms.

| Binary length | Exact shell log | Scalar-Chernoff log | Full bound log | Boundary bound log | Boundary coefficient loss |
|---:|---:|---:|---:|---:|---:|
| 48 | -3.2558 | -3.2558 | 1.4811 | 1.4594 | 4.7152 bits |
| 80 | -3.1172 | -3.1172 | 1.9828 | 1.9709 | 5.0881 bits |

Every one of the fifteen packet histograms has conditional minimum output
weight eight. There are no secretly impossible types for the coefficient
bound to discard. Across the types, the individual boundary-coefficient
losses range from 1.276 to 7.193 bits at length 48 and from 1.646 to 8.039
bits at length 80.

The exact shell values use integer accumulator enumeration and rational
outer weights. The scalar equality follows from (1). The multivariate values
use floating-point numerical optimization and are upper bounds, not proof
certificates against numerical error.

## Consequence for the proof

The current hierarchy of losses is now separated.

1. The hard cutoff (H\le2D) removes all higher input-weight shells exactly.
2. On the boundary (H=2D), scalar output-tail Chernoff is exact.
3. The remaining 4.7--5.1 bits arise because positive-fugacity evaluation
   omits the probability of landing on the requested packet histogram.

For fixed (N) and (H), the histogram obeys both
(\sum_kh_k=N) and (\sum_kkh_k=H). It consequently has only three free
coordinates at width four. The present positive-fugacity evaluation of (4)
also ranges over other values of (H), so it uses four effective packet
tilts after fixing one scale. A sharper route can first extract or condition
on (H), then apply a three-dimensional local coefficient or
anti-concentration correction. Equation (4) is the appropriate transfer
matrix for deriving that correction.

The small instances do not establish how the missing local factor scales
with (N). A saddle-point or local central-limit approximation may predict
it, but the proof needs a certified upper bound rather than an asymptotic
heuristic.

## Reproduction

Run the two probes sequentially:

```powershell
python scripts/probe_riffle_small_one_type_loss.py --data-blocks 4
python scripts/probe_riffle_small_one_type_loss.py --data-blocks 8
```

The receipts are
`receipts/goal21_one_type_loss_b4_p2_h16_d8.json` and
`receipts/goal21_one_type_loss_b8_p2_h16_d8.json`.

The next goal should derive and test a certified local coefficient factor for
(4). A useful first checkpoint is an exact-versus-saddle-point comparison on
all fifteen types at lengths 48 and 80, followed by a finite-size inequality
that exposes its dependence on the three free histogram coordinates.
