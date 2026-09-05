# Goal 20: compress typewise tilts into adaptive regions

## Result

A small number of tilt regions recovers most of the improvement from one
tilt per packet histogram. At output weight eight, eight regions recover
87--89% of the available improvement. Sixteen regions recover 98--99%.
The behavior is stable between binary lengths 48 and 80.

The experiment also reveals a sharper proof split. Most important packet
types lie in the minimum outer binary shell. After treating that shell
exactly, the remaining loss comes from coefficient domination within the
next shell. It does not come from the number of tilt regions.

## Region bound

Fix a partition \(\mathcal P\) of the outer packet histograms. For each
region \(S\in\mathcal P\), choose positive packet fugacities
\(\mathbf x_S\) and an output fugacity \(0<z_S<1\). Equation (2) of Goal 18
gives a valid bound for every histogram in \(S\). Summing within each region
and then across the partition gives

\[
\sum_{S\in\mathcal P}
z_S^{-D}e_0^{\mathsf T}M(\mathbf x_S,z_S)^N\mathbf1
\sum_{\mathbf h\in S}
\frac{C(\mathbf h)\mathbf x_S^{-\mathbf h}}
{T_{N,4}(\mathbf h)}.
\tag{1}
\]

Every fixed partition and every fixed set of fugacities makes (1) a valid
upper bound. Numerical optimization only searches for smaller valid values.

The probe first optimizes one tilt for each packet histogram. It clusters the
resulting five-dimensional parameter vectors. The clustering weights are
proportional to the optimized contribution of each type. The probe then
re-optimizes one common tilt for every cluster.

The learned partition is diagnostic. It is not yet a scalable rule because
constructing it requires the individual type optima.

## Region-count experiment

Both experiments use double parity and output threshold \(D=8\).

| Regions | \(n=48\) bound log | Recovered improvement | \(n=80\) bound log | Recovered improvement |
|---:|---:|---:|---:|---:|
| 1 | 10.2010 | 0.0% | 10.0751 | 0.0% |
| 2 | 9.0143 | 27.0% | 9.1776 | 20.2% |
| 4 | 7.4783 | 61.9% | 7.3252 | 61.8% |
| 8 | 6.2994 | 88.6% | 6.1893 | 87.3% |
| 16 | 5.8286 | 99.3% | 5.6920 | 98.5% |
| One per type | 5.7996 | 100% | 5.6253 | 100% |

Thus the four-dimensional histogram domain does not require one independently
optimized bound per integer point. A constant-size region family captures
most of the observed variation at these sizes.

## Deterministic input-weight cutoff

Let \(u_t\) be accumulator input packet \(t\). Let \(q_t\) be the state after
that packet, with \(q_0=0\). Define

\[
H:=\sum_{t=1}^N\operatorname{wt}(u_t),
\qquad
W:=\sum_{t=1}^N\operatorname{wt}(q_t).
\]

Since \(u_t=q_{t-1}+q_t\), the triangle inequality for Hamming weight gives

\[
\begin{split}
H
&\le\sum_{t=1}^N
\left(\operatorname{wt}(q_{t-1})+\operatorname{wt}(q_t)\right)\\
&=2W-\operatorname{wt}(q_N)
\le2W.
\end{split}
\tag{2}
\]

Therefore an output with \(W\le D\) must have outer binary weight
\(H\le2D\). This condition is exact and independent of the random
permutations.

At \(D=8\), the small double-parity construction has eligible outer weights
12 and 16 only. Every higher outer type has zero contribution to the bad
tail. The coefficient bound does not automatically enforce this hard zero,
so the proof should apply (2) before any tilt optimization.

## Exact minimum shell and clustered remainder

At \(n=48\), outer weight 12 contributes 83.1% of the exact bad tail. At
\(n=80\), it contributes 75.6%. The exact logs are:

| \(n\) | Complete exact tail | Exact weight-12 shell |
|---:|---:|---:|
| 48 | -0.6896 | -0.9564 |
| 80 | -1.0800 | -1.4828 |

After inserting the exact weight-12 contribution, the only possible
remainder at \(D=8\) has outer weight 16. There are fifteen such packet
histograms. Giving all fifteen separate tilts yields:

| \(n\) | Hybrid bound log | Loss over complete exact tail |
|---:|---:|---:|
| 48 | 1.7255 | 2.4151 bits |
| 80 | 2.1078 | 3.1878 bits |

The remaining loss persists after the region count reaches the number of
eligible remainder types. Hence more clustering cannot remove it. The next
analysis must strengthen the conditional bound for one fixed packet
histogram.

## Consequence for the target proof

The current evidence supports the following proof structure:

1. Apply the deterministic cutoff \(H\le2D\).
2. Evaluate the lowest outer shells exactly.
3. Partition the remaining eligible types into a small number of tilt
   regions.
4. Strengthen the one-type coefficient bound before scaling the region
   calculation to the target length.

The region experiment resolves the combinatorial \(N^4\) concern at the
numerical level. It does not yet provide a certified target-length bound.

## Reproduction

Run the learned-region probes with

```powershell
python scripts/probe_riffle_small_clustered_tilts.py --data-blocks 4 --parity-symbols 2 --distance 8
python scripts/probe_riffle_small_clustered_tilts.py --data-blocks 8 --parity-symbols 2 --distance 8
```

Run the sparse-isolation probes with

```powershell
python scripts/probe_riffle_small_sparse_isolation.py --data-blocks 4 --parity-symbols 2 --distance 8
python scripts/probe_riffle_small_sparse_isolation.py --data-blocks 8 --parity-symbols 2 --distance 8
```

Run the exact-shell hybrid probes with

```powershell
python scripts/probe_riffle_small_exact_shell_hybrid.py --data-blocks 4 --parity-symbols 2 --distance 8
python scripts/probe_riffle_small_exact_shell_hybrid.py --data-blocks 8 --parity-symbols 2 --distance 8
```

The next goal should decompose the remaining one-type loss into two parts:
the scalar output-tail Chernoff loss and the multivariate packet-coefficient
loss. That decomposition will identify whether a saddle-point prefactor or a
sharper output-tail calculation is the appropriate repair.
