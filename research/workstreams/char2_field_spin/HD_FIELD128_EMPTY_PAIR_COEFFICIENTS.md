# Exact empty-pair counts in dense regions

The 40% proof still has uncovered middle-weight classes. Reusing the
weight-4096 coefficient markers gives especially poor bounds above that
weight. At weight 5120, the largest saved term has no empty routed pairs.
That boundary admits a simpler coefficient calculation without changing
the construction.

This note uses the gamma30 inner and the complete paired polynomial from
[the all-pair composition argument](HD_FIELD128_ALL_PAIR_COMPOSITION.md).
Write that polynomial as \(Q(x,y)\), where \(x\) marks singleton pairs and
\(y\) marks full pairs. Each coefficient is a nonnegative two-by-two
matrix. A composed packet pair has 64 pair slots, so its ordinary total
degree is at most 64. Its weighted input degree can be 128; these are
different degree bounds.

There are \(P=4096\) pair slots in a region and 64 composed packet pairs.
Fix \(0\le r\le2P\) and
\(\max(0,r-P)\le f\le\lfloor r/2\rfloor\), both integers.
Here \(r\) is the original active-group weight and \(f\) is the full-pair count.
The empty and singleton counts are

\[
e=P-r+f,\qquad s=r-2f.
\]

The regional cap is

\[
R(s,f)=\frac{[x^s y^f]Q(x,y)^{64}}
 {H(e,s,f)},\qquad
H(e,s,f)=\frac{P!}{e!s!f!}.
\]

These coefficients bound the inner moment; they are not the exact
codeword enumerator. Matrix products remain chronological. In each
coefficient of \(Q\), the raw packet product precedes the entrywise
fractional power. Histogram multiplicities remain outside that power.

## No empty pairs

Suppose \(e=0\), so \(f=r-P\) and \(s=P-f\). Only the ordinary
degree-64 part of every factor can contribute. Define its univariate
form by

\[
B(t)=\sum_{i,j=0}^{32}\binom{32}{i}\binom{32}{j}t^{i+j}
       (M_{32+i}M_{32+j})^{\circ\alpha},
\]

where \(M_j\) is the existing raw packet cap at occupancy \(j\).
The exact regional coefficient identity is

\[
R(P-f,f)=\frac{[t^f]B(t)^{64}}{\binom{P}{f}}.
\]

There is no additional singleton-orientation factor. The regional route
already averages each fair swap; source-support multiplicities are counted
separately in the whole-code argument.

A finite positive marker \(t\) gives the usual coefficient upper bound
\(B(t)^{64}/(t^f\binom{P}{f})\). Direct extraction can be tighter.
The directed implementation multiplies polynomials with Arb enclosures,
rounds coefficients upward, and discards degrees above the requested
coefficient after each multiplication. Nonnegative degrees make that
truncation exact for every retained coefficient. Exact positive tilting
improves numerical conditioning without changing the coefficient.
Polynomial reversal may reduce the requested degree; it does not
complement the code, local matrices, input weight, or outer factor.

## A small number of empty pairs

Define

\[
C(u,t):=u^{64}Q(1/u,t/u).
\]

Every monomial has a nonnegative power of \(u\). That power is its empty
pair count, and \(t\) still marks full pairs. Therefore

\[
[x^s y^f]Q(x,y)^{64}=[u^e t^f]C(u,t)^{64}.
\]

For any fixed \(t>0\), positivity gives the entrywise bound

\[
R(s,f)\le
\frac{[u^e]C(u,t)^{64}}{t^f H(e,s,f)}.
\]

This counts empty pairs exactly while bounding the full-pair coefficient
by one positive marker. Only degrees up to \(e\) in \(u\) are needed.
It is useful near the boundary even when both original exponents are
large. The current implementation of this extension is a FLOAT diagnostic,
not directed arithmetic.

## Evidence and scope

At the original weight-4096 witness, the weight-5120 boundary profile has
the following diagnostic log bounds:

| Calculation for that profile | Base-two logarithm |
| --- | ---: |
| Saved two-variable alias bound | +195,673.94 |
| Fresh finite marker for \(B(t)\) | +1,302.15 |
| Direct positive convolution estimate | -264.17 |

The [convolution diagnostic](theory/hd_field128_gamma30_d40_r5120_empty_boundary_univariate_float_20261010.json)
flags underflow and has no directed correction. A fresh
[256-bit directed replay](theory/hd_field128_gamma30_d40_r5120_empty_boundary_directed256_20261010.json)
then proves a log upper bound of -264.166341358368 for this profile.
It recomputes the directed packet caps, polynomial coefficient, integer
support cap, and original-weight outer factor. It covers only the profile
\((f,f)=(1024,1024)\), not every support at weight 5120.
A [separate 384-bit replay](theory/hd_field128_gamma30_d40_r5120_empty_boundary_directed384_20261010.json)
gives the same rational ceiling. An
[independent audit](theory/hd_field128_gamma30_d40_r5120_empty_boundary_directed_dual_precision_independent_audit_20261010.json)
recomputes the count, outer factor, and final transfer from the saved
regional caps. It does not regenerate those caps. The boundary support
count also agrees with the exact cyclic count
\((8192/5120)\binom{5120}{3072}\).

Separately, a fresh 111-marker bank and bounded mark search reduce the
complete weight-5120 FLOAT sum from +195,673.94 to +26,384.16.
The largest remaining profile has \(f=1160\) in both alternating regions,
hence \(e=136\) and \(s=2800\). At those same searched marks, fitting
positive coefficient markers directly for this profile gives -38,463.26.
The exact-empty-count convolution estimate gives -39,711.28.
The [interior diagnostic](theory/hd_field128_gamma30_d40_r5120_full1160_empty_layers_float_20261010.json)
again covers only that profile and is not directed.

Adaptive fitting at the searched marks adds six finite markers, targeting
the current largest profile after each addition. The resulting
[117-marker complete FLOAT sum](theory/hd_field128_gamma30_d40_r5120_dense_adaptive_float_20261010.json)
has log bound -12,000.77997. All 2,362,369 profiles remain in each of the
seven completed sums. This result alone is a numerical screen, not a
directed weight certificate. Eighteen selected markers retain every
regional entrywise minimum.
A [fresh 18-marker replay](theory/hd_field128_gamma30_d40_r5120_winner18_replay_20261010/replay.json)
reproduces both the regional-array hash and the full FLOAT score exactly.
The subsequent directed band calculations use freshly recomputed
coefficients, not these FLOAT values; the [ladder index](HD_FIELD128_LADDER.md)
records their completed coverage separately.

These checks separate poor marker placement from losses in the local
moment caps. They do not establish low-weight codewords or require a new
encoder. The [partial-cover index](theory/hd_field128_gamma30_d40_partial_cover_v18_20261010.json)
contains only complete directed weight events. It includes the six new
bands from 4832 through 5215, not the isolated boundary-profile receipt.

The isolated helpers are
[boundary FLOAT extraction](theory/hd_field128_gamma_empty_boundary_screen.py),
[directed boundary extraction](theory/hd_field128_gamma_empty_boundary_certificate.py),
and [near-boundary FLOAT extraction](theory/hd_field128_gamma_empty_layers_screen.py).
Each has tiny independent polynomial tests and a separate source review.
