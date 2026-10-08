# Rank screens for GRS[33,22] and GRS[64,44]

Neither outer-only repair closes the present sparse proof gap. GRS[64,44]
clears the 40-bit target for eight active rows, but the same bound is
vacuous at sixteen rows and becomes much worse thereafter. GRS[33,22]
does not reach 40 bits even at eight rows. These are limitations of the
tested upper bound, not bad-code witnesses for the repaired geometries.

The objective remains relative distance 10% with setup failure at most
2^-40. This note does not lower that objective or revisit the already
established half-rate result. It reports bounded floating screens and a
diagnostic of their local return terms. No native code or benchmark was run.

## Candidates and unchanged proof interface

Both candidates retain the current width-sixteen, four-state Cauchy inner:
eta=2, independent nonzero input diagonals, and four fresh nonzero state
diagonals in each packet. The four fixed phases and one regional routing
stage are unchanged. Only the systematic GRS outer and number of regions
change.

| Outer | Rate | Parity p | K20 rows L | Encoded symbols N | Largest forbidden weight e |
|---|---:|---:|---:|---:|---:|
| GRS[33,22] | 2/3 | 11 | 47,680 | 1,573,440 | 157,343 |
| GRS[64,44] | 11/16 | 20 | 23,872 | 1,527,808 | 152,780 |

The source is padded to `L=64*ceil(2^20/(64*k))` rows and
`e=ceil(N/10)-1`. The script imports the existing two-state regional rank
bound without changing its transfer entries. States mean actual zero and
nonzero witness state. The one-, two-, and three-input return entries use
the audited joint rank-return bounds. No fixed-word probability is
substituted for a conditional rank bound.

For A active source rows, let F_h be the regional transfer summed over
h-element slot subsets and divided by C(L,h). The regional matrix and
source factor are

\[
 R_A=\sum_{h=0}^{A}\binom Ah\tau^hF_h,
 \qquad
 B_A=\binom LA\tau^{-pA-1}y^{-e}e_0^T R_A^b\mathbf1.
\]

Here b is the outer length, `0<y<=1`, and `1<=tau<=255`. The existing
logarithmic polynomial convolution accounts for routing without replacement.
The source relaxation includes infeasible supports conservatively.
Reported margin means `-log2(B_A)`, not a directed lower bound on that margin.

## Bounded sparse sweep

For each geometry and each listed A, the helper first evaluates the old
GRS[32,22] witness. It then runs one L-BFGS-B search, capped at 35 iterations.
There are no extra restarts. All twelve searches converged after 5--11
iterations. Thus a failed frozen witness, especially for the doubled
outer length, is not mistaken for the best result of this bounded search.

The marks are parameterized by
`y=exp(-exp(logtilt)*A/N)` and `tau=exp(logtau)`.

| A | GRS[33,22] frozen | GRS[33,22] optimized | GRS[64,44] frozen | GRS[64,44] optimized |
|---:|---:|---:|---:|---:|
| 4 | 95.741653 | 96.263963 | -259.085072 | 210.125209 |
| 8 | 6.230855 | 6.898277 | -112.530944 | 80.564676 |
| 16 | -105.226597 | -104.222350 | -370.780193 | -67.023181 |
| 32 | -296.732828 | -294.964151 | -878.418442 | -319.229832 |
| 64 | -626.983146 | -623.599590 | -1849.491657 | -755.359117 |
| 128 | -1213.423984 | -1208.666629 | -3747.918937 | -1559.681255 |

The optimized witnesses are:

| A | 33/22 logtilt | 33/22 logtau | 64/44 logtilt | 64/44 logtau |
|---:|---:|---:|---:|---:|
| 4 | 5.408110570083 | 5.141631096879 | 6.016522153875 | 4.955235100149 |
| 8 | 5.211180238308 | 2.981174404569 | 5.833848893086 | 2.892815739721 |
| 16 | 5.143314755516 | 2.402305424372 | 5.771379994298 | 2.332308214012 |
| 32 | 5.081133685645 | 2.181961383502 | 5.711679950113 | 2.116944864997 |
| 64 | 5.048849233719 | 2.094151235013 | 5.682437006272 | 2.035247900889 |
| 128 | 5.040511070629 | 2.067434889577 | 5.676299791410 | 2.014495743526 |

The GRS[64,44] eight- and sixteen-row results reproduce the earlier
checkpoints in `OUTER_LENGTH_TRADEOFF.md`. The extension through A=128
does not reveal a recovery of the sparse bound. No conclusion about all
unlisted row counts follows from these samples.

## Return diagnostic

After each optimization, the helper repeats the calculation at the same
marks while setting every nonzero-to-zero packet entry to zero. This
deletes real witness histories and is **invalid as an upper bound**.
Its only purpose is to quantify the contribution of the current return
cover. The original transfer is restored before any valid calculation.

| A | 33/22 with returns deleted, invalid | 64/44 with returns deleted, invalid |
|---:|---:|---:|
| 4 | 131.920521 | 274.011909 |
| 8 | 52.575105 | 167.291885 |
| 16 | 3.920869 | 141.367914 |
| 32 | -24.807120 | 201.294937 |
| 64 | -28.760686 | 394.997081 |
| 128 | 38.652492 | 873.129459 |

For GRS[64,44], deleting returns changes every tested A>=16 case from a
vacuous bound to more than 40 bits. Returns remain the dominant obstruction
to this proof calculation. GRS[33,22] also loses hundreds of bits to
returns, but deleting them is insufficient at the same marks for several
classes. These counterfactual values do not establish that a sharper
valid return bound can recover the difference.

## Lightweight input-density comparison

The existing global diagonal-refresh objective was also optimized at four
input densities alpha=W/N. It includes the return cover and the relaxed
outer row factor, but omits finite routing and boundary factors. Each
geometry uses one warm-started search per density, capped at 80 iterations.
All eight searches converged after 15--23 iterations.

The table gives the leading exponent `log2(upper envelope)/N`; negative
values are favorable. These are isolated points, not a certified interval
or an established density cutoff.

| alpha | GRS[33,22] exponent | GRS[64,44] exponent |
|---:|---:|---:|
| .12 | -.0452502361 | -.0385089537 |
| .11 | -.0307877723 | -.0252546584 |
| .10 | -.0178799540 | -.0134304311 |
| .08 | +.0028728278 | +.0055595289 |

The denser classes appear compatible with the target, as they did before
the structural rejection of GRS[32,22]. They do not resolve the missing
sparse classes. No partial directed-certificate campaign was started.

## Structural and performance context

Literal window witnesses reject some realized routes of the old
GRS[32,22], four-state geometry independently of inner coefficients. They
do not transfer automatically to either geometry in this note. Likewise,
the separate GRS[64,44] seed-one window search found D=0 against 256 exit
constraints because all initial states lay on the search's flat zero
plateau. That unsuccessful search provides no evidence of safety and no
distribution-level setup-failure bound.

The existing product-volume estimates are about 2.04545 GFNI512-equivalent
products per source record for GRS[33,22] and 2.02841 for GRS[64,44]. These
are respectively 6.51% and 5.62% above the parked 169/88 reference. The
GRS[64,44] outer also has a conservative Boolean count of 647 per four
rows, versus 170 for the old outer. These arithmetic estimates are not
timings; in particular, the longer outer's small product increase hides
a significant scheduling and Boolean-work risk.

An outer-only repair is therefore not ready for promotion on either proof
or performance grounds. The next useful gate is a structural/state-capacity
repair with a valid return law strong enough to pass an entire sparse
range, not another improvement confined to singletons or eight rows.
The separate state-eight work is not duplicated here. Any successor still
needs both an obstruction search that can escape its threshold plateau
and a complete proof before native performance work becomes justified.

## Replay

From `research/workstreams/char2_field_spin/theory`:

```text
python -B repaired_geometry_rank.py --check
python -B repaired_geometry_rank.py --active 4 8 16 32 64 128 --maxiter 35 --dense
```

The helper imports the existing regional and global-density objectives.
Its checks verify padded geometry, monotonicity of the invalid diagnostic,
and restoration of the valid transfer, including exceptional exits. They
are not a new proof of the imported transfer lemmas. The old scripts,
results, constructions, and certificates remain unchanged.
