# Goal 19: validate packet-type bounds against complete small instances

## Result

The small-instance oracle now computes the complete ensemble spectrum at
binary lengths 48 and 80. It also computes the exact low tail through weight
16 at length 112. The oracle includes every nonzero message, all local bit
permutations, the global packet permutation, and the parallel accumulator.

The comparison identifies two distinct losses. A common coefficient tilt
loses about 11 bits at output weight eight. Assigning a separate tilt to each
packet histogram recovers about 4.4 bits. The remaining coefficient-bound
loss is about 6.5--6.7 bits at the two completed typewise checkpoints.

Thus, adaptive packet-type regions have a measurable benefit. They do not
remove the complete loss from coefficient domination.

## Scaled construction

Fix \(B\in\{4,8,12\}\). A message contains \(B\) symbols in \(GF(16)\).
The field uses the polynomial \(X^4+X+1\). For double parity, define

\[
p_0:=\sum_{i=0}^{B-1}x_i,
\qquad
p_1:=\sum_{i=0}^{B-1}\gamma^{4+i}x_i,
\]

where \(\gamma=X\). The coefficients are distinct and nonzero because
\(B\le15\).

Each physical symbol is encoded by the canonical extended BCH code
\([8,4,4]\). The encoder independently permutes the eight bits in each BCH
word. It then forms two four-bit packets. One uniform permutation acts on all
packet positions before the four-lane accumulator.

The zero- and one-parity comparisons omit \(p_0,p_1\) and \(p_1\),
respectively. Every comparison retains the same data symbols and inner code.

## Exact oracle

The BCH code has spectrum

\[
A_0=1,\qquad A_4=14,\qquad A_8=1.
\]

The all-zero field symbol encodes to the zero word. The field symbol 15
encodes to the all-one word. Every other nonzero symbol has BCH weight four.

A permuted weight-four BCH word has one of three packet histograms. Their
counts among the \(\binom84=70\) binary words are

\[
2,\qquad32,\qquad36.
\]

The outer dynamic program tracks \((p_0,p_1)\in GF(16)^2\) and the three BCH
weight classes. At \(B=4\), direct enumeration of all \(2^{16}\) messages
matches this dynamic program at every parity level.

For each resulting packet histogram \(\mathbf h\), the accumulator dynamic
program computes \(A_{N,4}(\mathbf h,w)\) exactly. The oracle divides by the
exact packet-sequence count \(T_{N,4}(\mathbf h)\). It then sums the result
with rational arithmetic.

For \(B=12\), the dynamic program discards states after their output weight
exceeds 16. Accumulator output weight never decreases. Therefore the clipped
coefficients through weight 16 remain exact.

## Expected-count crossings

The table reports the first output weight whose cumulative expected number
of nonzero words is at least one.

| Data blocks | Parity symbols | Binary length | Outer histograms | Crossing |
|---:|---:|---:|---:|---:|
| 4 | 0 | 32 | 54 | 2 |
| 4 | 1 | 40 | 83 | 6 |
| 4 | 2 | 48 | 126 | 9 |
| 8 | 0 | 64 | 284 | 2 |
| 8 | 1 | 72 | 377 | 6 |
| 8 | 2 | 80 | 492 | 9 |
| 12 | 0 | 96 | 818 | 2 |
| 12 | 1 | 104 | 1,007 | 6 |
| 12 | 2 | 112 | 1,226 | 10 |

These crossings do not estimate the asymptotic distance of the target
construction. The scaled family fixes both the field size and local BCH
code. Its purpose is to validate enumeration and bounding methods.

## Common-tilt comparison

For double parity, the exact and bounded cumulative expected counts near the
crossing are as follows.

| \(B\) | \(n\) | \(D\) | Exact log count | Common-tilt log bound | Loss |
|---:|---:|---:|---:|---:|---:|
| 4 | 48 | 8 | -0.6896 | 10.2010 | 10.8906 bits |
| 8 | 80 | 8 | -1.0800 | 10.0751 | 11.1551 bits |
| 12 | 112 | 8 | -1.2385 | 10.0350 | 11.2735 bits |
| 4 | 48 | 9 | 0.6477 | 12.8291 | 12.1815 bits |
| 8 | 80 | 9 | 0.0641 | 12.8743 | 12.8102 bits |
| 12 | 112 | 9 | -0.2289 | 12.9243 | 13.1532 bits |

The common tilt does not certify the exact positive rows. Its loss grows
slowly with \(n\) at fixed \(D\) in these three instances.

## Typewise comparison

The typewise probe gives each outer packet histogram its own optimized
coefficient tilt. The probe then sums the resulting upper bounds.

| \(B\) | \(n\) | \(D\) | Common loss | Typewise loss | Improvement |
|---:|---:|---:|---:|---:|---:|
| 4 | 48 | 8 | 10.8906 | 6.4892 | 4.4014 bits |
| 8 | 80 | 8 | 11.1551 | 6.7053 | 4.4498 bits |

The recovered margin is stable across the two sizes. This evidence supports
adaptive boxes or clusters of nearby packet types. A single common tilt is
unnecessarily coarse.

The typewise optimizer reported two termination flags among 492 fits at
\(B=8\). Fallback optimization preserved the reported objective. Every
reported parameter still defines a valid coefficient upper bound. The flags
concern numerical convergence, not the direction of the inequality.

## libOTe ideas used

This implementation adopts four design choices from the libOTe enumeration
branch:

1. Keep an exact rational path as the reference calculation.
2. Label floating-point optimization separately from exact enumeration.
3. Clip the output tail only when monotonicity makes the clipping exact.
4. Prune partial states only when they cannot reach a requested final type.

No libOTe source was copied into this goal. The next optimized implementation
can reuse its multiprecision backends and discarded-mass accounting.

## Reproduction

Run the complete spectra with

```powershell
python scripts/analyze_riffle_small_exact_vs_type_bound.py --data-blocks 4
python scripts/analyze_riffle_small_exact_vs_type_bound.py --data-blocks 8
```

Run the length-112 low tail with

```powershell
python scripts/analyze_riffle_small_exact_vs_type_bound.py --data-blocks 12 --h-max 16
```

Run the typewise checkpoints with

```powershell
python scripts/probe_riffle_small_typewise_bound.py --data-blocks 4 --parity-symbols 2 --distance 8
python scripts/probe_riffle_small_typewise_bound.py --data-blocks 8 --parity-symbols 2 --distance 8
```

The next goal should replace one tilt per histogram by a small adaptive set of
tilt regions. It should reproduce most of the 4.4-bit typewise improvement
without enumerating every packet type independently.
