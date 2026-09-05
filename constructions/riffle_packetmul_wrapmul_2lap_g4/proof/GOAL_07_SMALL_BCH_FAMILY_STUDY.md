# Goal 07: smaller BCH family study

## Result

The smaller-BCH experiment found a faithful family and one important phase
change. The canonical primitive BCH construction reproduces the committed
64-bit parity map exactly. Its two smaller half-rate members are the extended
BCH codes ([8,4,4]) and ([32,16,8]).

An exhaustive census covered every nonzero lifted state for both smaller
members. Their 24-node distances are

\[
D_4(24)=24,
\qquad
D_{16}(24)=121.
\]

Both values exceed the scaled targets 9 and 36. The width-16 census examined
all (2^{32}-1) nonzero states. It also determined every minimum conditioned
on an anchor position and anchor weight.

These results do not support an induction to width 64. The width-64 map has
seven distinct primary components. The smaller maps have fewer components
and qualitatively different every-other-node observability. In particular,
the 24-node width-64 trace has a 31-dimensional subcode whose even outputs
all vanish. The corresponding subcodes are zero at widths 4 and 16.

This phase change explains why the smaller instances do not reproduce the
current obstruction. It also suggests a better finite target. At 38 nodes,
the even-zero subcode has dimension 10 and exact minimum weight 501. The
clean bound

\[
D_{64}(38)\ge228
\tag{A38}
\]

would close the current gap accounting. A counterexample to (A38) still has
an anchor of weight at most five. Thus the longer window removes most of the
alternating-zero freedom without enlarging the sparse-anchor set.

The study therefore supports a provisional move from the 24-node target to
the 38-node target. It does not prove (A38).

## The BCH family

For (m\in\{3,5,7\}), let (C_m) be the extended primitive narrow-sense BCH
code in the following table. The encoder uses the same cyclic coordinate
order, extension coordinate, and left-systematic convention as the committed
construction.

| (m) | Primitive BCH | Extended BCH | Width (b) | Generator polynomial |
|---:|:---|:---|---:|:---|
| 3 | ([7,4,3]) | ([8,4,4]) | 4 | `0xb` |
| 5 | ([31,16,7]) | ([32,16,8]) | 16 | `0x8faf` |
| 7 | ([127,64,21]) | ([128,64,22]) | 64 | `0xf4845518b9582a1f` |

Let (P_b:\mathbb F_2^b\to\mathbb F_2^b) be the systematic parity map of
the corresponding code. Let (A_b) be the (b)-bit prefix accumulator and
set (T_b=P_bA_b). The zero-input lifted recurrence is

\[
q_t=A_b^2a_t+A_bb_t,
\qquad
\begin{pmatrix}a_{t+1}\\b_{t+1}\end{pmatrix}
=
\begin{pmatrix}T_b&0\\T_bA_b&T_b\end{pmatrix}
\begin{pmatrix}a_t\\b_t\end{pmatrix}.
\]

For (n\ge1), define the state-based distance

\[
D_b(n)
:=
\min_{(a_0,b_0)\ne(0,0)}
\sum_{t=0}^{n-1}\operatorname{wt}(q_t).
\]

The generic parity-check construction at (m=7) produces all 64 committed
systematic parity columns without adjustment. Thus the small instances use
the same encoder convention, not merely equivalent BCH codes.

The three exponents do not form an inductive field tower. Neither 3 nor 5
divides 7. Consequently, a proof cannot embed the smaller recurrences into
the width-64 recurrence. A general result must instead use uniform properties
of the BCH encoders.

## Exact small-instance distances

The width-4 census enumerated 255 nonzero lifted states. The width-16 census
enumerated (4,294,967,295) nonzero lifted states. A Gray-code engine updated
the complete 24-node observation word after each state-bit change.

| Width (b) | Observation code | Scaled target | Exact (D_b(24)) | Minimum multiplicity |
|---:|:---|---:|---:|---:|
| 4 | ([96,8]) | 9 | 24 | 4 |
| 16 | ([384,32]) | 36 | 121 | 1 |
| 64 | ([1536,128]) | 144 | unknown | unknown |

One width-4 minimum has the periodic node-weight pattern

\[
(1,1,2,0)^6.
\]

The unique width-16 minimum starts from state `0xcd7e70d6` and has node
weights

\[
\begin{split}
(&7,4,4,5,5,2,2,7,3,6,3,4,\\
 &2,4,3,8,3,7,8,7,6,8,7,6).
\end{split}
\]

The exact width-16 zero-anchor minima range from 124 to 132 as the anchor
position varies. The weight-one-anchor minima range from 127 to 137. The
global minimum has no node of weight zero or one; its lightest nodes have
weight two.

The width-16 instance therefore validates the affine-anchor formulation. It
does not exhibit a difficult weight-one slice. This distinction is useful
for testing a decoder, but it weakens direct extrapolation to width 64.

## Why numerical extrapolation fails

The ordinary BCH minimum-distance ratios already vary substantially:

\[
\frac48=0.5,
\qquad
\frac8{32}=0.25,
\qquad
\frac{22}{128}=0.171875.
\]

The 24-node observation codes diverge more sharply. For comparison, let
(r(N,K)) be the first weight at which a random binary ([N,K]) code has
expected cumulative multiplicity at least one. The exact values used here
are:

| Width | Random crossing (r(24b,2b)) | Exact value or known upper bound |
|---:|---:|---:|
| 4 | 35 | 24 |
| 16 | 131 | 121 |
| 64 | 515 | at most 278 |

At width 16, the exact minimum is ten bits below the random crossing. At
width 64, the known 278-bit word is 237 bits below the random crossing. A
random code has expected log-count about (-365) through weight 278. Thus
the width-64 transient is a strong structural effect, not a finite-length
version of the width-16 minimum.

The autonomous maps expose the structural difference. Their minimum
polynomials factor as follows:

| Width | Factor degrees and multiplicities of the minimum polynomial of (T_b) |
|---:|:---|
| 4 | (1^4) |
| 16 | (1^7,9) |
| 64 | (1,2,4,9,10,18,20) |

The exact width-16 minimum uses both lifted primary components. The known
width-64 nine-node transient uses all seven primary components. Cancellation
between components becomes more expressive as the number of components
grows. An argument based only on the ordinary BCH weight spectrum cannot
control this cancellation.

## Every-other-node observability

The smaller family does reveal the mechanism behind the alternating-zero
obstruction. For (r\ge1), let (E_b(r)) be the rank of the map from the
lifted state to

\[
(q_0,q_2,\ldots,q_{2r-2}).
\]

The exact initial rank sequences are

\[
\begin{array}{c|l}
b& E_b(1),E_b(2),\ldots\\ \hline
4&4,7,8,8,\ldots\\
16&16,19,22,25,28,30,32,32,\ldots\\
64&64,67,70,73,\ldots,124,127,128.
\end{array}
\]

At width 64, the first 12 even samples have rank 97. This gives the
31-dimensional all-even-zero subcode in the 24-node analysis. The even
samples do not determine the lifted state until sample 23, which is output
node 44 in zero-based numbering.

Each additional even sample usually adds only three independent constraints.
This observation is exact data, not yet a proved rank formula. It persists
across all three BCH instances and is a plausible target for symbolic
analysis of the BCH generator family.

Longer windows shrink the width-64 even-zero subcode quickly:

| Window | Even samples | Kernel dimension | Exact nonzero minimum |
|---:|---:|---:|---:|
| 24 | 12 | 31 | greater than 143 in the existing certificate |
| 32 | 16 | 19 | 330 |
| 36 | 18 | 13 | 449 |
| 38 | 19 | 10 | 501 |
| 39 | 20 | 7 | 501 |
| 40 | 20 | 7 | 542 |
| 44 | 22 | 1 | 682 |
| 45 | 23 | 0 | no nonzero state |

The minima from 32 nodes onward are exhaustive Gray-code results over the
listed kernels. This table identifies the current obstruction more precisely:
24 nodes stop while the even subsequence retains 31 hidden state bits.

## Window choice

Suppose a window contains (n) nodes. The existing support-33 accounting
guarantees at least

\[
B_n
=
\left\lceil
\frac{32737-34(n-1)}{n}
\right\rceil
\]

complete windows. A clean six-bits-per-node bound is sufficient through
(n=39). Selected choices are:

| (n) | (B_n) | Proposed bound (6n) | Output margin | Even-zero kernel dimension |
|---:|---:|---:|---:|---:|
| 24 | 1332 | 144 | 3,043 | 31 |
| 32 | 991 | 192 | 1,507 | 19 |
| 36 | 877 | 216 | 667 | 13 |
| 38 | 829 | 228 | 247 | 10 |
| 39 | 807 | 234 | 73 | 7 |

For every row in this table, a violating word has total weight at most
(6n-1). Therefore some node has weight at most five. The set of possible
sparse anchor values remains

\[
\sum_{i=1}^{5}\binom{64}{i}=8{,}303{,}632.
\]

The 38-node target is the best provisional compromise. It retains the same
anchor cutoff as the 24-node target. It reduces the alternating-zero kernel
from dimension 31 to dimension 10. After fixing one anchor, it also leaves
2,368 observed coordinates rather than 1,472. Those coordinates may support
substantially more disjoint information sets.

The 39-node target gains little additional observability and leaves only a
73-bit accounting margin. At 40 nodes, the bound (6n=240) no longer
suffices. The forced anchor cutoff then increases to six if one uses the
minimum sufficient bound 241.

## Consequences for a general proof

The smaller BCH instances support three conclusions.

First, the recurrence family is well defined and exactly matches the
committed instance. The width-16 member supplies complete ground truth for
testing a sparse-anchor decoder.

Second, a theorem based only on BCH minimum distance or ordinary weight
spectrum is unlikely to close the width-64 case. Such a theorem would miss
the cancellation between primary components and the defect in even-sample
observability.

Third, a useful family theorem could instead have two parts:

1. determine the rank of each strided observation map from the BCH generator
   polynomial; and
2. bound the affine low-weight lists after a sparse output block is fixed.

The first part now has a concrete rank sequence to explain. The second part
is the existing joint sparse-anchor problem. The width-16 exact census can
validate an implementation of that decoder before it is applied at width 64.

## Recommended next goal

Replace the provisional 24-node decoder target with a 38-node preflight.
Pack disjoint information sets in the zero-anchor shortenings at positions
0, 19, and 37. Compare their counts and restriction radii with the existing
24-node values. Do not run the large enumerations during this preflight.

If the 38-node shortenings materially improve the packing, implement the
joint sparse-anchor decoder first at width 16. The exact width-16 census gives
the required oracle for anchors of weights zero and one. After the decoder
reproduces every exact conditional minimum, apply the same method to the
38-node width-64 code.

## Evidence and scope

The primary exhaustive engine is
`scripts/analyze_riffle_packetmul_wrapmul_2lap_bch32_family.cpp`. The audit is
`scripts/audit_riffle_packetmul_wrapmul_2lap_bch_family.py`. The audit
reconstructs each BCH code from its power-sum checks, verifies the committed
width-64 parity map, replays every width-16 minimum witness, and exhausts the
reported even-zero kernels.

The authenticated receipt is
`constructions/riffle_packetmul_wrapmul_2lap_g4/receipts/goal07_bch_family_audit.json`.
The raw width-16 census is
`constructions/riffle_packetmul_wrapmul_2lap_g4/receipts/goal07_bch32_family_exact_raw.json`.

The exact small-instance results do not prove a width-64 distance bound. The
38-node proposal remains a finite proof target, not an established property
of the construction.
