# RM(5,11) Q1 Envelope Status

Status: partial Q1 certificate closed; complete Q1 certificate open.

## Target

The candidate outer repeats one fixed
\(\operatorname{RM}(5,11)=[2048,1024,64]\) constituent 1024 times.  It has
no outer setup randomness.  The routing consists of independent uniform
2048-coordinate permutations in the rows and independent uniform
1024-position permutations in the transposed regions.  The inner is
RandomStepConv-M22.  The target is

\[
 k=2^{20},\qquad N=2^{21},\qquad D=228{,}590.
\]

The failure probability is over the routing and inner maps only.

## Proved low-shell contribution

The exact Kasami--Tokura formulas below twice the minimum distance give the
five nonzero shells 64, 96, 112, 120, and 124.  The formula implementation
also reproduces the authenticated RM(4,9) coefficients in the same range.

Outward arithmetic proves

\[
 \sum_{w<128} 1024 A_w T_1(w) < 2^{-83.3743314000},
\]

where \(A_w\) is the fixed constituent's weight enumerator and \(T_1(w)\)
is the occupation-one RandomStepConv-M22 transfer.  Weight 64 dominates the
sum.  This proves that the known low shells are not the finite-distance
obstruction.

## Missing finite envelope

All RM(5,11) weights are multiples of four.  A total-mass bound is sufficient
from weight 384 onward: assigning all \(2^{1024}\) constituent words to
weight 384 gives a binary64 diagnostic contribution below \(2^{-60.36}\).
Monotonicity of the inner tail must be included in any outward replay that
uses this cutoff.

Consequently, the unresolved Q1 input is an upper envelope for

\[
  A_{128},A_{132},\ldots,A_{380}.
\]

Representative per-shell caps sufficient for a 40-bit pointwise bound are:

| weight | required upper bound on \(\log_2 A_w\) |
|---:|---:|
| 128 | 231.18 |
| 160 | 312.55 |
| 192 | 398.32 |
| 224 | 488.81 |
| 256 | 584.78 |
| 320 | 796.48 |

These thresholds do not yet allocate a union budget among shells.  They are
diagnostic targets, not accepted caps.

## Inputs tested and rejected

The following statements are proved about the relaxations, not about the
true RM spectrum.

1. Type-II self-duality and nonnegativity alone are insufficient in the
   current Gleason linear program.  The relaxation permits nearly all code
   mass at weight 128.
2. Minimum distance 64 alone is insufficient.  Radius-15 Johnson-space
   packing exceeds the required caps by 227 to 328 bits at representative
   weights below 384.
3. The Kaufman--Lovett--Porat cumulative theorem is asymptotically suitable,
   but its explicit sampling proof has constants that are vacuous at eleven
   variables.
4. The transitive-code noise bounds contain an asymptotic finite-size error.
   No explicit value for that error has been established here, so those
   bounds are not certificate inputs.
5. Sampling estimates of the unknown RM(5,11) enumerator are evidence only.
   They cannot upper-bound a fixed deterministic constituent.

The 1976 Kasami--Tokura--Azumi formulas give exact coefficients below
\(2.5d=160\).  They have not been transcribed into this workstream.  Even
after that step, weights 160 through 380 still require an RM-specific bound.

## Next proof step

The best remaining route is a finite RM-specific cumulative bound for the
intermediate band.  A useful result must expose constants and produce an
integer or outward-rational envelope at \(m=11,r=5\).  Candidate mechanisms
are a recursion using the Plotkin decomposition and certified coset caps, or
a finite specialization of a modern RM weight-distribution inequality.

Do not begin Q2 or implementation benchmarking until complete Q1 closes.
