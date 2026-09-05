# Dense finite cover for EBCH32--ParityFanout--BA

## Purpose and status

The dense first-moment sum ranges over every occupation
(Q=65,\ldots,8192) and every three-band count
(q=(q_1,q_2,q_3)) with (q_1+q_2+q_3=Q).  A single global worst-term
bound charges all compositions to the weakest cell.  That relaxation became
artificially tight at 11%.  This document records the sharper cell-local
union rule used by the active diagnostic and required from the outward
verifier.

The current receipt is still incomplete and non-interval.  The lemma below
is exact; the recorded numerical cell bounds are diagnostic witnesses.

## Cell-local union lemma

Let

\[
 \mathcal Q=\{q\in\mathbb Z_{≥0}^3:
 65\le q_1+q_2+q_3\le8192\}.
\]

For each (q\in\mathcal Q), let (U(q)\) be an upper bound on the expected
number of bad words with band counts (q\).  Let
(T_1,\ldots,T_m\subseteq\mathbb R^3) cover (mathcal Q).  For each cell,
define its integer coordinate box

\[
 R_j=\prod_{i=1}^3
 \left[\left\lceil\min_{x\in T_j}x_i\right\rceil,
       \left\lfloor\max_{x\in T_j}x_i\right\rfloor\right]\cap\mathbb Z^3.
\]

If (U(q)\le 2^{-M_j}\) for every integer point in (T_j), then

\[
 \sum_{q\in\mathcal Q}U(q)
 \le \sum_{j=1}^m |R_j|2^{-M_j}.
 \tag{1}
\]

To prove (1), assign each (q\in\mathcal Q) to one covering cell.  The
assigned points in (T_j) form a subset of (R_j).  Summing the cell bounds
gives (1).  Counting points outside (T_j), counting inadmissible points, and
counting boundary points more than once can only increase the right-hand
side.

The active verifier target stores

\[
 e_j=\log_2|R_j|-M_j
\]

and must prove

\[
 \log_2\left(\sum_j2^{e_j}\right)\le-40.
 \tag{2}
\]

The diagnostic currently requires \(e_j\le-44\) as a refinement rule and
limits the receipt to 8,192 accepted cells as a resource guard. Neither
condition replaces (2). The explicit sum in (2) is the authoritative 40-bit
test.

## Continuous-cell witness

Each tetrahedron uses one fixed Renyi/change-of-measure witness.  For that
fixed witness, the raw reference exponent and the final exponent are convex
functions of the four type counts.  The verifier must establish the following
conditions at all four vertices:

1. every reference probability is positive;
2. the raw exponent is strictly negative;
3. the coefficient of the log-multinomial term has the sign required by the
   convexity argument; and
4. the final exponent is at most (-M_j\log 2).

Convexity then extends the bound throughout the tetrahedron.  The integer box
count converts that continuous bound into the finite union contribution in
(1).

## Small-box lattice witness

Continuous interpolation can be wasteful when a tetrahedron contains only a
few integer coordinate triples.  If its coordinate box contains at most four
admissible triples, the diagnostic evaluates each triple separately.  It
stores one witness per triple and sums those point bounds.  This path does not
use convex interpolation inside the cell.

The outward verifier must distinguish the two cell types.  It must not apply
a continuous-cell witness to an enumerated box or infer a pointwise bound from
the diagnostic optimizer.

## Delegated singleton refinement

The three-band relaxation fails at the external type

\[
 (q_1,q_2,q_3)=(3344,1,21).
\]

Seven depth-50 geometric leaves contain only this integer point. A stronger
local argument splits the first band into \([24,54]\) and \([55,100]\).
For each \(a=0,\ldots,3344\), the four active counts are

\[
 (a,3344-a,1,21).
\]

The diagnostic interval cover uses 17 consecutive intervals. Each interval
stores one fixed witness and checks both endpoints. The same convexity rule
then covers every integer \(a\) in the interval.

`certify_ebch32_parityfanout_split_low_interval_cover_outward.py` verifies
all 17 intervals with 256-bit Arb arithmetic. Their aggregate satisfies

\[
 \log_2 U_{(3344,1,21)}
 \le -198.0297879425064.
\]

The main partition may delegate any leaf whose integer coordinate box is
exactly \(\{(3344,1,21)\}\). The partition audit checks this geometry. The
final dense verifier must add the split-low bound once, even if several
geometric leaves delegate the same integer point.

## Current diagnostic checkpoint

`ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json` currently
records:

- 8,060 processed tetrahedra;
- 3,957 accepted cells;
- 149 pending cells;
- zero rejected cells;
- an exact recursive partition of the three explicit root tetrahedra; and
- diagnostic aggregate accepted-box margin 43.525570 bits.

The accepted-box sum omits the 149 pending cells. It is therefore not an
upper bound on the full dense contribution.  Completion requires an empty
pending worklist and zero rejected cells.

`audit_ebch32_parityfanout_dense_cover_partition.py` reconstructs the full
binary subdivision forest from the leaf geometries. The current audit
consumes all 4,106 leaves exactly once, verifies every recorded depth, and
checks every accepted box count. It also checks the exact root six-volume
identity

\[
 545393737728+4327464960+34336575
 =549755539263=8192^3-65^3.
\]

`certify_ebch32_parityfanout_dense_cells_outward.py` uses 256-bit Arb balls.
It has verified 3,716 nondelegated nonempty accepted cells. It also checks
that seven delegated cells contain only the external type \((3344,1,21)\).
The other 234 accepted cells have empty admissible lattice boxes. Every
nondelegated nonempty cell passes its configured \(2^{-44}\) threshold. The
outward aggregate, excluding the delegated point, satisfies

\[
 \log_2 U_{\rm accepted,ordinary}
 \le -43.525569895123915.
\]

The delegated point has the separate 198.029787-bit outward bound above.
These outward bounds still omit the 149 pending leaves. They are not a bound
on the complete dense range.

## Outward-verifier obligations

For every saved continuous witness, the verifier must outward-bound:

- the three conditioned EBCH32--ParityFanout--BA Renyi band moments;
- the reference type probability for each vertex;
- the finite RM2Sub-S19 transfer raised to 16,384 epochs;
- the Chernoff factor at (d=230686);
- the fixed-witness Renyi exponent; and
- the cell's contribution after the exact integer box count.

For every enumerated box, the verifier must perform the same checks at each
stored integer triple and outward-sum the point contributions.  Finally, it
must outward-sum all cell contributions and then add the independent
(Q=1) and (Q=2,\ldots,64) receipts.

The optimizer, its convergence flag, and its nearest-binary64 objective are
not trusted.  They select legal witnesses only.

## Open implementation issue

The current optimizer restarts each child cell from generic initial points.
This causes repeated multistart searches along a narrow low-band-heavy ridge
near (Q\approx3050).  A continuation should propagate a failed parent
witness to both children as an additional start.  That change affects witness
discovery only; it does not change the bound or the future verifier.
