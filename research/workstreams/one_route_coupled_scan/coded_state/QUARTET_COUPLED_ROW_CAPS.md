# Couple the population before coarsening input counts

The construction remains fixed. The fine joint-cohort calculation passes
the392..424 band numerically. The next calculation covers all input counts
0..512, but uses broad incoming count cells. Its first implementation is
too loose. This note identifies a precise loss and gives a local repair.
Its subsequent [regional composition](QUARTET_CAPPED_REGIONAL.md) covers all
counts at the selected r,c but still fails. This note preserves the local
derivation and the preceding checkpoint.

## What the coarse calculation covers

`quartet_mixed_counts.py` groups incoming counts into intervals0..63,
64..127,...,448..511, and the singleton512. Both age cohorts use width64
buckets. The resulting165 state cells have7725 outgoing edge labels.
Every edge retains its exact next input count. Every count path, regional
payment allocation, and terminal state is included at r=2100,c=111572.

The original normalization inherited from the narrow band gives the FLOAT
log2 expression+799982.94. Retuning the positive-vector count marker and
the payment marker gives+688051.01. Neither value is a certificate or a
lower bound on actual failure. The latter calculation places its maximum
regional row at incoming cell a in320..383,n1 in0..63,n2 in0..63.

The final replay reproduces+688051.0077459727 in116.60s. Its count marker
is log2(z)=-0.5626220703125 and payment marker is
log2(y)=-4.456752670299854. A common log scale is removed before exponentiating
the tilted matrix, and restored in the global factor. This avoids silently
losing extreme edges during normalization. The replay result is unchanged.

The scripts use upward positive floors rather than deliberately deleting
tiny coefficients. A separate bound places the added counter-floor mass
more than750bits below the regional norm in the retuned test. Thus that
particular numerical safeguard does not explain the remaining large loss.
Directed evaluation is still required. The count-marker proposal need not
converge: its maximum row ratio is used, not an asserted eigenvalue.

## Independent endpoint choices admit too many positions

Fix a cell with a in[a_lo,a_hi], and temporarily set every potential and
moment weight to1. The actual incoming histogram has w-a inactive positions
and a active positions, for a total of exactly w.

The factored coarse bound instead combines an inactive-coordinate factor
at a_lo with an active-coordinate factor that may use a_hi. When the
selected count k is at most a_lo, Vandermonde's identity gives the combined
selection multiplicity

    sum_u binom(a_hi,u) binom(w-a_lo,k-u)
      = binom(w+a_hi-a_lo,k).

The true population would give binom(w,k). Both endpoint choices are
individually upper bounds, so this does not invalidate the coarse bound.
It explains why the bound can be much too large.

For w=512,a_lo=320,a_hi=383, the relaxed population has575 positions.
Take outgoing count l=408 and selected count k=208. If J_(l,k) is the
exact quartet mask multiplicity, integer arithmetic proves

    J_(408,208) binom(575,208)
      / [binom(512,208) binom(512,408)] >= 2^40.

This is already a lower bound on one unweighted relaxed coefficient,
relative to the entire true unweighted row at l=408. It is an obstruction
to this factorization's tightness, not to the encoder's distance.

## A total-row bound keeps all512 positions coupled

Use the established positive potential phi(n)=product_i h_i^n_i and
weighted selection factors d_i. For an incoming cell i=(A,R,S), let
V_i contain the vertices of

    a in A, r in R, s in S, r+s<=a.

At z>0, put L_i(z)=(1+d_i z)/h_i and L_old=max(L_3,L_4). Define

    G_(i,k)(z) = z^(-k) max_((a,r,s) in V_i)
      L_0^(w-a) L_1^r L_2^s L_old^(a-r-s).

For every histogram n in cell i, its coefficient
[X^k] product_j(1+d_j X)^n_j / phi(n) is at most G_(i,k)(z).
Choose Gamma_(i,k) by evaluating this bound at any positive z, or by
minimizing over evaluation points. Take the vertex maximum first. This
retains n0+n1+...+n4=w at every vertex. A point-mass correction can sharpen
Gamma using a variance lower bound valid at every vertex and both old types.

Define the total-row cap at exact outgoing count l by

    g_(i,l) = sum_k J_(l,k)/binom(w,k) h1^(l-k) Gamma_(i,k).

For every actual incoming histogram in cell i, the sum of its weighted
transition coefficients to all edges labelled l is at most g_(i,l).
This is an additional constraint on the joint-cohort transfer; it does not
discard its individual outgoing-cohort bounds.

At the same incoming cell used by the full-count diagnostic, the full-width
FLOAT total-row comparison is:

| Outgoing count l | Uncoupled log2 row sum | Coupled log2 cap | Difference |
| --- | ---: | ---: | ---: |
| 320 | 422.64012 | 421.63103 | 1.00909 |
| 384 | 464.94364 | 430.74322 | 34.20042 |
| 392 | 468.91733 | 429.27558 | 39.64175 |
| 408 | 477.17620 | 428.98872 | 48.18748 |
| 424 | 493.25313 | 442.03174 | 51.22138 |
| 448 | 522.70843 | 459.14418 | 63.56425 |
| 512 | 525.83258 | 439.81983 | 86.01276 |

These are local total-row factors, not per-update security margins. They
cannot be multiplied by4096 and called a whole-code improvement.

## Combine total caps with the existing edge bounds

For one incoming cell and one exact count l, let m_e be the current
coefficientwise bounds on its outgoing edges. Let g be the coupled total
cap, optionally replaced by min(g,sum_e m_e). For nonnegative continuation
values f_e, define

    F(f) = max { sum_e k_e f_e : 0<=k_e<=m_e, sum_e k_e<=g }.

Every actual weighted transition row is feasible. Therefore F(f) bounds
its continuation value. The finite linear program has the equivalent
threshold expression

    F(f) = min_(theta>=0) [ theta*g + sum_e m_e max(f_e-theta,0) ].

An optimizer is among0 and the values f_e. Equivalently, fill the largest
f_e first, up to its cap, until the total allowance is exhausted. Any
chosen nonnegative threshold gives an upper bound; several thresholds
can be tested and their minimum retained.

Do not rescale every m_e to make their sum g. A total-mass bound does not
identify which outgoing edges carry that mass, and that rescaling need
not dominate the actual row.

The exact operator F is monotone, positively homogeneous, and subadditive.
These properties allow its use in a backward regional counter recurrence.
For each exact l, apply F_(i,l) to the continuation values at total
deviation s-(l-408), then sum over l. The same induction that validates
the linear bucket transfer validates this upper recurrence. The initial
row remains the true zero-state row; no edge returns to its dummy cell.
Terminal weights are still applied only once.

A positive vector can normalize this operator using its maximum component
ratio. After sixteen updates, apply the exact regional binomial factors
and the payment marker. This preserves all routing quotas. The threshold
operator and its numerical normalization must be checked before reusing
any error estimate derived specifically for the earlier linear operator.

## Checks and remaining implementation

`quartet_mixed_count_tests.py` checks8692 exact coarse-count edge bounds,
81 exact two-region quota pairs,10107 polynomial coefficients against the
interval FLOAT envelope, and all7725 full-width edge labels. It also checks
4805 coefficients against the coupled-population FLOAT bound and180 exact
vertex objectives. The threshold formula agrees with81 exact small linear
programs. The40-bit population-mismatch lower bound uses exact integers.

The full-count driver has a300s watchdog and a2GiB conservative array bound.
The total-row screen takes5.31s; the combined audit takes0.33s. Runs are
remote, serial, single-threaded, bounded, and file-free. A command-line
length failure occurred before remote launch; the loader was shortened to
send only required helper definitions. No stalled job was restarted.

The subsequent regional calculation implements this next step on the
existing all-count cells. Its FLOAT expression improves to+429554.70 but
still fails. See the linked checkpoint for the exact toy checks, coupled
outgoing-age refinement, and next discriminating test. The fine392..424
checkpoint is preserved. Other covered counts, row counts, and directed
whole-code verification remain open.
