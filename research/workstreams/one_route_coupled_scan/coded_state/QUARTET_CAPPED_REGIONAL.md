# Compose the coupled row caps through a region

The encoder is unchanged. This calculation tests whether the coupled
population cap repairs the coarse all-count bound at r=2100,c=111572.
It is a FLOAT diagnostic, not a minimum-distance certificate.
The earlier [population calculation](QUARTET_COUPLED_ROW_CAPS.md) supplies
both the individual outgoing-edge caps and their shared total cap.

## Retain both constraints

Fix a source cell i and an exact next input count l. Write m_e for its
individual edge bounds and g for its total bound. Replace g by
min(g,sum_e m_e), and put beta_e=min(1,m_e/g). For nonnegative continuation
values f_e, use

    F_(i,l)(f) = g min_(theta>=0)
      [theta + sum_e beta_e max(f_e-theta,0)].

This is the support function of 0<=k_e<=m_e, sum_e k_e<=g.
Every actual weighted transition row is feasible. The implementation sorts
f to propose a threshold, then evaluates the dual expression. Any proposed
nonnegative threshold gives an upper bound. It also takes the minimum with
g max_e f_e, which is a second valid bound. It never scales all edge caps
down proportionally.

Let v be any positive vector and z>0. Set

    rho = max_i [sum_l z^(l-408) F_(i,l)(v)] / v_i.

No convergence assumption about the proposed vector is needed. Normalize
each operator by z^(l-408)/(rho v_i), inserting v at the destination.
The sum of these normalized operators sends the all-one vector to at
most itself.

## Count before applying the regional quota

For a total input count b, initialize U_0(b,i)=1 if b=0 and zero otherwise.
For j=1,...,16, compute

    U_j(b,i) = sum_l normalized_F_(i,l)(U_(j-1)(b-l,.)),

with zero continuation for unreachable totals. Each edge retains its exact
l, even when several edges enter the same coarse source cell next time.
Monotonicity and positive homogeneity prove coefficientwise domination by
induction. Subadditivity permits the same construction at successive
regional boundaries. We do not claim that the sum over all b is at most1:
the maximizing row may differ for different requested totals.

After sixteen updates, weight total b=6528+s by the exact regional factor

    binom(6528+s,436+s) / binom(8192,2100),

and by the payment marker. The 256 regions together require sum s=-44.
The initial row is the actual all-zero-state row, represented by an extra
source cell with no incoming edges. The free-terminal factor is applied
once after the complete stream. Neither state age nor the terminal cost
is reset at region boundaries.

## Checks and numerical scope

`quartet_capped_regional_tests.py` passes81 exact dual evaluations and81
exact two-region quota comparisons on the small model. The production
counter agrees with170 exact capped coefficients, including a nonconstant
terminal factor. The replay took0.13s on Peach and wrote no files. It was
launched only after confirming that no previous remote job remained live;
the earlier tool result was unavailable.

The full-width driver uses165 source cells and7725 exact edge labels. It
includes every next count0..512 and every regional payment allocation for
the selected r,c. It has a1200s watchdog, a2GiB conservative overall array
budget, and a separate512MiB counter-workspace guard. Runs are serial,
single-threaded, and remote; source is loaded in memory and no files are
generated. Small positive floors prevent intentional deletion of tails.
These safeguards do not replace directed arithmetic.

The full-width run completes in750.39s. Its FLOAT log2 expression is
**+429554.70068170864**, compared with+688051.0077459727 for the preceding
uncapped calculation. It still fails. These are failed upper expressions,
not lower bounds on failure or security margins. Other r,c values, directed
verification, and the complete45%/40-bit certificate remain open.

The selected count marker is log2(z)=-0.17627457812105682 and the payment
marker is log2(y)=-4.124545331751425. The proposed vector has residual1.00668;
the maximum ratio, not an asserted eigenvalue, is used. The worst regional
source cell has a=512,n1 in128..191,n2 in64..127. The counter itself takes
622.07s. Observed process RSS during the run is about114MiB. No files are
written. The prior linear counter-floor estimate is not claimed for this
nonlinear calculation; its final expression includes the added floors.

## Keep the outgoing ages coupled as well

The total cap alone does not describe how mass is divided among outgoing
age cells. A joint moment can bound those cells without introducing more
states. Fix the actual incoming histogram n and a next input count l.
If v_j positions of type j are selected into fully active quartets, write
k=sum_j v_j. The next first and second cohort counts are

    r' = l-k+v_0,   s' = v_1.

Thus u=l-r'=sum_(j>=1) v_j. For positive markers x,y, replace the selection
weights by

    (d_0, x*y*d_1, x*d_2, x*d_3, x*d_4).

The same coupled-population coefficient calculation now bounds the entire
row moment sum_e M_e x^(l-r'_e) y^(s'_e). Call that uniform bound G_(i,l)(x,y).
For a destination cell R=[r_lo,r_hi], S=[s_lo,s_hi], its edge mass is at most

    G_(i,l)(x,y) x^(-l)
      max_(r in R) x^r max_(s in S) y^(-s).

This follows term by term from positivity. The endpoint maxima can include
infeasible pairs; that only weakens the upper bound. First maximize over
all actual incoming histograms, then minimize over chosen marker pairs.
At x=y=1 this recovers the total cap. Other pairs bound the allocation of
mass to age cells, while retaining the same population at every source
vertex. The resulting edge bounds can be intersected with the original
edge bounds and used together with the total cap.

`quartet_joint_moment_caps.py` implements a local test of this refinement.
It does not yet replace the edge bounds in the full regional calculation.
Tests pass3393 exact joint moment identities and3555 exact destination-cell
bounds. Another3393 comparisons test the modified coefficient weights
against the FLOAT coupled envelope. A25-marker local screen at the earlier
source cell a320..383,n1,n2 in0..63 gives:

| Next count | Old edge-sum log2 | Joint-moment edge-sum log2 |
| --- | ---: | ---: |
| 384 | 464.94364 | 431.74377 |
| 392 | 468.91733 | 430.40304 |
| 408 | 477.17620 | 430.57368 |
| 424 | 493.25313 | 443.03179 |
| 448 | 522.70843 | 460.14418 |

These sums include the intersection with the earlier edge bounds. They
are local factors, not end-to-end margins. The local test completes in
10.73s. It follows the full regional run serially.

An orthogonal possibility is to multiply the history potential by kappa^a,
where a is the preceding input count. A true transition acquires exactly
kappa^(l-a), so these factors telescope to a terminal correction. Changing
this potential can affect the slack of a coarse-cell bound without changing
the code or its actual path weights. It is not itself a distance gain.
The identity passes966 exact transition checks and1250 exact path checks.

## A bounded diagnostic of the remaining loss

`quartet_capped_cycle.py` checks a fixed count path inside the relaxed
operator. The first four counts are408. It then uses44 counts407 and4048
counts408, placing one407 in each of44 regions. The remaining212 regions
contain only408. Thus the path respects the regional payment counts.
After three transient transitions, it remains in a single coarse age cell.
For each edge it uses min(m_e,g), which is feasible in the row-cap program.

The largest tested FLOAT path expression is-27937.69012503454, at incoming
a384..447,n1,n2 in0..63. The check takes27.86s. This is a path in the relaxed
bound, not a codeword. Its negative value does NOT prove the whole bound
small. It also does not establish an intrinsic positive obstruction for
this coarse transfer or justify declaring its potential optimization
exhausted. The counter retains much more path entropy than this diagnostic.

The current potential is proposed from the unrestricted one-update
operator. That operator includes arbitrarily long full-input runs, although
the regional quota rules them out here. Before enlarging the grid, test a
potential proposed from the quota-weighted sixteen-update operator. The
joint moment edge caps offer a cheaper linear candidate for that proposal;
every proposed vector must still be checked against the actual upper
operator. Do not infer a certificate from proposal convergence. Preserve
the passing narrow-band checkpoint and the failed complete-count controls.
