# Retaining populations and run history

Checkpoint, 2026-10-30. The construction is unchanged. The full 45% distance,
40-bit failure-margin goal remains open. This note develops the
[closed-set transfer](QUARTET_CLOSED_TRANSFER.md); it is not a certificate.

The exact transfer tracks the allowed input set U and the canonical selected
state support S, with S contained in U. Its twelve pending types also retain
run ages. A useful bound must preserve enough history to prohibit repeated
births, premature deaths, and incompatible quartet requirements.

## Population-conditioned bounds

Initially retain just the counts of the three pending categories
(U,S)=(0,0),(1,0),(1,1). For a quartet population m, let F_m bound its
weighted transition sum, including a positive product potential. Given a
full population n of width512, the partition average is bounded by

    [z^n] (sum_m multinomial(4;m) F_m z^m)^128 / multinomial(512;n).

The coefficient notation selects precisely population n. The denominator
counts its assignments to labelled positions. Maximizing over run ages
within each quartet gives a valid upper relaxation, not the exact kernel.

`quartet_closed_population.py` extracts the population coefficient by a
positive array recurrence. It passed60 exact small coefficient checks and
60 independent-tilt replays. At the inherited product potential, one
width512 row decreases from the population-free bound668.44 to494.49 bits.
These are local FLOAT values, not failure margins.

`quartet_closed_count.py` instead extracts the original-inactive count i,
using positive evaluation for the source-population variables. It passed
480 exact small coefficient comparisons and480 tilt replays. Its outgoing
markers retain the allowed and selected populations. This original count
is distinct from the paid count; a region must have exactly2100 original
inactive positions in the current hard class.

`quartet_closed_cells.py` bounds all incoming populations in a specified
cell. For fixed positive source markers, the inverse multinomial factor
times the inverse marker monomial is log-convex in the population. Thus
its maximum over the cell occurs at a vertex. The implementation passed
585 comparisons with small integer coefficients. The full-size arithmetic
is FLOAT; no directed claim follows from these comparisons.

## The coarse regional shortcut fails

`quartet_closed_regional.py` constructs45 population cells, all513 original
inactive labels, and nine outgoing marker pairs. It uses no run-age
potential: an age-dependent potential can introduce substantial slack when
each transition independently maximizes over the forgotten ages.

A path that stays in one cell and uses twelve labels131 and four labels132
per region has the exact quota2100. For this fixed proposal, one such path
contributes a log2 value622862.10 to the relaxed final expression. The
initial factor and one terminal test are included. This rules out this
particular coarse transfer proposal. It does not lower-bound code failure
or the exact closed-set moment. The source cells and outgoing marker
bounds permit incompatible transitions.

The run took112.73 seconds, serially. Fifteen initial coefficient checks
and1215 outgoing-cell endpoint checks passed. Do not spend a large sweep
on this coarse matrix.

## Retain both endpoint populations

`quartet_closed_joint.py` retains incoming and outgoing category counts
jointly in the quartet polynomial. Positive evaluation bounds these
population coefficients; a positive recurrence extracts original
inactivity. Zero-population faces are removed exactly, not approximated
with extreme finite markers. Seventy-five comparisons with exact width8
coefficients pass. The routine also supports more than three categories.

There is a sharper identity when U=S at both endpoints. Write z and z'
for the incoming and outgoing zero counts. Every quartet containing a
selected coordinate forces all four next inputs to be allowed. Since the
outgoing allowed coordinates are all selected, outgoing zeros can occur
only in entirely unselected incoming quartets. Consequently z'<=z, and
there are no deaths.

For a fixed incoming zero set, let E be the number of entirely zero
quartets in a uniform partition. Its distribution is

    Pr(E=e) = binom(128,e)
      [x^(z-4e)] ((1+x)^4-x^4)^(128-e) / binom(512,z).

For each partition there are binom(4E,z') possible outgoing zero sets.
In particular, a transition with z'=z is impossible unless4 divides z.
When4 divides z, its routing factor is exactly

    binom(128,z/4) / binom(512,z).

With i original inactive positions, multiply this factor by

    binom(512-z,i-z) b^(lambda(512-z)-beta(i-z)),

where lambda=707/2560 and beta=6/5. This is the exact diagonal coefficient
for the unweighted pending-type kernel. `quartet_closed_diagonal.py`
checks275 exact occupancy coefficients. For example, the constant diagonal
population381 is impossible, whereas the coarse marker proposal assigned
it positive mass. These identities remove that apparent obstruction.

Mixed populations remain harder. At (allowed,selected)=(408,128), the
three-category conditional regional proposal is+1827.50 bits before the
final threshold factor. An inherited age potential reduces it to+1545.03,
but does not resolve the forgotten history.

## Exact seven-type interface

Several of the twelve pending types have identical future behavior. They
can be merged without maximizing over their histories. The resulting
types are:

| Type | Meaning |
| --- | --- |
| Z | Disallowed input; unselected state |
| V0 | Allowed, unselected; no preceding allowed run |
| V1 | Allowed, unselected; one preceding allowed position |
| B | Allowed, unselected; next quartet must not be full |
| L0 | Newly selected pending position |
| L1 | Next selected run stage |
| L2 | Selected run stage that can stop on the next exposure |

Here a full quartet means that all four next inputs belong to U. Any
L-type member forces its quartet to be full. A B member forbids fullness.
Thus a quartet containing both B and an L-type member has zero weight.

The old-to-new quotient is

    {0,3,6,15}->Z; 1->V0; 4->V1; {7,16}->B;
    2->L0; 11->L1; {14,17}->L2.

Conditioned on fullness, V0 advances to V1 and V1 advances to B.
Conditioned on non-fullness, each V-type or B member resets its age and
chooses Z, V0, or L0 according to its next input and selection bits.
A Z member also chooses among these three types. These choices must
agree with the shared quartet-fullness condition.

An L0 member advances to L1, and L1 advances to L2. An L2 member either
remains L2 or becomes B, recording a death. All L-type transitions require
fullness. Touched-input and field-score factors remain unchanged.

Initialize all positions at Z. At the single terminal test, accept only
Z and the three L types; short selected terminal runs remain allowed.
No terminal test or age reset occurs at regional boundaries.

`quartet_closed_quotient.py` verifies equality of all1365 original local
rows with the210 quotient rows, including their integer multiplicities
and marked degrees. All twelve terminal types agree. The quotient also
reproduces the eleven exact two-region coefficients, with131 maximum
states in that small replay. Its quartet polynomial has2245 marked
coefficients. This is an exact reduction, not a smaller-state code.

## Constraints that survive across updates

Write l0,l1,l2 for the selected type counts and D for current deaths.
The exact local identities imply

    l1_next = l0,
    l2_next = l1+l2-D,
    B_next >= D,
    selected_next-selected = l0_next-D.

Only L2 members can die, so D<=l2. All680 nonzero quotient rows satisfy
these identities. Independently, the original local test checks446
forbidden mixed quartets and6355 instances of the exact type16 death label.

For a constant-population trajectory, l0=l1=D. Therefore
selected>=3*l0 and B>=l0. Two previously tested candidate populations
violate this condition. Their large negative numerical proposals came
from an infeasible moment fit; they are not useful passing examples.

Keeping all seven types substantially changes the diagnostic bounds.
The following are FLOAT bounds for a fixed-population path with the same
twelve131/four132 schedule. They are normalized by the regional routing
denominator, but omit the final threshold, row union, initialization,
terminal factor, and all other trajectories.

| Population (Z,V0,V1,B,L0,L1,L2) | Conditional regional log2 bound |
| --- | ---: |
| (104,156,60,64,40,40,48) | -83.47 |
| (104,156,60,64,32,32,64) | +331.49 |
| (128,140,52,64,40,40,48) | -414.85 |
| (104,64,24,64,48,48,160) | -100.51 |
| (104,176,72,32,16,16,96) | +322.29 |

These values use a coefficient upper bound, not exact full-size
coefficients. A positive value does not prove that the corresponding
exact moment fails. A negative value does not certify the whole code.
The small tests and selected fits took under5 seconds each, serially.
No files or numerical dumps were generated on the remote host.

## Next discriminating test

Use the exact seven-type polynomial to retain the age-flow constraints
across updates. First test the two positive profiles consistent with age
flow using sharper joint coefficient extraction. Do not optimize profiles
that violate these exact constraints.
Then bound transitions between populations and compose the exact regional
quota, with one final terminal acceptance. Fixed-population checks alone
cannot replace that coverage. The other row counts, directed arithmetic,
and the final failure-budget sum also remain outstanding.
