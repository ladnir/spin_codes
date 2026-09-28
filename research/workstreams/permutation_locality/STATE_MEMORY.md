# Retain feedback history to bound cancellation

This note analyzes the unchanged one-column distribution from
[ONE_COLUMN.md](ONE_COLUMN.md), at K=2^20 with IMT(128,19). The expansion
and feedback maps are the actual selected production maps. Setup is ideal
uniform, not the heuristic permutation bank. The target bad event is an
output of weight at most 209715.

The preceding envelope treated every nonzero state as arbitrary after an
active epoch. That loses useful information: a single feedback input can
have an atom of 1/32, but the sum of two independent inputs has nonzero
atoms at most 1/512. The zero atom can still be 1/32. We retain zero
separately rather than applying the stronger bound to it.

## Condition on placement before using independence

Fix the message and every group's permuted column vectors. Within each
region, also condition on the assignment of the active labeled groups to
epochs, but not their window positions. The group's placement permutation
is uniform. Conditional on these epoch assignments, each epoch has an
independent uniform injection of its assigned groups into its 32 windows.
The lane shuffles are independent. Thus an epoch's input distribution is
independent of its incoming state, although multiple groups within that
epoch occupy distinct windows.

Every transition below bounds all fixed tuples of nonzero column weights.
Taking entrywise maxima removes those weights and the group labels from
the calculation. The exact occupancy polynomial in MULTI_GROUP.md then
averages the epoch assignments. No independent-window approximation is
made for groups sharing an epoch.

For epoch input x and entering state q, the output is x+Eq and the next
state is tau(q)+Bx. Each independently sampled transvection tau satisfies, for
fixed q!=0, the exact mixture: retain q with probability 1/2, otherwise
sample a uniform nonzero state. The output weight is evaluated before
this state update. Zero is fixed by the transvection.

## A nine-coordinate invariant

Let z=exp(-lambda), with lambda>0. The measure on the state is weighted by
z to the total weight emitted so far. Put m=2^19-1. The expansion classes
have weights v in {48,56,64,72,80} and sizes N_v.

For a in {1,2,3,4}, let p_a be the feedback distribution obtained from a
uniform window and a uniform weight-a mask in that window. The exact
single-window census constructs these distributions. They have disjoint
supports: each window map is injective, and distinct window images
intersect only at zero.

Decompose the weighted state measure into a zero component, a fresh
component, a mature component, and five expansion-class components. The
nine coordinates bound them as follows:

- Z bounds the mass at zero.
- F bounds the fresh component pointwise by sum_a beta_a p_a for some
  nonnegative beta_a with sum_a beta_a<=F.
- M bounds the total mass of the mature nonzero component.
- C bounds the mature component's mass at every individual nonzero state.
- L_v bounds the class-v component pointwise by L_v/N_v within that class,
  and by zero outside it.

The coefficients beta_a need not be stored. Bounds uniform over a suffice
for every mixture. M and C constrain the same component; they are not two
pieces of probability mass. The terminal functional is

    Z + F + M + sum_v L_v.

In particular, it excludes C. All coordinates and transfers are positive,
so componentwise upper bounds compose through epochs and regions.

## Preserve the invariant through an epoch

Fix the nonzero column weights in an active epoch. Their sum W is the
input weight because their windows are distinct. Write p for its feedback
distribution, p0=p(0), and alpha for an upper bound on its largest nonzero
atom. For every nonzero incoming state, a common moment bound is

    f = z^max(0,48-W).

For expansion class v, the pointwise bound is
tau_v=z^max(0,v-W). Let f_v bound its averaged moment and c_v bound its
weighted cancellation probability before the transvection coin. The
existing local enumerators give these quantities for one and two active
groups. Triangle bounds give them for three and four groups.

From zero, one input gives fresh mass z^W. Two or more inputs give zero
mass z^W p0, mature mass z^W(1-p0), and mature density at most z^W alpha.
The same-epoch pair and triple/quadruple censuses provide p0. Their maximum
atom bounds provide alpha.

For fresh mass, let c_F bound weighted lazy cancellation, and let d_F
bound the nonzero density after lazy addition, before the factor 1/2.
One active input uses exact independent feedback convolutions:

    d_F <= f max_a max_{y!=0} (p_a*p_b)(y).

Here b is the new input's column weight, and convolution is XOR addition.
The weighted c_F is counted exactly; cancellation is possible only when
the old and new single-input shapes agree. For two or more new inputs,
c_F and d_F are each at most f alpha. These follow by convolving a
probability distribution with p. No uniformity of the feedback is assumed.

The outgoing bounds from unit fresh mass are

    Z: c_F/2 + f/(2m),    M: f/2,    C: d_F/2,
    L_v: f N_v/(2m).

They hold for the fresh mixture by linearity and maximization over its
constituent shapes. The surviving lazy part becomes mature.

For the mature component, convolution cannot increase its density cap.
The incoming component has no mass at zero, so lazy cancellation is at
most f C(1-p0). Its outgoing bounds are

    Z: f C(1-p0)/2 + f M/(2m),
    M: f M/2,    C: f C/2,    L_v: f M N_v/(2m).

The measure can be counted twice by these upper bounds; this is harmless.
In particular, the nonzero-mass estimate need not subtract an upper bound
on the mass sent to zero.

For a unit class-v component, use the existing moment and cancellation
counts for Z and total lazy mass. The additional density bound follows
from its entering pointwise density 1/N_v and the pointwise tilt bound
tau_v:

    Z: c_v/2 + f_v/(2m),    M: f_v/2,
    C: tau_v/(2 N_v),    L_w: f_v N_w/(2m).

For all nonzero source components, the refreshed part has outgoing
pointwise density at most its moment divided by 2m. Feedback translation
can exclude a state but cannot increase this bound. This justifies both
the zero term and the expansion-class terms above.

An empty epoch leaves Z unchanged. It retains F with factor z^48/2,
and retains M and C with the same factor. Its refresh terms are the
corresponding mass times z^48 N_v/(2m). For a class-v component, its lazy
part stays in class v with factor z^v/2, and its refresh part contributes
z^v N_w/(2m) to each class w. Empty input cannot cancel a nonzero state.

These rules preserve all parts of the invariant. `occupancy_memory.py`
constructs them with outward Arb arithmetic and takes entrywise maxima
over input-weight tuples. Coefficient extraction uses the nine-coordinate
matrices without changing the placement distribution.

## Verification scope

`test_occupancy_memory.py` checks 98 exact rational density, convolution,
and fresh-mixture inequalities. It also compares 48 explicit two-step
calculations using the actual maps against the moment, zero, mature-mass,
and mature-density bounds. The latter use binary64 tolerances as regression
tests, not as outward certificates. The terminal-functional tests ensure
that C is excluded and that the previous all-mass default is preserved.

`occupancy_memory_verify.py` partitions every support from 38 through 256
into intervals, with extra boundaries at the rank minima. It covers every
ordered support box by summing sorted boxes with their exact multiplicity.
It multiplies by the number of group locations and the authenticated outer
CDF cap at each upper endpoint. Using a CDF cap for an interval overcounts;
it does not assume an exact support enumerator.

Each box selects a positive Cauchy witness from a finite list. The numerical
optimizer only proposes its probabilities. The replay uses their rounded
rational values, raises the positive nine-coordinate region matrix to 256,
and applies the terminal functional above. Binomial mass is log-concave,
so its minimum over each support interval occurs at an endpoint. Outward
lower bounds on those masses give a valid conditional bound for the whole
box. The final sum uses outward upper bounds throughout.

## Outward replay results

Both 192-bit and 384-bit replays pass on the complete support range, using
width-eight intervals with the extra rank boundaries:

| Active groups | Sorted support boxes | Failure upper, all locations | Margin (bits) |
|---|---:|---:|---:|
| Exactly 3 | 5456 | <2.265934e-22 | 71.9023 |
| Exactly 4 | 46376 | <1.967830e-23 | 75.4277 |

The 384-bit totals begin
2.26593383020843923996058205129427823563051652070702213649e-22
and 1.96782900032125412589089723758054356279867923558015515105e-23.
The 192-bit outward totals agree to the reported digits. Every rank tuple
is included. The weaker initial four-group width-twelve screen gave only
28.95 bits; this was interval-counting slack, removed by the finer cover.

Combined with the preceding one- and two-group bounds, the failure
contribution from messages occupying at most four groups is below
3.181390e-15, or about 48.1592 bits. This is not the full minimum-distance
certificate: messages in five or more groups remain uncovered.

The subsequent [adaptive extension](ADAPTIVE_OCCUPANCY.md) closes five and
six groups. The table above records this state-memory stage's own results.

Reproduce from the repository root:

    python -B research/workstreams/permutation_locality/test_occupancy_memory.py
    python -B research/workstreams/permutation_locality/occupancy_memory_verify.py --groups 3 --step 8 --precision 192
    python -B research/workstreams/permutation_locality/occupancy_memory_verify.py --groups 4 --step 8 --precision 192

Repeat the last two commands with `--precision 384`. The verifier
authenticates the existing 163 BCH dependencies and checks the 67 exact
shortening witnesses before using the outer caps. It checks the current
maps and local feedback counts; it does not rerun every historical proof.

The next step needs a scalable treatment of larger occupancies. As a first
probe, `occupancy_potential.py` searches for a common positive potential h
with R_r h<=a b^r h. This would reduce the outer sum to a one-dimensional
CDF evaluation. The tested four-group parameter grid is far too loose: its
best binary64 log2 upper is about +607.25. No potential witness is certified,
and that shortcut is not used in either replay above.

No encoder changes or timings are part of this refinement. The candidate
remains at 1.88--1.90x against the matched production controls, short of 2x.
