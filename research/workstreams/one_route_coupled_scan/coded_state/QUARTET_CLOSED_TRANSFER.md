# A local transfer for closed paid sets

The [closed-set reduction](QUARTET_CLOSED_PAID_SETS.md) requires the actual
touched-input count, not an independently chosen coverage bound. This note
gives an exact local transfer for that count. It preserves regional quotas
and the one free terminal state. The encoder is unchanged.

The transfer identity is valid at full width. Only its small instances have
been evaluated exactly; a useful full-width bound remains open.

## Recognize the canonical support with a finite automaton

Let z_t indicate that a coordinate is allowed as a state position. At the
current field weights, the canonical support selects every nonterminal
allowed run of length at least3, and every terminal allowed run. Instead
of waiting for two future positions, guess whether a run is selected and
validate that choice when it ends. There is exactly one accepting choice.

Use six automaton states:

| State | Meaning after the consumed state position |
| --- | --- |
| 0 | The position is disallowed |
| 1,2 | An unselected allowed run has length1 or2 |
| 3,4,5 | A selected allowed run has length1,2, or at least3 |

For z=0, selection must be zero. States0,1,2,5 may move to0; states3,4
cannot end because their selected nonterminal run is too short. For z=1,
an unselected run follows0->1->2 and cannot continue further. A selected
run follows0->3->4->5->5. A run cannot change its selection midway.
At the actual terminal boundary, accept0,3,4,5. Reject1,2: an allowed run
ending there must be selected, even if its length is only1 or2.

Thus nondeterministic run selection does not introduce multiplicity. Every
allowed bit stream has exactly one accepted support. The toy helper checks
all2046 nonempty streams of length at most10 against direct run selection.

## Carry a pending input and selection

Write U_t for the allowed-input bit and S_t for the selected-state bit.
The state-position indicator is

    z_(t-1) = U_(t-1) * G_t,

where G_t indicates that all four current inputs in the coordinate's
quartet are allowed. Thus an update at time t can finalize z_(t-1).
Carry the automaton state after z_(t-2), together with the pending pair
(U_(t-1),S_(t-1)). Selection implies U=1, so the pair has three choices:
(0,0),(1,0),(1,1).

Only twelve of the eighteen apparent types can accept a continuation:
states1,2 permit only pending S=0; states3,4 require pending(U,S)=(1,1);
states0,5 permit all three pairs. These restrictions follow from the
transition rules, not a heuristic truncation.

At an update, choose the current(U_t,S_t), compute the quartet's G_t,
and advance the old automaton using z_(t-1) and pending S_(t-1). The outgoing
type contains that updated automaton state and the current pending pair.
The initial type is(0,0,0), representing an empty dummy state. At the end,
finalize z_T=U_T and apply terminal acceptance once. No automaton state is
reset at a regional boundary.

## The quartet factor charges exactly the used payments

Fix the four incoming types and one valid current assignment. Let u be the
number of disallowed current inputs. If any incoming S_(t-1) is selected,
all four inputs are touched; otherwise exactly the current selected S_t
positions are touched. Let v be this touched count. Validity implies that
every touched input is allowed.

Introduce x for an original inactive input and y for a paid input. Every
disallowed input must be originally inactive and unpaid. An allowed input
that is not touched must be originally active. A touched input may be
originally active or paid. The exact factor is therefore

    x^u (1+x*y)^v.

If a current coordinates are selected and d incoming selected coordinates
die, multiply by b^(lambda*a-tau*d). Also mark each outgoing type. Summing
valid current assignments gives the quartet polynomial K_m, where m is
the incoming twelve-type histogram of total size4.

The code enumerates81 possible current(U,S) assignments, rejects invalid
ones, and combines identical outcomes. There are1365 incoming histograms;
919 have a nonzero kernel. Their marked moment expansion has20723 terms.
These are small local objects, not an enumeration of full-width histories.

Expanding (1+x*y)^v and extracting the original-inactive count across a
region reproduces binom(total touched,c). This is the closed-set regional
factor from the preceding note, without storing a separate coverage counter.

## Average the actual random partition with its finite population

Let n be the incoming twelve-type histogram of total size w=512. Introduce
formal incoming markers z_1,...,z_12 and outgoing type markers inside K_m.
Define

    Q(z) = sum_(|m|=4) multinomial(4;m) * z^m * K_m.

The exact whole-update kernel, including the original-inactive, paid,
outgoing-type, and moment markers, is

    [z^n] Q(z)^(w/4) / multinomial(w;n).

To see this, expose the quartets as ordered slots in labelled groups.
Conditional on n, every slot assignment of the incoming types is equally
likely. There are multinomial(w;n) such assignments. The coefficient counts
them with their local factors. Forgetting group labels and within-group
orders gives the same uniform quartet partition; those orderings have
constant multiplicity.

Apply this kernel4096 times. In each group of16 updates, extract exactly
r original-inactive inputs and divide by binom(8192,r). Carry the pending
types across regions; only the regional counter resets. A y coefficient
retains total paid count, or y=b^(-6/5) sums the relaxed all-count moment.
The latter is an upper relaxation if the surviving upper bound on paid
count is omitted. Its threshold is d=713983/5-111127. No lower paid-count
bound may be imposed after closure.

## What the tests establish

`quartet_closed_transfer.py` reproduces all11 paid-count/score coefficients
of the exhaustive width4, two-region calculation. Its maximum dynamic-
program size is131 states. The paid-count-zero total is binom(8,2)^2, as
required because the empty paid set is always closed.

`quartet_closed_partition.py` compares labelled enumeration of all35
width8 partitions with composition over70 ordered group splits. Four
incoming populations give126 exact coefficient comparisons, including a
population where selected and unselected run constraints forbid some
partitions. Another45 coefficients check the one-update free terminal
against binom(8,r)binom(r,c). This checks the population normalization and
terminal condition independently of the two-region calculation.

The runs take0.29s and2.39s, respectively, serially and remotely with90s
watchdogs, less than256MiB, and no generated files. These checks do not
turn the full-width transfer identity into an evaluated distance bound.

## Cheap scalar envelope: valid but not useful

A first screen chooses a positive product potential over the twelve types,
then bounds each quartet's weighted row ratio by the same scalar gamma.
For any realized partition, a whole update is bounded by gamma^128.
Evaluating x positively rather than extracting each quota gives

    log2 upper = 524288 log2(gamma) - 537600 log2(x)
                 - d log2(b) + 512 log2(terminal ratio)
                 - 255 log2 binom(8192,2100).

The last coefficient includes the row-set union. The initial potential is
normalized to1. This formula covers all histories and all paid counts for
the one row class, but loses the actual populations and exact regional
coefficient constraints.

`quartet_closed_envelope.py` evaluates this valid formula after proposing
weights with a small, bounded optimizer. The FLOAT result is+912692.56 at
log2(b)=2.72208. The maximizing local row has four initial-zero types. The
run takes0.69s. This is a failed upper candidate, not a proof that every
such potential fails or that the code fails. The optimizer needs no SciPy;
an unavailable import was replaced by a NumPy-only proposal routine. The
final maximum is evaluated directly, without assuming convergence.

Do not tune this scalar envelope indefinitely. The next calculation must
use the finite-population coefficient and regional quota, rather than
maximizing the incoming quartet independently at every update.

## Retune the moment before the next full-width attempt

`quartet_closed_burst.py` tests the previous boundary-crossing bursts with
the corrected payment count. For a burst occupying v prescribed inputs in
a region, choose paid positions only inside those guaranteed touched inputs.
Every resulting paid set is closed: the canonical support includes the
guaranteed selected run and cannot touch an unpaid input. The event factor is

    binom(v,c_j) binom(8192-v,2100-c_j) / binom(8192,2100).

This gives a lower bound on the exact closed moment, since the canonical
support may gain additional score outside the prescribed burst.

The seven-update halves lose over101214bits of redundant payment count;
their former positive b=16 lower diagnostic becomes negative. However,
nine-update halves still give a positive directed lower expression at
b=16. For the selected field weights, that witness implies the necessary
restriction log2(b)<2.88903170 for a40-bit target. Thus the old b=16 choice
cannot simply be carried into the repaired proof.

These are directed lower diagnostics of the moment method, not probability
upper bounds or low-distance witnesses. Negative diagnostics do not certify
the corresponding family. A smaller moment parameter is only a candidate;
the population-sensitive transfer still has to be bounded and verified.

The full45%/40-bit objective remains open. Prior certified exclusions,
the frozen proof helpers, and the construction remain unchanged.
