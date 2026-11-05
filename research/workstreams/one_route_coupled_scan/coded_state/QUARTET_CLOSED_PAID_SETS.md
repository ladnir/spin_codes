# Count payments used by the maximizing support

The [crossing-burst obstruction](QUARTET_EXACT_MOMENT_OBSTRUCTION.md)
rules out the current moment over all paid sets. This reduction changes
the counted objects, not the encoder, its routing, or its randomness.

Fix the original inactive set I and the quartet partitions. For P subset I,
write U=I-complement union P for the allowed-input set. Retain the current
lambda=707/2560 and tau=381/512. Their ratio lies strictly between2 and3.

## The canonical maximizing support

Let S(U) include every coordinate in each nonterminal allowed run of length
at least3, and every coordinate in each terminal allowed run. Other runs
contribute no selected coordinates. This is the unique maximizing support
for H(U)=max(lambda*A-tau*D): a nonterminal run contributes either zero or
lambda*l-tau, and no integer l makes the latter zero at these weights.

Let V(S) be the set of input positions touched by S. A state coordinate
touches its own input and all four inputs of its quartet at the next update;
at the final update it touches only its own input. Feasibility implies
V(S(U)) subset U.

Define

    P'=I intersect V(S(U)), U'=I-complement union P'.

Then P' subset P. Moreover, S(U) remains feasible under U', since all its
touched inputs remain allowed. Monotonicity gives H(U')<=H(U), while this
feasible support gives the reverse inequality. Thus H(U')=H(U). Uniqueness
implies S(U')=S(U), so applying the reduction again leaves P' unchanged.
Call a paid set closed when it is fixed by this reduction.

## The bad event still implies a closed witness

For the current multipliers, B(c)=(6/5)c+d with
d=713983/5-111127. A remaining bad codeword supplies an actual paid set P
with H(P)>B(|P|). Its reduction satisfies

    H(P')=H(P)>B(|P|)>=B(|P'|).

It is therefore sufficient to union-bound over closed paid sets. No
preimage multiplicity is needed for this implication about existence.

The original lower bound|P|>=27476 belongs to codeword histories and must
NOT be imposed on P'. Reduction may make P' smaller. The upper bound
|P'|<=255463 survives. Likewise, field constraints on the original activity
must not be silently imposed on the canonical maximizing support.

## The regional factor depends on actual coverage

Fix the entire allowed mask U and all quartet partitions, then compute
S(U). In region j let a_j count allowed inputs and v_j count touched inputs
V(S(U)). Every touched input is allowed. The original inactive set contains
all8192-a_j disallowed inputs, so its number of paid positions must be

    c_j=a_j-(8192-r).

For a closed paid set, these additional inactive positions must be chosen
from the v_j touched inputs. Therefore the exact number of compatible
original inactive sets in that region is binom(v_j,c_j), not binom(a_j,c_j).
The corresponding factor in the setup expectation is

    binom(v_j,c_j)/binom(8192,r),

with value zero unless0<=c_j<=min(r,v_j). Multiply these factors across
regions and sum the allowed masks and actual quartet partitions, retaining
the desired total payment count. There is no additional paid-set binomial
factor. Coverage and the canonical support must cross regional boundaries;
neither is reset there.

This formula targets the unused-payment multiplicity identified by the
obstruction. It is not yet an efficiently evaluated transfer or a small
probability bound. In particular, v_j depends on both the maximizing runs
and the actual next-update quartets; it cannot be replaced by a marginal
expectation or an independently chosen count.

## Exact finite checks and next test

`quartet_closed_paid_sets.py` exhausts12544 pairs(I,P) at width4 over four
updates with two regions. Every reduction preserves the maximizing support
and its score, and every reduced paid set is fixed. Of these pairs,4592
are closed. Separately summing all65536 allowed masks with the new regional
factor reproduces the direct closed-set count, coefficientwise in paid
count and score (11 coefficients).

These are exact checks of the reduction and counting identity, not a
full-size distance certificate. The current field weights make the
canonical support unique; changing them requires revisiting ties.

Next retain enough local run information to determine the canonical support
and touched-input count jointly. For the current onset, selection of a
nonterminal coordinate is determined by whether it belongs to a fully
allowed length-three window. This needs two updates of lookahead, with
special treatment at the one actual terminal boundary. Test the resulting
coverage count against the crossing-burst family before a large all-history
calculation. Preserve the existing proof exclusions and failure allocation.

Continuation: the [closed-set transfer](QUARTET_CLOSED_TRANSFER.md) now
implements this dependency with a six-state run automaton and twelve
pending-input types. Its exact regional and random-partition identities
pass toy checks. A scalar full-width envelope still fails, and directed
burst tests require a smaller moment parameter. The remaining task is a
population- and quota-sensitive full-width bound, not another closure claim.
