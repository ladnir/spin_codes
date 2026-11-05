# Episode geometry for the three-shear proof

Continuation: [UNIFORM_PAIR_DEFECTS.md](UNIFORM_PAIR_DEFECTS.md) now closes
all histories with at most four defective incoming batches at 46%, including
its route failure event. The zero-defect argument below is preserved.

The construction is fixed as in [the pin identity](UNIFORM_PAIR_PINS.md):
RS[256,128], one regional route, state512, rejected-uniform SL2 pairs over
GF(2^32), and Q_next=delta Q+Y. Setup is sampled once for every message.
No encoder or distribution change is made in this continuation.

The previously difficult no-birth/no-one-sided profile closes once its
temporal geometry is retained. More generally, the entire zero-defect
family below has a restricted 46% certificate for every route. This is
not a whole-code certificate: histories with births or one-sided pairs
remain outside that family.

## The zero-defect family

Fix a suffix of R regions, with zero incoming state and zero earlier
inputs/outputs. Each region has M=8192 coordinates. There are L=16R
batches of b=512 coordinates. Write Q_1,...,Q_(L-1) for the incoming
states after the initial update. Q_L is a free terminal state.

Require:

- No coordinate becomes active between these incoming states.
- Each incoming pair is either 00 or 11 under its current matching.

The first state Q_1 is unrestricted. The terminal update may create new
activity; its input/output must therefore be counted separately.
Let V be the total number of active coordinates in these incoming states.
It is even. Every live pair supplies two random inverse-input coordinates,
so their total number is exactly D=V.

Supports decrease over time. The backward-closure lemma from the previous
note implies that a proper nonempty support cannot persist through nine
distinct matching phases. Consequently the incoming states consist of:

    full prefix -> at most eight proper nonempty states -> zero suffix.

If a is the number of full states, then a<=floor(V/b). Therefore all
nonzero outputs, except possibly the terminal output, lie in the first

    p = min(L, floor(V/b)+9)

batches. The same prefix contains every random inverse-input coordinate.
The output support has at most U=min(RM,b(p+1)) allowed positions. The
extra batch allows the free terminal update even when the incoming state
there is zero.

## An inverse count that does not need random routing

Fix any route and any nonzero output word in this family. Its reconstructed
states determine the set F of random inverse-input coordinates. F lies in
at most r=ceil(p/16) regions. Every RS row restricted to those regions has
shortening dimension at most max(0,r-128). Thus

    d(F) <= min(V/2, 8192*max(0,r-128)).

The V/2 cap is the usual rate-half shortening cap, applied row by row.
For the fixed values on coordinates outside F, the number of compatible
outer words is either zero or at most q^d(F), where q=2^32. The pair-image
point bound therefore gives probability at most 2*q^(d(F)-V). The factor
two covers all pair-density losses together, not one factor per update.

Sum this bound over output words of weight at most W on the U allowed
positions. When q*(U-W+1)>=W and W<U, a convenient exact majorant is

    volume <= (W+1)*binom(U,W)*q^W.

Indeed, binom(U,j)*q^j increases through j=W under that condition, and
(q-1)^j<=q^j. Otherwise use volume<=q^U. The verifier uses integer
binomials and integer upper bounds on their logarithms.

There is no union over state masks in addition to this output count.
The inverse recurrence reconstructs the states uniquely from the output.
No route independence assumption or good-route conditioning is needed.

For the old failing profile R=256,V=1048320, the bound gives:

| Threshold | Log2 upper bound for this whole profile |
| --- | ---: |
| 40% | -5672925 |
| 46% | -1975691 |
| 49% | -227102 |
| 49.1% | -171302 |
| 49.5% | +48862, unresolved by this bound |

Here p=2056, U=1053184, r=129, and d(F)<=8192. All amplitudes and support
placements with the stated counts are included. A positive upper bound
does not establish a low-weight codeword or a distance limitation.

## Cover the whole zero-defect family at 46%

If the prefix plus the terminal batch touches at most128 regions, the
outer word must be zero. This follows from the RS distance129 and holds
for every route. V=0 is excluded for nonzero messages by the same argument.

In every remaining case,

    d(F)-V <= -K+24b,             K=1048576.

For r>128, the ceiling inequality 16r<=floor(V/b)+24 gives this directly
from d(F)<=8192(r-128). For r<=128, d(F)=0. A nonexcluded prefix must
touch at least128 regions before possibly adding the terminal region;
this implies floor(V/b)>=2024 and hence V>=K-24b.

There are at most 2^N possible output supports and at most q^W amplitudes
per support. Union over at most256 aligned suffix choices. The entire
family therefore has bad-word contribution at most

    2^[9+N+32*(W-K+24b)].

At W=floor(.46N), its log2 upper bound is **-194007**. This already sums
over all V: the bound is uniform for each output word before the output
sum. It does not need an additional activity-count union.

An independent interval replay checks all R=129,...,256 and all even
V=2,...,RM-b. Within each interval of fixed floor(V/b), U and the regional
cap are constant and d(F)-V is nonincreasing. Checking its left endpoint
covers the interval; the right endpoint is also checked as a control.
There are394240 intervals, of which259072 are excluded by the outer.
The largest retained coarse class bound is-217567. The simpler analytic
family bound above includes all suffix choices and is the reported result.

## Extend the geometry with defects

Call an incoming-state batch defective if it contains a one-sided pair or
if its following internal transition has any coordinate birth. At the last
incoming state, count only one-sided pairs; leave the terminal transition
free. Let h be the number of defective batches.

Removing these h indices leaves at most h+1 clean runs. Within each clean
run the supports decrease and are closed under the current matching.
The same backward argument bounds its partly active states by eight.
Including defective indices gives the deterministic bound

    number of partly active incoming states <= 9h+8.

For width2^d with d cycling matching bits, the corresponding bound is
d*h+d-1. An exhaustive support-path dynamic program verifies it on every
width-four, length20 support path, using4140 final dynamic states and
411280 transitions. This supplements the general clean-run argument.

The useful new parameter is the placement of defect batches, not just the
total number of births or one-sided pairs. Many births in one batch break
only one clean run. Conversely, spreading small defects across every batch
makes the bound vacuous. Such histories still need the weighted pin count
and the outer constraints; a small-defect theorem alone cannot cover them.

## Decision and next test

Continue the proof, rather than redesign now. Retaining interval geometry
changed the old large positive bound into a broad restricted certificate
without more arithmetic or randomness. That is concrete proof progress.

Next combine the defect positions with the actual regional support count.
Start with one defective batch, retaining its full/partial/zero intervals
and the free terminal update. Then separate clustered defects from widely
distributed defects, where the pin-transfer method should be tested.
Do not multiply the uniform-feedback certificate by an unexplained loss:
its forward reset bound need not hold for a partially active pair state.

The goal remains active. Stop and report a proposed redesign only if a
remaining family identifies a missing rank or probability property that
the fixed construction cannot supply cheaply, or if the bounded proof
experiments consistently cease to give useful information.

## Reproduction and preservation

Run serially:

    python -B diagnose_pair_episodes.py --test
    python -B diagnose_pair_episodes.py
    python -B diagnose_pair_episodes.py --cover
    python -B diagnose_pair_episodes.py --defect-test

Checks include972 exact Hamming-volume inequalities,274 exhaustive proper
support starts, the complete interval replay, and the defect-path test.
The code uses only the standard library and writes no data files.

C: had only about0.5GiB free, below the20GiB floor. Numerical work used the
configured Peach connection, with62GiB free on /tmp, a256MiB process limit,
and30/60-second timeouts. Every run finished in under two seconds excluding
connection overhead. Only source and this compact note are kept locally.
No native benchmark, library/paper change, commit, push, or promotion.
