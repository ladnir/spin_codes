# All pair types and verified count regions

Follow-up: the [location-aware count](LAYERED_OUTER_LOCATIONS.md) closes the
large-empty-pair profile left open below, over every active-outer-row count.
This note preserves the earlier intermediate bounds and coverage.

This extends the [joint component count](LAYERED_COMPONENT_COUNT.md) without
changing the outer, route, or inner. All 40 failures from the previous selected
mixed-history screen now pass exact rational bounds. The local count also
handles 00 pairs. Verified neighborhoods and one gapless interval are recorded
below; they are not full count-space coverage or a whole-code certificate.

## Local extension

Keep the three temporal states: inactive, active-pinned, and active-unpinned.
The component update and its terminal convention are unchanged. Let z,r,f,e
count 00 pairs, 01 pairs, full pairs, and exceptional full pairs. If P counts
all noninitial pairs, there are P-z-r-f type 10 pairs. For each actual pair,
use these fixed-input tags:

| Pair / branch | First coordinate | Second coordinate |
| --- | --- | --- |
| 00 | gamma=0 | gamma=0 |
| 10 | gamma=1 | free |
| 01 | free | free |
| 11 regular | free | free |
| 11 exceptional | gamma=0 | gamma=1 |

A fixed coordinate has X=Q_next+gamma*Q. A free coordinate retains the
existing inverse distribution. The 00 branch has no random free coordinate;
the 01 branch has two and cannot be exceptional. Endpoints of a full pair
share the same branch choice.

The aggregate dimensions are now

    V=P-z+f,       D=V+r-2e,       K_R=max(0,RM-N/2).

For a fixed history with Z internal active-to-active zero-output edges,
a=V+b-Z. The consistency lemma and its interpolation remain valid:

    eta * q^min(a+d_bar-D-e, kappa+K_R-D-e),    eta<2.

Use separate positive markers for z,r,f,e,B and counted output weight.
B counts newly active coordinates at internal transitions. The final output
and terminal births are still omitted; the terminal state remains free.
For theta=m/32, the history weights are 2^[(32-m)*kappa-m*Z].

When z=r=0, retain the forced-active first coordinate from the earlier
count. Otherwise that restriction is removed, including from the initial
sum and the outgoing state potentials. A zero prescribed count may exclude
that pair type entirely. Merely using the old six incoming rows for a
positive r or z count would omit histories and invalidate the bound.

There are now nine incoming coordinate-state pairs, 100 internal transition
terms, and 13 terminal branch terms. Exact tests enumerate 2,207 width-four
histories. They include every pair type, compare the component transfer with
the independent path oracle, and check a complete rational weighted sum
against the positive-potential bound.

## Replaying the previous failures

The original screen is reproduced, not replaced by a new sample: 55,161
selected count profiles, 9,194 excluded by necessary wiring constraints, and
45,967 retained profiles, of which 40 failed the older bounds. Of those 40,
27 have r=0 and 13 have r>0. All have z=0.

The original rational weights close 24 failures. A second set closes all 13
r-positive failures. Three r-zero failures at R=176 need further lookahead.
Two-step reweighting still fails the hardest of these, with an exact upper
bound +32,091. Three-step label lookahead closes it at −47,363.

The retained bank has three small rational inputs: two with interpolation
m=27 and depth 2, and one with m=25 and depth 3. Replaying that bank closes
all 40 selected failures; the largest conditional log2 upper bound is −47,363.
The exact rational norms, not the floating optimizer's status, establish
these results. Existing good-route conditioning and its failure budget remain
unchanged.

## Verified boxes, not just selected points

Fix R and W=838860. The route-dimension envelope is nonincreasing in z,e
and nondecreasing in r,f. For a box of integer counts with lower and upper
corners lo and hi, a valid dimension cap is therefore

    d_box = d_bar(R, z_lo, r_hi, f_hi, e_lo).

Choose local rows using the types permitted anywhere in the box: include 00
if z_hi>0 and 01 if r_hi>0. The local positive norms then apply throughout
the box, including at zero-count endpoints. Outside d_bar, the logarithm
of the bound is affine in z,r,f,e,B. Maximize each coefficient at the
appropriate endpoint, using integer upper logarithms. All counts are
nonnegative, so these coefficient upper bounds remain valid even when a
coefficient is negative.

This verifies 40 boxes around the retained failures. Each varying count has
radius 128, clipped at zero; z remains zero, and r remains zero when the
center has r=0. Other r counts vary. The worst box's conditional per-class
log2 upper bound is −36,212. Union over every integer count in that box gives
−36,187. An independent direct evaluation checks all 424 distinct box corners
against their box majorants. The monotonic/affine argument, not corner
sampling, certifies every interior count.

These boxes need not be disjoint. They cover neighborhoods of the selected
failures, not the entire no-00 domain. The verifier makes no completeness
claim about gaps between boxes or other suffix lengths.

## First verified interval with empty pairs

Retain R=192,r=24568,f=589632,e=368520,B=49136. The 01-positive marker set,
with a smaller 00 marker, closes every integer z from 0 through 8192.
The source partitions this interval into 33 gapless segments, verifies each
segment's monotonic/affine majorant, and includes the union over 8,193 counts.
The conditional log2 failure upper bound for the interval is −10,577.

A larger-empty-pair profile remains open. At z=49136 with the other counts
above, a bounded reweighting of the one-update operator gives floating
upper bound about +586,317 and exact rational upper bound +594,094. This is a
failed bound, not a low-weight codeword or a proof of poor distance.
The forced-active-coordinate restriction is no longer available; the count
must handle correlations between empty pairs under successive pairings.

### Two actual updates and forced outer zeros

We also composed two actual updates on four coordinate paths. The first
pairing is (0,1),(2,3); the second is (0,2),(1,3). The local operator therefore
retains correlations that separate worst-row bounds discard. A product
potential bounds its 81 incoming states. Four distinct consecutive label
bits give equal block populations; the internal transition count is even,
so no unmatched internal update is omitted. The original initial and
terminal factors remain in place.

An independent exact 81-by-81 composition agrees with the contracted local
operator. For the z=49136 profile, the new rational bound is +556,423. It
improves the one-update bound, but does not close the class. A first bounded
proposal expired; the retained fixed-interpolation proposal completed and
its norm was recomputed before exact replay. An optimizer iterate is never
accepted as a proof of a local inequality.

There is a further deterministic outer constraint. Apart from the terminal
batch, each 00 coordinate whose next state is inactive has X=0. At least

    J_0=max(0,2z-B-b)

input coordinates are therefore zero. The subtraction of b allows every
last-batch 00 coordinate to escape this count. For a shortened RS row with
g<=R remaining allowed positions,

    max(0,g-128) <= (R-128)*g/R.

Thus the outer space compatible with these forced zeros has dimension at most

    K_eff=floor((R-128)*(RM-J_0)/R).

These coordinates lie outside the inverse free set. Replacing K_R by K_eff
in the consistency branch is therefore valid; the shortening dimension on
the free set is unchanged. No generic-rank or route-independence assumption
is needed. For this profile, J_0=48624 and K_eff=508080, reducing the exact
two-update bound to +491,591. All 24,768 scalar inequalities for 129<=R<=256 and
0<=g<=R pass; the exhaustive history tests also check the forced-zero floor.

This does not identify a complete repair. At the current m=28 weights, even
the hypothetical full independent-rank gain from all 2z empty-pair coordinates
would leave a positive score of 163,335. That is an optimistic diagnostic,
not an available rank theorem. A further attempt must retain more joint
outer/history information or reweight the complete count, rather than assume
that the zero-count floor alone closes the gap.

## Reproduction and remaining coverage

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_general_components.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/layered_general_components.py --bank
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_component_boxes.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/layered_two_layer_components.py --verify

The marker bank is source, not a generated result dump. Optimizers only propose
finite rational inputs; replay uses exact fractions and integer-power log
bounds. Numerical jobs are bounded, serial, and write no data files.

Next extend verified coverage across the remaining count domain, especially
histories with substantial 00 populations. Retain actual forced-input-zero
placements jointly with the temporal components and outer constraints; first
test the strongest plausible gain before further subdivision. Two-update
composition alone was insufficient on the retained profile. Keep the GF32
sandwich as the whole-code control. No production or paper change, native
benchmark, or new near-half claim follows from these restricted 40% results.
