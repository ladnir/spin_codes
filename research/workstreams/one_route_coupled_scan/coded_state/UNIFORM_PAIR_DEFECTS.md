# Births, one-sided pairs, and full-state resets

Continuation: [UNIFORM_PAIR_REGIONAL_RANK.md](UNIFORM_PAIR_REGIONAL_RANK.md)
closes the five-defect family using actual regional occupancies and joint
episode geometry. This four-defect checkpoint and its verifier are preserved.

The fixed three-shear model now has a restricted 46% certificate for all
histories with at most four defective incoming batches. The bound is below
2^-149, including the route failure event. Histories with five or more
defective batches remain open. This is not a whole-code distance certificate.

The [zero-defect argument](UNIFORM_PAIR_EPISODES.md) is preserved. This
continuation keeps RS[256,128], the regional route, state width 512, and
the reset unchanged. It adds no mixing or arithmetic to the encoder.

## What counts as a defect

Work over J=GF(2^32), with N=2^21, K=2^20, b=512, and L=4096 updates.
There are 256 regions, each containing 16 batches and one coordinate from
each of the M=8192 RS rows. Setup is sampled once for every message.
The pair distribution and its point bound are specified in
[UNIFORM_PAIR_PINS.md](UNIFORM_PAIR_PINS.md).

Index updates by t=0,...,L-1, with Q_0=0 and free terminal state Q_L.
The output reconstructs every state through Q_(t+1)=delta Q_t+Y_t,
without revealing any pair maps. Count defects among Q_1,...,Q_(L-1).
Index t is defective if its current matching has a one-sided pair or
if Q_(t+1) has a coordinate birth. At t=L-1 count only one-sided pairs;
births into the terminal state are free.

Let h bound the number of defective indices. This counts batches, not
individual births or pairs: arbitrarily many defects within one batch
still count once. Unlike the earlier suffix convention, a first activation
after Q_1 consumes a defect. This proof counts the whole global timeline.

Let a indicate Q_1!=0. Let c be the number of intervals of nonzero incoming
states, and T their combined length. There is no nonzero message with c=0:
then only the terminal input can be nonzero, contrary to the outer distance.
Write P for the number of proper, nonempty incoming states.

## Geometry that survives defects

Every new nonzero interval after the first requires a birth from zero.
Therefore c<=h+a. Removing the h defective indices leaves clean runs with
decreasing supports, each closed under the cycling matchings. Each clean
run contains at most eight proper nonempty states, by backward matching
closure. If a=0, the leading clean run is zero. There are at most h+a
clean runs that can contain proper states.

Among defective indices, at least c-a have zero incoming state: they
activate a new interval. These indices cannot themselves be proper states.
Consequently

    P <= 8(h+a) + h-(c-a) = 9(h+a)-c.

Use P_max=min(T,9(h+a)-c). The random inverse-input set F contains both
coordinates of every live pair. Every full incoming state supplies b
coordinates, and every other nonzero state supplies at least two. Hence

    |F| >= D_0 = b(T-P_max)+2P_max.

An interval of ell batches intersects at most floor((ell+30)/16) regions.
If it starts at Q_1, the bound is floor((ell+16)/16). Summing over intervals
gives the regional cap

    r_F <= min(256, floor((T+30c-14a)/16)).

All possible inputs and outputs lie in the live intervals, their preceding
activation batches, and possibly the final batch. Call that set I. It has
at most U=min(N,b(T+c+1)) coordinates and intersects at most

    r_I <= min(256, floor((T+31c-15a)/16)+1)

regions. If r_I<=128, the RS outer excludes every nonzero input.
These are deterministic facts about the actual route, not independence
approximations.

## Lower-bound parity-check rank before counting outputs

For an RS row with f coordinates in F, the shortened dimension is
(f-128)_+. Thus the parity-check rank on F is sum_rows min(f,128).
If F intersects at most r regions, this rank is at least

    ceil(min(1,128/r)*|F|).

This retains more information than separately bounding dimension by K.
Apply it first with |F|>=D_0 and r<=r_F.

Full incoming states give a second bound. A transition from nonfull to
full requires a birth, so full states form at most h+1 intervals. If
there are h+1, their first interval starts at Q_1. Otherwise there are
at most h arbitrary intervals. Both cases give boundary allowance
E=30h+16.

There are at least A_0=T-P_max full states. Take the earliest A_0 of them
if there are more. This subset remains a union of the required number
of intervals; truncation only shortens the last retained interval.
The subset argument is important: the rounded numerical rank envelope
need not be monotone as A_0 increases.

Their bA_0 coordinates intersect at most min(256,floor((A_0+E)/16))
regions. They also contain at least max(0,ceil((A_0-E)/16)) complete
regions. Each complete region contributes one coordinate in every row.
Apply the previous rank inequality and, separately, count these complete
regions up to 128. Let R_0 be the largest of the three rank lower bounds.

For any fixed route and output word in the class, inverse reconstruction
has probability at most 2*q^-R_0 of landing in the outer code. The total
factor two is the proved density inflation over all live pairs. Therefore
the inverse bad-word contribution for fixed intervals is at most

    2 * Vol_q(U,W) * q^-R_0,       W=floor(.46N).

No extra state-mask union is needed: the output determines the states.

## Forward count for the same interval placement

The allowed input space supported in I has dimension at most

    d_I = floor(U*(r_I-128)/r_I)     when r_I>128.

Fix one such message and the route before exposing the pair maps. On
each nonzero incoming pair, its output distribution is bounded pointwise
by eta times the uniform distribution on J^2. Replace these outputs by
uniform draws adaptively; the product of all density factors is below two.
Zero incoming pairs remain deterministic in this comparison experiment.

If the actual history has at least D_0 random coordinates and weight at
most W, its first D_0 random coordinates have weight at most W. Pad the
comparison experiment with independent draws if it stops earlier. This
defines the first D_0 draws without conditioning on survival. Consequently
the forward contribution is at most

    2 * q^d_I * Vol_q(D_0,W) * q^-D_0.

Take the smaller complete forward or inverse bound for each class.
The count over messages is outside the probability; the message need not
be chosen independently of setup in the minimum-distance event.

The verifier bounds Hamming volumes with integer arithmetic. It uses
binom(n,z)<=(3n/z)^z and the ambient q^n bound. Selected borderline cells
use exact integer binomials through the preserved episode helper.

For fixed (a,c), there are at most (L+1)^(2c) interval placements, including
all lengths. With c<=h+1, the conservative logarithmic placement union is

    26(h+1) + ceil(log2(2(h+1))).

Do not additionally sum over T: the endpoint count already includes it.
All values of T are checked to obtain a uniform per-placement bound.

At h<=3, this argument alone gives log2 family upper bound -161687 at
46%, for every route. It checks 28656 count cells, excluding 13884
deterministically. Neither a route allowance nor the forward bound is
needed for that result; the forward bound becomes useful below.

## Separate a full reset from the free terminal interval

Call a history a full-reset history if some full incoming state is followed
by a zero incoming state before Q_L. If no full state ever occurs, T<=P_max.
Otherwise, in the absence of a full reset, the last live interval must
continue through Q_(L-1). Every preceding interval is entirely proper.
Their combined length is at most P_max.

This removes the extra terminal batch and sharpens r_I. The verifier
enumerates the combined length of the early intervals using a small
integer dynamic program. The first interval has an aligned start if a=1;
the last has an aligned end. Intermediate intervals have arbitrary offsets.
It maximizes regional cover before using the forward/inverse counts.

For h<=4 with no full reset, the complete result is -3723 at 46%.
All 36839 count cells are covered; 17971 are excluded by the outer.
The inverse-only bound is positive, so this result uses both arguments.

## The remaining four-defect classes use an actual route bound

For unrestricted h<=4, only 34 count cells fail the coarse -80-bit
family target. Exact Hamming volumes close four. The remaining 30 have

    |I| <= 994816,    I touches <=130 regions,
    I has <=11 partially occupied regions.

The last bound follows from at most five live intervals, their activation
batches, and the terminal batch. We now control the RS dimension supported
in every such I, rather than replacing it by its deterministic maximum.

Fix I. If it touches at most 128 regions, no nonzero RS row fits. Otherwise
let r be 129 or 130, and j<=11 the number of partial regions. An eligible
row must hit at least k=129-(r-j) partial regions. Write p_l for the allowed
fraction in each partial region. Their sum is at most
994816/8192-(r-j).

For one row, its eligibility indicator is bounded by the sum of products
over its k-subsets of partial regions. For any s specified rows, expand
the product of those positive sums. In region l, requiring m different
rows to land in I has probability

    (|I_l|)_m/(8192)_m <= p_l^m.

The regions are independent. Summing the expansion bounds joint eligibility
by e_k(p)^s, where e_k is the elementary symmetric polynomial. Pairwise
averaging of its arguments cannot decrease e_k at fixed sum: the affected
term has coefficient e_(k-2)>=0. Thus

    e_k(p) <= binom(j,k)*(sum(p)/j)^k.

The verifier checks every r,j possibility with exact fractions. The largest
upper bound is at r=130,j=11,k=10:

    beta = 11*(39/176)^10
         = 8140406085191601 / 2592590903942070665216.

The probability of at least 32 eligible rows is at most
binom(8192,32)*beta^32. Union over every relevant I uses the same
134-bit interval bound above. The resulting route-failure upper bound
is 2^-151. Outside that event, at most 31 rows are eligible. Each has
shortened dimension at most two, so d_I<=62 simultaneously for every I.

No assumption of independent row eligibility is used. The falling-factorial
calculation handles the dependencies within a regional permutation.

## Pay for the full reset under that route event

Fix a compatible message, a route, and an interval placement. If a full
state later becomes zero, every coordinate has an active-to-zero transition
before the terminal state. Choose one such transition for each coordinate.
There are at most L^b choices.

At each chosen transition the incoming pair is nonzero. Requiring one
specified outgoing coordinate to be zero costs at most eta/q; requiring
both costs eta/q^2. Reveal maps in time order. Incoming nonzero requirements
can only discard histories; they do not renormalize the maps. The total
density inflation remains below two. Thus a fixed message has probability
at most 2*L^b*q^-b of the selected full-reset event after summing choices.

With d_I<=62, the per-placement forward bound is therefore

    1 + 32*62 + 512*12 - 32*512 = -8255

in log2 units. No independent-state or uniform-whole-state claim is used.
This bound applies to all 30 remaining cells on the good-route event.
The other cells use their direct forward/inverse bounds.

The complete four-defect union is:

| Event | Log2 upper bound |
| --- | ---: |
| Route exception | -151 |
| Full reset and good route | -293 |
| No full reset | -3723 |
| Sum of the three, rounded upward | -149 |

The certificate concerns the actual rejected-uniform SL2 setup. It covers
all placements and amplitudes in the stated J-valued history family.
It does not imply whole-code closure over J or its extension field.

## Checks, preservation, and the next obstruction

Reproduce serially with the new standard-library-only helper:

    python -B diagnose_pair_defects.py --test
    python -B diagnose_pair_defects.py --test-route-reset
    python -B diagnose_pair_defects.py --ladder
    python -B diagnose_pair_defects.py --defects 4 --no-full-reset
    python -B diagnose_pair_defects.py --route-reset

Exact checks cover 5934 support-DP states and 121552 transitions,
32767 interval geometries, 16384 terminal-live geometries, 20099 binomial
majorants, 2108 Hamming volumes, 6912 adaptive-draw events, 1500 small-route
tail bounds, and 256 adaptive coordinate-death bounds. These finite tests
check the implementation; the arguments above justify the general bounds.
The full four-defect replay takes about 3.4--5.3 seconds locally. It has
a 30-second cap; other checks have 20-second caps and small bounded arrays.
No files, raw dumps, or caches are generated. C: had over 24 GiB free.

Five-defect histories remain open. Even the no-full-reset branch has a
positive coarse bound near T=2058, c=6. That frontier needs the actual
regional rank on full incoming intervals, not just the tiny-input-space
argument used for the four-defect tail. Defect-rich histories still need
the pin-transfer bounds or another joint count; do not extrapolate this
small-defect certificate to all histories.

Continue the proof before redesigning. Both the all-route rank calculation
and the route/reset combination remove identified losses without changing
the code. Next isolate the five-defect no-full-reset frontier, retain its
partial-region occupancies, and test whether the resulting outer rank
closes it. Preserve all earlier proofs and the certified GF32 control.
No native timing, production/paper edit, commit, push, or promotion occurred.
