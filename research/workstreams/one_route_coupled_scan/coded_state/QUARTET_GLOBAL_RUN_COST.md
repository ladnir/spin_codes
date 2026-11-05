# Global state-run closure cost

This checkpoint preserves the construction and the full 45%/40-bit goal.
It strengthens the grouped field count and records a corrected diagnosis
of the remaining global bound. The full certificate remains open.
All setup, field, footprint, rank and terminal conventions are those of
[the rank-gap checkpoint](QUARTET_GLOBAL_RANK_GAP.md).

## Correct the diagnostic before changing the proof

Finite differences at log-marker step 1/1024 crossed between different
dominant contributions. Their separate sensitivities must not be combined
into one trajectory. In particular, the earlier inference of thousands
of forced births per region in the 512-row class was incorrect. It is
not evidence of a quota obstruction. The coarse probe remains explicitly
labelled as a negative control.

Central differences at step 2^-20 are consistent across the coverage and
eligible/outside markers. Selected local sensitivities of the upper
envelope are approximately:

| h | covered F | rank R | available L_eff | covered in H | forced births |
|---:|---:|---:|---:|---:|---:|
| 512 | 129165 | 156362 | 107788 | 38339 | 6653 |
| 2048 | 1132763 | 903775 | 891035 | 491132 | 71426 |
| 4096 | 1631733 | 1083854 | 1281531 | 1072044 | 137634 |

These are sensitivities of a relaxed upper bound, not statistics of a
sampled code. They diagnose the missing lower eligible-degree constraint
at h=512 and upper rank/capacity constraints at h=4096. The h=2048 class
does not violate those aggregate constraints. Wider quota tilts and
joint degree markers do not close it.

Replacing eligible rows by actual active rows is not a useful repair by
itself. For a birth-free row with v>=129 covered positions, the shortened
RS space splits as

    q^(-128) = q^(-v) + (q^(v-128)-1) q^(-v).

Summing active-row subsets recovers the same rank moment. Exact nonzero
row enumeration saves less than a 2^-19 relative correction here. New
trajectory-dependent support constraints are needed for a material gain.

## Forced births start distinct availability runs

Recall A_t=V_t intersect the next fully covered quartets, and A_T=V_T.
Set A_0=empty. For a fixed coordinate, an availability run is a maximal
consecutive interval on which that coordinate belongs to A_t.

A canonical forced birth comes from a partial covered quartet, or a full
covered quartet with no old covered coordinate. In the first case,
A_(t-1) is empty throughout that quartet. In the second, the old covered
set is empty there. Thus every forced birth belongs to A_t\A_(t-1).
Each forced birth starts a distinct coordinate run. If f is their total
number, at least max(0,f-512) such runs end before the free terminal state.

This statement concerns canonical birth-forced positions. It does not
apply to every propagated nonzero-state label introduced below.

## A joint field-count bound

Fix 0<x<=1, and define

    H = [[1,qx],[x,1+qx]],
    lambda = 1+(q+1)x,    c=(q+1)x,
    e=min(1,lambda/(q+1)).

H upper-bounds one coordinate's state/output count while that coordinate
is allowed. A run of length ell that starts with a required birth has
terminal-free count at most c lambda^(ell-1). This follows from the
positive vector v=(1,1+1/q), since Hv<=lambda v and the birth contributes
qx in the nonzero state.

For a closed run, the last state must return to zero. Its count is

    qx e_1^T H^(ell-1) (1,x)^T = x (H^ell)_(0,1)
      <= x lambda^ell.

Combining this with the terminal-free bound gives
c lambda^(ell-1) e. Unforced runs contribute at most lambda^ell.
Multiplying over all coordinates and runs therefore yields

    state sum <= c^f lambda^(L-f) e^max(0,f-512)
              <= e^(-512) (c e)^f lambda^(L-f).       (1)

This is a single run-level argument; it does not multiply unrelated
upper bounds or charge a reset twice. It retains the free terminal.

After the fractional moment theta, the existing birth-only transfer can
use

    m_eff = (c e)^theta,
    b = lambda^theta,
    global factor e^(-512 theta).                    (2)

When lambda<q+1, (2) simplifies to m_eff=z b, where z=x^theta. The bulk
factor b is unchanged. Apply (1) to canonical footprints before adding
the positive extra annotations of the 129-state relaxation. Those extra
annotations need not themselves have the availability-run property.

Two independent audits support the argument. Exact rational checks cover
480 forced-start closed runs and 480 terminal runs, for q up to 2^32.

## A separate accepted-map propagation control

The grouped rank bound also omits a deterministic support property of
accepted matrices. In a full quartet with exactly one old covered
coordinate, any outside-H position forces that incoming state to be
nonzero. Nonzero columns of both A and A+delta I then force the outside-H
successor state AND output to be nonzero.

The implementation retains a guaranteed number of such singleton quartets.
Among k nonfresh full quartets and at most U old covered positions, at
least j=max(0,2k-U) are singletons. Select j actual singletons before
inspecting new routing labels. This auxiliary selection has total weight
one. Only afterward enlarge the old-partition count to a positive upper
kernel.

The reduced basis has 258 states for coverage-bin endpoints
64,128,192,255,512. In its outgoing polynomial, each selected position
contributes p y_O + out y_F, with out=(1-p)/a. The corresponding next-step
coefficient factor is p(1+bz)+out*mz. The propagated F position may use m
because both its state and output are forced nonzero; it is not merely
an active state. The terminal and zero operators remain exact.

Two audits and 34 exact factorization comparisons pass. This coarse
control still has positive middle-class bounds, and its weakened
partition count can outweigh its gain. It is a separate control, not
a claimed pointwise improvement over the 129-state transfer. Do not apply
the birth-only closure bonus to every merged propagated F position.

## Full-history numerical results

All numbers are FLOAT upper-bound diagnostics, not certified margins.
The closure-cost substitution, with nine quota tilts, gives:

| h | theta | log2 z | log2 eta | log2 upper bound |
|---:|---:|---:|---:|---:|
| 512 | 1/8 | -1/4 | 0 | +50908 |
| 2048 | 1/16 | -7/16 | 0 | +117834 |
| 4096 | 1/16 | -5/4 | 1/2 | +289604 |

The last two improve on the previous selected bounds of +136434 and
+411663. A bounded three-tilt grid finds +283412 at h=4096,
theta=3/32,log2 z=-3/2. The grid is not an optimality claim or a full
parameter cover. Every tested middle class still has a positive bound.

All previous checks also pass. Numerical jobs were serial, single-threaded,
file-free on Peach, below the 175/180-second and 2 GiB guards. No benchmark,
construction, production, paper, commit, or push changes were made.

## Remaining proof task

The useful next split retains joint rank and available-position ranges.
For fixed R,L, optimize the original grouped field bound before choosing
a common fractional moment:

    min(1,(gamma/q)q^(-R) inf_(0<x<=1)
          x^(-W)[1+(q+1)x]^L).

For L<=W the optimizer is x=1. For L>W it is
x=W/((q+1)(L-W)), capped at one. Count finite (R,L,F) boxes and apply their
appropriate field factors. Use pure geometric bounds only where the
conditional field factor saturates at one. Retain every earlier class
exclusion and the invariant-plateau control.

Alternatively, track more availability-run closures than the forced-birth
count alone guarantees. The relation D_A>=D_V-(F-L) is useless here:
A_t is a subset of V_t intersect V_(t+1), so D_V<=F-L already. A valid
closure count must retain overlap with the next availability holes.
Do not multiply (1) by another closure penalty spending the same runs.

The full fixed-construction 45%/40-bit certificate remains the goal.
