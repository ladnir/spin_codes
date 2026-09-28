# Track all-one columns in grouped BCH words

The one-column permutation preserves four-row column weights. A column of
weight four has only 32 possible feedback values, one per window. Columns
of weights one, two, and three have 128, 192, and 128 possibilities. Treating
every column as weight four can therefore exaggerate cancellation.

This refinement counts those columns explicitly. It leaves the encoder,
uniform setup distribution, K=2^20, IMT(128,19), and output threshold 209715
unchanged. It currently verifies one sixteen-group support class, not every
message with sixteen active groups or the full code.

## Count columns through a subcode

Fix four BCH words of rank h and union support u. Let J count coordinates
whose four bits are all one. Let H be the image of the even-parity coefficient
vectors under the linear map taking four coefficients to their combination
of these words. At a coordinate, all even combinations vanish exactly when
the column is zero or all one. Thus H has support v=u-J.

When J>0, evaluation at an all-one coordinate gives a nonzero functional on
the span of the four words. Its kernel is H, so dim(H)=h-1. Each of the four
words lies in the same nonzero coset of H. For a fixed H and extension coset,
the number of ordered spanning four-tuples is

    2^(h-1) product_{i=0}^{h-2}(8-2^i).

Indeed, choose the first word in that coset. Its differences with the other
three words must span H. For h=1, the tuple is simply (w,w,w,w), and J=u.
The ordinary weight-u spectrum cap handles that case directly.

For h>=2, divide the existing rank-(h-1) four-tuple support-CDF cap at v by
the number of spanning ordered four-tuples in an (h-1)-space. This bounds
the number of possible H. Let D(v) be the verified shortening-dimension cap.
For fixed H, extension cosets with J new support coordinates number at most

    binom(256-v,J) 2^D(v) / 2^(h-1).

Restriction outside the support of H has fibers of size at most 2^D(v).
Every coset has 2^(h-1) words in its fiber. Another bound is the number of
nonzero BCH words of weight at most floor(u-v/2): every such extension coset
contains a word at most its average weight u-v/2. Use the smaller bound.

Summing ranks gives an upper count C_u(J) for each J>0. At J=0 use the
existing total support cap T_u. Also impose sum_J A_u(J)<=T_u on the unknown
true counts. No differences of upper CDFs are interpreted as shell counts.

## Couple these counts to the inner bound

Choose a positive witness rho. For each epoch input shape, multiply its
nine-coordinate transfer by rho to the number of weight-four columns in
that epoch, **before** taking entrywise maxima. After composing all epochs,
the resulting moment bounds include rho^J_total, where J_total is fixed by
the message and is unchanged by setup permutations.

Remove this factor when summing messages. For a fixed support u, maximize

    sum_J A_u(J) rho^-J

subject to the nonnegative shell caps C_u(J) and total cap T_u. Integer
greedy allocation to the largest weights rho^-J solves this relaxation
exactly. Group counts multiply. The existing support-probability witness
and exact without-replacement placement recurrence remain unchanged.

The local high-occupancy shape compression now optionally retains the
number of weight-four columns, in addition to total weight and the largest
mask multiplicity. Defaults retain the previous unweighted operators.

## Evidence

The generalized flag count passes 88 exhaustive small-code shell checks.
The compressed shapes are checked against all ordered tuples through eight
inputs. Additional transfer tests check that both weight multipliers are
applied to each shape before aggregation. These component checks do not
replace the BCH certificate dependencies or the outward moment replay.

For sixteen active groups, **each with support exactly 80**, lambda=.0064
and rho=1/2 give the 192-bit outward upper

    9.65599027075553140709800769900548033730052992859946383003e-17,

or **53.2013533906 bits**. This includes every rank, every allowed column
composition, and all binom(2048,16) group locations in that support class.
It does not cover other support vectors. The verifier uses only the earlier
authenticated BCH spectrum/shortening caps; it does not need the optional
dual-distance/OA refinement.

The independent 384-bit replay passed with the same rational support witness
773484531/10^9 for each group. Its upper begins
9.6559902707555314070980076990054803373005299285994636935006869e-17,
agreeing with the lower-precision result to the displayed margin.

The preceding pair-conditioned numerical bound for this class had log2
contribution +57.94285. Thus this refinement removes about 111 bits of
slack at that point. Other tested equal-support cases remain open: the
refined binary64 log2 bounds are -26.69 at u=84, +27.06 at u=96, and +238.29
at u=128. These screens include the optional OA count caps and are not
outward certificates.

The earlier total-input-weight experiment, `occupancy_weighted.py`, passed
30 exact small-code count/allocation checks but did not improve these
difficult points over the unweighted witness. A coarser statistic alone
was insufficient; the all-one-column constraint captures a more specific
cancellation mechanism.

Reproduction:

    python -B research/workstreams/permutation_locality/test_occupancy_tail.py
    python -B research/workstreams/permutation_locality/occupancy_weighted.py
    python -B research/workstreams/permutation_locality/occupancy_allones.py
    python -B research/workstreams/permutation_locality/occupancy_allones.py --tilts .0055 .0064 .007 .008 --penalties .125 .25 .375 .5 .625 .75 --supports 80 84 96 128
    python -B research/workstreams/permutation_locality/occupancy_allones_verify.py
    python -B research/workstreams/permutation_locality/occupancy_allones_verify.py --precision 384

The follow-up combined this constraint with fresh-state and exhaustive
window averages; see [WINDOW_AVERAGES.md](WINDOW_AVERAGES.md). Improved
interval counting, joint witnesses, and mature-tail bounds extended full coverage through thirteen groups in
[CDF_COVER.md](CDF_COVER.md). Full support coverage beyond that occupancy,
larger occupancies, and the final 5-6% latency reduction needed for 2x are
still open. Production defaults and the paper remain unchanged.
