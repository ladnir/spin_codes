# A uniform bound using inactive coverage and support deaths

2026-10-27. The construction is unchanged. This lemma complements
[the low-activity certificate](QUARTET_LOW_ACTIVITY.md). Both statements
apply to the actual once-sampled accepted field-quartet setup over J=GF(2^32).

## Statement

Retain K=2^20, N=2^21, M=8192 outer RS[256,128] rows, and 4096 states of
width512. For a nonzero J-message, let n be its number of active outer rows.
Let A=sum_t wt(Q_t), and let D count the coordinates leaving state support
across all updates. Q_0=0, so the first update contributes no deaths.
Let C_T count covered input positions in inactive outer rows. A live
incoming quartet covers all four input positions. An empty incoming quartet
covers exactly its nonzero next-state positions. Define R=C_T+D.
The terminal state is unrestricted. Put W=floor(0.45N)=943718 and

    h(n)=max(0, W-128n+91522).

The probability that any nonzero J-message satisfies both wt(Y)<=W and
R>=h(n) is less than2^-63. The union includes all n, all row sets, and all
labeled support histories. It holds conditionally on any fixed route and
quartet partitions and hence under their sampled distributions.

Since C_T>=0, D>=h(n) is a sufficient condition for this exclusion.
In particular h(n)=0 for n>=8088. The lemma therefore covers every low-output
history in those row-count classes, not only a turnover-restricted subset.
It does not cover all nonzero messages of the full code. As with the
low-activity lemma, the restricted statement is explicitly over J.

At r=2100 inactive rows, n=6092 and h(n)=255464. Combined with the preceding
lemma, the remaining event has

    A>=694489,  C_T+D<=255463,  wt(Y)<=943718.

The two excluded events together have probability less than2^-44+2^-63,
which is below2^-43. These are global allowances, not costs to repeat for
each row count.

## Count supports together with selected zero outputs

Fix an active-row set, state supports, and a choice Z of selected zero
output coordinates in consecutive-support intersections. Write z=|Z|.
The established pin-and-shortening argument bounds the probability of a
compatible message by

    2 q^(A-128n-C_T-z),  q=2^32,

where C_T>=0 counts covered inactive-row positions. The factor two bounds
the global acceptance density for the actual field-multiplication setup.
This is the same conditional count used by the earlier quartet certificates.

The union of consecutive state supports has total size A+D. Thus a
low-output history admits a selection of

    z=max(0,A+D-W)

zero outputs in the intersections. If A+D>=W, the exponent above is
W-128n-C_T-D. Otherwise it is A-128n-C_T, which is no greater than that
expression. Retaining C_T gives the common bound

    2 q^(W-128n-R).

To count candidate support histories and selections, follow one fixed state
coordinate through time. When its incoming support bit is zero, its next
bit has two possibilities and cannot carry a selected intersection zero.
When its incoming bit is one, a next zero has one choice; a next one has
two choices according to whether its output is selected. The counting
matrix, indexed by incoming and outgoing support bits, is

    T = [[1,1],[1,2]].

Counting every selected subset, rather than only the required cardinality,
is an upper bound. With h=(1,13/8) and lambda=21/8, direct rational checks
give T h<=lambda h and (1,1)<=h. Starting from support bit zero, with a free
terminal bit, there are at most lambda^4096 histories for that coordinate.
The 512 coordinates therefore have at most lambda^N joint histories.
This is a count of all labels; it does not assume independent realized
state coordinates.

## Exact union calculation

The integer comparison

    21^1024 < 2^4498

proves lambda^1024<2^1426. Since N/1024=2048, the full history count is
less than2^2920448. Set

    L=(2920448+8192+64)/32=91522.

For R>=h(n), the exponent W-128n-R is at most -L, including the case where
h(n)=0. Sum the conditional bound over histories and active-row sets:

    Pr[excluded event]
      < 2 * 2^2920448 * q^(-L) * sum_(n=1)^M binom(M,n)
      < 2^(1+2920448+8192-32L)
      = 2^-63.

`quartet_high_turnover.py` replays the potential inequalities, acceptance
density, integer power comparison, threshold arithmetic, and final exponent
without floating-point arithmetic. The substantive pin-and-shortening
lemma remains part of the proof, not an assumption certified by this scalar
replay alone.

## Remaining-count interface

For n<8088, define a_min=b(n)n+1 from the low-activity certificate and
R_max=W-128n+91521. A remaining low-output history must satisfy A>=a_min
and R<=R_max. For eta>=1 and 0<zeta<=1,

    1[A>=a_min and R<=R_max]
      <= eta^(A-a_min) zeta^(R-R_max).

The implementation multiplies each local term of next weight v, deaths d,
and inactive coverage c by eta^v zeta^(c+d), then applies
eta^(-a_min) zeta^(-R_max) once globally. Omitting the covered marker gives
a valid weaker condition on D alone.
It does not impose a per-update activity floor. All states, turnovers,
activation times, and terminal weights remain in the transfer.

This marker inequality is valid but can be loose. The tests in
`quartet_remaining_activity.py` check its direction by exact short-path
enumeration and verify the modified activation/reset identities. Numerical
screens from that helper are not part of either exact certificate.

## Deterministic consequence of the routing quota

Let U count all uncovered input positions across the code. There are exactly
256r input positions from the r inactive rows. At most U can be uncovered,
so C_T>=256r-U for every route and history. The remaining event therefore
satisfies

    U-D >= 256r-R_max = 128r+K-W-91521.

At r=2100 this lower bound is282137. Because covered positions include every
nonzero next-state coordinate, U<=N-A. Consequently A+D<=1815015 in this
row-count class. This excludes the nearly full-state histories that dominate
some unconstrained numerical bounds. It does not assert that every remaining
history has low state weight at each update.

For gamma>=1, impose the first inequality by multiplying a local term by
gamma^(u-d), where u is its uncovered-position count, and dividing the final
sum by gamma^282137. In an empty incoming quartet u=4-v; in a live quartet
u=0. The marker acts on an exact necessary condition. The resulting FLOAT
transfer screen is still only diagnostic until verified with directed bounds.

## Bounded numerical checks

The following full 513-state screens use r=2100, ell=0.00377, and nine
Bernoulli proof parameters. Every turnover and terminal weight is retained.
Entries are FLOAT log2 upper bounds, not certificates. Positive entries
provide no nontrivial failure guarantee.

| theta | Activity eta | Marked constraint | Marker | log2 screen |
| ---: | ---: | --- | ---: | ---: |
| 0.175 | 1 | none (matched baseline) | 1 | +572899.41 |
| 0.175 | 1.02 | A lower bound only | 1 | +583912.00 |
| 0.175 | 1.1 | A lower bound only | 1 | +616789.56 |
| 0.175 | 1.5 | A lower bound only | 1 | +913374.68 |
| 0.175 | 1 | D upper bound | 0.9 | +597683.65 |
| 0.175 | 1.02 | D upper bound | 0.9 | +608800.82 |
| 0.175 | 1 | C_T+D upper bound | 0.9 | +592233.64 |
| 0.175 | 1.02 | C_T+D upper bound | 0.9 | +603602.36 |
| 0.175 | 1 | U-D lower bound | 1.1 | +622120.49 |
| 0.175 | 1.02 | U-D lower bound | 1.1 | +639531.92 |
| 0.175 | 1 | U-D lower bound | 4 | +3227199.07 |
| 0.5 | 1 | A lower bound only | 1 | +3238777.99 |
| 0.5 | 2 | A lower bound only | 1 | +4332932.84 |
| 0.5 | 8 | A lower bound only | 1 | +5935512.35 |
| 0.5 | 32 | A lower bound only | 1 | +7549276.22 |

These samples do not establish an optimum over markers. They show that the
new necessary conditions do not automatically improve the existing global
transfer bound. Penalizing one class can increase the charge for another.
They neither exhibit a bad codeword nor prove that the construction fails.

The exact replay checks eight short support/zero-selection automata.
The transfer helper checks nine rational marker inequalities, 54 boundary
cells, and 5460 exact coverage/quota inequalities. Numerical runs were serial,
file-free, and remote on Peach, with a90s watchdog and64MiB numeric-array
budget. The longest recorded run took about53s. The local20GiB disk guard
remains intact. An oversized command failed before launch; loading only the
required source functions into the remote process resolved that limit.

Next retain the joint remaining counts inside the routing calculation,
rather than rely on a single unconstrained positive weighting. The low-
activity and high-score lemmas are reusable exclusions. All other row
counts, all remaining histories, and a directed final union are still
required for the full45%/40-bit certificate.
