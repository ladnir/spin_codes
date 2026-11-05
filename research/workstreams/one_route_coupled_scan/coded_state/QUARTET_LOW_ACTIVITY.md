# A uniform bound for low total state activity

2026-10-27. The construction is unchanged. This lemma removes low-activity
histories from the remaining 45%-distance proof. It does not impose a
per-step weight band, a turnover bound, or a reset schedule.

## Statement

Use the actual accepted field-quartet construction over J=GF(2^32), with
K=2^20, N=2^21, M=8192 RS[256,128] outer rows, and 4096 updates of width512.
Setup is sampled once and fixed for all messages. Write n for the number
of nonzero outer rows of a nonzero J-message, and define its total state
activity by A=sum_(t=1)^4096 wt(Q_t). The initial state is zero; the terminal
state is unrestricted.

The probability that any such message has A<=b(n)*n is less than 2^-44,
where b(n) is given by this table:

| Active rows n | b(n) |
| --- | ---: |
| 1 | 84 |
| 2--3 | 86 |
| 4--7 | 88 |
| 8--15 | 90 |
| 16--31 | 92 |
| 32--63 | 95 |
| 64--127 | 97 |
| 128--255 | 100 |
| 256--511 | 102 |
| 512--1023 | 105 |
| 1024--2047 | 108 |
| 2048--4095 | 111 |
| 4096--8191 | 114 |
| 8192 | 118 |

The union covers every row set, every nonzero J-message, every labeled state
history in the stated activity range, and every first activation time.
No output-weight restriction is needed. The result holds conditionally on
any fixed route and quartet partitions, so averaging their setup choices
preserves it.

This is a base-field statement. A restricted activity property need not
transfer unchanged to an extension-field message. The eventual complete
minimum-distance statement over J does transfer to the extension encoder.

For example, r=2100 inactive rows means n=6092. The threshold is then
114*6092=694488, or 169.552734375 state coordinates per update on average.
The remaining proof may discard A<=694488 at this row count, while charging
the single global failure allowance above. It must still bound every
low-output-weight history with larger A.

## Conditional count without selected output zeros

Fix the route, quartet partitions, a set of n active outer rows, and all
state supports. A quartet is live if its incoming state is nonzero. Let L
count input coordinates in live quartets and B count nonzero next-state
coordinates in empty incoming quartets. These L+B positions form the
covered input set. Let C_T count its positions in inactive outer rows.

All input coordinates outside the covered set are zero. If an active RS
row has fewer than129 covered coordinates, no nonzero word can realize the
history. Otherwise, shortening each active RS row bounds the number of
outer messages by

    q^(L+B-C_T-128n),  q=2^32.

Fix one such message. In empty incoming quartets the next state equals the
input, so its B nonzero coordinates are already determined. The remaining
state amplitudes have at most q^(A-B) choices. For a fixed candidate state
trajectory, each raw live field-multiplication quartet has point probability
q^-4. Independence of the raw map samples gives q^-L for all live quartets.

The [actual field-quartet interface](INNER_MECHANISMS.md#exact-fixed-input-law-and-acceptance-bound)
has global acceptance density at most1024/989<2. Thus the probability that
any message realizes these supports under the accepted setup is at most

    2 q^(A-128n-C_T) <= 2 q^(A-128n).

This recovers the established pin-and-shortening count with no selected
output zeros. It counts candidate messages and amplitudes before exposing
the map values; no message is assumed independent of the realized setup.
The last inequality deliberately drops inactive-row coverage. It is useful
when A is small because its field exponent is already strongly negative.

## Sum all support histories

There are binom(N,A) labeled state-support histories with total activity A,
since the 4096 next states have N coordinates in total. For an integer
threshold a<=N/2, monotonicity of these binomial coefficients gives

    sum_(A=1)^a binom(N,A) q^A
      <= [q/(q-1)] binom(N,a) q^a.

After the row-set union, the contribution at row count n is therefore at most

    [2q/(q-1)] binom(M,n) binom(N,a) q^(a-128n).       (1)

Use binom(M,n)<=M^n and binom(N,a)<=(eN/a)^a. The elementary bound e<11/4
then turns (1), at a=bn, into

    [2q/(q-1)] [ M (11N/(4bn))^b q^(-(128-b)) ]^n.    (2)

The binomial bound follows from a!>=(a/e)^a. For example, integrating log x
from1 to a proves this factorial inequality. Bounding the exponential series
tail from1/4! geometrically proves e<11/4. All thresholds in the table have
bn<=N/2, as checked in the replay.

## Exact scalar certification

First consider b=84 for all n. Replace n in the denominator of (2) by one.
The resulting geometric ratio is

    R_0 = 2^33 (11/21)^84 < 4/(5*2^45).

The latter inequality is the integer comparison

    5*2^78*11^84 < 4*21^84.

Summing R_0^n for n>=1 already gives a uniform 84n threshold with failure
less than2^-44. The table strengthens the thresholds for larger n.

For each dyadic range with lower endpoint L=2^j, j=1,...,13, choose the
largest integer b in [84,127] satisfying

    R_j = M (11N/(4bL))^b q^(-(128-b))
          <= 2^(-ceil(60/L)).

These exact rational comparisons give the table. Within that range the
ratio in (2) is at most R_j. Since R_j<=1/2 and R_j^L<=2^-60, its sum over
n>=L is at most2^-59. The n=1 contribution uses R_0 directly. Consequently
the complete failure probability is at most

    [2q/(q-1)] [ R_0 + 13*2^-59 ] < 2^-44.

`quartet_low_activity.py` verifies these comparisons with integers and
rationals only. It also checks the simpler uniform84n bound and every
range endpoint. No floating-point logarithm, numerical optimizer, or
interpolation enters this certificate.

## Remaining obligation

This certificate removes low total state activity, not low output weight.
For each n, the remaining event has A>b(n)n and wt(Y)<=943718. Retain that
activity condition in the transfer or score-class sum. Adding an unconstrained
all-weight bound would pay again for the very histories eliminated here.
The remaining event, together with this single 2^-44 allowance, must still
fit the full 2^-40 failure target. No construction change is justified or
made by this lemma.

## Coverage-aware and trajectory-aware exclusions

The later [joint count](QUARTET_JOINT_PRIOR.md) also includes selected zero
outputs. It gives a sufficient weighted routing condition for the r=2100
distance class; that routing condition still needs an ensemble proof.

The basic field count retains more information than an activity bound. For
fixed geometry, row set, and exact support history S, it gives

    Pr_field[some message realizes S] <= 2 q^(A-128n-C_T).   (3)

The following base-field statements retain C_T and the shape of S. They
use the same once-sampled setup and free terminal state. No message is
assumed independent of the realized setup.

### Retaining inactive coverage

There are2^N state-support histories and fewer than2^8192 nonempty active-row
sets, including every cardinality. Thus the probability that any nonzero
message has A-C_T<=128n-65794 is less than

    2^(1+N+8192-32*65794) = 2^-63.                          (4)

This holds conditional on any fixed route and quartet partitions. There is
no extra factor8192: cardinalities are included in the power-set count.
At n=6092 the remaining histories have A-C_T>=713983. The helper
`quartet_rank_score.py` checks the exponent and64 exact power-set identities.

### A normalized prior on support trajectories

The power-set count treats a long nearly constant support as expensive as
an arbitrary changing one. To retain turnover, assign a probability mass
mu(S) as a counting device, not as a model of actual encoder states.

For each of the512 coordinates, begin at zero. Given beta,eta in(0,1),
generate independent binary chains with probabilities0->1 equal to beta and
1->0 equal to eta. Use the same parameters for all coordinates. Call this
normalized history prior mu_(beta,eta). Integrate both parameters independently
over[0,1] to obtain mu. Integrating normalized distributions gives sum_S mu(S)=1.

Let w be the final state weight. Transition counts are

    N00=N-A-D,    N01=D+w,
    N10=D,        N11=A-D-w.

Put n0=N-A+w and n1=A-w. Integrating x^k(1-x)^(m-k) gives
1/[(m+1)binom(m,k)], and hence the exact mass

    mu(S)=1/[(n0+1)binom(n0,D+w)(n1+1)binom(n1,D)].        (5)

Assign a normalized prior to active-row sets T as well:

    nu(T)=1/[(8193)binom(8192,|T|)].

Summing over all row sets, including the empty one, gives one. The empty
set contributes no nonzero-message event.

Fix the geometry. Sum(3) only over pairs(T,S) for which

    2 q^(A-128n-C_T) <= 2^-74 nu(T) mu(S).                 (6)

The sum is at most2^-74 by normalization. Except with probability at most
2^-74, no nonzero message realizes such a pair. Averaging geometry preserves
the bound. This covers every row count, every history, and the free terminal.

### Simultaneous linear consequences

For every beta,eta and every history,

    mu(S) >= mu_(beta,eta)(S)/[(n0+1)(n1+1)]
          >= mu_(beta,eta)(S)/2^44.                        (7)

The first inequality follows because each corresponding binomial probability
is at most one. The second uses(N+1)^2<2^44. Outside the single exceptional
event in(6), every pair of transition parameters gives the necessary condition

    32(A-128n-C_T) - log2 mu_(beta,eta)(S)
      + log2[(8193)binom(8192,n)] + 1+44+74 > 0.           (8)

There is no parameter union cost: (7) puts all these consequences on the
same good event, even if the proof selects parameters for a support class.

Write a=1-beta and c=1-eta. The negative log prior is L0+LA*A+LD*D+LT*w, where

    L0=-N log2(a),
    LA=log2(a)-log2(c),
    LD=log2(a)+log2(c)-log2(beta)-log2(eta),
    LT=log2(c)-log2(beta).

The helper rounds each coefficient upward using rational log intervals.
Since A,D,w are nonnegative, this weakens the necessary inequality safely.
A marker at least one raised to this upper score is at least one on the good
event. Retain the terminal coefficient; the final state is not zero.

Exact tests normalize three small shared-parameter chain families and eight
row priors, and check144 pointwise comparisons in(7). The default linear
consequence uses beta=1/9,eta=1/6. These are proof parameters, not changes
to the inner or its randomness.

### Current numerical scope

These exclusions are proved, but the output-distance calculation is open.
At r=2100 a stress test combined(8) with the unproved ensemble assumption
C_T+D>=0.24A-5000. Markers1.5 and4, respectively, give a FLOAT log2 screen
of+731706.98. The artificial flow shifts near weight371. This is not an
actual code-state distribution or a certificate. The final universal row
prior changes only the global constant from the earlier uniform row prior.

A separate test using(4), marker1.1, and cut marker2 gave+338842.54. These
specific trials are not optimizations or impossibility results. They do
not show that the proposed routing inequality would suffice. The best
unconditional full screen remains+271673.74.

All proved exclusions used so far, including the optional stronger single-
row activity checkpoint, sum to less than2^-43. A remaining all-history,
all-row-count distance contribution at most2^-41 would therefore suffice
for the original2^-40 target. `exclusion_budget()` checks this arithmetic;
it does not claim that the remaining contribution has been proved.
