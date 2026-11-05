# Global adaptive zero budgets

This note gives two all-history classes for the unchanged quartet code.
It bounds zeros created at live quartet coordinates by fixing the outer
message before revealing the random field maps. It retains every later
input activation, including activations after the state has become zero.
These classes do not by themselves prove the full 45%/40-bit distance
target.

Use q=2^32, N=2^21, M=8192 outer RS[256,128] rows, and the accepted-setup
density bound gamma<2. Fix a route and all quartet partitions. For an
actual history let H contain the rows with at least129 covered input
positions, and let h=|H|. Every nonzero input row lies in H: outside H,
the outer word is supported on at most128 positions and must be zero.

## Count fixed messages before adaptive zero events

Fix a candidate H and one outer message supported on its h rows. There
are at most q^(128h) such messages. The message and route fix the entire
input stream independently of the random quartet multipliers.

At a live quartet the incoming state is nonzero. Multiplication by a
fresh raw uniform GF(q^4) scalar makes its four next-state coordinates
jointly uniform over GF(q)^4 after adding the fixed input. Conditional
on all past choices, its four coordinate state-zero events are independent
and each has probability1/q.

A coordinate output is Y_new=Q_new+delta Q_old. Thus the event
Q_new=0 or Y_new=0 contains at most two possible values of Q_new. Fix k
prescribed coordinate positions, each required to be live and bad. The
probability that all these joint events occur is at most (2/q)^k. If any
prescribed position is not live, this event fails. Fresh-map conditioning
proves the bound successively across updates; it does not select positions
after observing zeros, or assume independence of the adaptively determined
state supports.

Let E count state zeros at live coordinates, and K_bad count live
coordinates satisfying either zero event, counting each coordinate once.
The elementary subset union gives

    Pr[E>=e | fixed message, route, partitions]
        <= binom(N,e) q^(-e),
    Pr[K_bad>=k | fixed message, route, partitions]
        <= binom(N,k) (2/q)^k.                         (1)

Taking the H/message union and paying the accepted density once gives

    Pr[an h-row history has E>=e]
        < 2 binom(M,h) q^(128h) binom(N,e) q^(-e),
    Pr[an h-row history has K_bad>=k]
        < 2 binom(M,h) q^(128h) binom(N,k) (2/q)^k.   (2)

No state-support, birth-mask, death-count or output-support union appears
in (2). The count covers the full time interval. In particular, replacing
extinction by an absorbing state would be an invalid shortcut here.

## Two uniform coarse classes

For integers 1<=k<=n, binom(n,k)<=(3n/k)^k. Since

    3M<2^15,
    3N/256=24576<2^15,
    3N/272<2^15,

the first bound in (2), with e=256h, is less than

    2^(1-241h) h^(-257h).

The second, with k=272h, is less than

    2^(1-241h) h^(-273h).

If a threshold exceeds N, its event is empty. Summing over all h gives,
separately for each class,

    Pr[there is a history with E>=256h] < 2^-239,
    Pr[there is a history with K_bad>=272h] < 2^-239.  (3)

These statements hold uniformly over the route and partitions before
their averaging. They use the outer dimension128h, not a count of
possible state trajectories.

## Tighter exact thresholds and profile consequences

[quartet_adaptive_zero_budget.py](quartet_adaptive_zero_budget.py) reuses
the directed rational log routines in
[quartet_selector_capacity.py](quartet_selector_capacity.py). For each h
it verifies a threshold with a 2^-92 family budget using the exponent

    4096h+log2 binom(M,h)+log2 binom(N,e)-c e+1,       (4)

where c=32 for state zeros and c=31 for joint bad coordinates. The first
binomial uses its exact integer and a directed logarithm. For e<=N/2,
the second is bounded by N H2(e/N), again with directed rational logs.
The function N H2(e/N)-c e is strictly decreasing for e>=1: its
derivative is log2((N-e)/e)-c<21-c. The directed entropy interval is
less than one bit wide, so this monotonicity also holds for the reported
upper exponent. Bisection finds its least qualifying integer on that
half interval. If none qualifies, the separate bound binom(N,e)<=2^N
provides a certified threshold in the upper half. That branch does not
claim an entropy-optimal threshold.

Every returned threshold is checked directly with scaled integers.
Summing one 2^-92 threshold over all8192 values of h costs less than
2^-79; using both zero budgets costs less than2^-78. The serial replay
passed 2079 entropy comparisons, six quartet zero-union checks, 13056
state/output overlap checks, and ten threshold-row checks. The full
18-h threshold table also passed. Examples are:

| h | State-zero threshold e_h | Joint bad-coordinate threshold k_h |
| --- | ---: | ---: |
| 1 | 241 | 254 |
| 64 | 11406 | 11891 |
| 1536 | 229434 | 237627 |

For actual support statistics, E=F-Ls exactly. Moreover,

    K_bad=E+Z_survivor,

where Z_survivor counts zero outputs whose old and new states are both
nonzero. A bad output of weight at most W therefore satisfies

    K_bad >= E+max(0,Ls+D-W) = max(F-Ls,F+D-W).        (5)

For a profile with a proved lower death bound D_min, compare F-Ls with
the state-zero threshold e_h and max(F-Ls,F+D_min-W) with the joint
threshold k_h. Outside these classes, the death interval also obeys

    D <= W-F+k_h-1.

The remaining small-h problem is now a bounded-zero growth problem with
all RS-forced restarts present. For example, the certified h=1 remainder
has at most240 state zeros and253 joint bad coordinates, less than the
512 state zeros required to erase a fully active state. This observation
alone does not prove that every message reaches full activity early
enough; the growth and sparse fringe argument remains to be supplied.
