# Fit the boundary weights to the regional quota

The construction is unchanged. This diagnostic covers every input count and
every regional payment allocation at r=2100,c=111572. Other row and covered
counts, directed arithmetic, and whole-code closure remain open.

The preceding [capped regional calculation](QUARTET_CAPPED_REGIONAL.md) fits
its positive vector to unrestricted updates. The regional quota excludes
some trajectories favored by that proposal. Here we fit a second vector
to the actual sixteen-update regional operator.

## A linear upper operator with coupled age bounds

Intersect each original edge bound with the joint-moment edge bounds from
the preceding checkpoint. Use all25 marker pairs with base-two logarithms
in {-1,-1/2,0,1/2,1}. This gives a linear upper transfer M_l for each exact
next count l. It does not enforce the additional total-row cap inside
each transition. Thus it is not the same operator as the earlier nonlinear
cap calculation. Both are valid upper transfers in exact arithmetic.

A positive one-update vector v, a count marker z, and a maximum row ratio
rho normalize these transfers to P_l. Every row sum of sum_l P_l is at
most1 in exact arithmetic. Include the true initial row as an extra source
with no incoming edges. The initial normalization factor is v_start.

For a positive boundary vector w, the exact-count recurrence computes

    U_s(w) = [X^(6528+s)] (sum_l X^l P_l)^16 w.

The initial source is evaluated by the same recurrence. Define

    Q_s = binom(6528+s,436+s) / binom(8192,2100),
    R_y(w) = sum_(s=-436)^1664 (Q_s/Q_0) (y/z)^s U_s(w).

Let gamma=max_i R_y(w)_i/w_i over ordinary boundary cells only. The first
region starts from the true initial row, not an arbitrary boundary cell.
Let T_i be the free-terminal factor after removing the one-update potential.
The complete marked expression is bounded by

    binom(8192,2100) b^(-165556)
      v_start rho^4096 Q_0^256 y^44
      R_y(w)_start gamma^255 max_i(T_i/w_i).

All factors are positive. This inequality holds for any positive w and y.
The proposal alternates a damped regional update of w with a scalar marker
search. Neither convergence nor a claim of global optimality is needed.
The terminal factor is applied once, and boundary ages persist across regions.

## Test whether this transfer itself is too loose

A poor maximum ratio can result from poor weights. To distinguish that
from an intrinsically large relaxed coefficient, retain only regional
totals s=0 and s=-1. Set

    a_s = min_i [(Q_s/Q_0) U_s(w)_i/w_i],
    t = min_i T_i/w_i.

Fix the first region to s=0. Among the remaining255 regions, choose44 with
s=-1 and211 with s=0. Positivity gives the following LOWER diagnostic for
the coefficient of this relaxed linear transfer:

    binom(8192,2100) b^(-165556)
      v_start rho^4096 Q_0^256 z^44
      U_0(w)_start binom(255,44) a_0^211 a_(-1)^44 t.

This statement concerns the relaxed transfer, not actual codewords or their
failure probability. The two regional operators need not commute: their
componentwise lower bounds multiply in either order. The extra source
remains first, and the terminal factor remains last. A positive lower
diagnostic would rule out repairing this transfer by changing only the
payment marker or the boundary vector.

## Current run and checks

The exact toy test checks nine three-region upper bounds, their numerical
replays, and the distinct initial/free-terminal handling for non-eigenvectors.
The balanced-coefficient lower bound has a separate exact toy check.
The balanced-coefficient test passes nine exact three-region comparisons.
On the width64 incoming-count cells, the full numerical results are:

| Boundary fitting | FLOAT upper log2 | FLOAT balanced lower diagnostic |
| --- | ---: | ---: |
| Initial boundary weights | +419478.71616 | +393110.89070 |
| After seven damped updates | +407961.58605 | +404094.84018 |

The final payment marker has log2(y)=-4.18191468458. The complete run takes
169.72s, and the replay with the balanced diagnostic takes167.46s. The
positive balanced lower diagnostic gives strong numerical evidence that this
relaxed linear transfer cannot be repaired just by changing boundary weights
or the payment marker. It is not a directed impossibility certificate and
does not say that the actual code has low distance. This diagnostic
concerns the new linear joint-moment transfer, not a lower bound on the
older nonlinear total-cap operator.

The counter-floor contribution is over713bits below every final regional
output in this coarse run. The current source uses a more conservative
32-floor allowance to cover the stated floating row-sum tolerance, rather
than the16-floor allowance used in the coarse replay. Neither calculation
is a directed certificate.

The next controlled comparison narrows only the incoming input-count cells
from width64 to16, keeping both age-cohort buckets at width64. It has525
states rather than165, with the same7725 exact outgoing count labels.
This isolates the effect of forgetting a=l_previous inside a coarse cell.

That width16 comparison completes in323.66s. After three damped boundary
updates, its FLOAT upper expression is+382506.65795879357 and its balanced
lower diagnostic is+377771.8287581421. At the same number of updates, the
width64 values are+407962.18495 and+404042.29409. The reduced input-count
uncertainty helps, but leaves most of the gap. Its counter-floor bound is
over639bits below every regional output. Further count refinement alone
is not supported as the main next step by this comparison.

A separate diagnostic retains exact counts
{0,312,384,392,407,408,424,512} and width16 age buckets. It includes extreme
counts outside the passing392..424 band and has2502 states. It is implemented
in `quartet_quota_sparse.py`. Its restricted count set is printed explicitly;
it cannot certify all input counts or replace the full target.

That restricted diagnostic completes in111.39s. Its FLOAT upper expression
is+263756.41563485947 and its balanced lower diagnostic is+221544.01083661325.
Thus exact selected counts and width16 age buckets do not by themselves
remove the problem introduced by these tails. This does not invalidate the
earlier passing392..424 calculation, which includes no such extreme counts.

## Test the remaining old-age relaxation before adding another state index

The two retained ages leave n3 and n4 merged. Their incoming coefficient
factor is bounded using max(L3(z),L4(z))^(n3+n4), where

    alpha=3*lambda-tau,
    L3(z)=(1+h4*b^alpha*z)/h3,
    L4(z)=(1+h4*b^lambda*z)/h4.

The actual n3 count is determined by the preceding selected type2 positions;
it is not an independent choice at every update. A bound that maximizes
this split repeatedly may handle bursts poorly even when its steady-count
diagnostic is tight.

A targeted potential test sets h3=h4*b^(alpha-lambda), with the corresponding
change d2=h3. This equalizes the selected-coordinate slopes d3/h3 and d4/h4.
All actual path weights still telescope with their correct terminal factors.
Only proof weights change. The original log2(h4/h3) is0.5814903353544948;
the aligned value is0.7671875. The selected-coordinate gap before alignment
is about0.18570bits. This is a local relaxation gap, not a security margin.
The aligned run is a single diagnostic candidate, not an assumed optimum
or a substitute for retaining the actual third cohort when needed.

The aligned diagnostic completes in110.24s. Its restricted upper expression
is+274462.68244578515 and its balanced lower diagnostic is+249420.5231063863.
The ordinary potential's corresponding values are+263756.42 and+221544.01.
Thus this single targeted potential change does not repair the loss; it is
not promoted as an improvement. Neither restricted run covers every count.

Next retain the actual third-age count through a two-update calculation.
If v2 type2 positions are selected at the first update, the next n3 is
exactly v2. The second update must use that same value. Test a tail
transition before enlarging the full state grid. This addresses a concrete
missing dependency; it does not yet establish that this dependency is the
only obstruction or that the construction has the target distance.

Runs are serial, remote, one-thread, bounded by1200s and2GiB of conservative
array storage, and file-free. Positive coefficient floors are retained in
upper diagnostics. Lower diagnostics subtract a bound on the explicitly
added counter floor. Floating arithmetic still prevents a directed claim.

## Extend a successful weight choice across covered counts

For the current multipliers mu=1/5,nu=eta=0, the necessary threshold is

    B(c)=(6/5)c+713983/5-111127.

The count marker contributes y^(111616-c). Holding all proof weights fixed,
the log2 bound therefore has slope -log2(y)-24/5 as a function of c. If that
slope is negative, the sum over a covered-count tail is bounded by a
geometric series. This uses the same validated transfer, not a new setup
distribution. Other row counts still require their own valid constraints.

There is also a direct all-c choice: set log2(y)=-24/5. Then the coefficient
of c cancels. More explicitly, every remaining bad message supplies its
actual covered set P and satisfies

    H(P)-(6/5)|P| > 713983/5-111127.

Union and Markov therefore bound the entire r=2100 class by

    binom(8192,2100) b^[-(713983/5-111127)]
      E[sum_(P subset I) b^[H(P)-(6/5)|P|]].

The regional transfer already sums every permitted local paid count.
At the stated marker it evaluates this generating function directly, not
one coefficient. No extra union factor over c is then needed. Including
P outside the necessary covered-count range can weaken this bound, so it
still needs an actual evaluation; cancellation alone does not close it.
