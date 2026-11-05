# Whole runs under exact regional routing quotas

2026-10-27. The construction is unchanged. The preceding
[joint count](QUARTET_JOINT_PRIOR.md) reduces the r=2100 distance class to
a sufficient geometric inequality. This note proves that inequality for
full-width bursts and extends the result deterministically to some partially
active histories. It does not cover arbitrary histories or other row counts.

## The geometric question

Write A for total state activity, D for support deaths, and C for covered
inactive input positions. Fix lambda=1/4 and tau=37/50, and define

    F(S)=lambda*A-C-tau*D

for a state-support history S. The existing directed field argument would
close the r=2100 class if F(S)<=18171 held for every history and row set,
apart from a sufficiently small setup exception. The initial state is empty;
the terminal state is unrestricted.

The failed numerical count alternated nearly full and nearly empty states.
First isolate full-width bursts: each state is either empty or all512
coordinates. Such a history is a binary sequence b_0,...,b_4096 with b_0=0.
Let a count its live states and d its 1->0 transitions. Then A=512a and
D=512d. An update covers all512 inputs exactly when b_(t-1) or b_t is one.
The quartet partitions are irrelevant to this special case.

## Exact quota moment and complete burst count

Fix an inactive-row set of size r=2100. In each region, its inputs occupy a
uniform r-subset of8192 positions, independently between regions. If h of
the region's16 updates are covered, put H=512h. For b>1, the exact moment is

    f_r(h)=E[b^(-C_region)]
          =sum_c binom(H,c)binom(8192-H,r-c)b^(-c)/binom(8192,r).   (1)

The sum ranges over feasible hypergeometric counts. It does not replace the
regional quota by independent labels or condition each update separately.

For each incoming/outgoing burst bit, enumerate all16-update binary paths,
retaining h until applying (1). Weight each path by

    b^(128*a - (9472/25)*d) f_r(h).

These four sums form a2-by-2 regional matrix G. Multiplying256 such matrices,
starting at bit zero and ending against(1,1), sums every burst pattern,
including arbitrary activation times, multiple runs, and a free terminal.
The moment for any fixed inactive-row set is at most (G^256 (1,1))[0].

For a threshold B, Markov's inequality and a union over all row sets give

    Pr_route[exists R,S_burst: F(S_burst)>B]
      <= binom(8192,2100)b^(-B)(G^256(1,1))[0].              (2)

This statement is uniform over all quartet partitions. The route is sampled
once. The quantifiers in(2) allow both the row set and burst pattern to be
chosen after that route is known.

## Directed finite result

Freeze log2(b)=125/64. Every live update contributes exactly2^250 and every
death contributes exactly2^-740. The helper upper-bounds b^-1 by an80-bit
dyadic number; integer128th-power comparisons verify its direction.
All hypergeometric recurrences, matrix entries, and vector updates are
rounded upward to bounded dyadic precision. The final logs use rational
intervals. No floating-point value is trusted.

The replay gives

    B=3702,
    log2 Pr_upper = -66035138609487/1099511627776 < -60.      (3)

Thus, except with probability at most2^-60, every full-width burst and every
inactive-row set of size2100 satisfies F<=3702. This is stronger than the
18171 sufficient threshold, but only for the stated geometric family.
It is not an end-to-end minimum-distance certificate.

Two other tested proof bases gave larger offsets27099 and53966. They are
not used. The frozen choice in(3) is independent of sampled setup; selecting
a proof parameter does not introduce a setup-dependent union cost.

## Deterministic extension to partially active states

Compare an arbitrary support history S with any full/empty template T.
Call a selected position in an empty template state an addition. Call a
missing position in a full template state a hole. Let E_plus and E_minus
count additions and holes. Let H_start count holes in the first state of
each full run. Let H_end count holes in the last state of each full run
that ends before the free terminal boundary.

The following coverage inequality holds for every route, row set, and
quartet partition:

    C(T)-C(S) <= E_minus + (4/5)H_start + (1/5)H_end.        (4)

To prove it, consider an inactive input covered by T but uncovered by S.
Every selected template neighbor of that input must be a hole. Charge one
unit equally among those neighbors. At an interior full-to-full update
there are five: four previous states and the current state. At a birth
there is one current neighbor; at a death there are four previous neighbors.
An interior hole receives at most one unit in total. A run-start hole may
receive another4/5, and a nonterminal run-end hole another1/5. Terminal holes
have no future input incidence, so the same upper bound remains valid.
Additions can only reduce the lost coverage and therefore do not weaken(4).

Because A(S)-A(T)=E_plus-E_minus, (4) implies

    F(S) <= F(T)+Delta(S,T),
    Delta = lambda E_plus +(1-lambda)E_minus +(4/5)H_start
          +(1/5)H_end - tau[D(S)-D(T)].                    (5)

In particular, extra support deaths reduce the correction. This is sharper
than treating every modified state position as an independent loss.

On the single good event in(3), (5) holds simultaneously for all supports
and all templates. Thus F(S)<=3702+min_T Delta(S,T). There is no extra union
over deviations or templates. Histories satisfying

    min_T Delta(S,T) <= 14469                             (6)

meet the sufficient routing threshold18171. Combining this fact with the
existing activity and joint field exclusions closes the r=2100 low-output
event restricted to(6), over the base field J. Unrestricted extension-field
distance still awaits the full base-field proof.

The minimum in(6) is a two-state dynamic program using only the weight
sequence and total deaths. A template-zero update costs lambda*w; a
template-one update costs(1-lambda)(512-w). A0->1 transition adds
(4/5)(512-w_current). A1->0 transition adds tau*512+(1/5)(512-w_previous).
Subtract tau*D(S) once globally and take a free terminal minimum.

Exhaustive C++ tests check(4)--(5) for1048576 support/template pairs. The
dynamic program agrees with all65536 independently enumerated support
cases. Python additionally checks1000 candidate templates.

## What remains in the failed numerical trajectories?

A rounded illustrative trace from the earlier geometry moment has4096
updates, A=1537024, and minimum possible D=67712. Its minimum template
charge is1072124/25, about42885, above14469. An additional1420798/37 deaths,
about38400, would suffice to bring that charge below the cutoff. These
numbers describe an illustrative weight sequence, not a codeword, sampled
encoder behavior, or a probability distribution on full trajectories.

The whole-run result excludes pure full-width bursts. The remaining
difficulty includes the partially active transition intervals and the
relation between their coverage and excess deaths.

## Keeping the reference under the same routing

We also tested a stronger use of(3). For a template T chosen from S and a
marker gamma>=1, the factor gamma^(3702-F(T)) is at least one on the good
event. Multiplying the geometric moment by this factor retains C(S) and
C(T) under the SAME route. It does not resample or independently bound them.

A hard threshold T_t=1[wt(S_t)>=256] fails numerically: the artificial count
can oscillate around the threshold and pay an excessive reference switch
factor. Its FLOAT offset screen is1477389.99, not a certificate.

### Averaging all thresholds

For each j in1..512, use T_t^(j)=1[wt(S_t)>=j]. Let u_t,v_t be the incoming
and outgoing state weights, and let L_t count inactive inputs at that update.
Averaging over j gives the exact identities

    average A(T^(j)) = A(S),
    average D(T^(j)) = D_min := sum_t max(u_t-v_t,0),
    average C(T^(j)) = sum_t L_t max(u_t,v_t)/512.

Hence (3) implies, simultaneously for every S and row set,

    F(S) <=3702+sum_t[L_t max(u_t,v_t)/512-C_t]
                  -tau[D(S)-D_min].                       (7)

This removes the threshold discontinuity and isolates the desired geometric
gain: coverage beyond a fractional full-width reference, together with
deaths beyond the minimum required by changing weights.

The first FLOAT moment test of(7), at log2(b)=8, uses exact positive-kernel
identities and upper chords for its log-convex label moment. Three chord
anchors give an offset screen329810.56, still far above18171. This does not
invalidate(7); it shows that this coarse moment calculation does not close.
The code labels these screens as untrusted and retains the separate2^-60
burst allowance. The derived full-width certificate remains directed.

## Validation and next step

The scalar replay checks64 exact hypergeometric cases and60 independently
enumerated regional matrix entries. Reference tests check50 exact label-tilt
identities,625 terminal-path identities,625 averaged-threshold identities,
and100 log-convexity cells. All numerical jobs are serial and remote, with
explicit time and memory limits. No raw data, build trees, or reports remain.

The next discriminating question is how much coverage and excess turnover
the partially active transition intervals must create. Use(7) to separate
that question from already-covered full-width switching. Avoid another
unrestricted search over scalar weights or treating a failed upper bound
as a low-distance witness. First compare exact partial-state moment cells
against the chord upper bounds, to distinguish numerical-envelope loss
from genuine loss in the combinatorial count.
