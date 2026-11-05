# Sustained one-sided states in the layered inner

This note keeps the outer, regional route, pair wiring, and pair-key distribution
of [the layered candidate](LAYERED_INNER_PROOF.md) unchanged. It closes a new
restricted family over J=GF(2^32): incoming states may remain half populated
throughout a regional suffix. It does not close the whole code. The GF32
sandwich remains the certified whole-code control.

## Statement and scope

Use K=2^20, N=2^21, M=8192 outer RS[256,128] blocks, batch width b=512,
and q=2^32. Each region independently permutes the M block labels uniformly.
The incoming state is paired along index bit t mod 9 in batch t, with the
bit-zero coordinate first. Keys are independent across pair locations.
Setup is sampled once; the bad message may be chosen after setup.

For any suffix of R complete regions, consider outputs that are zero before
that suffix and whose incoming state in every batch after its first batch
has exactly one nonzero coordinate in each pair. The suffix's first incoming
state is zero. There is no constraint on the terminal state. The theorem
unions all such outputs, all their amplitudes and orientations, and all R.

Except with setup probability below 2^-52, there is no nonzero outer message
whose output belongs to this family and has weight at most

    W=838860=floor(0.40*N).

The class requires immediate half-population after the first batch and no
empty or full pairs thereafter. It permits arbitrary changes of orientation
and amplitude, not just a stationary half-support. Mixed histories, delayed
activation within the first region, and multiple zero-state episodes are
not covered by this statement. Nor does this restricted J-valued result
alone imply the analogous claim for arbitrary GF(2^128)-valued histories.

## Keep the orientation: the two one-sided cases differ

Write an incoming pair as (x,y), with

    p=x+a*y,             a uniform in J*,
    k=c/a,              c uniform in S,
    A(x,y)=(p,y+k*p),    |S| >= q-4.

For fixed output Y, the identity Q_(t+1)=delta*Q_t+Y_t determines all states
without using keys. The inverse input is X_t=Y_t+(A_t+delta*I)Q_t.
The local bounds therefore apply to fixed state values in this direction.

- Type 10 (x nonzero, y=0): one inverse-input coordinate is fixed, and the
  second has point mass at most 1/(q-1).
- Type 01 (x=0, y nonzero): p is always nonzero. Both inverse-input
  coordinates are free, with joint point mass at most 1/((q-1)|S|).

Thus it is conservative, but wasteful, to discard both coordinates of every
nonfull pair. In particular, type 01 has no exceptional branch.

Put u=R*M, L=u/b, and P=(u-b)/2. There are P one-sided pairs. If r are type
01, the inverse free set F contains every second coordinate (a head) in the
live batches, and r first coordinates (tails). Its size is D=P+r.
The affine outer-intersection bound is

    Pr[X(Y) is an outer word] <= eta*q^(d_pi(F)-P-r),
    eta <= (q/(q-4))^N < 2.

Here d_pi(F)=sum_i(f_i-128)_+, where f_i counts free positions in block i.
The first batch contributes fixed inverse coordinates, not free ones.

## Bound all head and tail selections simultaneously

Pull back the physical pairs through each region's route. This gives an
oriented uniform matching of outer blocks. Its heads are a uniform M/2-subset.
Let H_R(S0) count heads in a fixed s-vertex subset S0 across R regions, and
let E_R(S0) count matching edges internal to it. The
[previous graph bound](LAYERED_ROUTE_GRAPH.md) supplies E_R(S0)<=e_R(s)
simultaneously for all subsets and suffixes, except with probability 2^-56.

For one region, any specified j distinct vertices are all heads with
probability (M/2)_j/(M)_j <= 2^-j. Expanding factorial moments gives an
upper exponential-moment bound by Bin(s,1/2); independent regions give
Bin(R*s,1/2). With c_s the bit length of binom(M,s), Hoeffding's bound gives
the following integer envelope:

    h_R(s)=min(R*s, R*M/2,
               ceil(R*s/2)+ceil(sqrt(7*R*s*(c_s+77)/20))).

The same subset/size/suffix union as before has failure at most 2^-56.
The empty and full vertex sets are deterministic and need no tail union.
Both envelopes therefore hold except with probability at most 2^-55;
their events need not be independent.

For the blocks S0 with f_i>128, their total free degree is bounded by each of

    P+r,        P+e_R(s),        h_R(s)+r,        R*s.

The second bound counts at most one selected head per live edge, plus at
most one extra contribution for internal edges. The third uses all regional
heads, a superset of the live heads. Hence, on this route event,

    d_pi(F) <= d_R(r)
      := min(floor((P+r)/2), max(0,u-K),
             max(0,max_s[min(P+r,P+e_R(s),h_R(s)+r,R*s)-128*s])).

This holds for every choice of the r tails. The function d_R is nondecreasing.
There is no union over tail masks in the routing failure bound.

## Forward bound: maintaining the half-support forces cancellations

Consecutive batch pairings use distinct bits i and j. Any support with
exactly one member of every j-pair has exactly b/4 zeros among coordinates
whose i-bit is one: that half is itself a union of j-pairs.

Fix an outer message and all previous keys. At a one-sided incoming pair,
p is nonzero for every a. Expose a first. For each required zero in the
second coordinate of the next state, the equation

    0=y+k*p+X_second

specifies at most one c in S. It costs at most 1/(q-4). The constraints
use distinct pair keys; the argument also holds conditionally on the past.
There are L-2 internal transitions, so the total charge is
(L-2)*b/4=u/4-b/2. The final state is not charged.

A message producing zero output before the suffix has zero input there.
The shortened outer has dimension max(0,u-K), for every route. Union over
its messages and at most 2^P state-orientation masks gives

    log2 Pr[bad output in the family]
      < P+32*(max(0,u-K)-u/4+b/2)+1.                 (F)

This is a bound on all messages, not merely a fixed-message probability.

## Inverse bound: output zeros impose independent temporal equations

Let S_t be the nonzero support of each of the L-1 live incoming states,
and let r_t be its number of type-01 pairs. Each |S_t|=b/2, and sum r_t=r.
For consecutive, differently paired supports,

    |S_t intersect S_(t+1)|
      <= b/4+min(r_t,r_(t+1),b/2-r_t,b/2-r_(t+1)).

For example, an all-type-10 support intersects the next support in exactly
b/4 points; changing r_t orientations can increase that by at most r_t.
Use also the all-type-01 reference and exchange the two supports.
Summing the inequality gives a total intersection at most

    (L-2)*b/4 + min(r,P-r).

The first output is the first live state. At an internal transition, an
output zero inside S_t union S_(t+1) is possible only in their intersection,
and imposes Q_(t+1,j)=delta*Q_(t,j). Consequently an output of weight <=W
requires at least

    C(r)=max(0,P-W,3*u/4-b-min(r,P-r)-W)

such equations. The P-W term also uses that the potential output support
is at least the total P live state variables. Each equation joins successive variables along a coordinate's
temporal path. Any selected collection is a forest, so its rank equals
its number of edges. This is exact rank, not independent-state randomness.

There are P live state variables. Choose orientation masks in at most 2^P
ways and zero-edge sets in at most 2^P ways. The latter bound holds because
there are at most P available temporal edges. For every set with at least
C(r) edges, the number of assignments is at most q^(P-C(r)). Allow q^b
choices for the final output, whose outgoing state is unconstrained.
Combining this count with the inverse point law yields

    log2 Pr[bad output in the class r, good route]
      < 2P+32*(b-r-C(r)+d_R(r))+1.                 (I)

Nonzero restrictions on state variables only reduce this count. We do not
multiply independent probabilities for equations sharing a state variable.

## Complete integer coverage and checks

Suffixes R<=128 cannot contain a nonzero shortened outer word. For every
R=129,...,256 and every integer r=0,...,P, take the smaller of (F) and (I).
The source verifier requires a bound <=-80. It uses (F) for entire suffixes
when possible. Otherwise it covers integer r-intervals: d_R(r)<=d_R(hi),
-r<=-lo, and C(r) is minimized nearest P/2. It bisects until each entire
interval passes and independently checks gapless endpoint coverage.

The 40% run covers 100892800 classes: 39 suffix lengths pass by (F), and
the remaining ones use 919 accepted intervals. The largest accepted
interval bound in the current replay is -9631, at R=204, r in [313249,365456].
The original replay without the extra P-W rank term gave -7487.
This is a conservative
interval bound, not an optimized worst-class margin. Union over fewer than
2^27 classes at 2^-80, then add both routing failures, to obtain <2^-52.

The exact tests include 371 head factorial comparisons, 6144 oriented
graph/free-tail fixtures, 768 consecutive-wiring support checks, and 290
GF4 temporal-rank checks. All certificate arithmetic is integer. NumPy is
used only for int64 vector operations on quantities far below overflow.

Serial reproduction (no generated data):

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/layered_one_sided_gate.py
```

The same bounds do not close at 41% or 45%. The 41% run first misses its
-80 per-class target at R=210,r=428520 (upper bound -63); 45% first misses
at R=168,r=75732. These are failed bounds, not codeword counterexamples.

## Useful extension and the remaining gap

The earlier dense-suffix argument also admits a cheaper defect measure.
Let z,l,r,f count pair types 00,10,01,11 after the first batch, and let e
be exceptional 11 pairs. The inverse free-coordinate count is

    D=l+2r+2f-2e,       D+e=2P-e-(2z+l).

Complete each free half-edge to a full edge, then use monotonicity of the
old graph dimension bound at P-e edges. The old dense-suffix certificate
therefore loses only 2z+l weight units, not two for every nonfull pair.
In particular, type 01 has no defect penalty. The old threshold W854892
extends to all histories with W+2z+l<=854892, with the same restricted
52-bit failure bound. No additional support-mask union is necessary:
that argument already unions all output vectors.

For a general mixed history with V=l+r+2f live state symbols, the temporal
rank method exposes a useful cancellation:

    V-(D+e)=e-r.

If its output-zero temporal constraints have rank C, its inverse amplitude
exponent is b+e-r-C+d_pi(F), before support/branch counting. This identity
does not close mixed histories by itself: both C and the free-set geometry
depend on where the empty, one-sided, and full pairs occur. The next step
is a joint bound retaining these quantities, rather than assuming dense
states or a uniform refresh. No construction change or new timing is claimed.
