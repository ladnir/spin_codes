# The complete permitted-pair term

The a=2 term of the [positive component expansion](HD_COMPONENT_RANK.md)
is at most **2^-41.11393101** at K=2^20 and forbidden weight838860.
This includes every permitted pair of source rows and every fixed regionwise
row order. It covers both post-inner distributions through their common
positive upper operator. It is not a whole-code distance certificate:
all terms with a>=3 remain unbounded.

Here a counts permitted rows in an auxiliary expansion. It does not count
the exact active outer rows of an actual codeword. Bounds for different
values of a must be added, not interpreted as independent failure events.

## Why the singleton bound needs refinement

Fix a region of R positions and condition on exactly h nonzero inputs.
The regional permutation gives a uniform h-subset. Let F_h(v) be the
normalized kernel of the shared zero/live upper operator

    T0 = [[1,0],[0,b]],   T1 = [[0,v],[v*p^2,b]],
    p=1/(q-1),           b=v+(1-v)*p,        0<v<1.

These are positive upper operators, not stochastic matrices. The existing
closed run sums compute F_h. For two permitted rows, each region has one
interior component and a component carried to the next region. Summing over
their field values can multiply a kernel entry by (q-1)^2. The loose uniform
return bound used for singletons therefore loses too much.

We tighten the live-to-live and zero-to-zero entries. We also retain a
minimum occupancy of two when two disjoint positive-length components are
both live. These changes do not modify the construction.

## Uniform positive-occupancy kernels

Write lambda=R*p^2 and eta=(1+p^2)^R-1. The application has lambda/3<1.
Let G_h be the initial zero gap of a uniform h-subset. The no-return birth
term is

    B_h = v*b^(R-1) E[b^(-G_h)].

The first selected position moves earlier as h increases. Thus, for every
h>=m>=1, B_h<=B_m. Return terms give the uniform bound

    F_h[0,1] <= B_m+eta,     F_h[1,0] = p^2 F_h[0,1].

The equality belongs to the two-state upper kernel. It is not an equality
between aggregated transitions of the actual four-state support process.
The helper evaluates B_m by the positive beta-binomial factorial series.

For F_h[1,1], the no-return term is b^R. With exactly one return, there are
h-1 possible intervening zero gaps. Summing that positive contribution over
all h bounds each individual occupancy. Vandermonde's identity gives

    sum_h h(h-1) C(R-h,g) = 2 C(R+1,g+3).

After normalizing the uniform-subset gap law, the resulting bound is

    S = sum_(g=0)^(R-2)
          b^(R-2-g) * 2(R+1)(R-g-1)/((g+1)(g+2)(g+3));

    F_h[1,1] <= b^R + v^2*p^2*S + lambda^2/(2(1-lambda/3)).

The last term bounds paths with at least two returns. A path with t returns
has total coefficient at most lambda^t/t!, since all output-mark factors
are at most one. The ratios from t=2 onward are at most lambda/3.

For F_h[0,0], one return requires h>=2 and spans the first and last live
inputs. Its output-mark factor is v^2*b^(range-1). Couple a uniform h-subset
with a uniform two-subset contained in it. Its range can only increase,
so the one-return moment is maximized at h=2. Direct summation gives

    E[b^(range(2)-1)]
      = 2((R-1)-R*b+b^R)/(R(R-1)(1-b)^2).

Consequently

    F_h[0,0] <= v^2*p^2*E[b^(range(2)-1)]
                + p^2*((1+p^2)^(R-1)-1).

The second term sums the remaining return paths. This bound also covers
h=1, whose exact zero-to-zero entry is zero.

Together these inequalities define an entrywise positive bound P_m for
every h in {m,...,R}. Let A bound every h in {0,...,R}: use P_1 except
replace its zero-to-zero entry by one. Keep F_0 exact and F_R bounded by
the directed closed run sum.

## Carrying the two permitted cuts

Fix two cuts 0<=j1<j2<R in a region. They form pieces of lengths
j1, j2-j1, and R-j2. The first component is carried from the previous
region; the middle component is local; the last component is carried on.
The middle and last lengths are positive.

Set Q=q-1. With block order given by the old/new carried component's
zero/live flag, the following four-by-four matrix dominates every pair
of cut positions. Each displayed block is two-by-two on post-state flags:

    M = [[F_0+Q*P_1,       Q*P_1+Q^2*P_2],
         [A+Q*P_1,         Q*P_1+Q^2*F_R]].

Each new live component has Q possible nonzero values. The upper-right
Q^2 branch has two live pieces of positive length, so its occupancy is
at least two. In the lower-right Q^2 branch every component is live,
giving occupancy R exactly. Bounding either branch as an arbitrary
positive occupancy would discard useful information.

The initial carried component and post state are both zero. There are
64 regions, with global carries. Therefore

    W <= e_0 M^64 1 - 1

bounds the summed output moment of nonzero vectors in each permitted-pair
component space. The subtracted one is the unique all-zero assignment;
its exact contribution remains one in the upper operator.

Each pair's factor in the positive expansion is
Q^-64*(1-Q^-32)^(R-2)/Q. Dropping the factor at most one and using a
class-specific mark yields

    pair-class contribution
      <= C(R,2) * Q^-65 * v^(-838860) * (e_0 M^64 1 - 1).

This is a sum over all permitted pairs, not a claim that their events are
disjoint. Independence of those events is unnecessary.

## Directed result and validation

At R=32768, q=65536, and the rational mark v=1579/1580, directed Arb
arithmetic gives

    log2(per-pair upper)  = -70.11388698359675...;
    log2(all-pairs upper) = -41.11393101182719... .

The common post operator covers the conditioned four-product inner and
the independent-coefficient [three-product inner](HD_THREE_PRODUCT.md).
No additional setup-conditioning factor is needed.

Five tests pass. They verify the Vandermonde sum and two-point range
formula exactly, compare every small occupancy with the uniform bounds,
check all small cut pairs, and replay the directed K20 gate. Independent
review added 636 exact/local comparisons and a 160-bit precision replay.

```text
python -B research/workstreams/char2_field_spin/theory/hd_component_pair_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_component_pair.py
```

The singleton and pair terms together use less than 47% of a 2^-40 total
failure budget. The remaining terms cannot be omitted. Next extend the
component transfer to a>=3, retaining how many positive-length components
are live. Larger permitted sets will also require their length structure;
a separate worst-case entry for every transition may lose too much.

This round changes only proof helpers and notes. Native timings and all
frozen certificates are unchanged.
