# Joint geometry of sparse component classes

Retaining component lengths gives a valid finite improvement, but it does
not yet bound the permitted-triple term. The fixed, zero-offset L4 placement
has a particularly difficult profile: rows (0,1,2) produce cuts (0,1,2) in
every region. At the tested mark 1149/1150, its common-control upper is
2^18.09317428, even without position bins. This is a failed upper bound,
not a code counterexample or a lower bound on failure probability.

Two reusable bounds are established below. The first reduces a relaxed
joint-length problem to finitely many simplex corners. The second bounds
the post kernel uniformly over an occupancy interval without replacing all
return paths by an occupancy-independent error. The latter is substantially
more useful for fixed placement and for averaging setup offsets.

The implementation is `theory/hd_component_geometry.py`; seven tests in
`theory/hd_component_geometry_tests.py` pass. Existing files remain unchanged.
No native build, remote command, performance benchmark, or commit was made.

## The contraction to be bounded

Use the positive expansion in [HD_COMPONENT_RANK.md](HD_COMPONENT_RANK.md).
Here a is the number of **permitted** source rows, not the number of active
rows in a codeword. Fix their set J. In each region, its sorted cuts divide
the prefix into lengths

    l_0 >= 0,    l_1,...,l_a >= 1,    sum_i l_i = R.

The first component is carried in; the last is carried out. The other a-1
components are local. A component's nonzero flag has Q=q-1 possible values.
For old/new carry flags s,t and local flags z_1,...,z_(a-1), the input
occupancy is

    h = s*l_0 + sum_i z_i*l_i + t*l_a.

Let F_h(v) be the normalized regional kernel of the common zero/live
post-state control. The component transfer is the four-by-four matrix

    K(l)[(s,i),(t,j)]
      = sum_z Q^(|z|+t) F_(s*l_0+z·l+t*l_a)(v)[i,j].

Its normalization is C(R,h), because routing is uniform within each region.
The component flag and post state both carry between regions. For the actual
fixed cut profile, its nonzero moment is e_0 K(l^(0))...K(l^(63)) 1-1.

At K20, R=32768, q=65536, and the forbidden weight is e=838860. The
contribution for this particular permitted triple, after dropping a factor
at most one, is Q^-97 v^-e times its nonzero moment. Summing such quantities
over all permitted triples is necessary; an isolated profile is not a class
bound.

## Convexity of the no-return birth term

Write p=1/Q and beta=v+(1-v)p. If h>=1 live input positions are uniformly
routed, the no-return zero-to-live term is

    B_h = v/C(R,h) sum_(t=h-1)^(R-1) C(t,h-1) beta^t.

Summation by parts gives the positive mixture

    B_h = v beta^(R-1)
          + v(1-beta) sum_(t=1)^(R-1)
              beta^(t-1) C(t,h)/C(R,h).

For fixed t, put f_h=C(t,h)/C(R,h), taking C(t,h)=0 when h>t.
The first and second forward differences satisfy

    Delta f_h = -(R-t)/(R-h) f_h,
    Delta^2 f_h = (R-t)(R-t-1)/((R-h)(R-h-1)) f_h.

Thus B_h is decreasing and discretely convex on 1,...,R. Its piecewise
linear interpolation is convex. Exact tests verify both formulas and the
resulting convexity at small R.

This property concerns the birth term, not every entry of F_h. In particular,
F_0[0,1]=0 must not be joined to positive-occupancy birth terms by an asserted
convexity across zero.

## A finite corner reduction, and its limitation

For variable positive occupancies, retain B_h+eta in the birth entry, where
eta=(1+p^2)^R-1. Use the existing minimum-span return bound for the zero-to-zero
entry and the uniform live-to-live bound from
[HD_COMPONENT_PAIRS.md](HD_COMPONENT_PAIRS.md). Retain F_0 and F_R whenever
the flag pattern forces occupancy zero or R. The reverse birth entry is
p^2 times the birth entry. These choices define an upper transfer K_bar(l).

Separate two faces: l_0=0, and l_0>=1. On either face, each flag pattern has
a fixed set of positive-length components. Its minimum live occupancy is
therefore fixed. Any nonconstant occupancy lies strictly between zero and R,
so the convexity just proved applies. Each entry of K_bar is convex jointly
in the lengths on that face.

Every feasible length vector is a convex combination of the face's simplex
vertices. Consequently the **entire matrix** K_bar(l) is bounded entrywise
by the same convex combination of its vertex matrices. The vertices put
all excess length into one component. There are a vertices on l_0=0 and
a+1 on l_0>=1, hence at most 2a+1 matrices.

For three cuts this gives seven matrices. A positive backward recursion
can bound every regional product: start from the all-one column vector,
and at each step take the maximum of K_vertex times that vector separately
in each row. Monotonicity proves the resulting upper for any fixed cut
profile. Exact tests check the convex-combination inequality for every
small cut set, and check the final two-region bound with global carries.

This relaxation is still too weak at K20. Its directed all-triples upper
logarithms at marks 799/800, 999/1000, and 1199/1200 are approximately
822.152, 568.116, and 410.350. The rowwise maximum may choose a different
geometry after each state branch. Thus this valid reduction does not yet
exploit the cross-region restrictions of the fixed affine placement.

## A post-kernel upper for a whole occupancy interval

The closed run sums in HD_COMPONENT_RANK express each entry using a
nonnegative coefficient depending on h and a selected-zero-gap moment.
For a fixed number L of gaps, define

    M_L(h) = E[beta^(-G_L)],

where G_L is the sum of L gaps of a uniform h-subset of the R positions.
The gaps are exchangeable, so take the first L. Their sum is the position
of the L-th selected element minus L-1. Couple uniform subsets by taking
successive prefixes of one random permutation. Adding a selected element
can only move that order statistic earlier. Therefore M_L(h) decreases
with h whenever the term is defined. The endpoint L=h+1 is also monotone:
there G_L=R-h is the maximum possible zero-gap sum.

Fix an occupancy interval l<=h<=u. In a retained run term, take its
binomial coefficient at u and its gap moment at the smallest feasible h:

| Entry | Run index | Coefficient upper | Moment occupancy |
|---|---:|---|---:|
| F01 | t>=0 | C(u-t-1,t) | max(l,1,2t+1) |
| F11 | t>=0 | C(u-t,t) | max(l,1,2t) |
| F00 | t>=1 | C(u-t-1,t-1) | max(l,1,2t) |

Only feasible terms are included. Their remaining factors are unchanged
from the closed run sums. F10=p^2 F01 for this common control. The existing
omitted-run tail is valid with lambda=u*p^2. If the interval includes zero,
also dominate the two diagonal entries of F_0.

This proves an entrywise upper for every h in the interval. It does not
assume that the full kernel is monotone in h. It also avoids charging a
uniform eta to the birth entry when the entire interval has large occupancy.
Every small occupancy interval is checked against exact polynomial kernels.

A useful special identity is

    one-return contribution to F_h[1,1]
      = (h-1) (v/beta) p^2 B_h.

Thus an interval upper for this term is (u-1)(v/beta)p^2 B_l for l>=1.
It can be capped by the earlier sum-over-all-occupancies upper.

## Retaining the actual L4 orbit

For a triple inside source group g, write its lanes as
0<=r_1<r_2<r_3<4. In region c its cuts are

    4*(a_c^-1*g mod8192) + (r_1,r_2,r_3).

Both interior lengths remain at most two, and their sum is at most three.
Only the complementary exterior lengths change with the affine orbit.
The helper bins the physical group position finely near either scan end
and geometrically in the interior. It uses the interval kernel above for
each live-component flag pattern.

When the old and new carry flags agree, the exterior contribution is
constant. The helper retains that equality instead of adding independent
extrema of the two complementary exterior intervals. It then multiplies
the bin matrices along the actual 64-region orbit of every source group.
Taking a positive entrywise maximum over the four lane triples gives an
upper for each of them; multiplying by four covers their sum.

Group zero is handled separately. Every multiplier fixes it, and the
triples with first cut zero and first cut one lie on different faces.
Their four full-profile moments are evaluated separately, avoiding an
extra state-dependent choice between those faces.

At v=1149/1150, 128-bit directed Arb arithmetic and 273 position bins give

    log2(upper for all 32768 within-group triples) = 18.0935120427... .

The four group-zero profile upper logarithms are

| Cuts in every region | Upper log2 contribution |
|---|---:|
| (0,1,2) | 18.0931742762 |
| (0,1,3) | -44.1335230066 |
| (0,2,3) | -44.1335230066 |
| (1,2,3) | -4.86484262949 |

The largest nonzero-group profile upper is approximately 2^-1.33078.
These results concern only triples contained in one L4 group. Triples
meeting two or three groups are not included. More importantly, even this
subclass upper fails, and the difficult group-zero term is evaluated without
position-bin relaxation. Refining those bins cannot remove that term.

## What to try next

Geometry alone has not closed the fixed-placement sparse envelope. A
separate total-degree screen reported by the parent analysis also failed
to improve the (0,1,2) profile at this mark: increasing its degree mark
from one only increased the upper. Full per-row degree filtering remains
distinct and unresolved. The more direct next comparison is independent
cyclic group offsets in setup, which change the repeated group-zero
profile. This recommendation does not treat a failed upper as a
counterexample.

Independent setup offsets have a clean interface to this contraction. For
each fixed J, average each regional transfer over its independent offset
before multiplying the 64 transfers. Linearity and independence give
E[K_0...K_63]=(E K_0)...(E K_63), with the carry states still present.
Offsets must be multiples of four to preserve grouped placement. A cyclic
shift can reorder the cuts and change the carry split, so its average must
include that wrap; merely shifting an exterior length is incorrect.

Reproduce the tests and retained finite diagnostics from the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/hd_component_geometry_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_component_geometry.py --denominators 800 1000 1200
python -B research/workstreams/char2_field_spin/theory/hd_component_geometry.py --grouped-triples --denominators 1150 --bin-subdivision 16
```

The interval-kernel and corner inequalities are proved upper bounds. None
of the reported failed gates is a lower bound on setup failure, and none
is a complete distance certificate.
