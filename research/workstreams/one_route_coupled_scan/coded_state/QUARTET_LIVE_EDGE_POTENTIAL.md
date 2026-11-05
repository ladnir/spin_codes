# A limit of the five-state temporal edge bound

The edge relaxation below cannot certify the remaining profile
(h,R,Ls,f,F)=(1536,753664,730726,4,851968), for any continuous
mu>=0, 0<=r<=4, or beta>1. Its exact conditional mean already exceeds
the required threshold throughout that parameter domain. This is a
limitation of this counting relaxation, not a bad-codeword construction.
Further marker grids for this profile are parked. The full-code goal
remains open.

The source is
[quartet_live_edge_potential.py](quartet_live_edge_potential.py).
It follows the actual-state and free-terminal conventions of the
[live-selector quotient](QUARTET_LIVE_SELECTOR_QUOTIENT.md).
The construction and setup distribution are unchanged.

## A setup statistic that bounds every selector history

Fix an eligible row set H of size h and every quartet partition. Let
j_q be the number of H-labelled slots in a quartet. Join quartets in
successive update layers by their common physical coordinate; each
quartet has four incoming and four outgoing coordinate edges, except
at the endpoints. There are n(T-1) edges, where n=512 and T=4096.

For mu>=0 and 0<=r<=4, define the edge potential

    g(k,j)=max_(a,b in {0,1})
             [ab+(mu/8)(a(k-r)+b(j-r))].                 (1)

The bits a,b indicate selected live quartets. Thus the maximum in (1)
bounds every selector assignment without a union over such assignments.
Put Z_(mu,r)(H)=sum_edges g(k,j). For an actual profile, let
G=F-f and C_out=R-128h. The necessary threshold is

    S_(mu,r)=Ls-512-f
       +mu*((1-r/4)G-C_out)
       -(mu/2)*128*(4-r).                               (2)

The first subtraction accounts for the free final state. The last term
accounts for the last selector layer, which has only one incident edge
layer. The first selector layer is empty. Any matching actual history
implies Z_(mu,r)(H)>=S_(mu,r).

Under independent slot labels with probability p, quartet counts have
p_j=binom(4,j)p^j(1-p)^(4-j). For beta>1, set

    L_(k,j)=p_j beta^(4g(k,j)).

If v is positive and lambda>=max_k (Lv)_k/v_k, Holder's inequality at
each new quartet gives

    E beta^Z <= lambda^(128(T-1))
                    * (sum_j p_j v_j / min_j v_j)^128.   (3)

The bound holds for every fixed partition sequence. Exact regional
quotas can be restored by dividing by
pi_h(p)^256, where pi_h(p)=binom(8192,h)p^h(1-p)^(8192-h).
A single binom(8192,h) factor unions the row sets H. Distinct
(h,mu,r) cut classes require a finite union if combined; choosing
p,beta,v for one fixed statistic requires no additional union.

## Test the exact conditional mean before optimizing markers

The true regional routing distribution is available without an iid
approximation. Write M=8192 and (a)_b=a(a-1)...(a-b+1). A quartet count
has marginal probability

    P_j=binom(4,j)(h)_j(M-h)_(4-j)/(M)_4.

For two quartets in the same region, their slot sets are disjoint, so
their joint eligible-count probability is

    P_same(j,k)=binom(4,j)binom(4,k)
                    (h)_(j+k)(M-h)_(8-j-k)/(M)_8.        (4)

Across a regional boundary it is P_j P_k. There are 3840*512 edges
within regions and 255*512 edges across their boundaries. Equations
(1) and (4) therefore give E[Z] exactly, conditional on the regional
quotas and on every fixed partition sequence.

If E[Z]>=S, Jensen's inequality gives

    E[beta^Z] / beta^S >= beta^(E[Z]-S) >= 1.             (5)

Consequently even the exact conditional Chernoff expression cannot
certify that cut. Replacing its moment by the upper bound (3), or
paying an H-set union, cannot improve this obstruction.

## Exhaust the continuous parameter domain

Put nu=mu*r and order the edge labels as j<=k. Then

    g(k,j)=max(0,(k*mu-nu)/8,
                1+(j+k)*mu/8-nu/4),

on the cone mu>=0, 0<=nu<=4mu. The required score becomes

    S=Ls-512-f +mu*(G-C_out-256)+nu*(64-G/4).             (6)

Thus E[Z]-S is a convex piecewise-affine function. Its break lines lie
in the three families

    nu=k*mu,    nu=(j+k)*mu/2+4,    nu=j*mu+8.

Intersecting these lines and the cone boundaries gives 21 feasible
vertices, including the origin. On each arrangement cell the function
is affine. A finite vertex check is complete once its recession slopes
are nonnegative. Along nu=r*mu, the recession slope is

    128(T-1) sum_j P_j (j-r)_+
      -(G-C_out-256)-r*(64-G/4).                         (7)

This is affine between integer r, so only r=0,1,2,3,4 need be checked.
A negative slope is reported as an unbounded direction; it must not be
discarded in favor of the finite vertex minimum.

The parent's exact rational evaluation found:

| Actual profile (h,R,Ls,f,F) | Minimum vertex gap E[Z]-S | Recession |
| --- | ---: | --- |
| (1536,753664,730726,4,851968) | 36014.3238686... at mu=4,r=2 | All nonnegative; minimum 15628.6... |
| (2048,819200,819200,4,950272) | -946.31256... at mu=4,r=2 | All nonnegative |
| (2560,786432,779879,4,950272) | No finite global minimum | A negative slope |

For the first row, the exact comparisons are
36014<minimum gap<36015 and 15628<minimum recession slope<15629.
Together with the exhaustive vertices and rays, they prove the
obstruction for all continuous parameters in the stated domain. The
second row leaves only a small mean gap for this method. Neither
conclusion extends to a different setup statistic.

Before this mean audit, a bounded FLOAT control for the third profile
gave about -939.98 bits at r=2, mu=8, beta=17/16 and p approximately
0.26876858. It remains a diagnostic requiring exact verification.
The negative recession in the table is consistent with possible
selector-capacity improvement; it is not itself a probability bound.

The parent-run checks passed 75 local maxima, 210 rational Holder
comparisons, 114 width-4/8 edge checks, and 1500 boundary checks.
They include physical slot enumeration and negative controls for
omitting either endpoint correction. The source additionally exposes
the exact 21-vertex/five-ray audit and its three profile controls.
No further marker search for the obstructed h=1536 profile is warranted
within this relaxation. The next step requires a statistic retaining
more than this local maximum over adjacent selector bits.
