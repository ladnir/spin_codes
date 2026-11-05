# Global rank: all full/empty bursts

2026-10-31. The fixed field-quartet construction is unchanged. This note
closes a complete support class, across every outer-row count: all
nonzero histories in which each 512-coordinate state is entirely full
or entirely empty. Their total probability of producing an output of
weight at most W=943718 is less than 2^-148. The full 45%/40-bit code
certificate remains open for other histories.

This is stronger in scope than the old burst-cut checkpoint, which
treated 2100 inactive outer rows. No such row-count restriction appears
here. The argument uses actual regional routing, accepted field maps,
and one free terminal state. The theorem is over J=GF(2^32); do not
transfer this restricted support class alone to the128-bit extension.

## Rank and selected output zeros

Use N=2^21, K=2^20, M=8192 outer rows, and T=4096 updates. Fix a state
support S and quartet partitions before exposing the route. Let A count
nonzero state coordinates over time, D state deaths, and I surviving
nonzero states. Let V be the covered set from the global-rank note, and
E=|V|-A. Write

    R(V)=sum_i min(128,v_i),
    dshort=|V|-R(V)=sum_i(v_i-128)_+,

where v_i is the number of covered coordinates in outer row i.

Every birth or death forces a nonzero output. Thus a low-weight output
needs at least z=max(0,A+D-W) zero outputs at surviving coordinates.
Choose z of those coordinates. Each chosen zero relates the new state
value to its preceding value within the same excursion. Its new-state
pivot lies in L. These pivots are distinct and independent of the
residual RS equations, whose variables are births outside L.

Normalize the first nonzero state value to1. The joint rank count gives

    Pr_field[bad output with state support S | route, partitions]
      <= min(1,(gamma/q) binom(I,z) q^(dshort-E-z)),           (1)

where q=2^32 and gamma<2 is the global setup acceptance-density factor.
If z>I, the event is impossible. The union in (1) chooses z zeros;
it leaves all other outputs unrestricted. No geometric tail of further
zero choices is needed. The statement covers all messages with this
state support, including ones selected after setup.

Choose t=1/10 and theta=t/ln(q). Taking the theta power of the conditional
bound before averaging the route yields

    Pr[bad output with S]
      <= binom(I,z)^theta exp[-t(E+z)] E_route exp(t*dshort). (2)

We dropped the factor (gamma/q)^theta<1. Quartets remain arbitrary in
this bound; no union over their realizations is required.

## Exact burst counts and regional variance

For a full/empty history, let r count updates that change a full state
to an empty state. Then

    D=E=512r,  |V|=A+E,
    number of histories with r deaths = H_r=binom(4097,2r+1).

The count includes either terminal state. It follows by selecting the
endpoints of the one-runs in the binary state sequence. There are
0<=r<=2048 possibilities for the death count.

Let p_j be the fraction covered in region j. In any one outer row,
the degree Z is a sum of independent Bernoulli(p_j) variables: the
regional permutations are independent. Its mean and variance satisfy

    mu=|V|/M,
    Var(Z)<=v_r:=min(64,(r+1)/2).                            (3)

To see (3), the binary covered stream is the union of the current and
preceding full-state indicators. It has at most r closed runs and one
terminal run. At most2r+1 regions are partially covered; every other
regional probability is0 or1. Each partial region contributes at most
1/4 to the variance. The looser (r+1)/2 is convenient below.

A compatible nonzero outer word needs at least129 occupied regions.
If c regions are completely covered, c<=floor(mu), so

    129<=number of occupied regions<=floor(mu)+2r+1.          (4)

In particular, when mu<124, r<=2 is impossible. For3<=r<=6,
compatibility requires mu>=128-2r.

## A row moment without independent-row assumptions

Row degrees are dependent because each region has an exact quota.
Nevertheless, for nonnegative row functions that are all increasing,
or all decreasing, their product expectation is at most the product
of their marginal expectations.

An elementary coupling proves the needed inequality. Condition on one
row's incidence vector. Increasing its incidence in a region removes
one slot from the remaining quota. Generate the other rows from one
random ordering and nested prefixes. Their conditional product
expectation then changes in the opposite direction to the chosen row
factor. The chosen row's regional incidences are independent Bernoulli
variables, so opposite-monotonicity correlation bounds the product by
the product of expectations. Remove that factor and repeat. Unused rows
remain in the sampling population throughout this argument.

Consequently,

    E_route exp[t*sum_i(Z_i-128)_+]
      <= [E exp(t(Z-128)_+)]^M,                              (5)

and the same statement holds for (128-Z_i)_+. This preserves the actual
sampling without replacement. It is not a replacement by independent rows.

For X=Z-mu and c(t)=e^t-1-t, independence within one row gives
E exp(tX)<=exp(v_r*c(t)). The same bound applies to -X. For h>0,

    E(X-h)_+ <= v_r/(4h),
    log E exp[t(X-h)_+]
      <= t*v_r/(4h)+e^(-th)(exp(v_r*c(t))-1).                (6)

For completeness, set R_t(x)=e^(tx)-1-tx. It is nonnegative, and
R_t((x-h)_+)<=e^(-th)R_t(x). Also
E(X-h)_+=(E|X-h|-h)/2 <=(sqrt(v_r+h^2)-h)/2. Taking expectations and
using log(1+x)<=x proves (6).

At h=0 we instead use the simpler bound

    log E exp(t X_+) <= v_r*c(t)+t*sqrt(v_r)/2.              (7)

The same formulas hold for -X. All following comparisons use

    c(1/10)<13/2500,
    exp(64*13/2500)<7/5,
    exp(-(1/10)*(64/5))<1/3,
    ln2<7/10.                                               (8)

The helper checks rational Taylor bounds for all four inequalities.

## Coverage below124M

Put mu0=576/5; W/M<mu0. Since I<=N,
binom(I,z)<=N^z, so its entropy in (2) costs at most21z/32 field
dimensions. Equations (2), (5), and (6) bound the natural logarithm of
the full r-class contribution by

    ln H_r -t[512r+(11/32)z]
      +M*v_r[t/(4(128-mu))
              +exp(-t(128-mu))*(13/2500)*(7/5)].             (9)

For mu<=mu0, drop z and use the first endpoint from (8):

    L_low(r)=(7/10)ceilLog2 H_r -t*512r
          +M*v_r[1/512+91/37500].                           (10)

For mu0<=mu<124, use z>=M(mu-mu0). The right side of (9), with this
substitution, is convex in mu. Its maximum lies at an endpoint. At124,
the valid coarse bound exp(-4t)<=1 gives

    L_end(r)=(7/10)ceilLog2 H_r
          -t[512r+(11/32)M(124-mu0)]
          +M*v_r[t/16+91/12500].                            (11)

For r>=7, (10) and (11) cover every mu<124. For3<=r<=6, condition (4)
restricts mu to [128-2r,124). Check that lower endpoint in (9), again
using exp(-t(128-mu))<=1, and check (11). There is no remaining low-r
case below124M.

Exact rational arithmetic verifies that every resulting r-class bound
is less than -112. The largest is at r=7 and equals

    -1095778/9375 = -116.8829866... .

This is a bound on the full r-class, not an estimate at one coverage
profile: H_r already counts every state sequence with that death count.

## Coverage at least124M

Here use binom(I,z)<=2^I. If124<=mu<=128, inequality (7) bounds the
positive part (Z-128)_+ by the centered positive part. The remaining
field deficit in (2) is at least

    E+|V|-W-I/32 >= E+(31/32)|V|-W >= E+40346.

If mu>=128, use the exact fixed-total identity

    dshort=|V|-K+sum_i(128-Z_i)_+.

After taking out that deterministic term, I<=N leaves a field deficit
of at least E+K-W-N/32=E+39322. Thus all mu>=124 are covered by

    L_high(r)=(7/10)ceilLog2 H_r -t(512r+39322)
             +M[(13/2500)v_r+(t/2)s_r],                    (12)

where s_r is sqrt(v_r) rounded upward to a multiple of1/16. This uses
only integer square roots. Every r=0,...,2048 passes L_high(r)<-112.
The largest bound is -1767647/625, at r=43.

There are2046 low-coverage r-classes and2049 high-coverage r-classes,
4095 in total. Each contributes less than exp(-112)<2^-160 by (8).
Their union is therefore less than2^-148. This proves the claimed
all-row full/empty-burst theorem.

## Regional operator diagnostic and remaining work

The same helper also implements an exact regional proof interface:
the513-state quartet kernel, composed over16 updates, retains boundary
state weight. A scalar renewal identity counts paths visiting at least
129 positive-footprint regions, removing the empty path exactly.
Twenty small renewal comparisons and106 kernel comparisons passed.

At full width, four median-coverage diagnostics remain positive, even
after choosing their output markers separately. They are not
certificates or obstructions to the construction. An early ordinary
coefficient calculation underflowed to negative infinity; it was rejected
and replaced by normalized vectors and logarithmic scalar convolutions.
The corrected diagnostics retain all destination weights and impose the
129-region condition, but their entrywise regional envelopes still lose
too much information to close the remaining histories.

The analytic theorem above does not use those floating diagnostics.
An independent read-only audit checked the selected-zero count, MGF,
branch coverage, free-terminal count, and failure allocation. Its4095
checks are exact rational comparisons in
`quartet_regional_rank_probe.py --burst-cert`.

Next extend the all-row burst argument to partially occupied states,
using the cost of suppressing quartet expansion and the edit distance
from an empty/full history. Keep the actual regional variance and
forced-input constraints. The full goal stays active; no construction,
paper, library, performance, or setup-distribution changes were made.
