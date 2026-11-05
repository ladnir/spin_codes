# Retain regional quotas, not only the global count

2026-10-30. The unchanged construction's full 45% distance / 40-bit margin
goal remains open. This note rules out the proposed global-subset moment,
then gives an exact transfer for the actual regional distribution. A block
relaxation reduces the remaining computational question to one region.
No new whole-code failure allowance is claimed.

Continuation: [the boundary-state audit](QUARTET_BOUNDARY_POTENTIAL.md)
proves that the one-region free-terminal product below cannot close the
middle class, even with exact coefficients. Two-region blocks also fail
at the selected weights. The exact regional transfer remains valid; its
new product-potential envelope retains the state across those boundaries.

## A forbidden suffix defeats even exact global coefficient extraction

Fix the middle class from [the preceding audit](QUARTET_WINDOW_QUOTA.md):
N=2097152, R=537600 inactive positions, c=111572 paid positions, and
u=R-c=426028 unpaid positions. Its field constraints give B/lambda<=825555
for every nonnegative multiplier in the current family. Recall that H is
the maximum of lambda*A-tau*D over allowed state supports. The terminal
state is free.

Suppose U were uniform among all u-subsets of the N input positions. Take
L=825856=1613*512. The event that U avoids the last L positions has probability

    binom(N-L,u) / binom(N,u).

On this event, choose the state support to consist of all coordinates in
those last1613 updates. Its own inputs and every subsequent quartet are
active. There are no deaths, including at the terminal boundary, so
H>=lambda*L>B. This is a geometric support, not a codeword claim.

Consequently, for EVERY b>=1 and every current field-multiplier choice,

    binom(R,c) b^-B E_uniform_u[b^H]
      >= binom(R,c) binom(N-L,u) / binom(N,u).

The directed log2 lower bound is

    42493049717272959/1099511627776 >38647.20358.

This is a lower bound on the exact fixed-count moment expression, not only
on a Cauchy or conditioning relaxation. Multiplying by the density factor
D(u)>=1 and the row-set union cannot repair it. Thus the proposed exact
global-count transfer cannot close this class, regardless of computational
resources. The preceding audit's open possibility of exact global
coefficient extraction is now ruled out for this particular comparison.

The event is impossible under the actual distribution. The suffix contains
100 complete regions, each with2100 original inactive positions. Paying
all of them would require at least210000 payments, but c=111572.

More precisely, let a region have m positions and r>0 original inactive
positions. For c<r*Breg payments over Breg regions, the longest possible
all-active suffix has length

    floor(c/r)*m + m-r + (c mod r).

Every complete terminal region costs r payments. In the preceding region,
the remaining r-(c mod r) unpaid positions must fit before the suffix.
Placing them there attains the bound. Our parameters give440540 positions,
well below825856. Retaining a global count alone loses this restriction.

## An exact transfer under the actual regional distribution

Fix the original inactive-row set. Its routed inactive positions I_j are
independent uniform r-subsets of the m positions in each region. The random
quartet partitions are independent of these sets. Define

    Z_c(b)=E_[I,partitions] sum_(P subset I, |P|=c) b^H(P).

For one region, fix the allowed-input mask after paying P. If it contains
a active positions, then its unpaid count is m-a. Therefore the number of
payments in this region must be c_j=a-(m-r). There are exactly
binom(a,c_j) original inactive sets consistent with that allowed mask:
choose the paid positions among its a active positions and add every
unpaid position. Its weight in Z is consequently

    binom(a,c_j) / binom(m,r),  0<=c_j<=r.               (1)

Use the five-type kernel from the preceding audit. Each transition records
the raw allowed-mask multiplicity and increments a regional active counter
by l. At a region boundary, discard a<m-r, multiply by (1), add c_j to a
global paid counter, and reset ONLY the regional counter. Carry the five
state-type counts unchanged into the next region. Apply the free-terminal
factor only after the final update. Extracting paid count c gives Z_c(b)
exactly. This does not replace the regional input distribution by fresh
independent bits.

Without moment weights the result is binom(Breg*r,c), for every c. This
normalization follows also by counting the subsets of I for each fixed I.
Do not multiply the result by that binomial coefficient a second time.

The middle-class union bound would be

    binom(8192,2100) b^-B(c) Z_c(b),

followed by a sum over the remaining covered counts. That expression is
not evaluated or certified here. The full-width five-type state space is
still too large for brute force.

## A regional block relaxation with an explicit boundary cost

There is a simpler upper-bound interface if carrying state between regions
is too expensive. Within each region, treat its final allowed-state bit as
its own active-input bit, and give the region a free terminal state. Start
the run calculation afresh in the next region. Let H_j be the resulting
regional maximum. For every fixed input mask and every partition sequence,

    H <= sum_j H_j <= H+(lambda+tau)*512*(Breg-1).       (2)

For the first inequality, restrict a globally feasible support to each
region. Removing cross-boundary death charges and the last-bit constraint
cannot lower its score.

For the second inequality, combine regional maximizing supports. At each
internal boundary, delete a selected last coordinate if it violates the
true next-quartet constraint. Each deletion loses at most lambda of activity
and introduces at most one death of cost tau. If the last coordinate is
retained, restoring its cross-boundary death costs at most tau instead.
These alternatives cost at most lambda+tau per coordinate per boundary.
The resulting support is globally feasible. This argument does not assume
independence between adjacent windows.

Define the one-region polynomial

    G(y)=E_[I_j,partitions_j] sum_(P_j subset I_j) y^|P_j| b^H_j(P_j).

The regions are independent after this relaxation, so positivity and (2)
give the COEFFICIENTWISE bound

    Z_c(b) <= [y^c] G(y)^256.                            (3)

Each coefficient of G retains the exact regional quota through (1). Unlike
an assertion obtained from Holder at positive markers, (3) follows from a
pointwise inequality for each P and actual independence between regions.
The paid coefficient must still be retained; no scalar marker optimization
is claimed sufficient.

At lambda=707/2560,tau=381/512, the worst-case additive cost in (2) is133212.
This is too large to declare the relaxation adequate without evaluating it.
It is a bounded next test, not another completed certificate.

## Checks and remaining task

`quartet_quota_audit.py` now checks the suffix lower bound with directed
logs, plus5632 exact small nested-quota/paid-set cases. Its complete run
took24.22 seconds remotely. `quartet_regional_transfer.py` checks12544
original-quota/paid-set cases against direct enumeration, all paid-count
coefficients, and the unweighted binomial normalization. Width-eight tests
include random quartet partitions and two regions, reaching805 states.
The block upper bound is checked coefficientwise at widths four and eight;
65536 direct masks also check (2). These are exact toy tests using
lambda=1,tau=5/2,b=4, not finite-distance certificates.

The transfer audit took0.34 seconds. All runs were serial, bounded, remote,
and file-free. The frozen defect helper, construction, paper, and library
were not changed. No full-width table was launched.

Next test a bounded full-width envelope for the sixteen-update regional
polynomial G, retaining its r=2100 quota and paid coefficients. First decide
whether (3) can close c=111572; if not, retain boundary state instead of
enlarging a knowingly inadequate relaxation. Other covered counts and row
counts remain obligations for the original end-to-end goal.
