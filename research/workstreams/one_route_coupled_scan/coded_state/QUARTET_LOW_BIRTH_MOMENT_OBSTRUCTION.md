# A low-birth obstruction to the current relaxed moment

2026-10-31. Exact forced-birth truncation alone cannot repair the current
joint support/routing moment. An explicit family with exactly four forced
births makes that moment bound exceed 2^16130 at the recorded marker choice.
The family violates the target eligible-row and footprint classes. This is
an obstruction to this bound, not a low-weight codeword or a failure of the
fixed construction. The full 45%/40-bit goal remains open.

The [joint actual-state interface](QUARTET_ACTUAL_MDS_JOINT.md) defines the
moment. Here n=512, T=4096, N=2^21, W=943718 and q=2^32. There are 256
regions of 16 updates each. Fix a candidate set H of h=4096 outer rows.
Its labels form an independent uniform 4096-subset of the 8192 slots in
each region. Quartet partitions are independent and uniform at each update.

For candidate state supports S_0=empty,S_1,...,S_T, use

    V_t = Gamma_t(S_(t-1)) union S_t,
    Ls = sum_t |S_t|,       D = sum_t |S_(t-1) minus S_t|,
    f = sum_t |S_t minus Gamma_t(S_(t-1))|,
    G = sum_t |Gamma_t(S_(t-1))|,       F=G+f.

The exact quota-conditioned relaxed moment sums

    s^Ls t^D z^G nu^f eta^Cout,                         (1)

where Cout counts covered positions outside the candidate H. It retains
the necessary MDS filter: an old-empty quartet requires every new state
position to lie in H; an old-live quartet of size r permits at most r-1
new zeros outside H. The final state is free. Actual eligible-row
classification and the previously proved high-F/high-R exclusions are
not retained in this moment.

Fix the tested markers

    s=2, z=1, nu=2^-18, eta=1/4,
    x=2^-16, theta=1/16, t=65537^(-1/16).                (2)

Only t>=1/4 is needed below; it follows from 65537<4^16. Thus the proof
also applies when t is rounded upward to 1/2 and s is rounded upward.

## The zero-to-zero diagnostic and its limitation

Let Z be the matrix with its only nonzero entry Z[0,0]=1. The regional
comparison has the form Z+Qmin, where
Qmin is the entrywise minimum of the positive quota-adjusted regional
operators. Its (0,0) entry counts nonempty regional excursions, so

    e_0 (Z+Qmin)^256 1 >= (1+Qmin[0,0])^256.             (3)

The root-run diagnostic reported log2 Qmin[0,0] about 2029.5566. Formula
(3) explains a large part of the positive comparison bound. This decimal
is a floating diagnostic; the exact argument below does not use it.

A genuine nonempty region that starts with an empty state requires at
least one forced birth. Repeating 256 such excursions requires f>=256.
An exact f=4 restriction removes that particular pattern. It does not
remove the persistent family constructed next.

More generally, define a state episode as a maximal interval with
S_t nonempty. Its first update has empty incoming state, so every state
position at that update is a forced birth. The episode count, and the
sum of its initial support sizes, are therefore at most f. Nonempty
regions with empty incoming state contain distinct episode starts.

## A persistent transition with no forced births

Fix any old support of size 448. Retain partitions with 64 quartets of
old size three and 64 quartets of old size four. Their probability is

    P = C(128,64) 4^64 / C(512,64).                     (4)

In each old type, choose 32 quartets to have new size four and the other
32 to have new size three. The total new size is again 448. Every quartet
is hit, so G gains 512 and f gains zero. Old and new sizes sum to at least
six, making the MDS zero filter automatic for every H-label allocation.

For one quartet, the weighted new-mask sums for old/new sizes
(3,3),(3,4),(4,3),(4,4) are respectively

    (1+3t)s^3,       s^4,       4t s^3,       s^4.

The selected new masks therefore contribute at least

    s^448 C(64,32)^2 [(1+3t)(4t)]^32.

Throughout a 16-update region of these transitions, every slot is
covered. The exact regional quota makes Cout=4096, so its label weight
is exactly eta^4096=2^-8192. Allocate a factor 2^-512 per update and
use s=2,t>=1/4. The resulting per-update contribution is at least

    B = 2^64 C(128,64) C(64,32)^2 (7/4)^32 / C(512,64)
      > 2^34.                                           (5)

All comparisons in (5) reduce to integers. Central binomial coefficients
exceed 2^128/129 and 2^64/65. The binomial theorem and

    C(512,64) 7^448 < 8^512,
    4*8^7 < 11*7^7

give C(512,64)<22^64. Consequently

    B > 2^192 7^32 / (129*65^2*11^64) > 2^34,

using 129*65^2<2^22, 11^2<2^7 and 7^4>2^11. Thus the exact
quota-conditioned regional transition 448->448, restricted to zero
forced births, exceeds 2^(16*34)=2^544. This lower bound does not rely
on the entrywise comparison over Bernoulli tilts.

## Enter the persistent family with exactly four births

Use the following support sizes in the first six updates:

    4, 16, 64, 256, 512, 448.

Fix the first four state positions. Choose partitions in updates 2--4
that place the old positions in distinct quartets, and take the entire
hit set as the new support. In update 5, choose a partition with two old
positions in every quartet and take the full new support. In update 6,
retain four positions in 64 quartets and three positions in the other
64. This loses exactly 64 state positions and satisfies the MDS filter
for every label allocation. All subsequent updates use the family above.

The only forced births are the first four positions. Five specified
partitions suffice. Each has probability greater than 2^-4608, because
the number of quartet partitions is less than 512^512. The first-region
quota places the four initial positions in H with probability

    (4096)_4/(8192)_4 > (1/4)^4 = 2^-8.

For this region, Cout<=4096, so eta^Cout>=2^-8192. The initial birth
weight is nu^4=2^-72; the 64 initial deaths cost at most 2^128. Ignoring
the other positive state weights in the prefix only lowers the count.
The remaining ten updates in the first region each contribute more than
2^(512+34)=2^546 before its label weight. Hence its 0->448 transition,
restricted to exactly four births, contributes more than

    2^[-5*4608-8-8192-72-128+10*546] = 2^-25980.         (6)

The next 255 regions contribute more than 2^(255*544). Combining (5)--(6),
the term with exactly f=4 in the moment evaluated at nu=2^-18 exceeds

    2^[-25980+255*544] = 2^112740.                       (7)

The bare coefficient [nu^4] is 72 bits larger. Distinguish that coefficient
from the evaluated term in (7) when implementing exact-f extraction.

## Consequence for the recorded profile bound

The tested profile is

    (h,R,Ls,f,F) = (4096,786432,786432,4,1048576).

Its joint-moment prefactor, with the existing gamma<2 upper bound, is

    C(8192,4096) 2^(W-256h-31/16+18f)
      * (65537^(1/16)/2)^Ls.                            (8)

There is no separate routing factor or external quota divisor. Since
C(8192,4096)>2^8192/8193>2^8178 and 65537>2^16, (8) exceeds 2^-96610.
Equations (7)--(8) therefore make this marker envelope exceed 2^16130,
even with exact regional quotas and exact restriction to f=4. The same
conclusion holds after removing the all-zero history or retaining at
least 129 occupied regions: this family occupies every region.

The family's footprint is

    F = 4+16+64+256+4092*512 = 2095444.

Every slot in its last 255 regions is covered. Thus every outer row has
at least 255 covered positions: its actual eligible set would contain
all 8192 rows, and its retained outer rank would be K=2^20. It violates
the target h=4096 class, the outside-row capacity for that class, and the
previously excluded high-footprint and high-rank ranges. No field
amplitudes, compatible outer message, or low-weight output are exhibited.

## Full row classification adds negligible gain at the exact target

The preceding family enters because the broad moment drops the target
statistics. Once the target (h,R,F)=(4096,786432,1048576) is imposed,
retaining every individual eligible-row degree would add almost nothing
to the present geometry count. This is a separate exact statement.

Fix the support history, partitions, and the H-labelled slot set in every
region. Only row identities within H and its complement remain random.
Those identities are independently permuted between regions. Let c_j
and e_j count covered H slots and covered outside slots in region j.
At this exact profile,

    sum_j c_j = C_H = 192h,
    sum_j e_j = C_out = 64(8192-h).

For one fixed H row, its covered indicators in the 256 regions are
independent Bernoulli variables with probabilities c_j/h. If Z_H is its
degree, the arithmetic-geometric mean inequality gives

    E[3^(-Z_H)] = product_j (1-2c_j/(3h)) <= 2^-256,
    Pr[Z_H<=128] <= (3/4)^128.

For an outside row Z_O, the same argument gives

    E[3^Z_O] = product_j (1+2e_j/(8192-h)) <= (3/2)^256,
    Pr[Z_O>=129] <= Pr[Z_O>=128] <= (3/4)^128.

No independence between different rows is needed. A union bound over all
8192 rows shows that H is the actual eligible set with conditional
probability at least 1-epsilon, uniformly over every fixed slot geometry,
where

    epsilon = 8192(3/4)^128 < 2^-40.

The last comparison is the already checked integer inequality
9^64<2^203. Consequently inserting actual row classification into any
nonnegative exact-profile geometry count whose weights depend only on
the supports, partitions and H labels retains at least a 1-epsilon
fraction of that count. Its possible improvement is less than 2^-39
bits: use -ln(1-epsilon)<=epsilon/(1-epsilon) and ln2>2/3.

This comparison does not apply after conditioning on additional RS
amplitude equations or a bad-codeword event. It limits the gain from
adding the row-degree filter to the current weighted geometry count.
The independent audit checked this uniform conditioning argument as
well as the persistent-family proof above.

## Consequence for the next proof step

The exact-f interface remains valid: retain the positive birth variable
through regional quota conditioning, truncate polynomial coefficients at
f_0, and compose coefficients before taking the final f coefficient or
sum. It would remove repeated independently restarted excursions. The
family above proves that this truncation alone is insufficient at (2).

Concretely, retain coefficients of degrees 0,...,f_0 in the regional
matrix polynomial R(y). For each degree, quota-adjust the 16-update
coefficient and then minimize its entries over the Bernoulli tilts in
the common support-count basis. Preserve the exact empty operator in
degree zero. The desired upper count is

    [y^f] e_0 R(y)^256 1,

or its sum over 1<=f<=f_0. Positive polynomial multiplication permits
discarding higher degrees after every multiplication. This needs f_0+1
coefficient matrices rather than a dense matrix on all birth counts.

The next bound must retain the excluded footprint/rank/row-class ranges
inside the positive sum, or change the field and coverage markers enough
to suppress them. One cannot remove their contribution merely by taking
the minimum with their separately proved failure probabilities.

This note uses derivation and integer inequalities only. An independent
audit checked the stationary count, prefix, quota conditioning, endpoint,
and all exponents. No numerical jobs, benchmark, or construction change
accompanied the proof. The joint-count note records subsequent stronger
field controls, including explicitly floored FLOAT comparisons. Their
positive values do not prove that the corresponding exact or unfloored
bounds fail; the obstruction proved here has the precise markers (2).
