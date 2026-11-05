# Width two fails the 35% target

This note gives a coefficient-independent obstruction for the width-two
[construction](CONSTRUCTION.md). It uses the actual paired regional route,
not a uniform global permutation or a fixed-message probability bound.

**Claim.** At K=2^20, a route sampled from the specified distribution has,
with probability greater than 0.841, a nonzero message whose output weight
is at most floor(0.35 N). This holds for every choice of outer mixers and
inner coefficients. Width two therefore cannot meet the proposed 35%
distance guarantee with setup-failure probability at most 2^-40.

The argument is analytic. It does not use a sampled route, a numerical
parameter search, or a low-weight message search.

## Fix a message space before sampling the route

Restrict the message to outer groups 0 through 31. These groups contain

    D = 32 * 128 = 4096

free F-valued message coordinates. All other groups are zero. Each of the
256 regions receives one coordinate from every group.

For an even region, the 32 selected groups occupy exactly 16 source pairs.
For an odd region, they occupy exactly 17: the two boundary pairs each
contain one selected group, and the other 15 pairs contain two. A pair is
called occupied here if it contains a selected group. Its actual input
values need not be nonzero.

At width two, each destination pair is exactly one input batch. Write

    M = 8192 / 2 = 4096

for the number of batches per region. The occupied batch positions form a
uniform p-subset of the M positions, with p=16 or p=17 according to region
parity. These subsets are independent across regions. The independent swap
bits change input lanes, but not the occupied positions.

Every unoccupied batch has input identically zero on the entire chosen
message space. This statement holds before choosing the message.

## Force seven empty gaps per region to emit zero

Fix the complete setup. In one region, order its occupied batch positions
and let Z_0,...,Z_p be the numbers of unoccupied batches in the gaps.
Here Z_0 is the region's initial gap, Z_i for 1<=i<p lies after occupied
batch i and before occupied batch i+1, and Z_p is the final gap inside
the region. Thus

    sum_(i=0)^p Z_i = M-p.

Do not use Z_0. Select seven largest gaps among Z_1,...,Z_p, with a fixed
tie-breaking rule. Let G be their total length in complete batches.

For each selected gap, impose that the two-coordinate inner state is zero
immediately after its preceding occupied batch. This gives at most two
homogeneous F-linear equations on the selected message space. All inputs
inside the gap are zero, so its state and outputs remain zero throughout.

Perform this selection in all 256 regions. The combined system has at most

    2 * 7 * 256 = 3584 < 4096 = D

equations. Consequently its kernel contains a nonzero message, for every
fixed setup. Dependence between equations can only enlarge that kernel.
The outer and the causal inner are injective, so this message has a nonzero
codeword.

Let G_r be the selected total in region r. The selected gaps are disjoint
and contain two output records per batch. The constructed message therefore
satisfies

    wt(output) <= N - 2 sum_r G_r
               = N (1 - (1/256) sum_r (G_r/M)).                 (1)

The final selected gap in a region ends at that region's boundary. No
constraint is imposed on a later region, unless it is separately selected.
There is no requirement that the terminal state of the full encoder vanish.
All constraints and the surviving message may depend on the full setup,
as permitted in a minimum-distance event.

## Expected size of the seven selected gaps

Uniform p-subsets of M positions are in bijection with weak compositions
of M-p into p+1 gaps. Hence (Z_0,...,Z_p) is uniform over those compositions.

One way to generate this distribution is to sample

    U = (U_0,...,U_p) from Dirichlet(1,...,1),
    Z conditional on U from Multinomial(M-p; U_0,...,U_p).

Indeed, integrating the multinomial probability gives the same value
for every nonnegative integer vector of total M-p.

The sum of the seven largest coordinates among coordinates 1 through p
is a convex function: it is the maximum of all seven-coordinate sums.
Conditional Jensen's inequality therefore gives

    E[G] >= (M-p) E[sum of the seven largest of U_1,...,U_p].    (2)

The sum S=U_1+...+U_p has mean p/(p+1). Independently of S, the normalized
vector (U_1/S,...,U_p/S) has distribution Dirichlet(1,...,1) on p coordinates.

Write H_j=sum_(i=1)^j 1/i. For a Dirichlet(1,...,1) vector on p coordinates,
the expected sum of its seven largest coordinates is

    (7/p) (1 + H_p - H_7).

For completeness, normalize p independent mean-one exponential variables.
Their sum is independent of the normalized vector and has mean p.
The expected j-th largest exponential is H_p-H_(j-1). Summing these
expectations for j=1,...,7 and dividing by p gives the displayed formula.

Combining this formula with (2) proves

    E[G/M] >= mu_p,
    mu_p = ((M-p)/M) * (7/(p+1)) * (1 + H_p - H_7).            (3)

For both possible values of p, this lower bound exceeds 0.71. The following
weaker rational estimates suffice:

    H_16-H_7 > 39/50,
    H_17-H_7 > 21/25,

    mu_16 > (4080/4096) * (7/17) * (89/50) > 71/100,
    mu_17 > (4079/4096) * (7/18) * (46/25) > 71/100.           (4)

These estimates include the discreteness of the route and exclude the
unusable initial gap. They make no exponential-spacing approximation.

## Probability of a low-weight codeword

The 256 random variables G_r/M lie in [0,1] and are independent, because
their regional pair permutations are independent. Their means all exceed
0.71 by (3)--(4). Hoeffding's inequality yields

    Pr[(1/256) sum_r (G_r/M) <= 0.65]
      <= exp(-2 * 256 * (0.71-0.65)^2)
       = exp(-1152/625)
       < 0.159.                                             (5)

Outside this event, (1) supplies a nonzero message with output weight
strictly below 0.35 N, and hence at most floor(0.35 N). The resulting
failure probability exceeds 0.841.

Only the route was random in (5). The conclusion holds after conditioning
on arbitrary outer mixers and inner coefficients. Averaging over their
setup distributions preserves the same bound.

## Consequence for the design

The obstruction is not a failure of coefficient mixing or a loose lower
bound. The 32-group message space has enough dimensions to reset a
two-coordinate state before many large input-free gaps. Paired routing
places two selected group coordinates in one batch, reducing the number
of distinct input batches and making this attack more effective.

Increasing coefficient entropy, changing the field, or improving the
width-two local rank calculation cannot remove this obstruction. A new
candidate must change the state dimension, the outer/route locality, or
another assumption used in the reset construction. Passing this particular
test at a larger width would still not establish a distance lower bound.
