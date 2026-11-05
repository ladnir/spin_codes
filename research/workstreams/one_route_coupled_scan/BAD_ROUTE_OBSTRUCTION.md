# Why the whole-code proof must retain route randomness

The deterministic regional rank cap is useful, but cannot replace the
distribution of routes in a whole-code certificate. Some valid routes
admit low-weight words for every choice of inner coefficients. The
construction below uses only linear dimension and bounded state size.
It is not a claim that a sampled route has a substantial failure probability.

## Parameters and the selected route

Let the outer consist of pG independent injective [B,k] constituent maps,
grouped into G groups of p rows. Its input dimension is pGk and its
output length is N=pGB. Region j contains position j from every row.
Each region independently permutes the G packets of p coordinates.

Let the inner process b coordinates per batch, with an s-dimensional
state over F. Assume that the inner is linear and invertible for every
setup, and that zero input and zero state produce zero output and zero
next state. These properties hold for the one-diagonal candidate, with
s=b. Assume that pG is divisible by b, so region boundaries align with
inner batches.

Choose g of the G outer groups. Consider any valid route placing their
packets in the first g packet positions of every region. The identity
route is one example when those groups are the first g groups. Their
relative packet order within that prefix does not matter.

Restrict the message to these groups. Its remaining dimension is gpk.
The resulting outer input to the inner is zero outside the first gp
coordinates of each region. Put

    ell=ceil(gp/b),

and assume ell*b<=pG. The possible nonzero inputs lie within the first
ell batches of each region.

## Reset constraints produce a low-weight codeword

Fix any realization of all inner coefficients. Impose that the inner
state is zero immediately after batch ell of each region, using
one-based local batch numbering. These are homogeneous linear equations
on the gpk selected message coordinates. Each region contributes at
most s equations, so their total rank is at most Bs.

If

    gpk > Bs,                                             (1)

a nonzero selected message satisfies every reset condition. The initial
state is zero. After each required reset, both input and state remain
zero through the rest of that region. Thus the next region also starts
at zero state. Induction over the regions shows that the output is
supported in at most

    B*b*ceil(gp/b)                                        (2)

coordinates. The output is nonzero because both the outer and inner
are injective. Therefore (2) is a deterministic upper bound on the code's
minimum distance for every inner setup under this route.

The message in this argument may depend on the realized inner setup.
That is the correct quantifier for a minimum-distance upper bound.
No claim is made about one message fixed before setup. Imposing a reset
in the final region is permitted even though the encoder's terminal
state is free; it only selects a subspace of messages for this example.

## The current packet-four parameters

Take the candidate parameters

    K=2^20,  N=2^21,  B=64,  k=32,
    p=4,     b=s=128, G=8192, g=65.

The selected message dimension and maximum constraint rank are

    gpk=65*4*32=8320,
    Bs=64*128=8192.

At least 128 dimensions remain after imposing all resets. The gp=260
possibly nonzero input coordinates per region occupy ell=3 batches.
Equation (2) gives

    d_min <= 64*3*128 = 24576 = (3/256)*N.

This is at most 1.171875% relative distance, far below the 35% target.
The argument does not use special diagonal values, a kernel of an outer
constituent, or the causal-minor obstruction in the separate note.

## Scope and consequence

For a specified set of g groups, the event that every region places
exactly those groups in its first g packet slots has probability

    binom(G,g)^(-B)

under the independent uniform regional permutations. Arbitrary packet
orders inside and outside the prefix are included in this event. This
formula describes one particular bad-route event; it is not the total
failure probability of the ensemble. It gives no practical contradiction
to a 40-bit certificate that averages over routes.

The deterministic cap

    rank(P_T) >= sum of the B-k largest regional occupancies of T

remains correct for the MDS outer; B-k=k at rate one half. What fails is
the goal of proving the target distance
uniformly over every allowed route. A complete certificate must retain
route averaging or first establish a suitable good-route event. The
guard D_out(U)>0 and the joint regional moment for T subset U are suitable
interfaces for retaining this distinction.

This note is an analytic obstruction. It requires no benchmark, search,
or numerical regression. Next combine the guarded route moment with a
complete support or excursion count, using deterministic regional caps
only as complementary bounds.
