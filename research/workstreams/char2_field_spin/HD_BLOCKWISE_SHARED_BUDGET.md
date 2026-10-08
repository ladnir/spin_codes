# A state budget for tile-local mixing

Local mixing removes the identity-feedthrough obstruction, but it does not
permit arbitrarily sparse state updates. Even a perfect eight-point MDS map
can only spread a quiet input inside its own eight-point tile. This gives a
new distance upper bound for every regional route, not just a random route.

At K20, eight-point tiles and width64/state1 cannot reach 40% distance:
their distance is at most22.439%. Width64/state2 has a40.342% ceiling, leaving
almost no room above the target. The active width64/state4 choice survives
this obstruction. That is a necessary gate, not a distance certificate.

## Construction covered by the bound

Let the outer be R independent injective linear maps from k to b field
symbols. Thus K=Rk and N=Rb. Each coordinate region contains one coordinate
from every row. Fix an arbitrary permutation within each region.

Partition every region into tiles of g positions. Require g to divide R
and the packet width t, and t to divide N. Each packet has an s-symbol
state and the linear recurrence

    Y_i = F_i X_i + A_i Q_i,
    Q_(i+1) = M_i Q_i + C_i X_i,       Q_0=0.

Each F_i is invertible and block diagonal on the g-position tiles. The
other maps may be arbitrary and may vary with i. There is no state flush;
the final state is free. In particular, neither MDS local maps nor random
mixers are assumptions of this bound.

## A quiet word with a small covered support

There are P=N/t packets. Define

    e = s(P-1),
    a = floor(e/k)+1,
    d = ak-e.

Assume a<=R. For each set J of a outer rows, restrict the message to J and
impose C_i X_i=0 in every preterminal packet. At most e equations leave a
subspace of dimension at least d>0. Every entering state is zero on this
subspace, so Y_i=F_i X_i.

Call a tile covered when it contains at least one routed coordinate from J.
All outputs of the quiet subspace lie in covered tiles. Fix the route and
average over the uniform choice of J. Every tile contains coordinates from
g distinct source rows. Its probability of being uncovered is

    C(R-g,a)/C(R,a) = C(R-a,g)/C(R,g).

The mean number of covered output positions is therefore

    U = N*(1-C(R-a,g)/C(R,g)).

Some J has at most floor(U) covered positions. Its quiet subcode has
dimension at least d, and the causal encoder is injective because every
F_i is invertible. The Singleton bound on that subcode gives

    minimum_distance <= floor(U)-d+1.

This proof averages a witness set after fixing the route. It does not
average code distance over sampled routes and has no failure probability.

One may use any larger a with ak>e. The expected covered support minus ak
is a concave sequence in a. Its minimum occurs at an endpoint. The endpoint
a=R gives N-K+e+1, weaker than the ordinary Singleton bound N-K+1. Thus
the smallest admissible a above suffices for this particular obstruction.

## Exact K20 consequences

Here the outer is [256,128], R=8192, and N=2097152. The reported quiet caps
come from integer arithmetic before conversion to percentages.

| Packet width t | State s | Tile g | Quiet weight cap | Relative cap |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 4 | 1 | 262141 | 12.4999% |
| 64 | 1 | 8 | 470570 | 22.4385% |
| 64 | 2 | 8 | 846024 | 40.3416% |
| 64 | 4 | 4 | 867965 | 41.3878% |
| 64 | 4 | 8 | 1376900 | 65.6557% |
| 128 | 4 | 8 | 846022 | 40.3415% |
| 256 | 4 | 8 | 470567 | 22.4384% |

The ordinary50.00005% Singleton cap is stronger than the active row's
65.6557% quiet cap. That row is not evidence of a65% achievable distance.

The budget couples three construction choices: tile mixing width, state
dimension, and update frequency. Reducing state to one symbol, or retaining
four state symbols but updating only once per256 inputs, is structurally
incompatible with40% distance at tile width8. Width64/state2 and
width128/state4 are marginal candidates, not harmless implementation changes.
The current width64/state4 design keeps more room for a finite proof.

## Checks and scope

[The helper](theory/hd_blockwise_shared_budget.py) uses exact rational
arithmetic. [Six tests](theory/hd_blockwise_shared_budget_tests.py) enumerate
every row subset in small fixed, correlated routes; verify the mean coverage
and concavity; and check the parameter gates and rounding. All six pass.

This result is an upper bound on distance, not a bound on low-weight
codeword multiplicity. It neither proves a full certificate for the active
candidate nor excludes improvements from a different tile map or route.
The bound does not apply unchanged when a tile combines repeated coordinates
from the same source row, or when feedthrough crosses the declared tiles.
