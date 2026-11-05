# One-diagonal feedback: next proof interface

Status: implemented and checked, but not distance-certified. The
[native gate](NATIVE_PERSISTENT.md) measures about 16--17 ms for the
packet-four B64/state128 profile at K20. This note does not change
the certified two-diagonal construction in [Additive MDS](ADDITIVE_MDS.md).
It records a possible reduction in encoder work and the proof obligations
that reduction introduces.

## Proposed recurrence and exact inverse law

Retain the field F, batch width b, globally disjoint additive cosets,
fixed superregular matrices D_t, outer code, and one route from
[Additive MDS](ADDITIVE_MDS.md). Write q=|F| and n=N/b.
Independently sample r_t from F^b for every batch, including zero entries.
The keys are independent of the route, sampled once, and fixed for all
messages. Starting with Z_0=0, compute

    V_t = D_t Z_t,
    Y_t = X_t + r_t .* V_t,
    Z_(t+1) = V_t + Y_t.                                  (1)

Here .* is coordinatewise multiplication in F. Emit all N coordinates
and discard the terminal state. Equation (1) uses one D application and
b variable field products per batch. This is an operation count, not a
timing or a claimed speedup.

For fixed output Y, reconstruct Z and V deterministically from the first
and third equations. Then recover X_t=Y_t+r_t.*V_t. Thus the inner is
F-linear and invertible for every setup. Define

    A(Y) = {(t,j) : V_t[j] != 0},       H(Y)=|A(Y)|.

Over the independent setup keys, the inverse input is uniform on the
affine space Y+F^A(Y). Outside A(Y), the inverse input equals Y exactly.
For a fixed route and its routed outer code C_out, put

    D_out(S)=dim_F(C_out intersect F^S).

The affine-fiber argument gives

    Pr_inner[I^-1(Y) in C_out]
      <= q^(D_out(A(Y))-H(Y)).                             (2)

This is the exact random-coordinate set. A nonzero state does not make
the entire batch uniform: D_t Z_t can have zero coordinates. In particular,
the uniform-batch and q^(-b) reset interfaces of
[Persistent MDS](PERSISTENT_MDS.md) do not transfer to (1).

## Two deterministic bounds from superregularity

For each batch write

    h_t=wt(V_t),       w_t=wt(Y_t),       w=sum_t w_t.

The weights count individual F coordinates, not packets or payload lanes.

### Adjacent support weights

For every live state Z_t!=0, the MDS property of [I;D_t] implies

    wt(Z_t)+h_t >= b+1.

Since Z_t=V_(t-1)+Y_(t-1), for t>=1 with Z_t!=0 this gives

    h_(t-1)+w_(t-1)+h_t >= b+1.                           (3)

No such lower bound is asserted for a zero state. If an output has c
excursions occupying v batches, exactly v-c batches enter with live state.
Summing (3) on those batches yields the conservative bound

    2H(Y)+w >= (b+1)(v-c).                                (4)

The first batch of each excursion has V=0, and each remaining h_t is
counted at most twice. Boundary terms can strengthen (4).

For activity spanning almost the whole stream and few excursions, (4)
provides roughly (N-w)/2 randomized coordinates. At w=0.35N, this is
roughly 0.325N, with a small b-dependent correction. That is insufficient
to reuse the old output-volume proof merely by replacing its random
coordinate count: even a q^w amplitude count can exceed this inverse
penalty before any outer dimension is charged. This identifies a loss
in that proof substitution, not an upper bound on the code's distance.

### Fresh output variables pay deterministic ranks

Fix complete support masks W_t,H_t for Y_t,V_t, respectively. Their
cardinalities are w_t,h_t. Introduce variables only for output coordinates
in W_t, temporarily allowing those variables to be zero. Impose the
homogeneous equations

    V_t[j]=0,       j not in H_t.

The variables in Y_(t-1) have not appeared in any earlier V equation.
Their coefficient matrix in the current equations is

    D_t[H_t^c,W_(t-1)].

Superregularity gives this matrix rank min(b-h_t,w_(t-1)). Earlier rows
are zero on these new variables. Consequently the rank of all imposed
zero equations is at least

    R = sum_(t=1)^(n-1) min(b-h_t,w_(t-1)).                 (5)

The initial V_0 is identically zero; it contributes no rank. The terminal
state is free, so (5) makes no charge after Y_(n-1).

There are therefore at most

    q^(w-R)                                               (6)

output assignments with the specified masks. Required nonzero coordinates
and other state conditions can only reduce this count. Equation (6) counts
amplitudes for fixed masks, not the number of masks. For prescribed weights,
there are at most

    product_t binom(b,w_t) binom(b,h_t)

mask pairs. Actual output and state relations exclude many such pairs.

Combining (2) and (6), a fixed mask pair has conditional existence bound

    min(1, q^(w-R+D_out(A)-sum_t h_t)),
    A={(t,j):j in H_t}.                                   (7)

The rank R in this expression is a deterministic loss of output dimension.
It is not an additional probability of a feedback coefficient vanishing.
Equation (7) is a local counting interface, not a completed sum over masks.

## A stronger bound on a live zero-output run

Use the conceptual syndrome from [Additive MDS](ADDITIVE_MDS.md):

    S_(t+1)=S_t+C_t Y_t,       V_t=C_t^(-1)S_t,

where C_t has Vandermonde columns at the globally distinct points of
the t-th additive coset. Suppose Y_t=0 throughout a run of ell batches,
and its entering syndrome S is nonzero. Then S remains fixed and live
throughout that run.

Let f_S be the unique degree<b polynomial whose values on the predecessor
coset V_h are C_(-1)^(-1)S. This polynomial is nonzero. The symmetric
evaluation-transfer identity proved in ADDITIVE_MDS.md gives

    C_t^(-1)S = (f_S(a_t+beta_j))_(j=0)^(b-1).

Every point encountered during the run is distinct. A nonzero degree<b
polynomial has at most b-1 roots, so

    sum_(t in run) h_t >= b ell-(b-1).                     (8)

The root budget is b-1 for the entire run, not for each batch. A quiet
zero-syndrome run is excluded from this assertion. Nonzero outputs change
S and hence f_S; (8) does not share one root budget across those changes.

## Remaining whole-code gate

The outer input is still confined to the union U of excursion batches.
Its shortening D_out(U), the global MDS reset constraints, and the
excursion representation count remain useful. The obstacle is the finer
set A in (2) and (7): it can omit individual coordinates within a live
batch and within a routed packet.

The exact p-row shortening identity in [Packet route](PACKET_ROUTE.md)
requires complete packets. It cannot be applied to this A without a new
argument. Replacing A by its complete-packet closure is safe for an upper
bound on D_out(A), but loses information. Replacing A by all live batches
also loses the actual number of randomized coordinates. Neither operation
justifies inheriting the old finite margin.

A simple valid fallback is available for every coordinate set:

    D_out(A) <= min(D_out(U), |A|/2).

Indeed each constituent contributes (t-B/2)_+<=t/2 to shortening, and
A is contained in U. This bound needs no packet-alignment assumption,
but its loss has not been evaluated in the support count.

The [delayed-rank refinement](SINGLE_DIAGONAL_RANK.md) now counts constraints
on older output variables, not only the preceding batch. The
[partial-packet moment](SHUFFLED_PACKET_ROUTE.md) also applies after summing
a complete invariant mask family, without adding shuffles to the encoder.
Neither result completes the whole-code count. Retain per-batch mask
multiplicities, activation batches, and the free terminal. The old
two-diagonal reset-probability bound is unavailable. No numerical margin
or whole-code distance claim follows from these local bounds.

The [outer-compatibility refinement](AFFINE_COMPATIBILITY.md) additionally
uses coordinates not randomized by the inverse. Inner constraints on
overlap columns contribute independent ranks to those outer equations.
Its complete mask grouping is valid, and 18 bounded profile diagnostics
now distinguish the useful H- and union-support bounds. A full sum is
still required. The [causal-minor obstruction](CAUSAL_MINOR_OBSTRUCTION.md)
also excludes generic multitime-minor arguments for both b64 and b128;
the valid single-observation and selected-time lemmas are unaffected.
