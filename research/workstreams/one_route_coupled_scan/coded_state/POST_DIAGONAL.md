# Move the random diagonal after expansion

This alternative addresses the two-tap proof's main obstruction: sparse
state supports produce low-dimensional, non-coordinate inverse spaces.
The [sparse proof](SPARSE_PROOF.md) now covers every message using at most31
outer blocks at >49.1% distance, with failure <2^-50. The
[regional episode bound](ROUTE_EPISODES.md) extends coverage to64 blocks;
the combined failure is <2^-46. These are complete source classes, not a
complete code certificate; larger supports remain open.
The separate [column-scaled variant](QUIET_BRANCH.md) preserves these
certificates and closes the entirely quiet branch at >49.1%/66 bits.
It adds one GF256 scaling per output; mixed live/quiet trajectories remain open.
Keep the same outer, regional route, A, and B. Change only the inner:

    Y_t = X_t + D_t A Q_{t-1},
    Q_t = Q_{t-1} + B Y_t,       Q_0=0.

Every entry of every D_t is sampled independently from all of F during
setup, then fixed for all messages. In particular, zero diagonal entries
are allowed. Final state is free. This is not the same distribution as
the one- or two-tap candidates, which remain preserved.

It uses128 full-field products per128 outputs, the same count as two taps.
It needs only the original A and B transforms: no second64-point inverse
transform is needed for feedback. Raw diagonal setup is32 MiB at K20.
These are operation/storage counts, not native performance measurements.

## Exact uniform state refresh from every nonzero state

Fix X and a nonzero incoming Q before exposing the current diagonal.
Set S=supp(AQ). Since im(A) is [128,64,65] MDS, |S|>=65. The vector
D A Q is uniform on the coordinate space F^S. Every64 columns of B are
independent, so B restricted to S is onto F^64. Consequently

    Q' = Q + BX + B D A Q

is exactly uniform in F^64. This holds for every nonzero Q, not only
full-support Q. It also holds for any fixed input X.

Unlike the two-tap design, an empty input can now reset a live state.
The exact probability is q^-64. A zero output still cannot reset a live
state because Q'=Q+B Y. The new randomness gains uniform refresh by
giving up deterministic empty-input persistence. A proof must count resets;
it cannot simply declare their fixed-input probability negligible before
taking the all-message union.

## Coordinate inverse law

For any fixed candidate output y, reconstruct Q from Q'=Q+B y independently
of setup. Let S(y) be the union of the coordinate supports of A Q_{t-1} in
their respective batches. The inverse input is exactly uniform in the
affine coordinate space y+F^S. Thus, for every fixed route,

    Pr_D[inverse(y) in C] <= q^(D_C(S)-|S|),

where D_C(S) is the routed outer's shortened dimension on S. In this
construction, D_C(S)=sum_i (t_i-128)_+, with t_i the number of selected
coordinates from outer block i. This is precisely the coordinate
shortening object used by the retained regional-moment proof.

This removes the non-coordinate intersection problem of the two-tap design.
It does not by itself give the enumeration of all possible S(y). Sets S
can have holes in every batch, so interval counting alone is insufficient.

## Two global classes that close immediately

First, put H(y)=N-|S(y)|. Using D_C(S)<=K gives
Pr[inverse(y) in C]<=q^(-K+H). There are at most2^N*q^W words of weight
at most W. At K=2^20, W=floor(.49N), all outputs with H<=4587 therefore
have total setup-failure probability at most2^-128. The first batch already
contributes128 holes. This is a class certificate, not full minimum distance.

Second, long low-dimensional state trajectories need not be controlled by
this hole budget. For any S, (t_i-128)_+<=t_i/2, so
D_C(S)<=floor(|S|/2), uniformly over the route. Suppose the incoming state
is nonzero at every batch after the first, and every output batch has
weight at most31. With n=N/128=16384,

    |S| >=65(n-1),
    #eligible outputs <=2^N*q^(31n).

The total failure bound is at most

    2^N*q^(31n-ceil(65(n-1)/2)) = 2^-1044480.

This covers some trajectories with enormous hole budgets. In particular,
it does not assume the state becomes generic. It does not cover arbitrary
zero-state intervals or batches with higher output weight. Both class
bounds are unions over fixed output candidates and therefore cover messages
chosen after setup.

## A remaining map alignment to examine, not silently ignore

The current product P=AB has strong structure. If the two halves of Y are
viewed as degree-<64 evaluations, its diagonal64-by64 blocks are
(1+alpha)I and alpha I. Its off-diagonal blocks are scaled extrapolation
matrices, which have no zero entries. Hence P has8064 zero entries.

For example, if outputs use only one coordinate j in the first half,
then Q is always a scalar multiple of B[:,j]. Its expanded state is a
multiple of P[:,j], with weight65 when that scalar is nonzero. Missing
expanded-state coordinates can persist for many batches without requiring
fresh zero equations. This is an output trajectory, not a codeword witness;
the preceding sparse-output bound already covers its full-scan live case.
It does refute an attempted proof that charges every such hole independently.

One possible countermeasure is to replace B by R B for a fixed invertible
diagonal R over GF256. This preserves B's MDS property and the exact refresh
and inverse proofs above, although BA=I is no longer claimed. It costs64
small-field scalings per batch. One fixed test choice reduces the zero-entry
count of A R B from8064 to65. That count is not a bound on its minors and is
not a security certificate; no scaled map has been selected or promoted.
Alternatively, retaining the original alignment may make a sharp structured
enumerator easier. Do not assume fewer zero entries necessarily means a
better end-to-end construction.

## Checks and next decision

`check_post_diagonal.py` exhausts65535 nonzero GF16 states, covering all93
resulting expanded supports, and verifies full compression rank on each.
It enumerates the next-state law for48 GF4 input/state fixtures and checks
the inverse. Over GF256 it constructs an expanded state of exact weight65,
verifies rank64 on that support, and checks32 recurrence steps.
The two class bounds use integer arithmetic only. No finite whole-code
distance certificate or native timing is claimed.

Next compare proof envelopes for this coordinate-support design against
two taps. Retain the same cheap maps initially. The main obligations are
to count zero-state intervals and correlated expanded-state holes, then
combine that count with actual regional shortened-dimension moments.
Only add the fixed row scaling if it removes a demonstrated global loss.
