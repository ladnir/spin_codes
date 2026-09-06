# Removing the parity cost without changing the encoder

The remaining question is a probability question, not a count question.
An enormous expected number of bad messages need not imply frequent setup
failure. This note defines a weighted count whose second moment is simpler.
It then selects existing message rows to enforce column parity.

These changes affect the proof witness only. The fixed BCH [256,128] outer,
RM2Sub t128_s15 inner, and setup distribution remain unchanged.
The full second moment is not yet bounded.

## The experiment and the zero-state event

Let C be the fixed BCH subcode defined by p37(c) in {0,...,31}.
It is an even [256,128] code containing the all-one vector, with minimum
distance at least 38. Its dual D has minimum distance at least 30.
Let T80 be the nonempty set of weight-80 words of C.

Write the setup as (rho,sigma,alpha). The 8192 row permutations rho_i
are independent uniform permutations of 256 coordinates. The 256 region
permutations sigma_j are independent uniform permutations of 8192 positions.
The nonzero state multipliers alpha are independent of these permutations.

Each region contains 64 epochs of 128 bits. Write B for the fixed
15-by-128 update map. If every epoch input X satisfies BX=0, the state
remains zero and the output equals the input. Call this event E_u for
candidate message u. It does not depend on alpha.

For an input region of weight x, let beta_x be the probability that all
64 epochs lie in ker B after its uniform region permutation. For two
regions with weights x,y and overlap k, let beta2(x,y,k) be the probability
that both pass under the same region permutation. These are exact
fixed-weight probabilities, not independent-bit replacements.

Fix rho. A finite candidate index set U may specify a message u(rho)
and a retention indicator G_u(rho). Neither may depend on sigma or alpha.
For every retained candidate, require

    p_u(rho) := product_j beta_{J_j(u(rho))} > 0,

where J_j is its input weight in region j. Define the nonnegative witness

    F(rho,sigma) := sum_{u in U: G_u(rho)=1} 1[E_u] / p_u(rho).

The positivity event F>0 means that at least one retained candidate has
zero state throughout. In the family below, every such candidate is a
nonzero bad message for the original code.

## Exact cancellation in the two moments

Conditional on rho, independence of the region permutations gives

    E_sigma F = sum_u G_u(rho),

    E_sigma F^2
      = sum_{u,v} G_u(rho) G_v(rho)
          product_j R(J_j(u), J_j(v), K_j(u,v)),

where K_j is the overlap of the two region inputs and

    R(x,y,k) := beta2(x,y,k) / (beta_x beta_y).

The second sum includes diagonal pairs. Both messages use the same rho
and sigma. Only different regions have independent sigma_j.
The factors beta_x beta_y cancel exactly; no approximation of their
variation across rho is needed.

For any finite nonnegative F, Cauchy--Schwarz gives

    Pr[F>0] >= (E F)^2 / E F^2.

`test_inverse_probability_witness.py` checks these identities in an
exhaustive toy model with shared permutations and nonconstant p_u.
The toy test supplements the derivation; it does not prove the BCH bound.

## Three existing rows enforce every column parity

Reserve zero-indexed row positions 8189,8190,8191. A core index specifies
2610 occupied positions among 0,...,8188 and one T80 word at each position.
All remaining core rows are zero. Thus

    M := |U| = binom(8189,2610) |T80|^2610.

For fixed reserved-row permutations, consider the linear map

    L_rho : C^3 -> V_even,
    (c1,c2,c3) |-> rho_8189(c1)+rho_8190(c2)+rho_8191(c3),

where V_even is the 255-dimensional space of even vectors in F2^256.
Let H be the event that L_rho is onto. On H, choose its right inverse
by fixed-order binary Gaussian elimination. Apply that right inverse to
the XOR of all permuted core rows, and place the resulting C words in
the three reserved positions. The XOR of all completed rows is then zero.

This selection is deterministic after rho is fixed. It does not depend
on the region permutations. Distinct core indices produce distinct
completed messages because their core rows remain unchanged.
Every completed message has occupancy between 2610 and 2613 and weight
at most

    2610*80 + 3*256 = 209568 < 209716.

Consequently, E_u implies a bad output for the fixed target. Selecting a
message after setup is permitted here because setup failure means that
there exists a bad nonzero message. No encoder step has been added.

### Probability that the correction map exists

The orthogonal complement of the image of L_rho is the intersection of
three independently permuted copies of D. Each contains the all-one vector.
The map fails to be onto precisely when this intersection has dimension
at least two. Such an intersection contains at least two vectors other
than zero and one.

Let B_w be the number of weight-w words of D. For a fixed vector of
weight w, membership in a uniformly permuted D has probability
B_w/binom(256,w). Hence

    Pr[not H] <= (1/2) sum_{w=1}^{255} B_w^3 / binom(256,w)^2.

D has orthogonal-array strength 37 because C has minimum distance at
least 38. Define K_j(w) as the binary Krawtchouk polynomial of degree j
at length 256. The degree-18 Christoffel bound gives

    B_w <= floor(2^128 / sum_{j=0}^{18} K_j(w)^2/binom(256,j)).

Indeed, the square of the degree-18 reproducing kernel has degree 36.
Its average on D equals its binomial average. Bounding the contribution
of weight w alone gives the displayed inequality. Also use the ambient
cap binom(256,w), minimum distance 30, evenness, and complement symmetry.

`certify_parity_corrected_witness.py` evaluates the resulting rational sum
exactly. It is below 2^-100 (diagnostic 1.812e-31). The accompanying sample
right inverse is checked on a basis of V_even; that sample is not used
as probability evidence.

## A first moment with no global-parity penalty

Retain a core index when H holds and every core region weight lies in
[740,897]. The three correction rows add between zero and three bits
per region. Every completed region therefore has even weight in [740,900].
All beta_x in this window are positive.

For any fixed core index, each individual core region weight has law
Bin(2610,5/16) over the independent core row permutations. Region weights
need not be independent. A union bound gives

    g := 1 - 256 Pr[Bin(2610,5/16) outside [740,897]].

H depends only on reserved-row permutations, independently of core rows.
The exact cancellation above therefore yields

    E F / M >= (1-Pr[not H]) g > 3/4.

The stored lower bound is approximately 0.7703919913290518; its exact
rational representation, rather than this decimal, is authoritative.
The producer uses 256-bit intervals; the 512-bit replay passed.
There is no 2^-255 parity loss in this first moment.

## What the pair certificates cover

Nine outward Fourier tiles cover all 793881 ordered types with even
x,y in [740,900] and integer k in [40,160]. They prove

    R(x,y,k) <= exp((3/25000) (k-xy/8192)^2).

The tiles use coefficient aliases with nonnegative coefficients, exact
dyadic tilts, and explicit outward errors for evaluation and radix-2
Fourier butterflies. Pair-swap symmetry covers the transposed cross tile.
`verify_overlap_atlas.py` checks the exact domain ledger; `--replay`
also repeats every tile with 512-bit root and final arithmetic.
Consult the replay receipt for completion status.

For independent uniform core indices u,v, the normalized second moment is

    E F^2/M^2
      = E_{rho,u,v}[1[H] Gcore_u Gcore_v product_j R(J_j(u),J_j(v),K_j)].

To refute the desired setup-failure bound, it would suffice to prove
this expectation below (9/16)*2^40. That would imply Pr[F>0]>2^-40.
It would not establish a positive all-occupancy guarantee.

The unproved obligation includes overlaps k<40 and k>160, diagonal pairs,
and correlations across all 256 regions. Correction rows can increase
the core overlap by zero through three. Their dependence on core parity
must be retained or covered by a uniform bound. Multiplying a per-type
maximum only over the certified window does not bound the full expectation.

## Retained alternative witness

`zero_state_polynomial_witness.py` improves the earlier, uncorrected T80
central-count witness using a degree-six lower envelope for log beta_x.
Exact additive moments and the existing weighted-parity bound give
an approximate log2(E Z_G) lower bound of 19633.21535078605, with the stored
rational bound authoritative.
Its 512-bit replay passed, gaining about 8.53 bits over the quadratic bound.
The corrected inverse-weight witness avoids this single-probability
approximation, but the earlier certificate remains valid and unchanged.
