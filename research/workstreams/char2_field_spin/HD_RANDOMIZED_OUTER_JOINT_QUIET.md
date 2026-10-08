# A joint quiet/output moment for the randomized outer

The joint moment below certifies every quiet low-weight word with
5940 through 8192 active outer rows. The total failure probability for this
class is below **2^-3832.4504**. This extends the earlier quiet dense range,
which began at 6566 rows. The proof checks two occupancy endpoints and
uses convexity to bound every intervening integer; it does not scan that range.

Together with the [randomized-outer proof](HD_RANDOMIZED_OUTER_PROOF.md),
the remaining quiet class has 954 through 5939 active rows. The subsequent
[weighted-state proof](HD_RANDOMIZED_OUTER_STATE_PROFILE.md) also closes
all one-active-row words. Nonquiet words with 2 through 7547 active rows
remain open. No complete distance certificate is claimed.

## The local quantity

Use the setup in the randomized-outer proof: q=65536, R=8192 outer rows,
256 independently permuted regions, 32768 packets of 64 positions, and
the shared beta64 maps F and C. Each packet comprises eight FFT8 tiles.
The output-weight limit is D=838860.

Fix one packet's potential input subset S, before drawing its artificial
uniform amplitudes in the outer comparison. Let X_S be uniform over
F_q^S and set all other input coordinates to zero. For 0<z<=1, define

    K_S(z) = E[1[CX=0] z^wt(FX)].

This expectation is joint. It does not separate the output weight from the
feedback syndrome. For a uniform j-subset S of the packet, let K_j(z)
be the average of K_S(z). The free final packet uses

    T_j(z) = E[z^wt(FX)],

with the same potential subset and amplitude laws, but no syndrome constraint.

For any output vector y of length 64,

    z^wt(y) = sum_(V subset [64])
               z^(64-|V|) (1-z)^|V| 1[y_V=0].

Every summand is nonnegative. We can therefore bound each simultaneous
output-zero and feedback-zero event before summing.

## Exact tile ranks

Fix one tile with input evaluation nodes alpha+{0,...,7} and output nodes
alpha+64+{0,...,7}. These sets are disjoint. Let a be the number of
potential input positions, and v the number of selected output-zero positions.

In the interpolating polynomial p of degree below eight, the 8-a excluded
input positions and v selected output positions are distinct roots. Thus
the output constraints have rank min(a,v). Their residual space has dimension

    d=max(a-v,0).

If d>0, write p=R_0 h, where R_0 is the monic polynomial with these 8-d
roots and deg(h)<d. The tile feedback consists of p's four highest
novel-basis coefficients. It vanishes exactly when deg(p)<4. The feedback
kernel therefore has dimension max(d-4,0), so its image rank is min(4,d).
The statement also holds at d=0.

The calculation is exact for every choice of masks, not just their sizes.
Degree-preserving changes between monomial and novel bases preserve the
argument. The ranks remain valid after extension from GF256 to F_q.

## A stronger rank bound across tiles

The total feedback is the sum of the eight tile feedback vectors. Its image
contains each tile's residual image, which gives rank at least
min(4,max_l d_l). There is also an inexpensive strengthening when several
tiles each leave one residual dimension.

For d=1, the residual polynomial is a scalar multiple of a monic polynomial
of degree seven. Put H=span_F2{1,2,4,64}. All seven roots lie in alpha+H,
so their sum lies in alpha+H. That sum is the polynomial's coefficient
of X^6. Its coefficient of X^7 is one.

Distinct tile offsets alpha in {0,8,...,56} lie in distinct cosets modulo H.
Consequently the normalized highest coefficients differ between tiles.
The same invertible triangular basis conversion acts on their high
coefficients, so their feedback lines are distinct. Two such lines have
joint rank two, also over F_q.

If b counts tiles with d_l>0, a valid combined rank lower bound is

    rho(d_1,...,d_8)
       = max(min(4,max_l d_l), min(2,b)).

This adds a genuine second feedback equation when all surviving tile images
are one-dimensional. No claim that the sum has rank min(4,sum_l d_l)
is made. Higher intersections require information about the actual subspaces.

## A small positive dynamic program

For potential counts a_l and selected-zero counts v_l in each tile,
put d_l=max(a_l-v_l,0). The probability of all output-zero equations is
q^(-sum_l min(a_l,v_l)). Conditioned on those equations, the surviving
input is uniform on the product of the residual tile spaces. Feedback zero
has probability at most q^(-rho). This conditional subspace calculation
justifies multiplying the two factors.

Averaging the input masks and summing output masks gives

    K_j(z) <= 1/binom(64,j)
       * sum_(a_1+...+a_8=j) sum_(v_1,...,v_8)
         product_l [binom(8,a_l) binom(8,v_l)
                    z^(8-v_l)(1-z)^v_l q^(-min(a_l,v_l))]
         * q^(-rho(d_1,...,d_8)).

All sums use 0<=a_l,v_l<=8. Omitting q^(-rho) gives the exact T_j(z).
The matrix-minor property of each FFT8 feedthrough makes the output ranks
independent of the particular masks at these fixed sizes.

The helper evaluates these sums with states (j,rho). There are only five
feedback states, numbered zero through four. For a new capped residual
dimension d, the update is

    rho'=2             if rho=d=1,
    rho'=max(rho,d)    otherwise.

State two includes either a rank-two tile or two distinct rank-one tile
images. Both interpretations preserve the lower bound under later updates.
Each tile coefficient is scaled by q^a; the final j-coefficient is divided
by q^j binom(64,j). Scaling prevents underflow in the floating screen.

There is an additional independent bound K_j(z)<=q^(-min(j,4)), because
z^wt(FX)<=1 and C_S has rank min(j,4). The helper takes the smaller of
these two complete moment bounds. At z=1 it uses this exact syndrome law.
It never combines favorable pieces of different moment decompositions.

## Regional composition and a whole interval

Fix an active-row set of size r. Conditional on its packet counts, average
the uniform amplitudes and within-packet subsets before applying the
fractional Markov bound. The [earlier route argument](HD_RANDOMIZED_OUTER_PROOF.md)
then applies with the computed K_j(z) and free-final T_j(z).

For rational marker p and exponent 0<alpha<=1, define

    S = sum_(j=0)^64 binom(64,j)p^j(1-p)^(64-j) K_j(z)^alpha,
    T = sum_(j=0)^64 binom(64,j)p^j(1-p)^(64-j) T_j(z)^alpha,
    c_r = binom(R,r)p^r(1-p)^(R-r),
    beta = q^256/(q-1)^128.

The upper bound for the complete r-row quiet class is

    B_r = binom(R,r) beta^(alpha r) z^(-alpha D)/(q-1)^alpha
             * S^(32767) T / c_r^256.

The outer binomial factor unions the active-row subsets. Its position
outside alpha is essential. The power 256 on c_r is the exact regional
normalization; p does not replace routing by independent Bernoulli choices.
The factor T accounts for the unconstrained final packet.

Use fixed witnesses

    z=1/16384,       p=127/128,       alpha=3/8.

At these fixed values, log2 B_r is an affine function of r minus
255 log2 binom(8192,r). Binomial coefficients are log-concave, so this
sequence is convex. Its maximum on a consecutive interval occurs at one
of the endpoints.

Directed 160-bit Arb evaluations give

    log2 B_5940 <= -3843.588047310521286548128518141267539159256,
    log2 B_8192 <= -532766.886869526820505645727438138583970150653.

There are 2253 occupancies in this interval. Convexity and a union bound yield

    log2 sum_(r=5940)^8192 B_r
       <= -3832.450415712285858591541173418694288524500.

This is a whole-class bound obtained from two directed endpoint calculations
and a proved convexity property. It is not interpolation between sampled
successful occupancies.

## Bounded screen and remaining work

A short witness grid was evaluated at nine selected occupancies. Directed
replay confirmed useful bounds at 6144 and 6565. The selected values
954, 1024, 1536, 2048, 4096, 5120, and 5632 were uninformative. These
failed witnesses do not establish a construction defect or an impossibility.
The separate endpoint calculation above selected the conservative threshold
5940 without a full occupancy scan.

The distinct-line refinement improves the new endpoint at r=6144 by about
0.54 bits. Its small numerical effect suggests that tracking more feedback
subspace information alone may not close the remaining middle range.
The next useful step is to compare that loss with the separate loss from
the regional coefficient marker, while preserving the joint local moment.

The new implementation is
[theory/hd_randomized_outer_joint_quiet.py](theory/hd_randomized_outer_joint_quiet.py).
[Seven tests](theory/hd_randomized_outer_joint_quiet_tests.py) pass. They check
1296 actual beta64 tile-mask cases, a bounded family of 1792 residual lines,
exact local MDS enumerators, the free-input moment, the zero-input mass,
the unweighted baseline, and directed point and interval endpoints.
The root-sum proof covers all residual-line cases beyond the tested family.

From the theory directory:

    python -B -m unittest hd_randomized_outer_joint_quiet_tests -v
    python -B hd_randomized_outer_joint_quiet.py

An independent mathematical audit checked the distinct-line argument and
the convex whole-interval step. No benchmark, native change, remote run,
or commit was made for this work.
