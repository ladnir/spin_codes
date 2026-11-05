# Rank-one refresh: sparse success and medium-weight proof obstructions

Replacing the current dense refresh could delete twenty fixed GFNI calls
per 64-record packet. The first proof gate does not justify that change.
Its sparse bounds pass, but its three-component majorant cannot close a
whole interval of medium weights, for any weight mark or fractional exponent.
This is a limitation of the majorant, not a low-distance codeword.
Retaining all block masks improves the tested medium bound, but does not
close it. A single controlled projection test also fails. No native change
or whole-code certificate follows from these diagnostics. Grouping masks,
a bounded joint witness search, and one two-refresh control also remain
positive at the tested medium profile.

## Construction and refresh law

Keep the Toeplitz outer, independently swapped alternating-pair route,
packet width64, state dimension16, and four RS[16,4] emission blocks.
Work over GF256. Write A for the emission map and C=A^T; CA=0.
The proposed packet step is

    Y = X + A Q,       W = Q + C X,       Q_next = M W.

The initial state is zero and the terminal state is unrestricted. For each
packet, independently sample u nonzero in F^16 and v uniformly subject to
v^T u != 1, then set M=I+uv^T. Sampling occurs at setup, once per packet.
This order differs from the earlier blockwise recurrence M Q+C X.

For nonzero W, the exact refresh kernel is

    P = a I + (1/q-a) L + (1-1/q) J,
    a = (q^15-1)/(q^16-1),       q=256.

Here L samples a nonzero scalar multiple of W, and J samples all nonzero
states. This is the law in [the rank-one analysis](HD_BLOCKWISE_INNER.md).
The combined identity/line mass is eta=1/q, not an independent uniform
refresh. Both these kernels preserve uniform-live density. Put gamma=1-eta.

## The tested three-component cap

Bound an incoming weighted state measure by Z at zero, arbitrary-live
mass B, and density U/d at each nonzero state, where d=q^16-1.
These are overlapping majorant components, not disjoint events.

Fix a potential input support with block occupancies l_i. Its values X
are uniform over that coordinate subspace, including zero values, as in
the established outer comparison. Define

    l = sum_i l_i,       rho = sum_i min(4,l_i),
    delta = max(0,13-max_i l_i),       t=z+(1-z)/q.

Use the existing ideal packet quantities K, T-K, Jbar/d, V from
[the common-refresh analysis](HD_COMMON_REFRESH.md), without its diagonal
density penalty or its inverse-refresh empty-packet bound. In particular,
T=t^l, K is the quiet marked mass, and V is the uniform-live marked mass.
The following additional caps hold for each fixed support:

    A_live = t^l z^delta,
    R_live = q^-rho min(t^(l-rho), z^max(1,delta))   if l>0,
    R_live = 0                                    if l=0,
    L_live = min(t^l,z^delta).

For arbitrary nonzero Q, active coordinates average to t^l; at least delta
inactive coordinates must remain nonzero. This proves A_live. Conditioning
CX=Q leaves an affine space of dimension l-rho. Its translated output has
that dimension because I+AC is invertible, proving the information-set
part of R_live. The return output is nonzero, giving its other cap.

For fixed nonzero W, the map X to (I+AC)X+AW is an affine injective image.
Excluding Q=0 can only decrease its weighted sum. Its information-set and
minimum-weight bounds give L_live, a pointwise pre-refresh density cap.

A valid matrix updating the row vector (Z,B,U) is

    [ K       eta(T-K)       gamma(T-K)       ]
    [ R_live  eta A_live     gamma A_live     ]
    [ Jbar/d  0              eta L_live+gamma V ].

The last row retains the identity/line branch as a density bound. Another
valid choice replaces only that row by (Jbar/d, eta V, gamma V), retaining
its nonuniform part as arbitrary mass. These choices are alternatives;
their entries cannot be minimized separately without further justification.
Each selected matrix is raised entrywise to alpha before averaging actual
pair placements. This yields a valid but potentially loose fractional bound.

## Bounded FLOAT gate

The new isolated helper is `theory/hd_rankone_common_screen.py`.
It changes no existing proof defaults or native code. At alpha=1, positive
truncated extraction includes every source class at weights one and two:

| Active groups | Weight mark | All-support log2 bound, FLOAT |
|---:|---:|---:|
| 1 | 0.9991588324466113 | -121.9270218 |
| 2 | 0.9983218568920199 | -236.6596796 |

For the single matching profile (r,a,b)=(1024,400,400), hold
z=0.4828132651072611 and alpha=0.10118724647150915 fixed. Optimizing only
two positive coefficient markers, with at most100 objective evaluations,
gives +125681.95 for the density row and +112117.03 for the mass row.
The same positive-marker engine with the old dense refresh gives +3978.93.
Thus this coarse extraction already loses some old-proof margin, but it
does not account for the much larger new loss. No z/alpha search was run
for this three-component cap.

These results use the sharper existing Toeplitz outer measure
beta=q^128(1-1/q)^-32. They are not directed receipts or whole-code bounds.

## Why further mark optimization cannot repair this cap

Let P=32768 be the packet count. The arbitrary-live self-loop satisfies

    eta A_live >= q^-1 z^(l+13).

If the first packet has an active potential coordinate, its birth entry
is at least q^-1 gamma z^l. For any fixed source matching profile of weight
r, the probability that the first routed potential coordinate is active is
r/8192. Across all packets, the potential degrees sum to256r.

Consequently, even charging only one source support, the resulting
fractional majorant has log2 value at least

    log2(r/8192)
    + alpha [r log2(beta) - log2(255) - 8P + log2(gamma)
             + (256r+13(P-1)-734003) log2(z)].

This is a lower bound on the evaluated majorant, not on setup failure.
For every integer 257<=r<=1203, the coefficient of log2(z) is nonpositive.
The remaining bracket exceeds 1024r-8P-9, which is positive throughout
that interval. Therefore this majorant exceeds r/8192>2^-5 for every
0<z<=1 and 0<alpha<=1. Exact coefficient extraction cannot repair it.
Neither choice of the flat row changes this obstruction.

The majorant loses the identity/line branch's persistent block structure.
It permits a fresh arbitrary state of emission weight13 at each packet.
The actual recurrence must also account for which blocks CX activates.

## The implemented seventeen-component cap

Retain zero, one arbitrary-live mass for each of the15 nonempty masks of
the four state blocks, and one uniform-live density component:17 components.
The identity/line refresh branch preserves the block mask of W=Q+CX.
Its mask transition must therefore use the actual feedback law, instead of
resetting to an arbitrary live state.

For one block, a zero incoming state has marked transitions K_l and
T_l-K_l. A live incoming state admits the fixed-state return cap R_l and
the total-emission cap A_l. Independent X blocks permit a tensor product
of these positive two-state caps. Multiply nonzero successor masks by eta;
send the gamma part to the flat component. A zero successor has probability
one under refresh and must not be multiplied by eta.

The isolated implementation is `theory/hd_rankone_mask_screen.py`. Its
sixteen mask components include mask0 for the zero state; component16 is
the uniform-live density. The tensor argument conditions on the entire
incoming state. It requires independence of the input blocks, not independence
of the incoming state blocks.

For an incoming mask I, let F_(I,S) be the tensor cap for successor mask S.
Let G_I be the product of T_l on zero incoming blocks and A_l on live blocks.
Put N_I=min(G_I,sum_(S!=0) F_(I,S)); for I=0 use the exact value T-K.
The outgoing caps are

    zero: F_(I,0),       mask S!=0: eta F_(I,S),       flat: gamma N_I.

The flat row uses a separate positive calculation. Put b=q^-4. In block i,
write K_i for its quiet marked mass, J_i for its nonquiet return cap, and
T_i=t^l_i. Write E_i for the excess over b in the uniform-state marked
moment on the 16-l_i inactive emission coordinates. These are the existing
local quantities in `hd_rankone_common_screen.local_values`.

For a uniformly sampled block state, including zero, the two successor
caps are b(K_i+J_i) for zero and b(T_i-K_i)+T_i E_i for nonzero.
The second cap subtracts the exact quiet contribution from the total cap
T_i(b+E_i). For a target mask S, define

    Z_i = K_i       and D_i = b J_i    when i is absent from S,
    Z_i = T_i-K_i   and D_i = T_i E_i  when i belongs to S.

The product of bZ_i+D_i bounds the full-uniform marked mass on S.
The excluded incoming all-zero state contributes exactly b^4 product_i Z_i.
Consequently a valid uniform-live cap on S is

    H_S = min(V, [product_i(b Z_i+D_i)-product_i(b Z_i)]/(1-b^4)).

The numerator is evaluated as a sum of nonnegative terms, not by subtracting
nearby floating-point values. The flat row sends H_0 to zero, eta H_S to
each nonempty mask, and gamma min(V,sum_(S!=0) H_S) to the flat component.
The zero entry agrees with the previous Jbar/d return bound.

The regional polynomial retains all four occupancy labels. Sorting
occupancies would change the state-mask operator. Each packet averages
entrywise alpha powers over its actual pair placements and fair swaps.
The implementation does not sample new source matching geometry between
regions.

## Mask gates and one controlled projection

At alpha=1, the seventeen-component operator gives the following complete
weight-one and weight-two FLOAT bounds. These use positive truncated
extraction, not an alias approximation.

| Active groups | Weight mark | All-support log2 bound, FLOAT |
|---:|---:|---:|
| 1 | 0.9991588324466113 | -121.9559023085 |
| 2 | 0.9983218568920199 | -240.9645891946 |

For the same medium profile and fixed z/alpha as above, the seventeen-component
bound is +99844.04545436724. The marker fit used64 objective evaluations and
returned log2 markers (-3.710160740811969,-2.660319071541368). As in the
three-component test, the support count is a single-matching upper bound.

The overlapping flat and arbitrary components can lose precision when their
matrix entries are separately raised to alpha. A controlled alternative
removes the flat component before that operation. For each nonempty mask S,
the uniform-live distribution has mask probability

    pi_S = (q^4-1)^|S|/(q^16-1),       pi_0=0.

Replace a row's nonempty-mask cap by its existing cap plus pi_S times its
flat cap. Leave the zero entry unchanged. This is a valid sixteen-component
measure majorant: it retains each mask's mass and discards its density bound.
The two contributions are added before applying alpha. In particular,
uniform-live mass is never projected onto mask0.

One evaluation at the unchanged z, alpha, and markers gives
+103577.27827627852. No projection-specific fitting was performed. The
projection therefore worsens this witness by about3733.23 bits. It removes
one fractionalization loss but discards useful density information; this
experiment does not determine how much of the remaining loss comes from
either effect.

These are FLOAT proof screens, not directed receipts. The sparse tests do
not cover other weights, and the medium failures are not construction
counterexamples. The seventeen-component and projected operators do not
inherit the universal three-component obstruction proved above.

## Grouping masks before fractionalization

A smaller cone combines destination masks before applying alpha.
This differs from grouping entries of the already fractionalized matrix.
It retains persistent activity while avoiding a separate alpha power for
every destination mask.

The actual pair law does not have full S4 block symmetry. Its pairs connect
blocks0 and1, or blocks2 and3, at matching coordinates. For example, a
packet containing one full pair can occupy blocks0 and1 but not blocks0
and2. The valid block symmetries swap within either paired half and exchange
the two halves. This order-eight group has six mask orbits:

    empty; one live block; two live blocks in one paired half;
    two live blocks in different halves; three live blocks; all four.

Keeping arbitrary mass in each orbit plus one flat component gives seven
components. A six-component version using only live-block cardinality is
also possible, but must explicitly maximize each row over both size-two
source orbits. Their equality cannot be inferred from the route law.

Fix a packet type (s,f), with
s singleton pairs and f full pairs. Conditional on all packet types in a
region, local pair placements and fair swaps are independent between
packets. Average the labelled caps over those local choices, then sum the
arbitrary-mass entries within each destination orbit. Keep zero and flat
components separate. For source orbit I and destination orbit J, the
implementation uses max_(i in I) sum_(j in J) E[M_ij | type]. This bounds
every incoming mixture on I without relying on numerically identical rows.
It does not resample the incoming state. Conditional on the entire type
sequence, future local placements are independent of incoming history.

Let H_(s,f) be this conditional mean cap. If Z is the complete marked
moment for a fixed sequence of local placements, concavity gives

    E[Z^alpha | all packet types]
      <= (E[Z | all packet types])^alpha
      <= (e_0 product_j H_(s_j,f_j) 1)^alpha
      <= e_0 product_j (H_(s_j,f_j) entrywise^alpha) 1.

This order first averages local placements, then groups destination masses,
then fractionalizes. Averaging introduces a Jensen loss, so the smaller
state space alone does not guarantee an improved bound.

The polynomial coefficient must retain the probability of each packet
type. With N_(s,f)=32!/((32-s-f)!s!f!), that coefficient is
N_(s,f) times H_(s,f) entrywise^alpha. Raising N_(s,f) H_(s,f) itself to
alpha would give the wrong placement normalization. The regional coefficient
extraction and global matching counts remain unchanged.

The isolated helper `theory/hd_rankone_orbit_screen.py` implements this
seven-component interface. One evaluation at the unchanged medium witness
and markers gives +101453.78397569829, versus +99844.04545436724 for the
seventeen-component cap. The grouped bound is worse by1609.73852133 bits.
There was no grouped-operator fit. This control does not show that grouping
can never help; it shows that this particular grouping loses at this witness.

## Bounded joint witness search

The previous z and alpha came from the dense-refresh proof. To test whether
that choice caused the medium failure, the seventeen-component cap received
one bounded joint search over z, alpha, and both coefficient markers.
The source profile remained (1024,400,400). The seven-component operator
and the obstructed three-component cap were not searched.

The new helper `theory/hd_rankone_mask_proposal.py` evaluates positive
markers without constructing every polynomial coefficient. For positive
markers x,y and formal variables u,v, the half-occupancy weight at (i,j) is

    [u^i v^j] (1+x(u+v)/2+yuv)^16.

The full packet weight is the product of its two ordered half weights.
Contracting these weights against the labelled packet caps, already raised
entrywise to alpha, gives exactly the previous marked polynomial evaluation.
This changes the screen computation, not the encoder or proof inequality.
At the frozen witness it reproduces +99844.0454543673, within6e-11 bits
of the previous evaluation. Small-degree tests independently compare every
matrix entry against explicit coefficient evaluation.

The budget comprised 47 search evaluations plus that initial agreement
control. Eight deterministic seeds preceded one bounded four-variable
Nelder-Mead search. Its best witness was

    z = 0.4353679232742893,
    alpha = 0.17088472633321333,
    log2 markers = (-3.3650867183232807,-1.996947766004256).

The FLOAT bound is +92363.8038867295. This improves the previous bound by
about7480.24 bits but does not approach closure. The bounded search does
not establish a global optimum or a construction counterexample.

A Perron diagnostic locates the long-run contribution within the marked
positive matrix. For nonnegative left and right Perron vectors ell,v,
normalize ell^T v=1 and report component weights ell_i v_i. At the best
witness, zero and masks of sizes1,2,3,4 contribute respectively
23.722%,20.859%,27.797%,19.032%,7.446%; the flat component contributes1.143%.
Thus the majorant's contribution is spread across mask classes rather than
concentrated in the flat component or one small mask. These are diagnostic
weights of the majorant, not probabilities of the actual state process.

## Two independent refreshes at one fixed witness

A separate construction applies two independent rank-one updates after
the same feedback addition. Its next state is M_2 M_1(Q+CX). Both updates
use the sampler defined above and are independent of all other packet keys.
Since L^2=L, LJ=JL=J, and J^2=J, the exact nonzero-state kernel is

    P^2 = a^2 I + (q^-2-a^2) L + (1-q^-2) J,

where a is the coefficient in the single-refresh law. The combined
identity/line probability becomes q^-2, not q^-1. Zero remains fixed.
In the seventeen-component cap, divide all nonempty-mask entries by q,
multiply flat entries by (1-q^-2)/(1-q^-1), and leave zero entries unchanged.
An exact GF4/state2 enumeration verifies this composition, including its
scalar-line term.

One evaluation at the best one-refresh witness gives +65853.25259219245.
No two-refresh fitting or additional profile screen followed. Two updates
would use sixteen variable GFNI products instead of the current refresh's
eight variable and twenty fixed products: twelve calls removed per packet.
That is an operation count, not a measured speedup. The failed bound does
not justify implementing this change in the encoder.

## Next proposal: alternating rank-one and dense refreshes

A different proposal would retain the current dense refresh after every
second packet. That refresh is D_2 H_dense D_1: H_dense is the fixed
interpolation map, and both diagonals have independent nonzero entries.
The intervening packet would use one rank-one refresh. The initial state
remains zero, and each packet applies its selected refresh to Q+CX.
All refresh keys would be sampled independently at setup as before.

At dense-refresh boundaries, retain only zero mass and uniform-live density.
The rank-one packet maps these two components into the seventeen-component
cone, using its zero and flat source rows. The following dense packet maps
those seventeen components back to two. Its zero entry uses the same return
cap. Its flat entry is at most b_D times the nonzero-emission cap, where
b_D=(q^16-1)/(q-1)^16 is the established diagonal-Cauchy density factor.
For a flat source, the total-emission cap V is also safe.

Compose this 2-by-17 and 17-by-2 pair before fractionalization. Dense refresh
then prevents a persistent arbitrary-mask component from crossing pair
boundaries. The proposal still needs a joint placement calculation with
the correct conditioning; no paired gate has been implemented or evaluated.
It would remove twenty fixed GFNI calls per two packets, or ten per packet
on average. No native timing or distance claim is attached to that count.

This paired gate, rather than another unrestricted one-refresh marker
search, is the recommended next structural proof test. All results above
remain FLOAT screens; none is a directed or whole-code certificate.

## Tests and reproduction

Six tests check the exact GF4/state2 refresh law, scalar-line probabilities,
the three-component cone for every live toy state, both flat-row choices,
ordered polynomial arithmetic, route normalization, and the path obstruction:

    python -B -m unittest hd_rankone_common_screen_tests

Seven further tests check the labelled mask cone against exact GF4 toy
probabilities, the flat return identity, block-permutation conjugacy,
exhaustive placements through weight two, ordered matrix products, resource
guards, and the projection's zero-mass and pre-alpha rules:

    python -B -m unittest hd_rankone_mask_screen_tests

Reproduce the single fixed projection from the theory directory:

    python -B -c "from hd_rankone_mask_screen import fixed_symmetric_gate; print(fixed_symmetric_gate(1024,400,.4828132651072611,.10118724647150915,[-3.710160740811969,-2.660319071541368],projection=True))"

Six orbit tests check the exact matching symmetry, row domination,
placement normalization, alpha-one equivalence, Jensen chronology, and
resource guards. Five proposal tests check direct contraction, the exact
two-refresh law, branch scaling, Perron diagnostics, and the hard budget.
The previously completed combined run passed all 24 tests:

    python -B -m unittest hd_rankone_common_screen_tests hd_rankone_mask_screen_tests hd_rankone_orbit_screen_tests hd_rankone_mask_proposal_tests

The bounded search and its fixed controls are reproducible from the theory
directory. These commands perform new FLOAT calculations when invoked:

    python -B hd_rankone_orbit_screen.py
    python -B hd_rankone_mask_proposal.py --search --max-evaluations 47
    python -B hd_rankone_mask_proposal.py --best
    python -B hd_rankone_mask_proposal.py --best --two-refresh
