# A finite rank continuation for the coupled scan

This note gives a composable upper bound for the construction in
[CONSTRUCTION.md](CONSTRUCTION.md). It is not a distance certificate. The
bound explicitly drops some deterministic rank gains; its usefulness must
be tested before committing to a large calculation.

Write q=256 and fix 1 <= tau <= q. The output marker z is positive; use
0 < z <= 1 when bounding output weight at most e. The state wiring R may be
the original cycle or a proposed pair-antipodal cycle. In either case it is
one fixed w-cycle, with w=4 or w=8. This note does not change the specified
construction or claim a certificate for either wiring.

## 1. A positive matrix uniform over the incoming quotient

Keep the common message-constraint space W and all previous nonzero
requirements from [PROOF.md](PROOF.md). Fix the exact input support I and
incoming state support H for one step. Let J and H' be proposed output and
outgoing state supports, and put K=R(H). All sets are subsets of the w
coordinate positions.

The two additions impose the necessary support conditions

    I symmetric_difference H subset J subset I union H,
    J symmetric_difference K subset H' subset J union K.

These conditions use nonzero diagonal coefficients. They do not assert
that every admitted support pattern has a witness.

Also exclude I=empty, H!=empty, H'=empty: the determinant lemma rules out
that whole-state extinction. In particular, I=empty forces J=H; and
H=empty forces J=H'=I.

For each admitted triple define

    p(H,J,H') = |H minus (J intersection H')|.

The fixed-mask lemma charges one fresh coefficient at each of these p
coordinates. There may also be deterministic rank gains d. Dropping d
gives the following 2^w by 2^w matrix:

    M_I(z)[H,H'] = sum_J z^|J| tau^(-p(H,J,H')),

where the sum includes only the admitted triples above. Every coefficient
is nonnegative.

For any incoming quotient, pool of input forms, and earlier nonzero
requirements, M_I bounds the expected sum of relative rank potentials over
the listed one-step branches. This is the fixed-mask lemma summed over J.
It is not a Markov model for a fixed message. Earlier equations remain in
W; the uniform bound allows the same matrix to apply to the new quotient
after each step. No independent-equation assumption is made.

If I is also to be enumerated, retain its physical positions rather than
only |I|. Pair-antipodal wiring acts on this geometry. Collapsing input
masks of equal weight can discard the reason for changing the wiring.

## 2. Average one actual region of the paired route

A region has M=4096 physical pairs. Each step consumes b=w/2 pairs, so
there are M/b steps per region. For a fixed exact source support S, let
n_0, n_1, n_2 count source pairs containing zero, one, and two positions of
S in that region. The matching used to compute these counts is the
specified even or odd source matching; n_0+n_1+n_2=M.

Introduce commuting markers x_0,x_1,x_2. For a sequence omega of b physical
pair masks, let I(omega) be the corresponding input-coordinate mask. Give
each physical pair the weight

    00: x_0,    01: x_1/2,    10: x_1/2,    11: x_2.

Define the polynomial matrix

    T(x_0,x_1,x_2;z)
      = sum_omega (product of its b pair weights) M_{I(omega)}(z).

Then a valid regional continuation matrix is

    B(n_0,n_1,n_2;z)
      = [x_0^n_0 x_1^n_1 x_2^n_2] T(x_0,x_1,x_2;z)^(M/b)
        / multinomial(M; n_0,n_1,n_2).

This normalization is specific to the actual paired route. A uniform
permutation of distinct source pairs gives a uniform word of the three
pair occupancy types. Each singleton is independently oriented by its
fair swap. Swapping a double-active pair changes the two nonzero forms,
but not its support; the uniform incoming-quotient bound covers either
order and their average. The same applies to the remaining label order
within each occupancy type.

The matrix power preserves the chronological order of steps. Pair types
are sampled without replacement through coefficient extraction and the
multinomial divisor; they are not replaced by independent per-step draws.

## 3. Compose regions and cover messages

Fix an outer realization and an exact, nonempty outer support S before
routing. Let V_S be the outer message space shortened to S, and let
D_S=dim_F V_S. Discard D_S=0. Also discard S if a coordinate required
nonzero is identically zero on V_S. Begin with W=0 inside V_S and require
every source coordinate in S to be nonzero.

For region r, compute its pair counts n(r,S). Let e_empty select the empty
state and let 1 be the vector with every entry equal to one. The latter
includes every terminal state; no final reset is required. A valid bound
for the existence of a nonzero message with this exact support and final
weight at most e is

    tau^(D_S-1) z^(-e)
      e_empty^T (product over r=0,...,255 of B(n(r,S);z)) 1.

The bound averages over the independent regional routes, swaps, and fresh
inner coefficients. It is conditional on the outer realization. Each
regional matrix is uniform over the quotient left by the preceding
regions, so multiplication does not reset W or assume independent rank
increments. The factor tau^(D_S-1) is the usual rank-event domination for
a nonzero witness. Finally, z^weight >= z^e on the forbidden event.

Sum this expression over exact supports S and average over the actual
outer setup. Supports that cannot occur may be included as an explicit
overcount, provided the corresponding source-dimension bound is valid.
An optimized coefficient-extraction bound must still cover every support
profile it claims to cover. No profile has been numerically evaluated here.

## 4. A valid bridge to the GF256 outer enumerator

The outer matrices have entries in E=GF256. For any fixed outer setup and
coordinate set S, the shortened dimension is the same over E and F:
matrix rank does not change under field extension. Consequently

    number of E-codewords supported within S = q^(D_S).

The word "within" is essential: it includes codewords with smaller
support. It does not claim uniform F-valued amplitudes or identify exact
support multiplicities over the two fields.

For d>=1 and 1<=tau<=q,

    tau^(d-1) <= (q^d-1)/(q-1).

For s active outer groups with shortened dimensions d_i>=1, the direct-sum
outer therefore gives

    tau^(sum_i d_i-1)
      <= tau^(s-1)/(q-1)^s product_i (q^(d_i)-1).

Each factor q^(d_i)-1 counts nonzero E-codewords supported within that
group's selected support. Thus actual GF256 support statistics can bound
the initial rank potential. The factor tau^(s-1) must not be omitted.
Independence of the outer setups permits averaging the product group by
group when the group supports are fixed.

For a fixed length-n E-code C, an illustrative generating identity is

    sum_{S subset [n]} u^|S| (q^(dim C_S)-1)
      = sum_{c in C, c!=0} u^wt(c) (1+u)^(n-wt(c)).

It follows by first fixing a codeword and summing its containing supports.
This identity is a support-containment transform, not an exact-support
identity. The global route depends on the placement of supports across
regions and neighboring groups. An ordinary weight enumerator alone
cannot replace those profiles in the regional product. Use a refined
support measure, or prove a domination that preserves the required route
weights. The supersets in this bridge are also a possible source of slack.

### Why averaging the real outer distribution is necessary

A conditional guarantee uniform over every allowed outer realization is
false. Consider the permitted realization in which every Toeplitz mixer
is the identity. In each group, use just one RS row and the same nonzero
degree-31 evaluation polynomial with 31 chosen evaluation roots, scaled
by one free F amplitude per group. Nonzero column multipliers in a GRS
variant do not change its support.

These 8192 amplitudes form a message subspace of dimension 8192. The outer
is supported in only 33 of the 256 scalar regions. For each of these live
regions impose Q=0 at its end. This uses at most 33w equations, fewer than
8192 for the candidate widths. A nonzero message survives, and every
other region emits zero: its input is zero and its incoming state is zero.
Thus this outer realization has

    d_min <= 33N/256,

regardless of inner coefficients. The outer realization may be extremely
unlikely; this does not refute an averaged setup-failure guarantee. It
does rule out proving the target by discarding outer randomness entirely.

## 5. Low-charge paths to retain in the audit

The natural cycle lets an input on one physical pair extinguish a singleton
state. Pair-antipodal wiring removes that particular path: a reset driven
by one physical pair requires at least w/2 active state coordinates, hence
at least w/2 fresh-coefficient charges. With two occupied physical pairs,
singleton resets can still occur. Their route cost must be counted rather
than ruled out by input weight alone.

There is a different path that exposes the cost of dropping deterministic
rank. Let the current state be supported only at j, and let k=R(j).
Choose the input only at k, with X[k]=Q[j]. The next state remains supported
only at j; the output is supported at j and k. No fresh coefficient is
charged, because the current live coordinate is nonzero in both output
and next state. The relation X[k]=Q[j] may add deterministic rank, but the
matrix M_I drops that rank.

Repeatedly admitting this path yields output density 2/w and input density
1/w. At w=8 its output density is 25%, below the target. This is not a
constructed codeword: future input labels cannot be chosen adaptively to
match the state. The actual outer, route normalization, and retained
deterministic equations must charge those matches. It is a specific
potential weakness of the relaxed continuation, not evidence of an actual
25% distance ceiling. At w=4 this particular path emits 50%.

The next section recovers a fresh-coefficient charge for this repeated
localization path. It changes the proof exposure, not the encoder. Other
low-charge paths may still require additional deterministic rank. Do not
respond to such a failure by expanding numerical precision or profile
searches without identifying the lost constraint.

## 6. One-step lookahead recovers localization rank

Fix the outer realization and the full route before exposing any inner
coefficients. Every input coordinate is then a known linear form on the
shortened message space. In particular, a future input form does not depend
on the current inner coefficients.

Fix a complete proposed sequence of support masks before taking the
expectation over inner setup. These masks specify the branch being counted;
they are not events on which to condition and renormalize the setup. The
proof may impose a necessary future zero equation early when that equation
uses only the current coefficient and already known forms.

Write rho(j) for the destination of state coordinate j under R, so
(R Q)[rho(j)]=Q[j]. There are T=N/w steps, indexed 0 through T-1. For a step
t<T-1 define the ordinary charged set and an additional lookahead set:

    P_t = H_t minus (J_t intersection H_(t+1)),

    L_t = {j in H_t intersection J_t intersection H_(t+1) :
             rho(j) is outside H_(t+1) union H_(t+2)}.

The two sets are disjoint. The definition of L_t requires the current live
coordinate to survive in both output and next state, so its coefficient
has not already received an ordinary charge.

For j in L_t, the proposed branch requires

    Q_(t+1)[rho(j)] = Q_(t+2)[rho(j)] = 0.

Applying the recurrence at step t+1 eliminates its diagonal coefficient:

    X_(t+1)[rho(j)] + Q_(t+1)[j] = 0.

Substituting the current update gives the necessary equation

    X_(t+1)[rho(j)] + (R Q_t)[j] + X_t[j]
      + a[t,j] Q_t[j] = 0.

Everything except a[t,j] is known before the current coefficients are
exposed. Also Q_t[j] is required nonzero. This is precisely the affine
fresh-coefficient equation of the fixed-mask rank lemma.

First impose the coefficient-independent current equations, with rank
increase d. Then expose the coefficients indexed by P_t union L_t, in any
fixed order, and impose their respective equations in the same W. Each
coefficient is charged once. At its exposure, either its equation adds one
rank, or survival requires its unique exceptional value with probability
at most 1/q. Earlier equations may change which case applies. If they
annihilate a previously required state form, discard the branch.

Consequently

    E[Phi_after] <= tau^(-d-|P_t|-|L_t|) Phi_before.

Append the current new nonzero requirements after their coefficients have
been exposed, as in the original lemma. An anticipated equation is already
in W when the next chronological step reaches it. It must not be charged
again as a new deterministic rank gain.

This argument never conditions on the proposed future masks. It counts
their necessary equations while the relevant coefficient is still fresh.
There is no lookahead charge at t=T-1, because the construction has no next
input step. The terminal state is still free.

### Consecutive singleton localization

Suppose L consecutive steps retain the singleton state support {j}, have
input support {rho(j)}, and output support {j,rho(j)}. The first match is

    X_t[rho(j)] = Q_t[j].

It is a coefficient-independent equation and may add zero or one rank.
Write its actual rank increment as d_0. The update then gives

    Q_(t+1)[j] = a[t,j] Q_t[j].

The next localization match requires

    X_(t+1)[rho(j)] = a[t,j] Q_t[j].

Charge this equation before exposing a[t,j]. Repeating the argument gives
L-1 different fresh-coefficient charges. Conditional on the preceding
exposure, the potential of the whole specified run is bounded by

    E[Phi_run_end] <= tau^(-d_0-(L-1)) Phi_run_start.

The last step's coefficient needs no charge when localization is not
required afterward. Thus the rank-free run admitted by the first-order
relaxation is not actually rank-free. No extra random diagonal or encoder
operation is needed to recover this particular cost.

### A finite augmented transfer

The extra charge uses three successive state masks. An explicit finite
transfer can therefore use augmented states (H,H_1). For fixed input mask
I, define a matrix whose only possible transitions are

    (H,H_1) -> (H_1,H_2).

For such a transition, sum over the admitted current output masks J and
give each the weight

    z^|J| tau^(-p(H,J,H_1)-ell(H,J,H_1,H_2)),

where ell is the cardinality of L_t above. This is a 2^(2w)-state
description, not a recommended dense storage layout. The matrix does not
need the next input mask in its weight: the rank bound is uniform over
the fixed future input forms, and the next chronological transition
checks that input mask's support compatibility.

Start with weight one on each pair (empty,H_1), and zero on other pairs.
For the last step use the terminal vector

    f_I(H,H_1) = M_I(z)[H,H_1],

which has only the ordinary charges. Multiplying the augmented transfers
for steps 0 through T-2 and then f for step T-1 sums every terminal state.
The first step's support compatibility enforces its actual H_1.

The paired-region coefficient extraction remains valid after this
augmentation. One may first bound inner expectations for each fixed full
route and then average those nonnegative mask-matrix products over the
route. All nonfinal regions use powers of the augmented pair-polynomial
matrix. The last region ends with its corresponding terminal pair-vector
instead of one augmented step. This preserves correlations across region
boundaries and does not insert artificial resets.
