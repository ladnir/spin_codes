# A whole-state bridge for the full-GF32 refresh

The full-GF32 diagonal refresh has a useful point-probability bound, but
the retained outer is only randomized over GF256. This note connects those
facts without asserting a uniform-GF32 outer density. It produces a positive
two-state operator for all actual source sizes. The operator is proved as an
upper bound; no complete numerical certificate or runtime follows from it.

The essential change is to count an entire live state, rather than its
coordinate support. A fourth root then converts a GF32 linear-space count
to a GF256 count. This root is applied to each requested-zero event before
the ordinary output-weight expansion. Auxiliary GF256 messages need not
have the same group support as the actual GF32 witness. Section 4 accounts
for that difference explicitly.

## 1. Construction and probability space

Use the exact construction in [SMALL_STATE_PROOF_GATE.md](SMALL_STATE_PROOF_GATE.md).
Write

    E=GF256,  q0=256,  F=GF(2^32),  q=q0^4,
    N=2^21,   n=N/64=32768,        s=16.

The outer has 8192 groups, each with 128 source symbols and 256 output
symbols. Retain its four RS[64,32] constituents, independent invertible
E-Toeplitz column mixers, and the alternating paired regional route.
For a fixed set S of actual source groups, I_t(S) is the set of coordinates
in batch t that can receive outputs from S. In particular,
sum_t |I_t(S)|=256|S|. This is an allowed set, not an exact nonzero support.

The inner is

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H16 D1_t (Q_t + C X_t),
    Q_0 = 0.

Here A, C, and H16 are the retained E-valued maps, C=A^T, and CA=0.
Each diagonal entry is independently sampled from F minus zero. Setup is
sampled once and fixed for all messages. The last state is discarded.

Let Bad_S be the existence of an F-valued message nonzero in every group
of S, zero outside S, and with output weight at most e=734003. The event
allows the message to depend on all sampled setup. Initially fix the outer
and route; only inner diagonal randomness remains.

For every fixed nonzero v and nonzero w in F^s, the sparse-gate lemma gives

    Pr[D2 H16 D1 v = w] <= gamma q^(-s),
    gamma = (q/(q-1))^s.                                      (1)

The map is invertible, so zero remains zero and a nonzero vector cannot
become zero during a refresh. Equation (1) holds conditionally on the past.

## 2. Fixed zero constraints and a fixed whole-state path

Fix requested output-zero sets Z_t and a path sigma_t in {0,1}, with
sigma_0=0. A value one denotes a nonzero incoming state in an actual
trajectory. No state after batch n-1 is specified. Put a=sum_t sigma_t.

Before relaxing this path, discard it if I_t(S) is empty and
sigma_(t+1) differs from sigma_t, for any t<n-1. Such a transition is
impossible for the actual code: X_t=0 and the refresh is invertible.

Introduce source variables m supported on S. Introduce an s-vector Q_t
when sigma_t=1, and set Q_t=0 otherwise. Impose only the linear equations

    (X_t + A Q_t) restricted to Z_t = 0,
    Q_t + C X_t = 0 whenever t<n-1 and sigma_(t+1)=0.             (2)

All equations have coefficients in E. Their coefficients depend on the
fixed outer and route but not on any inner diagonal. Denote their kernel
over E by V. Let d=dim_E V and let d0 be the dimension of its subspace
with m=0. Extension to F has the same two dimensions.

The relaxation drops Q_t!=0 on marked live states and drops the condition
Q_t+CX_t!=0 before a marked live next state. It retains every reset
equation in (2). These dropped inequalities can introduce artificial
trajectories, so the resulting operator is not the actual transition law.

For any fixed actual tuple (m,Q), each marked live state costs at most
gamma q^(-s), by (1). The fresh refreshes that generate those states are
distinct. Homogeneity lets us count projective tuples: multiplying m and
all Q by a common nonzero scalar describes the same existence event.
Thus the fixed-(Z,sigma) event has probability at most

    min(1, gamma^a q^(-sa) (q^d-q^d0)/(q-1)).                  (3)

If d=d0 the event is impossible. Otherwise set P_Q(k)=1+Q+...+Q^(k-1).
Since P_q0(k)^4 >= P_(q0^4)(k),

    [(q^d-q^d0)/(q-1)]^(1/4)
      <= (q0^d-q0^d0)/(q0-1).                                (4)

Applying min(1,x)<=x^(1/4) to (3) gives

    Pr[fixed Z,sigma event]
      <= gamma^(a/4) q0^(-sa) (q0^d-q0^d0)/(q0-1).             (5)

The numerator on the right is exactly the number of auxiliary E-valued
pairs in V with m!=0. It is not the number of actual F messages. A uniform
global factor is

    Gamma = gamma^((n-1)/4) < 2.                              (6)

For example, (1-1/q)^m>=1-m/q with m=4(n-1) proves (6). No charge uses
the unused final refresh.

## 3. Keep the fourth root outside the output-weight expansion

For 0<z<1 and any fixed output y,

    z^wt(y) = sum_(Z subset of zeros(y))
                 z^(N-|Z|) (1-z)^|Z|.                        (7)

If Bad_S occurs, select any bad witness and its actual state path.
Its contribution to the following sum is at least z^e:

    z^e 1[Bad_S]
      <= sum_Z z^(N-|Z|)(1-z)^|Z|
             sum_sigma 1[the fixed Z,sigma event occurs].     (8)

Apply (5) separately to each event in (8), then sum its nonnegative
weight. The weights in (7) are not raised to the fourth power or its
inverse. Concavity is not used to move a root through this weighted sum.

Expanding the auxiliary E-pair count in (5), and then summing Z, recovers
the ordinary factor z^wt(X+AQ). Each live state variable contributes
the normalization q0^(-s). This gives the local scalar operator in
Section 5 after the outer average.

## 4. Preserve the actual source class before choosing its marker

A single common marker z can be used directly in (8). To use different
markers for different actual sizes r=|S|, one must not identify S with
the support of an auxiliary E message.

For a retained fixed-(Z,sigma) event, every group projection of V is
nonzero. Otherwise every F solution has zero message in that group and
cannot belong to Bad_S. Fix a group whose projection has rank r_i>=1.
Among E pairs in V with m!=0, its missing probability is

    (q0^(d-r_i)-q0^d0)/(q0^d-q0^d0) < q0^(-r_i) <= 1/q0.       (9)

The zero-source subspace is contained in every missing-group subspace;
this is why the same d0 is subtracted in numerator and denominator.

Choose an integer h with 0<=h<r and

    c_r = 1-r/[q0(h+1)] > 0.                                 (10)

Markov's inequality shows that at least c_r of the auxiliary pairs miss
at most h groups of S. Therefore the entire pair count is at most c_r^-1
times the count with auxiliary support S' subset of S and |S\S'|<=h.
This comparison is made for each fixed kernel before averaging the outer.
The missing groups can differ between kernels; no common witness is chosen.

Use the retained sharper E-valued outer density from
[CHEAP_LANE_RANK_INTERFACE.md](../char2_field_spin/CHEAP_LANE_RANK_INTERFACE.md):

    beta = q0^160/(q0-1)^32.                                 (11)

For one nonzero E source group, its expected output counting measure is
dominated by beta times uniform E^256. Independent groups multiply this
bound. Equation (11) is not a density bound for F-valued inputs.

If S' misses j groups, apply (11) to its r-j nonzero groups. Complete its
allowed input set from I_t(S') to I_t(S) by fixing every added coordinate
to zero. Nonnegative sums imply a cost at most q0 per added coordinate.
Over the full word this costs q0^(256j). The constraints, weights, and
path-pruning rule are otherwise identical. In particular, pruning always
uses the actual I_t(S), not the possibly smaller auxiliary I_t(S').

Summing over the possible missing subsets gives the conservative factor

    beta^r C_r,
    C_r = sum_(j=0)^h binom(r,j) (q0^256/beta)^j.               (12)

Together, (10)--(12) justify markers chosen by actual r, with the explicit
penalty C_r/c_r. They do not justify applying old per-size markers to an
unlabelled auxiliary support.

There is a sharper bound when the missing-group strata are retained
through this density comparison. Put kappa=q0^256/beta>1. For one active
E group, let mu(v) be its expected output counting measure, summed over
nonzero messages. Injectivity gives mu(0)=0, and (11) gives
mu(v)<=beta q0^(-256). An auxiliary message missing exactly the groups
in J therefore has density at most beta^r kappa^|J| against the completed
uniform input measure, on the stratum

    v_g=0 exactly for g in J.

These strata are disjoint. After discarding the retained-kernel indicator,
all remaining equations, state variables, weights, and pruning depend
only on the completed input v and actual S, not on J. Thus their common
nonnegative reference function can be summed over these disjoint strata.
For |J|<=h the total is dominated by beta^r kappa^h times its unrestricted
uniform reference. In every subsequent bound one may consequently use

    C_r = kappa^h                                                (12a)

instead of (12). This is a sharper measure argument, not the false
numerical inequality sum_j binom(r,j) kappa^j <= kappa^h. It uses exact
group-zero strata and injectivity of the retained outer. If those facts
are discarded earlier, retain the conservative sum (12).

## 5. The two-state operator

Fix an allowed coordinate set I within one 64-symbol batch. Let x be
uniform in E^I and zero outside I; let u be independent uniform E^s.
Define four nonnegative moments:

    T_I(z) = E_x z^wt(x),
    K_I(z) = E_x z^wt(x) 1[Cx=0],
    J_I(z) = E_x z^wt((Id+AC)x),
    V_I(z) = E_(x,u) z^wt(x+Au).                              (13)

For nonempty actual I, the transition matrix and final vector are

                [ K_I          T_I ]
    M_I(z) =    [ q0^(-s) J_I   V_I ],
    f_I(z) = [ T_I, V_I ]^T.                                (14)

Rows specify the current state mark; columns specify the next mark.
The reset entry has a unique u=Cx and therefore carries q0^(-s).
The final vector imposes no reset or next-state condition.

For empty actual I, the structural pruning before (2) instead gives

    M_empty(z) = diag(1,V_empty(z)),
    f_empty(z) = [1,V_empty(z)]^T.                            (15)

For a fixed route and actual source set S, define

    P_S(z) = [1,0] M_(I_0(S))(z) ... M_(I_(n-2)(S))(z)
                         f_(I_(n-1)(S))(z).                 (16)

Every factor uses the actual allowed input set. The outer average and
the auxiliary support comparison above yield

    Pr[Bad_S | route]
      <= min(1, Gamma beta^r C_r z^(-e) P_S(z)
                       /[(q0-1)c_r]).                       (17)

The minimum in (17) is taken after averaging the outer; equivalently,
one may cap the fixed-outer bound and then apply Jensen. Both give the
fractional bound below.

For any 0<alpha<=1, define entrywise powers M_I^(o alpha) and
f_I^(o alpha). Positivity and subadditivity over state paths imply

    P_S(z)^alpha
      <= [1,0] product_(t=0)^(n-2) M_(I_t(S))(z)^(o alpha)
                                       f_(I_(n-1)(S))(z)^(o alpha).

Consequently the complete actual-size-r bound is

    Pr[there is a bad word on exactly r groups]
      <= [Gamma beta^r C_r z^(-e)/((q0-1)c_r)]^alpha
         * sum_(S:|S|=r) E_route[
             [1,0] product_t M_(I_t(S))(z)^(o alpha)
                              f_(I_(n-1)(S))(z)^(o alpha)].   (18)

The expectation is over the actual alternating pair permutations and
swaps. Their exact regional coefficient extraction remains necessary.
The choices z, alpha, and h may depend on r, but are fixed before setup.
Summing valid bounds (18) over r gives a whole-code bound. No such sum
has yet been certified for this operator.

## 6. What a bounded replay may reuse

The retained paired-route combinatorics can evaluate (18), provided its
local entries are replaced by (13)--(15). Its source prefactor must use
(10) and either (12) or (12a). Approximate moments need proved upward bounds.
Its exact terminal is f_I, not a transition with a charged last refresh.
For a simpler positive envelope, f_I<=M_I[1,1]^T entrywise, including
empty I. After entrywise alpha powers the same inequality follows from
subadditivity. A wrapper may therefore use n regular matrices and a final
all-ones vector. This adds a nonnegative envelope term; it does not assert
that an unused last refresh was sampled or contributes a probability cost.

The empty-input pruning removes artificial quiet-to-live transitions
across gaps that contain no coordinate from the actual source set.
Artificial transitions remain possible at a nonempty allowed set when
its auxiliary x happens to vanish. This loss is included in (14).
The current live-state moment also includes u=0. These are controlled
relaxations, not independence assertions about actual successive states.

The earlier positive operator has different activation and live entries.
Its numerical margins cannot be reused merely because the new refresh
has a smaller point probability. A pointwise operator comparison, or a
bounded directed replay of (18), is the next discriminating test.
The direct one-through-six-group bounds in
[SMALL_STATE_ROUTE_GATE.md](SMALL_STATE_ROUTE_GATE.md) remain available
independently and need not be reproved by this larger-source relaxation.

## 7. Penalize artificial zero states without coordinate masks

The auxiliary live state in (14) can equal zero. Repeated use of that
origin can enlarge a fractional bound. A second positive marker limits
this relaxation while retaining only two whole-state marks.

Fix a retained kernel V from Section 2. Let P be its E-valued pairs with
m!=0, and sample a pair uniformly from P. For every marked live Q_t,
the projection V -> Q_t has positive rank r_t. Otherwise an actual
F witness could not have Q_t!=0. Let j_t be the dimension of the intersection
of the zero-source fiber with ker Q_t. Rank-nullity gives j_t>=d0-r_t,
even though that fiber need not lie entirely in ker Q_t. Consequently

    Pr_P[Q_t=0]
      = (q0^(d-r_t)-q0^j_t)/(q0^d-q0^d0)
      <= q0^(-r_t) <= 1/q0 = 1/256.                         (19)

Let Z count zero Q_t among the a marked live states. Then E_P Z<=a/256.
Let G be the subset of P missing at most h source groups, as in Section 4.
Its relative size p satisfies p>=c_r. Therefore

    E[Z | G] <= a/(256 c_r).

For a real marker omega>=1, the function omega^(-x) is convex.
Jensen's inequality gives

    sum_(v in G) omega^(-Z(v))
      >= |P| p omega^(-E[Z|G])
      >= |P| c_r omega^(-a/(256 c_r)).

Equivalently,

    |P| <= c_r^(-1) omega^(eta_r a)
                    sum_(v in G) omega^(-Z(v)),
    eta_r = 1/(256 c_r).                                    (20)

No independence between missing source groups and zero states is used.
Equation (20) replaces the earlier unweighted source-support comparison;
it does not introduce a second factor c_r^(-1). The completed-input
function is still common to all exact missing-group strata. Thus the
sharper density factor beta^r kappa^h from (12a) still applies.

Put epsilon=q0^(-s). For nonempty actual I, replace (14) by

    Rbar_I = omega^eta_r epsilon (J_I-K_I+K_I/omega),
    Vbar_I = omega^eta_r
               (V_I-epsilon T_I+epsilon T_I/omega),
    Mbar_I = [[K_I,T_I],[Rbar_I,Vbar_I]],
    fbar_I = [T_I,Vbar_I]^T.                                (21)

In the reset entry, Q=Cx and Q=0 exactly when Cx=0. In the live entry,
the zero-Q contribution is epsilon T_I. These are the contributions
divided by omega in (21). The factor omega^eta_r is charged once per
marked incoming live state, including the last batch's incoming state.

For empty actual I the matrix remains diagonal:

    Mbar_empty = diag(1,omega^eta_r
                           (V_empty-epsilon+epsilon/omega)). (22)

Use fbar_empty=[1,Vbar_empty]^T. The final state is still free. As before,
fbar_I<=Mbar_I[1,1]^T is a valid positive terminal envelope.

Equation (18) remains valid with these new matrices and final vectors,
the same source prefactor, and C_r=kappa^h. Choose omega together with
z, alpha, and h by the actual source size r before setup. At omega=1,
the original operator is recovered. The modified entries need not be
at most one: they are positive counting envelopes, not transition
probabilities.

If a numerical implementation uses a majorant for J_I or V_I, the
subtractions in (21) must retain the exact K_I and T_I. Subtracting an
unrelated upper bound is not valid. The retained local provider has
exact K_I and T_I, so its proved return majorant can be used directly.
This refinement is an analytic bound only; no successful replay is
asserted here.

## 8. Joint weighting removes the missing-group cutoff

The same argument can weight missing source groups instead of discarding
a fraction of the auxiliary pairs. This removes h and c_r and keeps the
live-state exponent equal to 1/q0 for every actual source size.

Fix a retained kernel and its nonzero-source pair set P. Let J count
missing source groups and let Z count zero marked-live states. Equations
(9) and (19) give

    E_P J <= r/q0,             E_P Z <= a/q0.

For kappa=q0^256/beta>1 and omega>=1, joint convexity of the exponential
gives

    E_P[kappa^(-J) omega^(-Z)]
      >= kappa^(-E_P J) omega^(-E_P Z)
      >= kappa^(-r/q0) omega^(-a/q0).

Hence

    |P| <= kappa^(r/q0) omega^(a/q0)
                   sum_(v in P) kappa^(-J(v)) omega^(-Z(v)). (23)

This uses no independence between the two counts. As before, the actual
source set and actual state path determine which projections must be
nonzero before the auxiliary relaxation.

For an exact missing-group set of size j, its completed-input counting
density is at most beta^r kappa^j. The factor kappa^(-j) in (23) cancels
this excess pointwise. Distinct missing-group sets give disjoint strata,
so their total weighted measure is dominated by beta^r times the uniform
full-input measure. The source-zero stratum was excluded from P; adding
it to the uniform majorant only increases the bound. The remaining
function of completed input and auxiliary states is common across strata.

Thus a second complete interface is obtained from (18) by replacing its
source prefactor and live-state exponent as follows:

    Gamma beta^r C_r/[(q0-1)c_r]
        becomes Gamma beta^r kappa^(r/q0)/(q0-1),
    eta_r becomes 1/q0.                                    (24)

Use the weighted matrices and terminal in (21)--(22), with this new eta.
There is no missing-group cutoff and no retained-fraction denominator.
The parameters z, alpha, and omega can still be chosen by actual r before
setup. All source choices and the actual paired-route average remain
outside the corresponding fractional source prefactor, as in (18).

Equation (24) is especially useful when the minimal cutoff h made c_r
small. The earlier cutoff bound can nevertheless be better for some
sparse sizes, so both remain valid alternatives. A successful whole-code
claim still requires a complete directed replay of one of these bounds
for every remaining actual source class. The weighting argument itself
is not that replay.

## 9. A conditional route to intervals of source sizes

Exact regional coefficient extraction can be substantially sharper than
an independent-pair bound. Covering selected sizes does not, however,
cover the sizes between them. The following monotone envelope offers
one rigorous interpolation rule, at an unmeasured cost.

Set omega=1 and drop the special empty-input pruning. Use the raw matrix
from (14) also for empty I:

    Mraw_empty = [[1,1],[q0^(-s),V_empty]].                    (25)

This only enlarges the bound. Each exact entry of Mraw_I is nonincreasing
when the allowed input set I grows. For example, expand z^weight as a
positive sum of requested-zero indicators. Averaging each indicator over
uniform input gives q0 to the negative rank of its linear constraint
matrix. Adding an independent input coordinate adds a column, so that
rank cannot decrease. The same argument applies to K with its extra Cx=0
constraints, to J with output (Id+AC)x, and to V with its uniform state
columns. The local majorants need not themselves be monotone: the exact
moments at a larger I are still bounded by a majorant at a contained I.

Couple two source sets S0 subset S using the same complete sampled route.
Then I_t(S0) subset I_t(S) in every batch. Positivity and entrywise alpha
powers preserve the resulting matrix inequalities. Consequently, for
fixed z and alpha, the raw operator expectation for S is at most that
for S0. This is a comparison of actual nested source sets, not of
independently sampled regional histograms.

Suppose B_a(z,alpha) bounds the raw operator expectation for every source
set of size a, including all its realizable matching profiles. For
a<=r<=b, choose any a-subset of each actual r-set. The joint-source
correction (24) then gives the complete interval bound

    Pr[bad word on between a and b actual groups]
      <= B_a(z,alpha) [Gamma z^(-e)/(q0-1)]^alpha
         * sum_(r=a)^b binom(8192,r)
                         (beta kappa^(1/q0))^(alpha r).       (26)

All markers are fixed for the entire interval. The finite sum may be
bounded by its largest summand times b-a+1. For b<=4096 its summands
increase, so the endpoint b suffices. Above that range, either retain
the sum or verify its maximum. Do not silently use the lower endpoint's
source-density factor.

A bound on the average size-a operator also suffices in (26). Let O(S)
be the exact raw-operator expectation for source set S, and let Obar_a
be its average over all size-a source sets. For every size-r set S,
monotonicity gives

    O(S) <= average_(T subset S, |T|=a) O(T).

Average this inequality over all size-r sets. Every size-a set occurs
in exactly binom(8192-a,r-a) of them. Double counting therefore yields
Obar_r<=Obar_a. Thus (26) remains valid when B_a only bounds Obar_a.
The run census in Section 11 can give this averaged anchor by dividing
its route-only sum by binom(8192,a). All source factors in (26) still
use the actual r, and the raw-operator requirement is unchanged.

The pruned operator does not satisfy this monotonicity. Its empty-input
off-diagonal entries are zero, whereas adding an allowed coordinate
makes them positive. Therefore (26) cannot use an existing pruned
anchor without recomputing a valid raw-operator majorant.

No successful raw-operator anchor or interval receipt is asserted here.
The raw envelope permits artificial births throughout empty gaps; that
loss may remove the useful margin. Its source factor also grows by
roughly 1024 alpha bits per added group, before source-choice costs.
Thus an isolated margin near 1000 bits at alpha=1 does not even guarantee
coverage of the next size by this interpolation. A bounded comparison
of (25) against the pruned exact-coefficient operator is the prerequisite
to implementing an interval verifier, not a justified performance or
proof shortcut by itself.

## 10. Retain a source-zero cutoff in the two-state operator

A concentrated route profile can make the unrestricted uniform-input
reference count trajectories whose input constraints already force the
actual outer message to zero. The fixed-kernel proof excludes such
trajectories, but the later density majorant can put them back. A source
rank tag can preserve part of that exclusion without returning to a
coordinate-support state operator.

Fix actual source set S of size r, the outer and route, a whole-state path,
and requested output-zero sets Z as in Section 2. When the incoming state
is quiet, the output equations are pure input equations:

    X_t restricted to Z_t = 0.

If the next state is also quiet, add C X_t=0. Collect just these equations
into a source matrix B. Leave all live-input equations out of B; omitting
them only weakens its rank. This B is E-valued and is independent of the
inner refresh randomness.

Suppose a quiet batch has ell_b allowed input positions in emission block
b, of which a_b are requested zero. Its local input rank is exactly

    d_t = sum_b [a_b+min(4,ell_b-a_b)]        for 0 -> 0,
    d_t = sum_b a_b                         for 0 -> 1.        (27)

The requested coordinate zeros are independent. After imposing them,
the remaining columns of each RS feedback block have rank
min(4,ell_b-a_b). Set d_t=0 for live incoming transitions. These ranks
depend only on the allowed occupancies and the requested-zero counts.

Within each outer column, the corresponding input constraints involve
disjoint routed coordinates. Let d_l be the sum of its batch tags. The
invertible column mixers preserve rank, and any 32 RS columns are jointly
surjective onto their evaluations. Hence, for the original source space
of dimension D=128r,

    rank B >= top32(d_0,...,d_63),       0<=d_l<=4r.             (28)

If the right side reaches D, no nonzero source can satisfy the fixed
Z/path event. Discard this event before taking the fourth root in (4).
Equivalently, its kernel has d=d0 in the notation of Section 2. This
exclusion is uniform over every invertible outer-mixer realization; it
does not rely on a typical route or a density approximation.

The source cutoff is determined by the rank tags, so it can remain in
the positive zero-subset expansion and subsequent outer-density average.
It is common to all completed-input zero-group strata. The same source
corrections from Sections 4 or 8 therefore apply to the surviving terms.

The tagged quiet entries have explicit positive local polynomials. For
one emission block, let u record its input rank and define

    K_ell(z,u) = sum_(a=0)^ell binom(ell,a) z^(ell-a)(1-z)^a
                  q0^(-a-min(4,ell-a)) u^(a+min(4,ell-a)),
    T_ell(z,u) = sum_(a=0)^ell binom(ell,a) z^(ell-a)(1-z)^a
                  q0^(-a) u^a.                              (29)

Use product_b K_(ell_b) for 0->0 and product_b T_(ell_b) for 0->1.
The live row retains its existing positive majorants, tagged with rank
zero. Actual-empty pruning remains unchanged. At u=1 these are exactly
the prior quiet row, not new probability estimates. The zero-subset
choices outside the allowed input set sum to one and need no rank tag.
The exact terminal is the tagged vector [product_b T_(ell_b),V_I]^T:
its quiet component has output-zero ranks but no reset equations. A
regular-matrix terminal envelope remains possible only by explicitly
retaining this free branch. In a nonempty final batch it is the next-live
column; in an empty final batch it is the unchanged diagonal branch.
Removing full-source-rank branches must not turn the terminal into a
required reset.

The q0^(-d_t) coefficients in (29) already pay for the uniform-input
probability of the equations. Do not multiply by another q0^(-d_t) or
tau^(-d_t). The proposed improvement is solely the strict source cutoff
top32<D. No additional source-rank potential is being asserted.

For a first implementation take alpha=1. A much smaller counter than a
rank-degree polynomial suffices. Since d_l<=4r and D=32*4r,

    top32(d)<D  iff  at most31 columns have d_l=4r.             (30)

Each column contains exactly 4r allowed input coordinates. Thus d_l=4r
holds exactly when every local tag equals its allowed-input count.
Call such a local term clean; all other terms have a rank deficit.
Within one quiet-reset block, its clean coefficient is

    sum_(a=max(0,ell-4))^ell binom(ell,a)
            z^(ell-a)(1-z)^a q0^(-ell).                      (31)

Within a quiet-to-live block, only a=ell is clean, with coefficient
((1-z)/q0)^ell. Multiply the four block coefficients to obtain the
clean quiet entries. The remaining positive terms from (29) have a
deficit. A live incoming transition with nonempty I has a deficit,
because its retained tag is zero. Empty-I diagonal transitions are
clean; forbidden off-diagonal transitions remain zero.

Compose the deficit flag by Boolean OR through each batch, region,
and four-region column. Increment a global counter only on clean
columns; retain counts0,...,31 and discard count32. This implements
exactly the cutoff (30), with positive arithmetic and no union over
rank thresholds. Keep the free-terminal branch described above when
processing the last column.

Use one common output marker z before this expansion. Unlike the direct
rank gate, (27) depends on Z, not just the actual state path. Different
markers for clean and deficit terms are therefore not justified by the
present argument. Fractional alpha<1 requires applying subadditivity
to the tagged positive sum or another explicitly justified grouping;
it must not silently inherit an untagged fractional wrapper.

A simpler fallback tags only quiet-reset constraints: tag
c=sum_b min(4,ell_b) on 0->0, and zero elsewhere. Its cutoff depends only
on the state path and route. It loses the output-zero ranks in (27), but
does not have the same requested-zero-dependent marker restriction.

### Fractional grouping and the implemented control

Write C_I and D_I for the two-state clean and deficit matrices just
defined, and M_I=C_I+D_I. Extend the state by a flag recording whether
the current column already has a deficit. Its local matrix is

    W_I = [[C_I,D_I],[0,M_I]].                                (32)

The first flag value is clean. Once the flag is dirty, later clean and
dirty terms are grouped in M_I before taking a fractional power. Apply
entrywise alpha to (32) before route averaging. Subadditivity over paths
of this extended state gives the same type of bound as (18), with the
unchanged source prefactor. It would be incorrect to power M_I first
and then assign parts of that powered value to clean and dirty paths.

For one column, multiply its four regional matrices in even, odd, even,
odd order. Let A and B be the top-left and top-right two-state blocks of
this product. The column-count polynomial is B+u A. Reset the flag at
the next column and retain degrees zero through31 in (B+u A)^64. The
initial state is quiet and the final vector is all ones, with the
tag-preserving free branch described above. Every coefficient is
nonnegative. No coefficient is obtained by subtracting rounded totals.

The floating control is
[probe_small_state_column_cutoff.py](theory/probe_small_state_column_cutoff.py).
It computes exact rational quiet partitions, checks each regional integer
histogram mass, and uses log arithmetic for the clean-column counter.
At alpha=1, its regional clean-plus-dirty block agreed with the original
operator within 3.3e-16 relative error. An independent counter check
retaining all 65 possible counts agreed with the ordinary uncut matrix
power within 1.9e-11 in log2 units on the g17,(8,8) profile.

Selected concentrated profiles at g17,24,32 showed no useful gain from
this cutoff. Fractional splitting cost about 0.20 and 0.07 bits at the
latter two markers; the measured cutoff gain was below numerical noise.
Thus this refinement is not the basis of the directed closure through
g32. Those gates use the simpler untagged exact-route operator and adjusted
rational markers. Keep this source-cutoff implementation as a checked
diagnostic, not a certificate or a reason to add encoder operations.

## 11. Count source sets by their occupied runs

The worst matching profile can describe a rare arrangement of active
groups. Charging that profile for all binom(8192,g) source sets loses
information already available from the two fixed matchings. The following
refinement changes only the source-set count in (18).

Fix 1<=g<8192 and markers z, alpha, and a valid source correction.
The two matchings form one labeled cycle on L=8192 groups. For an actual
source set S of size g, let c be its number of nonempty occupied runs on
this cycle. The induced edges consist of g-c edges, so its ordered
matching counts satisfy

    f_even+f_odd=g-c.

The exact number of such source sets with c runs is

    N_g(c)=(L/c) binom(g-1,c-1) binom(L-g-1,c-1),
    1<=c<=min(g,L-g).                                        (33)

To count them, distinguish the start of one occupied run. Choose its
labeled vertex, then the positive lengths of the c occupied and c empty
runs. This counts every source set c times. Dividing by c gives (33).
In particular, these integer counts sum to binom(L,g).

Let R_(g,f)(z,alpha) be the exact regional matrix averaged over the
specified pair permutation and singleton orientations. Its local matrices
are powered entrywise before averaging. Define the free-terminal envelope

    O_g(f0,f1)=[1,0](R_(g,f0) R_(g,f1))^128 [1,1]^T.

For fixed S the route expectation depends only on its ordered pair counts.
The source prefactor in (18) depends only on g, not on S's locations.
For each c, take a maximum over all feasible ordered count pairs with
f0+f1=g-c. Including additional unrealizable profiles only increases it.
Then (18) gives

    Pr[bad word on exactly g actual groups]
      <= [Gamma beta^g C_g z^(-e)/((q0-1)c_g)]^alpha
         * sum_c N_g(c) max_(f0+f1=g-c) O_g(f0,f1).             (34)

The joint-source correction (24) can replace the prefactor in the same
way. Run counts remain outside alpha. There is no further binomial or
profile-count multiplier. Keep both orders of an asymmetric pair; matrix
multiplication need not commute.

Equation (34) applies to the pruned actual-empty matrices. It partitions
source sets directly and does not use the subset monotonicity from
Section 9. Markers are fixed for each g before setup. The full-support
endpoint g=L has one source set and the deterministic profile (4096,4096);
it is not covered by the proper-subset run formula.

## 12. Bound a complete interval of matching counts

For a fixed anchor g, large exact coefficient tables may be unnecessary.
One positive IID proposal can bound several full-pair counts at once.
This changes a route majorant, not the actual route distribution.

Fix positive rational probabilities p0,p1,p2 summing to one. In each
pair slot independently choose an empty, singleton, or full pair with
these probabilities. Orient singletons fairly in the actual striped
layout. Let L(p) be the expected entrywise-alpha raw matrix of one
32-pair batch. For full-pair count f, the required regional histogram is

    n(f)=(4096-g+f, g-2f, f).

Let H_f(p) be its exact multinomial probability under p. Conditioning
the IID placements on this histogram gives the actual region's pair
placement distribution. Positivity therefore gives the entrywise bound

    R_(g,f) <= L(p)^128 / H_f(p).                             (35)

For consecutive feasible counts, the histogram probabilities satisfy

    H_(f+1)/H_f = [(g-2f)(g-2f-1)/((4097-g+f)(f+1))]
                  * p0*p2/p1^2.

This ratio decreases with f. Thus log H_f is concave, and on any
integer interval [a,b] its minimum occurs at an endpoint. Replacing
the denominator in (35) by min(H_a,H_b) gives one regional matrix for
the entire interval. Partition all feasible f into such intervals.
The maximum over every ordered pair of interval matrices bounds every
source-set route operator at the anchor g. Including impossible interval
pairs is conservative. A collection of fixed proposals may be combined
by the entrywise minimum of their regional upper matrices, since each
bounds the same entries throughout that interval.

The [dense-anchor control](theory/probe_small_state_dense_anchor.py)
implements this interface with at most 128 intervals and five proposals.
It uses raw empty-input matrices, so Section 9 may extend a successful
anchor to larger source sizes with the same markers and their own source
factors. Its arithmetic is floating point. Complete matching-profile
coverage alone is not a directed certificate.

A common regional envelope can avoid enumerating all ordered interval
pairs. Let U be the entrywise maximum of valid regional upper matrices
over every interval. Every regional matrix is then at most U. Positivity
gives the source-set-uniform raw operator bound

    B_g(z,alpha) := [1,0] U^256 [1,1]^T.                     (36)

No commutation of different regional matrices is assumed. This replaces
each factor by the same entrywise upper matrix, in its original position.
For a directed implementation, take exact rational maxima of certified
upper endpoints. Within each interval, divide by a certified positive
lower bound on both endpoint histogram probabilities. The resulting
B_g can enter (26), with the actual source factors throughout that interval
of source sizes. Coverage of all matching intervals is still required.

## 13. Combine the cell census with averaged anchors

The common matrix in (36) can combine large entries from unrelated
matching-count cells. Keeping the ordered cell pair preserves that
information. The run census then charges clustered source sets for their
actual number, rather than for every source set of the anchor size.

Fix an anchor size 1<=a<L=8192 and theta=(z,alpha), with
0<z<1 and 0<alpha<=1. Throughout this section use omega=1, the exact
raw matrices from Section 9, and the positive free-terminal envelope.
Local entrywise powers precede route averaging. Write O_theta(S) for
the resulting exact operator expectation over the actual route for a
fixed source set S. Its monotonicity is the one proved in Section 9;
local numerical majorants need not themselves be monotone.

Partition every integer full-pair count in

    max(0,a-4096) <= f <= floor(a/2)

into disjoint cells C_i=[l_i,u_i]. For each cell, let E_i(theta) be a
nonnegative regional matrix that bounds R_(a,f)(theta) entrywise for
every f in C_i. Equation (35) supplies one such bound. A singleton cell
may instead use a regional coefficient bound at the same a and theta.
Define

    H_ij(theta) = [1,0](E_i(theta) E_j(theta))^128 [1,1]^T.    (37)

An integer pair in C_i times C_j satisfies f_even+f_odd=a-c exactly
when l_i+l_j<=a-c<=u_i+u_j. Thus this cell pair can contribute to the
run-count maximum precisely on the integer interval

    max(1,a-u_i-u_j) <= c <= min(a,L-a,a-l_i-l_j).

Additional cycle-realizability restrictions may exclude some count
pairs. Including them only enlarges the bound. For each
1<=c<=min(a,L-a), define B_c(theta) as the maximum of H_ij(theta)
over cell pairs whose interval contains c. Keep both orders: (37)
does not assume that the two regional matrices commute.

Every actual a-set with c runs has its ordered matching counts in one
of these cell pairs. Positivity therefore gives O_theta(S)<=B_c(theta).
Using the exact counts N_a(c) from (33), the averaged raw anchor obeys

    average_(S:|S|=a) O_theta(S)
      <= [sum_c N_a(c) B_c(theta)] / binom(L,a).              (38)

The numerator in (38) is route-only. It contains neither a source-density
factor nor an extra binomial or profile-count multiplier. All cells and
all run classes must be covered, and sum_c N_a(c)=binom(L,a) must hold.
Each ordered cell pair requires only one operator calculation; its value
can update every c in the interval above. This avoids a separate exact
regional coefficient table for every matching profile.

The FLOAT control
[probe_small_state_census_anchor.py](theory/probe_small_state_census_anchor.py)
implements this complete cell census with optional singleton tail cells.
Caller-supplied regional overrides must use the same raw operator,
anchor size, and markers. Such overrides are not independently certified
by the FLOAT interface. Directed use requires certified regional bounds,
upper maxima that cannot understate overlapping enclosures, and a complete
counted sum. The construction and its sampling distribution are unchanged.

### Markers chosen by the anchor subset's run class

Different run classes may favor different markers. They can still be
transported to larger source sizes, but their source factors must remain
inside the run-class sum.

Fix 1<=a<=b<=L with a<L. For each anchor run count c, choose
theta_c=(z_c,alpha_c) in the preceding domains before setup is sampled.
The choice depends on c and the fixed anchor, not on sampled setup.
Let B_c(theta_c) bound every a-set with c runs as above. Define the
per-source-set factor from the joint correction (24) by

    D_r(theta) = [Gamma beta^r kappa^(r/q0) z^(-e)
                            /(q0-1)]^alpha.                 (39)

This factor excludes the number of source sets. For an arbitrary actual
source set S of size r, the probability of Bad_S is over the ideal setup.
For every subset T of S with |T|=a, let c(T) count its occupied runs
and put theta_T=theta_(c(T)). The per-source-set bridge and raw
monotonicity give

    Pr_setup[Bad_S]
      <= D_r(theta_T) O_(theta_T)(S)
      <= D_r(theta_T) O_(theta_T)(T).

Average these valid upper bounds over all binom(r,a) choices of T.
Then sum over all size-r sets S. Every fixed size-a set T belongs to
binom(L-a,r-a) such sets S, and

    binom(L-a,r-a)/binom(r,a) = binom(L,r)/binom(L,a).

Grouping T by c(T), and then summing over r=a,...,b, proves

    Pr_setup[bad word on between a and b actual groups]
      <= (1/binom(L,a)) sum_c [N_a(c) B_c(theta_c)
                       * sum_(r=a)^b binom(L,r) D_r(theta_c)]. (40)

This is averaging over subsets of a fixed S, not selecting a favorable
marker after observing setup. If all theta_c agree, (40) reduces to the
averaged-anchor form of (26). When markers vary, do not factor a common
source term out of the c sum. Each D_r uses the actual r, never a in
place of r. The full-support endpoint r=L is allowed in (40); the run
census is taken only at the proper-subset anchor a<L.

### Optional deletion of singleton inputs

A regional bound may also discard singleton inputs before averaging.
This is a relaxation of the raw operator, not a change to the encoder
or to the actual source set.

Fix an actual source size g and one feasible regional histogram, with
g-2f singleton pairs and f full pairs. Couple it to the same placement
after deleting all singleton coordinates.
The remaining allowed set is contained in the original one in every
batch. Exact raw monotonicity and entrywise alpha powers therefore bound
the original regional operator by the operator after deletion. Under
the actual uniform pair permutation, the surviving f full-pair slots
are a uniform f-subset of the 4096 slots. Let Rfull_f(theta) denote
this full-pairs-only regional expectation. Then

    R_(g,f)(theta) <= Rfull_f(theta).                         (41)

For example, let Q_theta(y) be the one-batch matrix polynomial whose
coefficient of y^j sums the powered raw matrices over all j-full-pair
placements among its 32 slots. Exact univariate extraction gives

    Rfull_f(theta) = [y^f] Q_theta(y)^128 / binom(4096,f).     (42)

The same bound can be evaluated with a valid positive coefficient
majorant. Deletion is justified by the exact moments; no monotonicity
assumption on the numerical local caps is needed.

The retained full-pair count is not a new actual source size. A direct
size-g bound still uses D_g(theta) and the original source census.
Within (40), N_a(c) still counts the original anchor subsets and D_r
still uses the actual r. Deletion may be applied separately in each
region, so the retained coordinates need not form one common smaller
source set. Empty-input pruning cannot replace the raw operator in
this monotonicity argument.

Equations (38)--(42) are analytic bounds. This section asserts no new
numerical source-size interval or whole-code certificate. Such a claim
still requires complete directed evaluation at fixed witnesses.

## 14. Group two batches before the fractional power

Applying the fractional power separately to each batch splits the two
possible intermediate states. Retaining that sum for two consecutive
batches can reduce the bound without changing the encoder or route.
The following calculation evaluates this improvement without enumerating
all 4845 squared pairs of occupancy profiles.

Fix 0<z<1 and 0<alpha<=1. Use the exact raw matrices M_I from Section 9,
including raw empty-input transitions, and retain the batch order.
Write P^[alpha] for entrywise powering of a nonnegative matrix P.
For any two such matrices,

    (MN)^[alpha] <= M^[alpha] N^[alpha]                      (43)

entrywise, by subadditivity of the scalar alpha power. Group the 128
batches of each region into 64 consecutive pairs before applying that
power. Let R2_(g,f) be the resulting regional expectation under the
actual route with histogram n(f) from Section 12. This is a different,
no larger regional operator than the one-batch operator R_(g,f).
Subadditivity across the remaining group boundaries gives the same
source prefactor and free-terminal convention as before.

### Bound the two-path sum by bins of ratios

Let ell range over the 4845 sorted four-block occupancy profiles.
For each ell, fix a strictly positive matrix Mhat_ell that bounds the
exact raw matrix entrywise. The retained rational local majorants have
this property. Fix positive IID pair-type probabilities p as in Section 12,
and let pi_ell(p) be the complete striped occupancy probability of ell.
The two batches are independent under this unconditioned IID experiment.

For fixed output indices a,b in {0,1} and profiles ell,m, put

    x=(Mhat_ell)_(a,0) (Mhat_m)_(0,b),
    y=(Mhat_ell)_(a,1) (Mhat_m)_(1,b),
    lambda_a(ell)=(Mhat_ell)_(a,0)/(Mhat_ell)_(a,1),
    rho_b(m)=(Mhat_m)_(0,b)/(Mhat_m)_(1,b).

Thus x/y=lambda_a(ell)*rho_b(m). Define

    chi_alpha(r)=(1+r)^alpha/(1+r^alpha),
    (x+y)^alpha=chi_alpha(x/y)*(x^alpha+y^alpha).             (44)

The correction chi_alpha is at most one and satisfies
chi_alpha(r)=chi_alpha(1/r). Its logarithmic derivative is
alpha*(r/(1+r)-r^alpha/(1+r^alpha)), which is nonnegative for r>=1.
Consequently chi_alpha is nondecreasing in |log2 r|.

Fix a bin width w>0. Place a positive ratio r in bin
k=floor(log2(r)/w). If the left and right bins are h and k, their
product has log2 ratio in [(h+k)w,(h+k+2)w]. Set

    d_hk=w*max(|h+k|,|h+k+2|),
    b_hk=(1+2^(-d_hk))^alpha/(1+2^(-alpha*d_hk)).

Equation (44) therefore gives

    (x+y)^alpha <= b_hk*(x^alpha+y^alpha).

Using the closed upper endpoint is conservative, including when a ratio
lies exactly on a bin boundary.

For the fixed a,b, define the positive bin sums, for j in {0,1}, by

    A_(h,j)=sum_(ell: lambda_a(ell) in bin h)
                         pi_ell(p)*(Mhat_ell)_(a,j)^alpha,
    B_(k,j)=sum_(m: rho_b(m) in bin k)
                         pi_m(p)*(Mhat_m)_(j,b)^alpha.

Independence of the two IID batches factors their weights. Therefore
the matrix W(p) with entries

    W(p)_(a,b)=sum_(h,k) b_hk
                         [A_(h,0) B_(k,0)+A_(h,1) B_(k,1)]  (45)

bounds E_IID[(M_I M_J)^[alpha]] entrywise. Every occupancy pair appears
in exactly one bin pair. No profile or singleton orientation is omitted.
Since b_hk<=1, this bound is also at most the square of the one-batch
IID majorant formed from the same Mhat_ell. For alpha=1 the correction
is exactly one, and this comparison is an equality.

### Condition once, after the independent product

Under the unconditioned IID experiment, the 64 paired batches are
independent. Their expected ordered product is bounded by W(p)^64.
The full regional histogram event has probability H_f(p)>0 from
Section 12. Positivity and conditioning on that event give

    R2_(g,f) <= W(p)^64 / H_f(p).                            (46)

The batches are not independent after conditioning. Equation (46)
uses independence only before that step. Its denominator is the one
4096-pair histogram probability, without an alpha power or a separate
conditioning factor for each paired batch.

On a full-pair cell [l,u], replace the denominator by
min(H_l(p),H_u(p)), using the concavity proved in Section 12.
The resulting matrix can replace E_i in (37)--(40), with R2 as the
regional operator. The run census still uses the original source sets,
both orders of every compatible cell pair, and source counts outside
alpha. The source factor D_r still uses the actual source size r.

### Raw monotonicity survives grouping

Couple two source sets T contained in S by the same route realization.
The allowed input set from T is contained in that from S in each batch.
Section 9 gives M_(I(S))<=M_(I(T)) for the exact raw matrices.
Positive multiplication, entrywise alpha powers, and further positive
products preserve this order. Averaging over the shared route proves
the corresponding monotonicity for the exact grouped operator O2_theta:

    O2_theta(S) <= O2_theta(T).                              (47)

Thus the averaged-anchor argument in (38)--(40) also applies to grouped
operators. The numerical caps Mhat_ell and their bins need not be
monotone; evaluate them only to bound the operator at the anchor.
Empty-input pruning cannot be substituted in this argument.

The FLOAT control
[probe_small_state_grouped_route.py](theory/probe_small_state_grouped_route.py)
uses all 4845 raw local profiles, one positive rational proposal per
cell, and the complete census from Section 13. For w=1/v with
v in {1,2,4}, it computes each bin index as floor(log2(r^v)) using exact
rational comparisons. Powers, sums, and correction values remain floating
point. A directed version may instead start from certified positive
rational upper entries and bound every subsequent operation outward.
This section establishes the analytic majorant, not a numerical
source-size interval or a whole-code certificate.

## 15. Reuse regional numerators at each actual source size

The averaged-anchor bound keeps the smaller source set's route operator
while its source factor grows. Direct conditioning can retain the route
gain at each actual size without repeating the expensive IID calculation.
The reusable object is the unconditioned regional numerator, not a
conditional matrix or an averaged source-set bound.

Fix 1<=a<=b<L=8192 and common markers theta=(z,alpha) in the domains
of Section 14. Retain that section's grouped raw operator and free
terminal. Choose finitely many proposal families before setup. Within
each family, fix positive rational triples p_i summing to one, and define

    K_i = W_theta(p_i)^64,
    C_ij = [1,0](K_i K_j)^128 [1,1]^T.                      (48)

These matrices and scalars do not depend on the actual source size.
At each g in [a,b], use one predetermined family. Suppress its family
index below. Associate its active indices with nonempty integer cells
I_i(g)=[l_i(g),u_i(g)] that partition every feasible count

    max(0,g-4096) <= f <= floor(g/2).

The cells may move with g; the proposals p_i remain fixed within their
family. Let H_(g,f)(p) denote the multinomial probability of
n(g,f)=(4096-g+f,g-2f,f) under p. Choose positive numbers h_i(g) satisfying

    h_i(g) <= min(H_(g,l_i(g))(p_i), H_(g,u_i(g))(p_i)).

Section 12's concavity in f makes h_i(g) a lower bound throughout the
cell. Write O2_g(f0,f1) for the ordered free-terminal operator
[1,0](R2_(g,f0) R2_(g,f1))^128 [1,1]^T. Equation (46) and positive
multiplication therefore give

    R2_(g,f) <= K_i/h_i(g)                         for f in I_i(g),
    O2_g(f0,f1) <= C_ij/[h_i(g)h_j(g)]^128
                         for f0 in I_i(g), f1 in I_j(g).     (49)

Each scalar denominator occurs 128 times in that product.
No alpha power applies to the histogram probabilities,
and factoring these scalars does not commute the regional matrices.

For 1<=c<=min(g,L-g), define

    B_g(c) = max_(i,j: l_i(g)+l_j(g)<=g-c<=u_i(g)+u_j(g))
                         C_ij/[h_i(g)h_j(g)]^128.            (50)

Keep both orders of each compatible pair. Every size-g source set with
c occupied runs satisfies the matching-count relation in (50). The
per-source-set bridge bounds its bad-word probability by D_g(theta)B_g(c),
where D_g is (39) evaluated at the actual g. Summing the exact source
counts N_g(c) from (33), and then the source sizes, proves

    Pr_setup[bad word on between a and b actual groups]
      <= sum_(g=a)^b D_g(theta)
                         * sum_(c=1)^min(g,L-g) N_g(c)B_g(c). (51)

There is no additional binomial or profile-count multiplier: N_g(c)
already partitions all binom(L,g) source sets. There is also no anchor
binomial divisor. Unlike (40), this argument does not average over smaller
subsets or use inclusion monotonicity. The full-support endpoint g=L is
excluded here; it has the separate deterministic profile from Section 11.

### Moving cells and directed replay

The retained implementation uses one family at each parity's first size,
a or a+1 when present. At later g of the same parity, translate the
family's cells by (g-anchor)/2 full pairs. Set the lowest cell's lower
endpoint to max(0,g-4096), clip every cell to the feasible range, and discard empty
cells. Require exact adjacency and complete coverage after this operation.
A fixed positive proposal remains valid on the extended lowest cell;
only its endpoint histogram bound changes.

Away from that extension, translation keeps the deficit
j=floor(g/2)-f fixed. For a singleton tail cell, g-2f=(g mod 2)+2j is
therefore constant within its parity family. Separate families avoid
using a zero-singleton proposal's tiny pseudocount for a one-singleton
histogram. This is a numerical choice, not a mathematical requirement
of (48)--(51).

Cache all ordered C_ij, including pairs incompatible at the family anchor.
Recompute h_i(g), compatibility, and the maximum ordering at each actual
g. The scalar denominators can change that ordering. In directed replay,
use certified regional upper matrices and scalar upper endpoints in (48),
positive lower endpoints for h_i(g), and certified upper maxima in (50).
Check sum_c N_g(c)=binom(L,g) exactly for every g. Accept only after all
cells, run classes, and requested sizes are covered and the complete
sum in (51) is below the exact target.

The [FLOAT control](theory/probe_small_state_grouped_interval.py) and
[directed verifier](theory/verify_small_state_grouped_interval.py) implement
this reuse without changing the encoder. Fixed witnesses and numerical
coverage are recorded in the [exact-route checkpoint](SMALL_STATE_EXACT_ROUTE_GATE.md).

## 16. Retain longer windows with a projective grid

Longer windows retain more intermediate state sums before the fractional
power. A finite grid can bound their products without enumerating every
sequence of occupancy profiles.

Fix theta=(z,alpha) in the preceding domains and the strictly positive
rational cap matrices Mhat_ell from Section 14. Fix one positive IID
proposal p, with complete batch probabilities pi_ell(p). Choose a window
length k in {4,8,16,32,64,128}; thus k divides the 128 batches of a region.
All these choices precede setup.

Choose a finite rational grid G={xi_0<...<xi_m} of positive numbers.
Its first and last nodes bound every row ratio
(Mhat_ell)_(a,0)/(Mhat_ell)_(a,1), for a in {0,1}.
Require xi_(i+1)/xi_i<=rho for a fixed rho>1. Write d(r) for the least
index with xi_(d(r))>=r. A positive row has the form s(r,1). Its upward
rounding is s(xi_(d(r)),1), which enlarges it coordinatewise.

For a grid node xi_i and a profile ell, define

    b_(i,ell) = xi_i (Mhat_ell)_(0,1)+(Mhat_ell)_(1,1),
    r_(i,ell) = [xi_i (Mhat_ell)_(0,0)+(Mhat_ell)_(1,0)]
                                                       /b_(i,ell). (52)

Multiplication sends s(xi_i,1) to s b_(i,ell)(r_(i,ell),1).
For every nonnegative incoming ratio, the new ratio lies between the
two row ratios of Mhat_ell. Consequently every initial and subsequent
ratio has a ceiling in G.

Define the nonnegative matrices nu, T, and F by

    nu_(a,i) = sum_(ell: d((Mhat_ell)_(a,0)/(Mhat_ell)_(a,1))=i)
                            pi_ell (Mhat_ell)_(a,1)^alpha,
    T_(i,j) = sum_(ell: d(r_(i,ell))=j) pi_ell b_(i,ell)^alpha,
    F_(i,b) = sum_ell pi_ell
                   [xi_i (Mhat_ell)_(0,b)+(Mhat_ell)_(1,b)]^alpha.
                                                                  (53)

Here nu has two rows, T has one row and column per grid node, and F has
two columns. The first batch initializes a rounded row. The next k-2
batches update its scale and grid index. The last batch is not rounded.
Thus

    W_(k,G)(p) := nu T^(k-2) F
      >= E_IID[(Mhat_(ell_1)...Mhat_(ell_k))^[alpha]]        (54)

entrywise. To prove (54), fix an arbitrary sequence of k profiles.
Positive multiplication preserves the upper row, and each ceiling only
enlarges it. Taking alpha powers preserves the endpoint inequality.
Averaging these positive terms gives (53)--(54), with exactly k factors
of pi. No intermediate quiet/live path is split inside the window.

Each ceiling increases a row by at most a factor rho. There are k-1
ceilings, so the reverse comparison to the same unrounded cap product is

    W_(k,G)(p) <= rho^(alpha(k-1))
                   E_IID[(Mhat_(ell_1)...Mhat_(ell_k))^[alpha]]. (55)

Across n=32768 batches, this gives at most
alpha n(1-1/k) log2(rho) bits of rounding loss relative to that unrounded
IID grouping. It does not bound the separate losses from local caps or
histogram conditioning. Different grids need not give ordered bounds,
even when one refines another: the rounded trajectories can differ.

Let Rk_(g,f) be the exact raw regional operator obtained by grouping each
k consecutive batches before entrywise powering. Independent windows
under the IID proposal and one conditioning step give

    Rk_(g,f) <= W_(k,G)(p)^(128/k) / H_(g,f)(p).             (56)

The denominator is the full 4096-pair histogram probability, without an
alpha power. On a complete full-count cell, use its two endpoint minimum
as in Section 15. These regional matrices enter the cell census and
actual-size source accounting unchanged. The original refresh factor,
source-density correction, and free-terminal envelope also remain unchanged.
Exact raw inclusion monotonicity survives the fixed window partition;
the numerical caps and rounding map need not themselves be monotone.

The [projective-window control](theory/probe_small_state_window_route.py)
uses FLOAT arithmetic and FLOAT ceiling assignment. Directed replay must
certify the ceilings and all positive arithmetic, cover every profile
and cell, and retain one full histogram divisor per region.
The [diagnostic note](SMALL_STATE_BRIDGE_SCREEN.md#actual-size-reuse-and-longer-windows)
records the selected witnesses and checks; those FLOAT results are not
a directed source-size certificate.

## 17. A common linear functional as a diagnostic

A common positive functional offers a smaller alternative to the grid,
but takes a worst case over every incoming row direction at each batch.
This can lose the direction information retained in (53).

Fix eta>0 and put w_eta=(eta,1)^T. For the same cap matrices and an IID
proposal p, define

    A_ell = (Mhat_ell w_eta)_0/eta,
    B_ell = (Mhat_ell w_eta)_1,
    f_(eta,p)(t) = sum_ell pi_ell(p)
                              [t A_ell+(1-t)B_ell]^alpha,
    lambda_(eta,p) >= max_(0<=t<=1) f_(eta,p)(t).            (57)

For any nonzero nonnegative row v, set t=eta v_0/(v w_eta).
Then E_IID[(v Mhat_ell w_eta)^alpha] is at most
lambda_(eta,p)(v w_eta)^alpha. This bound is uniform in v, so it can be
iterated through independent batches without splitting state paths.

The function f_(eta,p) is concave. For any fixed t0 in [0,1], its tangent
therefore gives the global upper bound

    max_t f_(eta,p)(t)
      <= f_(eta,p)(t0)
           +max(-t0 f'_(eta,p)(t0), (1-t0)f'_(eta,p)(t0)).   (58)

Here f'(t)=alpha sum_ell pi_ell(A_ell-B_ell)
[t A_ell+(1-t)B_ell]^(alpha-1). The signed differences use the same
explicit cap matrices that define f; they do not subtract unrelated
upper bounds. If directed evaluation gives f(t0)<=U0 and f'(t0) in [L,U],
use U0+max(-t0 L,(1-t0)U,0) in (58).

For a fixed source set S, the bridge before state-path subadditivity
bounds its bad-word probability by D_|S|(theta) times

    E_route[([1,0] product_t Mraw_(I_t(S)) [1,1]^T)^alpha].

This follows from (17), the joint correction (24), and the raw terminal
envelope. For each regional cell i, let h_i be a positive histogram
lower bound and lambda_i=lambda_(eta,p_i). Conditioning the functional
bound once per region gives the factor lambda_i^128/h_i, uniformly in
the incoming row. The actual regions are independent for fixed S.
For an alternating cell pair i,j, composition therefore gives the scalar
route bound

    B_(eta,i,j) = max(eta,1)^alpha
                     [lambda_i^128 lambda_j^128/(h_i h_j)]^128. (59)

The endpoint factor follows from [1,1]^T<=max(1,1/eta)w_eta and
[1,0]w_eta=eta. Use one common eta throughout each composition in (59).
A fixed finite pool permits taking the minimum over eta for each cell
pair before the run-class maximum. Independently minimizing each cell's
lambda over eta would not justify that composition.

The [common-functional control](theory/probe_small_state_functional_route.py)
implements this valid alternative, but the retained complete diagnostic
did not close its source class. The
[diagnostic note](SMALL_STATE_BRIDGE_SCREEN.md#actual-size-reuse-and-longer-windows)
records that failure. This scalar relaxation is not a replacement for
the successful FLOAT projective-window witness or a directed certificate.
