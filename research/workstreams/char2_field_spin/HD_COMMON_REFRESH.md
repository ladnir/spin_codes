# Applying one refresh to state plus feedback

The 5.37-ms carried-state experiment suggests reusing the feedback circuit
for the incoming state. A small change makes its conditional distribution
easier to analyze: add the state before the first diagonal, rather than
after it. This note proves a local operator for that new construction.
It does not claim a whole-code distance certificate.

## Construction

Keep the GF65536 outer and regional route of the certified checkpoint.
Write q=65536 and s=16. Partition each 64-symbol input packet into four
blocks. The emission A is the direct sum of four normalized RS[16,4]
generators, and C=A^T. The fixed interpolation matrix H is unchanged.
For each packet, sample independent nonzero diagonal entries in D1,D2,
independently of all other setup coins. Setup is sampled once and fixed.
Starting with Q=0, compute

    Y = X + A Q,
    Q_next = T (Q + C X),       T = D2 H D1.

The terminal state is unrestricted. The identity feedthrough makes the
inner invertible. This recurrence contains no rank-one state mixer.
The scalar maps remain GF65536-linear.

For each local RS block, C A=0. Indeed, products of its polynomial basis
functions have degree at most six. Summing such a polynomial over the
sixteen-point additive subspace gives zero. Thus (I+A C)^2=I in
characteristic two. In particular, I+A C is invertible, both locally and
over the full packet.

## Local weighted operator

Fix the potential input positions U and draw X uniformly from F^U in
the existing outer comparison experiment. Define K,T_j,V as in the
[block-emission note](HD_BLOCK_EMISSION.md). Put d=q^s-1 and define

    J_U = E_X[1[C X != 0] z^wt((I+A C)X)].

The incoming weighted measure has mass at most Z at zero and density
at most B/d on each nonzero state. Conditional on Q and X, the argument
Q+CX of the refresh is fixed. If it is zero, the next state is zero.
Otherwise the diagonal-Cauchy lemma bounds every nonzero target's
probability by b/d, where b=d/(q-1)^s. The refresh is invertible, so it
cannot map a nonzero argument to zero.

For an incoming zero state, the outgoing zero mass is K. The outgoing
nonzero density coefficient is at most b(T_j-K). For a flat nonzero
incoming state, a return to zero requires exactly Q=CX!=0. Its weighted
mass is J_U/d. Removing these returning pairs from the uniform-state
sum leaves V-J_U/d. Hence the exact-J comparison operator is

    [[K,       b(T_j-K)],
     [J_U/d,   b(V-J_U/d)]].

There is no quiet identity branch. The price is a return moment different
from the rank-one construction's return moment. Applying an upper bound
to J_U in the subtracted entry would be invalid. For a computable cap
Jbar>=J_U, we instead use

    b * [[K,        T_j-K],
         [Jbar/d,   V]].

Positive products and the existing outer comparison preserve this bound.
The packetwise density cost is log2(b), not the additional conditioning
cost paid for the old rank-one mixer. Full row coverage is still required.

## A computable return cap

In a local block, let l be the potential input count and let k_l be its
quiet moment. The image of F^U under I+A C has dimension l. An information
set therefore bounds its total marked mass by t_l=(z+(1-z)/q)^l.
Quiet words are unchanged by this map and contribute k_l.

If CX!=0, the RS emission has at least thirteen nonzero coordinates.
Adding X can cancel at most l of them. The nonquiet probability is
1-q^-min(4,l). Consequently the local nonquiet return moment is at most

    j_l = min(t_l-k_l, (1-q^-min(4,l)) z^max(0,13-l)).

The blocks are independent in the comparison experiment. For occupancies
l_1,...,l_4, use

    Jbar = min(product_i(k_li+j_li)-product_i k_li, d V).

The implementation evaluates the product difference as a positive sum.
The second cap follows because the nonnegative uniform-state sum d V
contains every returning pair. Both bounds concern the weighted return
event, not merely the probability of CX=Q.

Exact GF4 enumeration checks every potential input support, both incoming
measure components, and the complete diagonal ensemble. Separate exact
checks confirm C A=0 for the concrete four- and eight-check maps. The
FLOAT paired-regional screen is not directed arithmetic and is not a
whole-code certificate.

## Refining the return moment

The minimum-weight cap loses most of the return map's typical weight.
For example, the sixteen one-coordinate input subspaces produce only
weight-fifteen or weight-sixteen nonzero words under I+A C. Their average
marked mass is therefore much smaller than a weight-twelve cap at useful
marks. A tighter enumerator is preferable to assuming uniform output.

For a local support U of size l, expand the output mark over sets of
coordinates constrained to zero. Split such a set into D inside U and
E outside U, with sizes c and e. Put rho=min(4,l). If e>=4, the MDS
emission forces CX=0, so its nonquiet contribution is zero.

For e<4, let nu_c be the maximum nullity of a c-by-c principal submatrix
of P=I+A C. A left-kernel vector of P[D,U] maps injectively to a vector
annihilating C[U minus D]. Thus its nullity is at most
min(nu_c,max(0,4-l+c)). This gives the probability cap

    q^-max(0,c-min(nu_c,max(0,4-l+c))).

There is a second cap. Set R=CX. The constraints on E put R in a space
of dimension at most min(rho,4-e). After fixing nonzero R, the active
output-zero constraints fix X[D]=A[D]R. The remaining equation has at
most q^max(0,l-c-4) solutions. Consequently use the additional cap

    q^(max(0,l-c-4)-l) (q^min(rho,4-e)-1).

Also cap by 1-q^-rho. Multiplying the minimum of these three caps by
binom(l,c) binom(16-l,e) z^(16-c-e)(1-z)^(c+e), then summing over c,e,
gives a valid local return bound. Every summand is nonnegative.

The principal-nullity scan covers all 65536 coordinate subsets exactly.
The determinant lemma reduces a principal nullity to that of I+C[D]A[D],
a four-by-four matrix. An independent test computes all minors of those
matrices and reproduces the Gaussian-elimination ranks. There are 704
singular principal restrictions; the maximum nullity is two. These are
finite facts about the specified maps, not sampled rank assumptions.

The resulting FLOAT screens at active rows256/512/2048 are approximately
+130/+298/+5185. They are failed upper bounds, not counterexamples.
The local proof therefore does not yet close the whole code.

An exact enumerator now considers the lattice of spans of rows of P[:,U].
Each span determines all output coordinates forced to zero and a kernel
size q^(l-rank). Descending subtraction through the lattice counts vectors
with each exact zero set. This computation uses ranks over GF256, which
remain valid over GF65536. Small-field exhaustive enumeration independently
checks the method.

The experiment averages these enumerators over all U of a fixed small
size. It may average the local operator before the paired fractional
moment: conditional on the two occupancies, their supports are independent,
so E[P_U P_V]=E[P_U]E[P_V], and entrywise concavity gives
E[(P_U P_V)^alpha]<=(E[P_U]E[P_V])^alpha for 0<alpha<=1. This Jensen step
must be retained when using average returns rather than pointwise caps.

## A possible fixed mixer between feedback blocks

One proposed modification replaces C by (B tensor I_4)C, where B is a
fixed four-by-four Cauchy matrix. Its columns mix the four feedback blocks.
The identities C A=0 and (I+A C)^2=I still hold for the modified C.
The implementation has not yet been changed or timed for this proposal.

If w feedback blocks are nonzero, at least 5-w mixed blocks are nonzero.
Apply the scalar MDS branch bound to any nonzero coordinate of those
four-dimensional blocks. This gives a first weighted return cap.

A sharper cap uses input ranks rho_i=min(4,l_i). For any set of k mixed
blocks required to be zero, the associated constraint map has rank at
least the sum of the k largest rho_i. Indeed, restrict to those k input
blocks; the corresponding Cauchy minor is invertible. Subtract the
all-quiet probability q^-sum(rho_i). Then sum over nonempty mixed support
sets W, with output-weight floor sum_(i in W) max(0,13-l_i).

This refinement gives FLOAT screens about -35.6/-34.6 at rows256/512,
but still fails badly at row2048. It is not a complete certificate and
does not justify reporting the unmixed kernel's time for the mixed design.

## Native checkpoint

The unmodified common-refresh mode is `bucket4-common`. All 35 native
research variants pass Release and ASAN/UBSAN checks. Two 501-call medians
are 5.318 and 5.328 ms, against a 5.777-ms frozen certified control.
The checksum is 15132924349672518052. The validated executable SHA256 is
`8fe2508c7b3e4d85cebe0ca288713c3613a00861636903d6d3246c5023efc65b`.
No production profile or paper result has changed.

## Scattering check coordinates without an extra mixer

Mode36, `bucket4-scatter`, changes the feedback map without adding a
field multiplication stage. Keep the block emission A. Evaluate each
local feedback polynomial at nodes16,...,31 instead of0,...,15. Then
transpose the four-by-four array of feedback coordinates: check j from
input block i becomes coordinate i of state block j. Denote this map by
C_sc. The recurrence is Q_next=D2 H D1(Q+C_sc X).

The node shift is an invertible local change of polynomial basis. The
coordinate transpose is a permutation. Thus ker(C_sc)=ker(C), and
C_sc A=0. The quiet moment K and the invertibility of I+A C_sc remain
unchanged. V depends only on A and is unchanged as well. The return
moment must be recomputed for P_sc=I+A C_sc.

Every shifted check is nonzero on a single input coordinate. The coordinate
transpose therefore sends that feedback into all four emission blocks.
Exact column enumeration gives return weights56,...,63 for a one-coordinate
input, instead of15 or16 in the unmodified common-refresh construction.
This improves a specific failure mode; it does not prove the full code.

For each subset Z of the four checks, the screen enumerates the row-space
lattice of their sixteen evaluation columns. It obtains the exact minimum
rank r_Z(l) over all l-column restrictions. If the mixed state blocks in
Z are zero, the independent input blocks impose total rank at least
sum_i r_Z(l_i). Subtracting the all-quiet probability gives an upper bound
on each nonempty mixed support W. Mark that support by the floor
sum_(i in W) max(0,13-l_i), then sum over W. The labelled occupancies are
averaged over their permutations before the paired Jensen comparison.

For one or two potential inputs in the whole packet, the screen instead
computes the exact average return enumerator. The two-input computation
counts projective directions: each nonzero output row excludes at most
one finite ratio, and all exceptional ratios lie in GF256. It counts the
remaining GF65536 ratios together. No random sampling is used for these
enumerators. Larger potential sets still use upper bounds.

The native circuit changes the constants in two fixed FFT stages and omits
the feedback record transpose. Release and ASAN/UBSAN checks pass for all
36 modes, including the independent dense-map digest 1225673016892679168.
Two 501-call timings are
5.324/5.349 ms, around a 5.369-ms common-refresh control and 5.797-ms
certified control. This is approximately the same cost as mode35, not a
new certified timing. The first FLOAT screens give -50.1/-39.5/+23267
at rows256/512/2048. They improve the sparse return analysis but do not
close the whole code. The medium-density return bound remains too loose.

## Retaining the ranks of inactive output-zero constraints

The return word is X+A C_sc X. At an input-inactive coordinate, only the
second term remains. The preceding bound retains which emission blocks
are nonzero, but loses the ranks of individual zero constraints within
those blocks. The following cap retains these ranks.

Fix a labelled input support U, with l_i potential positions in input
block i, for i=0,...,3. Put j=sum_i l_i and rho=sum_i min(4,l_i). For
X uniform on F^U, the syndrome R=C_sc X is uniform on a rho-dimensional
space. Conditional on R, the input is uniform on an affine space of
dimension j-rho. Its active-coordinate output is a translation of that
space. An information set bounds its marked expectation by
t_(j-rho)=(z+(1-z)/q)^(j-rho), for every R.

Choose e_i inactive output coordinates in emission block i and require
them to be zero. Set a_i=min(4,e_i). In the natural polynomial-degree
order, the first a_i columns of the local RS evaluation matrix have full
rank on these rows. Row elimination therefore gives pivots 0,...,a_i-1.
For p=0,...,3, define J_p={i:a_i>p}. The equations with pivot p have no
coefficients on earlier input blocks. Their coefficients on input block
p are the shifted checks in J_p. Consequently their combined rank is at
least

    R_lower(l,a) = sum_(p=0)^3 r_(J_p)(l_p).

The block-triangular structure justifies adding these four ranks. The
finite tables r_J are the exact minimum check ranks already defined.
The probability of these zero constraints and R!=0 is at most
q^(-R_lower)-q^(-rho). The subtraction removes the all-zero syndrome;
it is nonnegative because R_lower<=rho.

Expand the inactive output mark as a positive sum over zero sets. Weight
each count tuple e by product_i binom(16-l_i,e_i)
z^(16-l_i-e_i)(1-z)^e_i. Multiply the resulting sum by t_(j-rho).
This gives a return-moment cap for each labelled occupancy. Counts e_i>=4
may be grouped because their constraint space is the full local state.

The implementation provides two versions. The first uses the worst
labelling and the worst check subset of each size. The second averages
the labelled bounds over all permutations of the occupancy profile,
retaining the particular check sets J_p. This average precedes the
paired Jensen comparison; it is not asserted to hold pointwise for
every labelling. Both versions leave the encoder unchanged.

An initial prototype incorrectly used the best input-block order without
changing the emission basis. Ordinary RS evaluation does not guarantee
that arbitrary square minors are nonzero. For example, input coordinate
48 and output-zero coordinates {0,1} give rank zero in A C_sc, whereas
that prototype asserted rank one. Its approximately -310.4 screen at
512 active rows is invalid and discarded. The corrected worst-labelling
screen gives -133.99 there. Neither number is a directed certificate.

Five regression tests check the counterexample, 600 concrete restricted
matrix ranks, exact one/two-input return enumerators, full-rank syndrome
profiles, and domination of the averaged cap by the worst-labelled cap.
These finite checks supplement the rank argument; they do not replace
whole-code coverage or directed arithmetic.

## A rank-budget obstruction to further feedback pruning

Reducing the total feedback rank below thirteen cannot achieve the target
with this identity-feedthrough construction, regardless of state mixing.
To see this, let every packet feedback map have rank at most f. There
are P=N/64 packets. Restrict messages to any r outer blocks, giving a
message subspace of dimension128r. All packet feedback constraints
together have rank at most fP. If128r>fP, a nonzero message makes every
feedback zero. Starting at zero, the recursive state then stays zero.
Its output is the routed outer word, supported on at most256r positions.

At K20, P=32768. Taking r=256f+1 proves

    d_min/N <= f/32 + 1/8192.

For f=12 this is below 37.52%, ruling out 40% distance. This is a deterministic
obstruction, not a failed probability bound. It applies to the uniform
rank cap f on every packet, and does not exclude the existing sixteen-rank
block feedback or non-identity feedthrough designs. It prevents pursuing
a superficially cheap four- or eight-rank global feedback map at this target.
