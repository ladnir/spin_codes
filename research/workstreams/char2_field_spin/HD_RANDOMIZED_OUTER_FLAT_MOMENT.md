# Keeping state constraints inside occupied tiles

The three-component state bound previously used the state polynomial only
on untouched tiles. Its flat-state emission bound consequently lost useful
constraints when every tile had input. This note supplies a stronger local
bound that also retains constraints from occupied tiles.

The new bound is never weaker than the previous bound. A positive dynamic
program computes all 65 packet occupancies with at most 744 retained states.
For example, at weight mark z=1/256 and 32 potential input positions, its
base-two logarithm improves from about -347.00 to below -510.54.

These are local moment bounds for the existing construction. They are not a
whole-code distance certificate or a change to the encoder.

## The moment and its probability space

Use the shared beta64 maps over F=GF65536, with q=65536. A packet contains
eight tiles of width eight, and its old state has four coordinates. Denote
the feedthrough by F_64 and the state injection by A. The emitted packet is

    Y=F_64 X+A h.

Fix a potential input set U, sample its coordinates X_U uniformly in F^U,
and set all other input coordinates to zero. Independently sample h uniformly
from F^4 excluding zero. For 0<z<1, the desired flat-state moment is

    V_U(z)=E_(X_U,h)[z^wt(Y)].

Let V_j(z) average this moment over all j-subsets U of the packet. This is
the flat-component emission moment in
[HD_RANDOMIZED_OUTER_STATE_PROFILE.md](HD_RANDOMIZED_OUTER_STATE_PROFILE.md).
It concerns the artificial uniform input measure after outer domination,
before output weighting or state conditioning.

For an output-zero mask Z, let

    r_0=rank((F_64)_(Z,U)),
    r_h=rank([(F_64)_(Z,U),A_Z])-r_0.

The probability that every output in Z vanishes is exactly

    q^(-r_0) (q^(4-r_h)-1)/(q^4-1).

Indeed, eliminate the input variables. The remaining equations have rank
r_h on h, so q^(4-r_h)-1 nonzero states survive. Each surviving state has
q^(|U|-r_0) input solutions. This count removes the zero state exactly.

Expand the weight mark as

    z^wt(Y) = sum_(Z subset [64])
                  z^(64-|Z|)(1-z)^|Z| 1[Y_Z=0].

Every coefficient is nonnegative. Therefore any lower bound ell<=r_h gives
a valid upper moment by replacing r_h with ell in the displayed probability.
When ell=4, the contribution is zero. No subtraction of two upper bounds is
involved.

## An exact rank calculation within one tile

Fix one tile, with a potential input positions and v selected output zeros.
Its old evaluation nodes and new evaluation nodes lie in disjoint additive
cosets. On these two cosets, the feedthrough entries are a common nonzero
multiple of 1/(y-u), up to invertible column scaling.

The state polynomial g has degree below four. Thus the possible emitted
values at the new nodes are evaluations of rational functions

    g(X) + sum_(u in U_tile) c_u/(X-u).

Let D_U(X)=product_(u in U_tile)(X-u). The numerator of this rational function
ranges over every polynomial of degree below a+4. Partial fractions and
polynomial division establish this bijection. None of the selected output
nodes is a pole. Consequently

    rank([F_U,A]_Z)=min(v,a+4),
    rank((F_U)_Z)=min(a,v),
    dim(E_tile)=min(4,max(v-a,0)),

where E_tile is the space of linear constraints left on the common state
after eliminating that tile's input variables.

The tile input variables are disjoint. Across the packet,

    r_0=sum_i min(a_i,v_i),
    r_h=dimension of the sum of the eight spaces E_tile.

Taking the maximum local dimension is always a lower bound. Two inexpensive
improvements use the beta64 geometry.

## Distinct one- and three-dimensional constraint spaces

Write H for the binary span of {1,2,4,64} in GF256. All input and output nodes
of tile alpha lie in alpha+H. The eight tile offsets 0,8,...,56 are in distinct
cosets modulo H.

Suppose v=a+d with d=1 or d=3. Put

    t=sum_(u in U_tile) u + sum_(y in Z_tile) y.

There are 2a+d nodes in this sum, an odd number. Hence t belongs to alpha+H,
so its possible values are disjoint between distinct tiles.

For d=1, E_tile is a line of functionals on degree-below-four polynomials.
Normalize its nonzero functional Lambda by Lambda(1)=1. This normalization
is possible: the constant polynomial one cannot be canceled at a+1 output
nodes by a rational function with only a poles.

The polynomial X+t can be canceled. To see this, let R_Z be the monic
polynomial with roots Z_tile. The top two coefficients of D_U(X)(X+t) and
R_Z(X) agree, so their difference has degree below a and supplies the needed
partial-fraction numerator. Therefore Lambda(X)=t.

The first two entries of the normalized functional are (1,t). Distinct tiles
cannot give the same line. Two one-dimensional tile constraint spaces thus
span at least two dimensions.

For d=3, the states satisfying the tile's constraints form a line of cubics.
Every nonzero polynomial on that line has degree exactly three. Otherwise
its rational numerator would have degree below a+3 and could not vanish at
all a+3 distinct selected output nodes.

Normalize that cubic to be monic. Comparing the coefficient of degree a+2
in its numerator with the monic R_Z gives its degree-two coefficient t.
Distinct tiles therefore give distinct surviving state lines. Their
annihilators are distinct three-dimensional spaces, whose sum has dimension
four.

A safe accumulated rank rule is consequently

    1 joined with 1 -> 2,
    3 joined with 3 -> 4,
    otherwise take the maximum.

The accumulated value three under this rule arises from an actual
three-dimensional tile space. Separate information from untouched coordinates
must not be folded into this state before applying the 3+3 rule.

## Two further rank bounds

Selected output zeros in untouched tiles constrain evaluations of the state
polynomial directly. The state map is MDS[64,4], so these constraints have rank

    b=min(4, sum_(i:a_i=0) v_i).

This recovers the exact state information used by the previous moment bound.

A global bound is useful when nearly every output is constrained. Put j=|U|.
The identities F_64^(-1)A=A T and the MDS property imply

    rank([F_U,A])=min(64,j+4).

Its kernel consists of states whose word A T h is supported on U. The
64-j complementary evaluation constraints have rank min(4,64-j), which
proves the formula.

Deleting 64-sum_i v_i rows can reduce rank by at most that number. Define

    g=sum_i (8-max(a_i,v_i)).

After subtracting the known input rank r_0, this gives

    r_h >= max(0,min(4,64-j)-g).

The final state-rank lower bound ell is the maximum of the accumulated tile
bound, b, and this global bound. It always lies between zero and four.

## Positive dynamic program

For one tile with counts (a,v), its multiplicity and weight-mask coefficient
are

    binom(8,a) binom(8,v) z^(8-v)(1-z)^v.

The input-rank factor is q^(-min(a,v)). Scale the packet computation by q^j,
so the tile's coefficient instead contains the integer factor

    q^(a-min(a,v))=q^max(a-v,0).

Retain the state (j,r,b,g), where r is the accumulated tile rank bound and
b and g are capped at four. On adding a tile, increase j by a, update r by
the stated join rule, add v to b only when a=0, and add 8-max(a,v) to g.
States reaching r=4 or b=4 have zero nonzero-state contribution and may be
discarded.

At the end, set

    ell=max(r,b,max(0,min(4,64-j)-g)).

Multiply that state's coefficient by

    (q^(4-ell)-1)/(q^4-1),

sum the coefficients at j, and divide by q^j binom(64,j). The result bounds
V_j(z). All surviving coefficients are nonnegative. The actual computation
reaches at most 744 states after any tile.

Keeping only the bound b recovers the previous bound exactly: it uses the
zero-coset moments on occupied tiles and the nonzero-state MDS moment on
untouched coordinates. The new ell is at least b for every mask, so the new
bound is pointwise no larger before any averaging.

## Consequences for the other state components

The improved V_j can also tighten the residual and column bounds in the
three-component operator. Let D=q^4-1. For any fixed nonzero old state h,
nonnegativity gives

    E_(U,X)[z^wt(FX+A h)] <= D V_j(z).

For a fixed target and scalar lambda, the coupled column selects a state
h=lambda(b_0+CX) depending on X. At every fixed (U,X), its nonzero-state term
is still one term in the sum over all nonzero h. Thus the same bound D V_j
also applies to the normalized I/L column moment and its zero-target version.

Each of L_j, gamma_j, and gamma_j^0 may therefore be capped by D times the
new upper bound. This argument preserves the dependence of the selected state
on X; it does not treat the selected state as independent or uniform.

## Bounded checks and numerical gains

The following entries are directed upper logarithms for the new moment at
z=1/256. The final column compares them with the earlier floating bound.

| Potential inputs j | New log2 moment upper bound | Approximate gain in bits |
| --- | ---: | ---: |
| 16 | -468.0723795027 | 255.94 |
| 24 | -491.5738849949 | 208.49 |
| 32 | -510.5455403925 | 163.54 |
| 40 | -511.6413951347 | 114.15 |
| 48 | -511.6414322611 | 68.06 |

The 192-bit enclosures have radii below 8e-55 before decimal shortening.
At the weaker mark z=1/16, the gains at j=32 and j=40 are approximately
36.72 and 18.98 bits. Gains at z=1/2 are small. These are local improvements;
their effect on a whole-code bound requires composing the revised operator.

The implementation and tests are
[theory/hd_randomized_outer_flat_moment.py](theory/hd_randomized_outer_flat_moment.py)
and [theory/hd_randomized_outer_flat_moment_tests.py](theory/hd_randomized_outer_flat_moment_tests.py).
Five tests pass. They check 648 one-tile rank cases, both distinct-space joins,
48 packet masks on the actual beta64 maps, and every input/output-zero mask
of a four-position GF8 toy. The toy also checks the exact positive removal of
the zero state. Directed tests verify the middle-occupancy improvement and
the exact full-input endpoint.

The provider `flat_moments(z)` returns floating base-two upper logarithms.
With `directed=True`, it returns Arb enclosures of the moments themselves.
The directed interface accepts a rational mark or an Arb enclosure, including
an exact algebraic mark such as `arb(2)**(-arb(3)/256)`. The caller selects
the Arb precision.

    C:/Python314/python.exe -B -m unittest hd_randomized_outer_flat_moment_tests -v
    C:/Python314/python.exe -B hd_randomized_outer_flat_moment.py

The next step is to substitute these local bounds into the existing regional
composition, including the D V_j caps where useful. Keep exact regional
coefficients for sparse row counts. The present result supplies the local
rank proof and positive moment computation, not a new global range claim.
