# Weighted state density after an outer randomizer

A zero/nonzero state bound discards the flat part produced by each rank-one
refresh. A three-component bound retains it without treating births as
uniform. Combining this bound with the exact regional route certifies the
entire one-active-row class, including nonquiet words:

    Pr[there exists a nonzero one-active-row word of weight at most 838860]
        < 2^-391.9739840542.

The union over all 8192 choices of the active row is included. This is a
directed bound for the actual rank-one mixer. It does not cover every active
row count or constitute a complete distance certificate.

## State measures and the three components

Use the independently randomized GF65536 outer and shared beta64 inner from
[HD_RANDOMIZED_OUTER_PROOF.md](HD_RANDOMIZED_OUTER_PROOF.md). Put q=65536,
D=q^4-1, and let U_* be uniform on the nonzero four-coordinate states.
For a weight mark 0<z<1, weight every history by z to its emitted weight.
The resulting state measure is nonnegative; its mass need not be one.

Represent an upper bound on this measure as

    mu <= Z delta_0 + B U_* + nu,
    nu supported on nonzero states,       total_mass(nu) <= R.

The three stored numbers are Z, B, and R. B is a density cap relative to U_*.
The residual measure nu can have arbitrary shape. In particular, births from
zero enter R. This representation does not declare a weighted state uniform.

Write the actual rank-one refresh law on nonzero states as

    P = a I + b L + c J,
    a=(q^3-1)/D,       b=1/q-a,       c=1-1/q,
    epsilon=a+b=1/q.

Here L averages on a nonzero scalar line and J averages on all nonzero
states. Both I and L preserve pointwise density bounds. The J component
converts any incoming weighted live mass into a uniform nonzero measure,
before adding the current feedback.

## Row moments and column moments

Condition on j potential input positions U in the packet. Under the outer
uniform comparison, their amplitudes are independent uniform field elements,
including zero. After conditioning on packet counts, U is a uniform j-subset.
It is independent of the earlier weighted state measure.

Let a_i count the potential positions in tile i, and let u count occupied
tiles. A tile with a potential positions has a uniform MDS[8,a] feedthrough
image. Denote its weight moment by f_a(z), with f_0(z)=1, and put

    W_U = product_i f_(a_i)(z),
    ell_U = max(0,64-8u-3),
    L_U = W_U z^ell_U.

For a fixed nonzero old state h, the occupied-tile affine moments are at
most their zero-coset moments. This follows from the nonnegative Fourier
transform of the q-ary Hamming tilt. The untouched positions emit A h;
at most three of those values vanish. Hence L_U bounds the live row moment
uniformly over h.

A flat incoming component has a sharper average. Let v_m(z) be the exact
weight moment of A restricted to m coordinates when h is uniform nonzero.
The punctured state code is MDS when m>=4; for m<4, its evaluation map has
rank m and fibers of size q^(4-m). Its finite weight enumerator determines
v_m(z) in either case. Define

    V_U = W_U v_(64-8u)(z).

Then V_U bounds the averaged emission moment for the flat component.
The occupied-tile bound is uniform in h, so multiplying it by the untouched
state moment is valid. No independence between the two emissions is asserted.

Propagation of a density cap needs a column bound as well. Fix a successor
state b_0 and an I or L branch. Before refresh, the incoming state has the
form

    h=lambda (b_0+C X),       lambda!=0.

The emitted packet is consequently

    (F+lambda A C)X + lambda A b_0.

This expression retains the dependence of h on X. The fixed-map identities
F^2=I, F A=A T, and C A=0 imply

    C F^(-1) A=0,
    det(F+lambda A C)=det(F)!=0.

Therefore the restriction of F+lambda A C to j input columns has rank j.
Choose j output information coordinates. Their affine values are uniform
and independent, so their moment bounds the full output moment by t(z)^j,
where

    t(z) = (1+(q-1)z)/q.

The untouched-position bound still applies whenever h!=0. Thus every
normalized I or L column has the common upper bound

    gamma_U = min(t(z)^j, z^ell_U).

At successor zero, h=lambda C X. The zero input contributes q^-j to the
unrestricted moment and is excluded. Also C restricted to U has rank
r=min(j,4), so Pr[CX!=0]=1-q^-r. Consequently

    gamma_U^0 = min(t(z)^j-q^-j, z^ell_U (1-q^-r))

bounds the column at zero. For j=0 it is zero. These column bounds are not
obtained by substituting a dependent h into a bound proved for fixed h.

Average W_U, L_U, V_U, gamma_U, and gamma_U^0 over uniform j-subsets to
obtain W_j, L_j, V_j, gamma_j, and gamma_j^0. Tile counts suffice for all five
upper bounds. The helper uses a small dynamic program over (j,u).

## A positive local operator

For the zero-state component, use

    A_j = min(W_j, q^(-min(j,4))),
    B_j = W_j-q^-j.

A_j bounds the joint emission/zero-feedback moment. B_j bounds births into
the residual component by removing the zero input. Further nonzero quiet
inputs are retained conservatively in B_j.

For a fixed nonzero old state, an existing return bound gives, for j>0,

    kappa_j = min((q^3/D)L_j,
                  q^-r [epsilon+c(q^r-1)/D]),       r=min(j,4).

The first term uses the maximum atom of M h and the live emission moment.
For the second, the uniform syndrome hits M h with probability q^-r if
M h lies in its r-dimensional image. The nonuniform refresh contributes
at most epsilon to that event; its uniform part contributes (q^r-1)/D.
Define kappa_0=0 because an invertible mixer cannot kill a live state on
empty input.

With row and column order (Z,B,R), the following positive matrix propagates
the three-component upper bound for j>0:

    [ A_j                         0                    B_j          ]
    [ (c V_j+epsilon gamma_j^0)/D  c V_j+epsilon gamma_j 0            ]
    [ kappa_j                     c L_j                epsilon L_j  ].

At j=0, both live-to-zero entries are zero. All other entries have the same
form. Initial coefficients are (1,0,0), and terminal mass is bounded by
Z+B+R. The matrix is an upper operator; its rows need not sum to one.

To justify the flat-component row, first multiply the incoming uniform
measure by the emission tilt. Its mass is at most V_j. The J part, followed
by any feedback translation, has point mass at most c V_j/D at every target.
It therefore supplies flat coefficient c V_j, with a separate zero bound.
The I and L parts use the column caps just proved. Their flat coefficient
is at most epsilon gamma_j, and their zero mass is at most
epsilon gamma_j^0/D.

For a residual input of mass R, the weighted live mass is at most L_j R.
Its J part supplies flat coefficient c L_j R. The remaining nonzero part
has mass at most epsilon L_j R and remains residual. The separate zero
entry uses kappa_j. Possible overlap between these upper bounds only
increases the resulting measure.

## Empty gaps retain the refresh benefit

For an empty packet, L_0=z^61 and no new residual mass is born. After g
consecutive empty packets, the residual coefficient is at most

    R_after <= (z^61/q)^g R_before.

Meanwhile the refreshed mass remains in B. Its next return uses the
1/(q^4-1) density scale, with the emission and feedback dependence included
in the column norms. A zero/nonzero supremum would instead treat the entire
live measure as residual again at each packet.

This is a statement about unnormalized weighted mass. Conditioning on a
weighted history and renormalizing can increase a component's relative
share; no q^-g assertion about that conditional share is used.

## Exact regional routing at one active row

One active outer row supplies exactly one potential slot in each region.
Consequently, exactly one of the region's 128 packets uses the j=1 operator.
Every other packet uses j=0. If these two operators are P_0 and P_1, the
exact regional operator at fractional exponent one is

    R_1 = (1/128) sum_(i=0)^127 P_0^i P_1 P_0^(127-i).

The active packet is uniform; P_1 already averages its 64 possible positions.
The state carries between regions, so the complete bound uses R_1^256.
Outer domination and projective division give the failure upper bound

    8192 * beta/(q-1) * z^-838860 * e_0 R_1^256 1,
    beta=q^256/(q-1)^128.

This expression covers all 8192 active-row choices and all nonzero messages
in each row. It imposes no quietness condition. The artificial zero values
introduced by the outer envelope only enlarge the count.

At z=2^(-1/512), the 192-bit Arb enclosure of its base-two logarithm is

    [-391.973984054266951728929405218610708305145574594995312
      +/- 6.42e-52].

Only the j=0,1 operators are needed for this receipt. Their punctured state
moments have lengths 64 and 56 and positive integer MDS coefficients. The
helper computes the matrix coefficient of x in (P_0+x P_1)^128, divides by
128, and then raises the result to the 256th power. It uses no floating
acceptance test.

## What the bounded screens diagnose

The floating screen considers only active-row counts 1,16,128,512,1024,4096.
It includes fine tilts near one, rational fractional exponents, and several
Bernoulli markers. Its one- and two-packet marker bounds are:

| Active rows | One-packet log2 upper expression | Two-packet expression |
| --- | ---: | ---: |
| 1 | +482.06 | +478.57 |
| 16 | +2955.28 | +2917.58 |
| 128 | +18346.73 | +15704.56 |
| 512 | +46216.66 | +37509.82 |
| 1024 | +79827.67 | +74323.23 |
| 4096 | +209265.64 | +207475.71 |

These positive values are uninformative bounds, not evidence of distance
failure. Pairing steps helps but does not close any sampled marker bound.

An ablation identifies two different losses. Replacing the marker by the
exact one-row regional coefficient already gives about -204.44 bits with
the earlier two-state actual-mixer operator. Its ideal uniform-mixer control
gives about -207.40 bits. Keeping the flat density improves the actual-mixer
result to the directed -391.97-bit bound above. Thus regional relaxation is
the primary obstruction in this sparse example; density retention supplies
another approximately 187.5 bits of margin.

The next useful step is an exact regional coefficient calculation for active
row counts 2 through 16, combined with sharper occupied-tile state moments.
The present V_j bound gains information only from untouched
tiles; it loses that benefit when most tiles are occupied. A compact joint
rank calculation for [F_U,A] on selected output zeros can improve this part.
Increasing refresh strength alone does not address either of these losses.

## Reproduction

[theory/hd_randomized_outer_state_profile.py](theory/hd_randomized_outer_state_profile.py)
contains the positive operator, bounded floating screens, and directed
one-row receipt. Its five [tests](theory/hd_randomized_outer_state_profile_tests.py)
check the column bounds over a small exact field, perturbed-feedthrough
invertibility, punctured state enumerators, empty-gap contraction, and the
directed endpoint. All pass. The one- and two-packet exact regional
calculations agree at fractional exponent one.

    python -B -m unittest hd_randomized_outer_state_profile_tests -v
    python -B hd_randomized_outer_state_profile.py --receipt
    python -B hd_randomized_outer_state_profile.py --steps 1
    python -B hd_randomized_outer_state_profile.py --steps 2

No native implementation, benchmark, remote run, or prior proof file was
changed for this note.
