# Coupled scan construction

The candidate replaces the expensive inner maps by a small coupled state.
Its output can vanish without erasing the state, while a nonzero input always
activates an empty state. These properties address two different failure
mechanisms; neither alone establishes minimum distance.

## Field and target

Let E be the AES field GF(256), and use the existing tower representation

    F = E[T] / (T^16 + T^5 + T^3 + 0x32),
    theta = T mod (T^16 + T^5 + T^3 + 0x32).

The [prior construction](../char2_field_spin/CHEAP_LANE_MIXING_20261016.md)
specifies this representation and its irreducibility check. A record contains
the sixteen E coefficients of an F element. This is not an implicit conversion
from an application's arbitrary polynomial basis.

Fix K=2^20, N=2K, and g=K/128=8192 outer groups. Distance counts nonzero F
records, not bits. Widths w=2,4,8 define the initial comparison family.
Width two has since been [rejected](WIDTH_TWO_OBSTRUCTION.md) for 35% under
this route; its useful local algebra is retained below. All three widths
divide both N and the regional length g.

## Outer and one structured route

Retain the outer of the [certified baseline](../char2_field_spin/HD_SHIFTED_CLOSED_CERTIFICATE.md).
Split each group's 128 inputs into four rows of 32. Encode each row with
the same fixed RS[64,32] map over E, extended linearly to F. At each of the
64 output columns, apply an independently sampled invertible 4 by 4 Toeplitz
matrix over E to its four records. The distribution is uniform over
invertible Toeplitz matrices, not over all invertible matrices.

Write v[i,r] for group i's output coordinate r, where r=4j+a identifies
column j and component a. Coordinate r goes to region r, giving 256 regions
of g records each.

The following convention makes the route explicit. In region r, set e=r mod 2.
Its g/2 source pairs are

    ((2p+e) mod g, (2p+e+1) mod g),  0 <= p < g/2.

Sample an independent uniform permutation pi_r of the source pairs. For
each destination pair h, sample an independent fair swap bit sigma[r,h].
Read the two source indices from pair pi_r(h), swap their order if sigma is
one, and place their records v[i,r] in positions rg+2h and rg+2h+1.
These choices are independent across regions and of all outer and inner coins.

Concatenate the regions and divide the sequence into consecutive batches X_t
of w records. The definition uses one global permutation with internal
regional structure. It does not add an independent global interleaver.

## Inner setup and encoding

Let R rotate the w state coordinates by one:

    (R Q)[j] = Q[(j-1) mod w].

For every step t and coordinate j, independently sample c[t,j] from E and set

    a[t,j] = theta + c[t,j],
    D_t = diag(a[t,0], ..., a[t,w-1]).

All setup coins are sampled once and fixed for every message. The inner
coefficients are independent of the outer and route. Begin at Q_0=0 and emit

    Y_t = X_t + D_t Q_t,
    Q_(t+1) = R Q_t + Y_t.

Output all N records. State persists across regions. The terminal state is
discarded; no return-to-zero constraint is imposed on it. R is a fixed local
state shuffle, not a second global data permutation.

### Antipodal wiring used by the partial certificates

The width-four partial certificates use cycle order (0,2,1,3), rather than
the natural order (0,1,2,3). If this order is c, define
(R Q)[c[(i+1) mod 4]]=Q[c[i]]. The physical input pairs remain (0,1) and
(2,3); only the fixed local state shuffle changes. The determinant and
invertibility arguments below hold for either four-cycle.

The corresponding width-eight alternative has cycle order (0,2,4,6,1,3,5,7).
It is specified for comparison, not certified by the width-four calculations.

## Algebra valid for every setup

### Rate and inverse

From the output stream, reconstruct Q_0=0 and

    X_t = Y_t + D_t Q_t,
    Q_(t+1) = R Q_t + Y_t.

Thus the inner is an invertible F-linear map on F^N. The outer is injective,
so the full encoder is F-linear of rate 1/2.

### A live state cannot disappear on an empty input

Every a[t,j] is nonzero because theta is outside E. For a cyclic rotation,
in characteristic two,

    det(R+D_t) = 1 + product_j(theta+c[t,j]).

For 1 <= w < 16, the right side is a nonzero polynomial in theta of degree w.
It cannot vanish because theta has minimal polynomial of degree 16 over E.
Consequently R+D_t is invertible.

The recurrence therefore has four whole-vector properties:

- X_t=0 and Q_t!=0 imply Q_(t+1)!=0.
- Y_t=0 and Q_t!=0 imply Q_(t+1)=R Q_t!=0.
- Q_t=0 implies Y_t=Q_(t+1)=X_t.
- Q_(t+1)=0 with Q_t!=0 implies Y_t=R Q_t!=0 and
  X_t=(R+D_t)Q_t!=0.

Individual coordinates may still vanish. The minimal-polynomial proof does
not cover w=16; do not silently extend the width range.

### Width two on empty stretches

For w=2, R+D_t has nonzero entries and determinant. A nonzero state therefore
satisfies wt(Q_t)+wt(Q_(t+1)) >= 3 when X_t=0. Since D_t preserves support,
two consecutive empty-input steps emit at least three nonzero records out
of four if their initial state is nonzero.

There is also a sharper conditional statement for a fixed actual state
before exposing the current D_t. A state of weight one becomes weight two.
For weight two, the two cancellation targets are Q[1]/Q[0] and Q[0]/Q[1].
Both cannot lie in theta+E: their product is one, contradicting the degree-two
minimal-polynomial argument. At most one coordinate can therefore cancel,
with probability at most 1/256. The outgoing state has weight two otherwise.

This conditional probability does not cover a message chosen after setup.
The [rank argument](PROOF.md) addresses that quantifier order separately.

## Implementation hypothesis

Multiplication by theta+c is theta times the record XOR c times the record.
The former is a fixed byte shift and reduction; the latter uses byte-field
GFNI. It does not require a general 128-bit field multiplication.

The recurrence uses one such product and two record XORs per output record,
plus the state rotation. A later kernel must fuse regional pair reads into
the inner and write each output once, without an extra routed buffer.

This applies lane mixing four times as often as the previous 16-state,
64-output inner, but removes its A, C, and H circuits. Whether that tradeoff
wins is unknown. Width two has simple local algebra but fails the global
target. Widths four and eight may offer better SIMD scheduling. No timing
is assigned to these candidates.
