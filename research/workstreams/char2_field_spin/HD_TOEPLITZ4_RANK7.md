# Seven-product Toeplitz4: exact circuits, costly fixed transforms

An explicit seven-product circuit exists, but it does not yet justify native
work. Preserving the original Toeplitz map costs eight additional fixed byte
maps. A fixed-congruence variant reduces that overhead to six. The existing
Compact7 circuit already stores seven keys while retaining nine products.
Thus the new rank-seven circuit has no further key-stream saving over it.

No native source, provider, benchmark, or certificate changed.

## Field and setup

Work in the AES field GF256. It contains GF4 and GF16, but not GF8:
three does not divide eight. Seven small byte labels may be evaluation
points, but they are not an embedded GF8.

Fix the AES bytes

    omega=0xbc, omega^2=0xbd, gamma=0x5c.

They satisfy omega^2+omega=1 and gamma^2+gamma=omega. Use the seven
projective points (0,1,omega,omega^2,infinity,gamma,gamma+1).
For a finite point v, its evaluation row is (1,v,...,v^6); the infinity
row selects coefficient six. Let V be the resulting 7-by-7 matrix.
Given the original seven Toeplitz diagonals t, prepare lambda=V^(-T)t
in setup. Sample and reject the original Toeplitz matrix before this
deterministic key conversion.

Let E contain the corresponding evaluation rows through degree three.
Its infinity row selects coefficient three. Then

    H(t)[i,j]=t[i+j],    H(t)=E^T diag(lambda) E,
    T(t)x=H(t) reverse(x).

The seven variable products are lambda_i*(Ex)_i. Both E and E^T remain
runtime operations. Setup interpolation cannot remove the output transform.

## Exact original-map schedule

For monomial coefficients c, use novel coordinates

    u3=c3; u2=c2+c3; u1=c1+u2; u0=c0.

All plus signs here mean field addition, implemented by XOR. Evaluation is

    l=u1+u3; a=u0+u2
    y0=u0; y1=u0+u1
    y2=a+omega*l; y3=y2+l; y4=u3
    d=u1+omega*u3; b=u0+omega*u2
    y5=b+gamma*d; y6=y5+d.

This costs four fixed maps and nine XORs. Its explicit algebraic transpose
costs four fixed maps and twelve XORs. The two basis-conversion XORs and
their two transposed XORs make the complete original-map circuit
seven variable products, eight fixed maps, and 25 XORs. Input reversal is
register naming, not a shuffle of bytes within a record.

## A cheaper fixed congruence

Let C invert evaluation on the four GF4 points. Replacing the mixer by
C^T H(t) C makes the first four rows of E C the identity. Its evaluation is

    a=x0+x1; b=x2+x3; e=a+b
    d=omega*e+a; u=x0+x3; c=omega^2*u+x2
    v=gamma*d+c; w=v+d
    outputs=(x0,x1,x2,x3,e,v,w).

This costs three fixed maps and eight XORs. The transpose costs three
fixed maps and eleven XORs. Both schedules are explicit in the helper.

This is a new mixer map, not the old Toeplitz map in unchanged coordinates.
It preserves invertibility and the atom cap: for fixed nonzero x,
C x is nonzero, and C^T is a bijection on outputs. The uniform-measure
outer proof therefore transfers algebraically, with a new construction
binding. This argument does not assert that C preserves Hamming weight.

## Why GF4 plus infinity does not give seven GF4-only products

The Hankel space has dimension seven. In a seven-product field-bilinear
decomposition, the seven rank-one coefficient matrices must span that
space. Each must itself be Hankel. Over GF4, the projective rank-one
Hankel matrices are exactly the four evaluation directions and infinity.
Five directions cannot span a seven-dimensional space. This rules out
rank seven with both fixed transforms over GF4, even if setup key forms
use GF256 coefficients. It is not a lower bound on arbitrary bit-level
GFNI circuits.

Leaving the quadratic place un-diagonalized gives a concrete eight-product
alternative with GF4-only runtime transforms. Retain outputs
(x0,x1,x2,x3,e,c,d) from the congruent schedule, omitting the gamma map.
For the last two keys, put

    k=lambda5+lambda6
    ell=gamma*lambda5+(gamma+1)*lambda6.

Their bilinear block is

    [[k,ell],[ell,ell+omega*k]].

Apply it using three prepared products:

    shared=ell*(c+d)
    out_c=shared+(k+ell)*c
    out_d=shared+(omega*k)*d.

Together with five scalar products, this costs eight variable products,
four fixed maps, and 18 XORs. It realizes the same congruent mixer.
The quadratic polynomial splits over GF16; diagonalizing that place recovers
the seven-product circuit and its extra fixed gamma maps.

## Native cost and recommendation

A fixed multiplication by omega, omega^2, or gamma is one byte-linear
GFNI affine map. It is not a free XOR. These fixed constants need no
per-group key expansion, but their matrices consume registers or loads.
The counts below exclude loads of input/output vectors and setup work.
Each variable key still needs the existing four-byte expansion across
four 128-bit records.

| Circuit | Variable GFNI products | Fixed affine maps | XORs | Key expansions | RS + mixer GFNI sites |
|---|---:|---:|---:|---:|---:|
| Existing Toeplitz4 | 9 | 0 | 15 | 9 | 1092 |
| Existing Compact7 | 9 | 0 | 18 | 7 | 1092 |
| Exact original rank7 | 7 | 8 | 25 | 7 | 1476 |
| Congruent rank7 | 7 | 6 | 19 | 7 | 1348 |
| Congruent quadratic rank8 | 8 | 4 | 18 | 8 | 1284 |

The final column adds the existing 516 RS sites to 64 mixer calls per
four-group batch. Seven keys occupy 3.5 MiB at K20, as in Compact7;
eight keys occupy 4 MiB. GFNI multiplication and affine instructions may
have different port costs, so these totals are not runtime predictions.
However, neither tested circuit removes a pass or beats Compact7 on key
traffic. Six affine maps replace only two variable multiplications in
the best rank-seven variant.

Park this candidate. Reconsider native work only with a concrete way to
remove several fixed maps or verified machine-level evidence that their
port cost is sufficiently lower. No broader circuit search is proposed.

## Exact validation

The isolated helper is theory/hd_toeplitz4_rank7.py; its seven tests are in
hd_toeplitz4_rank7_tests.py. They passed in 0.300 seconds. For each of the
three circuits, all 56 key-bit by 32 input-bit basis pairs agree with the
specified dense map. Both transforms and their transposes are also checked
independently. Tests reconstruct every key direction, count the explicit
operation sites, enumerate the five GF4 rank-one Hankel directions, and
verify that the congruent map is not accidentally labelled as the old map.

Reproduce from the theory directory:

    C:/Python314/python.exe -m unittest hd_toeplitz4_rank7_tests -v
