# Exact local kernel for the XOR state refresh

The support-flat interface in [the refresh note](HD_XOR_STATE_REFRESH.md)
needs more information than a live-block mask. The local kernel is now
computed exactly for zero, one, or two potential input positions in one
current RS[16,4] block. This is a local result, not a distance certificate.
There is no 625-state operator, regional screen, or native change here.

## Local experiment

Let A be the current 16-by-4 normalized novel-basis emission matrix over
AES GF256, and let C=A^T. Fix an exact state-coordinate subset S and an
exact input-coordinate subset U with |U|<=2. Sample the coordinates Q_S
independently from F_q^*, and X_U independently from F_q. All other
coordinates are zero, and q=256. Compute

    Y=X+AQ,       W=Q+CX.

The local output is the joint support/weight enumerator

    E_(T,w) = #{(Q_S,X_U): support(W)=T, weight(Y)=w}.

Its denominator is D=(q-1)^|S| q^|U|. Thus the marked kernel entry is

    K_(U,S,T)(z) = sum_w E_(T,w) z^w / D.

The kernel precedes the XOR refresh. Fresh output diagonals restore uniform
nonzero amplitudes on each exact state support after that refresh, as proved
in the refresh note. This local experiment uses the uniform-field input law
from the Toeplitz comparison, not the forced-nonzero XOR-outer input law.

## Exact calculation with small ranks

Assume first that S is nonempty. Change variables to V=Y_U=X_U+A_(U,S)Q_S.
For each Q_S this is a translation of X_U. Therefore V is uniform on F_q^U
and independent of Q_S. Put R=I_4+C_U A_U. The maps become

    Y_outsideU=A_(outsideU,S)Q_S,
    W=R_(:,S)Q_S+C_U V.

The full RS code has minimum distance13. Consequently Y outside U has at
most three zero positions whenever Q is nonzero. This uses the full code's
distance; it makes no MDS claim about A restricted to selected columns.

To expand the output mark, take D_out outside U with |D_out|<=3 and any
D_in subset of U. These subsets specify output coordinates constrained to
zero. Fix also B, a set of feedback coordinates required to be zero.
For each J subset of S, force Q_J=0 and form the constraint matrix

    [ A_(D_out,S\J)       0                 ]
    [ R_(B,S\J)          C_(B,U\D_in)       ].

There are d_J=|S\J|+|U\D_in| remaining variables. If this actual matrix
has rank rho, it has q^(d_J-rho) solutions. Inclusion-exclusion over J
enforces every coordinate of Q_S to be nonzero. Write the resulting count
as G_(D_out,D_in)(B). The denominator remains D: deleting the forced-zero
V variables does not remove their original factors of q.

For exact feedback support T, use a second integer inclusion-exclusion:

    N_D(T)=sum_(L subset T) (-1)^|L| G_D(T^c union L).

All signed operations occur on exact integers. Summing N_D(T) over D with
|D_out|+|D_in|=k gives a nonnegative integer b_(T,k). The positive expansion is

    K_(U,S,T)(z)=sum_k b_(T,k) z^(16-k)(1-z)^k / D.

The identity z^weight(Y)=product_i[z+(1-z)1_(Y_i=0)] proves this formula.
All larger outside-zero subsets have zero contribution. There are at most
697,1152,1880 admissible subsets for |U|=0,1,2 respectively. Grouping equal
rowspaces and discarding spaces with no all-nonzero Q further reduces work.
The implementation also converts the result to ordinary nonnegative weight
counts E_(T,w), checks their signs, and checks that their total is D.

For S empty, Y=X and W=CX. Direct enumeration uses at most q^2 inputs and
avoids the nonzero-state distance argument entirely.

## Implementation and checks

The isolated files are `theory/hd_xor_state_local.py` and its test module.
`local_counts(U,S)` accepts S as a four-bit mask and returns exact integer
enumerators. `marked_kernel(data,z)` evaluates them at a rational z.
Neither routine replaces a source support by a fresh random support.

The production matrix comes from the existing feedback map's first block.
Tests reproduce every emission basis vector with the native four-point
butterfly order and verify A^T A=0. Every four-row restriction of the full
matrix is checked nonsingular. In contrast,

    A[rows(0,1), columns(0,2)] = [[1,0],[1,0]]

has rank one. The selected-column subsets are not assumed MDS. This also
gives a concrete count-lumping distinction: with U empty, the singleton
state support {0} emits weight16, whereas {2} emits weight14.

Exact batched elimination uses at most six variables. A conservative bound
on the absolute sum before feedback-support inversion is
697(q+1)^(|S|+|U|). Multiplying by16 still stays below2^63 at the supported
limits. Thus the signed int64 accumulation does not wrap; polynomial
conversion operates on smaller nonnegative exact counts.

Seven tests cover all 44 GF4 RS[4,2] support cases through two input
positions, independent rank elimination, the actual native basis, the
selected-column counterexample, and resource guards. They also compare an
actual GF256 one-input kernel with exhaustive sampling of its finite domain.
For full state support and U={0,5}, a separate enumeration of all 256^2
input vectors checks the exact W=0 fiber, where Q is uniquely fixed by CX.

## Bounded native-map examples

The following entries use z=1/2. Logarithms are display-only FLOAT values;
the underlying enumerators and marked moments are exact integers/rationals.

| U | State mask S | log2 total marked mass | log2 marked mass at W=0 |
|---|---:|---:|---:|
| empty | 15 | -15.91561134 | zero mass |
| {0} | 1 | -15.99437545 | -23 |
| {5} | 3 | -15.91726653 | zero mass |
| {0} | 15 | -15.90998679 | zero mass |
| {0,1} | 3 | -15.91164198 | -31.92289107 |
| {0,5} | 15 | -15.90998626 | -47.90624471 |
| {0,5} | 0 | -1.98875090 | -16 |

For S=15 and U={0,5}, all sixteen feedback supports occur. The exact
denominator is 255^4*256^2=277102632960000, and 65025 source pairs have W=0.
The local marked return is therefore much more informative than a uniform
arbitrary-state return cap. No global distance conclusion follows from it.

Run the local tests and the bounded seven-case example from the theory folder:

    python -B -m unittest hd_xor_state_local_tests
    python -B hd_xor_state_local.py --examples

The next proof gate should combine these exact local entries with the XOR
refresh's support law for a few controlled packet types. Incoming supports
with equal coordinate counts must be maximized over, not presumed equal.
Route averaging must keep its actual conditioning. A full 625-state sweep
or transfer of the existing certificate is not yet justified.
