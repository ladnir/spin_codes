# Small polynomial coefficients and invertible state mixing

The immediate two-scan construction can use cheap characteristic-two
coefficients without restricting its generator to an embedded GF(256).
The same coefficients also give a candidate for a cheaper batched SPIN inner.
This note establishes local algebra, not a complete distance certificate.
No performance measurement accompanies these arithmetic counts.

## Preserve the consumer's field operation

Fix F = GF(2^m). A field-linear encoder E satisfies E(a*x) = a*E(x)
for every a in F. That identity preserves field-valued linear combinations
in a PCS and the correlation q = t + Delta*x in field VOLE.
An arbitrary binary-linear map inside each field element need not satisfy it.
For example, squaring is binary-linear and invertible, but
(a*x)^2 = a^2*x^2 rather than a*x^2.

The existing binary SPIN matrix remains field-linear when its XOR circuit
acts coordinatewise on characteristic-two field records. Its binary
coefficients explain why scalar extension does not increase its distance.
Replacing those coefficients by genuine elements of F preserves field
linearity and changes the code family. Arbitrary binary mixing of the bits
inside a record has a different interface and cannot be substituted silently.

Field-valued coefficients also change the receiver's alphabet. They do not
automatically preserve binary choices required by a bit-COT consumer.
Field VOLE and field-linear PCS are the direct targets; binary OT needs its
own compatibility argument.

The new transposed encoder must transpose the matrix over F. A coefficient c
remains c under that transpose. The test identity uses the field pairing
sum_i x_i*y_i, not the parity of a polynomial-basis bitwise dot product.
Binary matrix adjoints of multiplication by c are generally different in
this basis. A standard bit-parity adjoint test would therefore reject a
correct field transpose unless it also accounts for the dual basis.

## A byte names a polynomial, not a subfield element

Represent F as GF(2)[X]/(p_m). Bit i names the coefficient of X^i, and let
theta denote the residue class of X. The exact Rabin test in `arithmetic.py`
confirms the following sparse irreducible polynomials:

| m | p_m(X) | Low reduction constant | theta^256, represented in this basis |
|---:|---|---:|---:|
| 8 | X^8 + X^4 + X^3 + X + 1 | 0x1b | 0x2 |
| 32 | X^32 + X^7 + X^3 + X^2 + 1 | 0x8d | 0x8c010001 |
| 64 | X^64 + X^4 + X^3 + X + 1 | 0x1b | 0x11011 |
| 128 | X^128 + X^7 + X^2 + X + 1 | 0x87 | 0x4015 |

Define S = {sum_(i=0)^7 c_i*theta^i : (c_0,...,c_7) != 0}.
The 255 byte patterns therefore name distinct nonzero field elements.
They are not small integers interpreted by repeated field addition.
In characteristic two, that interpretation would collapse to zero and one.

For m > 8, S together with zero is not the GF(256) subfield: theta belongs
to S but theta^256 != theta. More simply, 0x80*0x02 = 0x100 lies outside
the byte set. Because theta belongs to S and has degree m, S generates F.
Thus these coefficient choices do not confine the complete generator to
GF(256). This observation removes a subfield obstruction; it does not prove
that a particular construction achieves the larger-field GV distance.

Bytewise AES-field multiplication is not this operation in the chosen
polynomial basis. The tests exhibit a failure of field-linearity if those
bytewise products are substituted without a compatible field representation.

## Cheap multiplication by a member of S

For a fixed byte c, Horner evaluation uses at most seven reduced left shifts
and eight conditional XORs of the input. Its shifts multiply by theta;
an overflow bit selects the low reduction constant. This is an exact
reference algorithm, not a claim about optimal instruction scheduling.

For m >= 32, let p_m = X^m + r and compute the unreduced carryless product
x*c = low + X^m*high. Since c has degree at most seven, high has degree
at most six. Each selected r has degree at most seven, so

    x*c mod p_m = low XOR carryless_product(high, r).

The correction has degree at most thirteen, below m, and needs no second
fold. This exposes a short CLMUL implementation option. The GF256 case
needs its own reduction or the native byte-field multiplication operation.
Input and output retain the declared polynomial basis; no representation
conversion is hidden in this arithmetic description.

## What the small-set proof buys

Let r_1,...,r_w be fixed linear forms over F whose span has rank h.
Let ell be another fixed form. Sample a_j independently from a set of
Q distinct field elements. Then

    Pr[ell + sum_j a_j*r_j = 0] <= Q^(-h).

Condition on all coefficients except those of h independent pivot forms.
At most one assignment to the remaining h coefficients solves the equation.
This is the characteristic-independent pivot argument used by the companion's
[rank-mask transfer](../../../../../gen-BAA/field-universality/expander-lift/CONVOLUTION_RANK_MASK_TRANSFER.md).
Our polynomial set supplies Q=255. Smaller subsets remain available.
Nonzero coefficients prevent a scalar scan from dropping a nonzero state
when its next input is zero.

The forms must be fixed before these coefficients are exposed. Earlier
zero equations may be quotiented out, but the current coefficients must
remain independent of that earlier information. Uniform coordinate
marginals are insufficient: (a,a) has joint target probability 1/Q,
although each coordinate separately has atoms of size 1/Q.

## Three shears for a two-symbol state

For a later batched inner, sample independent a,b,c from S at each step.
For incoming state (x,y), perform

    u = x + a*y
    v = y + b*u
    output_state = (u + c*v, v).

Every stage is an invertible field-linear shear. Their product M has
determinant one, so a zero input cannot erase a nonzero state if the
state recurrence is Q_next = M*Q + C*X. This update costs three small
scalar products and three additions. It is not uniform GL(2,F).

Fix a nonzero incoming vector and an arbitrary target (z,w). There are
at most Q coefficient triples yielding that target, hence

    Pr[M*(x,y) = (z,w)] <= 1/Q^2.

If y=0, then u=x is nonzero. The second coordinate determines b and
the first determines c; a is unused. There are at most Q triples.
If y!=0 and w=y, then b*u=0. Since b!=0, this forces u=0,
fixes a and c, and leaves b free: at most Q triples. No u!=0
branch can reach this same target. If w differs from both zero and y,
each a fixes u, and the target determines b and c: at most Q triples.
If w=0, then u=z and the target fixes a and b while leaving c free.
That case also has at most Q triples.

Two shears are weaker. For (x,y)=(1,1) and a=1, they output (0,1)
for every b. Their joint atom is 1/Q, not 1/Q^2. The third shear
randomizes this exceptional branch instead of treating two marginals
as independent.

The block bound also has a rank-aware form. Fix a 2-by-D matrix H of
state forms and a target matrix T before sampling M. If rank(H)=1,
Pr[M*H=T] <= 1/Q^2 by the vector bound. If rank(H)=2, multiplication
by H is injective on 2-by-2 matrices. Moreover,

    M = [[1+b*c, a+c+a*b*c], [b, 1+a*b]].

Since b!=0, M determines b, then a, then c. Consequently
Pr[M*H=T] <= Q^(-3) at rank two. Rank zero is deterministic.
These bounds concern a specified complete block target. They do not
bound every subset of zero output coordinates, nor do they make entries
of M independent. A batched-inner certificate needs the corresponding
quotient transition analysis before using this replacement for uniform GL.

## Checks and next step

Twelve exact tests cover irreducibility, all byte coefficients against an
independent long-division oracle, subfield nonclosure, scalar-linearity,
the field transpose pairing, invertibility, shear atoms, and the
marginal-dependence counterexample.
The vector checks enumerate every projective state over GF4 and GF16,
including all 15 nonzero coefficients over GF16.

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/lane_mixing -p test_arithmetic.py -v
```

First use the polynomial coefficients in the companion's already-developed
two-scan proof. Keep three-shear state mixing as a separate next candidate
if a batched implementation can reduce routing or arithmetic cost.
