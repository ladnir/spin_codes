# The two halves of the state-32 return map

The state-32 return bound currently ignores where output-zero constraints
lie inside a packet. Its fixed map has more structure than an arbitrary
rank-32 matrix. This observation may support a sharper count; it does not
give a distance certificate or change the construction.

Use the AES representation of E=GF256 and its fixed embedding into
F=GF(2^128). Let U be the additive subspace with byte labels 0..31.
The state nodes are V=128+U, and the output nodes are U and 32+U.
Let A be the 64-by-32 Lagrange evaluation matrix from V to these outputs.
Put P=A A^T. All additions below are field additions, hence XOR on labels.

Define L_U(X)=product_(u in U)(X+u), L_V(X)=L_U(X+128),
and lambda=L_U'. The subspace polynomial is additive, and its derivative
lambda is a nonzero constant. For an output node x, put
c_x=L_V(x)/lambda. Lagrange interpolation gives

    A[x,j] = c_x / (x+128+j),       j in U.

The logarithmic derivative gives sum_(v in V) 1/(x+v)=1/c_x.
Thus P[x,x]=1. For distinct x and y, partial fractions give

    P[x,y] = (c_x+c_y)/(x+y).

The value c_x is constant on each U-coset. Therefore, when coordinates
are ordered by their two output halves,

    P = [ I  B   ],     B[x,y] = kappa/(x+y),
        [ B^T I  ]

where x is in U, y is in 32+U, and kappa=L_U(32)/lambda.
The two coordinate sets are disjoint, so B is a Cauchy matrix multiplied
by a nonzero scalar. Every square minor of B is nonsingular.
The identity A^T A=0 also gives P^2=0, hence B B^T=B^T B=I.

Consequently, the map relevant to a return to zero has the form

    I+gamma P = [ (1+gamma)I   gamma B     ].
                [ gamma B^T   (1+gamma)I  ]

The within-half off-diagonal entries vanish. It would therefore be wrong
to treat all off-diagonal submatrices as MDS. A future return calculation
could instead track active and zero coordinates in each half separately.
Such a calculation must preserve the actual paired-route distribution;
it cannot assume independent uniform support in the halves.

This identity is valid for the state-32 maps, not the state-30 certificate.
It does not remove the principal-nonsingularity hypothesis in the
[random-global-gamma construction](HD_FIELD128_RANDOM_GAMMA32.md), nor its
once-per-setup failure term.

An exact GF256 check on Peach rebuilt A from the existing generic Lagrange
routine, multiplied A A^T, and checked all 4096 entries on 2026-10-10.
Both diagonal blocks were identity; every cross entry satisfied
(x+y)B[x,y]=200; and B B^T=I. This finite check supports the algebra above,
not a numerical distance claim. No encoder change or timing was made.
