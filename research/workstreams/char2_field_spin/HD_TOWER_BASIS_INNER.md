# Factoring the inner maps through two quadratic extensions

The sixteen-state certificate uses the GF16-subspace ranks of its local
maps, not their particular polynomial bases. We can choose those bases
to simplify the native circuit without changing any local moment.
This note defines that refinement and explains why the certificate transfers.

Keep B=GF16 and F=GF65536. In the implementation, F is the AES GF256
extension with v^2+v+0x20=0; v is represented by 256. Put u=v^16+v,
represented by 79. Then u belongs to GF256, u^16=u+1, and v^256=v+1.
Thus (1,u,v,uv) is a B-basis of F: (1,u) is a basis of GF256 over B,
and (1,v) is a basis of F over GF256.

Use this basis in the sixteen-row LRS generator. For one group of four
rows, the Moore map evaluates the bilinear expression

    x0 + u_j*x1 + v_j*x2 + u_j*v_j*x3,

where (u_j,v_j) ranges over
(u,v), (u+1,v+u), (u,v+1), and (u+1,v+u+1).
It has the following circuit:

    a=x0+u*x1; b=x2+u*x3; c=a+x1; d=b+x3;
    h0=a+v*b; h1=c+(v+u)*d; h2=h0+b; h3=h1+d.

The transpose circuit reverses these linear stages. On four packed records,
each circuit requires one byte-field GFNI product and two GFNI products
for the full F multiplication. This replaces the general four-by-four map.
The LRS coefficient a_i^j remains a byte-subfield scalar; Frobenius rows
repeat after four. The native implementation packs four consecutive rows
and retains all sixteen state coordinates.

## Local polynomial basis

Let e0,e1,e2 be the first three elements of the original ordered B-basis.
Let V_j be their first-j-dimensional binary span. Define

    n_d(X) = product_(j: bit_j(d)=1)
                 product_(a in V_j)(X-a) / product_(a in V_j)(e_j-a).

The polynomials n_0,...,n_7 form a graded basis: n_d has degree d and
nonzero leading coefficient. Let E0 and E1 evaluate this basis on H and
delta+H, respectively. Set F0=E1 E0^-1, let A0 be E1's first four columns,
and let C0 be E0^-1's last four rows. This is the same local feedthrough
as before, with a different basis for emission and feedback.

Every polynomial of degree below four has zero high-four coordinates in
this basis. Translation by delta preserves that subspace. Hence C0 A0=0
and C0 F0^-1 A0=0 still hold. Compared with monomials, A0 undergoes an
invertible right B-basis change; C0 undergoes an invertible left B-basis
change on the degree-four-through-seven quotient.

For each potential support and output-zero mask, the residual local
subspaces therefore have exactly the same B-dimensions as before. They
remain B-defined. The LRS rank theorem applies to every such subspace and
to every B-basis of F, including (1,u,v,uv). Consequently r0, rho, eta,
the live-output floor, and all local moments in the existing certificate
are unchanged. Its paired and global witnesses transfer without loss.

Combining this basis choice with the
[diagonal-Cauchy feedback](HD_DIAGONAL_CAUCHY_REFRESH.md) therefore retains
the complete 40%-distance guarantee with more than 48.32 bits of margin.
This is a different sampled code ensemble, not a same-output optimization
of the earlier monomial-basis construction.

The exact Python map builder and native setup independently construct the
new matrices. They agree on the full F,A,C digest. Tests also check the
Moore circuit, nilpotence, and concrete output-mask ranks. The rank argument
above, rather than these finite tests, establishes the universal transfer.
