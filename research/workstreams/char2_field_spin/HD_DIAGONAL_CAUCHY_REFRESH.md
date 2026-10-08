# A cheaper feedback refresh with a complete distance bound

The [sixteen-state construction](HD_SUMRANK_FEEDBACK_CERTIFICATE.md) needs
nearly uniform feedback directions. Exact uniformity is not necessary.
Two random diagonal scalings around a fixed Cauchy map suffice, with a
loss that the existing distance certificate can absorb.

Keep F=GF65536, s=16, the same outer, route, fixed inner maps, and rank-one
state mixer. Fix an invertible s-by-s matrix H over the GF256 subfield of F
whose every square minor is nonzero. For each packet independently sample
diagonal matrices D1 and D2, with all 2s entries independent and uniform
in F excluding zero. Replace the feedback map T by D2 H D1. These coins
are independent of all other setup coins. Setup is fixed for every message.

The native H interpolates a polynomial of degree below 16 from the AES-field
points 0,...,15 and evaluates it at 16,...,31. Its entries are

    H[y,x] = P(y) / ((y-x) P'(x)),   P(X)=product_(a=0,...,15) (X-a).

All row and column factors are nonzero. The Cauchy determinant formula
therefore proves that every square minor is nonzero. The native circuit
uses the existing inverse and forward additive FFTs, not dense products.

## Pointwise distribution bound

Fix a nonzero x in F^s, with support S of size w. The vector D1 x has
independent uniform nonzero coordinates on S. Fix a nonzero target y,
and let Z be its zero coordinates, of size z. If z>=w, the event
(H D1 x)_Z=0 is impossible: any w of those rows restricted to S form an
invertible matrix. If z<w, choose z pivot columns in S. Each choice of
the other w-z nonzero inputs determines at most one choice of the pivots.
Thus

    Pr[(H D1 x)_Z=0] <= (q-1)^(-z).

Conditioned on a vector with exactly the support of y, D2 maps it to y
with probability (q-1)^(-(s-z)). Other supports cannot reach y. Hence

    Pr[D2 H D1 x=y] <= (q-1)^(-s) = b/(q^s-1),
    b := (q^s-1)/(q-1)^s.

The image cannot be zero because all three factors are invertible.
This is a pointwise bound for every fixed x, not an independence claim
between messages. The whole encoder remains F-linear and injective.

## Transferring the complete certificate

Let P_U be the certified two-component weighted-state operator for one
packet's potential input set U. On nonzero feedback, the new randomizer
has pointwise density at most b times the ideal nonzero-uniform law.
On zero feedback, both randomizers are zero. The new operator is therefore
entrywise bounded by b P_U. This argument also covers the return to zero:
the preceding rank-one state update is invertible, so its cancellation
requires a particular nonzero feedback value.

There are n_p=32768 packets. The paired-profile and packetwise witnesses
therefore increase by at most alpha*n_p*log2(b), where 0<alpha<=1 is the
witness exponent. Bounding this by n_p*log2(b) gives a uniform loss for
all receipts. Coefficient extraction, route averaging, and union bounds
preserve the positive majorization. The dense inverse-output bound remains
valid because the new inner is invertible and F-linear.

Consequently the existing complete certificate transfers after adding
n_p*log2(b) to its total logarithmic failure bound. This is a conservative
transfer; it does not require a new numerical search or pretend that the
diagonal-Cauchy distribution is exactly uniform.

The [directed transfer helper](theory/hd_diagonal_cauchy_refresh.py) evaluates
the loss, and can replay the complete underlying certificate before applying
it. Its quick mode uses only the conservative, already established bound
2^-59.87. The native encoder is a separate implementation obligation.
