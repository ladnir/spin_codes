# FFT circuits for the sixteen-state maps

The certified inner applies eight related maps to its state in each packet.
Previously, each map used a separate four-term weighted sum. Choosing the
LRS representatives on an additive coset lets these sums share FFT work.
This changes the sampled code, but preserves all moments used by the
existing certificate.

Keep the fields, local maps, and tower basis from
[the tower-basis construction](HD_TOWER_BASIS_INNER.md). Write B=GF16,
F=GF65536, and beta=(1,79,256,20224). Let h_0,...,h_15 be the ordered
elements of B used by the local FFT. Set a_i=h_(8+i), for 0<=i<8.
These representatives are nonzero, and their norms from F to B are
a_i^4. Those eight norms are distinct because fourth powering permutes B.
Thus the original LRS subspace-rank theorem applies to these representatives.

Let n_0,...,n_3 be the first four normalized additive-basis polynomials
defined in the tower-basis note. For 0<=g,r<4, define the new blocks by

    R'_i[4g+r,l] = a_i^r n_g(a_i)^4 beta_l^(16^r).

For each i, set C_i=R'_i C0 and A_i=A0 (R'_i)^T. The packet feedthrough
F0, randomized recurrence, route, and outer distribution remain unchanged.
In particular, retain the independent rank-one state mixer and independent
diagonal-Cauchy feedback randomizer in each packet.

## Certificate transfer

Write n_g(X)=sum_(h=0)^g c_(g,h) X^h. Its leading coefficient is nonzero.
Since beta_l^(16^4)=beta_l, the new row is the following combination of
rows of the original LRS block with the new representatives:

    R'_i[4g+r,l] = sum_(h=0)^g c_(g,h)^4 R_i[4h+r,l].

Consequently R'_i=S R_i for a single invertible B-matrix S, independent
of i. For each residue r, S is triangular with nonzero diagonal. Applying
S to the concatenated restricted blocks preserves their F-rank.
The transpose change on the emission map likewise preserves every rank
of the joint input/output constraint matrix. The local residual subspaces
are unchanged. Therefore r0, rho, eta, the live-output floor, and the exact
local moments are identical to those in the existing proof.

The feedback randomizer is defined in the new coordinates; no claim of
distributional conjugacy is needed. Its point-density bound holds for every
fixed nonzero input vector in those coordinates. The state transvection
has the same distributional property as before. Hence the complete paired
and global certificate transfers, including the diagonal-Cauchy loss:
at K=2^20, rate 1/2 and 40% F-symbol distance, the whole-code setup-failure
probability is below 2^(-48.3283), under the stated independent setup model.

## Native circuit and checks

Group state coordinates by g, packing the four residues r in one SIMD
vector. Evaluate the degree-below-four polynomial in the fourth-powered
additive basis on two four-point cosets. Butterfly coefficients are the
fourth powers of the existing GF16 FFT coefficients. At point a_i, scale
the four residues by (1,a_i,a_i^2,a_i^3), then apply the transposed tower
Moore circuit. Feedback reverses this circuit and uses the transposed FFT.
Both directions retain all sixteen state symbols.

The exact Python builder and the independent native dense-map builder give
the F,A,C digest 17847517350586406002. Tests check residual-mask ranks,
nilpotence, dense matrix agreement, all packet basis directions, and full
encodings. These checks guard implementation errors; the row-change
argument establishes the universal certificate transfer.
