# A diagonal-refresh packet with five-product maps

The diagonal-refresh candidate can use the same one-stage regional routing as the GRS[32,22] design without two extra dense changes of state basis per packet. Its two fixed maps admit direct four-point XOR-circulant circuits. Each map costs five GFNI512 products on sixteen 128-bit records.

Including input diagonals, state refresh, and a conservative emission scale gives \(169/88\approx1.92045\) GFNI512 processed-byte equivalents per source record. This is an arithmetic estimate, not an implementation or runtime result. The ensemble differs from the once-per-region construction.

## Fixed geometry and sampled maps

Use AES GF256 and its actual subfields

\[
 H=\mathbb F_{16},\qquad
 H_4=\mathbb F_4=\{0,1,188,189\}.
\]

The integers denote AES-field bytes, not polynomial-basis GF16 indices. Retain the four sixteen-point phases \(H+\delta\), with \(\delta\in\{0,2,6,4\}\). Their union has 64 distinct labels. Choose the four residue nodes

\[
 v_j\in16+H_4=\{16,17,172,173\},
\]

which are disjoint from all 64 labels. These nodes remain fixed through every packet and phase. Their monic vanishing polynomial is

\[
 Q(X)=X^4+X+c,\qquad c=16^4+16,\qquad Q'(X)=1.
\]

For a packet label \(u\), define the emission and injection matrices by

\[
 A_{u,j}=(u+v_j)^3+1,
 \qquad C_{j,u}=(u+v_j)^3.
\]

The packet has input \(X\in\mathbb F_{256}^{16}\) and state \(S\in\mathbb F_{256}^4\). Sample four independent nonzero diagonal entries, independently in every packet and independently of the current history. With \(D_s\) their diagonal matrix and fixed \(\eta=2\), apply

\[
 R=D_sS,\qquad
 Y=X+\eta AR,\qquad S'=R+CX.
\]

The input \(X\) already includes the independent nonzero coordinate diagonals of the outer construction. The outer remains fixed GRS[32,22], with independent regional permutations and one routing stage. Any once-per-region orientation can remain as a boundary change of state coordinates; the fresh diagonal is sampled in the residue coordinates defined here.

## Superregular emission and the old cubic pairing

Since each \(v_j\) is a root of \(Q\),

\[
 A_{u,j}=\frac{Q(u)}{u+v_j}.
\]

Every \(Q(u)\) is nonzero. Thus \(A\) is a row-scaled Cauchy matrix. Any square minor of size at most four is nonzero: the Cauchy determinant has distinct row labels, distinct residue nodes, and no zero cross-difference. This applies to the full 64-by-4 emission matrix, not only to each sixteen-row phase.

For comparison with the dual-MDS map, write

\[
 V_u=(1,u,u^2,u^3),\qquad
 M_{\cdot,j}=(v_j^3+1,v_j^2,v_j,1)^T,
\]

and let \(J\) reverse four coordinates. The columns of \(M\) are the coefficient vectors of \(Q(X)/(X+v_j)\). Evaluation at the four nodes gives the identity matrix because \(Q'(v_j)=1\). Consequently,

\[
 M^{-1}=(v_j^k)_{j,k=0}^3,
 \qquad A=V M,
 \qquad C=M^{-1}J V^T.
\]

The state-coordinate change \(B=M^{-1}J\) therefore transforms the old cubic pair into this pair. In particular,

\[
 A_uC_w=(u+w)^3.
\]

Also, \(CA=0\). One proof uses the displayed basis change and the vanishing power sums over an affine GF16 phase. A direct proof sums products of two degree-three polynomials over that phase; every power has degree at most six and its sum is zero.

For any allowed-output set \(J_0\) inside a packet,

\[
 I+\eta C_{J_0}A_{J_0}
\]

is similar to the old eta-scaled cubic update. The all-subsets invertibility theorem in [DUAL_MDS_OBSTRUCTIONS.md](DUAL_MDS_OBSTRUCTIONS.md) therefore carries over. The fact that \(\eta=2\notin H\) is essential to that theorem. There is no zero-overlap return.

These identities concern fixed maps. Sampling a new \(D_s\) every packet changes the encoder distribution and does not preserve the old constant-state run law.

## Direct circuit without state-basis conversions

Order \(H_4\) as \(0,1,\omega,1+\omega\), where \(\omega=188\). Partition each sixteen-point phase into four additive \(H_4\) cosets. For one coset, the injection block is the symmetric XOR-circulant matrix

\[
 G_{i,j}=g_{i\mathbin\oplus j},
 \qquad g_j=(a+h_j)^3,
\]

where \(a\) is its fixed offset from the residue coset. The emission block differs from \(G\) by the all-ones matrix.

For a four-vector \(x\), define its superset zeta transform

\[
 z=(x_0+x_1+x_2+x_3,\ x_1+x_3,\ x_2+x_3,\ x_3).
\]

This transform is its own inverse. Let \(\gamma\) be the same transform of \(g\). The cubic power sum over \(H_4\) gives \(\gamma_0=1\). XOR convolution becomes multiplication in
\(\mathbb F_{256}[e_1,e_2]/(e_1^2,e_2^2)\), so compute

\[
 \begin{aligned}
 w_0&=z_0,\\
 w_1&=z_1+\gamma_1z_0,\\
 w_2&=z_2+\gamma_2z_0,\\
 w_3&=z_3+\gamma_1z_2+\gamma_2z_1+\gamma_3z_0.
 \end{aligned}
\]

Apply the zeta transform again to obtain \(Gx\). There are exactly five scalar multiplication sites. The emission block's all-ones correction adds \(z_0\) to \(w_3\), so it needs no extra product.

Batch the four cosets into the four 128-bit record lanes of a ZMM register. Each multiplication site is then one GFNI512 instruction, for five per map. The state-to-output map shares its initial four-point transform across all four cosets. The input-to-state map sums their results. No general GF128 multiplication or dense state-basis conversion appears.

The scalar four-point circuit uses four XORs for each zeta transform and five XORs between them. The emission correction adds one XOR. A packed implementation can reduce the four cosets before the final injection transform. Record transposes, broadcasts, horizontal reductions, and XOR instructions still matter; the multiplication count does not account for their latency or register pressure.

## Complete arithmetic comparison

Use the existing outer cost of 41 GFNI512 products per four GRS[32,22] rows, hence \(41/88\) per source record. The following conservative packet count applies to 128-bit records:

| Operation on sixteen encoded records | GFNI512 products |
|---|---:|
| Fixed emission map | 5 |
| Fixed injection map | 5 |
| Scale all four emission-state records by \(\eta\) | 1 |
| Fresh four-coordinate state diagonal | 1 |
| Sixteen input diagonals | 4 |
| Total | 16 |

The phase offsets are already in the fixed circulant coefficients. There is no additional per-packet state translation. The two fixed maps share three packed coefficient vectors per phase; changing from injection to emission only adds the all-ones correction. Four phases require at most 768 bytes for those three-vector tables.

The new random state diagonal uses four coefficient bytes per packet, or \(N/4\) bytes. A full GF256 GL4 matrix per packet would use sixteen bytes before any packed expansion. Both estimates exclude the coordinate diagonals and regional route tables, which remain unchanged. Any retained region-entry matrices contribute only the same \(O(1/L)\) product term as in the regional comparison, where \(L\) is the number of outer rows.

At \(N/K=32/22\), the full product-volume estimate is

\[
 \frac{41}{88}+\frac{16}{16}\frac{32}{22}
 =\boxed{\frac{169}{88}\approx1.92045}.
\]

For comparison, the prior eta-scaled regional estimate is \(181/88\). A full GL4 refresh every packet raised that estimate to \(213/88\). The present count is 6.63% below the former and 20.66% below the latter. These comparisons concern product volume only. Packing the four state records, coefficient expansion, and register lifetimes still require a native implementation audit before any speed claim.

## Involution, adjoint, and proof scope

For fixed sampled state diagonals, the zero-initial-state inner map is an involution when final state is discarded. Indeed, \(CY=CX\) in each packet because \(CA=0\). A second pass therefore follows the same state trajectory, and cancels the same emitted vector. This statement concerns the inner map on \(X\); it does not say that the whole outer/route/input-diagonal encoder is an involution.

The extended adjoint is different. Given output adjoint \(z\) and terminal-state adjoint \(h\), it returns

\[
 x^*=z+C^Th,\qquad
 s^*=D_s(h+\eta A^Tz).
\]

The transpose applies the diagonal after accumulating the emission contribution. Reusing the forward packet as its adjoint is incorrect.

The new superregular emission permits the separate diagonal-refresh rank bound. That argument averages over the fresh state diagonal, not over a uniform GL4 matrix. The directed replay in [diagonal_refresh_certificate.py](theory/diagonal_refresh_certificate.py) currently covers \(0.11N\le W<0.24N\), with margins 202/1833/8535 bits at K16/K18/K20. The independent dense theorem covers \(W\ge0.24N\). Inputs with \(0<W<0.11N\) remain unresolved. None of these class margins is a whole-code distance certificate.

## Replay and next gate

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_diagonal_packet_reference.py --check
```

The reference checks all 256 scalar values at every input basis coordinate, for both maps and all four phases: 20,480 map comparisons. It also verifies the row-scaled Cauchy identity, \(CA=0\), every within-phase cubic pairing, 128 extended adjoints, and six stream involutions with matching final states. All passed. The tests establish circuit identities; the Cauchy determinant argument establishes superregularity for every relevant minor.

Next, sharpen and evaluate the remaining sparse class using the new diagonal-refresh law. A positive whole-proof gate should precede native implementation or runtime measurement. All earlier construction references remain unchanged.
