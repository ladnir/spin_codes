# Coverage audit for GRS[256,176], state eight, width sixteen

The integer cutoffs cover every nonzero word when the sparse and intermediate bounds hold for this same construction. The dense class has an independent exact certificate with at least 12,072 bits of margin. The other replays have now passed; their [aggregate certificate](FIELD256_K20_CERTIFICATE.md) gives 119 whole-code bits. This note isolates the coverage argument and dense replay.

The construction uses the fixed systematic GRS[256,176] outer over AES GF256, independent row permutations in its 256 regions, and independent nonzero input-coordinate diagonals. The inner is the symmetric state-eight, width-sixteen Cauchy recurrence, with independent nonzero state diagonals in every packet. Its initial state is zero and its final state is discarded. This note does not cover the proposed width-thirty-two variant.

## Integer geometry and the complete partition

Pad the requested `K=2^20` source length to a multiple of `176*64`. The resulting parameters are

\[
 L=6016,\qquad K_{\mathrm{pad}}=1{,}058{,}816,\qquad
 N=256L=1{,}540{,}096.
\]

The forbidden output weights are `1,...,e`, where

\[
 e=\lceil N/10\rceil-1=154{,}009.
\]

For a nonzero outer word, let `A` be its number of nonzero constituent rows. Let `W` be its weight after routing and the nonzero input diagonals. Each nonzero GRS[256,176] row has weight at least 81. Routing and diagonals preserve weight, so

\[
 W\ge81A.
\]

The three proposed class bounds have these exact domains:

| Bound | Integer domain | Status in this audit |
|---|---|---|
| Sparse active-row bound | `1 <= A <= 600`, every input weight | Separate directed replay passed, 119 bits |
| Intermediate weight bound | `46,203 <= W <= 369,623` | Separate directed replay passed, 344 bits |
| Dense weight bound | `369,624 <= W <= N` | Exact replay passed |

Indeed, `ceil(.03N)=46,203` and `ceil(.24N)=369,624`. If `A>=601`, then `W>=48,681`, which exceeds the intermediate lower endpoint by 2,478. Such a word lies in the intermediate or dense class. If `A<=600`, it lies in the sparse class regardless of weight. Finally, `A=0` is the zero word and cannot produce a nonzero output.

Thus there is no gap between the stated classes. Their overlap is harmless. If the corresponding failure bounds are `B_sparse`, `B_middle`, and `B_dense`, then

\[
 \Pr[\text{some nonzero output has weight at most }e]
 \le B_{\mathrm{sparse}}+B_{\mathrm{middle}}+B_{\mathrm{dense}}.
\]

The forty-bit goal requires this sum to be at most `2^-40`; separate forty-bit bounds would not suffice. The cutoffs alone do not establish the numerical bounds.

## Exact dense certificate and the length lift

Condition on the inner map `B`. Identity current-input feedthrough makes `B` invertible. The state diagonals determining `B` are independent of the outer route and input diagonals. Therefore conditioning on `B` preserves the randomness required by the [regional dense theorem](REGIONAL_GLOBAL_ROUTE.md).

For a fixed inverse input of density `alpha`, that theorem gives, for any positive tilt `lambda`,

\[
 \Pr[x\text{ belongs to the routed, diagonally scaled outer code}]
 \le (L+1)^b\left[
 \lambda^{-b\alpha}\left((1-\alpha)^b+
 256^{-(b-k)}(1-\alpha+\tfrac{256}{255}\lambda\alpha)^b\right)
 \right]^L.
\]

The bound is uniform over all regional weight profiles. There is no additional union over profiles or input weights.

The old GRS[32,22] rational witnesses cover every `alpha` in `[.24,1]`. To lift each witness, put

\[
 X=(1-\alpha)^{32},\qquad
 Y=256^{-10}(1-\alpha+\tfrac{256}{255}\lambda\alpha)^{32}.
\]

The GRS[256,176] bracket is `X^8+Y^8`, which is at most `(X+Y)^8` because `X,Y>=0`. The tilt power also increases by eight. Hence the same per-coordinate exponent applies at the longer outer length. This argument uses the common rate and the same coefficient field; it does not require a tensor-product code.

The replay checks all 152 closed density intervals `[j/200,(j+1)/200]`, `48<=j<200`. It verifies the small inequality, the eighth-power comparison, and the large inequality using exact rational arithmetic. Floating optimization only proposes rational tilts; it cannot accept an invalid inequality. Every density in `[.24,1]`, including both endpoints, is covered. The resulting membership bound is

\[
 (L+1)^{256}2^{-1.28N}.
\]

The GF256 ball of radius `e` has at most `2^(1.27N)` vectors. The exact comparison used is

\[
 (10/9)^{100}2295^{10}<2^{127}.
\]

Summing over all low-weight outputs, including more than the necessary projective representatives, gives

\[
 B_{\mathrm{dense}}\le6017^{256}2^{-N/100}
 \le2^{-(\lfloor N/100\rfloor-256\cdot13)}
 =2^{-12072}.
\]

The proof holds for every independent invertible inner, so it applies after averaging over the state diagonals. The original all-size helper is not used to certify this profile at K16 or K18; this audit is specifically for K20.

## Scalar extension without another word union

All rank-cover and counting proofs can be carried out over `F=GF256`. Their common source quotient, required-live forms, diagonal samples, and equation ranks use that field. In particular, the diagonal pivot range is `1<=tau<=255`. No diagonal is assumed uniform over the nonzero elements of a larger field.

For a fixed setup, let `G` be the complete generator matrix over `F`. Let `E/F` be any finite extension and expand a message as `m=sum_j beta_j m_j` in an `F`-basis of `E`, with each `m_j` over `F`. Then

\[
 Gm=\sum_j\beta_jGm_j.
\]

If `Gm` is nonzero, at least one component `Gm_j` is nonzero. Every such component has support contained in `supp(Gm)`. Conversely, every base-field codeword is also an extension-field codeword. The minimum symbol distances are therefore equal for each fixed setup, not merely equal in distribution.

Consequently the complete bad-setup event is identical over GF256 and over each compatible extension `GF(2^m)` for `m=8,32,64,128`. The same bound holds simultaneously for all these record fields. There is no additional factor for extension coordinates or fields, and no union over larger-field messages. A concrete implementation must realize the specified embedded GF256 scalar action. The same argument applies when records are explicitly stored as vectors over GF256.

This extension argument applies after covering all base-field words. A component can have smaller weight than the extension-field word, so a partial density-class certificate must not be lifted in isolation. Likewise, the GF256 state-image moment must not be relabeled as a moment over an arbitrary larger-field state image.

## Dense replay and complete aggregate

From this directory run:

```text
python -B theory/field256_coverage_audit.py
```

Add `--show-tilts` to print every accepted rational witness. The helper checks the finite geometry, cutoff overlap, dense interval inequalities, and 12,072-bit margin. It deliberately reports `whole_code_certified=false` and does not evaluate the other two class bounds.

The [aggregate helper](theory/field256_whole_certificate.py) additionally executes the successful sparse and intermediate replays and sums all three bounds exactly. Its result is below `2^-119`. Native performance remains a separate question.
