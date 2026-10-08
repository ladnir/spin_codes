# Higher rates with the cubic GF256 inner

The certified `[128,83]` construction is not close to an information-theoretic rate limit. At 10% relative distance, the GF256 Gilbert–Varshamov reference rate is about 0.8414. A rate-3/4 code therefore has room, but the current proof does not establish it.

This investigation keeps the outer length 128 and studies dimension 96. It uses requested `K=2^20`, padded dimension 1,050,624, and output length 1,400,832. The target excludes output weights through 140,083. The fixed inner has packet width eight, state dimension four, and the repaired cubic maps from `GF8_PACKET_MAPS.md`.

The setup distribution is unchanged: independent uniform regional permutations and independent nonzero GF256 input and state diagonals. All experiments here concern proof bounds. No native implementation or timing was changed.

## What the current bound loses

One active outer row is not the obstruction. At dimension 96, its existing bound still gives approximately 136.9 bits. The difficulty is intermediate numbers of active rows.

| Active outer rows | Existing two-state bound | Five-state bound at the same marks |
|---:|---:|---:|
| 16 | −325.4 bits | −307.8 bits |
| 64 | −2131.0 bits | −2057.8 bits |
| 128 | −4301.6 bits | −4158.9 bits |
| 256 | −8306.7 bits | −8032.0 bits |

Negative entries mean that this upper bound exceeds one. They do not provide a lower bound on failure probability or show that the code has low distance. The first column uses optimized floating-point marks. The second column replays their rational approximations with outward-rounded Arb arithmetic.

The five-state transfer retains a cap on the rank of the four incoming state forms. State zero denotes an actually zero witness state; caps one through four denote nonzero states. The transfer dominates nonnegative continuations that increase with the nonzero cap. It does not claim to sample a rank-transition distribution.

For a return with input weight `w`, selected inactive output zeros `z`, and incoming cap `c`, the canonical state-zero-flat rank is at most

\[
\min(w-1,c-1,3-z).
\]

The existing joint return equations and independently checked return minima still apply. Positive coefficient convolution and regional averaging preserve the continuation cone. This refinement removes some incompatible rank histories, but the table shows that it is not a substantial rate increase by itself.

## Returns matter more than state orientation

At the sparse optimum, the rank tilt is approximately six. The existing diagonal refresh already has nearly the same marginal root factors as a uniform `GL(4,256)` refresh:

| Selected roots | Diagonal bound, tilt 6 | Uniform-GL bound, tilt 6 |
|---:|---:|---:|
| 1 | 0.16669220 | 0.16667938 |
| 2 | 0.02778419 | 0.02777990 |
| 3 | 0.00462963 | 0.00462963 |

These factors alone give little reason to pay for a full state matrix. This comparison does not port the return proof to a changed refresh distribution. Joint return conditions would require a separate argument.

An optimistic diagnostic deletes every transition back to zero, including immediate zero-to-zero transitions with nonzero input. This deletion is not valid for the construction. It isolates how much the return terms cost the proof.

For 64 active rows, reoptimizing that diagnostic raises the score from approximately −2058 bits to +8863 bits. Its tilt moves from approximately 6.4 to almost 255. Thus, return paths prevent the sparse proof from using the field's full rank charge.

The diagnostic does not close every middle class automatically. At inverse-input density 0.03, its optimized middle bound still gives approximately −14,789 bits. At density 0.10, a bounded search found a positive score. Extending a sharper sparse cover may bridge this region; a failed middle relaxation alone does not rule that out.

The current evidence therefore supports improving short excursions and their proof. It does not support a claim that rate 3/4 is impossible, or that a more uniform state orientation alone will reach it.

## Two different local repairs

### Change the feedthrough map

Let `X` be the input after its diagonal, and let `R` be the refreshed state. The current maps satisfy `CA=0`. Instead of emitting `X+eta A R`, consider

\[
J=CX,\qquad Q'=R+J,\qquad Y=X+A(\eta R+\mu J).
\]

This replaces identity feedthrough by `I+mu AC`. The map remains invertible because `(AC)^2=0`. At a birth, the emitted word is `(I+mu AC)X`. At a return, it is `(I+(mu+eta)AC)X`.

For `eta=2` and `mu=4`, exact enumeration gives:

| Input weight | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Minimum birth weight | 8 | 6 | 5 | 5 | 3 | 2 | 2 | 1 |
| Minimum return weight | 8 | 6 | 6 | 4 | 4 | 2 | 2 | 1 |

These are local algebraic statements, not a distance certificate. In particular, active output equations no longer have the identity pivots used by the current proof.

The implementation must compute `C X` before the `A` pass. It therefore loses the existing paired five-GFNI512 map schedule. Two smaller passes may still be efficient, but they introduce a serial dependency. Native measurements are required before selecting this change.

The cheaper coefficient `mu=1` has birth minima `[8,4,3,2,4,3,2,1]`. It saves a multiplication but gives weaker birth spreading than `mu=4`.

### Mix pairs before the input diagonal

A separate candidate keeps the existing inner and its paired map schedule:

\[
\text{regional route}\;\longrightarrow\;
\text{fixed pairwise MDS mixing}\;\longrightarrow\;
\text{nonzero input diagonal}\;\longrightarrow\;
\text{current cubic inner}.
\]

There is no extra global permutation or pre-mix diagonal. A pair with one nonzero input produces two nonzero outputs. A pair with two nonzero inputs can produce one or two nonzero outputs, but cannot vanish entirely.

This directly improves singleton activation without serializing the two fixed inner maps. However, cancellation inside an occupied pair depends on the source forms. The proof must account for those mixed-coordinate zero subspaces. Treating each cancellation as an independent probability, or replacing every input weight `w` by `2w`, would be invalid.

The [pairwise pre-mix follow-up](GF8_PAIR_PREMIX.md) now has an audited local cover, but its first numerical gate is worse than the unmixed control. Keep the `mu=4` feedthrough change as a distinct unclosed candidate, rather than combining both changes before identifying their effects.

## Bounded follow-up: six state coordinates

A separate gate increased the state dimension to six while retaining packet width eight. This candidate uses

\[
A_{u,j}=(u+(128+j))^{-1},\qquad C=A^\mathsf{T},\qquad\eta=2,
\]

in the AES representation of GF256, with `u=0,...,7` and `j=0,...,5`. Addition here is field addition. Both matrices are superregular by the Cauchy determinant formula; the test also checked all 6004 square minors directly.

This candidate does not satisfy `CA=0`. That identity is not necessary for the local cover. The required replacement is invertibility of the extended packet map

\[
\begin{pmatrix}I_8&\eta A\\C&I_6\end{pmatrix},
\]

equivalently, invertibility of `I_6+eta CA`. The fixed candidate satisfies this condition. Its return map `I_8+eta AC` is also invertible and has exact return minima `[7,7,3,4,3,3,1,2]`.

The local proof ports for the following reasons. Superregular `A` gives the same diagonal root-potential bound, now with state dimension six. Superregular `C` makes a nonzero return require `h+w>6`, where `h` is the refreshed state's support size. Selecting `z` inactive output zeros requires `z<h`.

For `i` selected active output zeros, the joint return/output equation rank is at least

\[
\max\{6,\ i+z,\ h+w-v\},\qquad v=8-i-z.
\]

The last term follows by deleting the `v` unconstrained output rows from an invertible extended map. Substituting `h>=z+1` gives the same joint rank formula used by the existing cover. None of these steps requires `CA=0` once extended invertibility has been established.

At rate 3/4, the bounded screen returned:

| Class | Width-eight/state-six bound |
|---|---:|
| 64 active outer rows | −1313.7 bits |
| Inverse-input density 0.10 | −18241.9 bits |
| Inverse-input density 0.14 | −7926.7 bits |

All three numerical searches converged. Increasing the state helps, but this candidate still fails the gate. No broader search or native implementation was started. The files `theory/gf8_state6_packet_screen.py` and `theory/gf8_state6_packet_tests.py` retain the candidate and three passing tests.

A final bounded variant kept `A` and used `C[v,j]=(v+(64+j))^{-1}`, with `v=128,...,133` and `j=0,...,7`. It also satisfies extended invertibility at `eta=2`. Exact return minima become `[8,7,5,4,3,2,1,1]`: small input weights improve, but some larger weights weaken.

The asymmetric variant gives −1310.5 bits at 64 active rows, only three bits better. Its middle scores worsen to −19878.9 bits at density 0.10 and −12343.2 bits at density 0.14. All three searches converged. The tests check the algebra and all 6004 minors for both variants. This branch is parked; the `--asymmetric` option reproduces its gate.

The pairwise pre-mix also failed its first conservative gate; see its separate screen and audit. These negative results concern the tested covers, not impossibility results for either construction. The next useful design round should jointly strengthen activation and return behavior, rather than changing only outer size or state orientation. The feedthrough candidate above is one concrete starting point and still requires a new local proof.

## Reproduction and scope

The helper `theory/gf8_high_rate_inner_screen.py` implements the five-state transfer, sparse point replays, middle screens, diagnostic return deletion, and exact feedthrough branch profiles. Its output identifies diagnostics as non-rigorous and never claims a whole-code certificate.

From the `theory` directory:

```text
python -B -m unittest gf8_high_rate_inner_tests
python -B gf8_high_rate_inner_screen.py --dimension 96 --sparse 16 64 128 256 --optimize-control
python -B gf8_high_rate_inner_screen.py --dimension 96 --density .03 .10 .14 .20
python -B gf8_high_rate_inner_screen.py --feedthrough 1 4
python -B gf8_high_rate_inner_screen.py --gv
```

Seven tests check the continuation cone, comparison with the two-state cover, aggregate formulas, exact polynomial endpoints, branch minima, and diagnostic labels. Earlier whole-code certificates remain unchanged.
