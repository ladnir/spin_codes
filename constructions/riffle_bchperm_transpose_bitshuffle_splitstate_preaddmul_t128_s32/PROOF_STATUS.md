# Proof status

## Fixed-inner regular theorem

Let the message length be \(2^{20}\), and let the encoded length be
\(N=2^{21}\). Set \(d=188743=\lfloor0.09N\rfloor\).

Fix the maps \(A\) and \(B\) recorded in this folder. Assume the following
outer profile.

1. The binary [256,128] outer code has the recorded regular weight spectrum.
   Its regular density factor is at most two. Its only excluded codeword is
   the all-one word.

Setup samples the outer coordinate permutations, the region permutations,
and the epoch scalars independently. Under these assumptions, the expected
number of nonzero regular messages whose output weight is at most \(d\) is
at most

\[
2^{-80.8139}.
\]

Consequently, some fixed setup has no such regular message.

This theorem no longer assumes an ensemble of inner maps. The maps \(A\) and
\(B\) are fixed and exactly audited. The theorem remains conditional on the
modeled outer spectrum. Its numerical evaluation also uses nearest binary64
arithmetic.

## Fixed inner constituents

The fixed map \(B=[H\mid I_{32}]\) was selected with seed 3005876494. The 96
columns of \(H\) are distinct weight-three vectors. The complete kernel
spectrum was derived by enumerating all \(2^{32}\) dual words and applying
the MacWilliams transform. It has mass \(2^{96}\) and minimum distance four.
The exact nonactivation probability on each weight shell is no larger than
the ensemble target used by the regular recurrence. This holds for all
weights from one through 128. The largest fixed-to-target ratio is
0.8981027, at weight 38. Because the recurrence has nonnegative transfer
entries, substituting this fixed profile preserves the previous upper bound.

The fixed map \(A\) was selected with seed 2736783579. Exact enumeration of
all \(2^{32}\) image words gives minimum distance 28. All 128 coordinate
forms are nonzero and distinct, and no three coordinate forms are dependent.
The first three factorial moments of the enumerated spectrum agree exactly
with the values implied by these coordinate conditions. A direct matrix
audit also proves \(BA=0\). Thus the fixed map satisfies every property used
by the shell lemma, with distance 28 in place of the required distance 20.

## Fixed-code shell lemma

Fix \(0<z<1\). Condition one epoch input \(X\) to have weight \(h\).
The region permutation makes its support a uniform \(h\)-subset of the 128
epoch coordinates.

Fix a codeword \(c\in\operatorname{im}(A)\) of weight \(v\). Then
\(J=|\operatorname{supp}(X)\cap\operatorname{supp}(c)|\) is
hypergeometric, and

\[
K_h(v;z)
:=
\mathbb E_X[z^{\operatorname{wt}(X+c)}]
=
\sum_j
\frac{\binom vj\binom{128-v}{h-j}}{\binom{128}h}
z^{h+v-2j}.
\]

Let \(p_v\) be the weight distribution of a uniform nonzero word in
\(\operatorname{im}(A)\). The three coordinate conditions imply

\[
\sum_v p_v(v)_r
=
(128)_r\frac{2^{32-r}}{2^{32}-1},
\qquad r\in\{1,2,3\}.
\]

Minimum distance gives \(p_v=0\) for \(v<20\). Maximizing
\(\sum_v p_vK_h(v;z)\) subject to these constraints is a finite linear
program. Its optimum is a deterministic upper bound for the live-output
moment of every fixed map \(A\) with the stated profile.

This argument replaces the earlier random affine-coset assumption.

## State lemma for PreAddMul

Write \(D=2^{32}-1\). Suppose the incoming state is nonzero. Multiplication
by a fresh \(\alpha_i\) makes \(\alpha_iQ_i\) uniform on
\(\operatorname{GF}(2^{32})^*\), regardless of the distribution of
\(Q_i\).

If \(B(X_i)=0\), the next state is uniform and nonzero. If
\(B(X_i)\ne0\), termination has probability \(1/D\). Conditioned on no
termination, the next state is uniform on the nonzero field elements with
one value omitted.

For every nonnegative function on nonzero states, the punctured distribution
has expectation at most

\[
\frac{D}{D-1}
\]

times its uniform-nonzero expectation. The termination event depends on the
fresh scalar, while the emitted word does not. Therefore both outgoing
branches inherit the fixed-code shell moment. The termination branch pays
the additional factor \(1/D\).

The parent update does not have this property. There, termination fixes
\(Q_i=B(X_i)\), and the emitted word becomes the structured map
\((I+AB)X_i\).

## Numerical coverage

The regular recurrence conditions exactly on each occupation from 1 through
8192. It uses a separate Chernoff tilt for each occupation.

- The exact modeled one-active margin is 80.8139 bits.
- The two-active pointwise margin is 126.2135 bits.
- The 512-active pointwise margin is 25350.98 bits.
- The 1024-active pointwise margin is 39521.59 bits.
- After summing all regular occupations, the margin remains 80.8139 bits.

The all-one audits currently establish the following claims.

- Pure all-one checkpoints have at least 3472.65 bits of margin.
- Checkpoints with exactly one all-one block have at least 3472.65 bits.
- Every mixed occupation with total count at most 64 has 3985.99 bits of
  aggregate margin.

These calculations use nearest binary64 arithmetic. They are not
outward-rounded certificates.

## Remaining gates

The proof is structurally closed for regular outer words with fixed inner
maps. The remaining gates are:

1. replace the modeled outer spectrum with an explicit code and spectrum
   bound;
2. cover mixed occupations with at least two all-one words and total count
   above 64; and
3. rerun the final finite calculation with outward-rounded arithmetic.

The inner constituent, affine-coset, and termination gates are resolved for
PreAddMul. Using the exact distance-28 spectrum of \(A\) may improve the
parameters, but the current theorem does not need that improvement.

## Receipts

- `../riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/receipts/preaddmul_support_averaged_moment3_one_active.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/receipts/preaddmul_support_averaged_moment3_merged.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/receipts/preaddmul_all_one_only_checkpoints.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/receipts/preaddmul_one_all_one_checkpoints.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_t128_s32/receipts/preaddmul_mixed_low_total64.json`
- `receipts/fixed_b32_selection.json`
- `receipts/fixed_b32_dual_histogram.csv`
- `receipts/fixed_b32_kernel_spectrum.json`
- `receipts/fixed_b32_activation_domination.json`
- `receipts/fixed_b32_moment3_checkpoints.json`
- `receipts/fixed_a32_selection.json`
- `receipts/fixed_a32_histogram.csv`
- `receipts/fixed_a32_spectrum_audit.json`
