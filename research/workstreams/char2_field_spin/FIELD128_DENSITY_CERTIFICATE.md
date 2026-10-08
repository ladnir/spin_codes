# GRS[128,87]: directed middle-weight certificate

At K20, the unchanged width-sixteen, state-eight symmetric Cauchy inner with GRS[128,87] has a directed failure bound below `2^-408` for every input weight in `[.03N,.24N)`. This is a complete bound for that class, not a whole-code certificate. Matching sparse and dense bounds must be checked and summed separately.

The replay is [field128_density_certificate.py](theory/field128_density_certificate.py), with [eight tests](theory/field128_density_certificate_tests.py). It uses frozen rational witnesses and does not run an optimizer. It leaves the earlier GRS[256,176] certificate unchanged.

## Construction and exact interval

The field is GF256. Each packet has sixteen inputs and eight state coordinates. The symmetric maps, four packet phases, eta 2, and fresh independent nonzero input/state diagonals are those of [STATE8_DENSITY_CERTIFICATE.md](STATE8_DENSITY_CERTIFICATE.md). Routing consists of independent regional permutations applied in one global routing stage. The state starts at zero; its final value is discarded.

Only the systematic GRS outer changes. Its length, dimension, and redundancy are `b=128`, `k=87`, and `p=41`. Natural padding at requested dimension `2^20` gives

\[
 L=12096,\qquad K_{\rm pad}=87L=1052352,\qquad
 N=128L=1548288,\qquad e=\lceil N/10\rceil-1=154828.
\]

The padded rate is `87/128`, above `2/3`. Restricting the source to the requested dimension cannot introduce a new low-weight word.

Let `W` be the input support size before the inner encoder. The bounded event is existence of a nonzero source vector with

\[
 46449\le W\le371589
\]

and output weight at most `e`. Probability is over the ideal independent route and diagonal choices. The argument bounds existence in a source kernel; it is not a word first moment.

## Finite bound and directed replay

The packet matrix `T(x,y,tau)` is exactly the directed joint-return transfer from the earlier state-eight certificate. In particular, a return entry is the minimum of two complete covers. It does not multiply a return charge by a separately conditioned root charge. The local construction and all local assumptions are unchanged.

For `alpha=W/N`, the finite source/route comparison is now

\[
 B_W=\frac{(L+1)^{128}}{\tau(x\lambda)^W y^e}
 \left[(1-\alpha)^{128}
 +\tau^{-41}(1-\alpha+\tau\lambda\alpha)^{128}\right]^L
 e_0T(x,y,\tau)^{N/16}\mathbf1.
\]

The prefactor retains regional conditioning, and `tau^-1` is the nonzero-kernel rank threshold. The implementation explicitly uses outer length 128 and redundancy 41; it does not inherit the old values 256 and 80.

For each inclusive integer cell `[a,c]`, choose one frozen rational witness and retarget its source tilt at the cell midpoint. Retargeting preserves `lambda*d/(1-d)` from its anchor density `d`. It is an exact choice of a free parameter and does not interpolate probabilities.

Put `beta=tau*lambda`. The source bracket is uniformly bounded by

\[
 R_{a,c}=(1-a/N)^{128}
 +\tau^{-41}\bigl(1+(\beta-1)\alpha_*\bigr)^{128},
\]

where `alpha_*=c/N` if `beta>=1`, and `alpha_*=a/N` otherwise. Choose `W_*=a` if `x*lambda>=1`, and `W_*=c` otherwise. The sum over the cell is at most

\[
 (c-a+1)\frac{(L+1)^{128}R_{a,c}^{L}}
 {\tau(x\lambda)^{W_*}y^e}
 e_0T^{N/16}\mathbf1.
\]

Every operation uses outward-rounded eighty-digit decimal intervals. Exact `Fraction` values supply the witnesses, diagonal rank factors, and flat sums. Replay rejects floating witness fields. Explicit contiguity checks require the accepted cells to cover every integer weight exactly once, even with Python assertions disabled.

The eight rational anchors are at densities `.03,.06,.08,.10,.12,.14,.18,.24`; their quadruples are saved in `WITNESSES`. For each cell, the replay evaluates the three nearest anchors and takes the smallest directed upper endpoint. It starts with cells of at most 1549 weights and would bisect any cell whose upper endpoint exceeded `2^-256`.

## Result and checks

Run from the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/field128_density_certificate.py
python -B research/workstreams/char2_field_spin/theory/field128_density_certificate_tests.py
```

The complete replay accepted 210 cells, containing 1400 through 1549 weights each. It covered all 325141 weights and returned

\[
 \Pr[\text{bad word in the middle class}]
 \le 1.0889494960515043307856169234070111249533504489846018407873117305326180215123106
 \cdot10^{-123}<2^{-408}.
\]

The partition SHA-256 is `ea5af80dd068950da774d2daf8ed3cd5577588120652a17b7dc817952204fe83`. The receipt includes the outer parameters, field, width, state size, exact lengths, output threshold, covered endpoints, decimal bound, and dyadic exponent.

All eight tests pass. They independently check the local transfer with exact fractions; verify the length-128 conditioning prefactor and parity-41 source expression; exercise both endpoint choices; compare the original and generic floating objectives; verify rational retargeting; check gap-free refinement; and replay the whole interval.

The preceding floating search had positive margins at all eight anchor weights, ranging from 3517.94 to 363127.58 bits. Those samples only supplied rational witnesses. The directed interval replay, including cell variation and all integer-weight multiplicities, establishes the stated 408-bit class bound.

## Join required for a whole-code result

Every nonzero GRS[128,87] row has weight at least 42. Hence at least 1106 active rows imply `W>=1106*42=46452`, inside the certified interval or above it. A sparse certificate for all active-row counts 1 through 1105 therefore overlaps the middle class without a missing range.

The dense class must start at `ceil(.24N)=371590`, exactly one beyond this replay's upper endpoint. Its bound must use the same construction, padding, and threshold. Until those separate sparse and dense receipts are verified and summed, this file certifies only the middle class. No runtime claim follows from the proof.
