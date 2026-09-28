# A dense contribution without an inner-state envelope

Independent row shuffles permit a positive decomposition of the averaged
outer counting measure. Its dense part has an outward bound below
2^-125.9745 at K = 2^20 and output threshold 209715. This is a partial
first-moment bound, not coverage of an active-group occupancy range.

## Averaged row measure

There are L = 8192 BCH[256,128] rows. Let A_w be the number of codewords
of weight w. After an independent uniform coordinate permutation, the
expected number of messages yielding a specified row word x is

    mu(x) = A_w / binomial(256,w),   w = weight(x).

Thus mu is a counting measure of total mass 2^128, not a probability
distribution. Independence of the row permutations makes the complete
averaged outer measure mu^tensor L. The subsequent group placements,
lane permutations, and inner setup are independent of these row shuffles.

Fix those subsequent choices. Their complete linear map T is invertible:
the route is a permutation, and the inner is block triangular with
identity diagonal. The argument below holds for each fixed T.

Let Abar_w be the authenticated coefficientwise BCH caps. The code is
even. Let U be the uniform distribution on the 255-dimensional subspace
of even 256-bit words. Define c = 2^129 and a nonnegative measure D by

    D(x) = max(0, Abar_w/binomial(256,w) - 2^-126)

for even w, and D(x) = 0 for odd w. Then

    mu <= c U + D

pointwise. The residual mass is the exact rational number

    R = sum_even_w max(0, Abar_w - 2^-126 binomial(256,w)).

The authenticated caps give log2 R = 66.09646835664..., enclosed outward
by the verifier. The residual has weights 0, 38, 40, ..., 58, their
complements 198, 200, ..., 218, and 256. Thus most row-word mass is handled
by the uniform component; the residual retains the difficult low-weight
and nearly all-one words explicitly.

## Positive tensor expansion

Expand (c U + D)^tensor L. A term with j uniform factors has coefficient
c^j and L-j residual factors. Fix a choice of those row positions and fix
the residual words. The input is uniform on an affine subspace of
dimension 255j. Invertibility of T preserves that dimension.

Any binary affine subspace of dimension d has d coordinate positions
forming an information set. Its projection to these positions is uniform
on all d-bit strings, including when the subspace has a nonzero offset.
For 0 < z <= 1, dropping the other output coordinates therefore gives

    E[z^weight(T(input))] <= ((1+z)/2)^(255j).

Integrating the residual words multiplies this bound by R^(L-j).
Summing the binomial(L,j) position choices and applying the nonnegative
Chernoff bound at threshold H = 209715 gives

    U_j(z) = binomial(L,j) c^j R^(L-j)
             z^(-H) ((1+z)/2)^(255j).

For 255j > 2H, the verifier uses the exact rational witness
z = H/(255j-H). Otherwise it uses z = 1. Every arithmetic operation in
the final evaluation is enclosed by Arb; the accumulated upper is rounded
outward after each addition.

Here j counts factors in a positive measure expansion. It is neither
the number of active BCH rows nor the number q of active four-row groups.
The expansion terms overlap as measures on input words. They must not be
reported as disjoint classes of original messages.

## Verified dense contribution

Independent runs at 192-bit and 384-bit precision establish

    sum_(j=7878)^8192 U_j
      < 1.196428403684e-38
      < 2^(-125.9745).

The bound holds for every fixed subsequent route and invertible inner,
so averaging their setup randomness preserves it. Adding j = 7877 does
not close: the computed aggregate upper has log2 value +30.20254....
This is looseness of that upper bound, not evidence of a bad message.

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_uniform.py --precision 192
python -B research/workstreams/permutation_locality/independent_rows/dense_uniform.py --precision 384
```

Both runs authenticate the existing 163 BCH dependency files and pass
988 exact small-code checks. The checks include pointwise domination,
affine shifts, information-set bounds, and the positive tensor identity.
The argument above, rather than these small examples, establishes the
general inequalities.

The components with j < 7878 remain. The original all-zero message must
also be excluded when bounding that remaining contribution. The upper
measure contains the zero input, so blindly summing its entire low-j
contribution cannot establish a small nonzero-message failure bound.

## A reusable offset lemma

The dense remainder can use a stronger fact than the information-set
bound. For a binary linear subspace V, any offset b, and 0 <= z <= 1,

    E_(v uniform V)[z^weight(v+b)]
        <= E_(v uniform V)[z^weight(v)].

To see this, let f(x) = z^weight(x) on n bits. Its normalized Fourier
coefficient at a binary character a is

    fhat(a) = 2^-n (1+z)^(n-weight(a)) (1-z)^weight(a) >= 0.

Averaging f over V+b leaves only characters in V-perp and multiplies
their coefficients by (-1)^(a dot b). Dropping those signs gives the
average over V. This proof also includes z = 0 and z = 1 by direct
evaluation or continuity.

Consequently, for a fixed placement of uniform-row components, arbitrary
residual words can be replaced by zero when upper-bounding their tilted
output moment. This does not make the full residual contribution zero:
its mass R^(L-j), position choices, and nonzero-message restriction remain.

## Next step

The information-set bound discards every output coordinate outside the
255j selected coordinates. It is therefore weak when j is much smaller
than 8192, even if the recursive inner spreads those random bits widely.
The next useful refinement is an inner bound for uniform bits in selected
rows, using the offset lemma to remove residual-word values from that
inner calculation. Group the selected row positions by how many of their
four rows are uniform. For r such rows in a group, replacing uniform-even
rows by uniform unrestricted rows costs a factor 2^r. Each packet then
has independent uniform bits on r randomly selected lanes, with weight
distribution binomial(r,b)/2^r. These are actual averaged input laws,
not entrywise maxima over incompatible packet shapes.

That refinement still needs a proof and outward implementation. In
particular, column independence holds only after the stated replacement
of the even-row distribution. The factor for conditioning all selected
rows back to even parity must be retained. No existing shared-shuffle
dense certificate is automatically a certificate for this new ensemble.

One exact feedback fact supports this direction. The existing census of
all 2^19 dual characters shows that a nonzero character vanishes on at
most 13 complete four-bit windows. Therefore the feedback columns in
any 14 windows span all 19 state coordinates. For j uniformly selected
windows, let h_0(a) count the windows annihilated by character a. Then

    Pr[feedback rank < 19]
      <= sum_(a nonzero) binomial(h_0(a),j)/binomial(32,j).

The bound is zero for j >= 14. At j = 8 it is approximately 0.00033304,
and at j = 12 approximately 6.2004e-8. Character orthogonality also gives
the exact average probability of zero feedback under independent uniform
bits in the selected windows: add one to that character sum and divide
by 2^19.

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_uniform.py --self-test-only --feedback-rank
```

Surjectivity does not justify treating the state as uniform after weighting
paths by their output. The output tilt and feedback remain correlated.
A next transfer can bound that correlation with a character polynomial
for biased independent input bits, using the same 9218 window histograms.
